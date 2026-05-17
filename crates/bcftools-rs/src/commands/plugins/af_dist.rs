//! `bcftools +af-dist` (upstream `bcftools/plugins/af-dist.c` + `bin.c`).
//!
//! Collects the AF-deviation distribution and the HWE genotype-probability
//! distribution. Only non-reference genotypes with a known site allele
//! frequency are considered; per-genotype probabilities are `2*AF*(1-AF)`
//! (RA) and `AF**2` (AA). All binning arithmetic is done in `f32` to match
//! upstream's `float` precision exactly (the histogram edges are sensitive
//! to it).

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

/// Upstream default for both `-d`/`--dev-bins` and `-p`/`--prob-bins`.
pub const DEFAULT_BINS: &str = "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1";

/// Port of `bin.c`'s `_bin_t` for the `min..max` (here `0..1`) case.
struct Bins {
    bins: Vec<f32>,
}

impl Bins {
    /// Port of `bin_init(list_def, 0, 1)`.
    fn new(list_def: &str) -> Result<Bins, String> {
        let (min, max) = (0.0f32, 1.0f32);
        let list = if list_def.contains(',') {
            list_def.split(',').map(str::to_owned).collect()
        } else {
            read_bin_file(list_def)?
        };

        let mut bins: Vec<f32> = Vec::new();
        for tok in list {
            let v: f64 = tok
                .parse()
                .map_err(|_| format!("Could not parse {list_def}: {tok}"))?;
            let v = v as f32; // strtod -> float
            if min != max && (v < min || v > max) {
                return Err(format!(
                    "Expected values from the interval [{min},{max}], found {tok}"
                ));
            }
            bins.push(v);
        }
        if min != max {
            if bins.len() <= 1 {
                return Err(format!(
                    "Expected at least two bin boundaries in {list_def}"
                ));
            }
            let max_err = (bins[1] - bins[0]) * 1e-6;
            if (bins[0] - min).abs() > max_err {
                bins.insert(0, min);
            }
            let last = bins.len() - 1;
            if (bins[last] - max).abs() > max_err {
                bins.push(max);
            }
        }
        Ok(Bins { bins })
    }

    fn size(&self) -> usize {
        self.bins.len()
    }

    fn value(&self, idx: usize) -> f32 {
        self.bins[idx]
    }

    /// Faithful port of `bin_get_idx` (half-open `[)` binary search).
    fn idx(&self, value: f32) -> usize {
        let n = self.bins.len();
        if self.bins[n - 1] < value {
            return n - 1;
        }
        let mut imin: i64 = 0;
        let mut imax: i64 = n as i64 - 2;
        while imin < imax {
            let i = (imin + imax) / 2;
            let b = self.bins[i as usize];
            if value < b {
                imax = i - 1;
            } else if value > b {
                imin = i + 1;
            } else {
                return i as usize;
            }
        }
        if self.bins[imax as usize] <= value {
            imax as usize
        } else {
            (imin - 1) as usize
        }
    }
}

fn read_bin_file(list_def: &str) -> Result<Vec<String>, String> {
    let text =
        fs::read_to_string(list_def).map_err(|_| format!("Error: failed to read {list_def}"))?;
    Ok(text
        .lines()
        .filter(|line| !line.is_empty())
        .map(str::to_owned)
        .collect())
}

/// Reads the input VCF/BCF and returns the af-dist report text (including
/// the two `bcftools`-tagged provenance lines the harness strips with
/// `grep -v bcftools`).
pub fn run(
    input: &Path,
    af_tag: &str,
    dev_def: &str,
    prob_def: &str,
    list_def: Option<&str>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute_with_list(&text, af_tag, dev_def, prob_def, list_def).map_err(io::Error::other)
}

#[cfg(test)]
fn compute(text: &str, af_tag: &str, dev_def: &str, prob_def: &str) -> Result<String, String> {
    compute_with_list(text, af_tag, dev_def, prob_def, None)
}

fn compute_with_list(
    text: &str,
    af_tag: &str,
    dev_def: &str,
    prob_def: &str,
    list_def: Option<&str>,
) -> Result<String, String> {
    let dev_bins = Bins::new(dev_def)?;
    let prob_bins = Bins::new(prob_def)?;
    let list = list_def.map(parse_list_range).transpose()?;
    let mut dev_dist = vec![0u64; dev_bins.size()];
    let mut prob_dist = vec![0u64; prob_bins.size()];

    let mut out = String::new();
    // Provenance lines (stripped by the harness `grep -v bcftools`).
    out.push_str("# This file was produced by: bcftools +af-dist(bcftools-rs+htslib-rs)\n");
    out.push_str("# The command line was:\tbcftools +af-dist\n#\n");
    if let Some((min, max)) = list {
        out.push_str(&format!(
            "# GT, genotypes with P(AF) in [{:.6},{:.6}]; [2]Chromosome\t[3]Position[4]Sample\t[5]Genotype\t[6]AF-based probability\n",
            min as f64, max as f64
        ));
    }

    let mut nsmpl = 0usize;
    let mut samples_header: Vec<&str> = Vec::new();
    for line in text.lines() {
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            nsmpl = cols.len().saturating_sub(9);
            samples_header = if cols.len() > 9 {
                cols[9..].to_vec()
            } else {
                Vec::new()
            };
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if nsmpl == 0 {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let Some(af) = parse_af(f[7], af_tag) else {
            continue; // naf <= 0
        };

        let p_ra = 2.0f32 * af * (1.0f32 - af);
        let p_aa = af * af;
        let i_ra = prob_bins.idx(p_ra);
        let i_aa = prob_bins.idx(p_aa);
        let list_ra = list.is_some_and(|(min, max)| p_ra >= min && p_ra <= max);
        let list_aa = list.is_some_and(|(min, max)| p_aa >= min && p_aa <= max);

        let fmt = f[8];
        let gt_slot = fmt.split(':').position(|k| k == "GT");
        let samples = &f[9..];

        // ngt = max ploidy (vcf vector_end padding in upstream).
        let mut max_ploidy = 0usize;
        let parsed: Vec<Vec<&str>> = samples
            .iter()
            .map(|s| {
                let gt = match gt_slot {
                    Some(idx) => s.split(':').nth(idx).unwrap_or("."),
                    None => "",
                };
                let toks: Vec<&str> = gt.split(['/', '|']).collect();
                if toks.len() > max_ploidy {
                    max_ploidy = toks.len();
                }
                toks
            })
            .collect();

        let mut nals = 0i64;
        let mut nalt = 0i64;
        for (i, toks) in parsed.iter().enumerate() {
            let mut dosage = 0i64;
            let mut j = 0usize;
            while j < max_ploidy {
                let Some(&tok) = toks.get(j) else {
                    break; // vector_end
                };
                if tok == "." || tok.is_empty() {
                    break; // missing
                }
                let allele: i64 = match tok.parse() {
                    Ok(a) => a,
                    Err(_) => break,
                };
                if allele == 1 {
                    dosage += 1;
                }
                j += 1;
            }
            if j != max_ploidy {
                continue;
            }
            nals += j as i64;
            nalt += dosage;
            if dosage == 1 {
                prob_dist[i_ra] += 1;
                if list_ra {
                    out.push_str(&format!(
                        "GT\t{}\t{}\t{}\t1\t{:.6}\n",
                        f[0],
                        f[1],
                        samples_header.get(i).copied().unwrap_or("."),
                        p_ra as f64
                    ));
                }
            } else if dosage == 2 {
                prob_dist[i_aa] += 1;
                if list_aa {
                    out.push_str(&format!(
                        "GT\t{}\t{}\t{}\t2\t{:.6}\n",
                        f[0],
                        f[1],
                        samples_header.get(i).copied().unwrap_or("."),
                        p_aa as f64
                    ));
                }
            }
        }

        if nals != 0 && (nalt != 0 || af != 0.0) {
            let ratio = nalt as f32 / nals as f32;
            let af_dev = (af - ratio).abs();
            let i_af = dev_bins.idx(af_dev);
            dev_dist[i_af] += 1;
        }
    }

    out.push_str("# PROB_DIST, genotype probability distribution, assumes HWE\n");
    for (i, &count) in prob_dist.iter().enumerate().take(prob_bins.size() - 1) {
        out.push_str(&format!(
            "PROB_DIST\t{:.6}\t{:.6}\t{}\n",
            prob_bins.value(i) as f64,
            prob_bins.value(i + 1) as f64,
            count
        ));
    }
    out.push_str(&format!(
        "# DEV_DIST, distribution of AF deviation, based on {af_tag} and INFO/AN, AC calculated on the fly\n"
    ));
    for (i, &count) in dev_dist.iter().enumerate().take(dev_bins.size() - 1) {
        out.push_str(&format!(
            "DEV_DIST\t{:.6}\t{:.6}\t{}\n",
            dev_bins.value(i) as f64,
            dev_bins.value(i + 1) as f64,
            count
        ));
    }
    Ok(out)
}

fn parse_list_range(raw: &str) -> Result<(f32, f32), String> {
    let (min, max) = raw
        .split_once(',')
        .ok_or_else(|| format!("Could not parse: --list {raw}"))?;
    if max.contains(',') {
        return Err(format!("Could not parse: --list {raw}"));
    }
    let min = min
        .parse::<f32>()
        .map_err(|_| format!("Could not parse: --list {raw}"))?;
    let max = max
        .parse::<f32>()
        .map_err(|_| format!("Could not parse: --list {raw}"))?;
    Ok((min, max))
}

/// Returns the first value of INFO/`af_tag` as `f32`, or `None` when the tag
/// is absent (upstream `bcf_get_info_float` returning `naf <= 0`).
fn parse_af(info: &str, af_tag: &str) -> Option<f32> {
    if info == "." {
        return None;
    }
    for kv in info.split(';') {
        let mut it = kv.splitn(2, '=');
        let key = it.next()?;
        if key != af_tag {
            continue;
        }
        let val = it.next()?;
        let first = val.split(',').next()?;
        return first.parse::<f32>().ok();
    }
    None
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    if path == Path::new("-") {
        let tmp = stdin_tmp_path();
        let mut data = Vec::new();
        io::stdin().lock().read_to_end(&mut data)?;
        fs::write(&tmp, data)?;
        let result = read_vcf_text(&tmp);
        let _ = fs::remove_file(&tmp);
        return result;
    }

    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }

    let mut text = String::new();
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let file = File::open(path)?;
        let mut dec = MultiGzDecoder::new(file);
        dec.read_to_string(&mut text)?;
    } else {
        text = fs::read_to_string(path)?;
    }
    normalize_vcf_text(&mut text);
    Ok(text)
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-af-dist-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_bins_have_eleven_edges() {
        let b = Bins::new(DEFAULT_BINS).unwrap();
        assert_eq!(b.size(), 11);
        assert_eq!(b.value(0), 0.0);
        assert_eq!(b.value(10), 1.0);
    }

    #[test]
    fn bin_file_reads_one_boundary_per_line_and_inserts_extremes() {
        let tmp = tempfile::TempDir::new().unwrap();
        let path = tmp.path().join("bins.txt");
        fs::write(&path, "0.25\n0.5\n0.75\n\n").unwrap();

        let b = Bins::new(path.to_str().unwrap()).unwrap();
        assert_eq!(b.size(), 5);
        assert_eq!(b.value(0), 0.0);
        assert_eq!(b.value(1), 0.25);
        assert_eq!(b.value(2), 0.5);
        assert_eq!(b.value(3), 0.75);
        assert_eq!(b.value(4), 1.0);
    }

    #[test]
    fn bin_idx_matches_c_half_open_search() {
        let b = Bins::new(DEFAULT_BINS).unwrap();
        // exact edge returns that edge's index
        assert_eq!(b.idx(0.5), 5);
        // inside [0.4,0.5) -> index 4
        assert_eq!(b.idx(0.444444), 4);
        // 0 -> index 0
        assert_eq!(b.idx(0.0), 0);
        // above last edge -> n-1
        assert_eq!(b.idx(2.0), 10);
    }

    #[test]
    fn parse_af_first_value_or_none() {
        assert_eq!(parse_af("AN=6;AF=0.833333;AC=5", "AF"), Some(0.833333_f32));
        assert_eq!(
            parse_af("AN=6;AF=0.333333,0.333333", "AF"),
            Some(0.333333_f32)
        );
        assert_eq!(parse_af("AN=6", "AF"), None);
        assert_eq!(parse_af(".", "AF"), None);
    }

    #[test]
    fn small_record_prob_and_dev() {
        // af=0.5, two RA samples, one missing -> pRA=0.5 (bin 5), dev 0.
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
20\t326891\t.\tA\tAC\t999\t.\tAN=4;AF=0.5\tGT\t0|1\t0|1\t./.\n";
        let out = compute(vcf, "AF", DEFAULT_BINS, DEFAULT_BINS).unwrap();
        assert!(out.contains("PROB_DIST\t0.500000\t0.600000\t2\n"), "{out}");
        assert!(out.contains("DEV_DIST\t0.000000\t0.100000\t1\n"), "{out}");
    }

    #[test]
    fn list_range_emits_matching_genotypes_before_histograms() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
20\t326891\t.\tA\tAC\t999\t.\tAN=4;AF=0.5\tGT\t0|1\t1|1\t./.\n";
        let out =
            compute_with_list(vcf, "AF", DEFAULT_BINS, DEFAULT_BINS, Some("0.5,0.5")).unwrap();
        assert!(
            out.contains("# GT, genotypes with P(AF) in [0.500000,0.500000];"),
            "{out}"
        );
        assert!(out.contains("GT\t20\t326891\tA\t1\t0.500000\n"), "{out}");
        assert!(!out.contains("GT\t20\t326891\tB\t2\t0.250000\n"), "{out}");
        assert!(
            out.find("GT\t20\t326891\tA\t1\t0.500000\n").unwrap()
                < out.find("# PROB_DIST").unwrap(),
            "{out}"
        );
    }
}
