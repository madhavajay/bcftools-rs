//! `bcftools +fixploidy` (upstream `bcftools/plugins/fixploidy.c`).
//!
//! Rewrites FORMAT/GT so each sample's genotype has the ploidy implied
//! by a `CHROM FROM TO SEX PLOIDY` regions file (`-p`) and a
//! `NAME SEX` sample-sex file (`-s`), or a single forced ploidy (`-f`).
//! Mirrors the upstream `ploidy.c` query semantics. Only `-t GT` is
//! supported, matching upstream.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

/// One `CHROM FROM TO SEX PLOIDY` region (1-based inclusive).
struct PloidyRegion {
    chrom: String,
    from: u64,
    to: u64,
    sex: String,
    ploidy: i32,
}

/// Ported subset of `ploidy.c`: the regions plus the default ploidy.
struct Ploidy {
    regions: Vec<PloidyRegion>,
    dflt: i32,
    sex_defaults: HashMap<String, i32>,
}

impl Ploidy {
    fn from_file(text: &str, dflt: i32) -> Result<Ploidy, String> {
        let mut regions = Vec::new();
        let mut dflt = dflt;
        let mut sex_defaults = HashMap::new();
        for line in text.lines() {
            let s = line.trim();
            if s.is_empty() || s.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = s.split_whitespace().collect();
            if cols.len() < 5 {
                return Err(format!("Could not parse ploidy line: {line}"));
            }
            let ploidy: i32 = cols[4]
                .parse()
                .map_err(|_| format!("Could not parse: {line}"))?;
            if cols[0] == "*" && cols[1] == "*" && cols[2] == "*" {
                sex_defaults.insert(cols[3].to_string(), ploidy);
                continue;
            }
            let from: u64 = cols[1]
                .parse()
                .map_err(|_| format!("Could not parse: {line}"))?;
            let to: u64 = cols[2]
                .parse()
                .map_err(|_| format!("Could not parse: {line}"))?;
            regions.push(PloidyRegion {
                chrom: cols[0].to_string(),
                from,
                to,
                sex: cols[3].to_string(),
                ploidy,
            });
        }
        if let Some(star_default) = sex_defaults.get("*") {
            dflt = *star_default;
        }
        Ok(Ploidy {
            regions,
            dflt,
            sex_defaults,
        })
    }

    fn builtin_default(dflt: i32) -> Ploidy {
        let raw = "X 1 60000 M 1\n\
                   X 2699521 154931043 M 1\n\
                   Y 1 59373566 M 1\n\
                   Y 1 59373566 F 0\n\
                   MT 1 16569 M 1\n\
                   MT 1 16569 F 1\n";
        // The builtin string is well-formed by construction.
        Ploidy::from_file(raw, dflt).expect("builtin ploidy string")
    }

    /// Mirrors `ploidy_query`: returns `(sex_ploidy_for(sex), max_ploidy)`
    /// for `chrom:pos` (1-based), where `max_ploidy` only reflects the
    /// explicitly-listed non-default ploidies (upstream `_max`).
    fn query(&self, chrom: &str, pos: u64, sex: &str) -> (i32, i32) {
        let mut overlap = false;
        let mut sex_ploidy = self.dflt;
        let mut max_pl = -1i32;
        for r in &self.regions {
            if r.chrom == chrom && pos >= r.from && pos <= r.to {
                overlap = true;
                if r.ploidy != self.dflt {
                    if r.sex == sex {
                        sex_ploidy = r.ploidy;
                    }
                    if r.ploidy > max_pl {
                        max_pl = r.ploidy;
                    }
                }
            }
        }
        if !overlap {
            let sex_default = self.sex_defaults.get(sex).copied().unwrap_or(self.dflt);
            return (sex_default, self.dflt);
        }
        let max_ploidy = if max_pl == -1 { self.dflt } else { max_pl };
        (sex_ploidy, max_ploidy)
    }
}

/// Reads the input and returns the ploidy-fixed VCF text.
pub fn run(
    input: &Path,
    sex_file: Option<&Path>,
    ploidy_file: Option<&Path>,
    default_ploidy: i32,
    force_ploidy: Option<i32>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let sex_text = match sex_file {
        Some(p) => Some(fs::read_to_string(p)?),
        None => None,
    };
    let ploidy_text = match ploidy_file {
        Some(p) => Some(fs::read_to_string(p)?),
        None => None,
    };
    compute(
        &text,
        sex_text.as_deref(),
        ploidy_text.as_deref(),
        default_ploidy,
        force_ploidy,
    )
    .map_err(io::Error::other)
}

fn compute(
    text: &str,
    sex_text: Option<&str>,
    ploidy_text: Option<&str>,
    default_ploidy: i32,
    force_ploidy: Option<i32>,
) -> Result<String, String> {
    let ploidy = if force_ploidy.is_some() {
        None
    } else if let Some(pt) = ploidy_text {
        Some(Ploidy::from_file(pt, default_ploidy)?)
    } else {
        Some(Ploidy::builtin_default(default_ploidy))
    };

    // sample name -> sex; unlisted samples default to "F".
    let mut sample2sex: HashMap<String, String> = HashMap::new();
    if let Some(st) = sex_text {
        for line in st.lines() {
            let s = line.trim();
            if s.is_empty() || s.starts_with('#') {
                continue;
            }
            let mut it = s.split_whitespace();
            if let (Some(name), Some(sx)) = (it.next(), it.next()) {
                sample2sex.insert(name.to_string(), sx.to_string());
            }
        }
    }

    let lines: Vec<&str> = text.lines().collect();
    let samples: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();

    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));

    let mut out = String::new();
    for line in &lines {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            // bcftools emits a PASS filter header when writing.
            if !has_pass && line.starts_with("##fileformat=") {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let mut f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if f.len() < 10 {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let gt_slot = f[8].split(':').position(|k| k == "GT");
        let Some(gslot) = gt_slot else {
            // GT not present: emit unchanged (upstream `ngts<0`).
            out.push_str(line);
            out.push('\n');
            continue;
        };

        // Parse each sample's GT into alleles + separators.
        let mut parsed: Vec<(Vec<String>, Vec<char>)> = Vec::with_capacity(samples.len());
        let mut ngts = 0usize;
        for s in 0..samples.len() {
            let col = &f[9 + s];
            let gt = col.split(':').nth(gslot).unwrap_or(".");
            let (alleles, seps) = split_gt(gt);
            ngts = ngts.max(alleles.len());
            parsed.push((alleles, seps));
        }
        if ngts == 0 {
            out.push_str(line);
            out.push('\n');
            continue;
        }

        let pos: u64 = f[1].parse().unwrap_or(0);
        for (s, name) in samples.iter().enumerate() {
            let (per_ploidy, max_ploidy) = match (&ploidy, force_ploidy) {
                (_, Some(fp)) => (fp, fp),
                (Some(p), None) => {
                    let sx = sample2sex.get(*name).map(|x| x.as_str()).unwrap_or("F");
                    p.query(&f[0], pos, sx)
                }
                (None, None) => (default_ploidy, default_ploidy),
            };
            // Upstream leaves an already-haploid record untouched.
            if ngts >= max_ploidy as usize && ngts == 1 && max_ploidy == 1 {
                continue;
            }
            let (alleles, seps) = &parsed[s];
            let new_gt = rebuild_gt(alleles, seps, per_ploidy);
            let parts: Vec<&str> = f[9 + s].split(':').collect();
            let mut rebuilt: Vec<String> = parts.iter().map(|x| x.to_string()).collect();
            if gslot < rebuilt.len() {
                rebuilt[gslot] = new_gt;
            }
            f[9 + s] = rebuilt.join(":");
        }
        out.push_str(&f.join("\t"));
        out.push('\n');
    }
    Ok(out)
}

/// Splits a GT string into allele tokens and the separators preceding
/// alleles 1..n (`'/'` or `'|'`).
fn split_gt(gt: &str) -> (Vec<String>, Vec<char>) {
    let mut alleles = Vec::new();
    let mut seps = Vec::new();
    let mut cur = String::new();
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            alleles.push(std::mem::take(&mut cur));
            seps.push(ch);
        } else {
            cur.push(ch);
        }
    }
    alleles.push(cur);
    (alleles, seps)
}

/// Rebuilds a GT with exactly `ploidy` alleles (a single `.` for
/// ploidy 0), repeating the last allele/separator when extending.
fn rebuild_gt(alleles: &[String], seps: &[char], ploidy: i32) -> String {
    if ploidy <= 0 {
        return ".".to_string();
    }
    let p = ploidy as usize;
    let len = alleles.len();
    let mut out = String::new();
    out.push_str(alleles.first().map(|s| s.as_str()).unwrap_or("."));
    for k in 1..p {
        let sep = if k <= seps.len() {
            seps[k - 1]
        } else {
            *seps.last().unwrap_or(&'/')
        };
        let allele = if k < len {
            alleles[k].as_str()
        } else {
            alleles.last().map(|s| s.as_str()).unwrap_or(".")
        };
        out.push(sep);
        out.push_str(allele);
    }
    out
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
        ".bcftools-rs-fixploidy-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_split_join() {
        let (a, s) = split_gt("0|0");
        assert_eq!(a, vec!["0", "0"]);
        assert_eq!(s, vec!['|']);
        assert_eq!(rebuild_gt(&a, &s, 2), "0|0");
        assert_eq!(rebuild_gt(&a, &s, 1), "0");
        assert_eq!(rebuild_gt(&a, &s, 0), ".");
        assert_eq!(rebuild_gt(&a, &s, 5), "0|0|0|0|0");
    }

    #[test]
    fn gt_missing_unphased() {
        let (a, s) = split_gt("./.");
        assert_eq!(rebuild_gt(&a, &s, 4), "./././.");
    }

    #[test]
    fn ploidy_query_overlap_and_default() {
        let p = Ploidy::from_file("1 100 100 X 0\n1 100 100 Y 1\n1 100 100 Z 2\n", 2).unwrap();
        assert_eq!(p.query("1", 100, "X"), (0, 1));
        assert_eq!(p.query("1", 100, "Z"), (2, 1));
        // No overlap -> default for both.
        assert_eq!(p.query("1", 200, "X"), (2, 2));
    }

    #[test]
    fn ploidy_default_lines_apply_only_without_region_overlap() {
        let p = Ploidy::from_file("* * * M 1\n* * * F 2\n1 100 100 M 0\n", 2).unwrap();
        assert_eq!(p.query("1", 200, "M"), (1, 2));
        assert_eq!(p.query("1", 200, "F"), (2, 2));
        // Overlapping queries start from the global default, then apply
        // explicit region rows for matching sexes.
        assert_eq!(p.query("1", 100, "M"), (0, 0));
        assert_eq!(p.query("1", 100, "F"), (2, 0));
    }
}
