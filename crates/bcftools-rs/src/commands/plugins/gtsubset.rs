//! `bcftools +GTsubset` (upstream `bcftools/plugins/GTsubset.c`).
//!
//! Outputs only sites where every requested sample (`-s`) exclusively
//! shares a genotype: all selected samples must have the same GT and no
//! other sample may have it. Missing genotypes always pass. Phasing is
//! significant (the comparison is on the bcf-encoded allele/phase ints,
//! exactly as upstream).

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const MISSING: i64 = 0;
const VECTOR_END: i64 = -1;

/// Reads the input and returns the GTsubset-filtered VCF text.
pub fn run(input: &Path, samples: &str) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(&text, samples).map_err(io::Error::other)
}

/// Encodes one allele token to the htslib GT int: missing → 0, allele
/// `n` → `((n+1)<<1) | phase` (`phase` = 1 for a `|`-preceded allele).
fn encode(tok: &str, phased: bool) -> i64 {
    if tok == "." || tok.is_empty() {
        return MISSING;
    }
    match tok.parse::<i64>() {
        Ok(n) => ((n + 1) << 1) | (phased as i64),
        Err(_) => MISSING,
    }
}

/// Parses a sample GT subfield into encoded allele ints.
fn parse_gt(gt: &str) -> Vec<i64> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut phased = false; // separator preceding the current token
    let mut first = true;
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            out.push(encode(&cur, !first && phased));
            cur.clear();
            phased = ch == '|';
            first = false;
        } else {
            cur.push(ch);
        }
    }
    out.push(encode(&cur, !first && phased));
    out
}

fn compute(text: &str, samples_spec: &str) -> Result<String, String> {
    let lines: Vec<&str> = text.lines().collect();
    let samples: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();
    if samples.is_empty() {
        return Err("No samples in input file.".to_string());
    }
    let nsmp = samples.len();

    let sel_names: Vec<&str> = samples_spec.split(',').filter(|s| !s.is_empty()).collect();
    let mut selected = vec![false; nsmp];
    for name in &sel_names {
        let idx = samples
            .iter()
            .position(|s| s == name)
            .ok_or_else(|| format!("Sample '{name}' not in input vcf file."))?;
        selected[idx] = true;
    }

    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));

    let mut out = String::new();
    for line in &lines {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            if !has_pass && line.starts_with("##fileformat=") {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            }
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let gslot = f[8]
            .split(':')
            .position(|k| k == "GT")
            .ok_or_else(|| format!("GT not present at {}:{}", f[0], f[1]))?;

        // Per-sample encoded alleles; gte_smp is the max ploidy.
        let mut per: Vec<Vec<i64>> = Vec::with_capacity(nsmp);
        let mut gte_smp = 0usize;
        for s in 0..nsmp {
            let gt = f[9 + s].split(':').nth(gslot).unwrap_or(".");
            let enc = parse_gt(gt);
            gte_smp = gte_smp.max(enc.len());
            per.push(enc);
        }
        if gte_smp > 2 {
            return Err("GTsubset does not support ploidy higher than 2.".to_string());
        }
        // bcf_get_genotypes pads shorter samples with vector_end.
        let allele = |s: usize, j: usize| -> i64 {
            if j < per[s].len() {
                per[s][j]
            } else {
                VECTOR_END
            }
        };
        let second = |s: usize| -> i64 {
            if gte_smp == 2 {
                allele(s, 1)
            } else {
                VECTOR_END
            }
        };

        // a1/a2: GT of the first selected sample with both alleles set.
        let mut a1 = 0i64;
        let mut a2 = 0i64;
        let mut gt = -1i64;
        while a1 == 0 || a2 == 0 {
            gt += 1;
            if gt as usize == nsmp {
                break;
            }
            let g = gt as usize;
            if !selected[g] {
                continue;
            }
            a1 = allele(g, 0);
            a2 = second(g);
        }

        let mut pass = 0usize;
        for (i, &is_sel) in selected.iter().enumerate() {
            let b1 = allele(i, 0);
            let b2 = second(i);
            if b1 == MISSING || b2 == MISSING {
                pass += 1;
                continue;
            } else if is_sel {
                if b1 == a1 && b2 == a2 {
                    pass += 1;
                    continue;
                }
                break;
            } else {
                if b1 != a1 || b2 != a2 {
                    pass += 1;
                    continue;
                }
                break;
            }
        }
        if pass == nsmp {
            // bcftools expands a lone "." sample column to one missing
            // value per FORMAT subfield on write (e.g. "." -> ".:.:.").
            let nfmt = f[8].split(':').count();
            if nfmt > 1 && f[9..].contains(&".") {
                let dot = vec!["."; nfmt].join(":");
                let mut g: Vec<&str> = f.clone();
                for c in g.iter_mut().skip(9) {
                    if *c == "." {
                        *c = &dot;
                    }
                }
                out.push_str(&g.join("\t"));
            } else {
                out.push_str(line);
            }
            out.push('\n');
        }
    }
    Ok(out)
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
        ".bcftools-rs-gtsubset-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encoding_matches_htslib() {
        // "0/1" -> [unphased 0, unphased 1] = [2, 4]
        assert_eq!(parse_gt("0/1"), vec![2, 4]);
        // "0|1" -> [unphased 0, phased 1] = [2, 5]
        assert_eq!(parse_gt("0|1"), vec![2, 5]);
        // "./." -> [missing, missing]
        assert_eq!(parse_gt("./."), vec![0, 0]);
        // haploid "1" -> [4]
        assert_eq!(parse_gt("1"), vec![4]);
    }
}
