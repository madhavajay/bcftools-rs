//! `bcftools +dosage` (upstream `bcftools/plugins/dosage.c`).
//!
//! Prints per-sample allele dosage from `-t` tags (`PL`, `GL`, `GT`), trying
//! the handlers in order and using the first that applies. PL/GL dosages are
//! derived from the diploid GL-ordered genotype likelihoods
//! (`10^(-0.1*PL)` / `10^GL`, normalized) accumulated into per-allele
//! dosages; GT dosage is the alt-allele count. All likelihood arithmetic is
//! in `f32` to match upstream's `float` precision. Missing data yields `-1`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tag {
    Pl,
    Gl,
    Gt,
}

/// Reads the input VCF/BCF and returns the dosage table text.
pub fn run(input: &Path, tags: &[String]) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(&text, tags).map_err(io::Error::other)
}

fn compute(text: &str, tags: &[String]) -> Result<String, String> {
    let lines: Vec<&str> = text.lines().collect();

    // FORMAT tags declared in the header (for the PL/GL "present" gate).
    let mut has_pl = false;
    let mut has_gl = false;
    let mut samples: Vec<&str> = Vec::new();
    for l in &lines {
        if l.starts_with("##FORMAT=<ID=PL,") {
            has_pl = true;
        } else if l.starts_with("##FORMAT=<ID=GL,") {
            has_gl = true;
        } else if l.starts_with("#CHROM") {
            samples = l.split('\t').skip(9).collect();
        }
    }

    // Build the ordered handler list exactly like upstream `init`.
    let mut handlers: Vec<Tag> = Vec::new();
    for t in tags {
        match t.as_str() {
            "PL" => {
                if has_pl {
                    handlers.push(Tag::Pl);
                }
            }
            "GL" => {
                if has_gl {
                    handlers.push(Tag::Gl);
                }
            }
            "GT" => handlers.push(Tag::Gt),
            other => return Err(format!("No handler for tag \"{other}\"")),
        }
    }

    let mut out = String::new();
    out.push_str("#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT");
    for (i, s) in samples.iter().enumerate() {
        out.push_str(&format!("\t[{}]{s}", i + 5));
    }
    out.push('\n');

    for l in &lines {
        if l.starts_with('#') || l.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        let reference = f[3];
        let alts: Vec<&str> = if f[4] == "." {
            Vec::new()
        } else {
            f[4].split(',').collect()
        };
        let n_allele = 1 + alts.len();

        out.push_str(f[0]);
        out.push('\t');
        out.push_str(f[1]);
        out.push('\t');
        out.push_str(reference);
        if alts.is_empty() {
            out.push_str("\t.");
        } else {
            out.push('\t');
            out.push_str(&alts.join(","));
        }

        let nsmpl = if f.len() > 9 { f.len() - 9 } else { 0 };
        if n_allele == 1 {
            for _ in 0..nsmpl {
                out.push_str("\t0.0");
            }
            out.push('\n');
            continue;
        }

        let fmt = if f.len() > 8 { f[8] } else { "" };
        let sample_cols = if f.len() > 9 { &f[9..] } else { &[][..] };

        let mut printed = false;
        for &h in &handlers {
            if handler(h, fmt, sample_cols, n_allele, &mut out) {
                printed = true;
                break;
            }
        }
        if !printed {
            for _ in 0..nsmpl {
                out.push_str("\t-1.0");
            }
        }
        out.push('\n');
    }

    Ok(out)
}

/// Runs one handler. Returns `true` on success (dosages appended), `false`
/// if the tag does not apply to this record (try the next handler).
fn handler(tag: Tag, fmt: &str, samples: &[&str], n_allele: usize, out: &mut String) -> bool {
    match tag {
        Tag::Gt => {
            let Some(slot) = fmt.split(':').position(|k| k == "GT") else {
                return false;
            };
            for s in samples {
                let gt = s.split(':').nth(slot).unwrap_or(".");
                let mut dsg = vec![0.0f32; n_allele];
                let mut j = 0;
                for tok in gt.split(['/', '|']) {
                    if tok == "." || tok.is_empty() {
                        break; // missing
                    }
                    match tok.parse::<usize>() {
                        Ok(idx) if idx < n_allele => dsg[idx] += 1.0,
                        _ => break,
                    }
                    j += 1;
                }
                if j == 0 {
                    for d in dsg.iter_mut() {
                        *d = -1.0;
                    }
                }
                for (k, d) in dsg.iter().enumerate().skip(1) {
                    out.push(if k == 1 { '\t' } else { ',' });
                    out.push_str(&format!("{d:.1}"));
                }
            }
            true
        }
        Tag::Pl | Tag::Gl => {
            let key = if tag == Tag::Pl { "PL" } else { "GL" };
            let Some(slot) = fmt.split(':').position(|k| k == key) else {
                return false;
            };
            let expected = n_allele * (n_allele + 1) / 2;
            for s in samples {
                let sub = s.split(':').nth(slot).unwrap_or(".");
                let dsg = pl_gl_dosage(tag, sub, n_allele, expected);
                for (k, d) in dsg.iter().enumerate().skip(1) {
                    out.push(if k == 1 { '\t' } else { ',' });
                    out.push_str(&format!("{d:.6}"));
                }
            }
            true
        }
    }
}

/// One sample's PL/GL → per-allele dosage. Mirrors the upstream BRANCH:
/// any missing/short vector yields all `-1`.
fn pl_gl_dosage(tag: Tag, sub: &str, n_allele: usize, expected: usize) -> Vec<f32> {
    if sub == "." || sub.is_empty() {
        return vec![-1.0; n_allele];
    }
    let toks: Vec<&str> = sub.split(',').collect();
    let mut vals = vec![0.0f32; expected];
    let mut sum = 0.0f32;
    let mut j = 0;
    while j < expected {
        let Some(&t) = toks.get(j) else {
            break; // vector_end
        };
        if t == "." {
            break; // missing
        }
        let Ok(x) = t.parse::<f64>() else {
            break;
        };
        let v = match tag {
            Tag::Pl => 10f64.powf(-0.1 * x),
            _ => 10f64.powf(x),
        } as f32;
        vals[j] = v;
        sum += v;
        j += 1;
    }
    if j < expected {
        return vec![-1.0; n_allele];
    }
    if sum != 0.0 {
        for v in vals.iter_mut() {
            *v /= sum;
        }
    }
    vals[0] = 0.0;
    let mut dsg = vec![0.0f32; n_allele];
    let mut l = 0;
    for jj in 0..n_allele {
        for k in 0..=jj {
            dsg[jj] += vals[l];
            dsg[k] += vals[l];
            l += 1;
        }
    }
    dsg
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
        ".bcftools-rs-dosage-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.1\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"GL\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t3000150\t.\tC\tT\t59.2\tPASS\t.\tGT:PL:GL\t0/0:0,9,95:0,-0.9,-9.5\t.\n\
1\t3000152\t.\tC\tT,G\t59.2\tPASS\t.\tGT:PL:GL\t0/0:0,36,545,36,545,545:-0.00031839,-3.6121,-54.5012,-3.6121,-54.5012,-54.5012\t.\n";

    fn body(out: &str) -> Vec<&str> {
        out.lines().filter(|l| !l.starts_with('#')).collect()
    }

    #[test]
    fn pl_mode_matches_fixture_values() {
        let out = compute(VCF, &["PL".into()]).unwrap();
        assert_eq!(body(&out)[0], "1\t3000150\tC\tT\t0.111816\t-1.000000");
        assert_eq!(
            body(&out)[1],
            "1\t3000152\tC\tT,G\t0.000251,0.000251\t-1.000000,-1.000000"
        );
    }

    #[test]
    fn gl_mode_matches_fixture_values() {
        let out = compute(VCF, &["GL".into()]).unwrap();
        assert_eq!(body(&out)[0], "1\t3000150\tC\tT\t0.111816\t-1.000000");
    }

    #[test]
    fn gt_mode_matches_fixture_values() {
        let out = compute(VCF, &["GT".into()]).unwrap();
        assert_eq!(body(&out)[0], "1\t3000150\tC\tT\t0.0\t-1.0");
        assert_eq!(body(&out)[1], "1\t3000152\tC\tT,G\t0.0,0.0\t-1.0,-1.0");
    }

    #[test]
    fn header_lists_samples() {
        let out = compute(VCF, &["GT".into()]).unwrap();
        assert_eq!(
            out.lines().next().unwrap(),
            "#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT\t[5]A\t[6]B"
        );
    }
}
