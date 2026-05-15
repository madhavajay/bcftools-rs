//! `bcftools +ad-bias` (upstream `bcftools/plugins/ad-bias.c`), report mode.
//!
//! For each sample/control pair from the `-s` file, finds the two most
//! frequent alleles from FORMAT/AD (scanning the sample then the control,
//! exactly as upstream) and runs Fisher's exact test on the 2x2
//! ref/alt × sample/control table, emitting `FT` lines below the threshold
//! and an `SN` summary. Fisher's exact test routes through
//! `htslib_rs::math::kt_fisher_exact` (the HTSlib `kfunc.c` port). The
//! `-c`/`--clean-vcf` VCF allele-removal output and `-v` variant-type
//! filtering need infrastructure tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::math::kt_fisher_exact;

use crate::vcf_compat::normalize_vcf_text;

/// One AD entry: an integer count, a missing value, or past the end of the
/// vector (shorter ploidy / fewer values than the record max).
#[derive(Clone, Copy, PartialEq)]
enum Ad {
    Int(i32),
    Missing,
    End,
}

/// Reads the inputs and returns the ad-bias report (with the two
/// `bcftools`-tagged provenance lines the harness strips via
/// `grep -v bcftools`).
pub fn run(
    input: &Path,
    samples_file: &Path,
    threshold: f64,
    min_dp: i32,
    min_alt_dp: i32,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let pairs_raw = fs::read_to_string(samples_file)?;
    Ok(compute(&text, &pairs_raw, threshold, min_dp, min_alt_dp))
}

struct Pair {
    smpl: usize,
    ctrl: usize,
    smpl_name: String,
    ctrl_name: String,
}

/// C `printf("%e", x)`: one mantissa digit, 6 fractional, signed exponent
/// with at least two digits.
fn fmt_e(x: f64) -> String {
    let s = format!("{x:.6e}"); // e.g. "9.254363e-4", "0.000000e0"
    let (mant, exp) = s.split_once('e').unwrap();
    let e: i32 = exp.parse().unwrap();
    let sign = if e < 0 { '-' } else { '+' };
    format!("{mant}e{sign}{:02}", e.abs())
}

fn compute(text: &str, pairs_raw: &str, threshold: f64, min_dp: i32, min_alt_dp: i32) -> String {
    // Resolve pairs against the #CHROM sample order.
    let mut sample_idx: Vec<&str> = Vec::new();
    for line in text.lines() {
        if line.starts_with("#CHROM") {
            sample_idx = line.split('\t').skip(9).collect();
            break;
        }
    }
    let idx_of = |name: &str| sample_idx.iter().position(|s| *s == name);

    let mut pairs: Vec<Pair> = Vec::new();
    for line in pairs_raw.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 2 {
            // upstream errors; in this slice be lenient and skip
            continue;
        }
        let (Some(smpl), Some(ctrl)) = (idx_of(cols[0]), idx_of(cols[1])) else {
            continue; // sample not in VCF -> skipped, like upstream
        };
        pairs.push(Pair {
            smpl,
            ctrl,
            smpl_name: cols[0].to_owned(),
            ctrl_name: cols[1].to_owned(),
        });
    }

    let mut out = String::new();
    out.push_str("# This file was produced by: bcftools +ad-bias(bcftools-rs+htslib-rs)\n");
    out.push_str("# The command line was:\tbcftools +ad-bias\n#\n");
    out.push_str(
        "# FT, Fisher Test\t[2]Sample\t[3]Control\t[4]Chrom\t[5]Pos\t[6]REF\t[7]ALT\t\
[8]smpl.nREF\t[9]smpl.nALT\t[10]ctrl.nREF\t[11]ctrl.nALT\t[12]P-value\n",
    );

    let mut nsite: u64 = 0;
    let mut ncmp: u64 = 0;

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let mut alleles: Vec<&str> = vec![f[3]];
        if f[4] != "." {
            alleles.extend(f[4].split(','));
        }
        if alleles.len() < 2 {
            continue;
        }

        let fmt = f[8];
        let Some(ad_slot) = fmt.split(':').position(|k| k == "AD") else {
            continue; // bcf_get_format_int32 AD -> nad<0
        };
        let sample_cols = &f[9..];

        // Per-sample AD vectors; nad = max length across samples.
        let mut ads: Vec<Vec<Ad>> = Vec::with_capacity(sample_cols.len());
        let mut nad = 0usize;
        for s in sample_cols {
            let adstr = s.split(':').nth(ad_slot).unwrap_or(".");
            let v: Vec<Ad> = if adstr == "." {
                vec![Ad::Missing]
            } else {
                adstr
                    .split(',')
                    .map(|t| match t.parse::<i32>() {
                        Ok(n) => Ad::Int(n),
                        Err(_) => Ad::Missing,
                    })
                    .collect()
            };
            if v.len() > nad {
                nad = v.len();
            }
            ads.push(v);
        }

        nsite += 1;

        let chrom = f[0];
        let pos = f[1];

        for pair in &pairs {
            let aptr = ad_at(&ads, pair.smpl, nad);
            let bptr = ad_at(&ads, pair.ctrl, nad);

            let Some((iref, ialt, n11, n12, n21, n22)) =
                pick_alleles(&aptr, &bptr, min_dp, min_alt_dp)
            else {
                continue;
            };
            ncmp += 1;

            let fisher = kt_fisher_exact(n11, n12, n21, n22).two_tail;
            if fisher >= threshold {
                continue;
            }
            out.push_str(&format!(
                "FT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                pair.smpl_name,
                pair.ctrl_name,
                chrom,
                pos,
                alleles[iref],
                alleles[ialt],
                n11,
                n12,
                n21,
                n22,
                fmt_e(fisher),
            ));
        }
    }

    out.push_str(
        "# SN, Summary Numbers\t[2]Number of Pairs\t[3]Number of Sites\t\
[4]Number of comparisons\t[5]P-value output threshold\n",
    );
    out.push_str(&format!(
        "SN\t{}\t{}\t{}\t{}\n",
        pairs.len(),
        nsite,
        ncmp,
        fmt_e(threshold)
    ));
    out
}

fn ad_at(ads: &[Vec<Ad>], sample: usize, nad: usize) -> Vec<Ad> {
    let v = &ads[sample];
    (0..nad)
        .map(|j| v.get(j).copied().unwrap_or(Ad::End))
        .collect()
}

/// Port of the upstream two-most-frequent-allele search (over the sample
/// then the control AD vectors) plus the depth/guard checks. Returns
/// `(iref, ialt, n11, n12, n21, n22)` or `None` to skip the pair.
fn pick_alleles(
    aptr: &[Ad],
    bptr: &[Ad],
    min_dp: i32,
    min_alt_dp: i32,
) -> Option<(usize, usize, i32, i32, i32, i32)> {
    let mut nbig: i32 = -1;
    let mut nsmall: i32 = -1;
    let mut ibig: i32 = -1;
    let mut ismall: i32 = -1;

    for (j, av) in aptr.iter().enumerate() {
        match *av {
            Ad::Missing => continue,
            Ad::End => break,
            Ad::Int(v) => {
                if ibig == -1 {
                    ibig = j as i32;
                    nbig = v;
                    continue;
                }
                if nbig < v {
                    if ismall == -1 || nsmall < nbig {
                        ismall = ibig;
                        nsmall = nbig;
                    }
                    ibig = j as i32;
                    nbig = v;
                    continue;
                }
                if ismall == -1 || nsmall < v {
                    ismall = j as i32;
                    nsmall = v;
                }
            }
        }
    }
    for (j, bv) in bptr.iter().enumerate() {
        match *bv {
            Ad::Missing => continue,
            Ad::End => break,
            Ad::Int(v) => {
                if ibig == -1 {
                    ibig = j as i32;
                    nbig = v;
                    continue;
                }
                if ibig == j as i32 {
                    if nbig < v {
                        nbig = v;
                    }
                    continue;
                }
                if nbig < v {
                    if ismall == -1 || nsmall < nbig {
                        ismall = ibig;
                        nsmall = nbig;
                    }
                    ibig = j as i32;
                    nbig = v;
                    continue;
                }
                if ismall == -1 || nsmall < v {
                    ismall = j as i32;
                    nsmall = v;
                }
            }
        }
    }
    if ibig == -1 || ismall == -1 {
        return None;
    }
    if nbig + nsmall < min_dp {
        return None;
    }

    let int_at = |p: &[Ad], i: i32| -> Option<i32> {
        match p[i as usize] {
            Ad::Int(v) => Some(v),
            _ => None,
        }
    };
    let a_big = int_at(aptr, ibig)?;
    let b_big = int_at(bptr, ibig)?;
    let a_small = int_at(aptr, ismall)?;
    let b_small = int_at(bptr, ismall)?;

    let (iref, ialt, nalt) = if ibig > ismall {
        (ismall as usize, ibig as usize, nbig)
    } else {
        (ibig as usize, ismall as usize, nsmall)
    };
    if nalt < min_alt_dp {
        return None;
    }

    let (n11, n12) = if iref == ibig as usize {
        (a_big, a_small)
    } else {
        (a_small, a_big)
    };
    let (n21, n22) = if iref == ibig as usize {
        (b_big, b_small)
    } else {
        (b_small, b_big)
    };
    Some((iref, ialt, n11, n12, n21, n22))
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
        ".bcftools-rs-ad-bias-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c_style_e_format() {
        assert_eq!(fmt_e(0.0009254363), "9.254363e-04");
        assert_eq!(fmt_e(1.328836e-09), "1.328836e-09");
        assert_eq!(fmt_e(1e-3), "1.000000e-03");
        assert_eq!(fmt_e(0.0), "0.000000e+00");
        assert_eq!(fmt_e(1.5e10), "1.500000e+10");
    }

    #[test]
    fn picks_ref_alt_like_upstream() {
        // smpl AD [4,3], ctrl AD [55,0]: ibig=0 (T), ismall=1 (C).
        let a = vec![Ad::Int(4), Ad::Int(3)];
        let b = vec![Ad::Int(55), Ad::Int(0)];
        let got = pick_alleles(&a, &b, 0, 1).unwrap();
        // iref=0, ialt=1, n11=4 n12=3 n21=55 n22=0
        assert_eq!(got, (0, 1, 4, 3, 55, 0));
    }

    #[test]
    fn fisher_matches_first_fixture_row() {
        // ad-bias.out row 1: 4 3 55 0 -> 9.254363e-04
        let p = kt_fisher_exact(4, 3, 55, 0).two_tail;
        assert_eq!(fmt_e(p), "9.254363e-04");
    }
}
