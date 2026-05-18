//! `bcftools +setGT` (upstream `bcftools/plugins/setGT.c`).
//!
//! First slice: the filter-free target classes `-t .` / `-t ./.` /
//! `-t ./x` / `-t a` with the simple new-genotype modes `-n 0`
//! (reference) and `-n .` (missing). Per upstream `set_gt`, every
//! allele of a targeted sample genotype is replaced (ploidy preserved,
//! result unphased), so a partially-missing `./1` becomes `0/0` under
//! `-n 0`.
//!
//! Deferred to later slices (tracked in TODO.md): `-t q` (filter
//! engine), `-t X` random, `-n` major/minor allele inference (`m`/`M`),
//! custom `c:GT`, `-n i/p/u` phase ops, `-n X` (VAF), and the `binom()`
//! target. Those upstream rows stay deferred.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

/// Target-genotype class mask (upstream `tgt_mask`).
#[derive(Clone, Copy, Default)]
struct TgtMask {
    missing: bool, // `./.` — every allele missing
    partial: bool, // `./x` — at least one allele missing
    all: bool,     // `a`
}

#[derive(Clone, Copy)]
enum NewGt {
    Ref,
    Missing,
}

/// Parse `-t`; returns `None` for the still-deferred classes
/// (`q`/`X`/binom) so the caller can surface a clear error.
fn parse_target(spec: &str) -> Option<TgtMask> {
    match spec {
        "." => Some(TgtMask {
            missing: true,
            partial: true,
            all: false,
        }),
        "./x" => Some(TgtMask {
            partial: true,
            ..TgtMask::default()
        }),
        "./." => Some(TgtMask {
            missing: true,
            ..TgtMask::default()
        }),
        "a" => Some(TgtMask {
            all: true,
            ..TgtMask::default()
        }),
        _ => None,
    }
}

fn parse_new(spec: &str) -> Option<NewGt> {
    // Upstream: any '.' in the arg => GT_MISSING; "0" => GT_REF.
    if spec.contains('.') {
        return Some(NewGt::Missing);
    }
    match spec {
        "0" => Some(NewGt::Ref),
        _ => None,
    }
}

pub struct Options<'a> {
    pub target: &'a str,
    pub new_gt: &'a str,
}

/// Reads the input VCF/BCF and returns the genotype-rewritten VCF text.
pub fn run(input: &Path, opts: Options<'_>) -> io::Result<(String, u64)> {
    let tgt = parse_target(opts.target).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!(
                "setGT -t '{}' is not supported in this slice (only ., ./., ./x, a)",
                opts.target
            ),
        )
    })?;
    let new = parse_new(opts.new_gt).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!(
                "setGT -n '{}' is not supported in this slice (only 0, .)",
                opts.new_gt
            ),
        )
    })?;

    let text = read_vcf_text(input)?;
    let mut out = String::with_capacity(text.len());
    let mut nchanged: u64 = 0;

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let mut f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let gt_idx = f[8].split(':').position(|k| k == "GT");
        let mut rebuilt: Vec<String> = Vec::new();
        if let Some(gi) = gt_idx {
            for sample in &f[9..] {
                let mut sub: Vec<String> = sample.split(':').map(str::to_owned).collect();
                if gi < sub.len()
                    && let Some(newgt) = rewrite_gt(&sub[gi], tgt, new, &mut nchanged)
                {
                    sub[gi] = newgt;
                }
                rebuilt.push(sub.join(":"));
            }
            for (i, s) in rebuilt.iter().enumerate() {
                f[9 + i] = s;
            }
        }
        out.push_str(&f.join("\t"));
        out.push('\n');
    }
    Ok((out, nchanged))
}

/// Returns the rewritten GT string if this sample is targeted, else
/// `None`. Mirrors upstream `set_gt`: all alleles → the new allele,
/// ploidy preserved, output unphased.
fn rewrite_gt(gt: &str, tgt: TgtMask, new: NewGt, nchanged: &mut u64) -> Option<String> {
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    let ploidy = alleles.len();
    if ploidy == 0 {
        return None;
    }
    let nmiss = alleles
        .iter()
        .filter(|a| **a == "." || a.is_empty())
        .count();

    let do_set = tgt.all || (tgt.partial && nmiss > 0) || (tgt.missing && ploidy == nmiss);
    if !do_set {
        return None;
    }

    let new_allele = match new {
        NewGt::Ref => "0",
        NewGt::Missing => ".",
    };
    let rebuilt = vec![new_allele; ploidy].join("/");
    if rebuilt != gt {
        // upstream counts changed alleles, not samples
        *nchanged += alleles.iter().filter(|a| **a != new_allele).count() as u64;
    }
    Some(rebuilt)
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
        MultiGzDecoder::new(file).read_to_string(&mut text)?;
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
        ".bcftools-rs-setgt-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rw(gt: &str, t: &str, n: &str) -> Option<String> {
        let mut c = 0;
        rewrite_gt(gt, parse_target(t).unwrap(), parse_new(n).unwrap(), &mut c)
    }

    #[test]
    fn missing_to_ref() {
        assert_eq!(rw("./.", ".", "0").as_deref(), Some("0/0"));
        assert_eq!(rw(".", ".", "0").as_deref(), Some("0"));
        assert_eq!(rw(".|.", ".", "0").as_deref(), Some("0/0"));
    }

    #[test]
    fn present_gt_untouched() {
        assert_eq!(rw("0/1", ".", "0"), None);
        assert_eq!(rw("2", ".", "0"), None);
        assert_eq!(rw("1|1", ".", "0"), None);
    }

    #[test]
    fn partial_missing_targeted_by_dot_and_partial() {
        // `-t .` sets whenever any allele is missing -> all alleles reset.
        assert_eq!(rw("./1", ".", "0").as_deref(), Some("0/0"));
        // `-t ./x` (PARTIAL) fires on any missing allele (upstream
        // `GT_PARTIAL && nmiss`); a fully-present GT is untouched.
        assert_eq!(rw("0/1", "./x", "0"), None);
        assert_eq!(rw("./1", "./x", "0").as_deref(), Some("0/0"));
        // `-t ./.` only fires when fully missing.
        assert_eq!(rw("./1", "./.", "0"), None);
        assert_eq!(rw("./.", "./.", "0").as_deref(), Some("0/0"));
    }

    #[test]
    fn target_all_and_new_missing() {
        assert_eq!(rw("0/1", "a", "0").as_deref(), Some("0/0"));
        assert_eq!(rw("0/1", "a", ".").as_deref(), Some("./."));
        assert_eq!(rw("1", "a", ".").as_deref(), Some("."));
    }

    #[test]
    fn unsupported_modes_rejected() {
        assert!(parse_target("q").is_none());
        assert!(parse_new("pM").is_none());
        assert!(parse_new("c:1").is_none());
    }
}
