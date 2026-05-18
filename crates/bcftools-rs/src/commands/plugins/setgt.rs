//! `bcftools +setGT` (upstream `bcftools/plugins/setGT.c`).
//!
//! Target classes `-t .` / `-t ./.` / `-t ./x` / `-t a` (missing/all
//! masks) and `-t q` (samples selected by the `-i` filter, evaluated
//! per-sample through the shared filter engine), with new-genotype
//! modes `-n 0` (reference) and `-n .` (missing). Per upstream
//! `set_gt`, every allele of a targeted sample genotype is replaced
//! (ploidy preserved, result unphased), so a partially-missing `./1`
//! becomes `0/0` under `-n 0`.
//!
//! Deferred to later slices (tracked in TODO.md): `-t q` with `-e`
//! (per-sample exclude invert), the sample-subset `GT[@file]` /
//! `binom()` forms (`setGT.{2,3}.out`), `-t X` random, `-n` major/minor
//! allele inference (`m`/`M`), custom `c:GT`, `-n i/p/u` phase ops,
//! `-n X` (VAF), and the `binom()` target.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
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

/// Which genotypes `-t` targets.
#[derive(Clone, Copy)]
enum Target {
    /// `.` / `./.` / `./x` / `a` — a missing/all class mask.
    Mask(TgtMask),
    /// `q` — samples selected by the `-i`/`-e` filter expression.
    Query,
}

/// Parse `-t`; returns `None` for the still-deferred classes
/// (`X` random / `binom`) so the caller can surface a clear error.
fn parse_target(spec: &str) -> Option<Target> {
    match spec {
        "q" | "?" => Some(Target::Query),
        "." => Some(Target::Mask(TgtMask {
            missing: true,
            partial: true,
            all: false,
        })),
        "./x" => Some(Target::Mask(TgtMask {
            partial: true,
            ..TgtMask::default()
        })),
        "./." => Some(Target::Mask(TgtMask {
            missing: true,
            ..TgtMask::default()
        })),
        "a" => Some(Target::Mask(TgtMask {
            all: true,
            ..TgtMask::default()
        })),
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
    /// `-i`/`-e` expression as `(exclude, expr)`; required for `-t q`.
    pub filter: Option<(bool, &'a str)>,
}

/// Reads the input VCF/BCF and returns the genotype-rewritten VCF text.
pub fn run(input: &Path, opts: Options<'_>) -> io::Result<(String, u64)> {
    let target = parse_target(opts.target).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!(
                "setGT -t '{}' is not supported in this slice (only ., ./., ./x, a, q)",
                opts.target
            ),
        )
    })?;
    let filter = match (&target, opts.filter) {
        (Target::Query, Some((false, expr))) => Some(expr),
        (Target::Query, Some((true, _))) => {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "setGT -t q with -e (exclude) is not supported in this slice",
            ));
        }
        (Target::Query, None) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Expected -i/-e with -t q",
            ));
        }
        (Target::Mask(_), Some(_)) => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Expected -t q with -i/-e",
            ));
        }
        (Target::Mask(_), None) => None,
    };
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
            let owned: Vec<String> = f.iter().map(|s| (*s).to_owned()).collect();
            for (si, sample) in f[9..].iter().enumerate() {
                let mut sub: Vec<String> = sample.split(':').map(str::to_owned).collect();
                if gi < sub.len() {
                    let do_set = match target {
                        Target::Mask(m) => mask_targets(&sub[gi], m),
                        Target::Query => {
                            sample_passes(&owned, si, filter.expect("query needs filter"))?
                        }
                    };
                    if do_set {
                        sub[gi] = rewrite_allele(&sub[gi], new, &mut nchanged);
                    }
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

/// Whether a mask target (`.`/`./.`/`./x`/`a`) selects this genotype.
fn mask_targets(gt: &str, tgt: TgtMask) -> bool {
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    let ploidy = alleles.len();
    if ploidy == 0 {
        return false;
    }
    let nmiss = alleles
        .iter()
        .filter(|a| **a == "." || a.is_empty())
        .count();
    tgt.all || (tgt.partial && nmiss > 0) || (tgt.missing && ploidy == nmiss)
}

/// Upstream `set_gt`: replace every allele with the new allele, ploidy
/// preserved, output unphased. Updates the changed-allele counter.
fn rewrite_allele(gt: &str, new: NewGt, nchanged: &mut u64) -> String {
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    let ploidy = alleles.len().max(1);
    let new_allele = match new {
        NewGt::Ref => "0",
        NewGt::Missing => ".",
    };
    let rebuilt = vec![new_allele; ploidy].join("/");
    if rebuilt != gt {
        *nchanged += alleles.iter().filter(|a| **a != new_allele).count() as u64;
    }
    rebuilt
}

/// Evaluate the `-i` expression for sample `si` (upstream `-t q`
/// per-sample `smpl_pass`): a single-sample [`EvalContext`] over that
/// sample's FORMAT values, falling back to record-level lookups for
/// site fields (CHROM/POS/QUAL/INFO/…).
fn sample_passes(fields: &[String], si: usize, expr: &str) -> io::Result<bool> {
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let sample = &fields[9 + si];
    let values: Vec<&str> = sample.split(':').collect();
    let context = EvalContext::new().with_sample(
        format_keys
            .iter()
            .enumerate()
            .map(|(i, key)| {
                let raw = values.get(i).copied().unwrap_or(".");
                let value = if key.eq_ignore_ascii_case("GT") {
                    FilterValue::String(raw.to_owned())
                } else if raw == "." || raw.is_empty() {
                    FilterValue::Missing
                } else if let Ok(n) = raw.parse::<f64>() {
                    FilterValue::Number(n)
                } else {
                    FilterValue::String(raw.to_owned())
                };
                ((*key).to_owned(), value)
            })
            .collect::<Vec<_>>(),
    );
    Ok(
        bcffilter::eval_expression_with(expr, &context, |name, sample_index| {
            if sample_index.is_some() {
                return None;
            }
            crate::commands::filter::record_lookup(name, fields)
        })?
        .truthy(),
    )
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
        let Target::Mask(m) = parse_target(t).unwrap() else {
            panic!("rw helper expects a mask target");
        };
        if !mask_targets(gt, m) {
            return None;
        }
        let mut c = 0;
        Some(rewrite_allele(gt, parse_new(n).unwrap(), &mut c))
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
        assert!(parse_target("X").is_none());
        assert!(parse_new("pM").is_none());
        assert!(parse_new("c:1").is_none());
        assert!(matches!(parse_target("q"), Some(Target::Query)));
    }

    #[test]
    fn query_per_sample_filter() {
        // FORMAT GT:GQ:DP; only the sample matching the expression sets.
        let fields: Vec<String> = "1\t3177144\t.\tG\tT\t.\t.\t.\tGT:GQ:DP\t./.:150:30\t0/1:99:30"
            .split('\t')
            .map(str::to_owned)
            .collect();
        let expr = r#"GT~"." && FMT/DP=30 && GQ=150"#;
        assert!(sample_passes(&fields, 0, expr).unwrap());
        assert!(!sample_passes(&fields, 1, expr).unwrap());
    }
}
