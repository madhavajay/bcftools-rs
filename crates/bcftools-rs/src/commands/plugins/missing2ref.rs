//! `bcftools +missing2ref` (upstream `bcftools/plugins/missing2ref.c`).
//!
//! Sets missing genotypes to the reference or major allele. This local slice
//! rewrites missing `.` allele tokens inside the `GT` FORMAT subfield while
//! preserving all other FORMAT subfields byte-for-byte. Common `-i`/`-e`
//! record filters route through the shared text filter engine.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext};
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy)]
pub enum FilterMode {
    Include,
    Exclude,
}

#[derive(Clone, Copy)]
pub struct FilterSpec<'a> {
    pub mode: FilterMode,
    pub expr: &'a str,
}

/// Reads the input VCF/BCF and returns the missing-to-ref rewritten VCF text.
pub fn run(
    input: &Path,
    phased: bool,
    major: bool,
    filter: Option<FilterSpec<'_>>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    rewrite(&text, phased, major, filter)
}

fn rewrite(
    text: &str,
    phased: bool,
    major: bool,
    filter: Option<FilterSpec<'_>>,
) -> io::Result<String> {
    let mut out = String::with_capacity(text.len());
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if let Some(filter) = filter
            && !record_passes(&fields, filter)?
        {
            continue;
        }
        out.push_str(&rewrite_record(line, phased, major));
        out.push('\n');
    }
    Ok(out)
}

fn record_passes(fields: &[&str], filter: FilterSpec<'_>) -> io::Result<bool> {
    if fields.len() < 8 {
        return Ok(true);
    }
    let fields_owned: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
    let matched =
        bcffilter::eval_expression_with(filter.expr, &EvalContext::new(), |name, sample_index| {
            if sample_index.is_some() {
                return None;
            }
            crate::commands::filter::record_lookup(name, &fields_owned)
        })?
        .truthy();
    Ok(match filter.mode {
        FilterMode::Include => matched,
        FilterMode::Exclude => !matched,
    })
}

fn rewrite_record(line: &str, phased: bool, major: bool) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() <= 9 {
        // No samples (or malformed) — nothing to do.
        return line.to_owned();
    }
    let format = fields[8];
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else {
        return line.to_owned();
    };

    let replacement = if major {
        major_allele(&fields, gt_idx)
    } else {
        0
    };
    let mut rewritten: Vec<String> = fields.iter().take(9).map(|s| s.to_string()).collect();
    for sample in &fields[9..] {
        rewritten.push(rewrite_sample(sample, gt_idx, replacement, phased));
    }
    rewritten.join("\t")
}

fn rewrite_sample(sample: &str, gt_idx: usize, replacement: usize, phased: bool) -> String {
    let mut parts: Vec<String> = sample.split(':').map(str::to_owned).collect();
    if gt_idx >= parts.len() {
        return sample.to_owned();
    }
    parts[gt_idx] = fix_gt(&parts[gt_idx], replacement, phased);
    parts.join(":")
}

/// Replaces every missing allele (`.`) in a GT string. Phased mode mirrors
/// HTSlib's text writer by using `|` before alleles that were replaced.
fn fix_gt(gt: &str, replacement: usize, phased: bool) -> String {
    let mut alleles = Vec::new();
    let mut separators = Vec::new();
    let mut current = String::new();
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            alleles.push(std::mem::take(&mut current));
            separators.push(ch);
        } else {
            current.push(ch);
        }
    }
    alleles.push(current);

    let mut out = String::with_capacity(gt.len());
    for (i, allele) in alleles.iter().enumerate() {
        if i > 0 {
            let sep = if phased && allele == "." {
                '|'
            } else {
                separators.get(i - 1).copied().unwrap_or('/')
            };
            out.push(sep);
        }
        if allele == "." {
            out.push_str(&replacement.to_string());
        } else {
            out.push_str(allele);
        }
    }
    out
}

fn major_allele(fields: &[&str], gt_idx: usize) -> usize {
    let allele_count = if fields.get(4).is_none_or(|alt| *alt == ".") {
        1
    } else {
        fields[4].split(',').count() + 1
    };
    let mut counts = vec![0usize; allele_count];
    for sample in fields.iter().skip(9) {
        let Some(gt) = sample.split(':').nth(gt_idx) else {
            continue;
        };
        for allele in gt.split(['/', '|']) {
            if allele == "." || allele.is_empty() {
                continue;
            }
            if let Ok(idx) = allele.parse::<usize>()
                && let Some(count) = counts.get_mut(idx)
            {
                *count += 1;
            }
        }
    }
    counts
        .iter()
        .enumerate()
        .max_by(|(left_idx, left_count), (right_idx, right_count)| {
            left_count
                .cmp(right_count)
                .then_with(|| right_idx.cmp(left_idx))
        })
        .map(|(idx, _)| idx)
        .unwrap_or(0)
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
        ".bcftools-rs-missing2ref-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixes_diploid_and_haploid_missing() {
        assert_eq!(fix_gt("./.", 0, false), "0/0");
        assert_eq!(fix_gt(".|.", 0, false), "0|0");
        assert_eq!(fix_gt("./1", 0, false), "0/1");
        assert_eq!(fix_gt("1/.", 0, false), "1/0");
        assert_eq!(fix_gt(".", 0, false), "0");
        assert_eq!(fix_gt("0/1", 0, false), "0/1");
        assert_eq!(fix_gt("1|2", 0, false), "1|2");
    }

    #[test]
    fn phased_replacement_sets_separator_before_replaced_alleles() {
        assert_eq!(fix_gt("./.", 0, true), "0|0");
        assert_eq!(fix_gt("./1", 0, true), "0/1");
        assert_eq!(fix_gt("1/.", 0, true), "1|0");
        assert_eq!(fix_gt(".|1", 0, true), "0|1");
        assert_eq!(fix_gt("1/2", 0, true), "1/2");
    }

    #[test]
    fn rewrites_only_gt_subfield() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GQ:GL\t./.:245:-2.5,-5,-2.5";
        let got = rewrite_record(line, false, false);
        assert_eq!(
            got,
            "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GQ:GL\t0/0:245:-2.5,-5,-2.5"
        );
    }

    #[test]
    fn handles_gt_not_first_subfield() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGQ:GT\t245:./.";
        let got = rewrite_record(line, false, false);
        assert_eq!(got, "1\t10\t.\tC\tT\t.\tPASS\t.\tGQ:GT\t245:0/0");
    }

    #[test]
    fn preserves_headers_and_called_genotypes() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t1\t.\tC\tT\t.\tPASS\t.\tGT\t./.\t0/1\n";
        let out = rewrite(vcf, false, false, None).unwrap();
        assert!(out.contains("##fileformat=VCFv4.2\n"));
        assert!(out.contains("\tGT\t0/0\t0/1\n"), "{out}");
    }

    #[test]
    fn record_without_samples_is_untouched() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.";
        assert_eq!(rewrite_record(line, false, false), line);
    }

    #[test]
    fn major_mode_uses_most_common_called_allele_with_lowest_index_tie() {
        let line = "1\t10\t.\tC\tT,G\t.\tPASS\t.\tGT\t1/1\t0/1\t./.\t2/.";
        assert_eq!(
            rewrite_record(line, false, true),
            "1\t10\t.\tC\tT,G\t.\tPASS\t.\tGT\t1/1\t0/1\t1/1\t2/1"
        );
        let tie = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t./.";
        assert_eq!(
            rewrite_record(tie, false, true),
            "1\t10\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0"
        );
    }

    #[test]
    fn include_exclude_filters_select_records_before_rewrite() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t10\t.\tC\tT\t.\tPASS\tDP=5\tGT\t./.\n\
1\t20\t.\tC\tG\t.\tPASS\tDP=8\tGT\t./.\n";

        let out = rewrite(
            vcf,
            false,
            false,
            Some(FilterSpec {
                mode: FilterMode::Include,
                expr: "DP=5",
            }),
        )
        .unwrap();
        assert!(out.contains("1\t10\t.\tC\tT\t.\tPASS\tDP=5\tGT\t0/0\n"));
        assert!(!out.contains("1\t20\t"));

        let out = rewrite(
            vcf,
            false,
            false,
            Some(FilterSpec {
                mode: FilterMode::Exclude,
                expr: "DP=5",
            }),
        )
        .unwrap();
        assert!(!out.contains("1\t10\t"));
        assert!(out.contains("1\t20\t.\tC\tG\t.\tPASS\tDP=8\tGT\t0/0\n"));
    }
}
