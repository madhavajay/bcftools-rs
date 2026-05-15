//! `bcftools +missing2ref` (upstream `bcftools/plugins/missing2ref.c`).
//!
//! Sets missing genotypes to the reference allele. This local slice
//! implements the default behavior (missing allele -> `0`) over the text
//! VCF path: every `.` allele token inside the `GT` FORMAT subfield is
//! rewritten to `0` while phase/unphase separators and all other FORMAT
//! subfields are preserved byte-for-byte. The `-e`/major-allele modes
//! remain tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

/// Reads the input VCF/BCF and returns the missing-to-ref rewritten VCF text.
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(rewrite(&text))
}

fn rewrite(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        out.push_str(&rewrite_record(line));
        out.push('\n');
    }
    out
}

fn rewrite_record(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() <= 9 {
        // No samples (or malformed) — nothing to do.
        return line.to_owned();
    }
    let format = fields[8];
    let Some(gt_idx) = format.split(':').position(|k| k == "GT") else {
        return line.to_owned();
    };

    let mut rewritten: Vec<String> = fields.iter().take(9).map(|s| s.to_string()).collect();
    for sample in &fields[9..] {
        rewritten.push(rewrite_sample(sample, gt_idx));
    }
    rewritten.join("\t")
}

fn rewrite_sample(sample: &str, gt_idx: usize) -> String {
    let mut parts: Vec<String> = sample.split(':').map(str::to_owned).collect();
    if gt_idx >= parts.len() {
        return sample.to_owned();
    }
    parts[gt_idx] = fix_gt(&parts[gt_idx]);
    parts.join(":")
}

/// Replaces every missing allele (`.`) in a GT string with `0`, leaving
/// `/` and `|` separators untouched. `./.`→`0/0`, `.|1`→`0|1`, `.`→`0`.
fn fix_gt(gt: &str) -> String {
    let mut out = String::with_capacity(gt.len());
    for ch in gt.chars() {
        if ch == '.' {
            out.push('0');
        } else {
            out.push(ch);
        }
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
        ".bcftools-rs-missing2ref-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fixes_diploid_and_haploid_missing() {
        assert_eq!(fix_gt("./."), "0/0");
        assert_eq!(fix_gt(".|."), "0|0");
        assert_eq!(fix_gt("./1"), "0/1");
        assert_eq!(fix_gt("1/."), "1/0");
        assert_eq!(fix_gt("."), "0");
        assert_eq!(fix_gt("0/1"), "0/1");
        assert_eq!(fix_gt("1|2"), "1|2");
    }

    #[test]
    fn rewrites_only_gt_subfield() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GQ:GL\t./.:245:-2.5,-5,-2.5";
        let got = rewrite_record(line);
        assert_eq!(
            got,
            "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GQ:GL\t0/0:245:-2.5,-5,-2.5"
        );
    }

    #[test]
    fn handles_gt_not_first_subfield() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGQ:GT\t245:./.";
        let got = rewrite_record(line);
        assert_eq!(got, "1\t10\t.\tC\tT\t.\tPASS\t.\tGQ:GT\t245:0/0");
    }

    #[test]
    fn preserves_headers_and_called_genotypes() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t1\t.\tC\tT\t.\tPASS\t.\tGT\t./.\t0/1\n";
        let out = rewrite(vcf);
        assert!(out.contains("##fileformat=VCFv4.2\n"));
        assert!(out.contains("\tGT\t0/0\t0/1\n"), "{out}");
    }

    #[test]
    fn record_without_samples_is_untouched() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.";
        assert_eq!(rewrite_record(line), line);
    }
}
