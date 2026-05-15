//! `bcftools +fill-AN-AC` (upstream `bcftools/plugins/fill-AN-AC.c`).
//!
//! Fills `INFO/AN` (total called alleles) and `INFO/AC` (per-ALT allele
//! counts) from `FORMAT/GT`, over the text VCF path. Upstream is DEPRECATED
//! in favor of `+fill-tags` but is still exercised by `test_vcf_plugin`.
//! Header lines for `AC` and `AN` are inserted after the last existing
//! `##INFO` line, mirroring HTSlib `bcf_hdr_append` grouping.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const AC_HEADER: &str =
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">";
const AN_HEADER: &str = "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">";

/// Reads the input VCF/BCF and returns the AN/AC-filled VCF text.
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(fill(&text))
}

fn fill(text: &str) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let last_info = lines
        .iter()
        .rposition(|l| l.starts_with("##INFO="))
        .or_else(|| lines.iter().position(|l| l.starts_with("#CHROM")));

    let mut out = String::with_capacity(text.len() + 128);
    for (idx, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            if Some(idx) == last_info {
                if line.starts_with("##INFO=") {
                    out.push_str(AC_HEADER);
                    out.push('\n');
                    out.push_str(AN_HEADER);
                    out.push('\n');
                } else {
                    // Inserted right before #CHROM when no ##INFO lines exist;
                    // back out the just-written #CHROM, emit headers, re-emit.
                    out.truncate(out.len() - line.len() - 1);
                    out.push_str(AC_HEADER);
                    out.push('\n');
                    out.push_str(AN_HEADER);
                    out.push('\n');
                    out.push_str(line);
                    out.push('\n');
                }
            }
            continue;
        }
        if line.trim().is_empty() {
            out.push('\n');
            continue;
        }
        out.push_str(&fill_record(line));
        out.push('\n');
    }
    out
}

fn fill_record(line: &str) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return line.to_owned();
    }
    let alt = fields[4];
    let n_alt = if alt == "." || alt.is_empty() {
        0
    } else {
        alt.split(',').count()
    };

    let mut an: u64 = 0;
    let mut ac = vec![0u64; n_alt];
    if fields.len() > 9
        && let Some(gt_idx) = fields[8].split(':').position(|k| k == "GT")
    {
        for sample in &fields[9..] {
            let Some(gt) = sample.split(':').nth(gt_idx) else {
                continue;
            };
            for allele in gt.split(['/', '|']) {
                if allele == "." || allele.is_empty() {
                    continue;
                }
                let Ok(idx) = allele.parse::<usize>() else {
                    continue;
                };
                an += 1;
                if idx >= 1 && idx <= n_alt {
                    ac[idx - 1] += 1;
                }
            }
        }
    }

    let mut tags = format!("AN={an}");
    if n_alt > 0 {
        let ac_str: Vec<String> = ac.iter().map(|c| c.to_string()).collect();
        tags.push_str(&format!(";AC={}", ac_str.join(",")));
    }

    let info = strip_an_ac(fields[7]);
    let new_info = if info == "." || info.is_empty() {
        tags
    } else {
        format!("{info};{tags}")
    };

    let mut out: Vec<&str> = fields.clone();
    out[7] = new_info.as_str();
    out.join("\t")
}

fn strip_an_ac(info: &str) -> String {
    if info == "." {
        return info.to_owned();
    }
    info.split(';')
        .filter(|kv| {
            let key = kv.split_once('=').map(|(k, _)| k).unwrap_or(kv);
            key != "AN" && key != "AC"
        })
        .collect::<Vec<_>>()
        .join(";")
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
        ".bcftools-rs-fill-an-ac-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_an_ac_for_biallelic() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GQ\t0/1:9\t0/1:9";
        assert_eq!(
            fill_record(line),
            "1\t10\t.\tC\tT\t.\tPASS\tAN=4;AC=2\tGT:GQ\t0/1:9\t0/1:9"
        );
    }

    #[test]
    fn counts_per_alt_for_multiallelic_and_haploid() {
        let line = "1\t10\t.\tG\tT,C\t.\tPASS\tTEST=5\tGT:GQ\t0/1:9\t2:9";
        assert_eq!(
            fill_record(line),
            "1\t10\t.\tG\tT,C\t.\tPASS\tTEST=5;AN=3;AC=1,1\tGT:GQ\t0/1:9\t2:9"
        );
    }

    #[test]
    fn all_missing_yields_zero() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GQ\t./.:1\t./.:1";
        assert_eq!(
            fill_record(line),
            "1\t10\t.\tC\tT\t.\tPASS\tAN=0;AC=0\tGT:GQ\t./.:1\t./.:1"
        );
    }

    #[test]
    fn no_alt_emits_only_an() {
        let line = "1\t10\t.\tG\t.\t.\tPASS\t.\tGT:GQ\t./.:1\t./.:1";
        assert_eq!(
            fill_record(line),
            "1\t10\t.\tG\t.\t.\tPASS\tAN=0\tGT:GQ\t./.:1\t./.:1"
        );
    }

    #[test]
    fn existing_an_ac_are_replaced() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\tAN=99;DP=5;AC=99\tGT\t1/1";
        assert_eq!(
            fill_record(line),
            "1\t10\t.\tC\tT\t.\tPASS\tDP=5;AN=2;AC=2\tGT\t1/1"
        );
    }

    #[test]
    fn header_lines_inserted_after_last_info() {
        let vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##contig=<ID=1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t1\t.\tC\tT\t.\tPASS\t.\tGT\t1/1\n";
        let out = fill(vcf);
        let header_pos = out.find("##INFO=<ID=DP").unwrap();
        let ac_pos = out.find("##INFO=<ID=AC").unwrap();
        let an_pos = out.find("##INFO=<ID=AN").unwrap();
        let chrom_pos = out.find("#CHROM").unwrap();
        assert!(header_pos < ac_pos && ac_pos < an_pos && an_pos < chrom_pos);
        assert!(out.contains("\tAN=2;AC=2\tGT\t1/1\n"));
    }
}
