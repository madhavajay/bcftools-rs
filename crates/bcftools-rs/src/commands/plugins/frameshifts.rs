//! `bcftools +frameshifts` (upstream `bcftools/plugins/frameshifts.c`).
//!
//! Annotates indel alleles that overlap an exon with INFO/OOF values:
//! out-of-frame (1), in-frame (0), and not-applicable (-1). This local slice
//! supports plain BED-like exon files and text-backed VCF/BCF input.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const OOF_HEADER: &str = "##INFO=<ID=OOF,Number=A,Type=Integer,Description=\"Frameshift Indels: out-of-frame (1), in-frame (0), not-applicable (-1 or missing)\">";

#[derive(Clone, Debug)]
struct Exon {
    chrom: String,
    start: i64,
    end: i64,
}

/// Reads the input VCF/BCF and returns VCF text annotated with INFO/OOF.
pub fn run(input: &Path, exons_path: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let exons = read_exons(exons_path)?;
    compute(&text, &exons).map_err(io::Error::other)
}

fn compute(text: &str, exons: &[Exon]) -> Result<String, String> {
    let mut out = String::new();
    for line in text.lines() {
        if line.starts_with("##INFO=<ID=OOF,") {
            continue;
        }
        if line.starts_with("#CHROM") {
            out.push_str(OOF_HEADER);
            out.push('\n');
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields = annotate_record(line, exons)?;
        out.push_str(&fields);
        out.push('\n');
    }
    Ok(out)
}

fn annotate_record(line: &str, exons: &[Exon]) -> Result<String, String> {
    let mut fields = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
    if fields.len() < 8 {
        return Err(format!("Malformed VCF record: {line}"));
    }
    let pos0 = fields[1]
        .parse::<i64>()
        .map_err(|_| format!("Invalid VCF position: {}", fields[1]))?
        - 1;
    let ref_len = allele_len(&fields[3]);
    let alts = fields[4].split(',').collect::<Vec<_>>();
    if alts.is_empty() || fields[4] == "." {
        return Ok(line.to_string());
    }
    let diffs = alts
        .iter()
        .map(|alt| allele_len(alt) - ref_len)
        .collect::<Vec<_>>();
    if !diffs.iter().any(|diff| *diff != 0) {
        return Ok(line.to_string());
    }

    let most_negative = diffs.iter().copied().min().unwrap_or(0).min(0);
    let pos_to = if most_negative != 0 {
        pos0 - most_negative
    } else {
        pos0
    };
    let Some(exon) = first_overlap(exons, &fields[0], pos0, pos_to) else {
        return Ok(line.to_string());
    };

    let mut vals = Vec::with_capacity(diffs.len());
    for diff in diffs {
        if diff == 0 {
            vals.push(-1);
            continue;
        }
        let tlen = trimmed_exon_len(diff, pos0, exon);
        vals.push(if tlen == 0 {
            -1
        } else if tlen % 3 == 0 {
            0
        } else {
            1
        });
    }
    append_info(&mut fields[7], "OOF", &format_values(&vals));
    Ok(fields.join("\t"))
}

fn allele_len(raw: &str) -> i64 {
    if raw.starts_with('<') || raw == "*" {
        0
    } else {
        raw.len() as i64
    }
}

fn first_overlap<'a>(exons: &'a [Exon], chrom: &str, start: i64, end: i64) -> Option<&'a Exon> {
    exons
        .iter()
        .find(|exon| exon.chrom == chrom && exon.start <= end && exon.end > start)
}

fn trimmed_exon_len(diff: i64, pos0: i64, exon: &Exon) -> i64 {
    if diff > 0 {
        if exon.start <= pos0 && exon.end > pos0 {
            diff.abs()
        } else {
            0
        }
    } else if exon.start <= pos0 + diff.abs() {
        let mut tlen = diff.abs();
        if pos0 < exon.start {
            tlen -= exon.start - pos0 + 1;
        }
        if exon.end < pos0 + diff.abs() {
            tlen -= pos0 + diff.abs() - exon.end;
        }
        tlen.max(0)
    } else {
        0
    }
}

fn append_info(info: &mut String, tag: &str, value: &str) {
    let mut fields = info
        .split(';')
        .filter(|field| {
            if field.is_empty() || *field == "." {
                return false;
            }
            field.split_once('=').map(|(key, _)| key).unwrap_or(field) != tag
        })
        .map(str::to_owned)
        .collect::<Vec<_>>();
    fields.push(format!("{tag}={value}"));
    *info = fields.join(";");
}

fn format_values(values: &[i32]) -> String {
    values
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn read_exons(path: &Path) -> io::Result<Vec<Exon>> {
    let text = fs::read_to_string(path)?;
    let mut exons = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split_whitespace().collect::<Vec<_>>();
        if fields.len() < 3 {
            continue;
        }
        let start = fields[1]
            .parse::<i64>()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid exon start"))?;
        let end = fields[2]
            .parse::<i64>()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid exon end"))?;
        exons.push(Exon {
            chrom: fields[0].to_string(),
            start,
            end,
        });
    }
    Ok(exons)
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
        ".bcftools-rs-frameshifts-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotates_indel_frame_status() {
        let exons = vec![Exon {
            chrom: "1".to_string(),
            start: 9,
            end: 30,
        }];
        let text = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\t.\tA\tAT,ATGC,C\t.\t.\t.
1\t20\t.\tATGC\tA\t.\t.\tDP=1
1\t40\t.\tA\tAT\t.\t.\t.
";
        let out = compute(text, &exons).unwrap();
        assert!(out.contains("##INFO=<ID=OOF,"));
        assert!(out.contains("1\t10\t.\tA\tAT,ATGC,C\t.\t.\tOOF=1,0,-1\n"));
        assert!(out.contains("1\t20\t.\tATGC\tA\t.\t.\tDP=1;OOF=0\n"));
        assert!(out.contains("1\t40\t.\tA\tAT\t.\t.\t.\n"));
    }

    #[test]
    fn trims_deletions_to_exon_bounds() {
        let exon = Exon {
            chrom: "1".to_string(),
            start: 12,
            end: 15,
        };
        assert_eq!(trimmed_exon_len(-5, 9, &exon), 1);
    }
}
