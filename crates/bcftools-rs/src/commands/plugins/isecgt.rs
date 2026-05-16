//! `bcftools +isecGT` (upstream `bcftools/plugins/isecGT.c`).
//!
//! Compares two VCF/BCF inputs and writes the first file with genotypes set
//! to missing where the corresponding genotype in the second file differs.
//! This local slice is text-backed and pairs records by CHROM/POS/REF/ALT.
//! Full synced-reader region/target/index semantics are deferred to the
//! shared synced-reader work.

use std::collections::{HashMap, VecDeque};
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::smpl_ilist::{self, SMPL_STRICT};
use crate::vcf_compat::normalize_vcf_text;

const MISSING: i32 = 0;
const VECTOR_END: i32 = i32::MIN;

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct Key {
    chrom: String,
    pos: String,
    ref_allele: String,
    alt: String,
}

#[derive(Clone, Debug)]
struct RecordGt {
    genotypes: Vec<Vec<i32>>,
    ploidy: usize,
}

/// Reads two inputs and returns file A with non-identical genotypes set to
/// missing.
pub fn run(path_a: &Path, path_b: &Path) -> io::Result<String> {
    let text_a = read_vcf_text(path_a)?;
    let text_b = read_vcf_text(path_b)?;
    compute(&text_a, &text_b).map_err(io::Error::other)
}

fn compute(text_a: &str, text_b: &str) -> Result<String, String> {
    let samples_a =
        header_samples(text_a).ok_or_else(|| "No #CHROM header in first file".to_string())?;
    let samples_b =
        header_samples(text_b).ok_or_else(|| "No #CHROM header in second file".to_string())?;
    let smpl = smpl_ilist::map(&samples_a, &samples_b, SMPL_STRICT).map_err(|e| e.to_string())?;
    let mut b_records = collect_b_records(text_b, samples_b.len())?;

    let mut out = String::new();
    for line in text_a.lines() {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let mut fields = split_record(line);
        if fields.len() < 9 + samples_a.len() {
            return Err(format!("Malformed VCF record in first file: {line}"));
        }
        let key = record_key(&fields)?;
        let Some(queue) = b_records.get_mut(&key) else {
            out.push_str(line);
            out.push('\n');
            continue;
        };
        let Some(rec_b) = queue.pop_front() else {
            out.push_str(line);
            out.push('\n');
            continue;
        };
        let rec_a = parse_record_gt(&fields, samples_a.len())?;
        if rec_a.ploidy != rec_b.ploidy {
            return Err(format!(
                "Different genotype ploidy at {}:{}",
                fields[0], fields[1]
            ));
        }

        let gt_slot = gt_slot(&fields[8])
            .ok_or_else(|| format!("GT not present at {}:{}", fields[0], fields[1]))?;
        let mut dirty = false;
        for (sample_idx, b_idx) in smpl.idx.iter().copied().enumerate() {
            if rec_a.genotypes[sample_idx] != rec_b.genotypes[b_idx] {
                dirty = true;
                fields[9 + sample_idx] = replace_gt(&fields[9 + sample_idx], gt_slot, rec_a.ploidy);
            }
        }
        if dirty {
            out.push_str(&fields.join("\t"));
        } else {
            out.push_str(line);
        }
        out.push('\n');
    }
    Ok(out)
}

fn header_samples(text: &str) -> Option<Vec<String>> {
    text.lines()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| line.split('\t').skip(9).map(str::to_owned).collect())
}

fn collect_b_records(
    text: &str,
    samples_count: usize,
) -> Result<HashMap<Key, VecDeque<RecordGt>>, String> {
    let mut out: HashMap<Key, VecDeque<RecordGt>> = HashMap::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields = split_record(line);
        if fields.len() < 9 + samples_count {
            return Err(format!("Malformed VCF record in second file: {line}"));
        }
        let key = record_key(&fields)?;
        let gt = parse_record_gt(&fields, samples_count)?;
        out.entry(key).or_default().push_back(gt);
    }
    Ok(out)
}

fn split_record(line: &str) -> Vec<String> {
    line.split('\t').map(str::to_owned).collect()
}

fn record_key(fields: &[String]) -> Result<Key, String> {
    if fields.len() < 5 {
        return Err("Malformed VCF record".to_string());
    }
    Ok(Key {
        chrom: fields[0].clone(),
        pos: fields[1].clone(),
        ref_allele: fields[3].clone(),
        alt: fields[4].clone(),
    })
}

fn gt_slot(format: &str) -> Option<usize> {
    format.split(':').position(|field| field == "GT")
}

fn parse_record_gt(fields: &[String], samples_count: usize) -> Result<RecordGt, String> {
    let slot = gt_slot(&fields[8])
        .ok_or_else(|| format!("GT not present at {}:{}", fields[0], fields[1]))?;
    let mut raw = Vec::with_capacity(samples_count);
    let mut ploidy = 0usize;
    for sample in fields.iter().skip(9).take(samples_count) {
        let gt = sample.split(':').nth(slot).unwrap_or(".");
        let values = parse_gt(gt);
        ploidy = ploidy.max(values.len());
        raw.push(values);
    }
    if ploidy == 0 {
        ploidy = 1;
    }
    for gt in &mut raw {
        gt.resize(ploidy, VECTOR_END);
    }
    Ok(RecordGt {
        genotypes: raw,
        ploidy,
    })
}

fn parse_gt(gt: &str) -> Vec<i32> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut phased = false;
    let mut first = true;
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            out.push(encode_gt_allele(&cur, !first && phased));
            cur.clear();
            phased = ch == '|';
            first = false;
        } else {
            cur.push(ch);
        }
    }
    out.push(encode_gt_allele(&cur, !first && phased));
    out
}

fn encode_gt_allele(raw: &str, phased: bool) -> i32 {
    if raw == "." || raw.is_empty() {
        return MISSING;
    }
    raw.parse::<i32>()
        .map(|allele| ((allele + 1) << 1) | i32::from(phased))
        .unwrap_or(MISSING)
}

fn replace_gt(sample: &str, gt_slot: usize, ploidy: usize) -> String {
    let mut fields: Vec<String> = sample.split(':').map(str::to_owned).collect();
    while fields.len() <= gt_slot {
        fields.push(".".to_string());
    }
    fields[gt_slot] = missing_gt(ploidy);
    fields.join(":")
}

fn missing_gt(ploidy: usize) -> String {
    if ploidy <= 1 {
        ".".to_string()
    } else {
        vec!["."; ploidy].join("/")
    }
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
        ".bcftools-rs-isecgt-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_encoding_preserves_phase() {
        assert_eq!(parse_gt("0/1"), vec![2, 4]);
        assert_eq!(parse_gt("0|1"), vec![2, 5]);
        assert_eq!(parse_gt("1|0"), vec![4, 3]);
        assert_eq!(parse_gt("./."), vec![0, 0]);
    }

    #[test]
    fn core_comparison_sets_only_different_samples_to_missing() {
        let a = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:8\t0/0:9
1\t2\t.\tG\tT\t.\t.\t.\tGT:DP\t1|0:7\t0/0:6
1\t3\t.\tC\tG\t.\t.\t.\tGT:DP\t0/1:5\t0/1:4
";
        let b = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tA
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/1
1\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0|1
";
        let out = compute(a, b).unwrap();
        assert!(out.contains("1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:8\t0/0:9\n"));
        assert!(out.contains("1\t2\t.\tG\tT\t.\t.\t.\tGT:DP\t./.:7\t./.:6\n"));
        assert!(out.contains("1\t3\t.\tC\tG\t.\t.\t.\tGT:DP\t0/1:5\t0/1:4\n"));
    }
}
