//! `bcftools +check-sparsity` (upstream `bcftools/plugins/check-sparsity.c`).
//!
//! Prints samples without enough genotype calls per chromosome or requested
//! region. This local slice is text-backed and supports VCF/VCF.gz/BCF input,
//! `-n`, `-r`, and `-R`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Debug, Eq, PartialEq)]
struct Region {
    label: String,
    chrom: String,
    start: Option<i64>,
    end: Option<i64>,
}

/// Reads the input and returns the check-sparsity report.
pub fn run(
    input: &Path,
    min_sites: usize,
    region: Option<&str>,
    region_file: Option<&Path>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let regions = load_regions(region, region_file)?;
    compute(&text, min_sites, &regions).map_err(io::Error::other)
}

fn compute(text: &str, min_sites: usize, regions: &[Region]) -> Result<String, String> {
    let samples = header_samples(text).ok_or_else(|| "No #CHROM header in input".to_string())?;
    if samples.is_empty() {
        return Ok(String::new());
    }
    if regions.is_empty() {
        return compute_by_chrom(text, &samples, min_sites);
    }
    let records = parse_records(text)?;
    let mut out = String::new();
    for region in regions {
        let mut counts = vec![0usize; samples.len()];
        let mut seen = false;
        for rec in records
            .iter()
            .filter(|rec| region.contains(&rec.chrom, rec.pos))
        {
            seen = true;
            update_counts(rec, &mut counts);
        }
        if seen {
            report_sparse(&mut out, &region.label, &samples, &counts, min_sites);
        }
    }
    Ok(out)
}

fn compute_by_chrom(text: &str, samples: &[String], min_sites: usize) -> Result<String, String> {
    let mut out = String::new();
    let mut current_chrom: Option<String> = None;
    let mut counts = vec![0usize; samples.len()];

    for rec in parse_records(text)? {
        if current_chrom.as_deref() != Some(&rec.chrom) {
            if let Some(chrom) = current_chrom.take() {
                report_sparse(&mut out, &chrom, samples, &counts, min_sites);
            }
            current_chrom = Some(rec.chrom.clone());
            counts.fill(0);
        }
        update_counts(&rec, &mut counts);
    }

    if let Some(chrom) = current_chrom {
        report_sparse(&mut out, &chrom, samples, &counts, min_sites);
    }
    Ok(out)
}

#[derive(Clone, Debug)]
struct ParsedRecord {
    chrom: String,
    pos: i64,
    gts: Vec<String>,
}

fn parse_records(text: &str) -> Result<Vec<ParsedRecord>, String> {
    let nsamples = header_samples(text)
        .ok_or_else(|| "No #CHROM header in input".to_string())?
        .len();
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 9 + nsamples {
            return Err(format!("Malformed VCF record: {line}"));
        }
        let pos = fields[1]
            .parse::<i64>()
            .map_err(|_| format!("Invalid VCF position: {}", fields[1]))?;
        let Some(gt_slot) = fields[8].split(':').position(|field| field == "GT") else {
            continue;
        };
        let gts = fields
            .iter()
            .skip(9)
            .take(nsamples)
            .map(|sample| sample.split(':').nth(gt_slot).unwrap_or(".").to_string())
            .collect::<Vec<_>>();
        out.push(ParsedRecord {
            chrom: fields[0].to_string(),
            pos,
            gts,
        });
    }
    Ok(out)
}

fn update_counts(rec: &ParsedRecord, counts: &mut [usize]) {
    for (idx, gt) in rec.gts.iter().enumerate() {
        if gt_has_first_allele(gt) {
            counts[idx] += 1;
        }
    }
}

fn gt_has_first_allele(gt: &str) -> bool {
    let first = gt.split(['/', '|']).next().unwrap_or(".");
    first != "." && !first.is_empty()
}

fn report_sparse(
    out: &mut String,
    label: &str,
    samples: &[String],
    counts: &[usize],
    min_sites: usize,
) {
    for (sample, count) in samples.iter().zip(counts) {
        if *count < min_sites {
            out.push_str(label);
            out.push('\t');
            out.push_str(sample);
            out.push('\n');
        }
    }
}

fn header_samples(text: &str) -> Option<Vec<String>> {
    text.lines()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| line.split('\t').skip(9).map(str::to_owned).collect())
}

fn load_regions(region: Option<&str>, region_file: Option<&Path>) -> io::Result<Vec<Region>> {
    let mut out = Vec::new();
    if let Some(raw) = region {
        for item in raw.split(',').filter(|item| !item.is_empty()) {
            out.push(parse_region(item)?);
        }
    }
    if let Some(path) = region_file {
        for line in fs::read_to_string(path)?.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            out.push(parse_region(line)?);
        }
    }
    Ok(out)
}

fn parse_region(raw: &str) -> io::Result<Region> {
    if raw.split_whitespace().count() >= 3 {
        let fields = raw.split_whitespace().collect::<Vec<_>>();
        let start = fields[1]
            .parse::<i64>()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid region start"))?
            + 1;
        let end = fields[2]
            .parse::<i64>()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid region end"))?;
        return Ok(Region {
            label: raw.to_string(),
            chrom: fields[0].to_string(),
            start: Some(start),
            end: Some(end),
        });
    }

    let (chrom, range) = raw.split_once(':').unwrap_or((raw, ""));
    if range.is_empty() {
        return Ok(Region {
            label: raw.to_string(),
            chrom: chrom.to_string(),
            start: None,
            end: None,
        });
    }
    let (start, end) = range.split_once('-').unwrap_or((range, range));
    let start = parse_pos(start)?;
    let end = parse_pos(end)?;
    Ok(Region {
        label: raw.to_string(),
        chrom: chrom.to_string(),
        start: Some(start),
        end: Some(end),
    })
}

fn parse_pos(raw: &str) -> io::Result<i64> {
    raw.replace(',', "")
        .parse::<i64>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "invalid region position"))
}

impl Region {
    fn contains(&self, chrom: &str, pos: i64) -> bool {
        if self.chrom != chrom {
            return false;
        }
        if let Some(start) = self.start
            && pos < start
        {
            return false;
        }
        if let Some(end) = self.end
            && pos > end
        {
            return false;
        }
        true
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
        ".bcftools-rs-check-sparsity-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t./.\t.
1\t2\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:5\t0/.:7\t.:9
2\t1\t.\tG\tT\t.\t.\t.\tGT\t./.\t0/0\t0/1
";

    #[test]
    fn reports_sparse_samples_by_chromosome() {
        let out = compute(VCF, 1, &[]).unwrap();
        assert_eq!(out, "1\tC\n2\tA\n");
    }

    #[test]
    fn min_sites_threshold_and_region_filter() {
        let regions = vec![parse_region("1:1-2").unwrap()];
        let out = compute(VCF, 2, &regions).unwrap();
        assert_eq!(out, "1:1-2\tB\n1:1-2\tC\n");
    }

    #[test]
    fn first_allele_only_controls_missing_status() {
        assert!(gt_has_first_allele("0/."));
        assert!(!gt_has_first_allele("./0"));
        assert!(!gt_has_first_allele("."));
    }
}
