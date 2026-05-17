//! `bcftools +split` (upstream `bcftools/plugins/split.c`).
//!
//! This local slice covers the filter-free upstream fixtures: default
//! per-sample splitting, `-S/--samples-file`, and `-G/--groups-file`, writing
//! VCF text files. `-k/--keep-tags` is supported for INFO/FORMAT projection.
//! `-i`/`-e` filters route through the shared text filter engine after each
//! output record is projected, matching upstream's per-output filtering.
//! Region/target restriction, BCF output, and indexing are deferred to the
//! shared reader/writer work.

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
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

#[derive(Clone, Debug)]
struct Subset {
    samples: Vec<usize>,
    names: Vec<String>,
    fname: String,
}

#[derive(Default)]
struct UniqueNames {
    seen: HashSet<String>,
}

impl UniqueNames {
    fn create(&mut self, raw: &str) -> String {
        let base = sanitize_name(raw);
        let mut name = base.clone();
        let mut id = 0usize;
        while self.seen.contains(&name) {
            id += 1;
            name = format!("{base}-{id}");
        }
        self.seen.insert(name.clone());
        name
    }
}

pub fn run(
    input: &Path,
    output_dir: &Path,
    samples_file: Option<&Path>,
    groups_file: Option<&Path>,
    keep_tags: Option<&str>,
    compress_vcf: bool,
    filter: Option<FilterSpec<'_>>,
) -> io::Result<()> {
    if samples_file.is_some() && groups_file.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "split accepts only one of --samples-file or --groups-file",
        ));
    }
    let text = read_vcf_text(input)?;
    let files = compute(&text, samples_file, groups_file, keep_tags, filter)?;
    fs::create_dir_all(output_dir)?;
    for (fname, content) in files {
        let fh = File::create(output_dir.join(with_vcf_suffix(&fname, compress_vcf)))?;
        if compress_vcf {
            let mut bgzf = htslib_rs::bgzf::io::Writer::new(fh);
            bgzf.write_all(content.as_bytes())?;
            bgzf.finish()?;
        } else {
            let mut fh = fh;
            fh.write_all(content.as_bytes())?;
        }
    }
    Ok(())
}

fn compute(
    text: &str,
    samples_file: Option<&Path>,
    groups_file: Option<&Path>,
    keep_tags: Option<&str>,
    filter: Option<FilterSpec<'_>>,
) -> io::Result<Vec<(String, String)>> {
    let parsed = ParsedVcf::new(text)?;
    let keep = KeepTags::parse(keep_tags);
    let subsets = match (samples_file, groups_file) {
        (Some(path), None) => subsets_from_samples_file(path, &parsed.samples)?,
        (None, Some(path)) => subsets_from_groups_file(path, &parsed.samples)?,
        (None, None) => default_subsets(&parsed.samples),
        (Some(_), Some(_)) => unreachable!(),
    };

    let mut out = Vec::with_capacity(subsets.len());
    for subset in subsets {
        let mut buf = String::new();
        for line in keep.filter_meta(&parsed.meta) {
            buf.push_str(line);
            buf.push('\n');
        }
        buf.push_str(&parsed.fixed_header);
        for name in &subset.names {
            buf.push('\t');
            buf.push_str(name);
        }
        buf.push('\n');
        for record in &parsed.records {
            let projected = project_record(record, &subset.samples, &keep);
            if let Some(filter) = filter
                && !record_passes_filter(&projected, filter)?
            {
                continue;
            }
            buf.push_str(&projected);
            buf.push('\n');
        }
        out.push((subset.fname, buf));
    }
    Ok(out)
}

fn record_passes_filter(record: &str, filter: FilterSpec<'_>) -> io::Result<bool> {
    let fields = record.split('\t').map(str::to_owned).collect::<Vec<_>>();
    if fields.len() < 8 {
        return Ok(true);
    }
    let context = record_context(&fields);
    let matched = bcffilter::eval_expression_with(filter.expr, &context, |name, sample_index| {
        if sample_index.is_some() {
            return None;
        }
        crate::commands::filter::record_lookup(name, &fields)
    })?
    .truthy();
    Ok(match filter.mode {
        FilterMode::Include => matched,
        FilterMode::Exclude => !matched,
    })
}

fn record_context(fields: &[String]) -> EvalContext {
    if fields.len() <= 9 {
        return EvalContext::new();
    }

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    fields[9..]
        .iter()
        .fold(EvalContext::new(), |context, sample| {
            context.with_sample(sample_values(&format_keys, sample))
        })
}

fn sample_values(format_keys: &[&str], sample: &str) -> Vec<(String, FilterValue)> {
    let values: Vec<&str> = sample.split(':').collect();
    format_keys
        .iter()
        .enumerate()
        .map(|(i, key)| {
            let raw = values.get(i).copied().unwrap_or(".");
            let value = if key.eq_ignore_ascii_case("GT") {
                FilterValue::String(raw.to_owned())
            } else {
                format_value(raw)
            };
            ((*key).to_owned(), value)
        })
        .collect()
}

fn format_value(raw: &str) -> FilterValue {
    if raw == "." || raw.is_empty() {
        return FilterValue::Missing;
    }
    if raw.contains(',') {
        return FilterValue::List(raw.split(',').map(format_value).collect());
    }
    raw.parse::<f64>()
        .map(FilterValue::Number)
        .unwrap_or_else(|_| FilterValue::String(raw.to_owned()))
}

#[derive(Clone, Debug)]
struct KeepTags {
    keep_info: bool,
    keep_fmt: bool,
    info_tags: HashSet<String>,
    fmt_tags: HashSet<String>,
}

impl KeepTags {
    fn parse(raw: Option<&str>) -> Self {
        let Some(raw) = raw else {
            return Self::keep_all();
        };

        let mut keep = Self {
            keep_info: false,
            keep_fmt: false,
            info_tags: HashSet::new(),
            fmt_tags: HashSet::new(),
        };
        let mut mode = TagMode::None;
        for token in raw.split(',').filter(|s| !s.is_empty()) {
            let upper = token.to_ascii_uppercase();
            let (mode_for_tag, tag) = if upper == "INFO" {
                keep.keep_info = true;
                continue;
            } else if upper.starts_with("INFO/") {
                (TagMode::Info, &token[5..])
            } else if upper == "FMT" || upper == "FORMAT" {
                keep.keep_fmt = true;
                continue;
            } else if upper.starts_with("FMT/") {
                (TagMode::Fmt, &token[4..])
            } else if upper.starts_with("FORMAT/") {
                (TagMode::Fmt, &token[7..])
            } else {
                (mode, token)
            };

            mode = mode_for_tag;
            match mode_for_tag {
                TagMode::Info => {
                    keep.info_tags.insert(tag.to_string());
                }
                TagMode::Fmt => {
                    keep.fmt_tags.insert(tag.to_string());
                }
                TagMode::None => {}
            }
        }

        if !keep.keep_info
            && !keep.keep_fmt
            && keep.info_tags.is_empty()
            && keep.fmt_tags.is_empty()
        {
            keep.keep_info = true;
            keep.keep_fmt = true;
        }
        if !keep.keep_fmt && keep.fmt_tags.is_empty() {
            keep.keep_fmt = true;
        }
        keep
    }

    fn keep_all() -> Self {
        Self {
            keep_info: true,
            keep_fmt: true,
            info_tags: HashSet::new(),
            fmt_tags: HashSet::new(),
        }
    }

    fn filter_meta<'a>(&self, lines: &'a [String]) -> Vec<&'a str> {
        lines
            .iter()
            .filter(|line| self.keep_meta_line(line))
            .map(String::as_str)
            .collect()
    }

    fn keep_meta_line(&self, line: &str) -> bool {
        if line.starts_with("##INFO=<") {
            return self.keep_info || header_id(line).is_some_and(|id| self.info_tags.contains(id));
        }
        if line.starts_with("##FORMAT=<") {
            return self.keep_fmt || header_id(line).is_some_and(|id| self.fmt_tags.contains(id));
        }
        true
    }

    fn project_info(&self, raw: &str) -> String {
        if self.keep_info {
            return raw.to_string();
        }
        let kept = raw
            .split(';')
            .filter(|field| {
                let key = field.split_once('=').map(|(key, _)| key).unwrap_or(field);
                self.info_tags.contains(key)
            })
            .collect::<Vec<_>>();
        if kept.is_empty() {
            ".".to_string()
        } else {
            kept.join(";")
        }
    }

    fn project_format_and_samples(&self, format: &str, sample: &str) -> (String, String) {
        if self.keep_fmt {
            return (format.to_string(), sample.to_string());
        }
        let format_fields = format.split(':').collect::<Vec<_>>();
        let sample_fields = sample.split(':').collect::<Vec<_>>();
        let mut out_format = Vec::new();
        let mut out_sample = Vec::new();
        for (idx, key) in format_fields.iter().enumerate() {
            if self.fmt_tags.contains(*key) {
                out_format.push(*key);
                out_sample.push(sample_fields.get(idx).copied().unwrap_or("."));
            }
        }
        if out_format.is_empty() {
            (".".to_string(), ".".to_string())
        } else {
            (out_format.join(":"), out_sample.join(":"))
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum TagMode {
    None,
    Info,
    Fmt,
}

fn default_subsets(samples: &[String]) -> Vec<Subset> {
    let mut unique = UniqueNames::default();
    samples
        .iter()
        .enumerate()
        .map(|(idx, sample)| Subset {
            samples: vec![idx],
            names: vec![sample.clone()],
            fname: unique.create(sample),
        })
        .collect()
}

fn subsets_from_samples_file(path: &Path, samples: &[String]) -> io::Result<Vec<Subset>> {
    let text = fs::read_to_string(path)?;
    let sample_to_idx = sample_index(samples);
    let mut unique = UniqueNames::default();
    let mut out = Vec::new();

    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.is_empty() {
            continue;
        }

        let requested = split_csv(fields[0]);
        let mut idxs = Vec::new();
        let mut kept_samples = Vec::new();
        for sample in requested {
            if let Some(&idx) = sample_to_idx.get(sample) {
                idxs.push(idx);
                kept_samples.push(sample.to_string());
            } else {
                eprintln!("Warning: The sample \"{sample}\" is not present");
            }
        }
        if idxs.is_empty() {
            continue;
        }

        let names = if let Some(raw_renames) = fields.get(1) {
            if *raw_renames == "-" {
                kept_samples
            } else {
                split_csv(raw_renames)
                    .into_iter()
                    .map(ToOwned::to_owned)
                    .collect()
            }
        } else {
            kept_samples
        };
        if names.len() != idxs.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("sample and rename counts differ in {}", path.display()),
            ));
        }

        let fname_seed = fields
            .get(2)
            .copied()
            .filter(|s| *s != "-")
            .unwrap_or(&names[0]);
        let fname = unique.create(fname_seed);
        out.push(Subset {
            samples: idxs,
            names,
            fname,
        });
    }

    Ok(out)
}

fn subsets_from_groups_file(path: &Path, samples: &[String]) -> io::Result<Vec<Subset>> {
    let text = fs::read_to_string(path)?;
    let sample_to_idx = sample_index(samples);
    let mut unique = UniqueNames::default();
    let mut name_to_subset: HashMap<String, usize> = HashMap::new();
    let mut out: Vec<Subset> = Vec::new();

    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        let Some(sample) = fields.first().copied() else {
            continue;
        };
        let Some(&sample_idx) = sample_to_idx.get(sample) else {
            eprintln!("Warning: The sample \"{sample}\" is not present");
            continue;
        };
        let rename = fields
            .get(1)
            .copied()
            .filter(|s| *s != "-")
            .unwrap_or(sample);
        let group_list = fields.get(2).copied().unwrap_or(sample);

        for group in split_csv(group_list) {
            let subset_idx = if let Some(idx) = name_to_subset.get(group).copied() {
                idx
            } else {
                let idx = out.len();
                name_to_subset.insert(group.to_string(), idx);
                out.push(Subset {
                    samples: Vec::new(),
                    names: Vec::new(),
                    fname: unique.create(group),
                });
                idx
            };
            out[subset_idx].samples.push(sample_idx);
            out[subset_idx].names.push(rename.to_string());
        }
    }

    Ok(out)
}

struct ParsedVcf {
    meta: Vec<String>,
    fixed_header: String,
    samples: Vec<String>,
    records: Vec<String>,
}

impl ParsedVcf {
    fn new(text: &str) -> io::Result<Self> {
        let mut meta = Vec::new();
        let mut fixed_header = None;
        let mut samples = Vec::new();
        let mut records = Vec::new();

        for line in text.lines() {
            if line.starts_with("##") {
                meta.push(line.to_string());
            } else if line.starts_with("#CHROM\t") {
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() > 9 {
                    samples = cols[9..].iter().map(|s| (*s).to_string()).collect();
                }
                fixed_header = Some(cols[..cols.len().min(9)].join("\t"));
            } else if !line.trim().is_empty() {
                records.push(line.to_string());
            }
        }

        let fixed_header = fixed_header.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "VCF header line is missing")
        })?;
        if samples.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "No samples to split",
            ));
        }

        Ok(Self {
            meta,
            fixed_header,
            samples,
            records,
        })
    }
}

fn project_record(line: &str, sample_idxs: &[usize], keep: &KeepTags) -> String {
    let cols: Vec<&str> = line.split('\t').collect();
    let mut fixed = cols[..cols.len().min(9)]
        .iter()
        .map(|field| (*field).to_string())
        .collect::<Vec<_>>();
    if fixed.len() > 7 {
        fixed[7] = keep.project_info(&fixed[7]);
    }
    let format_raw = cols.get(8).copied().unwrap_or(".");
    let mut projected_format = None;
    let mut projected_samples = Vec::with_capacity(sample_idxs.len());
    for idx in sample_idxs {
        let sample = cols.get(9 + idx).copied().unwrap_or(".");
        if keep.keep_fmt {
            projected_samples.push(sample.to_string());
        } else {
            let (format, sample) = keep.project_format_and_samples(format_raw, sample);
            projected_format.get_or_insert(format);
            projected_samples.push(sample);
        }
    }
    if fixed.len() > 8 {
        fixed[8] = projected_format.unwrap_or_else(|| format_raw.to_string());
    }
    let mut out = fixed.join("\t");
    for sample in projected_samples {
        out.push('\t');
        out.push_str(&sample);
    }
    out
}

fn header_id(line: &str) -> Option<&str> {
    let body = line
        .strip_prefix("##")?
        .split_once('<')?
        .1
        .strip_suffix('>')?;
    for field in body.split(',') {
        if let Some(id) = field.strip_prefix("ID=") {
            return Some(id);
        }
    }
    None
}

fn sample_index(samples: &[String]) -> HashMap<&str, usize> {
    samples
        .iter()
        .enumerate()
        .map(|(idx, sample)| (sample.as_str(), idx))
        .collect()
}

fn split_csv(raw: &str) -> Vec<&str> {
    raw.split(',').filter(|s| !s.is_empty()).collect()
}

fn sanitize_name(raw: &str) -> String {
    raw.chars()
        .map(|c| match c {
            ':' | '\\' | '/' | ' ' | '\t' => '_',
            _ => c,
        })
        .collect()
}

fn with_vcf_suffix(fname: &str, compress_vcf: bool) -> String {
    let lower = fname.to_ascii_lowercase();
    if lower.ends_with(".bcf")
        || lower.ends_with(".vcf")
        || lower.ends_with(".vcf.gz")
        || lower.ends_with(".vcf.bgz")
    {
        fname.to_string()
    } else if compress_vcf {
        format!("{fname}.vcf.gz")
    } else {
        format!("{fname}.vcf")
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
        ".bcftools-rs-split-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unique_names_are_sanitized_with_suffixes() {
        let mut names = UniqueNames::default();
        assert_eq!(names.create("A/B"), "A_B");
        assert_eq!(names.create("A\\B"), "A_B-1");
        assert_eq!(names.create("A:B"), "A_B-2");
    }

    #[test]
    fn default_split_projects_samples() {
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/0\n";
        let out = compute(vcf, None, None, None, None).unwrap();
        assert_eq!(out[0].0, "A");
        assert!(
            out[0]
                .1
                .contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n")
        );
        assert!(out[0].1.contains("1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\n"));
    }

    #[test]
    fn keep_tags_projects_info_and_format() {
        let vcf = "\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Allele\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA
1\t1\t.\tA\tC\t.\t.\tDP=5;AA=T\tGT:AD:GQ\t0/1:2,3:9
";
        let out = compute(vcf, None, None, Some("INFO/DP,FMT/GT,AD"), None).unwrap();
        assert!(out[0].1.contains("##INFO=<ID=DP,"));
        assert!(!out[0].1.contains("##INFO=<ID=AA,"));
        assert!(out[0].1.contains("##FORMAT=<ID=GT,"));
        assert!(out[0].1.contains("##FORMAT=<ID=AD,"));
        assert!(!out[0].1.contains("##FORMAT=<ID=GQ,"));
        assert!(
            out[0]
                .1
                .contains("1\t1\t.\tA\tC\t.\t.\tDP=5\tGT:AD\t0/1:2,3\n")
        );
    }

    #[test]
    fn include_filter_is_applied_after_sample_projection() {
        let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/0
1\t2\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1
";
        let out = compute(
            vcf,
            None,
            None,
            None,
            Some(FilterSpec {
                mode: FilterMode::Include,
                expr: r#"GT="alt""#,
            }),
        )
        .unwrap();
        assert!(out[0].1.contains("1\t1\t.\tA\tC"));
        assert!(!out[0].1.contains("1\t2\t.\tA\tG"));
        assert!(!out[1].1.contains("1\t1\t.\tA\tC"));
        assert!(out[1].1.contains("1\t2\t.\tA\tG"));
    }
}
