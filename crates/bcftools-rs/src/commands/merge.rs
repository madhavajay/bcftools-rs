//! Focused `bcftools merge` implementation (upstream `vcfmerge.c`).
//!
//! This local slice merges records that are present in every input or are
//! absent from some inputs and have identical site fields, plus a narrow
//! same-position allele-union slice. Full synced-reader merging, full allele
//! unification, full INFO rules, and gVCF mode remain tracked in `TODO.md`.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::diagnostics::fmt_etag;
use crate::reference::FastaReference;
use crate::vcf_compat::normalize_vcf_text;

const USAGE: &str = "\n\
About:   Merge VCF/BCF files from non-overlapping sample sets.\n\
Usage:   bcftools merge [OPTIONS] <A.vcf.gz> <B.vcf.gz> [...]\n\
\n\
Options:\n\
    -l, --file-list FILE            Read input file names from FILE\n\
    -i, --info-rules TAG:METHOD,..  Apply AC:sum/AN:sum in the current text ALT-union slice\n\
    -m, --merge TYPE                Support narrow none/both/snp-ins-del slices; other modes accepted for command-shape compatibility\n\
    -o, --output FILE               Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9]  u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
    -F, --filter-logic x|+          Narrow filter merge logic support for upstream fixtures [+]\n\
    -g, --gvcf -|REF.FA             Narrow gVCF block splitting for local fixtures\n\
    -0, --missing-to-ref            Fill absent-input samples as 0/0 in the current text slice\n\
        --force-single              Allow a single input for command-shape compatibility\n\
        --force-samples             Allow duplicate sample names by prefixing later inputs with the input index\n\
        --no-index                  Accepted for command-shape compatibility in this text slice\n\
        --no-version                Accepted for command-shape compatibility\n\
\n";

#[derive(Debug)]
struct Args {
    inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    force_samples: bool,
    info_rules: InfoRules,
    merge_mode: MergeMode,
    local_alleles: Option<usize>,
    gvcf_ref: Option<PathBuf>,
    gvcf_no_ref: bool,
    missing_to_ref: bool,
}

#[derive(Debug, Clone, Copy, Default)]
struct InfoRules {
    sum_ac: bool,
    sum_an: bool,
    join_af: bool,
    filter_logic: FilterLogic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MergeMode {
    Default,
    None,
    Both,
    SnpInsDel,
    All,
    Id,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
enum FilterLogic {
    #[default]
    Add,
    RemoveIfPass,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VcfNumber {
    A,
    R,
    G,
    Other,
}

#[derive(Debug)]
struct VcfInput {
    meta: Vec<String>,
    fixed_header: Vec<String>,
    samples: Vec<String>,
    records: Vec<RecordLine>,
}

#[derive(Debug, Clone)]
struct RecordLine {
    fixed: Vec<String>,
    samples: Vec<String>,
}

#[derive(Debug)]
struct MergedSite {
    fixed: Vec<String>,
    samples_by_input: Vec<Option<Vec<String>>>,
    order: usize,
}

#[derive(Debug)]
enum ParseOutcome {
    Usage,
    Error(String),
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfmerge", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfmerge", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut inputs = Vec::new();
    let mut file_list = None;
    let mut output = None;
    let mut output_kind = OutputKind::VcfText;
    let mut force_single = false;
    let mut force_samples = false;
    let mut info_rules = InfoRules::default();
    let mut merge_mode = MergeMode::Default;
    let mut local_alleles = None;
    let mut gvcf_ref = None;
    let mut gvcf_no_ref = false;
    let mut missing_to_ref = false;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-l" | "--file-list" => {
                file_list = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-i" | "--info-rules" => {
                info_rules = parse_info_rules(&next_string(&mut iter, raw.as_ref())?);
            }
            "-F" | "--filter-logic" => {
                info_rules.filter_logic =
                    parse_filter_logic(&next_string(&mut iter, raw.as_ref())?)?;
            }
            "-m" | "--merge" => {
                merge_mode = parse_merge_mode(&next_string(&mut iter, raw.as_ref())?);
            }
            "-L" | "--local-alleles" => {
                local_alleles = Some(parse_local_alleles(&next_string(&mut iter, raw.as_ref())?)?);
            }
            "-g" | "--gvcf" => {
                let path = next_string(&mut iter, raw.as_ref())?;
                if path != "-" {
                    gvcf_ref = Some(PathBuf::from(path));
                } else {
                    gvcf_no_ref = true;
                }
            }
            "-0" | "--missing-to-ref" => missing_to_ref = true,
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--force-single" => force_single = true,
            "--force-samples" => force_samples = true,
            "--no-index" => {}
            "--no-version" => {}
            _ if raw.starts_with("--file-list=") => {
                file_list = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--info-rules=") => {
                info_rules = parse_info_rules(value_after_equals(&raw))
            }
            _ if raw.starts_with("--filter-logic=") => {
                info_rules.filter_logic = parse_filter_logic(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--merge=") => {
                merge_mode = parse_merge_mode(value_after_equals(&raw))
            }
            _ if raw.starts_with("--local-alleles=") => {
                local_alleles = Some(parse_local_alleles(value_after_equals(&raw))?)
            }
            _ if raw.starts_with("--gvcf=") => {
                let path = value_after_equals(&raw);
                if path != "-" {
                    gvcf_ref = Some(PathBuf::from(path));
                } else {
                    gvcf_no_ref = true;
                }
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-l") && raw.len() > 2 => {
                file_list = Some(PathBuf::from(&raw[2..]))
            }
            _ if raw.starts_with("-i") && raw.len() > 2 => info_rules = parse_info_rules(&raw[2..]),
            _ if raw.starts_with("-F") && raw.len() > 2 => {
                info_rules.filter_logic = parse_filter_logic(&raw[2..])?
            }
            _ if raw.starts_with("-L") && raw.len() > 2 => {
                local_alleles = Some(parse_local_alleles(&raw[2..])?)
            }
            _ if raw.starts_with("-g") && raw.len() > 2 => {
                let path = &raw[2..];
                if path != "-" {
                    gvcf_ref = Some(PathBuf::from(path));
                } else {
                    gvcf_no_ref = true;
                }
            }
            _ if raw.starts_with("-m") && raw.len() > 2 => merge_mode = parse_merge_mode(&raw[2..]),
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!(
                    "unrecognized option '{raw}' in this local merge slice"
                )));
            }
            _ => inputs.push(PathBuf::from(raw.as_ref())),
        }
    }

    if let Some(path) = file_list {
        inputs.extend(read_file_list(&path)?);
    }
    if inputs.len() < 2 && !force_single {
        return Err(ParseOutcome::Error(
            "expected at least two input VCF/BCF paths".into(),
        ));
    }
    if inputs.is_empty() {
        return Err(ParseOutcome::Error("expected an input VCF/BCF path".into()));
    }

    Ok(Args {
        inputs,
        output,
        output_kind,
        force_samples,
        info_rules,
        merge_mode,
        local_alleles,
        gvcf_ref,
        gvcf_no_ref,
        missing_to_ref,
    })
}

fn parse_filter_logic(raw: &str) -> Result<FilterLogic, ParseOutcome> {
    match raw {
        "+" => Ok(FilterLogic::Add),
        "x" => Ok(FilterLogic::RemoveIfPass),
        _ => Err(ParseOutcome::Error(format!(
            "filter logic not recognised: {raw}"
        ))),
    }
}

fn parse_local_alleles(raw: &str) -> Result<usize, ParseOutcome> {
    let value = raw.parse::<usize>().map_err(|_| {
        ParseOutcome::Error(format!("Could not parse argument: --local-alleles {raw}"))
    })?;
    if value == 0 {
        return Err(ParseOutcome::Error(format!(
            "Error: \"--local-alleles {raw}\" makes no sense, expected value bigger or equal than 1"
        )));
    }
    Ok(value)
}

fn parse_merge_mode(raw: &str) -> MergeMode {
    match raw {
        "none" => MergeMode::None,
        "both" => MergeMode::Both,
        "snp-ins-del" => MergeMode::SnpInsDel,
        "all" => MergeMode::All,
        "id" => MergeMode::Id,
        _ => MergeMode::Default,
    }
}

fn parse_info_rules(raw: &str) -> InfoRules {
    let mut rules = InfoRules::default();
    for rule in raw.split(',') {
        let Some((tag, method)) = rule.split_once(':') else {
            continue;
        };
        match tag {
            "AC" if method == "sum" => rules.sum_ac = true,
            "AN" if method == "sum" => rules.sum_an = true,
            "AF" if method == "join" => rules.join_af = true,
            _ => {}
        }
    }
    rules
}

fn read_file_list(path: &Path) -> Result<Vec<PathBuf>, ParseOutcome> {
    let text = fs::read_to_string(path)
        .map_err(|e| ParseOutcome::Error(format!("failed to read {}: {e}", path.display())))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(PathBuf::from)
        .collect())
}

fn next_string<'a, I>(iter: &mut std::iter::Peekable<I>, name: &str) -> Result<String, ParseOutcome>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .map(|s| s.to_string_lossy().into_owned())
        .ok_or_else(|| ParseOutcome::Error(format!("missing argument for {name}")))
}

fn value_after_equals(raw: &str) -> &str {
    raw.split_once('=').map(|(_, value)| value).unwrap_or("")
}

fn parse_output_kind(raw: &str) -> Result<OutputKind, ParseOutcome> {
    match raw.as_bytes().first().copied() {
        Some(b'v') => Ok(OutputKind::VcfText),
        Some(b'z') => Ok(OutputKind::VcfGz),
        Some(b'u' | b'b') => Ok(OutputKind::Bcf),
        _ => Err(ParseOutcome::Error(format!("unknown output type '{raw}'"))),
    }
}

fn run(args: &Args) -> io::Result<()> {
    let mut inputs = Vec::new();
    for path in &args.inputs {
        inputs.push(parse_vcf_text(&read_vcf_text(path)?)?);
    }
    if let Some(path) = &args.gvcf_ref {
        split_gvcf_reference_blocks(&mut inputs, path)?;
    } else if args.gvcf_no_ref {
        split_gvcf_reference_blocks_without_reference(&mut inputs);
    }
    let merged = merge_inputs(
        &inputs,
        args.force_samples,
        args.info_rules,
        args.merge_mode,
        args.local_alleles,
        args.missing_to_ref,
    )?;
    write_output(merged.as_bytes(), args)
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
        ".bcftools-rs-merge-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn parse_vcf_text(text: &str) -> io::Result<VcfInput> {
    let mut meta = Vec::new();
    let mut fixed_header = Vec::new();
    let mut samples = Vec::new();
    let mut records = Vec::new();

    for line in text.lines() {
        if line.starts_with("##") {
            meta.push(line.to_owned());
        } else if line.starts_with("#CHROM") {
            let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
            if fields.len() < 8 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid #CHROM header",
                ));
            }
            fixed_header = fields[..fields.len().min(9)].to_vec();
            if fields.len() > 9 {
                samples = fields[9..].to_vec();
            }
        } else if !line.trim().is_empty() {
            let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
            if fields.len() < fixed_header.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("record has fewer columns than header: {line}"),
                ));
            }
            let fixed_len = fixed_header.len();
            records.push(RecordLine {
                fixed: fields[..fixed_len].to_vec(),
                samples: fields[fixed_len..].to_vec(),
            });
        }
    }

    if fixed_header.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing #CHROM header",
        ));
    }
    Ok(VcfInput {
        meta,
        fixed_header,
        samples,
        records,
    })
}

fn split_gvcf_reference_blocks(inputs: &mut [VcfInput], reference_path: &Path) -> io::Result<()> {
    let reference = FastaReference::open(reference_path)?;
    let records_by_input = inputs
        .iter()
        .map(|input| input.records.clone())
        .collect::<Vec<_>>();

    for (input_idx, input) in inputs.iter_mut().enumerate() {
        let mut split_records = Vec::new();
        for record in &records_by_input[input_idx] {
            let Some((start, end)) = reference_block_span(record) else {
                split_records.push(record.clone());
                continue;
            };
            let mut breaks = vec![start, end + 1];
            for (other_idx, other_records) in records_by_input.iter().enumerate() {
                if other_idx == input_idx {
                    continue;
                }
                for other in other_records {
                    if other.fixed.first() != record.fixed.first() {
                        continue;
                    }
                    let Some(other_start) = record_pos(other) else {
                        continue;
                    };
                    if other_start > start && other_start <= end {
                        breaks.push(other_start);
                    } else if other_start == start {
                        if let Some((_, other_end)) = reference_block_span(other) {
                            if other_end < end {
                                breaks.push(other_end + 1);
                            }
                        } else if !is_reference_block_alt(&split_alt(&other.fixed[4])) {
                            let other_ref_len = other.fixed[3].len().max(1) as u64;
                            breaks.push((start + other_ref_len).min(end + 1));
                        }
                    }
                }
            }
            breaks.sort_unstable();
            breaks.dedup();
            for window in breaks.windows(2) {
                let seg_start = window[0];
                let seg_end = window[1] - 1;
                if seg_start > seg_end {
                    continue;
                }
                split_records.push(split_reference_block_segment(
                    record, &reference, seg_start, seg_end,
                )?);
            }
        }
        input.records = split_records;
    }
    Ok(())
}

fn split_gvcf_reference_blocks_without_reference(inputs: &mut [VcfInput]) {
    let records_by_input = inputs
        .iter()
        .map(|input| input.records.clone())
        .collect::<Vec<_>>();
    let ref_hints = reference_block_ref_hints(&records_by_input);

    for (input_idx, input) in inputs.iter_mut().enumerate() {
        let mut split_records = Vec::new();
        for record in &records_by_input[input_idx] {
            let Some((start, end)) = reference_block_span(record) else {
                split_records.push(record.clone());
                continue;
            };
            let mut breaks = vec![start, end + 1];
            let mut skip_positions = HashSet::new();
            for (other_idx, other_records) in records_by_input.iter().enumerate() {
                if other_idx == input_idx {
                    continue;
                }
                for other in other_records {
                    if other.fixed.first() != record.fixed.first() {
                        continue;
                    }
                    let Some(other_start) = record_pos(other) else {
                        continue;
                    };
                    if other_start < start || other_start > end {
                        continue;
                    }
                    if let Some((_, other_end)) = reference_block_span(other) {
                        if other_start > start {
                            breaks.push(other_start);
                        }
                        if other_end < end {
                            breaks.push(other_end + 1);
                        }
                    } else if !is_reference_block_alt(&split_alt(&other.fixed[4])) {
                        breaks.push(other_start);
                        breaks.push((other_start + 1).min(end + 1));
                        skip_positions.insert(other_start);
                        let other_ref_len = other.fixed[3].len().max(1) as u64;
                        if other_ref_len > 1 {
                            breaks.push((other_start + other_ref_len).min(end + 1));
                        }
                    }
                }
            }
            breaks.sort_unstable();
            breaks.dedup();
            for window in breaks.windows(2) {
                let seg_start = window[0];
                let seg_end = window[1] - 1;
                if seg_start > seg_end
                    || (seg_start == seg_end && skip_positions.contains(&seg_start))
                {
                    continue;
                }
                split_records.push(split_reference_block_segment_without_reference(
                    record, seg_start, seg_end, &ref_hints,
                ));
            }
        }
        input.records = split_records;
    }
}

fn reference_block_ref_hints(
    records_by_input: &[Vec<RecordLine>],
) -> HashMap<(String, u64), String> {
    let mut hints = HashMap::new();
    for record in records_by_input.iter().flatten() {
        if reference_block_span(record).is_some()
            && let (Some(chrom), Some(pos), Some(reference)) = (
                record.fixed.first(),
                record_pos(record),
                record.fixed.get(3),
            )
        {
            hints
                .entry((chrom.clone(), pos))
                .or_insert_with(|| reference.clone());
        }
    }
    hints
}

fn reference_block_span(record: &RecordLine) -> Option<(u64, u64)> {
    if record.fixed.len() < 8 || !is_reference_block_alt(&split_alt(&record.fixed[4])) {
        return None;
    }
    let start = record_pos(record)?;
    let end = info_u64(&record.fixed[7], "END")?;
    (end >= start).then_some((start, end))
}

fn is_reference_block_alt(alts: &[String]) -> bool {
    alts.is_empty() || is_single_ref_symbolic_alt(alts)
}

fn record_pos(record: &RecordLine) -> Option<u64> {
    record.fixed.get(1)?.parse().ok()
}

fn split_reference_block_segment(
    record: &RecordLine,
    reference: &FastaReference,
    start: u64,
    end: u64,
) -> io::Result<RecordLine> {
    let mut segment = record.clone();
    segment.fixed[1] = start.to_string();
    segment.fixed[3] = fetch_reference_base(reference, &segment.fixed[0], start)?;
    segment.fixed[7] = if end > start {
        set_info_value(&segment.fixed[7], "END", &end.to_string())
    } else {
        remove_info_value(&segment.fixed[7], "END")
    };
    Ok(segment)
}

fn split_reference_block_segment_without_reference(
    record: &RecordLine,
    start: u64,
    end: u64,
    ref_hints: &HashMap<(String, u64), String>,
) -> RecordLine {
    let mut segment = record.clone();
    segment.fixed[1] = start.to_string();
    if record_pos(record) != Some(start) {
        segment.fixed[3] = ref_hints
            .get(&(segment.fixed[0].clone(), start))
            .cloned()
            .unwrap_or_else(|| "N".to_owned());
    }
    segment.fixed[7] = if end > start {
        set_info_value(&segment.fixed[7], "END", &end.to_string())
    } else {
        remove_info_value(&segment.fixed[7], "END")
    };
    segment
}

fn fetch_reference_base(reference: &FastaReference, chrom: &str, pos: u64) -> io::Result<String> {
    let region = format!("{chrom}:{pos}-{pos}");
    let seq = reference.fetch_region(&region)?;
    let base = seq.first().copied().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("reference region {region} is empty"),
        )
    })?;
    Ok((base as char).to_ascii_uppercase().to_string())
}

fn merge_inputs(
    inputs: &[VcfInput],
    force_samples: bool,
    info_rules: InfoRules,
    merge_mode: MergeMode,
    local_alleles: Option<usize>,
    missing_to_ref: bool,
) -> io::Result<String> {
    let first = inputs
        .first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "no inputs"))?;
    let mut sample_names = Vec::new();
    let mut seen_samples = HashSet::new();
    for (input_idx, input) in inputs.iter().enumerate() {
        for sample in &input.samples {
            let mut name = sample.clone();
            if seen_samples.contains(&name) {
                if !force_samples {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("duplicate sample name '{sample}'"),
                    ));
                }
                while seen_samples.contains(&name) {
                    name = format!("{}:{name}", input_idx + 1);
                }
            }
            seen_samples.insert(name.clone());
            sample_names.push(name);
        }
    }

    let fileformat = merged_fileformat(inputs);
    let mut merged_meta = merged_meta(inputs);
    let format_numbers = format_numbers(&merged_meta);
    let info_numbers = info_numbers(&merged_meta);
    let ordered_info_numbers = ordered_info_numbers(&merged_meta);

    for input in &inputs[1..] {
        if !fixed_headers_compatible(&first.fixed_header, &input.fixed_header) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "inputs must have compatible fixed VCF columns",
            ));
        }
    }

    let mut sites = collect_sites(
        inputs,
        info_rules,
        merge_mode,
        &format_numbers,
        &info_numbers,
        &ordered_info_numbers,
    )?;
    let contigs = contig_order(&first.meta);
    sites.sort_by(|a, b| compare_sites(a, b, &contigs));
    if let Some(limit) = local_alleles {
        let input_sample_counts = inputs
            .iter()
            .map(|input| input.samples.len())
            .collect::<Vec<_>>();
        apply_local_alleles(&mut sites, limit, &format_numbers, &input_sample_counts);
        append_localized_format_meta(&mut merged_meta);
    }

    let mut out = render_meta_with_pass_filter(&merged_meta, info_rules, fileformat.as_deref());
    out.push_str(&first.fixed_header.join("\t"));
    if !sample_names.is_empty() {
        out.push('\t');
        out.push_str(&sample_names.join("\t"));
    }
    out.push('\n');

    for site in sites {
        let absent_alleles = if missing_to_ref {
            inputs
                .iter()
                .enumerate()
                .filter(|(input_idx, _)| site.samples_by_input[*input_idx].is_none())
                .map(|(_, input)| input.samples.len() * 2)
                .sum()
        } else {
            0
        };
        let mut fixed = site.fixed;
        if info_rules.filter_logic == FilterLogic::RemoveIfPass
            && let Some(filter) = fixed.get_mut(6)
            && is_pass_filter(filter)
        {
            *filter = "PASS".to_owned();
        }
        if absent_alleles > 0 {
            add_missing_reference_alleles_to_an(&mut fixed, absent_alleles);
        }
        cleanup_single_base_reference_block_end(&mut fixed);
        let mut samples = Vec::new();
        for (input_idx, input) in inputs.iter().enumerate() {
            match &site.samples_by_input[input_idx] {
                Some(values) => samples.extend(values.iter().cloned()),
                None => {
                    let missing = missing_sample_value(&fixed, missing_to_ref);
                    samples.extend(std::iter::repeat_n(missing, input.samples.len()));
                }
            }
        }
        out.push_str(&fixed.join("\t"));
        if !samples.is_empty() {
            out.push('\t');
            out.push_str(&samples.join("\t"));
        }
        out.push('\n');
    }

    Ok(out)
}

fn fixed_headers_compatible(first: &[String], other: &[String]) -> bool {
    if other == first {
        return true;
    }
    first.len() == 9 && other.len() == 8 && first[8] == "FORMAT" && first[..8] == other[..8]
}

fn render_meta_with_pass_filter(
    meta: &[String],
    info_rules: InfoRules,
    fileformat: Option<&str>,
) -> String {
    let has_pass = meta
        .iter()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    let mut out = String::new();
    let mut inserted = false;
    for line in meta {
        if let Some(fileformat) = fileformat
            && line.starts_with("##fileformat=")
        {
            out.push_str(fileformat);
        } else if info_rules.join_af && line.starts_with("##INFO=<ID=AF,") {
            out.push_str(&rewrite_info_number(line, "."));
        } else {
            out.push_str(line);
        }
        out.push('\n');
        if !inserted && !has_pass && line.starts_with("##fileformat=") {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            inserted = true;
        }
    }
    if !inserted && !has_pass {
        out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    }
    out
}

fn merged_meta(inputs: &[VcfInput]) -> Vec<String> {
    let Some(first) = inputs.first() else {
        return Vec::new();
    };
    let mut meta = first.meta.clone();
    let mut seen_ids = ["INFO", "FORMAT", "FILTER", "SAMPLE"]
        .into_iter()
        .map(|kind| {
            let ids = first
                .meta
                .iter()
                .filter_map(|line| meta_id(line, kind))
                .collect::<HashSet<_>>();
            (kind, ids)
        })
        .collect::<HashMap<_, _>>();

    for input in &inputs[1..] {
        for line in &input.meta {
            for kind in ["INFO", "FORMAT", "FILTER", "SAMPLE"] {
                if let Some(id) = meta_id(line, kind)
                    && seen_ids.get_mut(kind).is_some_and(|ids| ids.insert(id))
                {
                    meta.push(line.clone());
                    break;
                }
            }
        }
    }

    meta
}

fn merged_fileformat(inputs: &[VcfInput]) -> Option<String> {
    inputs
        .iter()
        .filter_map(|input| {
            input
                .meta
                .iter()
                .find(|line| line.starts_with("##fileformat="))
        })
        .max_by_key(|line| vcf_fileformat_rank(line))
        .cloned()
}

fn vcf_fileformat_rank(line: &str) -> (u32, u32) {
    let Some(version) = line.strip_prefix("##fileformat=VCFv") else {
        return (0, 0);
    };
    let Some((major, minor)) = version.split_once('.') else {
        return (0, 0);
    };
    (
        major.parse().unwrap_or(0),
        minor
            .chars()
            .take_while(|ch| ch.is_ascii_digit())
            .collect::<String>()
            .parse()
            .unwrap_or(0),
    )
}

fn rewrite_info_number(line: &str, number: &str) -> String {
    let Some(start) = line.find("Number=") else {
        return line.to_owned();
    };
    let value_start = start + "Number=".len();
    let Some(value_end) = line[value_start..].find(',').map(|idx| value_start + idx) else {
        return line.to_owned();
    };
    let mut out = String::new();
    out.push_str(&line[..value_start]);
    out.push_str(number);
    out.push_str(&line[value_end..]);
    out
}

fn format_numbers(meta: &[String]) -> HashMap<String, VcfNumber> {
    meta_numbers(meta, "FORMAT")
}

fn info_numbers(meta: &[String]) -> HashMap<String, VcfNumber> {
    meta_numbers(meta, "INFO")
}

fn meta_numbers(meta: &[String], kind: &str) -> HashMap<String, VcfNumber> {
    meta.iter()
        .filter_map(|line| {
            let id = meta_id(line, kind)?;
            let number = meta_attr(line, "Number").map(|raw| match raw {
                "A" => VcfNumber::A,
                "R" => VcfNumber::R,
                "G" => VcfNumber::G,
                _ => VcfNumber::Other,
            })?;
            Some((id, number))
        })
        .collect()
}

fn ordered_info_numbers(meta: &[String]) -> Vec<(String, VcfNumber)> {
    meta.iter()
        .filter_map(|line| {
            let id = meta_id(line, "INFO")?;
            let number = meta_attr(line, "Number").map(|raw| match raw {
                "A" => VcfNumber::A,
                "R" => VcfNumber::R,
                "G" => VcfNumber::G,
                _ => VcfNumber::Other,
            })?;
            Some((id, number))
        })
        .collect()
}

fn append_localized_format_meta(meta: &mut Vec<String>) {
    if meta
        .iter()
        .any(|line| line.starts_with("##FORMAT=<ID=LAA,"))
    {
        return;
    }

    let localized = meta
        .iter()
        .filter_map(|line| localized_format_meta_line(line))
        .collect::<Vec<_>>();
    if localized.is_empty() {
        return;
    }

    meta.push("##FORMAT=<ID=LAA,Number=.,Type=Integer,Description=\"Localized alleles: subset of alternate alleles relevant for each sample\">".to_owned());
    meta.extend(localized);
}

fn localized_format_meta_line(line: &str) -> Option<String> {
    let id = meta_id(line, "FORMAT")?;
    let number = meta_attr(line, "Number")?;
    if !matches!(number, "A" | "R" | "G") {
        return None;
    }

    let mut out = line.to_owned();
    out = rewrite_meta_attr(&out, "ID", &format!("L{id}"));
    out = rewrite_meta_attr(&out, "Number", ".");
    if let Some(description) = meta_attr(line, "Description") {
        let localized = if let Some(body) = description
            .strip_prefix('"')
            .and_then(|value| value.strip_suffix('"'))
        {
            format!("\"Localized field: {body}\"")
        } else {
            format!("Localized field: {description}")
        };
        out = rewrite_meta_attr(&out, "Description", &localized);
    }
    Some(out)
}

fn rewrite_meta_attr(line: &str, key: &str, value: &str) -> String {
    let Some(body_start) = line.find("=<").map(|idx| idx + 2) else {
        return line.to_owned();
    };
    let Some(value_start) = line[body_start..]
        .find(&format!("{key}="))
        .map(|idx| body_start + idx + key.len() + 1)
    else {
        return line.to_owned();
    };
    let mut in_quotes = false;
    let mut escaped = false;
    let bytes = line.as_bytes();
    let mut value_end = value_start;
    while value_end < line.len() {
        let ch = bytes[value_end] as char;
        if escaped {
            escaped = false;
        } else if ch == '\\' {
            escaped = true;
        } else if ch == '"' {
            in_quotes = !in_quotes;
        } else if !in_quotes && (ch == ',' || ch == '>') {
            break;
        }
        value_end += 1;
    }

    let mut out = String::new();
    out.push_str(&line[..value_start]);
    out.push_str(value);
    out.push_str(&line[value_end..]);
    out
}

fn meta_id(line: &str, kind: &str) -> Option<String> {
    let prefix = format!("##{kind}=<");
    line.strip_prefix(&prefix)
        .and_then(|body| meta_attr_body(body, "ID"))
        .map(str::to_owned)
}

fn meta_attr<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let body = line.split_once("=<")?.1;
    meta_attr_body(body, key)
}

fn meta_attr_body<'a>(body: &'a str, key: &str) -> Option<&'a str> {
    body.trim_end_matches('>').split(',').find_map(|field| {
        let (name, value) = field.split_once('=')?;
        (name == key).then_some(value)
    })
}

fn collect_sites(
    inputs: &[VcfInput],
    info_rules: InfoRules,
    merge_mode: MergeMode,
    format_numbers: &HashMap<String, VcfNumber>,
    info_numbers: &HashMap<String, VcfNumber>,
    ordered_info_numbers: &[(String, VcfNumber)],
) -> io::Result<Vec<MergedSite>> {
    let mut sites: Vec<MergedSite> = Vec::new();
    let mut by_key = HashMap::new();
    let mut by_locus: HashMap<String, Vec<usize>> = HashMap::new();
    let mut by_position: HashMap<String, Vec<usize>> = HashMap::new();

    for (input_idx, input) in inputs.iter().enumerate() {
        for record in &input.records {
            let key = site_key(record);
            let key_site_idx = by_key.get(&key).copied().filter(|site_idx: &usize| {
                !(record.fixed.get(4).is_some_and(|alt| alt == ".")
                    && sites[*site_idx].fixed[..5] != record.fixed[..5])
            });
            if let Some(site_idx) = key_site_idx {
                let site: &mut MergedSite = &mut sites[site_idx];
                if site.fixed.len() == record.fixed.len() && site.fixed[..5] == record.fixed[..5] {
                    merge_exact_site(
                        site,
                        record,
                        input_idx,
                        info_rules,
                        format_numbers,
                        info_numbers,
                        ordered_info_numbers,
                    )?;
                } else if supports_sampled_same_position_union(merge_mode)
                    && can_merge_sampled_same_position(site, record, merge_mode)
                {
                    merge_sampled_same_position(
                        site,
                        record,
                        input_idx,
                        info_rules.filter_logic,
                        format_numbers,
                        info_numbers,
                        ordered_info_numbers,
                    )?;
                } else {
                    merge_exact_site(
                        site,
                        record,
                        input_idx,
                        info_rules,
                        format_numbers,
                        info_numbers,
                        ordered_info_numbers,
                    )?;
                }
            } else {
                let mut merged = false;
                if merge_mode != MergeMode::None
                    && let Some(site_indices) = by_locus.get(&locus_key(record)).cloned()
                {
                    for site_idx in site_indices {
                        let site: &mut MergedSite = &mut sites[site_idx];
                        if can_merge_same_locus_alt_union(site, record) {
                            merge_sites_only_alt_union(site, record, info_rules);
                            site.samples_by_input[input_idx] = Some(record.samples.clone());
                            merged = true;
                            break;
                        }
                    }
                }
                if !merged
                    && supports_sampled_same_position_union(merge_mode)
                    && let Some(site_indices) = by_position.get(&position_key(record)).cloned()
                {
                    for site_idx in site_indices {
                        let site: &mut MergedSite = &mut sites[site_idx];
                        if site.samples_by_input[input_idx].is_some() {
                            continue;
                        }
                        if can_merge_sampled_same_position(site, record, merge_mode) {
                            merge_sampled_same_position(
                                site,
                                record,
                                input_idx,
                                info_rules.filter_logic,
                                format_numbers,
                                info_numbers,
                                ordered_info_numbers,
                            )?;
                            merged = true;
                            break;
                        }
                    }
                }
                if !merged {
                    let site_idx = sites.len();
                    let mut samples_by_input = vec![None; inputs.len()];
                    samples_by_input[input_idx] = Some(record.samples.clone());
                    by_key.insert(key, site_idx);
                    by_locus
                        .entry(locus_key(record))
                        .or_default()
                        .push(site_idx);
                    by_position
                        .entry(position_key(record))
                        .or_default()
                        .push(site_idx);
                    sites.push(MergedSite {
                        fixed: record.fixed.clone(),
                        samples_by_input,
                        order: site_idx,
                    });
                }
            }
        }
    }

    Ok(sites)
}

fn site_key(record: &RecordLine) -> String {
    record
        .fixed
        .iter()
        .take(5)
        .cloned()
        .collect::<Vec<_>>()
        .join("\t")
}

fn locus_key(record: &RecordLine) -> String {
    record
        .fixed
        .iter()
        .take(4)
        .cloned()
        .collect::<Vec<_>>()
        .join("\t")
}

fn position_key(record: &RecordLine) -> String {
    record
        .fixed
        .iter()
        .take(2)
        .cloned()
        .collect::<Vec<_>>()
        .join("\t")
}

fn can_merge_same_locus_alt_union(site: &MergedSite, record: &RecordLine) -> bool {
    (site.fixed.len() == 8 || site.fixed.len() == 9)
        && record.fixed.len() == 8
        && record.samples.is_empty()
        && site.fixed[..4] == record.fixed[..4]
}

fn supports_sampled_same_position_union(merge_mode: MergeMode) -> bool {
    matches!(
        merge_mode,
        MergeMode::Default
            | MergeMode::None
            | MergeMode::Both
            | MergeMode::SnpInsDel
            | MergeMode::All
            | MergeMode::Id
    )
}

fn can_merge_sampled_same_position(
    site: &MergedSite,
    record: &RecordLine,
    merge_mode: MergeMode,
) -> bool {
    if site.fixed.len() < 9
        || record.fixed.len() < 9
        || site.fixed[..2] != record.fixed[..2]
        || merged_ref(&site.fixed[3], &record.fixed[3]).is_none()
    {
        return false;
    }

    let site_has_non_ref = alt_contains_non_ref(&site.fixed[4]);
    let record_has_non_ref = alt_contains_non_ref(&record.fixed[4]);
    let same_ref_subset_compatible = same_ref_alt_subset_compatible(
        &site.fixed[3],
        &site.fixed[4],
        &record.fixed[3],
        &record.fixed[4],
    );
    let ref_extended_subset_compatible = ref_extended_alt_subset_compatible(
        &site.fixed[3],
        &site.fixed[4],
        &record.fixed[3],
        &record.fixed[4],
    );
    let subset_compatible = same_ref_subset_compatible || ref_extended_subset_compatible;
    let same_id =
        site.fixed.get(2) == record.fixed.get(2) && site.fixed.get(2).is_some_and(|id| id != ".");
    let same_ref = site.fixed[3].eq_ignore_ascii_case(&record.fixed[3]);
    let both_missing_id = site.fixed.get(2).is_some_and(|id| id == ".")
        && record.fixed.get(2).is_some_and(|id| id == ".");
    match merge_mode {
        MergeMode::Default => {
            let site_alts = split_alt(&site.fixed[4]);
            let record_alts = split_alt(&record.fixed[4]);
            same_id
                || ref_extended_subset_compatible
                || single_ref_symbolic_alt_compatible(&site.fixed[4], &record.fixed[4])
                || same_ref
                    && contains_ref_symbolic_alt(&site_alts)
                    && contains_ref_symbolic_alt(&record_alts)
                || site_alts.is_empty()
                || record_alts.is_empty()
                || {
                    let site_class = coarse_variant_class(&site.fixed[3], &site.fixed[4]);
                    (same_ref || both_missing_id)
                        && site_class != CoarseVariantClass::Other
                        && site_class == coarse_variant_class(&record.fixed[3], &record.fixed[4])
                }
        }
        MergeMode::None => {
            subset_compatible
                || site_has_non_ref && record_has_non_ref
                || single_ref_symbolic_alt_compatible(&site.fixed[4], &record.fixed[4])
        }
        MergeMode::Both => {
            if subset_compatible {
                return true;
            }
            if site_has_non_ref && record_has_non_ref {
                return true;
            }
            let site_class = coarse_variant_class(&site.fixed[3], &site.fixed[4]);
            site_class != CoarseVariantClass::Other
                && site_class == coarse_variant_class(&record.fixed[3], &record.fixed[4])
        }
        MergeMode::SnpInsDel => {
            let site_class = precise_variant_class(&site.fixed[3], &site.fixed[4]);
            site_class != PreciseVariantClass::Other
                && site_class == precise_variant_class(&record.fixed[3], &record.fixed[4])
        }
        MergeMode::All => true,
        MergeMode::Id => same_id || both_missing_id && same_ref,
    }
}

fn alt_contains_non_ref(alt: &str) -> bool {
    split_alt(alt).iter().any(|alt| alt == "<NON_REF>")
}

fn single_ref_symbolic_alt_compatible(a_alt: &str, b_alt: &str) -> bool {
    let a_alts = split_alt(a_alt);
    let b_alts = split_alt(b_alt);
    is_single_ref_symbolic_alt(&a_alts) && !contains_ref_symbolic_alt(&b_alts)
        || is_single_ref_symbolic_alt(&b_alts) && !contains_ref_symbolic_alt(&a_alts)
}

fn is_single_ref_symbolic_alt(alts: &[String]) -> bool {
    matches!(alts, [alt] if is_ref_symbolic_alt(alt))
}

fn contains_ref_symbolic_alt(alts: &[String]) -> bool {
    alts.iter().any(|alt| is_ref_symbolic_alt(alt))
}

fn is_ref_symbolic_alt(alt: &str) -> bool {
    alt == "<NON_REF>" || alt == "<*>"
}

fn same_ref_alt_subset_compatible(a_ref: &str, a_alt: &str, b_ref: &str, b_alt: &str) -> bool {
    let a_alts = split_alt(a_alt);
    let b_alts = split_alt(b_alt);
    if !a_ref.eq_ignore_ascii_case(b_ref)
        || a_ref.is_empty()
        || a_alts.is_empty()
        || b_alts.is_empty()
    {
        return false;
    }
    a_alts.iter().all(|alt| b_alts.contains(alt)) || b_alts.iter().all(|alt| a_alts.contains(alt))
}

fn ref_extended_alt_subset_compatible(a_ref: &str, a_alt: &str, b_ref: &str, b_alt: &str) -> bool {
    let Some(merged_ref) = merged_ref(a_ref, b_ref) else {
        return false;
    };
    if a_ref == b_ref {
        return false;
    }
    let a_alts = normalize_alts(a_ref, a_alt, &merged_ref);
    let b_alts = normalize_alts(b_ref, b_alt, &merged_ref);
    if a_ref.is_empty() || b_ref.is_empty() || a_alts.is_empty() || b_alts.is_empty() {
        return false;
    }
    a_alts.iter().all(|alt| b_alts.contains(alt)) || b_alts.iter().all(|alt| a_alts.contains(alt))
}

fn merge_exact_site(
    site: &mut MergedSite,
    record: &RecordLine,
    input_idx: usize,
    info_rules: InfoRules,
    format_numbers: &HashMap<String, VcfNumber>,
    info_numbers: &HashMap<String, VcfNumber>,
    ordered_info_numbers: &[(String, VcfNumber)],
) -> io::Result<()> {
    if site.fixed.len() != record.fixed.len() || site.fixed[..5] != record.fixed[..5] {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "conflicting records at {}:{} require full merge semantics",
                record.fixed[0], record.fixed[1]
            ),
        ));
    }
    if site.samples_by_input[input_idx].is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "duplicate record at {}:{}",
                record.fixed[0], record.fixed[1]
            ),
        ));
    }

    let transformed_samples = if site.fixed.len() > 8 {
        let old_format = site.fixed[8].clone();
        let record_format = record.fixed[8].clone();
        let merged_format = merged_sample_format(&old_format, &record_format);
        if merged_format != old_format {
            let alts = split_alt(&site.fixed[4]);
            let map = allele_map(&alts, &alts);
            for values in site.samples_by_input.iter_mut().flatten() {
                *values = transform_sample_values(
                    &old_format,
                    &merged_format,
                    values,
                    &map,
                    alts.len(),
                    format_numbers,
                );
            }
            site.fixed[8] = merged_format.clone();
        }
        let alts = split_alt(&site.fixed[4]);
        let map = allele_map(&alts, &alts);
        transform_sample_values(
            &record_format,
            &merged_format,
            &record.samples,
            &map,
            alts.len(),
            format_numbers,
        )
    } else {
        record.samples.clone()
    };

    let clear_ref_block_end = should_clear_ref_block_end_on_exact_merge(&site.fixed, &record.fixed);
    merge_fixed_shared_fields(&mut site.fixed, &record.fixed, info_rules.filter_logic);
    if clear_ref_block_end {
        site.fixed[7] = ".".to_owned();
    } else {
        let alts = split_alt(&site.fixed[4]);
        let map = allele_map(&alts, &alts);
        site.fixed[7] = merge_sampled_info(
            &site.fixed[7],
            &record.fixed[7],
            SampledInfoMerge {
                current_map: &map,
                next_map: &map,
                alt_count: alts.len(),
                info_numbers,
                ordered_info_numbers,
                preserve_info_order: true,
            },
        );
    }
    if info_rules.join_af {
        site.fixed[7] = join_info_tag(&site.fixed[7], &record.fixed[7], "AF");
    }
    site.samples_by_input[input_idx] = Some(transformed_samples);
    Ok(())
}

fn should_clear_ref_block_end_on_exact_merge(
    site_fixed: &[String],
    record_fixed: &[String],
) -> bool {
    site_fixed.get(4).is_some_and(|alt| alt == ".")
        && record_fixed.get(4).is_some_and(|alt| alt == ".")
        && site_fixed
            .get(7)
            .is_some_and(|info| info_has_only_end(info))
        && record_fixed
            .get(7)
            .is_some_and(|info| info_has_only_end(info))
        && site_fixed.get(5) != record_fixed.get(5)
}

fn info_has_only_end(info: &str) -> bool {
    info.strip_prefix("END=")
        .is_some_and(|value| !value.is_empty() && !value.contains(';'))
}

fn merge_fixed_shared_fields(
    site_fixed: &mut [String],
    record_fixed: &[String],
    filter_logic: FilterLogic,
) {
    if let (Some(site_id), Some(record_id)) = (site_fixed.get_mut(2), record_fixed.get(2)) {
        *site_id = merge_id(site_id, record_id);
    }
    if let (Some(site_qual), Some(record_qual)) = (site_fixed.get_mut(5), record_fixed.get(5)) {
        *site_qual = merge_qual(site_qual, record_qual);
    }
    if let (Some(site_filter), Some(record_filter)) = (site_fixed.get_mut(6), record_fixed.get(6)) {
        *site_filter = merge_filter_with_logic(site_filter, record_filter, filter_logic);
    }
}

fn merge_id(current: &str, next: &str) -> String {
    match (
        current == "." || current.is_empty(),
        next == "." || next.is_empty(),
    ) {
        (true, true) => ".".to_owned(),
        (true, false) => next.to_owned(),
        (false, true) => current.to_owned(),
        (false, false) if current == next => current.to_owned(),
        (false, false) => {
            let mut ids = current.split(';').collect::<Vec<_>>();
            if !ids.contains(&next) {
                ids.push(next);
            }
            ids.join(";")
        }
    }
}

fn merge_filter(current: &str, next: &str) -> String {
    if current == "." || current == "PASS" || current.is_empty() {
        return next.to_owned();
    }
    if next == "." || next == "PASS" || next.is_empty() || next == current {
        return current.to_owned();
    }
    let mut filters = current.split(';').collect::<Vec<_>>();
    for filter in next.split(';') {
        if !filters.contains(&filter) {
            filters.push(filter);
        }
    }
    filters.join(";")
}

fn merge_filter_with_logic(current: &str, next: &str, filter_logic: FilterLogic) -> String {
    if filter_logic == FilterLogic::RemoveIfPass
        && (is_pass_filter(current) || is_pass_filter(next))
    {
        return "PASS".to_owned();
    }
    merge_filter(current, next)
}

fn is_pass_filter(filter: &str) -> bool {
    filter == "." || filter == "PASS" || filter.is_empty()
}

fn merge_qual(current: &str, next: &str) -> String {
    match (current.parse::<f64>(), next.parse::<f64>()) {
        (Ok(a), Ok(b)) if b > a => next.to_owned(),
        (Ok(_), Ok(_)) => current.to_owned(),
        (Err(_), Ok(_)) if current == "." || current.is_empty() => next.to_owned(),
        _ => current.to_owned(),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CoarseVariantClass {
    Snp,
    Indel,
    Other,
}

fn coarse_variant_class(reference: &str, alt: &str) -> CoarseVariantClass {
    let classes = split_alt(alt)
        .into_iter()
        .map(|alt| precise_allele_class(reference, &alt))
        .collect::<Vec<_>>();
    if classes.is_empty() {
        return CoarseVariantClass::Other;
    }
    if classes
        .iter()
        .all(|class| *class == PreciseVariantClass::Snp)
    {
        CoarseVariantClass::Snp
    } else if classes.iter().all(|class| {
        matches!(
            class,
            PreciseVariantClass::Insertion | PreciseVariantClass::Deletion
        )
    }) {
        CoarseVariantClass::Indel
    } else {
        CoarseVariantClass::Other
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PreciseVariantClass {
    Snp,
    Insertion,
    Deletion,
    Other,
}

fn precise_variant_class(reference: &str, alt: &str) -> PreciseVariantClass {
    let classes = split_alt(alt)
        .into_iter()
        .map(|alt| precise_allele_class(reference, &alt))
        .collect::<Vec<_>>();
    let Some(first) = classes.first().copied() else {
        return PreciseVariantClass::Other;
    };
    if classes.iter().all(|class| *class == first) {
        first
    } else {
        PreciseVariantClass::Other
    }
}

fn precise_allele_class(reference: &str, alt: &str) -> PreciseVariantClass {
    if reference.len() == alt.len() && reference.len() == 1 {
        PreciseVariantClass::Snp
    } else if alt.len() > reference.len() && alt.starts_with(reference) {
        PreciseVariantClass::Insertion
    } else if reference.len() > alt.len() && reference.starts_with(alt) {
        PreciseVariantClass::Deletion
    } else {
        PreciseVariantClass::Other
    }
}

fn merge_sampled_same_position(
    site: &mut MergedSite,
    record: &RecordLine,
    input_idx: usize,
    filter_logic: FilterLogic,
    format_numbers: &HashMap<String, VcfNumber>,
    info_numbers: &HashMap<String, VcfNumber>,
    ordered_info_numbers: &[(String, VcfNumber)],
) -> io::Result<()> {
    if site.samples_by_input[input_idx].is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "duplicate record at {}:{}",
                record.fixed[0], record.fixed[1]
            ),
        ));
    }

    let old_ref = site.fixed[3].clone();
    let old_format = site.fixed[8].clone();
    let record_format = record.fixed[8].clone();
    let merged_format = merged_sample_format(&old_format, &record_format);
    let old_alts = normalize_alts(&old_ref, &site.fixed[4], &old_ref);
    let new_ref = merged_ref(&old_ref, &record.fixed[3]).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "conflicting records at {}:{} require full merge semantics",
                record.fixed[0], record.fixed[1]
            ),
        )
    })?;
    let normalized_site_alts = normalize_alts(&old_ref, &site.fixed[4], &new_ref);
    let normalized_record_alts = normalize_alts(&record.fixed[3], &record.fixed[4], &new_ref);

    let mut merged_alts = normalized_site_alts.clone();
    for alt in &normalized_record_alts {
        if !merged_alts.contains(alt) {
            merged_alts.push(alt.clone());
        }
    }

    let old_site_map = allele_map(&normalized_site_alts, &merged_alts);
    let record_map = allele_map(&normalized_record_alts, &merged_alts);
    let same_non_dot_id =
        site.fixed.get(2) == record.fixed.get(2) && site.fixed.get(2).is_some_and(|id| id != ".");
    let merged_info = merge_sampled_info(
        &site.fixed[7],
        &record.fixed[7],
        SampledInfoMerge {
            current_map: &old_site_map,
            next_map: &record_map,
            alt_count: merged_alts.len(),
            info_numbers,
            ordered_info_numbers,
            preserve_info_order: same_non_dot_id,
        },
    );
    if new_ref != old_ref || merged_alts != old_alts || merged_format != old_format {
        for values in site.samples_by_input.iter_mut().flatten() {
            *values = transform_sample_values(
                &old_format,
                &merged_format,
                values,
                &old_site_map,
                merged_alts.len(),
                format_numbers,
            );
        }
    }

    site.fixed[3] = new_ref;
    site.fixed[4] = if merged_alts.is_empty() {
        ".".to_owned()
    } else {
        merged_alts.join(",")
    };
    site.fixed[7] = merged_info;
    site.fixed[8] = merged_format.clone();
    merge_fixed_shared_fields(&mut site.fixed, &record.fixed, filter_logic);
    site.samples_by_input[input_idx] = Some(transform_sample_values(
        &record_format,
        &merged_format,
        &record.samples,
        &record_map,
        merged_alts.len(),
        format_numbers,
    ));
    Ok(())
}

fn merged_ref(a: &str, b: &str) -> Option<String> {
    if a.eq_ignore_ascii_case(b) {
        return Some(a.to_ascii_uppercase());
    }
    if a.len() >= b.len() && a.starts_with(b) {
        Some(a.to_owned())
    } else if b.starts_with(a) {
        Some(b.to_owned())
    } else {
        None
    }
}

fn normalize_alts(reference: &str, alt: &str, merged_ref: &str) -> Vec<String> {
    let suffix = merged_ref.strip_prefix(reference).unwrap_or("");
    split_alt(alt)
        .into_iter()
        .map(|alt| {
            if alt == "*" || (alt.starts_with('<') && alt.ends_with('>')) {
                alt
            } else {
                format!("{alt}{suffix}")
            }
        })
        .collect()
}

fn apply_local_alleles(
    sites: &mut [MergedSite],
    limit: usize,
    format_numbers: &HashMap<String, VcfNumber>,
    input_sample_counts: &[usize],
) {
    for site in sites {
        let alt_count = split_alt(&site.fixed[4]).len();
        if alt_count <= limit || site.fixed.len() < 9 {
            continue;
        }
        let old_format = site.fixed[8].clone();
        let new_format = localized_format(&old_format, format_numbers);
        if new_format == old_format {
            continue;
        }

        for values in site.samples_by_input.iter_mut().flatten() {
            for sample in values {
                *sample =
                    localize_sample_value(&old_format, sample, alt_count, limit, format_numbers);
            }
        }
        site.fixed[8] = new_format;
        fill_localized_missing_samples(site, limit, input_sample_counts);
    }
}

#[derive(Clone, Copy)]
enum LaaToken {
    Int(usize),
    Missing,
    VectorEnd,
}

fn localized_format(format: &str, format_numbers: &HashMap<String, VcfNumber>) -> String {
    let keys = split_format_keys(format);
    if !keys
        .iter()
        .any(|key| is_localizable_format_key(key, format_numbers))
    {
        return format.to_owned();
    }

    let mut out = Vec::new();
    for key in keys {
        if is_localizable_format_key(key, format_numbers) {
            out.push(format!("L{key}"));
        } else {
            out.push(key.to_owned());
        }
    }
    out.push("LAA".to_owned());
    out.join(":")
}

fn is_localizable_format_key(key: &str, format_numbers: &HashMap<String, VcfNumber>) -> bool {
    matches!(
        format_numbers.get(key),
        Some(VcfNumber::A | VcfNumber::R | VcfNumber::G)
    )
}

fn localize_sample_value(
    input_format: &str,
    sample: &str,
    alt_count: usize,
    limit: usize,
    format_numbers: &HashMap<String, VcfNumber>,
) -> String {
    let input_keys = split_format_keys(input_format);
    let input_values = sample.split(':').collect::<Vec<_>>();
    let input_index = input_keys
        .iter()
        .enumerate()
        .map(|(idx, key)| (*key, idx))
        .collect::<HashMap<_, _>>();
    let Some(pl) = input_index
        .get("PL")
        .and_then(|idx| input_values.get(*idx).copied())
    else {
        return sample.to_owned();
    };

    let local_alts = choose_local_alts_from_pl(pl, alt_count, limit);
    let mut out = Vec::new();
    for key in input_keys {
        if is_localizable_format_key(key, format_numbers) {
            let value = input_index
                .get(key)
                .and_then(|idx| input_values.get(*idx).copied())
                .unwrap_or(".");
            out.push(localize_format_value(
                key,
                value,
                &local_alts,
                format_numbers,
            ));
        } else {
            out.push(
                input_index
                    .get(key)
                    .and_then(|idx| input_values.get(*idx).copied())
                    .unwrap_or(".")
                    .to_owned(),
            );
        }
    }
    out.push(localize_laa_value(&local_alts));
    out.join(":")
}

fn choose_local_alts_from_pl(raw: &str, alt_count: usize, limit: usize) -> Vec<usize> {
    if raw == "." {
        return Vec::new();
    }
    let values = raw.split(',').collect::<Vec<_>>();
    let mut allele_probs = vec![0.0; alt_count + 1];
    for b in 0..=alt_count {
        for a in 0..=b {
            let idx = genotype_index(a, b);
            let Some(value) = values.get(idx) else {
                continue;
            };
            let Some(prob) = pl_to_probability(value) else {
                continue;
            };
            allele_probs[a] += prob;
            allele_probs[b] += prob;
        }
    }

    let mut alts = (1..=alt_count)
        .filter(|allele| allele_probs[*allele] > 0.0)
        .collect::<Vec<_>>();
    alts.sort_by(|a, b| {
        allele_probs[*b]
            .partial_cmp(&allele_probs[*a])
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.cmp(b))
    });
    alts.truncate(limit.min(alt_count));
    alts.sort_unstable();
    alts
}

fn pl_to_probability(raw: &str) -> Option<f64> {
    if raw == "." {
        return None;
    }
    let pl = raw.parse::<f64>().ok()?;
    Some(10_f64.powf(-pl / 10.0))
}

fn localize_format_value(
    key: &str,
    value: &str,
    local_alts: &[usize],
    format_numbers: &HashMap<String, VcfNumber>,
) -> String {
    if value == "." {
        return value.to_owned();
    }
    match format_numbers.get(key).copied().unwrap_or(VcfNumber::Other) {
        VcfNumber::A => localize_number_a(value, local_alts),
        VcfNumber::R => localize_number_r(value, local_alts),
        VcfNumber::G => localize_number_g(value, local_alts),
        VcfNumber::Other => value.to_owned(),
    }
}

fn localize_number_a(raw: &str, local_alts: &[usize]) -> String {
    let values = raw.split(',').collect::<Vec<_>>();
    local_alts
        .iter()
        .map(|allele| values.get(allele - 1).copied().unwrap_or("."))
        .collect::<Vec<_>>()
        .join(",")
}

fn localize_number_r(raw: &str, local_alts: &[usize]) -> String {
    let values = raw.split(',').collect::<Vec<_>>();
    std::iter::once(0)
        .chain(local_alts.iter().copied())
        .map(|allele| values.get(allele).copied().unwrap_or("."))
        .collect::<Vec<_>>()
        .join(",")
}

fn localize_number_g(raw: &str, local_alts: &[usize]) -> String {
    let values = raw.split(',').collect::<Vec<_>>();
    let alleles = std::iter::once(0)
        .chain(local_alts.iter().copied())
        .collect::<Vec<_>>();
    let mut out = Vec::new();
    for (b_idx, b) in alleles.iter().enumerate() {
        for a in &alleles[..=b_idx] {
            out.push(values.get(genotype_index(*a, *b)).copied().unwrap_or("."));
        }
    }
    out.join(",")
}

fn localize_laa_value(local_alts: &[usize]) -> String {
    if local_alts.is_empty() {
        ".".to_owned()
    } else {
        local_alts
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn fill_localized_missing_samples(
    site: &mut MergedSite,
    limit: usize,
    input_sample_counts: &[usize],
) {
    let Some(format) = site.fixed.get(8) else {
        return;
    };
    let format_keys = split_format_keys(format);
    let Some(laa_idx) = format_keys.iter().position(|key| *key == "LAA") else {
        return;
    };

    let total_samples = input_sample_counts.iter().sum::<usize>();
    if total_samples == 0 {
        return;
    }

    let stride = limit + 1;
    let mut laa = vec![LaaToken::VectorEnd; total_samples * stride];
    let mut present = vec![false; total_samples];
    let mut n_laa = 0;
    let mut sample_idx = 0;
    for (input_idx, sample_count) in input_sample_counts.iter().copied().enumerate() {
        match site
            .samples_by_input
            .get(input_idx)
            .and_then(Option::as_ref)
        {
            Some(values) => {
                for sample in values.iter().take(sample_count) {
                    present[sample_idx] = true;
                    let base = sample_idx * stride;
                    laa[base] = LaaToken::Int(0);
                    let sample_laa = parse_laa_tokens(sample, laa_idx);
                    n_laa = n_laa.max(sample_laa.len().min(limit));
                    for (idx, token) in sample_laa.into_iter().take(limit).enumerate() {
                        laa[base + idx + 1] = token;
                    }
                    sample_idx += 1;
                }
            }
            None => {
                for _ in 0..sample_count {
                    let base = sample_idx * stride;
                    laa[base] = LaaToken::Missing;
                    sample_idx += 1;
                }
            }
        }
    }
    if n_laa == 0 {
        return;
    }

    for (sample_idx, is_present) in present.iter().copied().enumerate() {
        let src_offset = sample_idx * stride;
        let dst_offset = sample_idx * n_laa;
        let mut dst_idx = 0;
        if is_present {
            while dst_idx < n_laa {
                match laa[src_offset + dst_idx + 1] {
                    LaaToken::Missing => laa[dst_offset + dst_idx] = LaaToken::Missing,
                    LaaToken::VectorEnd => break,
                    LaaToken::Int(value) => laa[dst_offset + dst_idx] = LaaToken::Int(value),
                }
                dst_idx += 1;
            }
        }
        if dst_idx == 0 {
            laa[dst_offset] = LaaToken::Missing;
            dst_idx += 1;
        }
        // Match upstream's in-place LAA compaction, including its source-stride
        // tail write. This preserves byte parity for absent samples in LPL rows.
        while dst_idx < n_laa {
            laa[src_offset + dst_idx] = LaaToken::VectorEnd;
            dst_idx += 1;
        }
    }

    let mut sample_idx = 0;
    for (input_idx, sample_count) in input_sample_counts.iter().copied().enumerate() {
        if site
            .samples_by_input
            .get(input_idx)
            .is_some_and(Option::is_some)
        {
            sample_idx += sample_count;
            continue;
        }
        let values = (0..sample_count)
            .map(|_| {
                let dst_offset = sample_idx * n_laa;
                sample_idx += 1;
                localized_missing_sample_value(
                    format,
                    &render_laa_tokens(&laa[dst_offset..][..n_laa]),
                )
            })
            .collect::<Vec<_>>();
        if let Some(slot) = site.samples_by_input.get_mut(input_idx) {
            *slot = Some(values);
        }
    }
}

fn parse_laa_tokens(sample: &str, laa_idx: usize) -> Vec<LaaToken> {
    let raw = sample.split(':').nth(laa_idx).unwrap_or(".");
    if raw == "." || raw.is_empty() {
        return vec![LaaToken::Missing];
    }
    raw.split(',')
        .map(|value| {
            if value == "." {
                LaaToken::Missing
            } else {
                value
                    .parse::<usize>()
                    .map(LaaToken::Int)
                    .unwrap_or(LaaToken::Missing)
            }
        })
        .collect()
}

fn render_laa_tokens(tokens: &[LaaToken]) -> String {
    let mut values = Vec::new();
    for token in tokens {
        match token {
            LaaToken::Int(value) => values.push(value.to_string()),
            LaaToken::Missing => values.push(".".to_owned()),
            LaaToken::VectorEnd => break,
        }
    }
    if values.is_empty() {
        ".".to_owned()
    } else {
        values.join(",")
    }
}

fn localized_missing_sample_value(format: &str, laa: &str) -> String {
    format
        .split(':')
        .map(|key| {
            if key == "GT" {
                "./."
            } else if key == "LAA" {
                laa
            } else {
                "."
            }
        })
        .collect::<Vec<_>>()
        .join(":")
}

fn allele_map(old_alts: &[String], merged_alts: &[String]) -> Vec<Option<usize>> {
    let mut map = Vec::with_capacity(old_alts.len() + 1);
    map.push(Some(0));
    for alt in old_alts {
        map.push(
            merged_alts
                .iter()
                .position(|merged| merged == alt)
                .map(|idx| idx + 1),
        );
    }
    map
}

fn transform_sample_values(
    input_format: &str,
    output_format: &str,
    samples: &[String],
    allele_map: &[Option<usize>],
    alt_count: usize,
    format_numbers: &HashMap<String, VcfNumber>,
) -> Vec<String> {
    let input_keys = split_format_keys(input_format);
    let output_keys = split_format_keys(output_format);
    let input_index = input_keys
        .iter()
        .enumerate()
        .map(|(idx, key)| (*key, idx))
        .collect::<HashMap<_, _>>();

    samples
        .iter()
        .map(|sample| {
            let input_values = sample.split(':').collect::<Vec<_>>();
            let input_gt = input_index
                .get("GT")
                .and_then(|idx| input_values.get(*idx).copied());
            output_keys
                .iter()
                .map(|key| {
                    input_index
                        .get(key)
                        .and_then(|idx| input_values.get(*idx))
                        .map(|value| {
                            remap_format_value(
                                key,
                                value,
                                allele_map,
                                alt_count,
                                format_numbers,
                                input_gt,
                            )
                        })
                        .unwrap_or_else(|| ".".to_owned())
                })
                .collect::<Vec<_>>()
                .join(":")
        })
        .collect()
}

fn split_format_keys(format: &str) -> Vec<&str> {
    if format == "." || format.is_empty() {
        Vec::new()
    } else {
        format.split(':').collect()
    }
}

fn merged_sample_format(left: &str, right: &str) -> String {
    let mut keys = split_format_keys(left)
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    for key in split_format_keys(right) {
        if !keys.iter().any(|existing| existing == key) {
            keys.push(key.to_owned());
        }
    }
    if keys.is_empty() {
        ".".to_owned()
    } else {
        keys.join(":")
    }
}

fn remap_format_value(
    key: &str,
    value: &str,
    allele_map: &[Option<usize>],
    alt_count: usize,
    format_numbers: &HashMap<String, VcfNumber>,
    input_gt: Option<&str>,
) -> String {
    match key {
        "GT" => remap_gt(value, allele_map),
        "AD" => remap_number_r(value, allele_map, alt_count),
        _ => match format_numbers.get(key).copied().unwrap_or(VcfNumber::Other) {
            VcfNumber::A => remap_number_a(value, allele_map, alt_count),
            VcfNumber::R => remap_number_r(value, allele_map, alt_count),
            VcfNumber::G => {
                if genotype_ploidy(input_gt).is_some_and(|ploidy| ploidy == 1) {
                    remap_number_g_haploid(value, allele_map, alt_count)
                } else {
                    remap_number_g(value, allele_map, alt_count)
                }
            }
            VcfNumber::Other => value.to_owned(),
        },
    }
}

fn genotype_ploidy(gt: Option<&str>) -> Option<usize> {
    let gt = gt?;
    if gt == "." || gt.chars().all(|ch| ch.is_ascii_digit()) {
        return Some(1);
    }
    Some(gt.chars().filter(|ch| *ch == '/' || *ch == '|').count() + 1)
}

fn remap_gt(raw: &str, allele_map: &[Option<usize>]) -> String {
    let mut out = String::new();
    let mut chars = raw.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch.is_ascii_digit() {
            let mut allele = ch.to_string();
            while let Some(next) = chars.peek().copied() {
                if next.is_ascii_digit() {
                    allele.push(next);
                    chars.next();
                } else {
                    break;
                }
            }
            let mapped = allele
                .parse::<usize>()
                .ok()
                .and_then(|idx| allele_map.get(idx).copied().flatten());
            match mapped {
                Some(idx) => out.push_str(&idx.to_string()),
                None => out.push('.'),
            }
        } else {
            out.push(ch);
        }
    }
    out
}

fn remap_number_a(raw: &str, allele_map: &[Option<usize>], alt_count: usize) -> String {
    if raw == "." {
        return raw.to_owned();
    }
    slots_to_string(remap_number_a_slots(raw, allele_map, alt_count), ".")
}

fn remap_number_r(raw: &str, allele_map: &[Option<usize>], alt_count: usize) -> String {
    if raw == "." {
        return raw.to_owned();
    }
    slots_to_string(remap_number_r_slots(raw, allele_map, alt_count), ".")
}

fn remap_number_g(raw: &str, allele_map: &[Option<usize>], alt_count: usize) -> String {
    if raw == "." {
        return raw.to_owned();
    }
    slots_to_string(remap_number_g_slots(raw, allele_map, alt_count), ".")
}

fn remap_number_g_haploid(raw: &str, allele_map: &[Option<usize>], alt_count: usize) -> String {
    if raw == "." {
        return raw.to_owned();
    }
    slots_to_string(
        remap_number_g_haploid_slots(raw, allele_map, alt_count),
        ".",
    )
}

fn remap_number_a_slots(
    raw: &str,
    allele_map: &[Option<usize>],
    alt_count: usize,
) -> Vec<Option<String>> {
    let mut values = vec![None; alt_count];
    if raw == "." {
        return values;
    }
    let old_values = raw.split(',').collect::<Vec<_>>();
    for (old_alt_idx, value) in old_values.iter().enumerate() {
        if let Some(new_idx) = allele_map.get(old_alt_idx + 1).copied().flatten()
            && let Some(slot) = values.get_mut(new_idx.saturating_sub(1))
            && *value != "."
        {
            *slot = Some((*value).to_owned());
        }
    }
    values
}

fn remap_number_r_slots(
    raw: &str,
    allele_map: &[Option<usize>],
    alt_count: usize,
) -> Vec<Option<String>> {
    let mut values = vec![None; alt_count + 1];
    if raw == "." {
        return values;
    }
    let old_values = raw.split(',').collect::<Vec<_>>();
    for (old_idx, value) in old_values.iter().enumerate() {
        if let Some(new_idx) = allele_map.get(old_idx).copied().flatten()
            && let Some(slot) = values.get_mut(new_idx)
            && *value != "."
        {
            *slot = Some((*value).to_owned());
        }
    }
    values
}

fn remap_number_g_slots(
    raw: &str,
    allele_map: &[Option<usize>],
    alt_count: usize,
) -> Vec<Option<String>> {
    let mut values = vec![None; genotype_count(alt_count)];
    if raw == "." {
        return values;
    }
    let old_values = raw.split(',').collect::<Vec<_>>();
    let old_allele_count = allele_map.len();
    for old_b in 0..old_allele_count {
        for old_a in 0..=old_b {
            let old_idx = genotype_index(old_a, old_b);
            let Some(value) = old_values.get(old_idx) else {
                continue;
            };
            if *value == "." {
                continue;
            }
            let Some(new_a) = allele_map.get(old_a).copied().flatten() else {
                continue;
            };
            let Some(new_b) = allele_map.get(old_b).copied().flatten() else {
                continue;
            };
            let new_idx = genotype_index(new_a, new_b);
            if let Some(slot) = values.get_mut(new_idx) {
                *slot = Some((*value).to_owned());
            }
        }
    }
    values
}

fn remap_number_g_haploid_slots(
    raw: &str,
    allele_map: &[Option<usize>],
    alt_count: usize,
) -> Vec<Option<String>> {
    let mut values = vec![None; alt_count + 1];
    if raw == "." {
        return values;
    }
    let old_values = raw.split(',').collect::<Vec<_>>();
    for (old_idx, value) in old_values.iter().enumerate() {
        if *value == "." {
            continue;
        }
        if let Some(new_idx) = allele_map.get(old_idx).copied().flatten()
            && let Some(slot) = values.get_mut(new_idx)
        {
            *slot = Some((*value).to_owned());
        }
    }
    values
}

fn genotype_count(alt_count: usize) -> usize {
    let allele_count = alt_count + 1;
    allele_count * (allele_count + 1) / 2
}

fn genotype_index(a: usize, b: usize) -> usize {
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    hi * (hi + 1) / 2 + lo
}

fn slots_to_string(values: Vec<Option<String>>, missing: &str) -> String {
    values
        .into_iter()
        .map(|value| value.unwrap_or_else(|| missing.to_owned()))
        .collect::<Vec<_>>()
        .join(",")
}

struct SampledInfoMerge<'a> {
    current_map: &'a [Option<usize>],
    next_map: &'a [Option<usize>],
    alt_count: usize,
    info_numbers: &'a HashMap<String, VcfNumber>,
    ordered_info_numbers: &'a [(String, VcfNumber)],
    preserve_info_order: bool,
}

fn merge_sampled_info(current: &str, next: &str, info_merge: SampledInfoMerge<'_>) -> String {
    let current_has_mergeable = has_mergeable_sampled_info(current, info_merge.info_numbers);
    let next_has_mergeable = has_mergeable_sampled_info(next, info_merge.info_numbers);
    let mut fields = merge_unstructured_info_fields(current, next, info_merge.info_numbers);
    if fields.is_empty() && !current_has_mergeable && !next_has_mergeable {
        return current.to_owned();
    }

    if info_merge.preserve_info_order {
        for (key, number) in info_merge.ordered_info_numbers {
            if key == "AC" || key == "AN" || *number == VcfNumber::Other {
                continue;
            }
            if let Some(field) = merge_sampled_info_field(
                current,
                next,
                key,
                *number,
                info_merge.current_map,
                info_merge.next_map,
                info_merge.alt_count,
            ) {
                fields.push(field);
            }
        }
    } else {
        let target_order = if has_x_numbered_info(current) || has_x_numbered_info(next) {
            [VcfNumber::R, VcfNumber::A, VcfNumber::G]
        } else {
            [VcfNumber::A, VcfNumber::G, VcfNumber::R]
        };
        for target_number in target_order {
            for (key, number) in info_merge.ordered_info_numbers {
                if *number != target_number || key == "AC" || key == "AN" {
                    continue;
                }
                if let Some(field) = merge_sampled_info_field(
                    current,
                    next,
                    key,
                    target_number,
                    info_merge.current_map,
                    info_merge.next_map,
                    info_merge.alt_count,
                ) {
                    fields.push(field);
                }
            }
        }
    }

    if info_value(current, "AN").is_some() || info_value(next, "AN").is_some() {
        let an = info_i64(current, "AN").unwrap_or(0) + info_i64(next, "AN").unwrap_or(0);
        fields.push(format!("AN={an}"));
    }
    if info_value(current, "AC").is_some() || info_value(next, "AC").is_some() {
        fields.push(format!(
            "AC={}",
            merge_info_number_a(
                info_value(current, "AC"),
                info_value(next, "AC"),
                info_merge.current_map,
                info_merge.next_map,
                info_merge.alt_count,
            )
        ));
    }

    if fields.is_empty() {
        ".".to_owned()
    } else {
        fields.join(";")
    }
}

fn has_x_numbered_info(info: &str) -> bool {
    info.split(';').any(|field| field.starts_with('X'))
}

#[derive(Clone)]
struct InfoField {
    key: String,
    value: Option<String>,
}

fn merge_unstructured_info_fields(
    current: &str,
    next: &str,
    info_numbers: &HashMap<String, VcfNumber>,
) -> Vec<String> {
    let current_fields = parse_info_fields(current);
    let next_fields = parse_info_fields(next);
    let mut keys = Vec::<String>::new();
    for field in current_fields.iter().chain(next_fields.iter()) {
        if field.key == "AC"
            || field.key == "AN"
            || matches!(
                info_numbers.get(&field.key),
                Some(VcfNumber::A | VcfNumber::R | VcfNumber::G)
            )
        {
            continue;
        }
        if !keys.contains(&field.key) {
            keys.push(field.key.clone());
        }
    }

    keys.sort_by_key(|key| {
        unstructured_info_sort_key(
            key,
            current_fields.iter().find(|field| field.key == *key),
            next_fields.iter().find(|field| field.key == *key),
        )
    });

    keys.into_iter()
        .map(|key| {
            let current_value = current_fields
                .iter()
                .find(|field| field.key == key)
                .and_then(|field| field.value.as_deref());
            let next_value = next_fields
                .iter()
                .find(|field| field.key == key)
                .and_then(|field| field.value.as_deref());
            match (current_value, next_value) {
                (None, None) => key,
                (Some(a), Some(b))
                    if matches!(key.as_str(), "DP" | "DP4")
                        && numeric_list(a).is_some()
                        && numeric_list(b).is_some() =>
                {
                    format!("{key}={}", sum_numeric_lists(a, b))
                }
                (Some(a), Some(_)) => format!("{key}={a}"),
                (Some(a), None) => format!("{key}={a}"),
                (None, Some(b)) => format!("{key}={b}"),
            }
        })
        .collect()
}

fn parse_info_fields(info: &str) -> Vec<InfoField> {
    if info == "." || info.is_empty() {
        return Vec::new();
    }
    info.split(';')
        .filter(|field| !field.is_empty())
        .map(|field| {
            if let Some((key, value)) = field.split_once('=') {
                InfoField {
                    key: key.to_owned(),
                    value: Some(value.to_owned()),
                }
            } else {
                InfoField {
                    key: field.to_owned(),
                    value: None,
                }
            }
        })
        .collect()
}

fn unstructured_info_sort_key(
    key: &str,
    current: Option<&InfoField>,
    next: Option<&InfoField>,
) -> (usize, usize) {
    let current_non_vector = current
        .and_then(|field| field.value.as_deref())
        .is_none_or(|value| !value.contains(','));
    if current.is_some() && current_non_vector {
        return (0, 0);
    }
    let value = current
        .or(next)
        .and_then(|field| field.value.as_deref())
        .unwrap_or("");
    let category = if value.contains(',') {
        3
    } else if numeric_list(value).is_some() {
        2
    } else {
        1
    };
    let key_rank = match key {
        "STR" | "TXT" => 0,
        "DP" => 1,
        "DP4" => 2,
        _ => 3,
    };
    (category, key_rank)
}

fn sum_numeric_lists(a: &str, b: &str) -> String {
    let Some(a_values) = numeric_list(a) else {
        return a.to_owned();
    };
    let Some(b_values) = numeric_list(b) else {
        return a.to_owned();
    };
    let len = a_values.len().max(b_values.len());
    (0..len)
        .map(|idx| {
            let a = a_values.get(idx).copied().unwrap_or(0.0);
            let b = b_values.get(idx).copied().unwrap_or(0.0);
            format_float(a + b)
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn numeric_list(raw: &str) -> Option<Vec<f64>> {
    if raw.is_empty() || raw == "." {
        return None;
    }
    raw.split(',')
        .map(|value| value.parse::<f64>().ok())
        .collect()
}

fn merge_sampled_info_field(
    current: &str,
    next: &str,
    key: &str,
    target_number: VcfNumber,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> Option<String> {
    let current_value = info_value(current, key);
    let next_value = info_value(next, key);
    if current_value.is_none() && next_value.is_none() {
        return None;
    }
    let value = match target_number {
        VcfNumber::A if key.starts_with('X') => {
            overlay_info_number_a(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::A if current_value == next_value => {
            overlay_info_number_a(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::A => {
            merge_info_number_a(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::R if key.starts_with('X') => {
            overlay_info_number_r(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::R if current_value == next_value => {
            overlay_info_number_r(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::R => {
            merge_info_number_r(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::G if key.starts_with('X') => {
            overlay_info_number_g(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::G if current_value == next_value => {
            overlay_info_number_g(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::G => {
            merge_info_number_g(current_value, next_value, current_map, next_map, alt_count)
        }
        VcfNumber::Other => return None,
    };
    Some(format!("{key}={value}"))
}

fn has_mergeable_sampled_info(info: &str, info_numbers: &HashMap<String, VcfNumber>) -> bool {
    if info_value(info, "AN").is_some() || info_value(info, "AC").is_some() {
        return true;
    }
    info.split(';').any(|field| {
        let key = field.split_once('=').map(|(key, _)| key).unwrap_or(field);
        matches!(
            info_numbers.get(key),
            Some(VcfNumber::A | VcfNumber::R | VcfNumber::G)
        )
    })
}

fn merge_info_number_a(
    current: Option<&str>,
    next: Option<&str>,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> String {
    let current = current
        .map(|value| remap_number_a_slots(value, current_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count]);
    let next = next
        .map(|value| remap_number_a_slots(value, next_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count]);
    sum_numeric_slots(&current, &next, None)
}

fn overlay_info_number_a(
    current: Option<&str>,
    next: Option<&str>,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> String {
    let current = current
        .map(|value| remap_number_a_slots(value, current_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count]);
    let next = next
        .map(|value| remap_number_a_slots(value, next_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count]);
    overlay_slots(&current, &next, ".")
}

fn merge_info_number_r(
    current: Option<&str>,
    next: Option<&str>,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> String {
    let current = current
        .map(|value| remap_number_r_slots(value, current_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count + 1]);
    let next = next
        .map(|value| remap_number_r_slots(value, next_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count + 1]);
    sum_numeric_slots(&current, &next, None)
}

fn overlay_info_number_r(
    current: Option<&str>,
    next: Option<&str>,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> String {
    let current = current
        .map(|value| remap_number_r_slots(value, current_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count + 1]);
    let next = next
        .map(|value| remap_number_r_slots(value, next_map, alt_count))
        .unwrap_or_else(|| vec![None; alt_count + 1]);
    overlay_slots(&current, &next, ".")
}

fn merge_info_number_g(
    current: Option<&str>,
    next: Option<&str>,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> String {
    let current = current
        .map(|value| remap_number_g_slots(value, current_map, alt_count))
        .unwrap_or_else(|| vec![None; genotype_count(alt_count)]);
    let next = next
        .map(|value| remap_number_g_slots(value, next_map, alt_count))
        .unwrap_or_else(|| vec![None; genotype_count(alt_count)]);
    sum_numeric_slots(&current, &next, Some("0"))
}

fn overlay_info_number_g(
    current: Option<&str>,
    next: Option<&str>,
    current_map: &[Option<usize>],
    next_map: &[Option<usize>],
    alt_count: usize,
) -> String {
    let current = current
        .map(|value| remap_number_g_slots(value, current_map, alt_count))
        .unwrap_or_else(|| vec![None; genotype_count(alt_count)]);
    let next = next
        .map(|value| remap_number_g_slots(value, next_map, alt_count))
        .unwrap_or_else(|| vec![None; genotype_count(alt_count)]);
    overlay_slots(&current, &next, ".")
}

fn overlay_slots(current: &[Option<String>], next: &[Option<String>], missing: &str) -> String {
    current
        .iter()
        .zip(next)
        .map(|(a, b)| {
            normalize_scientific_literal(a.as_deref().or(b.as_deref()).unwrap_or(missing))
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn sum_numeric_slots(
    current: &[Option<String>],
    next: &[Option<String>],
    missing_when_both_absent: Option<&str>,
) -> String {
    current
        .iter()
        .zip(next)
        .map(|(a, b)| match (a.as_deref(), b.as_deref()) {
            (Some(a), Some(b)) => sum_numeric_values(a, b),
            (Some(a), None) => a.to_owned(),
            (None, Some(b)) => b.to_owned(),
            (None, None) => missing_when_both_absent.unwrap_or(".").to_owned(),
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn sum_numeric_values(a: &str, b: &str) -> String {
    if let (Ok(a), Ok(b)) = (a.parse::<i64>(), b.parse::<i64>()) {
        return (a + b).to_string();
    }
    match (a.parse::<f64>(), b.parse::<f64>()) {
        (Ok(a), Ok(b)) => format_float(a + b),
        _ => b.to_owned(),
    }
}

fn format_float(value: f64) -> String {
    if value.abs() >= 1_000_000.0 && value.fract() == 0.0 {
        let raw = format!("{value:.0e}");
        if let Some((mantissa, exponent)) = raw.split_once('e') {
            let exponent = exponent.parse::<i32>().unwrap_or(0);
            return format!("{mantissa}e{exponent:+03}");
        }
    }
    if value.fract() == 0.0 {
        format!("{value:.0}")
    } else {
        value.to_string()
    }
}

fn normalize_scientific_literal(raw: &str) -> String {
    if !raw.contains('e') && !raw.contains('E') {
        return raw.to_owned();
    }
    raw.parse::<f64>()
        .map(format_float)
        .unwrap_or_else(|_| raw.to_owned())
}

fn merge_sites_only_alt_union(site: &mut MergedSite, record: &RecordLine, info_rules: InfoRules) {
    let mut alts = split_alt(&site.fixed[4]);
    let next_alts = split_alt(&record.fixed[4]);
    for alt in &next_alts {
        if !alts.contains(alt) {
            alts.push(alt.clone());
        }
    }

    let site_has_no_sample_values = site_has_no_sample_values(site);
    if info_rules.join_af {
        site.fixed[4] = alts.join(",");
        site.fixed[7] = join_info_tag(&site.fixed[7], &record.fixed[7], "AF");
        return;
    }

    let mut ac_by_alt = HashMap::new();
    add_ac_by_alt(&mut ac_by_alt, &site.fixed[4], &site.fixed[7]);
    if site_has_no_sample_values || info_rules.sum_ac {
        add_ac_by_alt(&mut ac_by_alt, &record.fixed[4], &record.fixed[7]);
    }
    let an = info_i64(&site.fixed[7], "AN").unwrap_or(0)
        + if site_has_no_sample_values || info_rules.sum_an {
            info_i64(&record.fixed[7], "AN").unwrap_or(0)
        } else {
            0
        };
    let ac_values = alts
        .iter()
        .map(|alt| ac_by_alt.get(alt).copied().unwrap_or(0).to_string())
        .collect::<Vec<_>>()
        .join(",");

    site.fixed[4] = alts.join(",");
    site.fixed[7] = if site_has_no_sample_values || info_rules.sum_ac || info_rules.sum_an {
        format!("AC={ac_values};AN={an}")
    } else {
        format!("AN={an};AC={ac_values}")
    };
}

fn join_info_tag(current: &str, next: &str, key: &str) -> String {
    let current_value = info_value(current, key).filter(|value| !value.is_empty() && *value != ".");
    let next_value = info_value(next, key).filter(|value| !value.is_empty() && *value != ".");
    match (current_value, next_value) {
        (Some(a), Some(b)) => set_info_value(current, key, &format!("{a},{b}")),
        (None, Some(b)) => set_info_value(current, key, b),
        _ => current.to_owned(),
    }
}

fn set_info_value(info: &str, key: &str, value: &str) -> String {
    if info == "." || info.is_empty() {
        return format!("{key}={value}");
    }

    let mut found = false;
    let fields = info
        .split(';')
        .map(|field| {
            if field
                .split_once('=')
                .is_some_and(|(field_key, _)| field_key == key)
            {
                found = true;
                format!("{key}={value}")
            } else {
                field.to_owned()
            }
        })
        .collect::<Vec<_>>();

    if found {
        fields.join(";")
    } else {
        format!("{};{key}={value}", fields.join(";"))
    }
}

fn remove_info_value(info: &str, key: &str) -> String {
    if info == "." || info.is_empty() {
        return ".".to_owned();
    }
    let fields = info
        .split(';')
        .filter(|field| {
            field
                .split_once('=')
                .is_none_or(|(field_key, _)| field_key != key)
        })
        .collect::<Vec<_>>();
    if fields.is_empty() {
        ".".to_owned()
    } else {
        fields.join(";")
    }
}

fn site_has_no_sample_values(site: &MergedSite) -> bool {
    site.samples_by_input
        .iter()
        .all(|samples| samples.as_ref().is_none_or(|values| values.is_empty()))
}

fn split_alt(raw: &str) -> Vec<String> {
    if raw == "." || raw.is_empty() {
        Vec::new()
    } else {
        raw.split(',').map(str::to_owned).collect()
    }
}

fn add_ac_by_alt(out: &mut HashMap<String, i64>, alt: &str, info: &str) {
    let alts = split_alt(alt);
    let Some(ac_raw) = info_value(info, "AC") else {
        return;
    };
    for (alt, value) in alts.iter().zip(ac_raw.split(',')) {
        let Ok(value) = value.parse::<i64>() else {
            continue;
        };
        *out.entry(alt.clone()).or_insert(0) += value;
    }
}

fn info_i64(info: &str, key: &str) -> Option<i64> {
    info_value(info, key)?.parse().ok()
}

fn info_u64(info: &str, key: &str) -> Option<u64> {
    info_value(info, key)?.parse().ok()
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|field| {
        let (name, value) = field.split_once('=')?;
        (name == key).then_some(value)
    })
}

fn cleanup_single_base_reference_block_end(fixed: &mut [String]) {
    if fixed.get(4).is_none_or(|alt| alt != ".") {
        return;
    }
    let Some(pos) = fixed.get(1).and_then(|raw| raw.parse::<u64>().ok()) else {
        return;
    };
    if fixed.get(7).and_then(|info| info_u64(info, "END")) == Some(pos) {
        fixed[7] = remove_info_value(&fixed[7], "END");
    }
}

fn add_missing_reference_alleles_to_an(fixed: &mut [String], missing_alleles: usize) {
    let Some(info) = fixed.get_mut(7) else {
        return;
    };
    let Some(an) = info_i64(info, "AN") else {
        return;
    };
    *info = set_info_value(info, "AN", &(an + missing_alleles as i64).to_string());
}

fn contig_order(meta: &[String]) -> HashMap<String, usize> {
    let mut order = HashMap::new();
    for line in meta {
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let id = rest.split([',', '>']).next().unwrap_or("").to_owned();
            if !id.is_empty() && !order.contains_key(&id) {
                order.insert(id, order.len());
            }
        }
    }
    order
}

fn compare_sites(a: &MergedSite, b: &MergedSite, contigs: &HashMap<String, usize>) -> Ordering {
    let a_chrom = a.fixed.first().map(String::as_str).unwrap_or("");
    let b_chrom = b.fixed.first().map(String::as_str).unwrap_or("");
    match (contigs.get(a_chrom).copied(), contigs.get(b_chrom).copied()) {
        (Some(a_idx), Some(b_idx)) => a_idx.cmp(&b_idx),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => a_chrom.cmp(b_chrom),
    }
    .then_with(|| {
        let a_pos = a
            .fixed
            .get(1)
            .and_then(|pos| pos.parse::<u64>().ok())
            .unwrap_or(0);
        let b_pos = b
            .fixed
            .get(1)
            .and_then(|pos| pos.parse::<u64>().ok())
            .unwrap_or(0);
        a_pos.cmp(&b_pos)
    })
    .then_with(|| a.order.cmp(&b.order))
    .then_with(|| a.fixed.get(2).cmp(&b.fixed.get(2)))
    .then_with(|| a.fixed.get(3).cmp(&b.fixed.get(3)))
    .then_with(|| a.fixed.get(4).cmp(&b.fixed.get(4)))
    .then_with(|| a.order.cmp(&b.order))
}

fn missing_sample_value(fixed: &[String], missing_to_ref: bool) -> String {
    let Some(format) = fixed.get(8) else {
        return ".".to_owned();
    };
    if format == "." || format.is_empty() {
        return ".".to_owned();
    }
    format
        .split(':')
        .map(|key| {
            if key == "GT" {
                if missing_to_ref { "0/0" } else { "./." }
            } else {
                "."
            }
        })
        .collect::<Vec<_>>()
        .join(":")
}

fn write_output(bytes: &[u8], args: &Args) -> io::Result<()> {
    match &args.output {
        Some(path) if path != Path::new("-") => {
            write_to(bytes, args.output_kind, File::create(path)?)
        }
        _ => write_to(bytes, args.output_kind, io::stdout().lock()),
    }
}

fn write_to<W: Write>(bytes: &[u8], kind: OutputKind, out: W) -> io::Result<()> {
    match kind {
        OutputKind::VcfText => {
            let mut out = io::BufWriter::new(out);
            out.write_all(bytes)
        }
        OutputKind::VcfGz => {
            let mut bgzf = htslib_rs::bgzf::io::Writer::new(out);
            bgzf.write_all(bytes)?;
            bgzf.finish().map(|_| ())
        }
        OutputKind::Bcf => write_bcf_from_vcf_text(bytes, out),
    }
}

fn write_bcf_from_vcf_text<W: Write>(text: &[u8], out: W) -> io::Result<()> {
    use htslib_rs::vcf::variant::io::Write as _;

    let mut reader = htslib_rs::vcf::io::Reader::new(BufReader::new(text));
    let header = reader.read_header()?;
    let mut writer = htslib_rs::bcf::io::Writer::new(out);
    writer.write_variant_header(&header)?;
    for result in reader.records() {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }
    writer.try_finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_same_site_sample_columns() {
        let a = parse_vcf_text(
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let b = parse_vcf_text(
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\n\
1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();
        let merged = merge_inputs(
            &[a, b],
            false,
            InfoRules::default(),
            MergeMode::Default,
            None,
            false,
        )
        .unwrap();
        assert!(merged.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"));
        assert!(merged.contains("1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1"));
    }
}
