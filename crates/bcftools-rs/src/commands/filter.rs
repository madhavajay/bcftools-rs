//! Port of `bcftools filter` (upstream `vcffilter.c`).
//!
//! MVP behavior
//! ------------
//!
//! - Reads VCF/VCF.gz/BCF and writes VCF/VCF.gz/BCF (`-o`/`-O v|z|b|u`).
//! - `-i EXPR` / `-e EXPR` evaluate against text VCF records using the shared
//!   [`crate::filter`] expression engine. Records that don't match (`-i`) or
//!   that match (`-e`) are either dropped (hard filter) or labeled (soft).
//! - `-s STRING` annotates the FILTER column of failed records with `STRING`
//!   instead of dropping them. Combined with `-m +`, the new tag is appended
//!   to existing FILTER values. With no mode flag the existing FILTER value
//!   is replaced.
//! - `-m x` resets PASS at sites that pass; `-m +` adds the soft-filter tag
//!   to existing FILTER values rather than replacing.
//! - `--no-version` suppresses the `##bcftools_filter{Version,Command}` header
//!   lines.
//! - `-r/-R/-t/-T` perform simple POS-based region/target restriction
//!   (no overlap-aware semantics).
//! - `--mask`/`-M` soft-filter records by inline/file regions. `^` negates the
//!   mask. `--mask-overlap 0|1|2` supports POS and REF-span matching.
//! - `-S`/`--set-GTs .|0` rewrites GT fields at failed sites to missing or
//!   reference, preserving ploidy separators and recalculating existing
//!   INFO/AC and INFO/AN tags. Simple FORMAT-scoped expressions rewrite only
//!   failed samples.
//! - `-g`/`--SnpGap` and `-G`/`--IndelGap` add local gap filters in text mode.
//!
//! Deferred: exact buffered gap-filter tie-breaking parity and full
//! filter-expression sample-vector semantics. These are tracked in TODO.md for
//! the filter engine.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    build_bcf_csi_with_min_shift, build_vcf_csi_from_path_with_min_shift, build_vcf_tbi_from_path,
    write_csi, write_tbi,
};

use crate::diagnostics::fmt_etag;
use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::header_version::{build_lines, command_time};
use crate::io::{VariantOutputFormat, apply_verbosity, init_index2, write_index_parse};

const USAGE: &str = "\n\
About: Apply fixed-threshold filters.\n\
Usage: bcftools filter [options] <in.vcf>\n\
\n\
Options:\n\
   -e, --exclude EXPR             Exclude sites for which the expression is true\n\
   -g, --SnpGap INT[:TYPE]        Filter SNPs within INT bp of indel/mnp/bnd/other/overlap [indel]\n\
   -G, --IndelGap INT             Filter clustered indels separated by INT or fewer bp\n\
   -i, --include EXPR             Include only sites for which the expression is true\n\
       --mask [^]REGION           Soft filter regions, \"^\" to negate\n\
   -M, --mask-file [^]FILE        Soft filter regions listed in a file, \"^\" to negate\n\
       --mask-overlap 0|1|2       Mask if POS in region (0), record overlaps (1), variant overlaps (2) [1]\n\
   -m, --mode [+x]                \"+\": add to existing FILTER; \"x\": reset filters at sites which pass\n\
       --no-version               Do not append version and command line to the header\n\
   -o, --output FILE              Write output to a file [standard output]\n\
   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF [v]\n\
   -r, --regions REGION           Restrict to comma-separated list of regions\n\
   -R, --regions-file FILE        Restrict to regions listed in a file\n\
   -s, --soft-filter STRING       Annotate FILTER column with <string>\n\
   -S, --set-GTs .|0              Set GTs at failed sites to missing (.) or ref (0)\n\
   -t, --targets REGION           Similar to -r but streams\n\
   -T, --targets-file FILE        Similar to -R but streams\n\
       --threads INT              Use multithreaded BGZF compression for compressed output\n\
   -v, --verbosity INT            Verbosity level\n\
   -W, --write-index[=FMT]        Automatically index the output files [off]\n\
\n";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

impl OutputKind {
    fn parse(raw: &str) -> Option<Self> {
        match raw.chars().next()? {
            'v' => Some(Self::VcfText),
            'z' => Some(Self::VcfGz),
            'b' | 'u' => Some(Self::Bcf),
            '0'..='9' => Some(Self::VcfGz),
            _ => None,
        }
    }
}

#[derive(Debug, Default, Clone)]
struct ModeFlags {
    additive: bool,
    reset_pass: bool,
}

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    include_expr: Option<String>,
    exclude_expr: Option<String>,
    soft_filter: Option<String>,
    set_gts: Option<SetGts>,
    mode: ModeFlags,
    regions: Vec<RegionSpec>,
    targets: Vec<RegionSpec>,
    mask: Vec<RegionSpec>,
    mask_negate: bool,
    mask_overlap: OverlapMode,
    snp_gap: Option<SnpGap>,
    indel_gap: Option<i64>,
    no_version: bool,
    write_index: Option<i32>,
    thread_count: Option<NonZero<usize>>,
}

#[derive(Debug, Clone)]
struct RegionSpec {
    contig: String,
    start: Option<i64>,
    end: Option<i64>,
}

#[derive(Debug, Clone)]
struct SnpGap {
    distance: i64,
    kinds: Vec<VariantKind>,
    label: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantKind {
    Snp,
    Indel,
    Mnp,
    Bnd,
    Other,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SetGts {
    Missing,
    Ref,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OverlapMode {
    Pos,
    Record,
    Variant,
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args, argv) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcffilter", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcffilter", &message));
            ExitCode::FAILURE
        }
    }
}

#[derive(Debug)]
enum ParseOutcome {
    Usage,
    Error(String),
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut input: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut output_kind = OutputKind::VcfText;
    let mut include_expr = None;
    let mut exclude_expr = None;
    let mut soft_filter = None;
    let mut set_gts = None;
    let mut mode = ModeFlags::default();
    let mut regions = Vec::new();
    let mut targets = Vec::new();
    let mut mask = Vec::new();
    let mut mask_negate = false;
    let mut mask_overlap = OverlapMode::Record;
    let mut snp_gap = None;
    let mut indel_gap = None;
    let mut no_version = false;
    let mut write_index = None;
    let mut thread_count = None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "--no-version" => no_version = true,
            "-i" | "--include" => include_expr = Some(next_string(&mut iter, "--include")?),
            "-e" | "--exclude" => exclude_expr = Some(next_string(&mut iter, "--exclude")?),
            "-g" | "--SnpGap" => {
                snp_gap = Some(parse_snp_gap(&next_string(&mut iter, "--SnpGap")?)?);
            }
            "-G" | "--IndelGap" => {
                indel_gap = Some(parse_indel_gap(&next_string(&mut iter, "--IndelGap")?)?);
            }
            "-s" | "--soft-filter" => {
                soft_filter = Some(next_string(&mut iter, "--soft-filter")?);
            }
            "-S" | "--set-GTs" => {
                set_gts = Some(parse_set_gts(&next_string(&mut iter, "--set-GTs")?)?);
            }
            "--mask" => {
                let raw_mask = next_string(&mut iter, "--mask")?;
                parse_mask_regions(&mut mask, &mut mask_negate, &raw_mask)?;
            }
            "-M" | "--mask-file" => {
                let raw_mask = next_string(&mut iter, "--mask-file")?;
                load_mask_file(&mut mask, &mut mask_negate, &raw_mask)?;
            }
            "--mask-overlap" => {
                mask_overlap = parse_overlap_mode(&next_string(&mut iter, "--mask-overlap")?)?;
            }
            "-m" | "--mode" => {
                set_mode(&mut mode, &next_string(&mut iter, "--mode")?)?;
            }
            "-W" | "--write-index" => {
                write_index = parse_write_index(None)?;
            }
            "-r" | "--regions" => {
                parse_region_list(&mut regions, &next_string(&mut iter, "--regions")?)?;
            }
            "-R" | "--regions-file" => {
                load_region_file(&mut regions, &next_string(&mut iter, "--regions-file")?)?;
            }
            "-t" | "--targets" => {
                parse_region_list(&mut targets, &next_string(&mut iter, "--targets")?)?;
            }
            "-T" | "--targets-file" => {
                load_region_file(&mut targets, &next_string(&mut iter, "--targets-file")?)?;
            }
            "--threads" => {
                thread_count = parse_threads(&next_string(&mut iter, "--threads")?)?;
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, "--output")?));
            }
            "-O" | "--output-type" => {
                let value = next_string(&mut iter, "--output-type")?;
                output_kind = parse_output_kind(&value)?;
            }
            "-v" | "--verbosity" => {
                let value = next_string(&mut iter, "--verbosity")?;
                if apply_verbosity(&value).is_err() {
                    return Err(ParseOutcome::Error(format!(
                        "Could not parse argument: --verbosity {value}"
                    )));
                }
            }
            _ if raw.starts_with("--include=") => {
                include_expr = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--exclude=") => {
                exclude_expr = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--SnpGap=") => {
                snp_gap = Some(parse_snp_gap(value_after_equals(&raw))?)
            }
            _ if raw.starts_with("--IndelGap=") => {
                indel_gap = Some(parse_indel_gap(value_after_equals(&raw))?)
            }
            _ if raw.starts_with("--soft-filter=") => {
                soft_filter = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--set-GTs=") => {
                set_gts = Some(parse_set_gts(value_after_equals(&raw))?)
            }
            _ if raw.starts_with("--mask=") => {
                parse_mask_regions(&mut mask, &mut mask_negate, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--mask-file=") => {
                load_mask_file(&mut mask, &mut mask_negate, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--mask-overlap=") => {
                mask_overlap = parse_overlap_mode(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--mode=") => set_mode(&mut mode, value_after_equals(&raw))?,
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--write-index=") => {
                write_index = parse_write_index(Some(value_after_equals(&raw)))?
            }
            _ if raw.starts_with("--regions=") => {
                parse_region_list(&mut regions, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--regions-file=") => {
                load_region_file(&mut regions, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--targets=") => {
                parse_region_list(&mut targets, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--targets-file=") => {
                load_region_file(&mut targets, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--threads=") => {
                thread_count = parse_threads(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-i") && raw.len() > 2 => include_expr = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-e") && raw.len() > 2 => exclude_expr = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-g") && raw.len() > 2 => {
                snp_gap = Some(parse_snp_gap(&raw[2..])?)
            }
            _ if raw.starts_with("-G") && raw.len() > 2 => {
                indel_gap = Some(parse_indel_gap(&raw[2..])?)
            }
            _ if raw.starts_with("-s") && raw.len() > 2 => soft_filter = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-S") && raw.len() > 2 => {
                set_gts = Some(parse_set_gts(&raw[2..])?)
            }
            _ if raw.starts_with("-m") && raw.len() > 2 => set_mode(&mut mode, &raw[2..])?,
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-W=") => write_index = parse_write_index(Some(&raw[3..]))?,
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!("Unrecognized option: {raw}")));
            }
            _ => {
                if input.is_some() {
                    return Err(ParseOutcome::Error(format!(
                        "Multiple positional inputs are not yet supported: {raw}"
                    )));
                }
                input = Some(PathBuf::from(arg));
            }
        }
    }

    let input = input.ok_or(ParseOutcome::Usage)?;
    if !mask.is_empty() && soft_filter.is_none() {
        return Err(ParseOutcome::Error(
            "The option --soft-filter is required with --mask and --mask-file options".into(),
        ));
    }

    Ok(Args {
        input,
        output,
        output_kind,
        include_expr,
        exclude_expr,
        soft_filter,
        set_gts,
        mode,
        regions,
        targets,
        mask,
        mask_negate,
        mask_overlap,
        snp_gap,
        indel_gap,
        no_version,
        write_index,
        thread_count,
    })
}

fn parse_snp_gap(raw: &str) -> Result<SnpGap, ParseOutcome> {
    let (distance, type_part) = raw.split_once(':').unwrap_or((raw, "indel"));
    let distance = distance
        .parse::<i64>()
        .map_err(|_| ParseOutcome::Error(format!("Could not parse argument: --SnpGap {raw}")))?;
    if distance < 0 {
        return Err(ParseOutcome::Error(format!(
            "Could not parse argument: --SnpGap {raw}"
        )));
    }
    let mut kinds = Vec::new();
    for token in type_part.split(',') {
        let kind = match token.to_ascii_lowercase().as_str() {
            "indel" => VariantKind::Indel,
            "mnp" => VariantKind::Mnp,
            "bnd" => VariantKind::Bnd,
            "other" => VariantKind::Other,
            "overlap" => VariantKind::Other,
            "" => continue,
            _ => {
                return Err(ParseOutcome::Error(format!(
                    "Could not parse \"{token}\" in \"--SnpGap {raw}\""
                )));
            }
        };
        kinds.push(kind);
    }
    if kinds.is_empty() {
        kinds.push(VariantKind::Indel);
    }
    Ok(SnpGap {
        distance,
        kinds,
        label: type_part.to_owned(),
    })
}

fn parse_indel_gap(raw: &str) -> Result<i64, ParseOutcome> {
    let distance = raw
        .parse::<i64>()
        .map_err(|_| ParseOutcome::Error(format!("Could not parse argument: --IndelGap {raw}")))?;
    if distance < 0 {
        return Err(ParseOutcome::Error(format!(
            "Could not parse argument: --IndelGap {raw}"
        )));
    }
    Ok(distance)
}

fn parse_set_gts(raw: &str) -> Result<SetGts, ParseOutcome> {
    match raw {
        "." => Ok(SetGts::Missing),
        "0" => Ok(SetGts::Ref),
        _ => Err(ParseOutcome::Error(format!(
            "The argument to -S not recognised: {raw}"
        ))),
    }
}

fn parse_overlap_mode(raw: &str) -> Result<OverlapMode, ParseOutcome> {
    match raw {
        "0" => Ok(OverlapMode::Pos),
        "1" => Ok(OverlapMode::Record),
        "2" => Ok(OverlapMode::Variant),
        _ => Err(ParseOutcome::Error(format!(
            "Could not parse: --mask-overlap {raw}"
        ))),
    }
}

fn set_mode(mode: &mut ModeFlags, raw: &str) -> Result<(), ParseOutcome> {
    for ch in raw.chars() {
        match ch {
            '+' => mode.additive = true,
            'x' | 'X' => mode.reset_pass = true,
            _ => {
                return Err(ParseOutcome::Error(format!(
                    "Unrecognized --mode flag '{ch}' (expected '+' or 'x')"
                )));
            }
        }
    }
    Ok(())
}

fn parse_mask_regions(
    out: &mut Vec<RegionSpec>,
    negate: &mut bool,
    raw: &str,
) -> Result<(), ParseOutcome> {
    let raw = parse_negated_mask(raw, negate);
    parse_region_list(out, raw)
}

fn load_mask_file(
    out: &mut Vec<RegionSpec>,
    negate: &mut bool,
    raw: &str,
) -> Result<(), ParseOutcome> {
    let raw = parse_negated_mask(raw, negate);
    load_region_file(out, raw)
}

fn parse_negated_mask<'a>(raw: &'a str, negate: &mut bool) -> &'a str {
    if let Some(rest) = raw.strip_prefix('^') {
        *negate = true;
        rest
    } else {
        raw
    }
}

fn parse_region_list(out: &mut Vec<RegionSpec>, raw: &str) -> Result<(), ParseOutcome> {
    for token in raw.split(',') {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        out.push(parse_region(token)?);
    }
    Ok(())
}

fn parse_region(token: &str) -> Result<RegionSpec, ParseOutcome> {
    let (contig, range) = match token.rfind(':') {
        Some(i) => (&token[..i], Some(&token[i + 1..])),
        None => (token, None),
    };
    let (start, end) = match range {
        None => (None, None),
        Some(r) => match r.split_once('-') {
            None => {
                let pos = parse_pos(r)?;
                (Some(pos), Some(pos))
            }
            Some((s, e)) => {
                let s = if s.is_empty() {
                    None
                } else {
                    Some(parse_pos(s)?)
                };
                let e = if e.is_empty() {
                    None
                } else {
                    Some(parse_pos(e)?)
                };
                (s, e)
            }
        },
    };
    Ok(RegionSpec {
        contig: contig.to_owned(),
        start,
        end,
    })
}

fn parse_pos(raw: &str) -> Result<i64, ParseOutcome> {
    raw.replace(',', "")
        .parse::<i64>()
        .map_err(|e| ParseOutcome::Error(format!("Could not parse region position '{raw}': {e}")))
}

fn load_region_file(out: &mut Vec<RegionSpec>, path: &str) -> Result<(), ParseOutcome> {
    let f = File::open(path)
        .map_err(|e| ParseOutcome::Error(format!("Could not read regions file {path}: {e}")))?;
    let is_bed = path.to_ascii_lowercase().ends_with(".bed")
        || path.to_ascii_lowercase().ends_with(".bed.gz");
    let reader: Box<dyn BufRead> = if path.ends_with(".gz") {
        Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(f)))
    } else {
        Box::new(BufReader::new(f))
    };
    for line in reader.lines() {
        let line = line.map_err(|e| ParseOutcome::Error(e.to_string()))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        if fields.len() >= 3 {
            let chrom = fields[0].to_owned();
            let mut start = parse_pos(fields[1])?;
            let end = parse_pos(fields[2])?;
            if is_bed {
                start += 1; // BED is 0-based half-open; convert to 1-based inclusive.
            }
            out.push(RegionSpec {
                contig: chrom,
                start: Some(start),
                end: Some(end),
            });
        } else {
            out.push(parse_region(trimmed)?);
        }
    }
    Ok(())
}

fn run(args: &Args, argv: &[OsString]) -> io::Result<()> {
    if args.write_index.is_some() && args.output.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W requires an output file",
        ));
    }

    let text = read_vcf_text(&args.input)?;
    let original_header = extract_header(&text);
    let body = &text[original_header.len()..];
    let header_text = ensure_filter_header(
        &original_header,
        args.soft_filter.as_deref(),
        args.include_expr.as_deref(),
        args.exclude_expr.as_deref(),
    );
    let header_text = ensure_gap_headers(&header_text, args);
    let header_text = if args.no_version {
        header_text
    } else {
        add_version_header(&header_text, argv)
    };

    let mut buffer = Vec::new();
    {
        let mut out: Box<dyn Write> = Box::new(&mut buffer);
        out.write_all(header_text.as_bytes())?;
        let mut records = collect_body_records(body);
        apply_gap_filters(&mut records, args);
        for trimmed_line in records {
            if trimmed_line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = trimmed_line.split('\t').collect();
            if fields.len() < 8 {
                continue;
            }
            if !record_in_regions(&fields, &args.regions, &args.targets) {
                continue;
            }
            let sample_passes = if args.set_gts.is_some() {
                sample_passes_for_set_gts(&fields, args)?
            } else {
                None
            };
            let mut pass = match (&sample_passes, args.soft_filter.as_deref()) {
                (Some(passes), Some(_)) => passes.iter().all(|passed| *passed),
                (Some(passes), None) => passes.iter().any(|passed| *passed),
                (None, _) => evaluate(&fields, args)?,
            };
            if !args.mask.is_empty() {
                let masked = record_overlaps_mask(&fields, &args.mask, args.mask_overlap);
                let mask_pass = if args.mask_negate { masked } else { !masked };
                pass &= mask_pass;
            }
            let gap_failed = has_gap_filter(&fields);
            if gap_failed {
                pass = false;
            }
            match (pass, args.soft_filter.as_deref()) {
                (true, _) => {
                    let mut rendered = if args.mode.reset_pass {
                        replace_filter(&fields, "PASS")
                    } else {
                        trimmed_line
                    };
                    if let (Some(target), Some(passes)) = (args.set_gts, sample_passes.as_deref())
                        && passes.iter().any(|passed| !*passed)
                    {
                        rendered = set_genotypes(&rendered, target, Some(passes));
                    }
                    out.write_all(rendered.as_bytes())?;
                    out.write_all(b"\n")?;
                }
                (false, Some(_)) if gap_failed => {
                    let mut rendered = trimmed_line;
                    if let Some(target) = args.set_gts {
                        rendered = set_genotypes(&rendered, target, sample_passes.as_deref());
                    }
                    out.write_all(rendered.as_bytes())?;
                    out.write_all(b"\n")?;
                }
                (false, Some(tag)) => {
                    let mut rendered = if args.mode.additive {
                        append_filter(&fields, tag)
                    } else {
                        replace_filter(&fields, tag)
                    };
                    if let Some(target) = args.set_gts {
                        rendered = set_genotypes(&rendered, target, sample_passes.as_deref());
                    }
                    out.write_all(rendered.as_bytes())?;
                    out.write_all(b"\n")?;
                }
                (false, None) if args.set_gts.is_some() => {
                    let mut rendered = trimmed_line;
                    if let Some(target) = args.set_gts {
                        rendered = set_genotypes(&rendered, target, sample_passes.as_deref());
                    }
                    out.write_all(rendered.as_bytes())?;
                    out.write_all(b"\n")?;
                }
                (false, None) => {
                    // Hard filter: drop the record.
                }
            }
        }
    }

    write_output(args, &buffer)?;
    if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
        write_index(path, args.output_kind, index_format)?;
    }
    Ok(())
}

fn collect_body_records(body: &str) -> Vec<String> {
    body.split_inclusive('\n')
        .filter(|line| !line.starts_with('#'))
        .map(|line| line.trim_end_matches('\n').trim_end_matches('\r'))
        .filter(|line| !line.is_empty())
        .map(ToOwned::to_owned)
        .collect()
}

fn ensure_gap_headers(header: &str, args: &Args) -> String {
    let mut out = header.to_owned();
    if let Some(spec) = &args.snp_gap {
        out = ensure_named_filter_header(
            &out,
            "SnpGap",
            &format!("SNP within {} bp of {}", spec.distance, spec.label),
        );
    }
    if let Some(distance) = args.indel_gap {
        out = ensure_named_filter_header(
            &out,
            "IndelGap",
            &format!("Indel within {distance} bp of an indel"),
        );
    }
    out
}

fn ensure_named_filter_header(header: &str, id: &str, description: &str) -> String {
    let needle = format!("##FILTER=<ID={id},");
    if header.lines().any(|line| line.starts_with(&needle)) {
        return header.to_owned();
    }
    let mut out = String::new();
    let mut inserted = false;
    for line in header.split_inclusive('\n') {
        if !inserted && line.starts_with("#CHROM\t") {
            out.push_str(&format!(
                "##FILTER=<ID={id},Description=\"{description}\">\n"
            ));
            inserted = true;
        }
        out.push_str(line);
    }
    out
}

fn evaluate(fields: &[&str], args: &Args) -> io::Result<bool> {
    if let Some(expr) = &args.include_expr {
        let value = evaluate_expression(expr, fields)?;
        if !value.truthy() {
            return Ok(false);
        }
    }
    if let Some(expr) = &args.exclude_expr {
        let value = evaluate_expression(expr, fields)?;
        if value.truthy() {
            return Ok(false);
        }
    }
    Ok(true)
}

fn evaluate_expression(expr: &str, fields: &[&str]) -> io::Result<FilterValue> {
    let context = record_context(fields);
    let fields_owned: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
    bcffilter::eval_expression_with(expr, &context, |name, sample_index| {
        if sample_index.is_some() {
            return None;
        }
        record_lookup(name, &fields_owned)
    })
}

fn record_context(fields: &[&str]) -> EvalContext {
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

fn sample_passes_for_set_gts(fields: &[&str], args: &Args) -> io::Result<Option<Vec<bool>>> {
    let expression = match (&args.include_expr, &args.exclude_expr) {
        (Some(expr), _) if expression_is_format_scoped(expr) => Some((expr.as_str(), true)),
        (_, Some(expr)) if expression_is_format_scoped(expr) => Some((expr.as_str(), false)),
        _ => None,
    };
    let Some((expr, is_include)) = expression else {
        return Ok(None);
    };
    if fields.len() <= 9 {
        return Ok(None);
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let wrapped_expr = format!("N_PASS({expr}) > 0");
    let fields_owned: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
    let mut passes = Vec::with_capacity(fields.len().saturating_sub(9));
    for sample in &fields[9..] {
        let context = EvalContext::new().with_sample(sample_values(&format_keys, sample));
        let matched =
            bcffilter::eval_expression_with(&wrapped_expr, &context, |name, sample_i| {
                if sample_i.is_some() {
                    return None;
                }
                record_lookup(name, &fields_owned)
            })?
            .truthy();
        passes.push(if is_include { matched } else { !matched });
    }
    Ok(Some(passes))
}

fn expression_is_format_scoped(expr: &str) -> bool {
    expr.contains("FMT/") || expr.contains("FORMAT/")
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

/// Look up a record-level identifier (CHROM/POS/ID/REF/ALT/QUAL/FILTER/INFO
/// tag) using the raw VCF text fields. Used by `filter` and other commands
/// that share the same text-based filter expression engine.
pub(crate) fn record_lookup(token: &str, fields: &[String]) -> Option<FilterValue> {
    let upper = token.to_ascii_uppercase();
    match upper.as_str() {
        "CHROM" => Some(FilterValue::String(fields[0].clone())),
        "POS" => fields[1].parse::<f64>().ok().map(FilterValue::Number),
        "ID" => Some(FilterValue::String(fields[2].clone())),
        "REF" => Some(FilterValue::String(fields[3].clone())),
        "ALT" => Some(FilterValue::List(
            fields[4]
                .split(',')
                .map(|s| FilterValue::String(s.to_owned()))
                .collect(),
        )),
        "QUAL" => match fields[5].as_str() {
            "." => Some(FilterValue::Missing),
            other => other
                .parse::<f64>()
                .ok()
                .map(FilterValue::Number)
                .or(Some(FilterValue::Missing)),
        },
        "FILTER" => Some(FilterValue::String(fields[6].clone())),
        "INDEL" => match info_value(&upper, fields.get(7)?.as_str()) {
            Some(FilterValue::Missing) => Some(FilterValue::Bool(false)),
            value => value,
        },
        _ => info_value(&upper, fields.get(7)?.as_str()),
    }
}

fn info_value(name: &str, info: &str) -> Option<FilterValue> {
    if info == "." {
        return Some(FilterValue::Missing);
    }
    for entry in info.split(';') {
        let (key, value) = match entry.split_once('=') {
            Some((k, v)) => (k, Some(v)),
            None => (entry, None),
        };
        if !key.eq_ignore_ascii_case(name) {
            continue;
        }
        return Some(match value {
            Some(v) => match v.parse::<f64>() {
                Ok(f) => FilterValue::Number(f),
                Err(_) => FilterValue::String(v.to_owned()),
            },
            None => FilterValue::Bool(true),
        });
    }
    Some(FilterValue::Missing)
}

fn record_in_regions(fields: &[&str], regions: &[RegionSpec], targets: &[RegionSpec]) -> bool {
    let chrom = fields[0];
    let pos = fields[1].parse::<i64>().unwrap_or(-1);
    let in_regions = regions.is_empty() || matches_any(regions, chrom, pos);
    let in_targets = targets.is_empty() || matches_any(targets, chrom, pos);
    in_regions && in_targets
}

fn record_overlaps_mask(fields: &[&str], mask: &[RegionSpec], mode: OverlapMode) -> bool {
    let chrom = fields[0];
    let pos = fields[1].parse::<i64>().unwrap_or(-1);
    let (beg, end) = match mode {
        OverlapMode::Pos => (pos, pos),
        OverlapMode::Record | OverlapMode::Variant => {
            let span_end = pos + fields[3].len().saturating_sub(1) as i64;
            (pos, span_end.max(pos))
        }
    };
    mask.iter().any(|spec| {
        if spec.contig != chrom {
            return false;
        }
        let spec_start = spec.start.unwrap_or(i64::MIN);
        let spec_end = spec.end.unwrap_or(i64::MAX);
        end >= spec_start && beg <= spec_end
    })
}

fn apply_gap_filters(records: &mut [String], args: &Args) {
    if args.snp_gap.is_none() && args.indel_gap.is_none() {
        return;
    }
    let metadata: Vec<RecordMeta> = records
        .iter()
        .map(|record| RecordMeta::from_record(record))
        .collect();
    let mut tags: Vec<Vec<&'static str>> = vec![Vec::new(); records.len()];
    if let Some(spec) = &args.snp_gap {
        mark_snp_gap(&metadata, spec, &mut tags);
    }
    if let Some(distance) = args.indel_gap {
        mark_indel_gap(&metadata, distance, &mut tags);
    }
    for (record, tags) in records.iter_mut().zip(tags) {
        for tag in tags {
            *record = append_filter_to_record(record, tag);
        }
    }
}

#[derive(Debug, Clone)]
struct RecordMeta {
    chrom: String,
    pos: i64,
    end: i64,
    qual: Option<f64>,
    first_alt_ac: Option<i64>,
    kind: VariantKind,
}

impl RecordMeta {
    fn from_record(record: &str) -> Self {
        let fields: Vec<&str> = record.split('\t').collect();
        let chrom = fields.first().copied().unwrap_or("").to_owned();
        let pos = fields
            .get(1)
            .and_then(|value| value.parse::<i64>().ok())
            .unwrap_or(-1);
        let reference = fields.get(3).copied().unwrap_or("");
        let alts: Vec<&str> = fields
            .get(4)
            .copied()
            .unwrap_or(".")
            .split(',')
            .filter(|alt| *alt != ".")
            .collect();
        let end = pos + reference.len().saturating_sub(1) as i64;
        let qual = fields.get(5).and_then(|value| value.parse::<f64>().ok());
        let first_alt_ac = first_alt_allele_count(&fields);
        let kind = classify_variant(reference, &alts);
        Self {
            chrom,
            pos,
            end: end.max(pos),
            qual,
            first_alt_ac,
            kind,
        }
    }
}

fn mark_snp_gap(metadata: &[RecordMeta], spec: &SnpGap, tags: &mut [Vec<&'static str>]) {
    for (i, rec) in metadata.iter().enumerate() {
        if rec.kind != VariantKind::Snp {
            continue;
        }
        let near_gap_variant = metadata.iter().enumerate().any(|(j, other)| {
            i != j
                && rec.chrom == other.chrom
                && spec.kinds.contains(&other.kind)
                && intervals_within(rec.pos, rec.end, other.pos, other.end, spec.distance)
        });
        if near_gap_variant && !tags[i].contains(&"SnpGap") {
            tags[i].push("SnpGap");
        }
    }
}

fn mark_indel_gap(metadata: &[RecordMeta], distance: i64, tags: &mut [Vec<&'static str>]) {
    let mut i = 0;
    while i < metadata.len() {
        if metadata[i].kind != VariantKind::Indel {
            i += 1;
            continue;
        }
        let chrom = metadata[i].chrom.clone();
        let mut cluster = vec![i];
        let mut cluster_end = metadata[i].end + distance;
        let mut j = i + 1;
        while j < metadata.len() {
            let rec = &metadata[j];
            if rec.chrom != chrom || rec.pos > cluster_end {
                break;
            }
            if rec.kind == VariantKind::Indel {
                cluster.push(j);
                cluster_end = cluster_end.max(rec.end + distance);
            }
            j += 1;
        }
        if cluster.len() > 1 {
            let best = best_indel_gap_record(&cluster, metadata);
            for idx in cluster.iter().copied().filter(|idx| *idx != best) {
                if !tags[idx].contains(&"IndelGap") {
                    tags[idx].push("IndelGap");
                }
            }
        }
        i = j.max(i + 1);
    }
}

fn best_indel_gap_record(cluster: &[usize], metadata: &[RecordMeta]) -> usize {
    let mut best_qual_idx = None;
    let mut best_qual = f64::NEG_INFINITY;
    for &idx in cluster {
        let Some(qual) = metadata[idx].qual.filter(|qual| *qual > 0.0) else {
            continue;
        };
        if best_qual_idx.is_none() || qual > best_qual {
            best_qual_idx = Some(idx);
            best_qual = qual;
        }
    }
    if let Some(idx) = best_qual_idx {
        return idx;
    }

    let mut best_ac_idx = None;
    let mut best_ac = i64::MIN;
    for &idx in cluster {
        let Some(ac) = metadata[idx].first_alt_ac else {
            continue;
        };
        if best_ac_idx.is_none() || ac > best_ac {
            best_ac_idx = Some(idx);
            best_ac = ac;
        }
    }
    best_ac_idx.unwrap_or(cluster[0])
}

fn intervals_within(a_start: i64, a_end: i64, b_start: i64, b_end: i64, distance: i64) -> bool {
    if a_end < b_start {
        b_start - a_end <= distance
    } else if b_end < a_start {
        a_start - b_end <= distance
    } else {
        true
    }
}

fn classify_variant(reference: &str, alts: &[&str]) -> VariantKind {
    let mut kind = None;
    for alt in alts {
        let alt_kind = if alt.starts_with('<') || alt.contains('[') || alt.contains(']') {
            VariantKind::Bnd
        } else if reference.len() == 1 && alt.len() == 1 {
            VariantKind::Snp
        } else if reference.len() == alt.len() {
            VariantKind::Mnp
        } else if reference.len() != alt.len() {
            VariantKind::Indel
        } else {
            VariantKind::Other
        };
        kind = match kind {
            None => Some(alt_kind),
            Some(prev) if prev == alt_kind => Some(prev),
            Some(_) => Some(VariantKind::Other),
        };
    }
    kind.unwrap_or(VariantKind::Other)
}

fn first_alt_allele_count(fields: &[&str]) -> Option<i64> {
    if fields.len() <= 9 {
        return None;
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let gt_index = format_keys.iter().position(|key| *key == "GT")?;
    let mut count = 0_i64;
    for sample in &fields[9..] {
        let Some(gt) = sample.split(':').nth(gt_index) else {
            continue;
        };
        for allele in gt.split(['/', '|']) {
            if allele == "1" {
                count += 1;
            }
        }
    }
    Some(count)
}

fn append_filter_to_record(record: &str, tag: &str) -> String {
    let fields: Vec<&str> = record.split('\t').collect();
    append_filter(&fields, tag)
}

fn has_gap_filter(fields: &[&str]) -> bool {
    fields.get(6).is_some_and(|filter| {
        filter
            .split(';')
            .any(|tag| tag == "SnpGap" || tag == "IndelGap")
    })
}

fn matches_any(specs: &[RegionSpec], chrom: &str, pos: i64) -> bool {
    specs.iter().any(|spec| {
        spec.contig == chrom
            && spec.start.map(|s| pos >= s).unwrap_or(true)
            && spec.end.map(|e| pos <= e).unwrap_or(true)
    })
}

fn replace_filter(fields: &[&str], new_value: &str) -> String {
    let mut out = String::new();
    for (i, field) in fields.iter().enumerate() {
        if i > 0 {
            out.push('\t');
        }
        if i == 6 {
            out.push_str(new_value);
        } else {
            out.push_str(field);
        }
    }
    out
}

fn append_filter(fields: &[&str], tag: &str) -> String {
    let existing = fields[6];
    let new_value = if existing == "." || existing == "PASS" || existing.is_empty() {
        tag.to_owned()
    } else if existing.split(';').any(|t| t == tag) {
        existing.to_owned()
    } else {
        format!("{existing};{tag}")
    };
    replace_filter(fields, &new_value)
}

fn set_genotypes(record: &str, target: SetGts, sample_passes: Option<&[bool]>) -> String {
    let mut fields: Vec<String> = record.split('\t').map(ToOwned::to_owned).collect();
    if fields.len() <= 9 {
        return record.to_owned();
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let Some(gt_index) = format_keys.iter().position(|key| *key == "GT") else {
        return record.to_owned();
    };
    for (sample_i, sample) in fields[9..].iter_mut().enumerate() {
        if sample_passes
            .and_then(|passes| passes.get(sample_i))
            .copied()
            .unwrap_or(false)
        {
            continue;
        }
        let mut parts: Vec<String> = sample.split(':').map(ToOwned::to_owned).collect();
        if let Some(gt) = parts.get_mut(gt_index) {
            *gt = rewrite_gt(gt, target);
            *sample = parts.join(":");
        }
    }
    recalculate_ac_an(&mut fields, gt_index);
    fields.join("\t")
}

fn recalculate_ac_an(fields: &mut [String], gt_index: usize) {
    if fields.len() <= 9 {
        return;
    }
    let alt_count = fields[4]
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .count();
    let mut ac = vec![0usize; alt_count];
    let mut an = 0usize;
    for sample in &fields[9..] {
        let Some(gt) = sample.split(':').nth(gt_index) else {
            continue;
        };
        for allele in gt.split(['/', '|']) {
            if allele == "." || allele.is_empty() {
                continue;
            }
            let Ok(index) = allele.parse::<usize>() else {
                continue;
            };
            an += 1;
            if index > 0
                && let Some(count) = ac.get_mut(index - 1)
            {
                *count += 1;
            }
        }
    }
    update_info_ac_an(&mut fields[7], &ac, an);
}

fn update_info_ac_an(info: &mut String, ac: &[usize], an: usize) {
    let mut saw_ac = false;
    let mut saw_an = false;
    let new_ac = ac
        .iter()
        .map(|count| count.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let new_an = an.to_string();
    let mut updated = Vec::new();

    if info != "." && !info.is_empty() {
        for item in info.split(';') {
            if item == "AC" || item.starts_with("AC=") {
                saw_ac = true;
                updated.push(format!("AC={new_ac}"));
            } else if item == "AN" || item.starts_with("AN=") {
                saw_an = true;
                updated.push(format!("AN={new_an}"));
            } else {
                updated.push(item.to_owned());
            }
        }
    }

    if saw_ac || saw_an {
        if updated.is_empty() {
            *info = ".".to_owned();
        } else {
            *info = updated.join(";");
        }
    }
}

fn rewrite_gt(gt: &str, target: SetGts) -> String {
    let replacement = match target {
        SetGts::Missing => ".",
        SetGts::Ref => "0",
    };
    let mut out = String::new();
    let mut allele = String::new();
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            if !allele.is_empty() {
                out.push_str(replacement);
                allele.clear();
            }
            out.push(ch);
        } else {
            allele.push(ch);
        }
    }
    if !allele.is_empty() {
        out.push_str(replacement);
    }
    if out.is_empty() {
        replacement.to_owned()
    } else {
        out
    }
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }
    let mut text = String::new();
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let mut dec = flate2::read::MultiGzDecoder::new(f);
        dec.read_to_string(&mut text)?;
    } else {
        text = fs::read_to_string(path)?;
    }
    crate::vcf_compat::normalize_vcf_text(&mut text);
    Ok(text)
}

fn extract_header(text: &str) -> String {
    let mut out = String::new();
    for line in text.split_inclusive('\n') {
        if !line.starts_with('#') {
            break;
        }
        out.push_str(line);
        if line.starts_with("#CHROM\t") {
            break;
        }
    }
    out
}

fn ensure_filter_header(
    header: &str,
    tag: Option<&str>,
    include_expr: Option<&str>,
    exclude_expr: Option<&str>,
) -> String {
    let Some(tag) = tag else {
        return header.to_owned();
    };
    let needle = format!("##FILTER=<ID={tag},");
    if header.lines().any(|l| l.starts_with(&needle)) {
        return header.to_owned();
    }
    let description = if let Some(expr) = exclude_expr {
        format!("Set if true: {}", escape_header_description(expr))
    } else if let Some(expr) = include_expr {
        format!("Set if not true: {}", escape_header_description(expr))
    } else {
        "Set if not true: filter expression".to_owned()
    };
    let mut out = String::new();
    let mut inserted = false;
    let mut pass_inserted = header
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    for line in header.split_inclusive('\n') {
        out.push_str(line);
        if !pass_inserted && line.starts_with("##fileformat=") {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            pass_inserted = true;
        }
        if !inserted && line.starts_with("#CHROM\t") {
            let chrom_line = out.split_off(out.len() - line.len());
            out.push_str(&format!(
                "##FILTER=<ID={tag},Description=\"{description}\">\n"
            ));
            out.push_str(&chrom_line);
            inserted = true;
        }
    }
    out
}

fn escape_header_description(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}

fn add_version_header(header: &str, argv: &[OsString]) -> String {
    let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
    prog_argv.extend(argv.iter().cloned());
    let lines = build_lines("bcftools_filter", &prog_argv, command_time());
    let mut out = String::new();
    let mut inserted = false;
    for line in header.split_inclusive('\n') {
        if !inserted && line.starts_with("#CHROM\t") {
            out.push_str(&lines.version_line);
            out.push('\n');
            out.push_str(&lines.command_line);
            out.push('\n');
            inserted = true;
        }
        out.push_str(line);
    }
    out
}

fn write_output(args: &Args, buffer: &[u8]) -> io::Result<()> {
    match args.output_kind {
        OutputKind::VcfText => match &args.output {
            Some(path) => fs::write(path, buffer),
            None => io::stdout().lock().write_all(buffer),
        },
        OutputKind::VcfGz => {
            if let (Some(path), Some(thread_count)) = (&args.output, args.thread_count) {
                let file = File::create(path)?;
                let mut bgzf =
                    htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
                bgzf.write_all(buffer)?;
                let _file = bgzf.finish()?;
                return Ok(());
            }
            let inner: Box<dyn Write> = match &args.output {
                Some(path) => Box::new(File::create(path)?),
                None => Box::new(io::stdout().lock()),
            };
            let mut bgzf = htslib_rs::bgzf::io::Writer::new(inner);
            bgzf.write_all(buffer)?;
            bgzf.finish()?;
            Ok(())
        }
        OutputKind::Bcf => {
            // Re-parse the buffered VCF text into BCF.
            use htslib_rs::vcf::variant::io::Write as _;
            let mut reader = htslib_rs::vcf::io::Reader::new(BufReader::new(buffer));
            let header = reader.read_header()?;
            if let (Some(path), Some(thread_count)) = (&args.output, args.thread_count) {
                let file = File::create(path)?;
                let bgzf =
                    htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
                let mut writer = htslib_rs::bcf::io::Writer::from(bgzf);
                writer.write_variant_header(&header)?;
                for result in reader.records() {
                    let record = result?;
                    writer.write_variant_record(&header, &record)?;
                }
                let mut bgzf = writer.into_inner();
                let _file = bgzf.finish()?;
                return Ok(());
            }
            let dst: Box<dyn Write> = match &args.output {
                Some(path) => Box::new(File::create(path)?),
                None => Box::new(io::stdout().lock()),
            };
            let mut writer = htslib_rs::bcf::io::Writer::new(dst);
            writer.write_variant_header(&header)?;
            for result in reader.records() {
                let record = result?;
                writer.write_variant_record(&header, &record)?;
            }
            writer.try_finish()
        }
    }
}

fn parse_output_kind(raw: &str) -> Result<OutputKind, ParseOutcome> {
    OutputKind::parse(raw)
        .ok_or_else(|| ParseOutcome::Error(format!("The output type \"{raw}\" not recognised")))
}

fn parse_write_index(raw: Option<&str>) -> Result<Option<i32>, ParseOutcome> {
    write_index_parse(raw).map(Some).ok_or_else(|| {
        ParseOutcome::Error(format!(
            "The index format \"{}\" not recognised",
            raw.unwrap_or("")
        ))
    })
}

fn parse_threads(raw: &str) -> Result<Option<NonZero<usize>>, ParseOutcome> {
    raw.parse::<usize>()
        .map(NonZero::new)
        .map_err(|_| ParseOutcome::Error(format!("Could not parse argument: --threads {raw}")))
}

fn write_index(path: &Path, output_kind: OutputKind, index_format: i32) -> io::Result<()> {
    let path_str = path
        .to_str()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "non-UTF-8 output path"))?;
    let output_format = match output_kind {
        OutputKind::VcfText | OutputKind::VcfGz => VariantOutputFormat::Vcf,
        OutputKind::Bcf => VariantOutputFormat::Bcf,
    };
    let Some(plan) =
        init_index2(Some(path_str), index_format, output_format).map_err(io::Error::other)?
    else {
        return Ok(());
    };

    if output_format == VariantOutputFormat::Bcf {
        if plan.min_shift == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BCF requires CSI (min_shift > 0)",
            ));
        }
        let index = build_bcf_csi_with_min_shift(path, plan.min_shift)?;
        write_csi(plan.index_path, &index)
    } else if plan.min_shift == 0 {
        let index = build_vcf_tbi_from_path(path)?;
        write_tbi(plan.index_path, &index)
    } else {
        let index = build_vcf_csi_from_path_with_min_shift(path, plan.min_shift)?;
        write_csi(plan.index_path, &index)
    }
}

fn next_string<'a, I>(iter: &mut I, name: &str) -> Result<String, ParseOutcome>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .map(|v| v.to_string_lossy().into_owned())
        .ok_or_else(|| ParseOutcome::Error(format!("missing argument for {name}")))
}

fn value_after_equals(raw: &str) -> &str {
    raw.split_once('=').map(|(_, v)| v).unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_region_handles_chrom_only() {
        let r = parse_region("chr1").unwrap();
        assert_eq!(r.contig, "chr1");
        assert_eq!(r.start, None);
        assert_eq!(r.end, None);
    }

    #[test]
    fn parse_region_handles_chrom_pos() {
        let r = parse_region("chr1:100").unwrap();
        assert_eq!(r.contig, "chr1");
        assert_eq!(r.start, Some(100));
        assert_eq!(r.end, Some(100));
    }

    #[test]
    fn parse_region_handles_chrom_range() {
        let r = parse_region("chr1:100-200").unwrap();
        assert_eq!(r.contig, "chr1");
        assert_eq!(r.start, Some(100));
        assert_eq!(r.end, Some(200));
    }

    #[test]
    fn append_filter_dedupes_existing_tag() {
        let fields = &["1", "100", ".", "A", "C", "100", "LowQual", "."];
        assert_eq!(
            append_filter(fields, "LowQual"),
            "1\t100\t.\tA\tC\t100\tLowQual\t."
        );
        assert_eq!(
            append_filter(fields, "Other"),
            "1\t100\t.\tA\tC\t100\tLowQual;Other\t."
        );
    }

    #[test]
    fn replace_filter_overwrites_filter_column() {
        let fields = &["1", "100", ".", "A", "C", "100", "LowQual", "."];
        assert_eq!(
            replace_filter(fields, "PASS"),
            "1\t100\t.\tA\tC\t100\tPASS\t."
        );
    }

    #[test]
    fn truthy_handles_basic_kinds() {
        assert!(FilterValue::Bool(true).truthy());
        assert!(!FilterValue::Bool(false).truthy());
        assert!(FilterValue::Number(5.0).truthy());
        assert!(!FilterValue::Number(0.0).truthy());
        assert!(!FilterValue::Missing.truthy());
    }

    #[test]
    fn record_lookup_pulls_basic_fields() {
        let fields =
            ["1", "100", "rsX", "A", "C", "37.5", "PASS", "DP=12;AF=0.25"].map(String::from);
        let v = record_lookup("DP", &fields);
        match v {
            Some(FilterValue::Number(n)) => assert_eq!(n, 12.0),
            other => panic!("DP did not resolve to Number: {other:?}"),
        }
        let v = record_lookup("AF", &fields);
        match v {
            Some(FilterValue::Number(n)) => assert_eq!(n, 0.25),
            other => panic!("AF did not resolve to Number: {other:?}"),
        }
    }
}
