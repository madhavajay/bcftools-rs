//! Port of `bcftools stats` (upstream `vcfstats.c`).
//!
//! MVP behavior
//! ------------
//!
//! - Single-input statistics and text-backed two-input set statistics. Full
//!   indexed synced-reader parity remains deferred in TODO.md.
//! - Emits the canonical `# SN` (summary numbers) and `# TSTV`
//!   (transitions/transversions) sections, plus a `# AF` allele-frequency
//!   distribution and a `# ST` substitution-types breakdown. Output text is
//!   structured to match upstream's tab-delimited shape so downstream
//!   `plot-vcfstats` can still parse it.
//! - Honors `-f`/`--apply-filters`, `-r`/`-R`/`-t`/`-T` (POS-based),
//!   `-i`/`-e` expression filters via the shared filter engine,
//!   `--af-bins`, `--af-tag`, `-u`/`--user-tstv`, `-d`/`--depth`,
//!   `-E`/`--exons` for local indel frame-shift classification,
//!   `-F`/`--fasta-ref` for basic indel-context sections,
//!   `-1`/`--1st-allele-only`, `-I`/`--split-by-ID`, `-s`/`-S` sample
//!   selection with core PSC, PSI, SiS, and VAF sections.
//! - Reads VCF/VCF.gz/BCF inputs through the existing readers and the
//!   Kestrel-tolerant `vcf_compat::NormalizeFileformat` adapter for plain
//!   text.
//!
//! Deferred (tracked in TODO.md): full indexed two-input synced-reader
//! parity and exact PSC/PSI/indel-context parity.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{self, Compression, Exact};

use crate::diagnostics::fmt_etag;
use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::io::apply_verbosity;
use crate::reference::FastaReference;
use crate::synced::CollapseMode;

const INDEL_CONTEXT_REPEAT_LEN: usize = 10;
const SUBSTITUTION_TYPES: [&str; 12] = [
    "A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G",
];

const USAGE: &str = "\n\
About: Parses VCF or BCF and produces stats.\n\
Usage: bcftools stats [options] <A.vcf>\n\
\n\
Options:\n\
       --af-bins LIST           Allele frequency bins (default 0.1,0.5,1)\n\
       --af-tag TAG             INFO tag to use for allele frequency bins [AF]\n\
   -1, --1st-allele-only        Include only 1st allele at multiallelic sites\n\
   -c, --collapse STRING        Treat compatible records as identical [none]\n\
   -d, --depth INT,INT,INT      Depth distribution: min,max,bin size [0,500,1]\n\
   -e, --exclude EXPR           Exclude sites for which the expression is true\n\
   -E, --exons FILE             Tab-delimited exons file for indel frameshifts\n\
   -f, --apply-filters LIST     Require at least one of the listed FILTER strings\n\
   -F, --fasta-ref FILE         Reference FASTA for indel context stats\n\
   -i, --include EXPR           Select sites for which the expression is true\n\
   -I, --split-by-ID            Collect stats for sites with ID separately\n\
   -r, --regions REGION         Restrict to comma-separated list of regions\n\
   -R, --regions-file FILE      Restrict to regions listed in a file\n\
   -s, --samples LIST           Comma-separated samples, ^ to exclude, - for all\n\
   -S, --samples-file FILE      File of samples, ^FILE to exclude\n\
   -t, --targets REGION         Similar to -r but streams\n\
   -T, --targets-file FILE      Similar to -R but streams\n\
   -u, --user-tstv TAG[:min:max:n] Collect Ts/Tv stats for an INFO tag [0:1:100]\n\
   -v, --verbosity INT          Verbosity level\n\
\n";

#[derive(Debug)]
struct Args {
    inputs: Vec<PathBuf>,
    include_expr: Option<String>,
    exclude_expr: Option<String>,
    apply_filters: Option<Vec<String>>,
    af_bins: Vec<f64>,
    af_tag: String,
    depth: DepthSpec,
    first_allele_only: bool,
    split_by_id: bool,
    regions: Vec<RegionSpec>,
    targets: Vec<RegionSpec>,
    exons: Vec<RegionSpec>,
    fasta_ref: Option<PathBuf>,
    user_tstv: Vec<UserTstvSpec>,
    sample_list: Option<String>,
    samples_is_file: bool,
    collapse: CollapseMode,
}

#[derive(Debug, Clone)]
struct RegionSpec {
    contig: String,
    start: Option<i64>,
    end: Option<i64>,
}

#[derive(Debug, Clone)]
struct UserTstvSpec {
    tag: String,
    index: usize,
    min: f64,
    max: f64,
    nbins: usize,
}

#[derive(Debug, Clone)]
struct SelectedSamples {
    names: Vec<String>,
    indices: Vec<usize>,
}

#[derive(Debug, Clone)]
struct DepthSpec {
    min: i64,
    max: i64,
    step: i64,
    nbins: usize,
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args, argv) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfstats", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfstats", &message));
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
    let mut inputs: Vec<PathBuf> = Vec::new();
    let mut include_expr = None;
    let mut exclude_expr = None;
    let mut apply_filters: Option<Vec<String>> = None;
    let mut af_bins: Vec<f64> = vec![0.1, 0.5, 1.0];
    let mut af_tag = "AF".to_owned();
    let mut depth = DepthSpec::default();
    let mut first_allele_only = false;
    let mut split_by_id = false;
    let mut regions = Vec::new();
    let mut targets = Vec::new();
    let mut exons = Vec::new();
    let mut fasta_ref = None;
    let mut user_tstv = Vec::new();
    let mut sample_list = None;
    let mut samples_is_file = false;
    let mut collapse = CollapseMode::None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-1" | "--1st-allele-only" => first_allele_only = true,
            "-c" | "--collapse" => {
                collapse = parse_collapse(&next_string(&mut iter, raw.as_ref())?)?
            }
            "-I" | "--split-by-ID" => split_by_id = true,
            "-i" | "--include" => include_expr = Some(next_string(&mut iter, "--include")?),
            "-e" | "--exclude" => exclude_expr = Some(next_string(&mut iter, "--exclude")?),
            "-E" | "--exons" => load_region_file(&mut exons, &next_string(&mut iter, "--exons")?)?,
            "-F" | "--fasta-ref" => {
                fasta_ref = Some(PathBuf::from(next_string(&mut iter, "--fasta-ref")?))
            }
            "-d" | "--depth" => {
                depth = parse_depth(&next_string(&mut iter, "--depth")?)?;
            }
            "-f" | "--apply-filters" => {
                apply_filters = Some(parse_filter_list(&next_string(&mut iter, raw.as_ref())?));
            }
            "--af-bins" => {
                af_bins = parse_af_bins(&next_string(&mut iter, "--af-bins")?)?;
            }
            "--af-tag" => {
                af_tag = next_string(&mut iter, "--af-tag")?;
            }
            "-r" | "--regions" => {
                parse_region_list(&mut regions, &next_string(&mut iter, "--regions")?)?;
            }
            "-R" | "--regions-file" => {
                load_region_file(&mut regions, &next_string(&mut iter, "--regions-file")?)?;
            }
            "-s" | "--samples" => {
                sample_list = Some(next_string(&mut iter, "--samples")?);
                samples_is_file = false;
            }
            "-S" | "--samples-file" => {
                sample_list = Some(next_string(&mut iter, "--samples-file")?);
                samples_is_file = true;
            }
            "-t" | "--targets" => {
                parse_region_list(&mut targets, &next_string(&mut iter, "--targets")?)?;
            }
            "-T" | "--targets-file" => {
                load_region_file(&mut targets, &next_string(&mut iter, "--targets-file")?)?;
            }
            "-u" | "--user-tstv" => {
                user_tstv.push(parse_user_tstv(&next_string(&mut iter, "--user-tstv")?)?);
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
            _ if raw.starts_with("--exons=") => {
                load_region_file(&mut exons, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--fasta-ref=") => {
                fasta_ref = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--depth=") => depth = parse_depth(value_after_equals(&raw))?,
            _ if raw.starts_with("--collapse=") => {
                collapse = parse_collapse(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--apply-filters=") => {
                apply_filters = Some(parse_filter_list(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--af-bins=") => {
                af_bins = parse_af_bins(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--af-tag=") => {
                af_tag = value_after_equals(&raw).to_owned();
            }
            _ if raw.starts_with("--regions=") => {
                parse_region_list(&mut regions, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--regions-file=") => {
                load_region_file(&mut regions, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--samples=") => {
                sample_list = Some(value_after_equals(&raw).to_owned());
                samples_is_file = false;
            }
            _ if raw.starts_with("--samples-file=") => {
                sample_list = Some(value_after_equals(&raw).to_owned());
                samples_is_file = true;
            }
            _ if raw.starts_with("--targets=") => {
                parse_region_list(&mut targets, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--targets-file=") => {
                load_region_file(&mut targets, value_after_equals(&raw))?
            }
            _ if raw.starts_with("--user-tstv=") => {
                user_tstv.push(parse_user_tstv(value_after_equals(&raw))?);
            }
            _ if raw.starts_with("-c") && raw.len() > 2 => {
                collapse = parse_collapse(&raw[2..])?;
            }
            _ if raw.starts_with("-E") && raw.len() > 2 => {
                load_region_file(&mut exons, &raw[2..])?;
            }
            _ if raw.starts_with("-F") && raw.len() > 2 => {
                fasta_ref = Some(PathBuf::from(&raw[2..]));
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!("Unrecognized option: {raw}")));
            }
            _ => inputs.push(PathBuf::from(arg)),
        }
    }

    if inputs.is_empty() {
        return Err(ParseOutcome::Usage);
    }
    if inputs.len() > 2 {
        return Err(ParseOutcome::Usage);
    }

    Ok(Args {
        inputs,
        include_expr,
        exclude_expr,
        apply_filters,
        af_bins,
        af_tag,
        depth,
        first_allele_only,
        split_by_id,
        regions,
        targets,
        exons,
        fasta_ref,
        user_tstv,
        sample_list,
        samples_is_file,
        collapse,
    })
}

fn parse_collapse(raw: &str) -> Result<CollapseMode, ParseOutcome> {
    match raw {
        "none" | "snps" | "indels" | "both" | "any" | "all" | "some" => {
            raw.parse::<CollapseMode>().map_err(|_| {
                ParseOutcome::Error(format!("The --collapse string \"{raw}\" not recognised."))
            })
        }
        _ => Err(ParseOutcome::Error(format!(
            "The --collapse string \"{raw}\" not recognised."
        ))),
    }
}

fn parse_filter_list(raw: &str) -> Vec<String> {
    raw.split(',').map(|s| s.trim().to_owned()).collect()
}

fn parse_af_bins(raw: &str) -> Result<Vec<f64>, ParseOutcome> {
    raw.split(',')
        .map(|s| {
            s.trim().parse::<f64>().map_err(|e| {
                ParseOutcome::Error(format!("Could not parse --af-bins value '{s}': {e}"))
            })
        })
        .collect()
}

impl Default for DepthSpec {
    fn default() -> Self {
        Self {
            min: 0,
            max: 500,
            step: 1,
            nbins: 4 + (500usize),
        }
    }
}

fn parse_depth(raw: &str) -> Result<DepthSpec, ParseOutcome> {
    let fields: Vec<&str> = raw.split(',').collect();
    if fields.len() != 3 {
        return Err(ParseOutcome::Error(format!(
            "Could not parse --depth {raw}"
        )));
    }
    let min = fields[0]
        .parse::<i64>()
        .map_err(|_| ParseOutcome::Error(format!("Could not parse --depth {raw}")))?;
    let max = fields[1]
        .parse::<i64>()
        .map_err(|_| ParseOutcome::Error(format!("Could not parse --depth {raw}")))?;
    let step = fields[2]
        .parse::<i64>()
        .map_err(|_| ParseOutcome::Error(format!("Could not parse --depth {raw}")))?;
    if min < 0 || min >= max || step <= 0 || step > max - min + 1 {
        return Err(ParseOutcome::Error(format!(
            "Is this a typo? --depth {raw}"
        )));
    }
    Ok(DepthSpec {
        min,
        max,
        step,
        nbins: 4 + ((max - min) / step) as usize,
    })
}

fn parse_user_tstv(raw: &str) -> Result<UserTstvSpec, ParseOutcome> {
    let mut parts = raw.split(':');
    let tag_part = parts.next().unwrap_or("");
    if tag_part.is_empty() {
        return Err(ParseOutcome::Error("empty --user-tstv tag".into()));
    }
    let (tag, index) = parse_tag_index(tag_part)?;
    let min = match parts.next() {
        Some(v) if !v.is_empty() => parse_f64(v, "--user-tstv min")?,
        _ => 0.0,
    };
    let max = match parts.next() {
        Some(v) if !v.is_empty() => parse_f64(v, "--user-tstv max")?,
        _ => 1.0,
    };
    let nbins = match parts.next() {
        Some(v) if !v.is_empty() => v.parse::<usize>().map_err(|e| {
            ParseOutcome::Error(format!("Could not parse --user-tstv bin count '{v}': {e}"))
        })?,
        _ => 100,
    };
    if parts.next().is_some() {
        return Err(ParseOutcome::Error(format!(
            "Too many fields in --user-tstv {raw}"
        )));
    }
    if nbins == 0 {
        return Err(ParseOutcome::Error(
            "--user-tstv bin count must be positive".into(),
        ));
    }
    if max <= min {
        return Err(ParseOutcome::Error(format!(
            "--user-tstv max must be greater than min ({min}:{max})"
        )));
    }
    Ok(UserTstvSpec {
        tag,
        index,
        min,
        max,
        nbins,
    })
}

fn parse_tag_index(raw: &str) -> Result<(String, usize), ParseOutcome> {
    if let Some(prefix) = raw.strip_suffix(']')
        && let Some((tag, index)) = prefix.rsplit_once('[')
    {
        let index = index.parse::<usize>().map_err(|e| {
            ParseOutcome::Error(format!("Could not parse --user-tstv index in '{raw}': {e}"))
        })?;
        if tag.is_empty() {
            return Err(ParseOutcome::Error(format!(
                "empty --user-tstv tag in '{raw}'"
            )));
        }
        return Ok((tag.to_owned(), index));
    }
    Ok((raw.to_owned(), 0))
}

fn parse_f64(raw: &str, name: &str) -> Result<f64, ParseOutcome> {
    raw.parse::<f64>()
        .map_err(|e| ParseOutcome::Error(format!("Could not parse {name} '{raw}': {e}")))
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
        .map_err(|e| ParseOutcome::Error(format!("Could not parse position '{raw}': {e}")))
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
                start += 1;
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

#[derive(Debug, Default, Clone)]
struct Counters {
    n_records: u64,
    n_no_alts: u64,
    n_snps: u64,
    n_mnps: u64,
    n_indels: u64,
    n_others: u64,
    n_multiallelic: u64,
    n_multiallelic_snp: u64,
    in_frame: u64,
    out_frame: u64,
    na_frame: u64,
    in_frame_alt1: u64,
    out_frame_alt1: u64,
    na_frame_alt1: u64,
    repeat_context: [[u64; 4]; INDEL_CONTEXT_REPEAT_LEN],
    repeat_na: u64,
    singleton_snps: u64,
    singleton_ts: u64,
    singleton_tv: u64,
    singleton_indels: u64,
    singleton_repeats: [u64; 3],
    ts: u64,
    tv: u64,
    ts_alt1: u64,
    tv_alt1: u64,
    quality_counts: std::collections::BTreeMap<String, QualityCounters>,
    indel_distribution: std::collections::BTreeMap<i64, IndelDistributionCounters>,
    hwe: std::collections::BTreeMap<String, Vec<f64>>,
    af_counts: Vec<AfBinCounters>,
    depth_genotypes: Vec<u64>,
    depth_sites: Vec<u64>,
    substitutions: std::collections::BTreeMap<String, u64>,
    user_tstv: Vec<UserTstvCounters>,
    sample_stats: Vec<SampleCounters>,
}

#[derive(Debug, Default, Clone, Copy)]
struct AfBinCounters {
    snps: u64,
    indels: u64,
    others: u64,
}

#[derive(Debug, Default, Clone, Copy)]
struct QualityCounters {
    snps: u64,
    ts_alt1: u64,
    tv_alt1: u64,
    indels: u64,
}

#[derive(Debug, Default, Clone, Copy)]
struct IndelDistributionCounters {
    sites: u64,
    genotypes: u64,
    vaf_sum: f64,
    vaf_observations: u64,
}

#[derive(Debug, Default, Clone)]
struct UserTstvCounters {
    ts: Vec<u64>,
    tv: Vec<u64>,
}

#[derive(Debug, Default, Clone)]
struct SampleCounters {
    ref_hom: u64,
    nonref_hom: u64,
    hets: u64,
    ts: u64,
    tv: u64,
    indels: u64,
    depth_sum: u64,
    depth_count: u64,
    singletons: u64,
    hap_ref: u64,
    hap_alt: u64,
    missing: u64,
    frame_na: u64,
    frame_in: u64,
    frame_out: u64,
    ins_hets: u64,
    del_hets: u64,
    ins_alt_homs: u64,
    del_alt_homs: u64,
    vaf_snv: [u64; 21],
    vaf_indel: [u64; 21],
}

fn run(args: &Args, _argv: &[OsString]) -> io::Result<()> {
    let reference = open_fasta_reference(args)?;
    if args.inputs.len() == 2 {
        return run_pairwise(args, reference.as_ref());
    }

    let path = &args.inputs[0];
    let text = read_vcf_text(path)?;
    let body = body_after_header(&text);

    let header_samples = sample_names_from_header(&text);
    let selected_samples = select_samples(&header_samples, args)?;
    let n_samples = selected_samples.names.len();

    let mut total = new_counters(args, n_samples);
    let mut known = new_counters(args, n_samples);
    let mut novel = new_counters(args, n_samples);

    for line in body.split_inclusive('\n') {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }
        if !record_in_regions(&fields, &args.regions, &args.targets) {
            continue;
        }
        if !apply_filters_pass(&fields, args.apply_filters.as_deref()) {
            continue;
        }
        if !expression_pass(&fields, args)? {
            continue;
        }
        accumulate(
            &mut total,
            &fields,
            args,
            &selected_samples.indices,
            reference.as_ref(),
        )?;
        if args.split_by_id {
            if fields[2] == "." || fields[2].is_empty() {
                accumulate(
                    &mut novel,
                    &fields,
                    args,
                    &selected_samples.indices,
                    reference.as_ref(),
                )?;
            } else {
                accumulate(
                    &mut known,
                    &fields,
                    args,
                    &selected_samples.indices,
                    reference.as_ref(),
                )?;
            }
        }
    }

    print_report(args, &selected_samples, &total, &known, &novel)
}

fn open_fasta_reference(args: &Args) -> io::Result<Option<FastaReference>> {
    args.fasta_ref
        .as_deref()
        .map(FastaReference::open)
        .transpose()
}

fn run_pairwise(args: &Args, reference: Option<&FastaReference>) -> io::Result<()> {
    if args.split_by_id {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "Only one file can be given with -I.",
        ));
    }

    let left_text = read_vcf_text(&args.inputs[0])?;
    let right_text = read_vcf_text(&args.inputs[1])?;
    let left_samples = sample_names_from_header(&left_text);
    let selected_samples = select_samples(&left_samples, args)?;
    let n_samples = selected_samples.names.len();

    let left = collect_pairwise_records(&left_text, args)?;
    let right = collect_pairwise_records(&right_text, args)?;

    let mut left_only = new_counters(args, n_samples);
    let mut right_only = new_counters(args, n_samples);
    let mut shared = new_counters(args, n_samples);

    for (key, left_records) in &left {
        if let Some(right_records) = right.get(key) {
            for fields in left_records.iter().take(right_records.len()) {
                let refs: Vec<&str> = fields.iter().map(String::as_str).collect();
                accumulate(
                    &mut shared,
                    &refs,
                    args,
                    &selected_samples.indices,
                    reference,
                )?;
            }
            if left_records.len() > right_records.len() {
                for fields in &left_records[right_records.len()..] {
                    let refs: Vec<&str> = fields.iter().map(String::as_str).collect();
                    accumulate(
                        &mut left_only,
                        &refs,
                        args,
                        &selected_samples.indices,
                        reference,
                    )?;
                }
            }
        } else {
            for fields in left_records {
                let refs: Vec<&str> = fields.iter().map(String::as_str).collect();
                accumulate(
                    &mut left_only,
                    &refs,
                    args,
                    &selected_samples.indices,
                    reference,
                )?;
            }
        }
    }

    for (key, right_records) in &right {
        let matched = left.get(key).map_or(0, Vec::len);
        if matched >= right_records.len() {
            continue;
        }
        for fields in &right_records[matched..] {
            let refs: Vec<&str> = fields.iter().map(String::as_str).collect();
            accumulate(
                &mut right_only,
                &refs,
                args,
                &selected_samples.indices,
                reference,
            )?;
        }
    }

    print_pairwise_report(args, &selected_samples, &left_only, &right_only, &shared)
}

fn print_pairwise_report(
    args: &Args,
    selected_samples: &SelectedSamples,
    left_only: &Counters,
    right_only: &Counters,
    shared: &Counters,
) -> io::Result<()> {
    let mut out = io::stdout().lock();
    writeln!(out, "# This file was produced by bcftools stats")?;
    writeln!(out, "# The command line was: bcftools stats")?;
    writeln!(out, "#")?;
    writeln!(out, "# Definition of sets:")?;
    writeln!(out, "# ID\t[2]id\t[3]tab-separated file names")?;
    writeln!(out, "ID\t0\t{}", args.inputs[0].display())?;
    writeln!(out, "ID\t1\t{}", args.inputs[1].display())?;
    writeln!(
        out,
        "ID\t2\t{}\t{}",
        args.inputs[0].display(),
        args.inputs[1].display()
    )?;

    print_sn_sets(
        &mut out,
        selected_samples.names.len(),
        &[(0, left_only), (1, right_only), (2, shared)],
    )?;
    print_tstv_sets(&mut out, &[(0, left_only), (1, right_only), (2, shared)])?;
    print_st_sets(&mut out, &[(0, left_only), (1, right_only), (2, shared)])?;
    print_af_sets(
        &mut out,
        args,
        &[(0, left_only), (1, right_only), (2, shared)],
    )?;
    print_quality_sets(&mut out, &[(0, left_only), (1, right_only), (2, shared)])?;
    print_indel_distribution_sets(&mut out, &[(0, left_only), (1, right_only), (2, shared)])?;
    print_hwe_sets(&mut out, &[(0, left_only), (1, right_only), (2, shared)])?;
    print_depth_sets(
        &mut out,
        args,
        &[(0, left_only), (1, right_only), (2, shared)],
    )?;
    print_fs_sets(
        &mut out,
        args,
        &[(0, left_only), (1, right_only), (2, shared)],
    )?;
    print_indel_context_sets(
        &mut out,
        args,
        &[(0, left_only), (1, right_only), (2, shared)],
    )?;
    print_sis_sets(&mut out, &[(0, left_only), (1, right_only), (2, shared)])?;
    if args.sample_list.is_some() {
        print_psc_sets(
            &mut out,
            selected_samples,
            &[(0, left_only), (1, right_only), (2, shared)],
        )?;
        print_psi_sets(
            &mut out,
            selected_samples,
            &[(0, left_only), (1, right_only), (2, shared)],
        )?;
        print_vaf_sets(
            &mut out,
            selected_samples,
            &[(0, left_only), (1, right_only), (2, shared)],
        )?;
    }
    print_user_tstv_sets(
        &mut out,
        args,
        &[(0, left_only), (1, right_only), (2, shared)],
    )?;
    Ok(())
}

fn collect_pairwise_records(
    text: &str,
    args: &Args,
) -> io::Result<std::collections::BTreeMap<PairwiseKey, Vec<Vec<String>>>> {
    let mut records: std::collections::BTreeMap<PairwiseKey, Vec<Vec<String>>> =
        std::collections::BTreeMap::new();
    for line in body_after_header(text).split_inclusive('\n') {
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r');
        if trimmed.is_empty() {
            continue;
        }
        let fields: Vec<&str> = trimmed.split('\t').collect();
        if fields.len() < 8 {
            continue;
        }
        if !record_in_regions(&fields, &args.regions, &args.targets) {
            continue;
        }
        if !apply_filters_pass(&fields, args.apply_filters.as_deref()) {
            continue;
        }
        if !expression_pass(&fields, args)? {
            continue;
        }
        records
            .entry(pairwise_key(&fields, args.collapse))
            .or_default()
            .push(fields.iter().map(|field| (*field).to_owned()).collect());
    }
    Ok(records)
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct PairwiseKey {
    chrom: String,
    pos: i64,
    class: String,
}

fn pairwise_key(fields: &[&str], collapse: CollapseMode) -> PairwiseKey {
    let chrom = fields[0].to_owned();
    let pos = fields[1].parse::<i64>().unwrap_or(-1);
    let reference = fields[3];
    let alts: Vec<&str> = if fields[4] == "." {
        Vec::new()
    } else {
        fields[4].split(',').collect()
    };
    let class = match collapse {
        CollapseMode::None => format!("{}>{}", reference, fields[4]),
        CollapseMode::Id => fields.get(2).copied().unwrap_or(".").to_owned(),
        CollapseMode::Any | CollapseMode::All | CollapseMode::Some => "site".to_owned(),
        CollapseMode::Snps => {
            if alts
                .iter()
                .any(|alt| classify(reference, alt) == VariantKind::Snp)
            {
                "snp".to_owned()
            } else {
                format!("{}>{}", reference, fields[4])
            }
        }
        CollapseMode::Indels => {
            if alts
                .iter()
                .any(|alt| classify(reference, alt) == VariantKind::Indel)
            {
                "indel".to_owned()
            } else {
                format!("{}>{}", reference, fields[4])
            }
        }
        CollapseMode::Both => {
            if alts.iter().any(|alt| {
                matches!(
                    classify(reference, alt),
                    VariantKind::Snp | VariantKind::Indel
                )
            }) {
                "variant".to_owned()
            } else {
                format!("{}>{}", reference, fields[4])
            }
        }
    };
    PairwiseKey { chrom, pos, class }
}

fn print_report(
    args: &Args,
    selected_samples: &SelectedSamples,
    total: &Counters,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let mut out = io::stdout().lock();
    writeln!(out, "# This file was produced by bcftools stats")?;
    writeln!(out, "# The command line was: bcftools stats")?;
    writeln!(out, "#")?;
    writeln!(out, "# Definition of sets:")?;
    writeln!(out, "# ID\t[2]id\t[3]tab-separated file names")?;
    writeln!(out, "ID\t0\t{}", args.inputs[0].display())?;

    print_sn(
        &mut out,
        selected_samples.names.len(),
        total,
        args.split_by_id,
        known,
        novel,
    )?;
    print_tstv(&mut out, total, args.split_by_id, known, novel)?;
    print_st(&mut out, total)?;
    print_af(&mut out, args, total)?;
    print_quality(&mut out, total, args.split_by_id, known, novel)?;
    print_indel_distribution(&mut out, total, args.split_by_id, known, novel)?;
    print_hwe(&mut out, total, args.split_by_id, known, novel)?;
    print_depth(&mut out, args, total, args.split_by_id, known, novel)?;
    print_fs(&mut out, args, total, args.split_by_id, known, novel)?;
    print_indel_context(&mut out, args, total, args.split_by_id, known, novel)?;
    print_sis(&mut out, total, args.split_by_id, known, novel)?;
    if args.sample_list.is_some() {
        print_psc(
            &mut out,
            selected_samples,
            total,
            args.split_by_id,
            known,
            novel,
        )?;
        print_psi(
            &mut out,
            selected_samples,
            total,
            args.split_by_id,
            known,
            novel,
        )?;
        print_vaf(
            &mut out,
            selected_samples,
            total,
            args.split_by_id,
            known,
            novel,
        )?;
    }
    print_user_tstv(&mut out, args, total, args.split_by_id, known, novel)?;
    Ok(())
}

fn print_sn<W: Write>(
    out: &mut W,
    n_samples: usize,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters, &str)> = if split {
        vec![(0, total, "all"), (1, known, "known"), (2, novel, "novel")]
    } else {
        vec![(0, total, "all")]
    };
    let sets: Vec<(usize, &Counters)> = sets.into_iter().map(|(id, c, _)| (id, c)).collect();
    print_sn_sets(out, n_samples, &sets)
}

fn print_sn_sets<W: Write>(
    out: &mut W,
    n_samples: usize,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    writeln!(out, "# SN, Summary numbers:")?;
    writeln!(out, "# SN\t[2]id\t[3]key\t[4]value")?;
    for &(id, c) in sets {
        writeln!(out, "SN\t{id}\tnumber of samples:\t{n_samples}")?;
        writeln!(out, "SN\t{id}\tnumber of records:\t{}", c.n_records)?;
        writeln!(out, "SN\t{id}\tnumber of no-ALTs:\t{}", c.n_no_alts)?;
        writeln!(out, "SN\t{id}\tnumber of SNPs:\t{}", c.n_snps)?;
        writeln!(out, "SN\t{id}\tnumber of MNPs:\t{}", c.n_mnps)?;
        writeln!(out, "SN\t{id}\tnumber of indels:\t{}", c.n_indels)?;
        writeln!(out, "SN\t{id}\tnumber of others:\t{}", c.n_others)?;
        writeln!(
            out,
            "SN\t{id}\tnumber of multiallelic sites:\t{}",
            c.n_multiallelic
        )?;
        writeln!(
            out,
            "SN\t{id}\tnumber of multiallelic SNP sites:\t{}",
            c.n_multiallelic_snp
        )?;
    }
    Ok(())
}

fn print_tstv<W: Write>(
    out: &mut W,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_tstv_sets(out, &sets)
}

fn print_tstv_sets<W: Write>(out: &mut W, sets: &[(usize, &Counters)]) -> io::Result<()> {
    writeln!(out, "# TSTV, transitions/transversions:")?;
    writeln!(
        out,
        "# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)"
    )?;
    for &(id, c) in sets {
        let ratio = if c.tv != 0 {
            c.ts as f64 / c.tv as f64
        } else {
            0.0
        };
        let ratio_alt1 = if c.tv_alt1 != 0 {
            c.ts_alt1 as f64 / c.tv_alt1 as f64
        } else {
            0.0
        };
        writeln!(
            out,
            "TSTV\t{id}\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}",
            c.ts, c.tv, ratio, c.ts_alt1, c.tv_alt1, ratio_alt1
        )?;
    }
    Ok(())
}

fn print_st<W: Write>(out: &mut W, total: &Counters) -> io::Result<()> {
    print_st_sets(out, &[(0, total)])
}

fn print_st_sets<W: Write>(out: &mut W, sets: &[(usize, &Counters)]) -> io::Result<()> {
    writeln!(out, "# ST, Substitution types:")?;
    writeln!(out, "# ST\t[2]id\t[3]type\t[4]count")?;
    for &(id, counters) in sets {
        for key in SUBSTITUTION_TYPES {
            let count = counters.substitutions.get(key).copied().unwrap_or(0);
            writeln!(out, "ST\t{id}\t{key}\t{count}")?;
        }
    }
    Ok(())
}

fn print_af<W: Write>(out: &mut W, args: &Args, total: &Counters) -> io::Result<()> {
    print_af_sets(out, args, &[(0, total)])
}

fn print_af_sets<W: Write>(
    out: &mut W,
    args: &Args,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    writeln!(out, "# AF, Stats by non-reference allele frequency:")?;
    writeln!(
        out,
        "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of indels\t[6]number of others"
    )?;
    for &(id, total) in sets {
        for (i, bin) in args.af_bins.iter().enumerate() {
            let counters = &total.af_counts[i];
            writeln!(
                out,
                "AF\t{id}\t{:.6}\t{}\t{}\t{}",
                bin, counters.snps, counters.indels, counters.others
            )?;
        }
        let overflow = &total.af_counts[args.af_bins.len()];
        if overflow.snps + overflow.indels + overflow.others > 0 {
            writeln!(
                out,
                "AF\t{id}\t{:.6}\t{}\t{}\t{}",
                1.0, overflow.snps, overflow.indels, overflow.others
            )?;
        }
    }
    Ok(())
}

fn print_quality<W: Write>(
    out: &mut W,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_quality_sets(out, &sets)
}

fn print_quality_sets<W: Write>(out: &mut W, sets: &[(usize, &Counters)]) -> io::Result<()> {
    writeln!(out, "# QUAL, Stats by quality")?;
    writeln!(
        out,
        "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels"
    )?;
    for &(id, counters) in sets {
        let mut qualities: Vec<(&String, &QualityCounters)> =
            counters.quality_counts.iter().collect();
        qualities.sort_by(|(left, _), (right, _)| quality_label_cmp(left, right));
        for (quality, counts) in qualities {
            writeln!(
                out,
                "QUAL\t{id}\t{quality}\t{}\t{}\t{}\t{}",
                counts.snps, counts.ts_alt1, counts.tv_alt1, counts.indels
            )?;
        }
    }
    Ok(())
}

fn print_indel_distribution<W: Write>(
    out: &mut W,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_indel_distribution_sets(out, &sets)
}

fn print_indel_distribution_sets<W: Write>(
    out: &mut W,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    writeln!(out, "# IDD, InDel distribution:")?;
    writeln!(
        out,
        "# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF"
    )?;
    for &(id, counters) in sets {
        for (length, counts) in &counters.indel_distribution {
            let mean_vaf = if counts.vaf_observations == 0 {
                ".".to_string()
            } else {
                format!("{:.2}", counts.vaf_sum / counts.vaf_observations as f64)
            };
            writeln!(
                out,
                "IDD\t{id}\t{length}\t{}\t{}\t{mean_vaf}",
                counts.sites, counts.genotypes
            )?;
        }
    }
    Ok(())
}

fn print_hwe<W: Write>(
    out: &mut W,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_hwe_sets(out, &sets)
}

fn print_hwe_sets<W: Write>(out: &mut W, sets: &[(usize, &Counters)]) -> io::Result<()> {
    writeln!(out, "# HWE")?;
    writeln!(
        out,
        "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile"
    )?;
    for &(id, counters) in sets {
        for (af, values) in &counters.hwe {
            if values.is_empty() {
                continue;
            }
            let mut sorted = values.clone();
            sorted.sort_by(|left, right| left.total_cmp(right));
            let q25 = percentile_value(&sorted, 0.25);
            let median = percentile_value(&sorted, 0.50);
            let q75 = percentile_value(&sorted, 0.75);
            writeln!(
                out,
                "HWE\t{id}\t{af}\t{}\t{q25:.6}\t{median:.6}\t{q75:.6}",
                sorted.len()
            )?;
        }
    }
    Ok(())
}

fn print_depth<W: Write>(
    out: &mut W,
    args: &Args,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_depth_sets(out, args, &sets)
}

fn print_depth_sets<W: Write>(
    out: &mut W,
    args: &Args,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    writeln!(out, "# DP, depth:")?;
    writeln!(out, "#   - set id, see above")?;
    writeln!(
        out,
        "#   - the depth bin, corresponds to the depth (unless --depth was given)"
    )?;
    writeln!(
        out,
        "#   - number of genotypes with this depth (zero unless sample columns are present)"
    )?;
    writeln!(out, "#   - fraction of genotypes with this depth")?;
    writeln!(out, "#   - number of sites with this depth")?;
    writeln!(out, "#   - fraction of sites with this depth")?;
    writeln!(out, "# DP, Depth distribution")?;
    writeln!(
        out,
        "# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)"
    )?;
    for &(id, counters) in sets {
        let genotype_sum: u64 = counters.depth_genotypes.iter().sum();
        let site_sum: u64 = counters.depth_sites.iter().sum();
        for i in 0..args.depth.nbins {
            let genotypes = counters.depth_genotypes[i];
            let sites = counters.depth_sites[i];
            if genotypes == 0 && sites == 0 {
                continue;
            }
            let genotype_pct = if genotype_sum == 0 {
                0.0
            } else {
                genotypes as f64 * 100.0 / genotype_sum as f64
            };
            let site_pct = if site_sum == 0 {
                0.0
            } else {
                sites as f64 * 100.0 / site_sum as f64
            };
            writeln!(
                out,
                "DP\t{id}\t{}\t{}\t{:.6}\t{}\t{:.6}",
                depth_bin_label(&args.depth, i),
                genotypes,
                genotype_pct,
                sites,
                site_pct
            )?;
        }
    }
    Ok(())
}

fn print_fs<W: Write>(
    out: &mut W,
    args: &Args,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_fs_sets(out, args, &sets)
}

fn print_fs_sets<W: Write>(
    out: &mut W,
    args: &Args,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    if args.exons.is_empty() {
        return Ok(());
    }
    writeln!(out, "# FS, Indel frameshifts:")?;
    writeln!(
        out,
        "# FS\t[2]id\t[3]in-frame\t[4]out-frame\t[5]not applicable\t[6]out/(in+out) ratio\t[7]in-frame (1st ALT)\t[8]out-frame (1st ALT)\t[9]not applicable (1st ALT)\t[10]out/(in+out) ratio (1st ALT)"
    )?;
    for &(id, counters) in sets {
        let ratio = ratio_or_zero(counters.out_frame, counters.in_frame);
        let ratio_alt1 = ratio_or_zero(counters.out_frame_alt1, counters.in_frame_alt1);
        writeln!(
            out,
            "FS\t{id}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{:.2}",
            counters.in_frame,
            counters.out_frame,
            counters.na_frame,
            ratio,
            counters.in_frame_alt1,
            counters.out_frame_alt1,
            counters.na_frame_alt1,
            ratio_alt1
        )?;
    }
    Ok(())
}

fn print_indel_context<W: Write>(
    out: &mut W,
    args: &Args,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_indel_context_sets(out, args, &sets)
}

fn print_indel_context_sets<W: Write>(
    out: &mut W,
    args: &Args,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    if args.fasta_ref.is_none() {
        return Ok(());
    }
    writeln!(out, "# ICS, Indel context:")?;
    writeln!(
        out,
        "#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]"
    )?;
    writeln!(
        out,
        "# ICS\t[2]id\t[3]repeat-consistent\t[4]repeat-inconsistent\t[5]not applicable\t[6]c/(c+i) ratio"
    )?;
    for &(id, counters) in sets {
        let consistent: u64 = counters
            .repeat_context
            .iter()
            .map(|counts| counts[0] + counts[2])
            .sum();
        let inconsistent: u64 = counters
            .repeat_context
            .iter()
            .map(|counts| counts[1] + counts[3])
            .sum();
        let ratio = if consistent + inconsistent == 0 {
            0.0
        } else {
            consistent as f64 / (consistent + inconsistent) as f64
        };
        writeln!(
            out,
            "ICS\t{id}\t{consistent}\t{inconsistent}\t{}\t{ratio:.4}",
            counters.repeat_na
        )?;
    }
    writeln!(out, "# ICL, Indel context by length:")?;
    writeln!(
        out,
        "#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]"
    )?;
    writeln!(
        out,
        "# ICL\t[2]id\t[3]length of repeat element\t[4]repeat-consistent deletions)\t[5]repeat-inconsistent deletions\t[6]consistent insertions\t[7]inconsistent insertions\t[8]c/(c+i) ratio"
    )?;
    for &(id, counters) in sets {
        for i in 1..INDEL_CONTEXT_REPEAT_LEN {
            let counts = counters.repeat_context[i];
            let consistent = counts[0] + counts[2];
            let inconsistent = counts[1] + counts[3];
            let ratio = if consistent + inconsistent == 0 {
                0.0
            } else {
                consistent as f64 / (consistent + inconsistent) as f64
            };
            writeln!(
                out,
                "ICL\t{id}\t{}\t{}\t{}\t{}\t{}\t{ratio:.4}",
                i + 1,
                counts[0],
                counts[1],
                counts[2],
                counts[3]
            )?;
        }
    }
    Ok(())
}

fn print_sis<W: Write>(
    out: &mut W,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_sis_sets(out, &sets)
}

fn print_sis_sets<W: Write>(out: &mut W, sets: &[(usize, &Counters)]) -> io::Result<()> {
    writeln!(out, "# SiS, Singleton stats:")?;
    writeln!(
        out,
        "#   - allele count, i.e. the number of singleton genotypes (AC=1)"
    )?;
    writeln!(out, "#   - number of transitions, see above")?;
    writeln!(out, "#   - number of transversions, see above")?;
    writeln!(
        out,
        "#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]"
    )?;
    writeln!(
        out,
        "# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable"
    )?;
    for &(id, counters) in sets {
        writeln!(
            out,
            "SiS\t{id}\t1\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            counters.singleton_snps,
            counters.singleton_ts,
            counters.singleton_tv,
            counters.singleton_indels,
            counters.singleton_repeats[0],
            counters.singleton_repeats[1],
            counters.singleton_repeats[2]
        )?;
    }
    Ok(())
}

fn ratio_or_zero(numerator: u64, denominator_other: u64) -> f64 {
    let denominator = numerator + denominator_other;
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn quality_label_cmp(left: &str, right: &str) -> std::cmp::Ordering {
    match (left.parse::<f64>(), right.parse::<f64>()) {
        (Ok(l), Ok(r)) => l.partial_cmp(&r).unwrap_or(std::cmp::Ordering::Equal),
        (Ok(_), Err(_)) => std::cmp::Ordering::Less,
        (Err(_), Ok(_)) => std::cmp::Ordering::Greater,
        (Err(_), Err(_)) => left.cmp(right),
    }
}

fn percentile_value(sorted: &[f64], percentile: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let index = ((sorted.len() - 1) as f64 * percentile).round() as usize;
    sorted[index]
}

fn hwe_af_bin_label(af: f64) -> String {
    let minor_af = af.min(1.0 - af);
    let binned = if minor_af == 0.0 {
        0.0
    } else if minor_af >= 0.5 {
        0.49
    } else {
        (minor_af * 100.0).floor() / 100.0
    };
    format!("{binned:.6}")
}

fn hardy_weinberg_exact(obs_hets: u64, obs_hom1: u64, obs_hom2: u64) -> f64 {
    let obs_homc = obs_hom1.min(obs_hom2) as usize;
    let obs_homr = obs_hom1.max(obs_hom2) as usize;
    let obs_hets = obs_hets as usize;
    let rare_copies = 2 * obs_homc + obs_hets;
    let genotypes = obs_hets + obs_homc + obs_homr;
    if genotypes == 0 {
        return 1.0;
    }
    let mut probs = vec![0.0; rare_copies + 1];
    let mut mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
    if (rare_copies & 1) != (mid & 1) {
        mid += 1;
    }
    probs[mid] = 1.0;
    let mut sum = probs[mid];

    let mut curr_hets = mid;
    let mut curr_homr = (rare_copies - mid) / 2;
    let mut curr_homc = genotypes - curr_hets - curr_homr;
    while curr_hets >= 2 {
        let next = curr_hets - 2;
        probs[next] = probs[curr_hets] * curr_hets as f64 * (curr_hets - 1) as f64
            / (4.0 * (curr_homr + 1) as f64 * (curr_homc + 1) as f64);
        sum += probs[next];
        curr_hets = next;
        curr_homr += 1;
        curr_homc += 1;
    }

    curr_hets = mid;
    curr_homr = (rare_copies - mid) / 2;
    curr_homc = genotypes - curr_hets - curr_homr;
    while curr_hets + 2 <= rare_copies {
        let next = curr_hets + 2;
        probs[next] = probs[curr_hets] * 4.0 * curr_homr as f64 * curr_homc as f64
            / ((curr_hets + 2) as f64 * (curr_hets + 1) as f64);
        sum += probs[next];
        curr_hets = next;
        curr_homr = curr_homr.saturating_sub(1);
        curr_homc = curr_homc.saturating_sub(1);
    }

    for probability in &mut probs {
        *probability /= sum;
    }
    let observed = probs[obs_hets];
    probs
        .into_iter()
        .filter(|probability| *probability <= observed + 1e-12)
        .sum::<f64>()
        .min(1.0)
}

fn print_psc<W: Write>(
    out: &mut W,
    selected_samples: &SelectedSamples,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_psc_sets(out, selected_samples, &sets)
}

fn print_psc_sets<W: Write>(
    out: &mut W,
    selected_samples: &SelectedSamples,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    writeln!(
        out,
        "# PSC, Per-sample counts. Ref/het/hom counts include SNPs; indels are counted separately."
    )?;
    writeln!(
        out,
        "# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing"
    )?;
    for &(id, counters) in sets {
        for (sample_i, name) in selected_samples.names.iter().enumerate() {
            let c = &counters.sample_stats[sample_i];
            let avg_depth = if c.depth_count == 0 {
                0.0
            } else {
                c.depth_sum as f64 / c.depth_count as f64
            };
            writeln!(
                out,
                "PSC\t{id}\t{name}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}",
                c.ref_hom,
                c.nonref_hom,
                c.hets,
                c.ts,
                c.tv,
                c.indels,
                avg_depth,
                c.singletons,
                c.hap_ref,
                c.hap_alt,
                c.missing
            )?;
        }
    }
    Ok(())
}

fn print_psi<W: Write>(
    out: &mut W,
    selected_samples: &SelectedSamples,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_psi_sets(out, selected_samples, &sets)
}

fn print_psi_sets<W: Write>(
    out: &mut W,
    selected_samples: &SelectedSamples,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    writeln!(
        out,
        "# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets."
    )?;
    writeln!(
        out,
        "# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms"
    )?;
    for &(id, counters) in sets {
        for (sample_i, name) in selected_samples.names.iter().enumerate() {
            let c = &counters.sample_stats[sample_i];
            let frame_ratio = ratio_or_zero(c.frame_out, c.frame_in);
            writeln!(
                out,
                "PSI\t{id}\t{name}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}",
                c.frame_in,
                c.frame_out,
                c.frame_na,
                frame_ratio,
                c.ins_hets,
                c.del_hets,
                c.ins_alt_homs,
                c.del_alt_homs
            )?;
        }
    }
    Ok(())
}

fn print_vaf<W: Write>(
    out: &mut W,
    selected_samples: &SelectedSamples,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_vaf_sets(out, selected_samples, &sets)
}

fn print_vaf_sets<W: Write>(
    out: &mut W,
    selected_samples: &SelectedSamples,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    if selected_samples.names.is_empty() {
        return Ok(());
    }
    writeln!(
        out,
        "# VAF, Variant Allele Frequency determined as fraction of alternate reads in FORMAT/AD"
    )?;
    writeln!(
        out,
        "# VAF\t[2]id\t[3]sample\t[4]SNV VAF distribution\t[5]indel VAF distribution"
    )?;
    for &(id, counters) in sets {
        for (sample_i, name) in selected_samples.names.iter().enumerate() {
            let c = &counters.sample_stats[sample_i];
            writeln!(
                out,
                "VAF\t{id}\t{name}\t{}\t{}",
                join_counts(&c.vaf_snv),
                join_counts(&c.vaf_indel)
            )?;
        }
    }
    Ok(())
}

fn join_counts<const N: usize>(counts: &[u64; N]) -> String {
    counts
        .iter()
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn print_user_tstv<W: Write>(
    out: &mut W,
    args: &Args,
    total: &Counters,
    split: bool,
    known: &Counters,
    novel: &Counters,
) -> io::Result<()> {
    let sets: Vec<(usize, &Counters)> = if split {
        vec![(0, total), (1, known), (2, novel)]
    } else {
        vec![(0, total)]
    };
    print_user_tstv_sets(out, args, &sets)
}

fn print_user_tstv_sets<W: Write>(
    out: &mut W,
    args: &Args,
    sets: &[(usize, &Counters)],
) -> io::Result<()> {
    for (spec_i, spec) in args.user_tstv.iter().enumerate() {
        writeln!(
            out,
            "# USR:{}/{}\t[2]id\t[3]{}/{}\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)",
            spec.tag, spec.index, spec.tag, spec.index
        )?;
        for &(id, counters) in sets {
            let user = &counters.user_tstv[spec_i];
            for bin in 0..spec.nbins {
                let ts = user.ts[bin];
                let tv = user.tv[bin];
                if ts + tv == 0 {
                    continue;
                }
                let value = user_tstv_bin_value(spec, bin);
                writeln!(
                    out,
                    "USR:{}/{}\t{}\t{:.6}\t{}\t{}\t{}",
                    spec.tag,
                    spec.index,
                    id,
                    value,
                    ts + tv,
                    ts,
                    tv
                )?;
            }
        }
    }
    Ok(())
}

fn new_counters(args: &Args, n_samples: usize) -> Counters {
    Counters {
        af_counts: vec![AfBinCounters::default(); args.af_bins.len() + 1],
        depth_genotypes: vec![0; args.depth.nbins],
        depth_sites: vec![0; args.depth.nbins],
        sample_stats: vec![SampleCounters::default(); n_samples],
        user_tstv: args
            .user_tstv
            .iter()
            .map(|spec| UserTstvCounters {
                ts: vec![0; spec.nbins],
                tv: vec![0; spec.nbins],
            })
            .collect(),
        ..Counters::default()
    }
}

fn accumulate(
    c: &mut Counters,
    fields: &[&str],
    args: &Args,
    selected_sample_indices: &[usize],
    fasta_ref: Option<&FastaReference>,
) -> io::Result<()> {
    c.n_records += 1;
    let reference = fields[3];
    let alts: Vec<&str> = if fields[4] == "." {
        Vec::new()
    } else {
        fields[4].split(',').collect()
    };
    if alts.is_empty() {
        c.n_no_alts += 1;
        return Ok(());
    }
    if alts.len() > 1 {
        c.n_multiallelic += 1;
        if alts.iter().all(|a| a.len() == 1 && reference.len() == 1) {
            c.n_multiallelic_snp += 1;
        }
    }
    let mut saw_snp = false;
    let mut saw_mnp = false;
    let mut saw_indel = false;
    let mut saw_other = false;
    let mut first_alt_tstv: Option<bool> = None;
    let scan_alts: &[&str] = if args.first_allele_only && !alts.is_empty() {
        &alts[..1]
    } else {
        &alts
    };
    for (i, alt) in scan_alts.iter().enumerate() {
        match classify(reference, alt) {
            VariantKind::Snp => {
                saw_snp = true;
                if let Some((ts, tv)) = ts_tv_count(reference, alt) {
                    c.ts += ts;
                    c.tv += tv;
                    if i == 0 {
                        c.ts_alt1 += ts;
                        c.tv_alt1 += tv;
                        first_alt_tstv = Some(ts > 0);
                    }
                }
                let key = format!("{}>{}", reference, alt).to_ascii_uppercase();
                *c.substitutions.entry(key).or_default() += 1;
            }
            VariantKind::Mnp => saw_mnp = true,
            VariantKind::Indel => saw_indel = true,
            VariantKind::Other => saw_other = true,
        }
        if !args.exons.is_empty() && classify(reference, alt) == VariantKind::Indel {
            let frame = indel_frame(fields[0], fields[1], reference, alt, &args.exons);
            accumulate_frame_counts(c, i == 0, frame);
        }
        if let Some(fasta) = fasta_ref
            && classify(reference, alt) == VariantKind::Indel
        {
            accumulate_indel_context(c, fields, reference, alt, fasta)?;
        }
    }
    if saw_snp {
        c.n_snps += 1;
    }
    if saw_mnp {
        c.n_mnps += 1;
    }
    if saw_indel {
        c.n_indels += 1;
    }
    if saw_other {
        c.n_others += 1;
    }
    accumulate_quality_counts(c, fields[5], saw_snp, saw_indel, first_alt_tstv);
    accumulate_indel_distribution(c, fields, reference, scan_alts, selected_sample_indices);
    accumulate_hwe(c, fields, reference, scan_alts, selected_sample_indices);
    if let Some(is_ts) = first_alt_tstv {
        accumulate_user_tstv(c, fields[7], args, is_ts);
    }
    accumulate_singleton_stats(
        c,
        fields,
        reference,
        &alts,
        selected_sample_indices,
        fasta_ref,
    )?;
    if let Some(dp) = info_int_tag(fields[7], "DP") {
        let bin = depth_bin_index(&args.depth, dp);
        c.depth_sites[bin] += 1;
    }
    for dp in sample_depths(fields) {
        if dp > 0 {
            let bin = depth_bin_index(&args.depth, dp);
            c.depth_genotypes[bin] += 1;
        }
    }
    accumulate_sample_stats(
        c,
        fields,
        reference,
        &alts,
        selected_sample_indices,
        &args.exons,
    );
    let af = info_numeric_tag(fields[7], &args.af_tag).unwrap_or(0.0);
    let bin_index = af_bin_index(af, &args.af_bins);
    let bucket = &mut c.af_counts[bin_index];
    if saw_snp {
        bucket.snps += 1;
    }
    if saw_indel {
        bucket.indels += 1;
    }
    if saw_other || saw_mnp {
        bucket.others += 1;
    }
    Ok(())
}

fn accumulate_quality_counts(
    counters: &mut Counters,
    quality: &str,
    saw_snp: bool,
    saw_indel: bool,
    first_alt_tstv: Option<bool>,
) {
    let bucket = counters
        .quality_counts
        .entry(quality.to_string())
        .or_default();
    if saw_snp {
        bucket.snps += 1;
    }
    if saw_indel {
        bucket.indels += 1;
    }
    match first_alt_tstv {
        Some(true) => bucket.ts_alt1 += 1,
        Some(false) => bucket.tv_alt1 += 1,
        None => {}
    }
}

fn accumulate_indel_distribution(
    counters: &mut Counters,
    fields: &[&str],
    reference: &str,
    alts: &[&str],
    selected_sample_indices: &[usize],
) {
    let indel_lengths: Vec<(usize, i64)> = alts
        .iter()
        .enumerate()
        .filter_map(|(i, alt)| indel_length(reference, alt).map(|length| (i + 1, length)))
        .collect();
    if indel_lengths.is_empty() {
        return;
    }
    for &(_, length) in &indel_lengths {
        counters.indel_distribution.entry(length).or_default().sites += 1;
    }
    if selected_sample_indices.is_empty() || fields.len() <= 9 {
        return;
    }

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let gt_index = format_keys.iter().position(|key| *key == "GT");
    let ad_index = format_keys.iter().position(|key| *key == "AD");
    let Some(gt_index) = gt_index else {
        return;
    };
    for &sample_i in selected_sample_indices {
        let Some(sample) = fields.get(9 + sample_i) else {
            continue;
        };
        let values: Vec<&str> = sample.split(':').collect();
        let Some(gt) = values.get(gt_index) else {
            continue;
        };
        let called: Vec<usize> = parse_gt_alleles(gt).into_iter().flatten().collect();
        if called.is_empty() {
            continue;
        }
        let ad_values = ad_index.and_then(|i| values.get(i).and_then(|raw| parse_ad_values(raw)));
        let ad_total = ad_values
            .as_ref()
            .map(|values| values.iter().sum::<i64>())
            .filter(|total| *total > 0);
        let mut lengths_for_sample: std::collections::BTreeMap<
            i64,
            std::collections::BTreeSet<usize>,
        > = std::collections::BTreeMap::new();
        for allele in called.into_iter().filter(|&allele| allele > 0) {
            if let Some((_, length)) = indel_lengths
                .iter()
                .find(|(indel_allele, _)| *indel_allele == allele)
            {
                lengths_for_sample
                    .entry(*length)
                    .or_default()
                    .insert(allele);
            }
        }
        for (length, alleles) in lengths_for_sample {
            let bucket = counters.indel_distribution.entry(length).or_default();
            bucket.genotypes += 1;
            if let (Some(ad_values), Some(ad_total)) = (&ad_values, ad_total) {
                let alt_depth: i64 = alleles
                    .iter()
                    .filter_map(|&allele| ad_values.get(allele).copied())
                    .filter(|depth| *depth > 0)
                    .sum();
                if alt_depth > 0 {
                    bucket.vaf_sum += alt_depth as f64 / ad_total as f64;
                    bucket.vaf_observations += 1;
                }
            }
        }
    }
}

fn accumulate_hwe(
    counters: &mut Counters,
    fields: &[&str],
    reference: &str,
    alts: &[&str],
    selected_sample_indices: &[usize],
) {
    if selected_sample_indices.is_empty()
        || fields.len() <= 9
        || alts
            .first()
            .is_none_or(|alt| classify(reference, alt) != VariantKind::Snp)
    {
        return;
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let Some(gt_index) = format_keys.iter().position(|key| *key == "GT") else {
        return;
    };
    let mut ref_hom = 0_u64;
    let mut het = 0_u64;
    let mut alt_hom = 0_u64;
    for &sample_i in selected_sample_indices {
        let Some(sample) = fields.get(9 + sample_i) else {
            continue;
        };
        let values: Vec<&str> = sample.split(':').collect();
        let Some(gt) = values.get(gt_index) else {
            continue;
        };
        let alleles: Vec<usize> = parse_gt_alleles(gt).into_iter().flatten().collect();
        if alleles.len() != 2 {
            continue;
        }
        match (alleles[0], alleles[1]) {
            (0, 0) => ref_hom += 1,
            (1, 1) => alt_hom += 1,
            (0, 1) | (1, 0) => het += 1,
            _ => {}
        }
    }
    let called = ref_hom + het + alt_hom;
    if called == 0 {
        return;
    }
    let alt_alleles = het + 2 * alt_hom;
    let af = alt_alleles as f64 / (2 * called) as f64;
    let p_value = hardy_weinberg_exact(het, ref_hom, alt_hom);
    counters
        .hwe
        .entry(hwe_af_bin_label(af))
        .or_default()
        .push(p_value);
}

fn accumulate_sample_stats(
    c: &mut Counters,
    fields: &[&str],
    reference: &str,
    alts: &[&str],
    selected_sample_indices: &[usize],
    exons: &[RegionSpec],
) {
    if selected_sample_indices.is_empty() || fields.len() <= 9 {
        return;
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let gt_index = format_keys.iter().position(|key| *key == "GT");
    let site_is_snp = alts
        .iter()
        .any(|alt| classify(reference, alt) == VariantKind::Snp);
    let selected_alt_counts =
        selected_alt_allele_counts(fields, gt_index, selected_sample_indices, alts.len());
    for (out_i, &sample_i) in selected_sample_indices.iter().enumerate() {
        let Some(sample) = fields.get(9 + sample_i) else {
            continue;
        };
        if let Some(dp) = sample_depth(fields, sample_i)
            && dp > 0
        {
            c.sample_stats[out_i].depth_sum += dp as u64;
            c.sample_stats[out_i].depth_count += 1;
        }
        let Some(gt_index) = gt_index else {
            continue;
        };
        let values: Vec<&str> = sample.split(':').collect();
        let Some(gt) = values.get(gt_index) else {
            continue;
        };
        let alleles = parse_gt_alleles(gt);
        if alleles.is_empty() || alleles.iter().any(Option::is_none) {
            c.sample_stats[out_i].missing += 1;
            continue;
        }
        let called: Vec<usize> = alleles.into_iter().flatten().collect();
        if called.len() == 1 {
            if called[0] == 0 {
                c.sample_stats[out_i].hap_ref += 1;
            } else {
                c.sample_stats[out_i].hap_alt += 1;
            }
        }
        let mut singleton_alleles = std::collections::BTreeSet::new();
        for allele in called.iter().copied().filter(|&a| a > 0) {
            if selected_alt_counts
                .get(allele - 1)
                .is_some_and(|count| *count == 1)
            {
                singleton_alleles.insert(allele);
            }
        }
        c.sample_stats[out_i].singletons += singleton_alleles.len() as u64;
        if called.iter().all(|&a| a == 0) {
            if site_is_snp {
                c.sample_stats[out_i].ref_hom += 1;
            }
            continue;
        }
        accumulate_sample_vaf(
            &mut c.sample_stats[out_i],
            reference,
            alts,
            &called,
            &values,
            &format_keys,
        );
        let has_indel = called.iter().copied().filter(|&a| a > 0).any(|a| {
            alts.get(a - 1)
                .is_some_and(|alt| classify(reference, alt) == VariantKind::Indel)
        });
        if has_indel {
            c.sample_stats[out_i].indels += 1;
            accumulate_sample_indel_shape(
                &mut c.sample_stats[out_i],
                fields,
                &called,
                reference,
                alts,
                exons,
            );
        }
        let all_same = called.windows(2).all(|w| w[0] == w[1]);
        if site_is_snp && !called.is_empty() {
            if all_same && called[0] > 0 {
                c.sample_stats[out_i].nonref_hom += 1;
            } else if !all_same {
                c.sample_stats[out_i].hets += 1;
            }
            let mut counted = std::collections::BTreeSet::new();
            for allele in called.iter().copied().filter(|&a| a > 0) {
                if !counted.insert(allele) {
                    continue;
                }
                let Some(alt) = alts.get(allele - 1) else {
                    continue;
                };
                if classify(reference, alt) == VariantKind::Snp
                    && let Some((ts, tv)) = ts_tv_count(reference, alt)
                {
                    c.sample_stats[out_i].ts += ts;
                    c.sample_stats[out_i].tv += tv;
                }
            }
        }
    }
}

fn accumulate_sample_vaf(
    counters: &mut SampleCounters,
    reference: &str,
    alts: &[&str],
    called: &[usize],
    sample_values: &[&str],
    format_keys: &[&str],
) {
    let Some(ad_index) = format_keys.iter().position(|key| *key == "AD") else {
        return;
    };
    let Some(ad_values) = sample_values
        .get(ad_index)
        .and_then(|raw| parse_ad_values(raw))
    else {
        return;
    };
    let total: i64 = ad_values.iter().sum();
    if total <= 0 {
        return;
    }
    let mut counted = std::collections::BTreeSet::new();
    for allele in called.iter().copied().filter(|&allele| allele > 0) {
        if !counted.insert(allele) {
            continue;
        }
        let Some(depth) = ad_values.get(allele).copied() else {
            continue;
        };
        if depth < 0 {
            continue;
        }
        let Some(alt) = alts.get(allele - 1).copied() else {
            continue;
        };
        let bin = vaf_bin(depth as f64 / total as f64);
        match classify(reference, alt) {
            VariantKind::Snp => counters.vaf_snv[bin] += 1,
            VariantKind::Indel => counters.vaf_indel[bin] += 1,
            VariantKind::Mnp | VariantKind::Other => {}
        }
    }
}

fn parse_ad_values(raw: &str) -> Option<Vec<i64>> {
    let mut values = Vec::new();
    for part in raw.split(',') {
        if part == "." || part.is_empty() {
            return None;
        }
        values.push(part.parse::<i64>().ok()?);
    }
    Some(values)
}

fn vaf_bin(vaf: f64) -> usize {
    (vaf.clamp(0.0, 1.0) / 0.05).round() as usize
}

fn accumulate_sample_indel_shape(
    counters: &mut SampleCounters,
    fields: &[&str],
    called: &[usize],
    reference: &str,
    alts: &[&str],
    exons: &[RegionSpec],
) {
    if called.len() < 2 || called.iter().all(|&allele| allele == 0) {
        return;
    }
    let all_same = called.windows(2).all(|window| window[0] == window[1]);
    if all_same {
        let allele = called[0];
        if allele == 0 {
            return;
        }
        for allele in called.iter().copied().filter(|&allele| allele > 0) {
            accumulate_sample_frame(counters, fields, reference, alts, allele, exons);
        }
        match indel_direction(reference, alts.get(allele - 1).copied()) {
            Some(IndelDirection::Insertion) => counters.ins_alt_homs += 1,
            Some(IndelDirection::Deletion) => counters.del_alt_homs += 1,
            None => {}
        }
        return;
    }

    let mut has_insertion = false;
    let mut has_deletion = false;
    for allele in called.iter().copied().filter(|&allele| allele > 0) {
        accumulate_sample_frame(counters, fields, reference, alts, allele, exons);
        match indel_direction(reference, alts.get(allele - 1).copied()) {
            Some(IndelDirection::Insertion) => has_insertion = true,
            Some(IndelDirection::Deletion) => has_deletion = true,
            None => {}
        }
    }
    if has_insertion {
        counters.ins_hets += 1;
    }
    if has_deletion {
        counters.del_hets += 1;
    }
}

fn accumulate_sample_frame(
    counters: &mut SampleCounters,
    fields: &[&str],
    reference: &str,
    alts: &[&str],
    allele: usize,
    exons: &[RegionSpec],
) {
    let Some(alt) = alts.get(allele - 1) else {
        return;
    };
    if exons.is_empty() {
        return;
    }
    match indel_frame(fields[0], fields[1], reference, alt, exons) {
        FrameShift::NotApplicable => counters.frame_na += 1,
        FrameShift::InFrame => counters.frame_in += 1,
        FrameShift::OutFrame => counters.frame_out += 1,
    }
}

fn accumulate_frame_counts(counters: &mut Counters, is_alt1: bool, frame: FrameShift) {
    match frame {
        FrameShift::NotApplicable => {
            counters.na_frame += 1;
            if is_alt1 {
                counters.na_frame_alt1 += 1;
            }
        }
        FrameShift::InFrame => {
            counters.in_frame += 1;
            if is_alt1 {
                counters.in_frame_alt1 += 1;
            }
        }
        FrameShift::OutFrame => {
            counters.out_frame += 1;
            if is_alt1 {
                counters.out_frame_alt1 += 1;
            }
        }
    }
}

fn accumulate_indel_context(
    counters: &mut Counters,
    fields: &[&str],
    reference: &str,
    alt: &str,
    fasta: &FastaReference,
) -> io::Result<()> {
    match indel_context_class(fields, reference, alt, fasta)? {
        Some(IndelContextClass::ConsistentDeletion { repeat_len }) => {
            counters.repeat_context[repeat_len - 1][0] += 1;
        }
        Some(IndelContextClass::InconsistentDeletion { repeat_len }) => {
            counters.repeat_context[repeat_len - 1][1] += 1;
        }
        Some(IndelContextClass::ConsistentInsertion { repeat_len }) => {
            counters.repeat_context[repeat_len - 1][2] += 1;
        }
        Some(IndelContextClass::InconsistentInsertion { repeat_len }) => {
            counters.repeat_context[repeat_len - 1][3] += 1;
        }
        Some(IndelContextClass::NotApplicable) | None => counters.repeat_na += 1,
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IndelContextClass {
    ConsistentDeletion { repeat_len: usize },
    InconsistentDeletion { repeat_len: usize },
    ConsistentInsertion { repeat_len: usize },
    InconsistentInsertion { repeat_len: usize },
    NotApplicable,
}

impl IndelContextClass {
    fn singleton_bucket(self) -> usize {
        match self {
            Self::ConsistentDeletion { .. } | Self::ConsistentInsertion { .. } => 0,
            Self::InconsistentDeletion { .. } | Self::InconsistentInsertion { .. } => 1,
            Self::NotApplicable => 2,
        }
    }
}

fn indel_context_class(
    fields: &[&str],
    reference: &str,
    alt: &str,
    fasta: &FastaReference,
) -> io::Result<Option<IndelContextClass>> {
    let delta = alt.len() as i64 - reference.len() as i64;
    if delta == 0 {
        return Ok(None);
    }
    let Some((repeat_count, repeat_len)) = indel_repeat_context(fields[0], fields[1], fasta)?
    else {
        return Ok(Some(IndelContextClass::NotApplicable));
    };
    if repeat_count <= 1 || !(2..=INDEL_CONTEXT_REPEAT_LEN).contains(&repeat_len) {
        return Ok(Some(IndelContextClass::NotApplicable));
    }
    let consistent = delta.unsigned_abs().is_multiple_of(repeat_len as u64);
    let class = match (delta < 0, consistent) {
        (true, true) => IndelContextClass::ConsistentDeletion { repeat_len },
        (true, false) => IndelContextClass::InconsistentDeletion { repeat_len },
        (false, true) => IndelContextClass::ConsistentInsertion { repeat_len },
        (false, false) => IndelContextClass::InconsistentInsertion { repeat_len },
    };
    Ok(Some(class))
}

fn indel_repeat_context(
    chrom: &str,
    pos_raw: &str,
    fasta: &FastaReference,
) -> io::Result<Option<(usize, usize)>> {
    let Ok(pos) = pos_raw.parse::<i64>() else {
        return Ok(None);
    };
    let Some(sequence_len) = fasta.sequence_len(chrom) else {
        return Ok(None);
    };
    if pos < 1 || pos as u64 > sequence_len {
        return Ok(None);
    }
    let end = (pos + 50).min(sequence_len as i64);
    let region = format!("{chrom}:{pos}-{end}");
    let seq = fasta.fetch_region(&region)?;
    Ok(best_tandem_repeat(&seq))
}

fn best_tandem_repeat(seq: &[u8]) -> Option<(usize, usize)> {
    if seq.is_empty() {
        return None;
    }
    let mut best_count = 1;
    let mut best_len = 1;
    for unit_len in 1..=INDEL_CONTEXT_REPEAT_LEN.min(seq.len()) {
        let mut offset = 0;
        while offset + unit_len <= seq.len() {
            let unit = &seq[offset..offset + unit_len];
            let mut count = 1;
            let mut next = offset + unit_len;
            while next + unit_len <= seq.len()
                && seq[next..next + unit_len].eq_ignore_ascii_case(unit)
            {
                count += 1;
                next += unit_len;
            }
            if count > best_count || (count == best_count && unit_len > best_len) {
                best_count = count;
                best_len = unit_len;
            }
            offset += 1;
        }
    }
    Some((best_count, best_len))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FrameShift {
    NotApplicable,
    InFrame,
    OutFrame,
}

fn indel_frame(
    chrom: &str,
    pos_raw: &str,
    reference: &str,
    alt: &str,
    exons: &[RegionSpec],
) -> FrameShift {
    if exons.is_empty() || classify(reference, alt) != VariantKind::Indel {
        return FrameShift::NotApplicable;
    }
    let Ok(pos) = pos_raw.parse::<i64>() else {
        return FrameShift::NotApplicable;
    };
    let Some(affected_len) = exon_affected_indel_len(chrom, pos, reference, alt, exons) else {
        return FrameShift::NotApplicable;
    };
    if affected_len % 3 == 0 {
        FrameShift::InFrame
    } else {
        FrameShift::OutFrame
    }
}

fn exon_affected_indel_len(
    chrom: &str,
    pos: i64,
    reference: &str,
    alt: &str,
    exons: &[RegionSpec],
) -> Option<usize> {
    let len_delta = alt.len() as i64 - reference.len() as i64;
    if len_delta == 0 {
        return None;
    }
    if len_delta > 0 {
        return exons
            .iter()
            .any(|exon| exon_contains_pos(exon, chrom, pos))
            .then_some(len_delta.unsigned_abs() as usize);
    }

    let deleted_start = pos + 1;
    let deleted_end = pos + len_delta.unsigned_abs() as i64;
    let mut affected = 0usize;
    for exon in exons.iter().filter(|exon| exon.contig == chrom) {
        let start = exon.start.unwrap_or(i64::MIN).max(deleted_start);
        let end = exon.end.unwrap_or(i64::MAX).min(deleted_end);
        if start <= end {
            affected += (end - start + 1) as usize;
        }
    }
    (affected > 0).then_some(affected)
}

fn exon_contains_pos(exon: &RegionSpec, chrom: &str, pos: i64) -> bool {
    exon.contig == chrom
        && exon.start.map(|start| pos >= start).unwrap_or(true)
        && exon.end.map(|end| pos <= end).unwrap_or(true)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum IndelDirection {
    Insertion,
    Deletion,
}

fn indel_direction(reference: &str, alt: Option<&str>) -> Option<IndelDirection> {
    let alt = alt?;
    if classify(reference, alt) != VariantKind::Indel {
        return None;
    }
    if alt.len() >= reference.len() {
        Some(IndelDirection::Insertion)
    } else {
        Some(IndelDirection::Deletion)
    }
}

fn indel_length(reference: &str, alt: &str) -> Option<i64> {
    (classify(reference, alt) == VariantKind::Indel)
        .then_some(alt.len() as i64 - reference.len() as i64)
}

fn accumulate_singleton_stats(
    counters: &mut Counters,
    fields: &[&str],
    reference: &str,
    alts: &[&str],
    selected_sample_indices: &[usize],
    fasta_ref: Option<&FastaReference>,
) -> io::Result<()> {
    if selected_sample_indices.is_empty() || fields.len() <= 9 {
        return Ok(());
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let gt_index = format_keys.iter().position(|key| *key == "GT");
    let alt_counts =
        selected_alt_allele_counts(fields, gt_index, selected_sample_indices, alts.len());
    for (alt_i, count) in alt_counts.into_iter().enumerate() {
        if count != 1 {
            continue;
        }
        let Some(alt) = alts.get(alt_i).copied() else {
            continue;
        };
        match classify(reference, alt) {
            VariantKind::Snp => {
                counters.singleton_snps += 1;
                if let Some((ts, tv)) = ts_tv_count(reference, alt) {
                    counters.singleton_ts += ts;
                    counters.singleton_tv += tv;
                }
            }
            VariantKind::Indel => {
                counters.singleton_indels += 1;
                let bucket = match fasta_ref
                    .map(|fasta| indel_context_class(fields, reference, alt, fasta))
                    .transpose()?
                    .flatten()
                {
                    Some(class) => class.singleton_bucket(),
                    None => 2,
                };
                counters.singleton_repeats[bucket] += 1;
            }
            VariantKind::Mnp | VariantKind::Other => {}
        }
    }
    Ok(())
}

fn selected_alt_allele_counts(
    fields: &[&str],
    gt_index: Option<usize>,
    selected_sample_indices: &[usize],
    n_alts: usize,
) -> Vec<u64> {
    let mut counts = vec![0; n_alts];
    let Some(gt_index) = gt_index else {
        return counts;
    };
    for &sample_i in selected_sample_indices {
        let Some(sample) = fields.get(9 + sample_i) else {
            continue;
        };
        let values: Vec<&str> = sample.split(':').collect();
        let Some(gt) = values.get(gt_index) else {
            continue;
        };
        for allele in parse_gt_alleles(gt).into_iter().flatten() {
            if allele > 0
                && let Some(count) = counts.get_mut(allele - 1)
            {
                *count += 1;
            }
        }
    }
    counts
}

fn accumulate_user_tstv(c: &mut Counters, info: &str, args: &Args, is_ts: bool) {
    for (i, spec) in args.user_tstv.iter().enumerate() {
        let Some(value) = info_numeric_tag_index(info, &spec.tag, spec.index) else {
            continue;
        };
        let bin = user_tstv_bin_index(spec, value);
        if is_ts {
            c.user_tstv[i].ts[bin] += 1;
        } else {
            c.user_tstv[i].tv[bin] += 1;
        }
    }
}

fn user_tstv_bin_index(spec: &UserTstvSpec, value: f64) -> usize {
    if value <= spec.min {
        0
    } else if value >= spec.max {
        spec.nbins - 1
    } else {
        let frac = (value - spec.min) / (spec.max - spec.min);
        (frac * (spec.nbins - 1) as f64) as usize
    }
}

fn user_tstv_bin_value(spec: &UserTstvSpec, bin: usize) -> f64 {
    if spec.nbins <= 1 {
        spec.min
    } else {
        spec.min + (spec.max - spec.min) * bin as f64 / (spec.nbins - 1) as f64
    }
}

fn af_bin_index(af: f64, bins: &[f64]) -> usize {
    for (i, b) in bins.iter().enumerate() {
        if af <= *b {
            return i;
        }
    }
    bins.len()
}

fn info_numeric_tag(info: &str, tag: &str) -> Option<f64> {
    info_numeric_tag_index(info, tag, 0)
}

fn info_numeric_tag_index(info: &str, tag: &str, index: usize) -> Option<f64> {
    if info == "." {
        return None;
    }
    for entry in info.split(';') {
        let Some((key, value)) = entry.split_once('=') else {
            continue;
        };
        if key == tag {
            return value.split(',').nth(index)?.parse::<f64>().ok();
        }
    }
    None
}

fn info_int_tag(info: &str, tag: &str) -> Option<i64> {
    info_numeric_tag(info, tag).map(|value| value as i64)
}

fn sample_depths(fields: &[&str]) -> Vec<i64> {
    if fields.len() <= 9 {
        return Vec::new();
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let dp_index = format_keys.iter().position(|key| *key == "DP");
    let ad_index = format_keys.iter().position(|key| *key == "AD");
    let mut depths = Vec::new();
    for sample in &fields[9..] {
        let values: Vec<&str> = sample.split(':').collect();
        if let Some(i) = dp_index
            && let Some(dp) = values.get(i).and_then(|value| parse_sample_i64(value))
        {
            depths.push(dp);
            continue;
        }
        if let Some(i) = ad_index
            && let Some(ad) = values.get(i).and_then(|value| parse_ad_depth(value))
        {
            depths.push(ad);
        }
    }
    depths
}

fn sample_depth(fields: &[&str], sample_index: usize) -> Option<i64> {
    if fields.len() <= 9 {
        return None;
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let values: Vec<&str> = fields.get(9 + sample_index)?.split(':').collect();
    if let Some(i) = format_keys.iter().position(|key| *key == "DP")
        && let Some(dp) = values.get(i).and_then(|value| parse_sample_i64(value))
    {
        return Some(dp);
    }
    if let Some(i) = format_keys.iter().position(|key| *key == "AD") {
        return values.get(i).and_then(|value| parse_ad_depth(value));
    }
    None
}

fn parse_gt_alleles(gt: &str) -> Vec<Option<usize>> {
    gt.split(['/', '|'])
        .map(|allele| {
            if allele == "." || allele.is_empty() {
                None
            } else {
                allele.parse::<usize>().ok()
            }
        })
        .collect()
}

fn parse_sample_i64(raw: &str) -> Option<i64> {
    let first = raw.split(',').next().unwrap_or(raw);
    if first == "." || first.is_empty() {
        None
    } else {
        first.parse::<i64>().ok()
    }
}

fn parse_ad_depth(raw: &str) -> Option<i64> {
    let mut sum = 0;
    let mut has_value = false;
    for value in raw.split(',') {
        if value == "." || value.is_empty() {
            continue;
        }
        sum += value.parse::<i64>().ok()?;
        has_value = true;
    }
    has_value.then_some(sum)
}

fn depth_bin_index(spec: &DepthSpec, value: i64) -> usize {
    if value < spec.min {
        0
    } else if value > spec.max {
        spec.nbins - 1
    } else {
        1 + ((value - spec.min) / spec.step) as usize
    }
}

fn depth_bin_label(spec: &DepthSpec, index: usize) -> String {
    if index == 0 {
        format!("<{}", spec.min)
    } else if index + 1 == spec.nbins {
        format!(">{}", spec.max)
    } else {
        (spec.min + index as i64 - 1).to_string()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantKind {
    Snp,
    Mnp,
    Indel,
    Other,
}

fn classify(reference: &str, alt: &str) -> VariantKind {
    if alt == "." || alt == reference {
        return VariantKind::Other;
    }
    if alt.starts_with('<') || alt.contains('[') || alt.contains(']') {
        return VariantKind::Other;
    }
    if reference.len() == 1 && alt.len() == 1 {
        if reference.chars().all(is_dna) && alt.chars().all(is_dna) {
            return VariantKind::Snp;
        }
        return VariantKind::Other;
    }
    if reference.len() == alt.len() {
        return VariantKind::Mnp;
    }
    VariantKind::Indel
}

fn is_dna(b: char) -> bool {
    matches!(b, 'A' | 'C' | 'G' | 'T' | 'N' | 'a' | 'c' | 'g' | 't' | 'n')
}

fn ts_tv_count(reference: &str, alt: &str) -> Option<(u64, u64)> {
    let r = reference.chars().next()?.to_ascii_uppercase();
    let a = alt.chars().next()?.to_ascii_uppercase();
    let pair = (r, a);
    let is_transition = matches!(pair, ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C'));
    if is_transition {
        Some((1, 0))
    } else if matches!(r, 'A' | 'C' | 'G' | 'T') && matches!(a, 'A' | 'C' | 'G' | 'T') {
        Some((0, 1))
    } else {
        None
    }
}

fn record_in_regions(fields: &[&str], regions: &[RegionSpec], targets: &[RegionSpec]) -> bool {
    let chrom = fields[0];
    let pos = fields[1].parse::<i64>().unwrap_or(-1);
    let in_regions = regions.is_empty() || matches_any(regions, chrom, pos);
    let in_targets = targets.is_empty() || matches_any(targets, chrom, pos);
    in_regions && in_targets
}

fn matches_any(specs: &[RegionSpec], chrom: &str, pos: i64) -> bool {
    specs.iter().any(|spec| {
        spec.contig == chrom
            && spec.start.map(|s| pos >= s).unwrap_or(true)
            && spec.end.map(|e| pos <= e).unwrap_or(true)
    })
}

fn apply_filters_pass(fields: &[&str], filters: Option<&[String]>) -> bool {
    let Some(filters) = filters else {
        return true;
    };
    let f = fields[6];
    if filters.iter().any(|t| t == ".") && f == "." {
        return true;
    }
    let parts: Vec<&str> = if f == "." {
        Vec::new()
    } else {
        f.split(';').collect()
    };
    parts.iter().any(|p| filters.iter().any(|t| t == p))
}

fn expression_pass(fields: &[&str], args: &Args) -> io::Result<bool> {
    if let Some(expr) = &args.include_expr
        && !evaluate_expression(expr, fields)?.truthy()
    {
        return Ok(false);
    }
    if let Some(expr) = &args.exclude_expr
        && evaluate_expression(expr, fields)?.truthy()
    {
        return Ok(false);
    }
    Ok(true)
}

fn evaluate_expression(expr: &str, fields: &[&str]) -> io::Result<FilterValue> {
    let context = EvalContext::new();
    let owned: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
    bcffilter::eval_expression_with(expr, &context, |name, sample_index| {
        if sample_index.is_some() {
            return None;
        }
        stats_record_lookup(name, &owned)
    })
}

fn stats_record_lookup(token: &str, fields: &[String]) -> Option<FilterValue> {
    if token.eq_ignore_ascii_case("TYPE") {
        return Some(FilterValue::String(record_type_label(fields)));
    }
    if let Some(info_token) = token
        .strip_prefix("INFO/")
        .or_else(|| token.strip_prefix("Info/"))
        .or_else(|| token.strip_prefix("info/"))
    {
        return super::filter::record_lookup(info_token, fields);
    }
    super::filter::record_lookup(token, fields)
}

fn record_type_label(fields: &[String]) -> String {
    let reference = fields.get(3).map(String::as_str).unwrap_or(".");
    let alt = fields.get(4).map(String::as_str).unwrap_or(".");
    let mut saw_snp = false;
    let mut saw_mnp = false;
    let mut saw_indel = false;
    let mut saw_other = false;

    for allele in alt.split(',').filter(|allele| !allele.is_empty()) {
        match classify(reference, allele) {
            VariantKind::Snp => saw_snp = true,
            VariantKind::Mnp => saw_mnp = true,
            VariantKind::Indel => saw_indel = true,
            VariantKind::Other => saw_other = true,
        }
    }

    let mut labels = Vec::new();
    if saw_snp {
        labels.push("snp");
    }
    if saw_mnp {
        labels.push("mnp");
    }
    if saw_indel {
        labels.push("indel");
    }
    if saw_other {
        labels.push("other");
    }

    if labels.is_empty() {
        "ref".into()
    } else {
        labels.join(",")
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

fn body_after_header(text: &str) -> &str {
    let header_len = extract_header_len(text);
    &text[header_len..]
}

fn extract_header_len(text: &str) -> usize {
    let mut len = 0;
    for line in text.split_inclusive('\n') {
        if !line.starts_with('#') {
            break;
        }
        len += line.len();
        if line.starts_with("#CHROM\t") {
            break;
        }
    }
    len
}

fn sample_names_from_header(text: &str) -> Vec<String> {
    for line in text.split_inclusive('\n') {
        if line.starts_with("#CHROM\t") {
            let fields: Vec<&str> = line.trim_end().split('\t').collect();
            return if fields.len() > 9 {
                fields[9..].iter().map(|name| (*name).to_owned()).collect()
            } else {
                Vec::new()
            };
        }
        if !line.starts_with('#') {
            break;
        }
    }
    Vec::new()
}

fn select_samples(header_samples: &[String], args: &Args) -> io::Result<SelectedSamples> {
    let sample_list = match args.sample_list.as_deref() {
        Some("-") | None => None,
        other => other,
    };
    let selected = crate::smpl_ilist::init(
        header_samples,
        sample_list,
        args.samples_is_file,
        crate::smpl_ilist::SMPL_STRICT,
    )?;
    let names = selected
        .idx
        .iter()
        .filter_map(|&idx| header_samples.get(idx).cloned())
        .collect();
    Ok(SelectedSamples {
        names,
        indices: selected.idx,
    })
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
    fn classify_distinguishes_kinds() {
        assert_eq!(classify("A", "C"), VariantKind::Snp);
        assert_eq!(classify("AT", "CG"), VariantKind::Mnp);
        assert_eq!(classify("A", "AT"), VariantKind::Indel);
        assert_eq!(classify("A", "<INS>"), VariantKind::Other);
        assert_eq!(classify("A", "."), VariantKind::Other);
    }

    #[test]
    fn ts_tv_recognizes_transitions() {
        assert_eq!(ts_tv_count("A", "G"), Some((1, 0)));
        assert_eq!(ts_tv_count("C", "T"), Some((1, 0)));
        assert_eq!(ts_tv_count("A", "C"), Some((0, 1)));
        assert_eq!(ts_tv_count("G", "T"), Some((0, 1)));
    }

    #[test]
    fn af_bin_index_picks_first_le() {
        let bins = vec![0.1, 0.5, 1.0];
        assert_eq!(af_bin_index(0.05, &bins), 0);
        assert_eq!(af_bin_index(0.1, &bins), 0);
        assert_eq!(af_bin_index(0.3, &bins), 1);
        assert_eq!(af_bin_index(0.99, &bins), 2);
        assert_eq!(af_bin_index(1.5, &bins), 3);
    }

    #[test]
    fn depth_bins_match_upstream_bucket_edges() {
        let spec = parse_depth("10,20,5").unwrap();
        assert_eq!(depth_bin_index(&spec, 9), 0);
        assert_eq!(depth_bin_index(&spec, 10), 1);
        assert_eq!(depth_bin_index(&spec, 14), 1);
        assert_eq!(depth_bin_index(&spec, 15), 2);
        assert_eq!(depth_bin_index(&spec, 21), spec.nbins - 1);
        assert_eq!(depth_bin_label(&spec, 0), "<10");
        assert_eq!(depth_bin_label(&spec, spec.nbins - 1), ">20");
    }
}
