//! Focused `bcftools norm` implementation (upstream `vcfnorm.c`).
//!
//! This first local slice supports duplicate-record removal with
//! `-d/--rm-dup`, a narrow `-c s` reference-swap path, and simple
//! multiallelic splitting. Full left alignment, atomization, join
//! multiallelics, and INFO/FORMAT projection remain tracked in `TODO.md`.

use std::collections::HashSet;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::diagnostics::fmt_etag;
use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::reference::FastaReference;
use crate::vcf_compat::normalize_vcf_text;

const USAGE: &str = "\n\
About:   Left-align and normalize indels; this local slice supports duplicate removal.\n\
Usage:   bcftools norm [OPTIONS] <in.vcf.gz>\n\
\n\
Options:\n\
    -d, --rm-dup TYPE              Remove duplicate records: snps|indels|both|all|exact|none|any\n\
    -c, --check-ref MODE           Reference check mode; this slice supports 's' swap\n\
    -f, --fasta-ref FILE           Accepted by the duplicate-removal slice for command compatibility\n\
    -i, --include EXPR             Include records matching EXPR in duplicate-removal decisions\n\
    -m, --multiallelics MODE       Split multiallelic records; this slice supports '-'\n\
    -S, --sort MODE                Sort split records; this slice supports 'lex'\n\
    -o, --output FILE              Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9] u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
        --no-version               Accepted for command-shape compatibility\n\
\n";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    fasta_ref: Option<PathBuf>,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    rm_dup: DupMode,
    include_expr: Option<String>,
    check_ref: Option<CheckRefMode>,
    split_multiallelic: bool,
    split_sort: SplitSortMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DupMode {
    None,
    Exact,
    Snps,
    Indels,
    Both,
    All,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CheckRefMode {
    Swap,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SplitSortMode {
    Input,
    Lex,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum VariantKind {
    Snp,
    Indel,
    Other,
}

#[derive(Debug, Hash, PartialEq, Eq)]
struct DupKey {
    chrom: String,
    pos: String,
    rest: String,
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
                eprintln!("{}", fmt_etag("main_vcfnorm", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfnorm", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut input = None;
    let mut fasta_ref = None;
    let mut output = None;
    let mut output_kind = OutputKind::VcfText;
    let mut rm_dup = DupMode::None;
    let mut include_expr = None;
    let mut check_ref = None;
    let mut split_multiallelic = false;
    let mut split_sort = SplitSortMode::Input;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-d" | "--rm-dup" | "--rm-dups" => {
                rm_dup = parse_dup_mode(&next_string(&mut iter, raw.as_ref())?)?
            }
            "-c" | "--check-ref" => {
                check_ref = Some(parse_check_ref_mode(&next_string(
                    &mut iter,
                    raw.as_ref(),
                )?)?)
            }
            "-f" | "--fasta-ref" => {
                fasta_ref = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?));
            }
            "-i" | "--include" => include_expr = Some(next_string(&mut iter, raw.as_ref())?),
            "-m" | "--multiallelics" => {
                split_multiallelic =
                    parse_multiallelic_mode(&next_string(&mut iter, raw.as_ref())?)?
            }
            "-S" | "--sort" => {
                split_sort = parse_split_sort_mode(&next_string(&mut iter, raw.as_ref())?)?
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--no-version" => {}
            _ if raw.starts_with("--rm-dup=") || raw.starts_with("--rm-dups=") => {
                rm_dup = parse_dup_mode(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--check-ref=") => {
                check_ref = Some(parse_check_ref_mode(value_after_equals(&raw))?);
            }
            _ if raw.starts_with("--fasta-ref=") => {
                fasta_ref = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--include=") => {
                include_expr = Some(value_after_equals(&raw).to_owned());
            }
            _ if raw.starts_with("--multiallelics=") => {
                split_multiallelic = parse_multiallelic_mode(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--sort=") => {
                split_sort = parse_split_sort_mode(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-f") && raw.len() > 2 => {
                fasta_ref = Some(PathBuf::from(&raw[2..]));
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-c") && raw.len() > 2 => {
                check_ref = Some(parse_check_ref_mode(&raw[2..])?)
            }
            _ if raw.starts_with("-d") && raw.len() > 2 => rm_dup = parse_dup_mode(&raw[2..])?,
            _ if raw.starts_with("-i") && raw.len() > 2 => include_expr = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-m") && raw.len() > 2 => {
                split_multiallelic = parse_multiallelic_mode(&raw[2..])?
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with("-S") && raw.len() > 2 => {
                split_sort = parse_split_sort_mode(&raw[2..])?
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!(
                    "unrecognized option '{raw}' in this local norm slice"
                )));
            }
            _ if input.is_none() => input = Some(PathBuf::from(raw.as_ref())),
            _ => {
                return Err(ParseOutcome::Error(format!(
                    "unexpected extra input '{raw}'"
                )));
            }
        }
    }

    let input =
        input.ok_or_else(|| ParseOutcome::Error("expected one input VCF/BCF path".into()))?;
    Ok(Args {
        input,
        fasta_ref,
        output,
        output_kind,
        rm_dup,
        include_expr,
        check_ref,
        split_multiallelic,
        split_sort,
    })
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

fn parse_dup_mode(raw: &str) -> Result<DupMode, ParseOutcome> {
    match raw {
        "none" | "exact" => Ok(DupMode::Exact),
        "snps" => Ok(DupMode::Snps),
        "indels" => Ok(DupMode::Indels),
        "both" | "any" => Ok(DupMode::Both),
        "all" => Ok(DupMode::All),
        _ => Err(ParseOutcome::Error(format!(
            "unknown duplicate mode '{raw}'"
        ))),
    }
}

fn parse_check_ref_mode(raw: &str) -> Result<CheckRefMode, ParseOutcome> {
    match raw {
        "s" => Ok(CheckRefMode::Swap),
        _ => Err(ParseOutcome::Error(format!(
            "unsupported check-ref mode '{raw}' in this local norm slice"
        ))),
    }
}

fn parse_multiallelic_mode(raw: &str) -> Result<bool, ParseOutcome> {
    match raw {
        "-" => Ok(true),
        _ => Err(ParseOutcome::Error(format!(
            "unsupported multiallelic mode '{raw}' in this local norm slice"
        ))),
    }
}

fn parse_split_sort_mode(raw: &str) -> Result<SplitSortMode, ParseOutcome> {
    match raw {
        "lex" => Ok(SplitSortMode::Lex),
        _ => Err(ParseOutcome::Error(format!(
            "unsupported split sort mode '{raw}' in this local norm slice"
        ))),
    }
}

fn run(args: &Args) -> io::Result<()> {
    let input = read_vcf_text(&args.input)?;
    let reference = args
        .fasta_ref
        .as_ref()
        .map(FastaReference::open)
        .transpose()?;
    let output = if let Some(CheckRefMode::Swap) = args.check_ref {
        let reference = reference.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "-c s requires -f/--fasta-ref in this local norm slice",
            )
        })?;
        check_reference_swap(&input, reference)?
    } else if args.split_multiallelic {
        split_multiallelic_records(&input, args.split_sort)?
    } else {
        let include_filter = args
            .include_expr
            .as_deref()
            .map(IncludeFilter::from_expr)
            .transpose()?;
        remove_duplicates(
            &input,
            args.rm_dup,
            reference.as_ref(),
            include_filter.as_ref(),
        )?
    };
    write_output(output.as_bytes(), args)
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
        ".bcftools-rs-norm-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn remove_duplicates(
    input: &str,
    mode: DupMode,
    reference: Option<&FastaReference>,
    include_filter: Option<&IncludeFilter>,
) -> io::Result<String> {
    let mut out = String::with_capacity(input.len());
    let mut seen = HashSet::new();
    let filter_header_seen = input
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    let mut filter_header_inserted = false;

    for line in input.lines() {
        if line.starts_with("##fileformat=") {
            out.push_str(line);
            out.push('\n');
            if !filter_header_seen && !filter_header_inserted {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
                filter_header_inserted = true;
            }
            continue;
        }
        if line.starts_with("#CHROM") {
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
        let mut fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if fields.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid VCF record with fewer than 8 columns: {line}"),
            ));
        }
        if let Some(reference) = reference {
            normalize_duplicate_record(&mut fields, reference)?;
        }
        if include_filter
            .map(|filter| filter.matches(&fields))
            .transpose()?
            == Some(false)
        {
            out.push_str(&fields.join("\t"));
            out.push('\n');
            continue;
        }
        let field_refs: Vec<&str> = fields.iter().map(String::as_str).collect();
        let keys = duplicate_keys(&field_refs, mode);
        let is_dup = keys.iter().any(|key| seen.contains(key));
        if !is_dup {
            seen.extend(keys);
            out.push_str(&fields.join("\t"));
            out.push('\n');
        }
    }

    Ok(out)
}

fn split_multiallelic_records(input: &str, sort_mode: SplitSortMode) -> io::Result<String> {
    let mut out = String::with_capacity(input.len());
    let filter_header_seen = input
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    let mut filter_header_inserted = false;

    for line in input.lines() {
        if line.starts_with("##fileformat=") {
            out.push_str(line);
            out.push('\n');
            if !filter_header_seen && !filter_header_inserted {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
                filter_header_inserted = true;
            }
            continue;
        }
        if line.starts_with("#CHROM") {
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

        let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if fields.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid VCF record with fewer than 8 columns: {line}"),
            ));
        }
        let alts: Vec<&str> = fields[4].split(',').collect();
        if alts.len() <= 1 {
            out.push_str(line);
            out.push('\n');
            continue;
        }

        let mut split_rows: Vec<Vec<String>> = alts
            .into_iter()
            .map(|alt| {
                let mut row = fields.clone();
                row[4] = alt.to_owned();
                row
            })
            .collect();
        if sort_mode == SplitSortMode::Lex {
            split_rows.sort_by(|left, right| left[4].cmp(&right[4]));
        }
        for row in split_rows {
            out.push_str(&row.join("\t"));
            out.push('\n');
        }
    }

    Ok(out)
}

#[derive(Debug)]
enum IncludeFilter {
    Expression(String),
    FileMembership {
        token: String,
        negate: bool,
        values: HashSet<String>,
    },
}

impl IncludeFilter {
    fn from_expr(raw: &str) -> io::Result<Self> {
        if let Some((token, negate, path)) = parse_file_membership(raw) {
            let text = fs::read_to_string(path)?;
            let values = text
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(ToOwned::to_owned)
                .collect();
            return Ok(Self::FileMembership {
                token,
                negate,
                values,
            });
        }
        Ok(Self::Expression(raw.to_owned()))
    }

    fn matches(&self, fields: &[String]) -> io::Result<bool> {
        match self {
            Self::Expression(expr) => evaluate_expression(expr, fields),
            Self::FileMembership {
                token,
                negate,
                values,
            } => {
                let matched = super::filter::record_lookup(token, fields)
                    .is_some_and(|value| value_matches_any(&value, values));
                Ok(matched != *negate)
            }
        }
    }
}

fn parse_file_membership(raw: &str) -> Option<(String, bool, &str)> {
    for (needle, negate) in [("!=", true), ("==", false), ("=", false)] {
        let Some((lhs, rhs)) = raw.split_once(needle) else {
            continue;
        };
        let path = rhs.trim().strip_prefix('@')?;
        return Some((
            lhs.trim().to_owned(),
            negate,
            path.trim_matches(|c| c == '"' || c == '\''),
        ));
    }
    None
}

fn evaluate_expression(expr: &str, fields: &[String]) -> io::Result<bool> {
    bcffilter::eval_expression_with(expr, &EvalContext::new(), |name, sample_index| {
        if sample_index.is_some() {
            return None;
        }
        super::filter::record_lookup(name, fields)
    })
    .map(|value| value.truthy())
}

fn value_matches_any(value: &FilterValue, values: &HashSet<String>) -> bool {
    match value {
        FilterValue::Missing => values.contains("."),
        FilterValue::Bool(value) => values.contains(if *value { "true" } else { "false" }),
        FilterValue::Number(value) => {
            if value.fract() == 0.0 {
                values.contains(&format!("{value:.0}"))
            } else {
                values.contains(&value.to_string())
            }
        }
        FilterValue::String(value) => values.contains(value),
        FilterValue::List(list) => list.iter().any(|value| value_matches_any(value, values)),
    }
}

fn check_reference_swap(input: &str, reference: &FastaReference) -> io::Result<String> {
    let mut out = String::with_capacity(input.len());
    let filter_header_seen = input
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    let mut filter_header_inserted = false;

    for line in input.lines() {
        if line.starts_with("##fileformat=") {
            out.push_str(line);
            out.push('\n');
            if !filter_header_seen && !filter_header_inserted {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
                filter_header_inserted = true;
            }
            continue;
        }
        if line.starts_with("#CHROM") {
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

        let mut fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if fields.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid VCF record with fewer than 8 columns: {line}"),
            ));
        }
        swap_ref_alt_if_needed(&mut fields, reference)?;
        out.push_str(&fields.join("\t"));
        out.push('\n');
    }

    Ok(out)
}

fn swap_ref_alt_if_needed(fields: &mut [String], reference: &FastaReference) -> io::Result<()> {
    if fields[3].len() != 1 {
        return Ok(());
    }
    let pos = parse_pos(&fields[1])?;
    let Some(reference_base) = fetch_base(reference, &fields[0], pos)? else {
        return Ok(());
    };
    if fields[3].as_bytes()[0].to_ascii_uppercase() == reference_base {
        return Ok(());
    }
    let Some(alt_index) = fields[4]
        .split(',')
        .position(|alt| alt.len() == 1 && alt.as_bytes()[0].to_ascii_uppercase() == reference_base)
    else {
        return Ok(());
    };

    let old_ref = fields[3].clone();
    let mut alts: Vec<String> = fields[4].split(',').map(str::to_owned).collect();
    alts[alt_index] = old_ref;
    fields[3] = char::from(reference_base).to_string();
    fields[4] = alts.join(",");
    swap_genotypes(fields, alt_index + 1);
    Ok(())
}

fn swap_genotypes(fields: &mut [String], alt_allele: usize) {
    if fields.len() <= 9 {
        return;
    }
    let gt_index = fields[8].split(':').position(|key| key == "GT");
    let Some(gt_index) = gt_index else {
        return;
    };
    for sample in &mut fields[9..] {
        let mut parts: Vec<String> = sample.split(':').map(str::to_owned).collect();
        if let Some(gt) = parts.get_mut(gt_index) {
            *gt = swap_gt_alleles(gt, alt_allele);
            *sample = parts.join(":");
        }
    }
}

fn swap_gt_alleles(gt: &str, alt_allele: usize) -> String {
    let mut out = String::with_capacity(gt.len());
    let mut allele = String::new();
    for ch in gt.chars() {
        if ch == '/' || ch == '|' {
            out.push_str(&swapped_allele(&allele, alt_allele));
            allele.clear();
            out.push(ch);
        } else {
            allele.push(ch);
        }
    }
    out.push_str(&swapped_allele(&allele, alt_allele));
    out
}

fn swapped_allele(raw: &str, alt_allele: usize) -> String {
    match raw.parse::<usize>() {
        Ok(0) => alt_allele.to_string(),
        Ok(value) if value == alt_allele => "0".to_owned(),
        _ => raw.to_owned(),
    }
}

fn normalize_duplicate_record(fields: &mut [String], reference: &FastaReference) -> io::Result<()> {
    if fields.len() < 8 || fields[4].contains(',') {
        return Ok(());
    }

    if fields[4] == "<DEL>" {
        normalize_symbolic_deletion(fields, reference)
    } else {
        normalize_simple_deletion(fields, reference)
    }
}

fn normalize_symbolic_deletion(
    fields: &mut [String],
    reference: &FastaReference,
) -> io::Result<()> {
    if !info_has_negative_svlen(&fields[7]) {
        return Ok(());
    }
    let pos = parse_pos(&fields[1])?;
    let Some(first_deleted) = fetch_base(reference, &fields[0], pos)? else {
        return Ok(());
    };
    let deletion_start = leftmost_repeat_start(reference, &fields[0], pos, first_deleted)?;
    if deletion_start <= 1 || deletion_start == pos {
        return Ok(());
    }
    let anchor_pos = deletion_start - 1;
    if let Some(anchor) = fetch_base(reference, &fields[0], anchor_pos)? {
        fields[1] = anchor_pos.to_string();
        fields[3] = char::from(anchor).to_string();
    }
    Ok(())
}

fn normalize_simple_deletion(fields: &mut [String], reference: &FastaReference) -> io::Result<()> {
    let reference_allele = fields[3].as_bytes();
    let alternate = fields[4].as_bytes();
    if reference_allele.len() <= alternate.len() || !reference_allele.starts_with(alternate) {
        return Ok(());
    }
    let deleted = &reference_allele[alternate.len()..];
    if deleted.len() != 1 {
        return Ok(());
    }
    let deleted_base = deleted[0];

    let pos = parse_pos(&fields[1])?;
    let deletion_start = pos + alternate.len();
    let deletion_start = leftmost_repeat_start(
        reference,
        &fields[0],
        deletion_start,
        deleted_base.to_ascii_uppercase(),
    )?;
    if deletion_start <= 1 {
        return Ok(());
    }
    let anchor_pos = deletion_start - 1;
    if let Some(anchor) = fetch_base(reference, &fields[0], anchor_pos)? {
        fields[1] = anchor_pos.to_string();
        fields[3] = format!("{}{}", char::from(anchor), char::from(deleted_base));
        fields[4] = char::from(anchor).to_string();
    }
    Ok(())
}

fn info_has_negative_svlen(info: &str) -> bool {
    info.split(';').any(|entry| {
        entry
            .strip_prefix("SVLEN=")
            .and_then(|value| value.split(',').next())
            .is_some_and(|value| value.starts_with('-'))
    })
}

fn parse_pos(raw: &str) -> io::Result<usize> {
    raw.parse::<usize>()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, format!("invalid POS '{raw}'")))
}

fn leftmost_repeat_start(
    reference: &FastaReference,
    chrom: &str,
    mut pos: usize,
    base: u8,
) -> io::Result<usize> {
    while pos > 1 && fetch_base(reference, chrom, pos - 1)? == Some(base.to_ascii_uppercase()) {
        pos -= 1;
    }
    Ok(pos)
}

fn fetch_base(reference: &FastaReference, chrom: &str, pos: usize) -> io::Result<Option<u8>> {
    let region = format!("{chrom}:{pos}-{pos}");
    match reference.fetch_region(&region) {
        Ok(sequence) => Ok(sequence.first().map(|base| base.to_ascii_uppercase())),
        Err(error) if error.kind() == io::ErrorKind::NotFound => Ok(None),
        Err(error) => Err(error),
    }
}

fn duplicate_keys(fields: &[&str], mode: DupMode) -> Vec<DupKey> {
    let mut keys = Vec::new();
    match mode {
        DupMode::All => keys.push(dup_key(fields, "")),
        DupMode::Exact | DupMode::None => keys.push(dup_key(fields, &exact_rest(fields))),
        DupMode::Snps => {
            if allele_class(fields).has_snp {
                keys.push(dup_key(fields, "Snp"));
            }
        }
        DupMode::Indels => {
            if allele_class(fields).has_indel {
                keys.push(dup_key(fields, "Indel"));
            }
        }
        DupMode::Both => {
            let class = allele_class(fields);
            if class.has_snp {
                keys.push(dup_key(fields, "Snp"));
            }
            if class.has_indel {
                keys.push(dup_key(fields, "Indel"));
            }
        }
    }
    keys
}

fn dup_key(fields: &[&str], rest: &str) -> DupKey {
    DupKey {
        chrom: fields[0].to_owned(),
        pos: fields[1].to_owned(),
        rest: rest.to_owned(),
    }
}

fn exact_rest(fields: &[&str]) -> String {
    let mut alts: Vec<&str> = fields[4].split(',').collect();
    alts.sort_unstable();
    if alts.iter().any(|alt| alt.starts_with('<')) {
        format!("{}:{}:{}", fields[3], alts.join(","), fields[7])
    } else {
        format!("{}:{}", fields[3], alts.join(","))
    }
}

#[derive(Clone, Copy)]
struct AlleleClass {
    has_snp: bool,
    has_indel: bool,
}

fn allele_class(fields: &[&str]) -> AlleleClass {
    let alts: Vec<&str> = fields[4].split(',').collect();
    let kind = classify_kind(fields[3], &alts);
    AlleleClass {
        has_snp: matches!(kind, VariantKind::Snp | VariantKind::Other)
            && alts
                .iter()
                .any(|alt| alt.len() == 1 && fields[3].len() == 1 && *alt != "." && *alt != "*"),
        has_indel: matches!(kind, VariantKind::Indel | VariantKind::Other)
            && alts
                .iter()
                .any(|alt| alt.len() != fields[3].len() && *alt != "." && *alt != "*"),
    }
}

fn classify_kind(reference: &str, alts: &[&str]) -> VariantKind {
    let mut has_snp = false;
    let mut has_indel = false;
    for alt in alts {
        if *alt == "." || *alt == "*" {
            continue;
        }
        if alt.len() == 1 && reference.len() == 1 {
            has_snp = true;
        } else if alt.len() != reference.len() {
            has_indel = true;
        }
    }
    if has_snp && !has_indel {
        VariantKind::Snp
    } else if has_indel && !has_snp {
        VariantKind::Indel
    } else {
        VariantKind::Other
    }
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
    fn classify_kind_distinguishes_snp_and_indel() {
        assert_eq!(classify_kind("A", &["C"]), VariantKind::Snp);
        assert_eq!(classify_kind("A", &["AT"]), VariantKind::Indel);
        assert_eq!(classify_kind("AT", &["A"]), VariantKind::Indel);
        assert_eq!(classify_kind("A", &["C", "AT"]), VariantKind::Other);
    }

    #[test]
    fn rmdup_all_keeps_first_record_per_position() {
        let input = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t2\t.\tA\tC\t.\t.\t.\n\
1\t2\t.\tA\tG\t.\t.\t.\n";
        let out = remove_duplicates(input, DupMode::All, None, None).unwrap();
        assert!(out.contains("1\t2\t.\tA\tC\t.\t.\t."));
        assert!(!out.contains("1\t2\t.\tA\tG\t.\t.\t."));
    }
}
