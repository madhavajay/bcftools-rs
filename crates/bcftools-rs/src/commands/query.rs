//! Partial port of `bcftools query` (upstream `vcfquery.c`).
//!
//! This lands `-l/--list-samples` plus a small text-backed `-f` formatter for
//! core fields and simple sample loops. The full upstream formatter still
//! depends on the Phase 1 `convert` formatter.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Read as _, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::expr::{Filter, Value};
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};

const USAGE: &str = "\n\
About:   Extract fields from VCF/BCF files and print sample lists.\n\
Usage:   bcftools query [OPTIONS] <in.vcf.gz>|<in.bcf>\n\
\n\
Options:\n\
    -f, --format STR                 format string\n\
    -H, --print-header               print output header, -HH omits column indices\n\
    -i, --include EXPR               include only records matching expression\n\
    -e, --exclude EXPR               exclude records matching expression\n\
    -l, --list-samples               print sample names and exit\n\
    -r, --regions LIST               comma-separated regions\n\
    -R, --regions-file FILE          restrict to regions in FILE\n\
    -s, --samples LIST               comma-separated sample list\n\
    -S, --samples-file FILE          file of samples, optionally prefixed with ^\n\
    -t, --targets LIST               comma-separated targets, optionally prefixed with ^\n\
    -T, --targets-file FILE          restrict to targets in FILE, optionally prefixed with ^\n\
\n";

#[derive(Debug, Clone, PartialEq, Eq)]
struct Args {
    list_samples: bool,
    format: Option<String>,
    header_level: u8,
    samples: Option<String>,
    samples_is_file: bool,
    regions: Option<RegionFilterSpec>,
    filter: Option<FilterSpec>,
    input: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RegionFilterSpec {
    raw: String,
    is_file: bool,
    exclude: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FilterSpec {
    raw: String,
    exclude: bool,
}

/// Subcommand entry point. `argv[0]` is `"query"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args, io::stdout().lock()) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfquery", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ParseOutcome {
    Usage,
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let long_opts = [
        LongOpt::new("format", HasArg::Required, b'f' as i32),
        LongOpt::new("print-header", HasArg::None, b'H' as i32),
        LongOpt::new("include", HasArg::Required, b'i' as i32),
        LongOpt::new("exclude", HasArg::Required, b'e' as i32),
        LongOpt::new("list-samples", HasArg::None, b'l' as i32),
        LongOpt::new("regions", HasArg::Required, b'r' as i32),
        LongOpt::new("regions-file", HasArg::Required, b'R' as i32),
        LongOpt::new("samples", HasArg::Required, b's' as i32),
        LongOpt::new("samples-file", HasArg::Required, b'S' as i32),
        LongOpt::new("targets", HasArg::Required, b't' as i32),
        LongOpt::new("targets-file", HasArg::Required, b'T' as i32),
    ];

    let mut list_samples = false;
    let mut format = None;
    let mut header_level = 0u8;
    let mut samples = None;
    let mut samples_is_file = false;
    let mut regions = None;
    let mut filter = None;

    let mut g = Getopt::new("e:f:Hi:lR:r:s:S:T:t:", &long_opts, argv);
    loop {
        match g.next() {
            Ok(Some(m)) => match m.code {
                v if v == b'l' as i32 => list_samples = true,
                v if v == b'e' as i32 => {
                    if let Some(value) = m.value {
                        filter = Some(FilterSpec {
                            raw: value,
                            exclude: true,
                        });
                    }
                }
                v if v == b'f' as i32 => format = m.value,
                v if v == b'H' as i32 => header_level = header_level.saturating_add(1),
                v if v == b'i' as i32 => {
                    if let Some(value) = m.value {
                        filter = Some(FilterSpec {
                            raw: value,
                            exclude: false,
                        });
                    }
                }
                v if v == b'r' as i32 => {
                    if let Some(value) = m.value {
                        regions = Some(RegionFilterSpec {
                            raw: value,
                            is_file: false,
                            exclude: false,
                        });
                    }
                }
                v if v == b'R' as i32 => {
                    if let Some(value) = m.value {
                        regions = Some(RegionFilterSpec {
                            raw: value,
                            is_file: true,
                            exclude: false,
                        });
                    }
                }
                v if v == b's' as i32 => {
                    samples = m.value;
                    samples_is_file = false;
                }
                v if v == b'S' as i32 => {
                    samples = m.value;
                    samples_is_file = true;
                }
                v if v == b't' as i32 => {
                    if let Some(value) = m.value {
                        let (exclude, raw) = strip_exclusion_prefix(value);
                        regions = Some(RegionFilterSpec {
                            raw,
                            is_file: false,
                            exclude,
                        });
                    }
                }
                v if v == b'T' as i32 => {
                    if let Some(value) = m.value {
                        let (exclude, raw) = strip_exclusion_prefix(value);
                        regions = Some(RegionFilterSpec {
                            raw,
                            is_file: true,
                            exclude,
                        });
                    }
                }
                _ => return Err(ParseOutcome::Usage),
            },
            Ok(None) => break,
            Err(_) => return Err(ParseOutcome::Usage),
        }
    }

    if !list_samples && format.is_none() {
        return Err(ParseOutcome::Usage);
    }

    let rest = g.rest();
    if rest.len() > 1 {
        return Err(ParseOutcome::Usage);
    }
    let input = rest
        .first()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "-".into());

    Ok(Args {
        list_samples,
        format,
        header_level,
        samples,
        samples_is_file,
        regions,
        filter,
        input,
    })
}

fn strip_exclusion_prefix(value: String) -> (bool, String) {
    value
        .strip_prefix('^')
        .map(|s| (true, s.to_string()))
        .unwrap_or((false, value))
}

fn run<W: Write>(args: &Args, mut out: W) -> io::Result<()> {
    let input = materialize_input(&args.input)?;
    if args.list_samples {
        for sample in sample_names_from_path(&input, args.samples.as_deref(), args.samples_is_file)?
        {
            writeln!(out, "{sample}")?;
        }
    }
    if let Some(format) = &args.format {
        let options = QueryFormatOptions {
            sample_list: args.samples.as_deref(),
            sample_list_is_file: args.samples_is_file,
            header_level: args.header_level,
            region_spec: args.regions.as_ref(),
            filter_spec: args.filter.as_ref(),
        };
        query_format_from_path(&input, format, &options, &mut out)?;
    }
    Ok(())
}

fn materialize_input(input: &str) -> io::Result<PathBuf> {
    if input != "-" {
        return Ok(PathBuf::from(input));
    }

    let tmp = stdin_tmp_path();
    let mut data = Vec::new();
    io::stdin().lock().read_to_end(&mut data)?;
    fs::write(&tmp, data)?;
    Ok(tmp)
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-query-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn sample_names_from_path<P>(
    path: P,
    sample_list: Option<&str>,
    sample_list_is_file: bool,
) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    let samples = header_sample_names_from_path(path)?;
    let selected = crate::smpl_ilist::init(
        &samples,
        sample_list,
        sample_list_is_file,
        crate::smpl_ilist::SMPL_STRICT,
    )?;
    Ok(selected
        .idx
        .into_iter()
        .map(|idx| samples[idx].clone())
        .collect())
}

fn header_sample_names_from_path<P>(path: P) -> io::Result<Vec<String>>
where
    P: AsRef<Path>,
{
    use htslib_rs::variant_io_compat::{
        read_bcf_header_from_path, read_vcf_header, read_vcf_header_from_path,
    };

    let path = path.as_ref();
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    let header = if fmt.exact == Exact::Bcf {
        read_bcf_header_from_path(path)?
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        read_vcf_header(BufReader::new(dec))?
    } else {
        read_vcf_header_from_path(path)?
    };

    Ok(header
        .sample_names()
        .iter()
        .map(ToString::to_string)
        .collect())
}

#[derive(Debug, Clone, Copy)]
struct QueryFormatOptions<'a> {
    sample_list: Option<&'a str>,
    sample_list_is_file: bool,
    header_level: u8,
    region_spec: Option<&'a RegionFilterSpec>,
    filter_spec: Option<&'a FilterSpec>,
}

fn query_format_from_path<W: Write>(
    path: &Path,
    format: &str,
    options: &QueryFormatOptions<'_>,
    out: &mut W,
) -> io::Result<()> {
    let text = vcf_text_from_path(path)?;
    query_format_text(text.as_str(), format, options, out)
}

fn vcf_text_from_path(path: &Path) -> io::Result<String> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let mut dec = flate2::read::MultiGzDecoder::new(f);
        let mut text = String::new();
        dec.read_to_string(&mut text)?;
        return Ok(text);
    }
    fs::read_to_string(path)
}

fn query_format_text<W: Write>(
    text: &str,
    format: &str,
    options: &QueryFormatOptions<'_>,
    out: &mut W,
) -> io::Result<()> {
    let mut samples = Vec::new();
    let mut sample_indices = Vec::new();
    let mut region_filter: Option<RegionFilter> = None;
    let query_filter = options
        .filter_spec
        .map(QueryFilter::from_spec)
        .transpose()?;
    let mut wrote_header = false;
    for line in text.lines() {
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with("#CHROM\t") {
            samples = line.split('\t').skip(9).map(ToOwned::to_owned).collect();
            sample_indices =
                query_sample_indices(&samples, options.sample_list, options.sample_list_is_file)?;
            region_filter = options
                .region_spec
                .map(RegionFilter::from_spec)
                .transpose()?;
            if options.header_level > 0 {
                out.write_all(
                    render_format_header(format, &samples, &sample_indices, options.header_level)
                        .as_bytes(),
                )?;
                wrote_header = true;
            }
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        if let Some(filter) = &region_filter
            && !filter.matches(line)?
        {
            continue;
        }
        let record = TextRecord::parse(line, &samples, &sample_indices);
        if let Some(filter) = &query_filter
            && !filter.matches(&record)?
        {
            continue;
        }
        let rendered = render_format(format, &record);
        out.write_all(rendered.as_bytes())?;
    }
    if options.header_level > 0 && !wrote_header {
        out.write_all(render_format_header(format, &[], &[], options.header_level).as_bytes())?;
    }
    Ok(())
}

#[derive(Debug, Clone)]
struct RegionFilter {
    regions: Vec<QueryRegion>,
    exclude: bool,
}

#[derive(Debug, Clone)]
struct QueryRegion {
    chrom: String,
    start: Option<i64>,
    end: Option<i64>,
}

impl RegionFilter {
    fn from_spec(spec: &RegionFilterSpec) -> io::Result<Self> {
        let regions = if spec.is_file {
            read_region_file(&spec.raw)?
        } else {
            parse_region_list(&spec.raw)?
        };
        Ok(Self {
            regions,
            exclude: spec.exclude,
        })
    }

    fn matches(&self, record_line: &str) -> io::Result<bool> {
        let mut fields = record_line.split('\t');
        let chrom = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing CHROM field"))?;
        let pos = fields
            .next()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing POS field"))?
            .parse::<i64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        let matched = self.regions.iter().any(|region| {
            region.chrom == chrom
                && region.start.is_none_or(|start| pos >= start)
                && region.end.is_none_or(|end| pos <= end)
        });
        Ok(matched != self.exclude)
    }
}

fn parse_region_list(raw: &str) -> io::Result<Vec<QueryRegion>> {
    raw.split(',')
        .filter(|item| !item.trim().is_empty())
        .map(|item| parse_region_item(item.trim()))
        .collect()
}

fn parse_region_item(raw: &str) -> io::Result<QueryRegion> {
    let (chrom, coordinates) = raw.split_once(':').unwrap_or((raw, ""));
    if chrom.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "empty region chromosome",
        ));
    }
    if coordinates.is_empty() {
        return Ok(QueryRegion {
            chrom: chrom.to_string(),
            start: None,
            end: None,
        });
    }
    let (start, end) = coordinates
        .split_once('-')
        .unwrap_or((coordinates, coordinates));
    let start = parse_region_position(start)?;
    let end = parse_region_position(end)?;
    Ok(QueryRegion {
        chrom: chrom.to_string(),
        start: Some(start),
        end: Some(end),
    })
}

fn parse_region_position(raw: &str) -> io::Result<i64> {
    raw.replace(',', "")
        .parse::<i64>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn read_region_file(path: &str) -> io::Result<Vec<QueryRegion>> {
    let text = if path.ends_with(".gz") {
        let file = File::open(path)?;
        let mut dec = flate2::read::MultiGzDecoder::new(file);
        let mut text = String::new();
        dec.read_to_string(&mut text)?;
        text
    } else {
        fs::read_to_string(path)?
    };
    let is_bed = Path::new(path)
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bed"));
    text.lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(|line| parse_region_file_line(line, is_bed))
        .collect()
}

fn parse_region_file_line(line: &str, is_bed: bool) -> io::Result<QueryRegion> {
    let fields = line.split_whitespace().collect::<Vec<_>>();
    if fields.len() < 3 {
        return parse_region_item(line.trim());
    }
    let chrom = fields[0].to_string();
    let raw_start = parse_region_position(fields[1])?;
    let raw_end = parse_region_position(fields[2])?;
    let (start, end) = if is_bed {
        (raw_start + 1, raw_end)
    } else {
        (raw_start, raw_end)
    };
    Ok(QueryRegion {
        chrom,
        start: Some(start),
        end: Some(end),
    })
}

#[derive(Debug, Clone)]
struct QueryFilter {
    kind: QueryFilterKind,
    exclude: bool,
}

#[derive(Debug, Clone)]
enum QueryFilterKind {
    Expr(Filter),
    FilterIdMatch { id: String, negate: bool },
    PredicateGroups(Vec<Vec<SimplePredicate>>),
}

#[derive(Debug, Clone)]
struct SimplePredicate {
    lhs: String,
    vector_any: bool,
    op: PredicateOp,
    rhs: String,
}

#[derive(Debug, Clone, Copy)]
enum PredicateOp {
    Eq,
    Ne,
    Regex,
    NotRegex,
}

impl QueryFilter {
    fn from_spec(spec: &FilterSpec) -> io::Result<Self> {
        let kind = parse_filter_id_match(&spec.raw)
            .map(|(id, negate)| QueryFilterKind::FilterIdMatch { id, negate })
            .or_else(|| {
                parse_simple_predicate_groups(&spec.raw).map(QueryFilterKind::PredicateGroups)
            })
            .unwrap_or_else(|| {
                QueryFilterKind::Expr(Filter::new(normalize_filter_expr(&spec.raw)))
            });
        Ok(Self {
            kind,
            exclude: spec.exclude,
        })
    }

    fn matches(&self, record: &TextRecord<'_>) -> io::Result<bool> {
        let matched = match &self.kind {
            QueryFilterKind::Expr(filter) => {
                let value = filter
                    .eval_with(|src| lookup_filter_symbol(src, record))
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
                value.truth()
            }
            QueryFilterKind::FilterIdMatch { id, negate } => {
                record.filter_has_id(id.as_str()) != *negate
            }
            QueryFilterKind::PredicateGroups(groups) => groups
                .iter()
                .any(|predicates| predicates.iter().all(|predicate| predicate.matches(record))),
        };
        Ok(matched != self.exclude)
    }
}

impl SimplePredicate {
    fn matches(&self, record: &TextRecord<'_>) -> bool {
        let values = record.filter_values(&self.lhs, self.vector_any);
        match self.op {
            PredicateOp::Eq => values.iter().any(|value| value == &self.rhs),
            PredicateOp::Ne => values.iter().all(|value| value != &self.rhs),
            PredicateOp::Regex => values
                .iter()
                .any(|value| regex::Regex::new(&self.rhs).is_ok_and(|re| re.is_match(value))),
            PredicateOp::NotRegex => values
                .iter()
                .all(|value| regex::Regex::new(&self.rhs).is_ok_and(|re| !re.is_match(value))),
        }
    }
}

fn parse_simple_predicate_groups(raw: &str) -> Option<Vec<Vec<SimplePredicate>>> {
    let groups = split_simple_or(raw)
        .into_iter()
        .map(|term| {
            split_simple_and(term)
                .into_iter()
                .map(parse_simple_predicate)
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()?;
    (!groups.is_empty() && groups.iter().all(|group| !group.is_empty())).then_some(groups)
}

fn split_simple_or(raw: &str) -> Vec<&str> {
    split_simple_binary(raw, b'|')
}

fn split_simple_and(raw: &str) -> Vec<&str> {
    split_simple_binary(raw, b'&')
}

fn split_simple_binary(raw: &str, delimiter: u8) -> Vec<&str> {
    let mut parts = Vec::new();
    let mut start = 0usize;
    let mut in_string = false;
    let bytes = raw.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        match bytes[i] {
            b'"' => {
                in_string = !in_string;
                i += 1;
            }
            ch if ch == delimiter && !in_string => {
                parts.push(raw[start..i].trim());
                i += usize::from(i + 1 < bytes.len() && bytes[i + 1] == delimiter) + 1;
                start = i;
            }
            _ => i += 1,
        }
    }
    parts.push(raw[start..].trim());
    parts
}

fn parse_simple_predicate(raw: &str) -> Option<SimplePredicate> {
    let (lhs, op, rhs) = split_simple_predicate(raw)?;
    let lhs = lhs.trim();
    let (lhs, vector_any) = lhs
        .strip_suffix("[*]")
        .map(|lhs| (lhs, true))
        .unwrap_or((lhs, false));
    let rhs = parse_quoted_rhs(rhs.trim())?;
    Some(SimplePredicate {
        lhs: lhs.trim().to_string(),
        vector_any,
        op,
        rhs: rhs.to_string(),
    })
}

fn split_simple_predicate(raw: &str) -> Option<(&str, PredicateOp, &str)> {
    for (needle, op) in [
        ("!~", PredicateOp::NotRegex),
        ("!=", PredicateOp::Ne),
        ("==", PredicateOp::Eq),
        ("=", PredicateOp::Eq),
        ("~", PredicateOp::Regex),
    ] {
        if let Some((lhs, rhs)) = raw.split_once(needle) {
            return Some((lhs, op, rhs));
        }
    }
    None
}

fn parse_quoted_rhs(raw: &str) -> Option<&str> {
    let rest = raw.strip_prefix('"')?;
    let end = rest.find('"')?;
    rest[end + 1..].trim().is_empty().then_some(&rest[..end])
}

fn parse_filter_id_match(raw: &str) -> Option<(String, bool)> {
    let raw = raw.trim();
    let (lhs, op, rhs) = if let Some((lhs, rhs)) = raw.split_once("!~") {
        (lhs, "!~", rhs)
    } else if let Some((lhs, rhs)) = raw.split_once('~') {
        (lhs, "~", rhs)
    } else {
        return None;
    };
    if lhs.trim() != "FILTER" {
        return None;
    }
    let rhs = rhs.trim();
    let id = rhs.strip_prefix('"')?.strip_suffix('"')?;
    Some((id.to_string(), op == "!~"))
}

fn normalize_filter_expr(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    let mut chars = raw.chars().peekable();
    let mut in_string = false;
    let mut prev_non_ws = None;
    while let Some(ch) = chars.next() {
        if ch == '"' {
            in_string = !in_string;
            out.push(ch);
            prev_non_ws = Some(ch);
            continue;
        }
        if in_string {
            out.push(ch);
            continue;
        }
        let next = chars.peek().copied();
        match ch {
            '=' if !matches!(next, Some('=') | Some('~'))
                && !matches!(prev_non_ws, Some('!' | '<' | '>' | '=')) =>
            {
                out.push_str("==")
            }
            '&' if next == Some('&') => {
                out.push_str("&&");
                chars.next();
            }
            '&' => {
                out.push_str("&&");
            }
            '|' if next == Some('|') => {
                out.push_str("||");
                chars.next();
            }
            '|' => {
                out.push_str("||");
            }
            '~' if !matches!(prev_non_ws, Some('!' | '=')) => {
                out.push_str("=~");
            }
            _ => out.push(ch),
        }
        if !ch.is_whitespace() {
            prev_non_ws = Some(ch);
        }
    }
    out
}

fn lookup_filter_symbol(src: &str, record: &TextRecord<'_>) -> Option<(Value, usize)> {
    let had_percent = src.starts_with('%');
    let token_src = src.strip_prefix('%').unwrap_or(src);
    let len = token_src
        .char_indices()
        .take_while(|(_, ch)| ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | '/'))
        .map(|(idx, ch)| idx + ch.len_utf8())
        .last()?;
    if token_src[len..].starts_with('(') {
        return None;
    }
    let token = &token_src[..len];
    let value = filter_token_value(token, had_percent, record);
    let consumed = len + usize::from(had_percent);
    Some((value, consumed))
}

fn filter_token_value(
    token: &str,
    explicit_formatter_token: bool,
    record: &TextRecord<'_>,
) -> Value {
    let raw = if token.eq_ignore_ascii_case("type") {
        record.core("TYPE")
    } else if explicit_formatter_token && token == "ILEN" {
        record.computed_ilen().to_string()
    } else {
        render_token(token, record, None)
    };
    filter_value(raw)
}

fn filter_value(raw: String) -> Value {
    if raw == "." {
        return Value::string(raw);
    }
    raw.parse::<f64>()
        .map(Value::number)
        .unwrap_or_else(|_| Value::string(raw))
}

fn query_sample_indices(
    samples: &[String],
    sample_list: Option<&str>,
    sample_list_is_file: bool,
) -> io::Result<Vec<usize>> {
    let flags = crate::smpl_ilist::SMPL_STRICT | crate::smpl_ilist::SMPL_REORDER;
    Ok(crate::smpl_ilist::init(samples, sample_list, sample_list_is_file, flags)?.idx)
}

#[derive(Debug)]
struct TextRecord<'a> {
    fields: Vec<&'a str>,
    samples: &'a [String],
    sample_indices: &'a [usize],
}

impl<'a> TextRecord<'a> {
    fn parse(line: &'a str, samples: &'a [String], sample_indices: &'a [usize]) -> Self {
        Self {
            fields: line.split('\t').collect(),
            samples,
            sample_indices,
        }
    }

    fn core(&self, key: &str) -> String {
        match key {
            "CHROM" => self.fields.first().copied().unwrap_or(".").to_string(),
            "POS" => self.fields.get(1).copied().unwrap_or(".").to_string(),
            "ID" => self.fields.get(2).copied().unwrap_or(".").to_string(),
            "REF" => self.fields.get(3).copied().unwrap_or(".").to_string(),
            "ALT" => self.fields.get(4).copied().unwrap_or(".").to_string(),
            "QUAL" => self.fields.get(5).copied().unwrap_or(".").to_string(),
            "FILTER" => self.fields.get(6).copied().unwrap_or(".").to_string(),
            "INFO" => self.fields.get(7).copied().unwrap_or(".").to_string(),
            "FORMAT" => self.fields.get(8).copied().unwrap_or(".").to_string(),
            "N_ALT" => self.n_alt().to_string(),
            "N_SAMPLES" => self.samples.len().to_string(),
            "TYPE" => self.variant_type_label(),
            "LINE" => {
                let mut line = self.fields.join("\t");
                line.push('\n');
                line
            }
            _ => ".".to_string(),
        }
    }

    fn n_alt(&self) -> usize {
        match self.fields.get(4).copied().unwrap_or(".") {
            "." => 0,
            alt => alt.split(',').filter(|allele| !allele.is_empty()).count(),
        }
    }

    fn format_with_selected_samples(&self) -> String {
        let mut out = self.fields.get(8).copied().unwrap_or(".").to_string();
        for &sample_index in self.sample_indices {
            if let Some(sample) = self.fields.get(9 + sample_index) {
                out.push('\t');
                out.push_str(sample);
            }
        }
        out
    }

    fn variant_type_label(&self) -> String {
        let ref_allele = self.fields.get(3).copied().unwrap_or(".");
        let alt = self.fields.get(4).copied().unwrap_or(".");
        let mut variant_type = VariantType::REF;
        for alt_allele in alt.split(',').filter(|allele| !allele.is_empty()) {
            variant_type |= classify_variant(ref_allele, alt_allele).variant_type;
        }
        variant_type_label(variant_type)
    }

    fn computed_ilen(&self) -> i32 {
        let ref_len = self.fields.get(3).copied().unwrap_or(".").len() as i32;
        self.fields
            .get(4)
            .copied()
            .unwrap_or(".")
            .split(',')
            .filter(|allele| !allele.is_empty() && *allele != ".")
            .map(|allele| (allele.len() as i32 - ref_len).abs())
            .max()
            .unwrap_or(0)
    }

    fn info(&self, key: &str) -> String {
        let Some(info) = self.fields.get(7) else {
            return ".".into();
        };
        for field in info.split(';') {
            let (name, value) = field.split_once('=').unwrap_or((field, "1"));
            if name == key {
                return value.to_string();
            }
        }
        ".".into()
    }

    fn filter_has_id(&self, id: &str) -> bool {
        self.fields
            .get(6)
            .copied()
            .unwrap_or(".")
            .split(';')
            .any(|filter| filter == id)
    }

    fn filter_values(&self, key: &str, vector_any: bool) -> Vec<String> {
        let value = if key.eq_ignore_ascii_case("type") {
            self.core("TYPE")
        } else {
            render_token(key, self, None)
        };
        if vector_any || key == "ALT" {
            value.split(',').map(ToOwned::to_owned).collect()
        } else {
            vec![value]
        }
    }

    fn format_value(&self, sample_index: usize, key: &str) -> String {
        let Some(format) = self.fields.get(8) else {
            return ".".into();
        };
        let Some(&header_sample_index) = self.sample_indices.get(sample_index) else {
            return ".".into();
        };
        let Some(sample) = self.fields.get(9 + header_sample_index) else {
            return ".".into();
        };
        for (idx, name) in format.split(':').enumerate() {
            if name == key {
                return sample.split(':').nth(idx).unwrap_or(".").to_string();
            }
        }
        ".".into()
    }
}

fn render_format(format: &str, record: &TextRecord<'_>) -> String {
    let mut out = String::new();
    let mut rest = format;
    while let Some(start) = rest.find('[') {
        render_segment(&rest[..start], record, None, &mut out);
        let after_start = &rest[start + 1..];
        let Some(end) = after_start.find(']') else {
            render_segment(&rest[start..], record, None, &mut out);
            return out;
        };
        let block = &after_start[..end];
        for sample_index in 0..record.sample_indices.len() {
            render_segment(block, record, Some(sample_index), &mut out);
        }
        rest = &after_start[end + 1..];
    }
    render_segment(rest, record, None, &mut out);
    out
}

fn render_format_header(
    format: &str,
    samples: &[String],
    sample_indices: &[usize],
    header_level: u8,
) -> String {
    let mut out = String::from("#");
    let mut counter = 1usize;
    let mut rest = format;
    let indexed = header_level == 1;

    while let Some(start) = rest.find('[') {
        render_header_segment(&rest[..start], None, indexed, &mut counter, &mut out);
        let after_start = &rest[start + 1..];
        let Some(end) = after_start.find(']') else {
            render_header_segment(&rest[start..], None, indexed, &mut counter, &mut out);
            return finish_header(out);
        };
        let block = &after_start[..end];
        if segment_has_newline(block) {
            render_header_segment(block, None, indexed, &mut counter, &mut out);
        } else {
            for &sample_index in sample_indices {
                let sample = samples.get(sample_index).map(String::as_str);
                render_header_segment(block, sample, indexed, &mut counter, &mut out);
            }
        }
        rest = &after_start[end + 1..];
    }

    render_header_segment(rest, None, indexed, &mut counter, &mut out);
    finish_header(out)
}

fn finish_header(mut out: String) -> String {
    if !out.ends_with('\n') {
        out.push('\n');
    }
    out
}

fn segment_has_newline(segment: &str) -> bool {
    segment.contains('\n') || segment.contains("\\n")
}

fn render_header_segment(
    segment: &str,
    sample_prefix: Option<&str>,
    indexed: bool,
    counter: &mut usize,
    out: &mut String,
) {
    let mut chars = segment.char_indices().peekable();
    while let Some((idx, ch)) = chars.next() {
        match ch {
            '\\' => {
                if let Some((_, next)) = chars.next() {
                    match next {
                        'n' => out.push('\n'),
                        't' => out.push('\t'),
                        _ => out.push(next),
                    }
                } else {
                    out.push(ch);
                }
            }
            '%' => {
                let token_start = idx + ch.len_utf8();
                let mut token_end = token_start;
                while let Some(&(next_idx, next)) = chars.peek() {
                    if next.is_ascii_alphanumeric() || matches!(next, '_' | '.' | '/') {
                        token_end = next_idx + next.len_utf8();
                        chars.next();
                    } else {
                        break;
                    }
                }
                let label = header_token_label(&segment[token_start..token_end]);
                if indexed {
                    out.push_str(&format!("[{}]", *counter));
                }
                if let Some(sample) = sample_prefix {
                    out.push_str(sample);
                    out.push(':');
                }
                out.push_str(&label);
                *counter += 1;
            }
            _ => out.push(ch),
        }
    }
}

fn header_token_label(token: &str) -> String {
    let token = token.strip_prefix('/').unwrap_or(token);
    token
        .strip_prefix("INFO/")
        .or_else(|| token.strip_prefix("FMT/"))
        .or_else(|| token.strip_prefix("FORMAT/"))
        .unwrap_or(token)
        .to_string()
}

fn render_segment(
    segment: &str,
    record: &TextRecord<'_>,
    sample_index: Option<usize>,
    out: &mut String,
) {
    let mut chars = segment.char_indices().peekable();
    while let Some((idx, ch)) = chars.next() {
        match ch {
            '\\' => {
                if let Some((_, next)) = chars.next() {
                    match next {
                        'n' => out.push('\n'),
                        't' => out.push('\t'),
                        _ => out.push(next),
                    }
                } else {
                    out.push(ch);
                }
            }
            '%' => {
                let token_start = idx + ch.len_utf8();
                let mut token_end = token_start;
                while let Some(&(next_idx, next)) = chars.peek() {
                    if next.is_ascii_alphanumeric() || matches!(next, '_' | '.' | '/') {
                        token_end = next_idx + next.len_utf8();
                        chars.next();
                    } else {
                        break;
                    }
                }
                let token = &segment[token_start..token_end];
                out.push_str(&render_token(token, record, sample_index));
            }
            _ => out.push(ch),
        }
    }
}

fn render_token(token: &str, record: &TextRecord<'_>, sample_index: Option<usize>) -> String {
    let force_record_namespace = token.starts_with('/');
    let token = token.strip_prefix('/').unwrap_or(token);
    if token == "SAMPLE" {
        return sample_index
            .and_then(|i| record.sample_indices.get(i))
            .and_then(|&i| record.samples.get(i))
            .cloned()
            .unwrap_or_else(|| ".".into());
    }
    if let Some(key) = token.strip_prefix("INFO/") {
        return record.info(key);
    }
    if token == "FORMAT" && sample_index.is_none() {
        return record.format_with_selected_samples();
    }
    if let Some(key) = token
        .strip_prefix("FMT/")
        .or_else(|| token.strip_prefix("FORMAT/"))
    {
        return sample_index
            .map(|i| record.format_value(i, key))
            .unwrap_or_else(|| ".".into());
    }
    if !force_record_namespace && let Some(i) = sample_index {
        let value = record.format_value(i, token);
        if value != "." {
            return value;
        }
    }
    match token {
        "CHROM" | "POS" | "ID" | "REF" | "ALT" | "QUAL" | "FILTER" | "INFO" | "FORMAT"
        | "N_ALT" | "N_SAMPLES" | "TYPE" | "LINE" => record.core(token),
        _ => record.info(token),
    }
}

fn variant_type_label(variant_type: VariantType) -> String {
    if variant_type == VariantType::REF {
        return "ref".into();
    }
    let mut labels = Vec::new();
    if variant_type.contains(VariantType::SNP) {
        labels.push("snp");
    }
    if variant_type.contains(VariantType::MNP) {
        labels.push("mnp");
    }
    if variant_type.contains(VariantType::INDEL) {
        labels.push("indel");
    }
    if variant_type.contains(VariantType::BND) {
        labels.push("bnd");
    }
    if variant_type.contains(VariantType::OVERLAP) {
        labels.push("overlap");
    }
    if variant_type.contains(VariantType::OTHER) {
        labels.push("other");
    }
    labels.join(",")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_query_options<'a>() -> QueryFormatOptions<'a> {
        QueryFormatOptions {
            sample_list: None,
            sample_list_is_file: false,
            header_level: 0,
            region_spec: None,
            filter_spec: None,
        }
    }

    #[test]
    fn parses_list_samples_mode() {
        let argv = [
            OsString::from("query"),
            OsString::from("-l"),
            OsString::from("in.vcf"),
        ];
        assert_eq!(
            parse_args(&argv).unwrap(),
            Args {
                list_samples: true,
                format: None,
                header_level: 0,
                samples: None,
                samples_is_file: false,
                regions: None,
                filter: None,
                input: "in.vcf".into()
            }
        );
    }

    #[test]
    fn parses_format_mode() {
        let argv = [
            OsString::from("query"),
            OsString::from("-f"),
            OsString::from("%CHROM\n"),
            OsString::from("in.vcf"),
        ];
        assert_eq!(
            parse_args(&argv).unwrap(),
            Args {
                list_samples: false,
                format: Some("%CHROM\n".into()),
                header_level: 0,
                samples: None,
                samples_is_file: false,
                regions: None,
                filter: None,
                input: "in.vcf".into()
            }
        );
    }

    #[test]
    fn renders_core_and_sample_fields() {
        let text = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t2\trs1\tA\tC\t.\tPASS\tDP=7\tGT:DP\t0/1:3\t0/0:4\n";
        let mut out = Vec::new();
        query_format_text(
            text,
            "%CHROM\\t%POS\\t%DP[\\t%SAMPLE=%GT:%DP]\\n",
            &default_query_options(),
            &mut out,
        )
        .unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "1\t2\t7\tA=0/1:3\tB=0/0:4\n"
        );
    }

    #[test]
    fn sample_selection_reorders_format_loops() {
        let text = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t00\t11\n\
chr1\t10000\t.\tA\tC\t.\t.\t.\tGT\t0/0\t1/1\n";
        let mut out = Vec::new();
        let options = QueryFormatOptions {
            sample_list: Some("11,00"),
            ..default_query_options()
        };
        query_format_text(text, "[%SAMPLE %GT\\n]", &options, &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "11 1/1\n00 0/0\n");
    }

    #[test]
    fn renders_indexed_headers_for_sample_blocks() {
        let samples = vec!["C".to_string(), "D".to_string()];
        let indices = vec![0, 1];
        assert_eq!(
            render_format_header("%CHROM %POS[ %SAMPLE %DP %GT]\\n", &samples, &indices, 1),
            "#[1]CHROM [2]POS [3]C:SAMPLE [4]C:DP [5]C:GT [6]D:SAMPLE [7]D:DP [8]D:GT\n"
        );
        assert_eq!(
            render_format_header("%CHROM %POS[ %SAMPLE][ %DP][ %GT]", &samples, &indices, 2),
            "#CHROM POS C:SAMPLE D:SAMPLE C:DP D:DP C:GT D:GT\n"
        );
    }

    #[test]
    fn region_filter_supports_inline_and_exclusion() {
        let filter = RegionFilter::from_spec(&RegionFilterSpec {
            raw: "1:10-20".into(),
            is_file: false,
            exclude: false,
        })
        .unwrap();
        assert!(filter.matches("1\t10\t.\tA\tC\t.\t.\t.").unwrap());
        assert!(!filter.matches("1\t21\t.\tA\tC\t.\t.\t.").unwrap());

        let filter = RegionFilter::from_spec(&RegionFilterSpec {
            raw: "1:10-20".into(),
            is_file: false,
            exclude: true,
        })
        .unwrap();
        assert!(!filter.matches("1\t10\t.\tA\tC\t.\t.\t.").unwrap());
        assert!(filter.matches("1\t21\t.\tA\tC\t.\t.\t.").unwrap());
    }

    #[test]
    fn query_filter_supports_core_info_and_exclusion() {
        let text = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\trs1\tA\tC\t.\tPASS\tDP=7;STR=abc\n\
1\t20\trs2\tG\tT\t.\tq10\tDP=2;STR=xyz\n";
        let mut out = Vec::new();
        let include_spec = FilterSpec {
            raw: r#"DP>=7 && STR="abc""#.into(),
            exclude: false,
        };
        let options = QueryFormatOptions {
            filter_spec: Some(&include_spec),
            ..default_query_options()
        };
        query_format_text(text, "%CHROM:%POS:%STR\\n", &options, &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "1:10:abc\n");

        let mut out = Vec::new();
        let exclude_spec = FilterSpec {
            raw: r#"FILTER="PASS""#.into(),
            exclude: true,
        };
        let options = QueryFormatOptions {
            filter_spec: Some(&exclude_spec),
            ..default_query_options()
        };
        query_format_text(text, "%ID\\n", &options, &mut out).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "rs2\n");
    }
}
