//! Port of `bcftools sort` (upstream `vcfsort.c`).
//!
//! This implements the command shape VNtyper uses after Kestrel:
//!
//! ```text
//! bcftools sort input.vcf -o output.vcf.gz -W -O z
//! ```
//!
//! The Rust port preserves the same coordinate/ref/alt ordering, writes
//! VCF/VCF.gz, supports automatic CSI/TBI indexing with `-W`, and spills sorted
//! runs to temporary files when `--max-mem` is exceeded.

use std::cmp::Ordering;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    build_vcf_csi_from_path_with_min_shift, build_vcf_tbi_from_path, write_csi, write_tbi,
};
use htslib_rs::vcf;
use htslib_rs::vcf::variant::io::Write as _;

use crate::diagnostics::fmt_etag;
use crate::io::{VariantOutputFormat, apply_verbosity, init_index2, write_index_parse};
use crate::vcf_compat::NormalizeFileformat;

const USAGE: &str = "\n\
About:   Sort VCF/BCF file.\n\
Usage:   bcftools sort [OPTIONS] <FILE.vcf>\n\
\n\
Options:\n\
    -m, --max-mem FLOAT[kMG]       Maximum memory to use [768M]\n\
    -o, --output FILE              Output file name [stdout]\n\
    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n\
    -T, --temp-dir DIR             Temporary files [/tmp/bcftools.XXXXXX]\n\
        --threads INT              Use multithreaded BGZF compression for compressed output\n\
    -v, --verbosity INT            Verbosity level\n\
    -W, --write-index[=FMT]        Automatically index the output files [off]\n\
\n";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
}

impl OutputKind {
    fn parse(raw: &str) -> Option<Self> {
        let ty = raw.chars().next()?;
        match ty {
            'v' => Some(Self::VcfText),
            'z' => Some(Self::VcfGz),
            '0'..='9' => Some(Self::VcfGz),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    write_index: Option<i32>,
    max_mem: usize,
    temp_dir: PathBuf,
    thread_count: Option<NonZero<usize>>,
}

/// Subcommand entry point. `argv[0]` is `"sort"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfsort", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfsort", &message));
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
    let mut write_index = None;
    let mut max_mem = 768_000_000usize;
    let mut temp_dir = std::env::temp_dir();
    let mut thread_count = None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-o" | "--output" | "--output-file" => {
                output = Some(next_path(&mut iter, "--output")?);
            }
            "-O" | "--output-type" => {
                let value = next_string(&mut iter, "--output-type")?;
                output_kind = parse_output_kind(&value)?;
            }
            "-W" | "--write-index" => {
                write_index = parse_write_index(None)?;
            }
            "-m" | "--max-mem" => {
                max_mem = parse_max_mem(&next_string(&mut iter, raw.as_ref())?)?;
            }
            "-T" | "--temp-dir" => {
                temp_dir = PathBuf::from(next_string(&mut iter, raw.as_ref())?);
            }
            "--threads" => {
                thread_count = parse_threads(&next_string(&mut iter, "--threads")?)?;
            }
            "-v" | "--verbosity" => {
                let value = next_string(&mut iter, "--verbosity")?;
                if apply_verbosity(&value).is_err() {
                    return Err(ParseOutcome::Error(format!(
                        "Could not parse argument: --verbosity {value}"
                    )));
                }
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--output-file=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--write-index=") => {
                write_index = parse_write_index(Some(value_after_equals(&raw)))?;
            }
            _ if raw.starts_with("--max-mem=") => {
                max_mem = parse_max_mem(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--temp-dir=") => {
                temp_dir = PathBuf::from(value_after_equals(&raw));
            }
            _ if raw.starts_with("--threads=") => {
                thread_count = parse_threads(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("-m") && raw.len() > 2 => {
                max_mem = parse_max_mem(&raw[2..])?;
            }
            _ if raw.starts_with("-T") && raw.len() > 2 => {
                temp_dir = PathBuf::from(&raw[2..]);
            }
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?;
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => {
                output = Some(PathBuf::from(&raw[2..]));
            }
            _ if raw.starts_with("-W=") => {
                write_index = parse_write_index(Some(&raw[3..]))?;
            }
            _ if raw.starts_with("-W") && raw.len() > 2 => {
                write_index = parse_write_index(Some(&raw[2..]))?;
            }
            _ if raw.starts_with('-') => return Err(ParseOutcome::Usage),
            _ => {
                if input.is_some() {
                    return Err(ParseOutcome::Error(format!(
                        "multiple input files are not yet supported: {raw}"
                    )));
                }
                input = Some(PathBuf::from(arg));
            }
        }
    }

    let input = input.ok_or(ParseOutcome::Usage)?;
    Ok(Args {
        input,
        output,
        output_kind,
        write_index,
        max_mem,
        temp_dir,
        thread_count,
    })
}

fn run(args: &Args) -> io::Result<()> {
    if args.write_index.is_some() && args.output.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W requires an output file",
        ));
    }

    let fmt = format::detect_path(&args.input).map_err(|e| io::Error::other(e.to_string()))?;
    let (header, mut records) = read_records(&args.input, fmt)?;
    records = if estimated_records_bytes(&records) > args.max_mem {
        external_sort_records(&header, records, &args.temp_dir, args.max_mem)?
    } else {
        records.sort_by(|a, b| compare_records(&header, a, b));
        records
    };

    match args.output.as_deref() {
        Some(path) => write_to_path(path, args.output_kind, &header, &records, args.thread_count)?,
        None => write_vcf(io::stdout().lock(), &header, &records)?,
    }

    if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
        write_index(path, index_format)?;
    }

    Ok(())
}

fn estimated_records_bytes(records: &[vcf::variant::RecordBuf]) -> usize {
    records
        .iter()
        .map(|record| {
            128 + record.reference_sequence_name().len()
                + record.reference_bases().len()
                + alternate_bases_key(record).len()
        })
        .sum()
}

fn external_sort_records(
    header: &vcf::Header,
    records: Vec<vcf::variant::RecordBuf>,
    temp_dir: &Path,
    max_mem: usize,
) -> io::Result<Vec<vcf::variant::RecordBuf>> {
    fs::create_dir_all(temp_dir)?;
    let estimated_bytes = estimated_records_bytes(&records);
    let average_record_bytes = 1.max(estimated_bytes / records.len().max(1));
    let chunk_size = 1.max(max_mem / average_record_bytes);
    let mut chunk_paths = Vec::new();
    for (i, chunk) in records.chunks(chunk_size).enumerate() {
        let mut chunk = chunk.to_vec();
        chunk.sort_by(|a, b| compare_records(header, a, b));
        let path = temp_dir.join(format!("bcftools-rs-sort.{}.{}.vcf", std::process::id(), i));
        write_to_path(&path, OutputKind::VcfText, header, &chunk, None)?;
        chunk_paths.push(path);
    }

    let mut runs = Vec::new();
    for path in &chunk_paths {
        let (_, records) = read_records(
            path,
            format::Format {
                compression: Compression::None,
                exact: Exact::Vcf,
                category: format::Category::VariantData,
            },
        )?;
        runs.push(records);
    }

    let sorted = merge_sorted_runs(header, runs);
    for path in chunk_paths {
        let _ = fs::remove_file(path);
    }
    Ok(sorted)
}

fn merge_sorted_runs(
    header: &vcf::Header,
    runs: Vec<Vec<vcf::variant::RecordBuf>>,
) -> Vec<vcf::variant::RecordBuf> {
    let mut positions = vec![0usize; runs.len()];
    let mut out = Vec::new();

    loop {
        let mut best: Option<usize> = None;
        for (i, run) in runs.iter().enumerate() {
            if positions[i] >= run.len() {
                continue;
            }
            if best.is_none_or(|j| {
                compare_records(header, &run[positions[i]], &runs[j][positions[j]])
                    == Ordering::Less
            }) {
                best = Some(i);
            }
        }
        let Some(i) = best else {
            break;
        };
        out.push(runs[i][positions[i]].clone());
        positions[i] += 1;
    }

    out
}

fn read_records(
    path: &Path,
    fmt: format::Format,
) -> io::Result<(vcf::Header, Vec<vcf::variant::RecordBuf>)> {
    use htslib_rs::bcf;

    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let header = reader.read_header()?;
        let records = reader
            .record_bufs(&header)
            .collect::<io::Result<Vec<_>>>()?;
        return Ok((header, records));
    }

    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let normalized = NormalizeFileformat::new(BufReader::new(dec))?;
        let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
        let header = reader.read_header()?;
        let records = reader
            .records()
            .map(|result| {
                let record = result?;
                vcf::variant::RecordBuf::try_from_variant_record(&header, &record)
            })
            .collect::<io::Result<Vec<_>>>()?;
        return Ok((header, records));
    }

    let file = File::open(path)?;
    let normalized = NormalizeFileformat::new(BufReader::new(file))?;
    let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
    let header = reader.read_header()?;
    let records = reader
        .records()
        .map(|result| {
            let record = result?;
            vcf::variant::RecordBuf::try_from_variant_record(&header, &record)
        })
        .collect::<io::Result<Vec<_>>>()?;
    Ok((header, records))
}

fn compare_records(
    header: &vcf::Header,
    a: &vcf::variant::RecordBuf,
    b: &vcf::variant::RecordBuf,
) -> Ordering {
    contig_order(header, a.reference_sequence_name())
        .cmp(&contig_order(header, b.reference_sequence_name()))
        .then_with(|| a.variant_start().cmp(&b.variant_start()))
        .then_with(|| {
            a.reference_bases()
                .to_ascii_lowercase()
                .cmp(&b.reference_bases().to_ascii_lowercase())
        })
        .then_with(|| alternate_bases_key(a).cmp(&alternate_bases_key(b)))
}

fn contig_order(header: &vcf::Header, name: &str) -> usize {
    header.contigs().get_index_of(name).unwrap_or(usize::MAX)
}

fn alternate_bases_key(record: &vcf::variant::RecordBuf) -> String {
    record
        .alternate_bases()
        .as_ref()
        .join(",")
        .to_ascii_lowercase()
}

fn write_to_path(
    path: &Path,
    output_kind: OutputKind,
    header: &vcf::Header,
    records: &[vcf::variant::RecordBuf],
    thread_count: Option<NonZero<usize>>,
) -> io::Result<()> {
    let file = File::create(path)?;
    match output_kind {
        OutputKind::VcfText => write_vcf(file, header, records),
        OutputKind::VcfGz => {
            if let Some(thread_count) = thread_count {
                let bgzf =
                    htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
                let mut writer = vcf::io::Writer::new(bgzf);
                write_records(&mut writer, header, records)?;
                let mut bgzf = writer.into_inner();
                let _file = bgzf.finish()?;
            } else {
                let bgzf = htslib_rs::bgzf::io::Writer::new(file);
                let mut writer = vcf::io::Writer::new(bgzf);
                write_records(&mut writer, header, records)?;
                let bgzf = writer.into_inner();
                let _file = bgzf.finish()?;
            }
            Ok(())
        }
    }
}

fn write_vcf<W: Write>(
    out: W,
    header: &vcf::Header,
    records: &[vcf::variant::RecordBuf],
) -> io::Result<()> {
    let mut writer = vcf::io::Writer::new(out);
    write_records(&mut writer, header, records)
}

fn write_records<W: Write>(
    writer: &mut vcf::io::Writer<W>,
    header: &vcf::Header,
    records: &[vcf::variant::RecordBuf],
) -> io::Result<()> {
    writer.write_header(header)?;
    for record in records {
        writer.write_variant_record(header, record)?;
    }
    Ok(())
}

fn write_index(path: &Path, index_format: i32) -> io::Result<()> {
    let path = path
        .to_str()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "non-UTF-8 output path"))?;
    let Some(plan) = init_index2(Some(path), index_format, VariantOutputFormat::Vcf)
        .map_err(io::Error::other)?
    else {
        return Ok(());
    };

    if plan.min_shift == 0 {
        let index = build_vcf_tbi_from_path(path)?;
        write_tbi(plan.index_path, &index)
    } else {
        let index = build_vcf_csi_from_path_with_min_shift(path, plan.min_shift)?;
        write_csi(plan.index_path, &index)
    }
}

fn parse_output_kind(raw: &str) -> Result<OutputKind, ParseOutcome> {
    OutputKind::parse(raw)
        .ok_or_else(|| ParseOutcome::Error(format!("The output type \"{raw}\" not recognised")))
}

fn parse_max_mem(raw: &str) -> Result<usize, ParseOutcome> {
    let (number, multiplier) = match raw.as_bytes().last().copied() {
        Some(b'k' | b'K') => (&raw[..raw.len() - 1], 1_000f64),
        Some(b'm' | b'M') => (&raw[..raw.len() - 1], 1_000_000f64),
        Some(b'g' | b'G') => (&raw[..raw.len() - 1], 1_000_000_000f64),
        _ => (raw, 1f64),
    };
    let value = number
        .parse::<f64>()
        .map_err(|_| ParseOutcome::Error(format!("Could not parse --max-mem {raw}")))?;
    if value <= 0.0 {
        return Err(ParseOutcome::Error(format!(
            "--max-mem must be positive: {raw}"
        )));
    }
    Ok((value * multiplier) as usize)
}

fn parse_threads(raw: &str) -> Result<Option<NonZero<usize>>, ParseOutcome> {
    raw.parse::<usize>()
        .map(NonZero::new)
        .map_err(|_| ParseOutcome::Error(format!("Could not parse argument: --threads {raw}")))
}

fn parse_write_index(raw: Option<&str>) -> Result<Option<i32>, ParseOutcome> {
    write_index_parse(raw).map(Some).ok_or_else(|| {
        ParseOutcome::Error(format!("Unsupported index format '{}'", raw.unwrap_or("")))
    })
}

fn next_path<'a, I>(iter: &mut I, name: &str) -> Result<PathBuf, ParseOutcome>
where
    I: Iterator<Item = &'a OsString>,
{
    next_string(iter, name).map(PathBuf::from)
}

fn next_string<'a, I>(iter: &mut I, name: &str) -> Result<String, ParseOutcome>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .map(|value| value.to_string_lossy().into_owned())
        .ok_or_else(|| ParseOutcome::Error(format!("missing argument for {name}")))
}

fn value_after_equals(raw: &str) -> &str {
    raw.split_once('=').map(|(_, value)| value).unwrap_or("")
}
