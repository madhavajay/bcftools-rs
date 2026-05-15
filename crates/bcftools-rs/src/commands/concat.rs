//! Port of `bcftools concat` (upstream `vcfconcat.c`).
//!
//! This MVP supports vertical concatenation of VCF/BCF files that share the
//! same sample columns in the same order. Records from each input are emitted
//! in the order they appear, preserving the first input's header. Output type
//! follows `-O v|z|b|u` and the file extension when `-o` is given.
//!
//! Implemented options
//!
//! - `-o`/`--output FILE`
//! - `-O`/`--output-type u|b|v|z[0-9]`
//! - `-f`/`--file-list FILE` — read input list from a file (one per line).
//! - `-G`/`--drop-genotypes` — strip FORMAT and sample columns.
//! - `-D`/`--remove-duplicates` — alias for `-d exact`.
//! - `-d`/`--rm-dups STRING` — drop duplicate records: `snps|indels|both|all|exact`.
//! - `-a`/`--allow-overlaps` — allow adjacent input files to overlap.
//! - `-n`/`--naive` — fast text VCF concatenation preserving the first header.
//! - `--naive-force` — skip header equality checks in naive mode.
//! - `--no-version` — suppress the per-command header lines.
//!
//! Deferred (intentional gaps tracked in TODO.md): `-l/--ligate`,
//! `--ligate-force`, `--ligate-warn`, `-c/--compact-PS`, `-q/--min-PQ`.
//! Some of these depend on synced reader parity in `htslib-rs` and are tracked
//! at the bottom of TODO.md.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    build_bcf_csi_with_min_shift, build_vcf_csi_from_path_with_min_shift, build_vcf_tbi_from_path,
    write_csi, write_tbi,
};
use htslib_rs::vcf;
use htslib_rs::vcf::variant::io::Write as _;

use crate::diagnostics::fmt_etag;
use crate::header_version::{build_lines, command_time};
use crate::io::{
    VariantOutputFormat, apply_verbosity, init_index2, parse_overlap_option, write_index_parse,
};
use crate::vcf_compat::NormalizeFileformat;

const USAGE: &str = "\n\
About: Concatenate or combine VCF/BCF files.\n\
Usage: bcftools concat [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n\
\n\
Options:\n\
   -a, --allow-overlaps           Allow records from adjacent input files to overlap\n\
   -d, --rm-dups STRING           Output duplicate records present in multiple files only once: <snps|indels|both|all|exact>\n\
   -D, --remove-duplicates        Alias for -d exact\n\
   -f, --file-list FILE           Read the list of files from a file.\n\
   -G, --drop-genotypes           Drop individual genotype information.\n\
   -n, --naive                    Concatenate VCF bodies without re-encoding\n\
       --naive-force              Skip header equality checks in naive mode\n\
       --no-version               Do not append version and command line to the header\n\
   -o, --output FILE              Write output to a file [standard output]\n\
   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF [v]\n\
   -r, --regions REGION           Restrict to comma-separated list of regions\n\
   -R, --regions-file FILE        Restrict to regions listed in a file\n\
       --regions-overlap 0|1|2    Include if POS in region (0) or record overlaps (1/2) [0]\n\
       --threads INT              Use multithreaded BGZF compression for compressed output\n\
   -v, --verbosity INT            Set verbosity level\n\
   -W, --write-index[=FMT]        Automatically index the output files [off]\n\
\n";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
    BcfUncompressed,
}

impl OutputKind {
    fn parse(raw: &str) -> Option<Self> {
        match raw.chars().next()? {
            'v' => Some(Self::VcfText),
            'z' => Some(Self::VcfGz),
            'b' => Some(Self::Bcf),
            'u' => Some(Self::BcfUncompressed),
            '0'..='9' => Some(Self::VcfGz),
            _ => None,
        }
    }

    fn for_path(path: &Path) -> Self {
        let lower = path.as_os_str().to_string_lossy().to_ascii_lowercase();
        if lower.ends_with(".bcf") {
            Self::Bcf
        } else if lower.ends_with(".vcf.gz") || lower.ends_with(".vcf.bgz") {
            Self::VcfGz
        } else {
            Self::VcfText
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DupMode {
    Exact,
    Snps,
    Indels,
    Both,
    All,
}

impl DupMode {
    fn parse(raw: &str) -> Option<Self> {
        match raw {
            "exact" => Some(Self::Exact),
            "snps" => Some(Self::Snps),
            "indels" => Some(Self::Indels),
            "both" => Some(Self::Both),
            "all" => Some(Self::All),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct Args {
    inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    drop_genotypes: bool,
    rm_dups: Option<DupMode>,
    allow_overlaps: bool,
    naive: bool,
    naive_force: bool,
    no_version: bool,
    write_index: Option<i32>,
    regions: Vec<RegionSpec>,
    regions_overlap: u8,
    thread_count: Option<NonZero<usize>>,
}

#[derive(Debug, Clone)]
struct RegionSpec {
    contig: String,
    start: Option<i64>,
    end: Option<i64>,
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args, argv) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfconcat", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfconcat", &message));
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
    let mut file_list: Option<PathBuf> = None;
    let mut output: Option<PathBuf> = None;
    let mut explicit_kind: Option<OutputKind> = None;
    let mut drop_genotypes = false;
    let mut rm_dups: Option<DupMode> = None;
    let mut allow_overlaps = false;
    let mut naive = false;
    let mut naive_force = false;
    let mut no_version = false;
    let mut write_index = None;
    let mut regions = Vec::new();
    let mut regions_overlap = 0;
    let mut thread_count = None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "--no-version" => no_version = true,
            "-G" | "--drop-genotypes" => drop_genotypes = true,
            "-a" | "--allow-overlaps" => allow_overlaps = true,
            "-n" | "--naive" => naive = true,
            "--naive-force" => {
                naive = true;
                naive_force = true;
            }
            "-D" | "--remove-duplicates" => rm_dups = Some(DupMode::Exact),
            "-W" | "--write-index" => {
                write_index = parse_write_index(None)?;
            }
            "-d" | "--rm-dups" => {
                let value = next_string(&mut iter, raw.as_ref())?;
                rm_dups = Some(parse_dup_mode(&value)?);
            }
            "-f" | "--file-list" => {
                file_list = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?));
            }
            "-r" | "--regions" => {
                parse_region_list(&mut regions, &next_string(&mut iter, "--regions")?)?;
            }
            "-R" | "--regions-file" => {
                load_region_file(&mut regions, &next_string(&mut iter, "--regions-file")?)?;
            }
            "--regions-overlap" => {
                regions_overlap =
                    parse_regions_overlap(&next_string(&mut iter, "--regions-overlap")?)?;
            }
            "--threads" => {
                thread_count = parse_threads(&next_string(&mut iter, "--threads")?)?;
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?));
            }
            "-O" | "--output-type" => {
                let value = next_string(&mut iter, "--output-type")?;
                explicit_kind = Some(parse_output_kind(&value)?);
            }
            "-v" | "--verbosity" => {
                let value = next_string(&mut iter, "--verbosity")?;
                if apply_verbosity(&value).is_err() {
                    return Err(ParseOutcome::Error(format!(
                        "Could not parse argument: --verbosity {value}"
                    )));
                }
            }
            _ if raw.starts_with("--rm-dups=") => {
                rm_dups = Some(parse_dup_mode(value_after_equals(&raw))?);
            }
            _ if raw.starts_with("--file-list=") => {
                file_list = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--output-type=") => {
                explicit_kind = Some(parse_output_kind(value_after_equals(&raw))?);
            }
            _ if raw.starts_with("--write-index=") => {
                write_index = parse_write_index(Some(value_after_equals(&raw)))?;
            }
            _ if raw.starts_with("--regions=") => {
                parse_region_list(&mut regions, value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--regions-file=") => {
                load_region_file(&mut regions, value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--regions-overlap=") => {
                regions_overlap = parse_regions_overlap(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--threads=") => {
                thread_count = parse_threads(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                explicit_kind = Some(parse_output_kind(&raw[2..])?);
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => {
                output = Some(PathBuf::from(&raw[2..]));
            }
            _ if raw.starts_with("-W=") => {
                write_index = parse_write_index(Some(&raw[3..]))?;
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!("Unrecognized option: {raw}")));
            }
            _ => inputs.push(PathBuf::from(arg)),
        }
    }

    if let Some(list) = file_list {
        let f = File::open(&list).map_err(|e| {
            ParseOutcome::Error(format!("Could not read file list {}: {e}", list.display()))
        })?;
        for line in BufReader::new(f).lines() {
            let line = line.map_err(|e| ParseOutcome::Error(e.to_string()))?;
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            inputs.push(PathBuf::from(trimmed));
        }
    }

    if inputs.is_empty() {
        return Err(ParseOutcome::Usage);
    }

    let output_kind = explicit_kind.unwrap_or_else(|| {
        output
            .as_deref()
            .map(OutputKind::for_path)
            .unwrap_or(OutputKind::VcfText)
    });

    Ok(Args {
        inputs,
        output,
        output_kind,
        drop_genotypes,
        rm_dups,
        allow_overlaps,
        naive,
        naive_force,
        no_version,
        write_index,
        regions,
        regions_overlap,
        thread_count,
    })
}

fn run(args: &Args, argv: &[OsString]) -> io::Result<()> {
    if args.write_index.is_some() && args.output.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W requires an output file",
        ));
    }
    if args.naive {
        return run_naive(args);
    }

    let mut header = read_header(&args.inputs[0])?;
    if !args.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_concat", &prog_argv, command_time());
        for line in [&lines.version_line, &lines.command_line] {
            htslib_rs::header_compat::append_line(&mut header, line)?;
        }
    }
    let header = if args.drop_genotypes {
        drop_genotype_columns(&header)
    } else {
        header
    };

    let writer: Box<dyn ConcatWriter> = match args.output_kind {
        OutputKind::VcfText => match &args.output {
            Some(path) => Box::new(VcfTextWriter::new_file(path, &header, args.no_version)?),
            None => Box::new(VcfTextWriter::new_stdout(&header, args.no_version)?),
        },
        OutputKind::VcfGz => Box::new(VcfGzWriter::new(
            args.output.as_deref(),
            &header,
            args.no_version,
            args.thread_count,
        )?),
        OutputKind::Bcf | OutputKind::BcfUncompressed => Box::new(BcfWriter::new(
            args.output.as_deref(),
            &header,
            args.thread_count,
        )?),
    };

    let mut writer = writer;
    let mut seen: Option<DupSet> = args.rm_dups.map(|_| DupSet::default());
    let mut previous_file_last: Option<RecordSpan> = None;

    for (i, input) in args.inputs.iter().enumerate() {
        let input_header = read_header(input)?;
        if i > 0 {
            check_sample_columns(&header, &input_header, input)?;
        }
        let mut current_file_last = None;
        let mut checked_first_emitted = false;
        for_each_record(input, &input_header, |rec| {
            if !record_in_regions(rec, &args.regions, args.regions_overlap) {
                return Ok(());
            }
            let span = RecordSpan::from(rec);
            if args.drop_genotypes {
                let projected = strip_format_and_samples(rec);
                if let Some(seen) = seen.as_mut()
                    && let Some(mode) = args.rm_dups
                    && seen.is_dup(&projected, mode)
                {
                    return Ok(());
                }
                check_file_overlap(
                    args.allow_overlaps,
                    &mut checked_first_emitted,
                    previous_file_last.as_ref(),
                    &span,
                )?;
                current_file_last = Some(span);
                writer.write_record(&projected)
            } else {
                if let Some(seen) = seen.as_mut()
                    && let Some(mode) = args.rm_dups
                    && seen.is_dup(rec, mode)
                {
                    return Ok(());
                }
                check_file_overlap(
                    args.allow_overlaps,
                    &mut checked_first_emitted,
                    previous_file_last.as_ref(),
                    &span,
                )?;
                current_file_last = Some(span);
                writer.write_record(rec)
            }
        })?;
        if current_file_last.is_some() {
            previous_file_last = current_file_last;
        }
    }

    writer.finish()?;

    if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
        write_index(path, args.output_kind, index_format)?;
    }

    Ok(())
}

fn check_file_overlap(
    allow_overlaps: bool,
    checked_first_emitted: &mut bool,
    previous_file_last: Option<&RecordSpan>,
    span: &RecordSpan,
) -> io::Result<()> {
    if allow_overlaps || *checked_first_emitted {
        return Ok(());
    }
    *checked_first_emitted = true;
    if let Some(previous) = previous_file_last
        && previous.overlaps(span)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Input files overlap at {}:{}; use -a/--allow-overlaps to concatenate overlapping files",
                span.chrom, span.start
            ),
        ));
    }
    Ok(())
}

fn run_naive(args: &Args) -> io::Result<()> {
    if args.output_kind == OutputKind::Bcf || args.output_kind == OutputKind::BcfUncompressed {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "naive concat currently supports VCF/VCF.gz output only",
        ));
    }
    if args.drop_genotypes || args.rm_dups.is_some() || !args.regions.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "naive concat cannot be combined with record-transforming options",
        ));
    }

    let mut buffer = Vec::new();
    let mut expected_header: Option<String> = None;
    for (i, input) in args.inputs.iter().enumerate() {
        let text = read_vcf_text(input)?;
        let header_end = header_text_len(&text);
        let header = &text[..header_end];
        let comparable_header = comparable_naive_header(header);
        if i == 0 {
            expected_header = Some(comparable_header);
            buffer.extend_from_slice(header.as_bytes());
        } else if !args.naive_force && expected_header.as_deref() != Some(&comparable_header) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Different headers in {} (use --naive-force to override)",
                    input.display()
                ),
            ));
        }
        buffer.extend_from_slice(&text.as_bytes()[header_end..]);
    }

    write_naive_output(args, &buffer)?;
    if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
        write_index(path, args.output_kind, index_format)?;
    }
    Ok(())
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return read_bcf_as_vcf_text(path);
    }
    let mut text = String::new();
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let mut dec = flate2::read::MultiGzDecoder::new(f);
        dec.read_to_string(&mut text)?;
    } else {
        File::open(path)?.read_to_string(&mut text)?;
    }
    crate::vcf_compat::normalize_vcf_text(&mut text);
    Ok(text)
}

fn read_bcf_as_vcf_text(path: &Path) -> io::Result<String> {
    let mut reader = File::open(path).map(htslib_rs::bcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let mut out = Vec::new();
    {
        let mut writer = vcf::io::Writer::new(&mut out);
        writer.write_header(&header)?;
        for result in reader.record_bufs(&header) {
            let record = result?;
            writer.write_variant_record(&header, &record)?;
        }
    }
    String::from_utf8(out).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn header_text_len(text: &str) -> usize {
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

fn comparable_naive_header(header: &str) -> String {
    header
        .lines()
        .filter(|line| {
            !line.starts_with("##bcftools_")
                || !(line.contains("Version=") || line.contains("Command="))
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn write_naive_output(args: &Args, buffer: &[u8]) -> io::Result<()> {
    match args.output_kind {
        OutputKind::VcfText => match &args.output {
            Some(path) => std::fs::write(path, buffer),
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
        OutputKind::Bcf | OutputKind::BcfUncompressed => unreachable!("validated by run_naive"),
    }
}

fn read_header(path: &Path) -> io::Result<vcf::Header> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(htslib_rs::bcf::io::Reader::new)?;
        return reader.read_header();
    }
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let normalized = NormalizeFileformat::new(BufReader::new(dec))?;
        let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
        return reader.read_header();
    }
    let f = File::open(path)?;
    let normalized = NormalizeFileformat::new(BufReader::new(f))?;
    let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
    reader.read_header()
}

fn for_each_record<F>(path: &Path, header: &vcf::Header, mut visit: F) -> io::Result<()>
where
    F: FnMut(&vcf::variant::RecordBuf) -> io::Result<()>,
{
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(htslib_rs::bcf::io::Reader::new)?;
        let _ = reader.read_header()?;
        for result in reader.record_bufs(header) {
            visit(&result?)?;
        }
        return Ok(());
    }
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let normalized = NormalizeFileformat::new(BufReader::new(dec))?;
        let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
        let _ = reader.read_header()?;
        for result in reader.records() {
            let record = result?;
            let buf = vcf::variant::RecordBuf::try_from_variant_record(header, &record)?;
            visit(&buf)?;
        }
        return Ok(());
    }
    let f = File::open(path)?;
    let normalized = NormalizeFileformat::new(BufReader::new(f))?;
    let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
    let _ = reader.read_header()?;
    for result in reader.records() {
        let record = result?;
        let buf = vcf::variant::RecordBuf::try_from_variant_record(header, &record)?;
        visit(&buf)?;
    }
    Ok(())
}

fn check_sample_columns(
    expected: &vcf::Header,
    actual: &vcf::Header,
    path: &Path,
) -> io::Result<()> {
    let exp: Vec<&str> = expected.sample_names().iter().map(|s| s.as_str()).collect();
    let act: Vec<&str> = actual.sample_names().iter().map(|s| s.as_str()).collect();
    if exp != act {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Different sample columns in {}: expected {:?}, got {:?}",
                path.display(),
                exp,
                act
            ),
        ));
    }
    Ok(())
}

fn drop_genotype_columns(header: &vcf::Header) -> vcf::Header {
    let mut text = Vec::new();
    vcf::io::Writer::new(&mut text)
        .write_header(header)
        .expect("serialize header");
    let text = String::from_utf8(text).expect("header is utf-8");
    let mut out = String::new();
    for line in text.split_inclusive('\n') {
        if line.starts_with("##FORMAT=") {
            continue;
        }
        if line.starts_with("#CHROM\t") {
            let trimmed = line.trim_end_matches('\n');
            let parts: Vec<&str> = trimmed.split('\t').collect();
            let kept: Vec<&str> = parts.iter().take(8).copied().collect();
            out.push_str(&kept.join("\t"));
            out.push('\n');
            continue;
        }
        out.push_str(line);
    }
    out.parse().expect("re-parse drop-genotypes header")
}

fn strip_format_and_samples(record: &vcf::variant::RecordBuf) -> vcf::variant::RecordBuf {
    let mut new = record.clone();
    *new.samples_mut() = vcf::variant::record_buf::Samples::default();
    new
}

trait ConcatWriter {
    fn write_record(&mut self, record: &vcf::variant::RecordBuf) -> io::Result<()>;
    fn finish(self: Box<Self>) -> io::Result<()>;
}

struct VcfTextWriter {
    writer: vcf::io::Writer<Box<dyn Write>>,
    header: vcf::Header,
}

impl VcfTextWriter {
    fn new_file(path: &Path, header: &vcf::Header, no_version: bool) -> io::Result<Self> {
        let inner: Box<dyn Write> = Box::new(File::create(path)?);
        Self::new_inner(inner, header, no_version)
    }

    fn new_stdout(header: &vcf::Header, no_version: bool) -> io::Result<Self> {
        let inner: Box<dyn Write> = Box::new(io::stdout().lock());
        Self::new_inner(inner, header, no_version)
    }

    fn new_inner(
        inner: Box<dyn Write>,
        header: &vcf::Header,
        no_version: bool,
    ) -> io::Result<Self> {
        let _ = no_version; // header-version line is not yet appended.
        let mut writer = vcf::io::Writer::new(inner);
        writer.write_header(header)?;
        Ok(Self {
            writer,
            header: header.clone(),
        })
    }
}

impl ConcatWriter for VcfTextWriter {
    fn write_record(&mut self, record: &vcf::variant::RecordBuf) -> io::Result<()> {
        self.writer.write_variant_record(&self.header, record)
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        Ok(())
    }
}

enum VcfGzWriter {
    Single {
        writer: vcf::io::Writer<htslib_rs::bgzf::io::Writer<Box<dyn Write>>>,
        header: vcf::Header,
    },
    ThreadedFile {
        writer: vcf::io::Writer<htslib_rs::bgzf::io::MultithreadedWriter<File>>,
        header: vcf::Header,
    },
}

impl VcfGzWriter {
    fn new(
        output: Option<&Path>,
        header: &vcf::Header,
        no_version: bool,
        thread_count: Option<NonZero<usize>>,
    ) -> io::Result<Self> {
        let _ = no_version;
        if let (Some(path), Some(thread_count)) = (output, thread_count) {
            let file = File::create(path)?;
            let bgzf =
                htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
            let mut writer = vcf::io::Writer::new(bgzf);
            writer.write_header(header)?;
            return Ok(Self::ThreadedFile {
                writer,
                header: header.clone(),
            });
        }
        let inner: Box<dyn Write> = match output {
            Some(path) => Box::new(File::create(path)?),
            None => Box::new(io::stdout().lock()),
        };
        let bgzf = htslib_rs::bgzf::io::Writer::new(inner);
        let mut writer = vcf::io::Writer::new(bgzf);
        writer.write_header(header)?;
        Ok(Self::Single {
            writer,
            header: header.clone(),
        })
    }
}

impl ConcatWriter for VcfGzWriter {
    fn write_record(&mut self, record: &vcf::variant::RecordBuf) -> io::Result<()> {
        match self {
            Self::Single { writer, header } => writer.write_variant_record(header, record),
            Self::ThreadedFile { writer, header } => writer.write_variant_record(header, record),
        }
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        match *self {
            Self::Single { writer, .. } => {
                let bgzf = writer.into_inner();
                bgzf.finish()?;
            }
            Self::ThreadedFile { writer, .. } => {
                let mut bgzf = writer.into_inner();
                let _file = bgzf.finish()?;
            }
        }
        Ok(())
    }
}

enum BcfWriter {
    Single {
        writer: htslib_rs::bcf::io::Writer<htslib_rs::bgzf::io::Writer<Box<dyn Write>>>,
        header: vcf::Header,
    },
    ThreadedFile {
        writer: htslib_rs::bcf::io::Writer<htslib_rs::bgzf::io::MultithreadedWriter<File>>,
        header: vcf::Header,
    },
}

impl BcfWriter {
    fn new(
        output: Option<&Path>,
        header: &vcf::Header,
        thread_count: Option<NonZero<usize>>,
    ) -> io::Result<Self> {
        if let (Some(path), Some(thread_count)) = (output, thread_count) {
            let file = File::create(path)?;
            let bgzf =
                htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
            let mut writer = htslib_rs::bcf::io::Writer::from(bgzf);
            writer.write_variant_header(header)?;
            return Ok(Self::ThreadedFile {
                writer,
                header: header.clone(),
            });
        }
        let inner: Box<dyn Write> = match output {
            Some(path) => Box::new(File::create(path)?),
            None => Box::new(io::stdout().lock()),
        };
        let mut writer = htslib_rs::bcf::io::Writer::new(inner);
        writer.write_variant_header(header)?;
        Ok(Self::Single {
            writer,
            header: header.clone(),
        })
    }
}

impl ConcatWriter for BcfWriter {
    fn write_record(&mut self, record: &vcf::variant::RecordBuf) -> io::Result<()> {
        match self {
            Self::Single { writer, header } => writer.write_variant_record(header, record),
            Self::ThreadedFile { writer, header } => writer.write_variant_record(header, record),
        }
    }

    fn finish(self: Box<Self>) -> io::Result<()> {
        match *self {
            Self::Single { mut writer, .. } => writer.try_finish(),
            Self::ThreadedFile { writer, .. } => {
                let mut bgzf = writer.into_inner();
                let _file = bgzf.finish()?;
                Ok(())
            }
        }
    }
}

#[derive(Default)]
struct DupSet {
    seen: std::collections::HashSet<DupKey>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RecordSpan {
    chrom: String,
    start: i64,
    end: i64,
}

impl RecordSpan {
    fn from(record: &vcf::variant::RecordBuf) -> Self {
        let start = record
            .variant_start()
            .map(|p| {
                let raw: usize = p.get();
                raw as i64
            })
            .unwrap_or(0);
        let ref_len = record.reference_bases().len().max(1) as i64;
        Self {
            chrom: record.reference_sequence_name().to_owned(),
            start,
            end: start + ref_len - 1,
        }
    }

    fn overlaps(&self, other: &Self) -> bool {
        self.chrom == other.chrom && other.start <= self.end
    }
}

impl DupSet {
    fn is_dup(&mut self, record: &vcf::variant::RecordBuf, mode: DupMode) -> bool {
        let key = DupKey::from(record, mode);
        !self.seen.insert(key)
    }
}

#[derive(Hash, PartialEq, Eq, Clone)]
struct DupKey {
    chrom: String,
    pos: i64,
    rest: String,
}

impl DupKey {
    fn from(record: &vcf::variant::RecordBuf, mode: DupMode) -> Self {
        let chrom = record.reference_sequence_name().to_owned();
        let pos = record
            .variant_start()
            .map(|p| {
                let raw: usize = p.get();
                raw as i64
            })
            .unwrap_or(0);
        let reference = record.reference_bases();
        let alts: Vec<&str> = record
            .alternate_bases()
            .as_ref()
            .iter()
            .map(|a| a.as_str())
            .collect();
        let rest = match mode {
            DupMode::Exact => format!("{}:{}", reference, alts.join(",")),
            DupMode::All => String::new(),
            DupMode::Snps | DupMode::Indels | DupMode::Both => {
                let kind = classify_kind(reference, &alts);
                format!("{kind:?}")
            }
        };
        Self { chrom, pos, rest }
    }
}

#[derive(Debug, PartialEq, Eq, Hash, Clone, Copy)]
enum VariantKind {
    Snp,
    Indel,
    Other,
}

fn classify_kind(reference: &str, alts: &[&str]) -> VariantKind {
    let mut has_snp = false;
    let mut has_indel = false;
    for alt in alts {
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

fn record_in_regions(
    record: &vcf::variant::RecordBuf,
    regions: &[RegionSpec],
    overlap_mode: u8,
) -> bool {
    if regions.is_empty() {
        return true;
    }
    let chrom = record.reference_sequence_name();
    let pos = record
        .variant_start()
        .map(|p| {
            let raw: usize = p.get();
            raw as i64
        })
        .unwrap_or(-1);
    let end = if overlap_mode == 0 {
        pos
    } else {
        let ref_len = record.reference_bases().len().max(1) as i64;
        pos + ref_len - 1
    };
    regions.iter().any(|spec| {
        spec.contig == chrom
            && match overlap_mode {
                0 => {
                    spec.start.map(|s| pos >= s).unwrap_or(true)
                        && spec.end.map(|e| pos <= e).unwrap_or(true)
                }
                _ => {
                    let region_start = spec.start.unwrap_or(i64::MIN);
                    let region_end = spec.end.unwrap_or(i64::MAX);
                    pos <= region_end && end >= region_start
                }
            }
    })
}

fn parse_regions_overlap(raw: &str) -> Result<u8, ParseOutcome> {
    parse_overlap_option(raw).ok_or_else(|| {
        ParseOutcome::Error(format!(
            "Could not parse --regions-overlap {raw} (expected 0|1|2, pos|record|variant)"
        ))
    })
}

fn parse_output_kind(raw: &str) -> Result<OutputKind, ParseOutcome> {
    OutputKind::parse(raw)
        .ok_or_else(|| ParseOutcome::Error(format!("The output type \"{raw}\" not recognised")))
}

fn parse_dup_mode(raw: &str) -> Result<DupMode, ParseOutcome> {
    DupMode::parse(raw).ok_or_else(|| {
        ParseOutcome::Error(format!(
            "Unrecognized --rm-dups value: {raw} (use snps|indels|both|all|exact)"
        ))
    })
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
        OutputKind::Bcf | OutputKind::BcfUncompressed => VariantOutputFormat::Bcf,
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
        .map(|value| value.to_string_lossy().into_owned())
        .ok_or_else(|| ParseOutcome::Error(format!("missing argument for {name}")))
}

fn value_after_equals(raw: &str) -> &str {
    raw.split_once('=').map(|(_, v)| v).unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_dup_mode_recognizes_all_aliases() {
        assert_eq!(DupMode::parse("snps"), Some(DupMode::Snps));
        assert_eq!(DupMode::parse("indels"), Some(DupMode::Indels));
        assert_eq!(DupMode::parse("both"), Some(DupMode::Both));
        assert_eq!(DupMode::parse("all"), Some(DupMode::All));
        assert_eq!(DupMode::parse("exact"), Some(DupMode::Exact));
        assert_eq!(DupMode::parse("nope"), None);
    }

    #[test]
    fn output_kind_for_path_recognizes_extensions() {
        assert_eq!(
            OutputKind::for_path(Path::new("a.vcf")),
            OutputKind::VcfText
        );
        assert_eq!(
            OutputKind::for_path(Path::new("a.vcf.gz")),
            OutputKind::VcfGz
        );
        assert_eq!(OutputKind::for_path(Path::new("a.bcf")), OutputKind::Bcf);
        assert_eq!(
            OutputKind::for_path(Path::new("a.unknown")),
            OutputKind::VcfText
        );
    }

    #[test]
    fn classify_kind_distinguishes_snp_indel() {
        assert_eq!(classify_kind("A", &["C"]), VariantKind::Snp);
        assert_eq!(classify_kind("A", &["AT"]), VariantKind::Indel);
        assert_eq!(classify_kind("AT", &["A"]), VariantKind::Indel);
        assert_eq!(classify_kind("A", &["C", "AT"]), VariantKind::Other);
    }
}
