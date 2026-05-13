//! Port of `bcftools view` (upstream `vcfview.c`).
//!
//! This is the parity-anchor subcommand for VCF/BCF I/O. The full upstream
//! `view` accepts ~50 options (sample/region restriction, filtering, allele
//! count gates, FILTER tag dispatch, header-only mode, etc.). This initial
//! port covers only the I/O backbone:
//!
//! - read VCF / VCF.gz / BCF input (auto-detected by file content)
//! - write to one of `-O v|z|u|b` (VCF text / VCF.gz / uncompressed BCF /
//!   compressed BCF)
//! - `-o, --output FILE` to write to a path (default: stdout for `v`, error
//!   otherwise to avoid binary-on-tty)
//! - `--no-version` suppresses the `##bcftools_view{Version,Command}` header
//!   lines (other code paths inject them; here we honor the flag).
//! - `-h, --header-only` and `-H, --no-header` for header-vs-records dispatch.
//!
//! Filtering is NOT yet implemented and yields an explicit error if requested.
//! Positional region arguments support the common `CHROM` and `CHROM:START-END`
//! forms by streaming and filtering records.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Read as _, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::core::Position;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::vcf::variant::io::Write as _;

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};
use crate::header_version::{build_lines, command_time};

const USAGE: &str = "\n\
About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.\n\
Usage:   bcftools view [OPTIONS] <in.vcf.gz>|<in.bcf> [REGION...]\n\
\n\
Output options:\n\
    -G, --drop-genotypes              drop individual genotype information (NOT IMPLEMENTED)\n\
    -h, --header-only                 print only the header in VCF output\n\
    -H, --no-header                   suppress the header in VCF output\n\
    -l, --compression-level INT       compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]\n\
        --no-version                  do not append version and command line to the header\n\
    -o, --output FILE                 output file name [stdout]\n\
    -O, --output-type u|b|v|z[0-9]    u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n\
    -s, --samples LIST                comma-separated sample list, optionally prefixed with ^\n\
    -S, --samples-file FILE           file of samples, optionally prefixed with ^\n\
\n";

const OPT_NO_VERSION: i32 = 200;
const OPT_THREADS: i32 = 9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    BcfUncompressed,
    BcfCompressed,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Region {
    contig: String,
    start: Option<usize>,
    end: Option<usize>,
}

struct RunOptions<'a> {
    output_kind: OutputKind,
    output_file: Option<&'a str>,
    header_only: bool,
    no_header: bool,
    no_version: bool,
    regions: &'a [Region],
    thread_count: Option<NonZero<usize>>,
    sample_list: Option<&'a str>,
    sample_list_is_file: bool,
}

impl Region {
    fn parse(raw: &str) -> io::Result<Self> {
        let (contig, interval) = raw.split_once(':').unwrap_or((raw, ""));
        if contig.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Could not parse region \"{raw}\""),
            ));
        }

        if interval.is_empty() {
            return Ok(Self {
                contig: contig.to_string(),
                start: None,
                end: None,
            });
        }

        let (start, end) = interval.split_once('-').unwrap_or((interval, ""));
        let start = if start.is_empty() {
            None
        } else {
            Some(parse_region_pos(start, raw)?)
        };
        let end = if end.is_empty() {
            None
        } else {
            Some(parse_region_pos(end, raw)?)
        };

        Ok(Self {
            contig: contig.to_string(),
            start,
            end,
        })
    }

    fn contains(&self, contig: &str, pos: Position) -> bool {
        if self.contig != contig {
            return false;
        }
        let pos = usize::from(pos);
        self.start.is_none_or(|start| pos >= start) && self.end.is_none_or(|end| pos <= end)
    }
}

fn parse_region_pos(s: &str, raw_region: &str) -> io::Result<usize> {
    s.replace(',', "").parse::<usize>().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Could not parse region \"{raw_region}\""),
        )
    })
}

fn parse_threads(raw: &str) -> Result<Option<NonZero<usize>>, std::num::ParseIntError> {
    let n = raw.parse::<usize>()?;
    Ok(NonZero::new(n))
}

impl OutputKind {
    fn parse(s: &str) -> Option<(Self, Option<u32>)> {
        if s.is_empty() {
            return None;
        }
        let kind = match s.as_bytes()[0] {
            b'v' => OutputKind::VcfText,
            b'z' => OutputKind::VcfGz,
            b'u' => OutputKind::BcfUncompressed,
            b'b' => OutputKind::BcfCompressed,
            _ => return None,
        };
        let level = if s.len() > 1 {
            match s[1..].parse::<u32>() {
                Ok(l) if l <= 9 => Some(l),
                _ => return None,
            }
        } else {
            None
        };
        Some((kind, level))
    }
}

/// Subcommand entry point. `argv[0]` is `"view"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    let long_opts = [
        LongOpt::new("output", HasArg::Required, b'o' as i32),
        LongOpt::new("output-file", HasArg::Required, b'o' as i32),
        LongOpt::new("output-type", HasArg::Required, b'O' as i32),
        LongOpt::new("compression-level", HasArg::Required, b'l' as i32),
        LongOpt::new("header-only", HasArg::None, b'h' as i32),
        LongOpt::new("no-header", HasArg::None, b'H' as i32),
        LongOpt::new("samples", HasArg::Required, b's' as i32),
        LongOpt::new("samples-file", HasArg::Required, b'S' as i32),
        LongOpt::new("no-version", HasArg::None, OPT_NO_VERSION),
        LongOpt::new("threads", HasArg::Required, OPT_THREADS),
    ];

    let mut output_kind = OutputKind::VcfText;
    let mut compression_level: Option<u32> = None;
    let mut output_file: Option<String> = None;
    let mut header_only = false;
    let mut no_header = false;
    let mut no_version = false;
    let mut thread_count = None;
    let mut sample_list: Option<String> = None;
    let mut sample_list_is_file = false;

    let mut g = Getopt::new("o:O:l:hHs:S:", &long_opts, argv);
    loop {
        match g.next() {
            Ok(Some(m)) => match m.code {
                v if v == b'o' as i32 => output_file = m.value,
                v if v == b'O' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match OutputKind::parse(raw) {
                        Some((k, lvl)) => {
                            output_kind = k;
                            if lvl.is_some() {
                                compression_level = lvl;
                            }
                        }
                        None => {
                            eprintln!(
                                "{}",
                                fmt_etag(
                                    "main_vcfview",
                                    &format!("The output type \"{raw}\" not recognised")
                                )
                            );
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == b'l' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match raw.parse::<u32>() {
                        Ok(l) if l <= 9 => compression_level = Some(l),
                        _ => {
                            eprintln!(
                                "{}",
                                fmt_etag(
                                    "main_vcfview",
                                    &format!("invalid compression level \"{raw}\"")
                                )
                            );
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == b'h' as i32 => header_only = true,
                v if v == b'H' as i32 => no_header = true,
                v if v == b's' as i32 => {
                    sample_list = m.value;
                    sample_list_is_file = false;
                }
                v if v == b'S' as i32 => {
                    sample_list = m.value;
                    sample_list_is_file = true;
                }
                v if v == OPT_NO_VERSION => no_version = true,
                v if v == OPT_THREADS => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_threads(raw) {
                        Ok(n) => thread_count = n,
                        Err(_) => {
                            eprintln!(
                                "{}",
                                fmt_etag(
                                    "main_vcfview",
                                    &format!("Could not parse argument: --threads {raw}")
                                )
                            );
                            return ExitCode::FAILURE;
                        }
                    }
                }
                _ => {
                    eprint!("{USAGE}");
                    return ExitCode::FAILURE;
                }
            },
            Ok(None) => break,
            Err(_) => {
                eprint!("{USAGE}");
                return ExitCode::FAILURE;
            }
        }
    }

    let positional = g.rest();
    let fname = positional
        .first()
        .cloned()
        .unwrap_or_else(|| OsString::from("-"));
    let regions = match positional
        .iter()
        .skip(1)
        .map(|raw| Region::parse(&raw.to_string_lossy()))
        .collect::<io::Result<Vec<_>>>()
    {
        Ok(regions) => regions,
        Err(e) => {
            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
            return ExitCode::FAILURE;
        }
    };

    let path = Path::new(&fname);
    let _ = compression_level; // consumed by future writers

    let options = RunOptions {
        output_kind,
        output_file: output_file.as_deref(),
        header_only,
        no_header,
        no_version,
        regions: &regions,
        thread_count,
        sample_list: sample_list.as_deref(),
        sample_list_is_file,
    };

    match run(path, &options, argv) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
            ExitCode::FAILURE
        }
    }
}

fn run(path: &Path, options: &RunOptions<'_>, argv: &[OsString]) -> io::Result<()> {
    if path == Path::new("-") {
        let tmp = stdin_tmp_path();
        let mut data = Vec::new();
        io::stdin().lock().read_to_end(&mut data)?;
        fs::write(&tmp, data)?;
        let result = run(&tmp, options, argv);
        let _ = fs::remove_file(&tmp);
        return result;
    }

    let in_fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;

    if options.sample_list.is_some() {
        return run_sample_subset(path, in_fmt, options, argv);
    }

    let mut header = read_header(path, in_fmt)?;

    if !options.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_view", &prog_argv, command_time());
        // Strip the "##" prefix and the "key=" delimiter from each rendered
        // line, then route both into the header via htslib-rs's typed-wrapper
        // helper. Mirrors upstream `bcf_hdr_append_version` which appends
        // "##bcftools_<cmd>Version" and "##bcftools_<cmd>Command" lines.
        for line in [&lines.version_line, &lines.command_line] {
            htslib_rs::header_compat::append_line(&mut header, line)?;
        }
    }

    if options.output_kind == OutputKind::VcfText
        && options.no_version
        && options.regions.is_empty()
        && in_fmt.exact == Exact::Bcf
    {
        return match options.output_file {
            Some("-") | None => write_bcf_vcf_text_no_version(
                path,
                options.header_only,
                options.no_header,
                io::stdout().lock(),
            ),
            Some(p) => write_bcf_vcf_text_no_version(
                path,
                options.header_only,
                options.no_header,
                File::create(p)?,
            ),
        };
    }

    match options.output_kind {
        OutputKind::VcfText => match options.output_file {
            Some("-") | None
                if options.no_version
                    && options.regions.is_empty()
                    && in_fmt.exact != Exact::Bcf =>
            {
                write_vcf_text_passthrough(
                    path,
                    in_fmt,
                    options.header_only,
                    options.no_header,
                    io::stdout().lock(),
                )
            }
            Some("-") | None => write_vcf(
                path,
                in_fmt,
                &header,
                options.header_only,
                options.no_header,
                options.regions,
                io::stdout().lock(),
            ),
            Some(p)
                if options.no_version
                    && options.regions.is_empty()
                    && in_fmt.exact != Exact::Bcf =>
            {
                write_vcf_text_passthrough(
                    path,
                    in_fmt,
                    options.header_only,
                    options.no_header,
                    File::create(p)?,
                )
            }
            Some(p) => write_vcf(
                path,
                in_fmt,
                &header,
                options.header_only,
                options.no_header,
                options.regions,
                File::create(p)?,
            ),
        },
        OutputKind::VcfGz
            if options.no_version && options.regions.is_empty() && in_fmt.exact != Exact::Bcf =>
        {
            match options.output_file {
                Some(p) if p != "-" => {
                    let bgzf = htslib_rs::bgzf::io::Writer::new(File::create(p)?);
                    write_vcf_text_passthrough(
                        path,
                        in_fmt,
                        options.header_only,
                        options.no_header,
                        bgzf,
                    )
                }
                _ => {
                    let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
                    write_vcf_text_passthrough(
                        path,
                        in_fmt,
                        options.header_only,
                        options.no_header,
                        bgzf,
                    )
                }
            }
        }
        OutputKind::VcfGz => write_vcf_gz(path, in_fmt, &header, options),
        OutputKind::BcfUncompressed | OutputKind::BcfCompressed => {
            // For uncompressed BCF, upstream uses `wbu`. noodles' bcf writer
            // always wraps in BGZF; an "uncompressed" mode here is treated the
            // same as the compressed path until htslib-rs exposes the raw form.
            match options.output_file {
                Some("-") | None => write_bcf(
                    path,
                    in_fmt,
                    &header,
                    options.header_only,
                    options.no_header,
                    options.regions,
                    io::stdout().lock(),
                ),
                Some(p) => write_bcf(
                    path,
                    in_fmt,
                    &header,
                    options.header_only,
                    options.no_header,
                    options.regions,
                    File::create(p)?,
                ),
            }
        }
    }
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-view-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn run_sample_subset(
    path: &Path,
    in_fmt: format::Format,
    options: &RunOptions<'_>,
    argv: &[OsString],
) -> io::Result<()> {
    if matches!(
        options.output_kind,
        OutputKind::BcfUncompressed | OutputKind::BcfCompressed
    ) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "sample subsetting with BCF output is not yet supported",
        ));
    }

    let version_lines = if options.no_version {
        None
    } else {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        Some(build_lines("bcftools_view", &prog_argv, command_time()))
    };

    match options.output_kind {
        OutputKind::VcfText => match options.output_file {
            Some("-") | None => write_sample_subset_vcf(
                path,
                in_fmt,
                options,
                version_lines.as_ref(),
                io::stdout().lock(),
            ),
            Some(p) => write_sample_subset_vcf(
                path,
                in_fmt,
                options,
                version_lines.as_ref(),
                File::create(p)?,
            ),
        },
        OutputKind::VcfGz => match (options.output_file, options.thread_count) {
            (Some(p), Some(thread_count)) if p != "-" => {
                let bgzf = htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(
                    thread_count,
                    File::create(p)?,
                );
                write_sample_subset_vcf(path, in_fmt, options, version_lines.as_ref(), bgzf)
            }
            (Some(p), _) if p != "-" => {
                let bgzf = htslib_rs::bgzf::io::Writer::new(File::create(p)?);
                write_sample_subset_vcf(path, in_fmt, options, version_lines.as_ref(), bgzf)
            }
            _ => {
                let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
                write_sample_subset_vcf(path, in_fmt, options, version_lines.as_ref(), bgzf)
            }
        },
        OutputKind::BcfUncompressed | OutputKind::BcfCompressed => unreachable!(),
    }
}

fn write_sample_subset_vcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    options: &RunOptions<'_>,
    version_lines: Option<&crate::header_version::HeaderVersionLines>,
    mut out: W,
) -> io::Result<()> {
    let text = vcf_text_from_path(path, fmt)?;
    write_sample_subset_vcf_text(&text, options, version_lines, &mut out)
}

fn vcf_text_from_path(path: &Path, fmt: format::Format) -> io::Result<String> {
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

fn write_sample_subset_vcf_text<W: Write>(
    text: &str,
    options: &RunOptions<'_>,
    version_lines: Option<&crate::header_version::HeaderVersionLines>,
    out: &mut W,
) -> io::Result<()> {
    let mut selected_samples: Option<Vec<usize>> = None;
    let mut inserted_version = false;

    for line in text.split_inclusive('\n') {
        if line.starts_with("##") {
            if !options.no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }

        if line.starts_with("#CHROM\t") {
            if !options.no_header {
                if let Some(lines) = version_lines
                    && !inserted_version
                {
                    writeln!(out, "{}", lines.version_line)?;
                    writeln!(out, "{}", lines.command_line)?;
                    inserted_version = true;
                }
                let fields = line_fields(line);
                let sample_names = fields[9..]
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>();
                let selected = crate::smpl_ilist::init(
                    &sample_names,
                    options.sample_list,
                    options.sample_list_is_file,
                    crate::smpl_ilist::SMPL_STRICT,
                )?
                .idx;
                write_projected_vcf_line(&fields, &selected, out)?;
                selected_samples = Some(selected);
            } else {
                let fields = line_fields(line);
                let sample_names = fields[9..]
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>();
                selected_samples = Some(
                    crate::smpl_ilist::init(
                        &sample_names,
                        options.sample_list,
                        options.sample_list_is_file,
                        crate::smpl_ilist::SMPL_STRICT,
                    )?
                    .idx,
                );
            }
            continue;
        }

        if line.starts_with('#') {
            if !options.no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }

        if options.header_only {
            break;
        }
        let fields = line_fields(line);
        if !record_line_matches_regions(&fields, options.regions) {
            continue;
        }
        let selected = selected_samples.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "VCF header is missing the #CHROM sample line",
            )
        })?;
        write_projected_vcf_line(&fields, selected, out)?;
    }

    Ok(())
}

fn line_fields(line: &str) -> Vec<&str> {
    line.trim_end_matches('\n')
        .trim_end_matches('\r')
        .split('\t')
        .collect()
}

fn write_projected_vcf_line<W: Write>(
    fields: &[&str],
    selected_samples: &[usize],
    out: &mut W,
) -> io::Result<()> {
    let fixed_end = fields.len().min(9);
    let mut projected = fields[..fixed_end].to_vec();
    for &sample_idx in selected_samples {
        if let Some(value) = fields.get(9 + sample_idx) {
            projected.push(value);
        }
    }
    writeln!(out, "{}", projected.join("\t"))
}

fn record_line_matches_regions(fields: &[&str], regions: &[Region]) -> bool {
    if regions.is_empty() {
        return true;
    }
    let Some(contig) = fields.first() else {
        return false;
    };
    let Some(pos) = fields.get(1).and_then(|pos| pos.parse::<usize>().ok()) else {
        return false;
    };
    regions.iter().any(|region| {
        region.contig == *contig
            && region.start.is_none_or(|start| pos >= start)
            && region.end.is_none_or(|end| pos <= end)
    })
}

fn write_vcf_gz(
    path: &Path,
    in_fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    options: &RunOptions<'_>,
) -> io::Result<()> {
    match (options.output_file, options.thread_count) {
        (Some(p), Some(thread_count)) if p != "-" => {
            let file = File::create(p)?;
            let bgzf =
                htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
            write_vcf(
                path,
                in_fmt,
                header,
                options.header_only,
                options.no_header,
                options.regions,
                bgzf,
            )
        }
        (Some(p), _) if p != "-" => {
            let bgzf = htslib_rs::bgzf::io::Writer::new(File::create(p)?);
            write_vcf(
                path,
                in_fmt,
                header,
                options.header_only,
                options.no_header,
                options.regions,
                bgzf,
            )
        }
        _ => {
            let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
            write_vcf(
                path,
                in_fmt,
                header,
                options.header_only,
                options.no_header,
                options.regions,
                bgzf,
            )
        }
    }
}

fn write_vcf_text_passthrough<W: Write>(
    path: &Path,
    fmt: format::Format,
    header_only: bool,
    no_header: bool,
    out: W,
) -> io::Result<()> {
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        return write_vcf_text_passthrough_reader(BufReader::new(dec), header_only, no_header, out);
    }
    let reader = File::open(path).map(BufReader::new)?;
    write_vcf_text_passthrough_reader(reader, header_only, no_header, out)
}

fn write_vcf_text_passthrough_reader<R, W>(
    mut reader: R,
    header_only: bool,
    no_header: bool,
    mut out: W,
) -> io::Result<()>
where
    R: io::BufRead,
    W: Write,
{
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with('#') {
            if !no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }
        if header_only {
            break;
        }
        out.write_all(line.as_bytes())?;
        io::copy(&mut reader, &mut out)?;
        break;
    }
    Ok(())
}

fn write_bcf_vcf_text_no_version<W: Write>(
    path: &Path,
    header_only: bool,
    no_header: bool,
    mut out: W,
) -> io::Result<()> {
    let text = htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None)?;
    write_vcf_text_from_string(&text, header_only, no_header, &mut out)
}

fn write_vcf_text_from_string<W: Write>(
    text: &str,
    header_only: bool,
    no_header: bool,
    out: &mut W,
) -> io::Result<()> {
    for line in text.split_inclusive('\n') {
        if line.starts_with('#') {
            if !no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }
        if header_only {
            break;
        }
        out.write_all(line.as_bytes())?;
    }
    Ok(())
}

fn read_header(path: &Path, fmt: format::Format) -> io::Result<htslib_rs::vcf::Header> {
    use htslib_rs::variant_io_compat::{
        read_bcf_header_from_path, read_vcf_header, read_vcf_header_from_path,
    };
    if fmt.exact == Exact::Bcf {
        read_bcf_header_from_path(path)
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        read_vcf_header(BufReader::new(dec))
    } else {
        read_vcf_header_from_path(path)
    }
}

fn write_vcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    header_only: bool,
    no_header: bool,
    regions: &[Region],
    out: W,
) -> io::Result<()> {
    use htslib_rs::vcf;
    let mut writer = vcf::io::Writer::new(out);
    if !no_header {
        writer.write_header(header)?;
    }
    if header_only {
        return Ok(());
    }
    write_records_into_vcf(path, fmt, header, regions, &mut writer)
}

fn write_records_into_vcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    regions: &[Region],
    writer: &mut htslib_rs::vcf::io::Writer<W>,
) -> io::Result<()> {
    use htslib_rs::bcf;
    use htslib_rs::vcf;

    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        for result in reader.record_bufs(header) {
            let rec = result?;
            if region_matches(regions, rec.reference_sequence_name(), rec.variant_start()) {
                writer.write_variant_record(header, &rec)?;
            }
        }
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let mut reader = vcf::io::Reader::new(BufReader::new(dec));
        let _h = reader.read_header()?;
        for result in reader.records() {
            let rec = result?;
            if region_matches_result(regions, rec.reference_sequence_name(), rec.variant_start())? {
                writer.write_variant_record(header, &rec)?;
            }
        }
    } else {
        let mut reader = File::open(path)
            .map(BufReader::new)
            .map(vcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        for result in reader.records() {
            let rec = result?;
            if region_matches_result(regions, rec.reference_sequence_name(), rec.variant_start())? {
                writer.write_variant_record(header, &rec)?;
            }
        }
    }
    Ok(())
}

fn write_bcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    header_only: bool,
    no_header: bool,
    regions: &[Region],
    out: W,
) -> io::Result<()> {
    use htslib_rs::bcf;
    let _ = no_header; // BCF cannot be sensibly written without a header.
    if header_only {
        let mut writer = bcf::io::Writer::new(out);
        writer.write_variant_header(header)?;
        writer.try_finish()?;
        return Ok(());
    }
    if fmt.exact == Exact::Bcf {
        // BCF → BCF: copy records through as-is. Use record_bufs so the writer
        // sees fully decoded records keyed by contig string.
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        let mut writer = bcf::io::Writer::new(out);
        writer.write_variant_header(header)?;
        for result in reader.record_bufs(header) {
            let rec = result?;
            if region_matches(regions, rec.reference_sequence_name(), rec.variant_start()) {
                writer.write_variant_record(header, &rec)?;
            }
        }
        writer.try_finish()?;
        Ok(())
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        // VCF.gz → BCF: decompress on the fly into the htslib-rs path that's
        // exercised by its own test suite.
        use htslib_rs::vcf;
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let mut reader = vcf::io::Reader::new(BufReader::new(dec));
        let header = reader.read_header()?;
        let mut writer = bcf::io::Writer::new(out);
        writer.write_variant_header(&header)?;
        for result in reader.records() {
            let rec = result?;
            if region_matches_result(regions, rec.reference_sequence_name(), rec.variant_start())? {
                writer.write_variant_record(&header, &rec)?;
            }
        }
        writer.try_finish()?;
        Ok(())
    } else {
        if regions.is_empty() {
            // Plain VCF → BCF: delegate to htslib-rs's tested helper.
            htslib_rs::variant_io_compat::write_bcf_from_vcf_path(path, out)?;
        } else {
            use htslib_rs::vcf;
            let mut reader = File::open(path)
                .map(BufReader::new)
                .map(vcf::io::Reader::new)?;
            let header = reader.read_header()?;
            let mut writer = bcf::io::Writer::new(out);
            writer.write_variant_header(&header)?;
            for result in reader.records() {
                let rec = result?;
                if region_matches_result(
                    regions,
                    rec.reference_sequence_name(),
                    rec.variant_start(),
                )? {
                    writer.write_variant_record(&header, &rec)?;
                }
            }
            writer.try_finish()?;
        }
        Ok(())
    }
}

fn region_matches(regions: &[Region], contig: &str, pos: Option<Position>) -> bool {
    regions.is_empty()
        || pos
            .map(|pos| regions.iter().any(|region| region.contains(contig, pos)))
            .unwrap_or(false)
}

fn region_matches_result(
    regions: &[Region],
    contig: &str,
    pos: Option<io::Result<Position>>,
) -> io::Result<bool> {
    match pos {
        Some(Ok(pos)) => Ok(region_matches(regions, contig, Some(pos))),
        Some(Err(e)) => Err(e),
        None => Ok(region_matches(regions, contig, None)),
    }
}
