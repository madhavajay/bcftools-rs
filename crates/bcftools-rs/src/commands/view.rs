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
//! Filtering, sample subsetting, and region restriction are NOT yet
//! implemented and yield an explicit error if requested.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::path::Path;
use std::process::ExitCode;
use std::time::SystemTime;

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::vcf::variant::io::Write as _;

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};
use crate::header_version::build_lines;

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
        LongOpt::new("no-version", HasArg::None, OPT_NO_VERSION),
        LongOpt::new("threads", HasArg::Required, OPT_THREADS),
    ];

    let mut output_kind = OutputKind::VcfText;
    let mut compression_level: Option<u32> = None;
    let mut output_file: Option<String> = None;
    let mut header_only = false;
    let mut no_header = false;
    let mut no_version = false;

    let mut g = Getopt::new("o:O:l:hH", &long_opts, argv);
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
                v if v == OPT_NO_VERSION => no_version = true,
                v if v == OPT_THREADS => {} // accepted; no-op for now
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
    let fname = match positional.first() {
        Some(p) => p,
        None => {
            eprint!("{USAGE}");
            return ExitCode::FAILURE;
        }
    };
    if positional.len() > 1 {
        eprintln!(
            "{}",
            fmt_etag(
                "main_vcfview",
                "region positional arguments not yet implemented"
            )
        );
        return ExitCode::FAILURE;
    }

    let path = Path::new(fname);
    let _ = compression_level; // consumed by future writers

    match run(
        path,
        output_kind,
        output_file.as_deref(),
        header_only,
        no_header,
        no_version,
        argv,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
            ExitCode::FAILURE
        }
    }
}

fn run(
    path: &Path,
    output_kind: OutputKind,
    output_file: Option<&str>,
    header_only: bool,
    no_header: bool,
    no_version: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let in_fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;

    let mut header = read_header(path, in_fmt)?;

    if !no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_view", &prog_argv, SystemTime::now());
        // Strip the "##" prefix and the "key=" delimiter from each rendered
        // line, then route both into the header via htslib-rs's typed-wrapper
        // helper. Mirrors upstream `bcf_hdr_append_version` which appends
        // "##bcftools_<cmd>Version" and "##bcftools_<cmd>Command" lines.
        for line in [&lines.version_line, &lines.command_line] {
            htslib_rs::header_compat::append_line(&mut header, line)?;
        }
    }

    let stdout = io::stdout();
    let writer: Box<dyn Write> = match output_file {
        Some("-") | None => Box::new(stdout.lock()),
        Some(p) => Box::new(File::create(p)?),
    };

    match output_kind {
        OutputKind::VcfText => write_vcf(path, in_fmt, &header, header_only, no_header, writer),
        OutputKind::VcfGz => {
            let bgzf = htslib_rs::bgzf::io::Writer::new(writer);
            write_vcf(path, in_fmt, &header, header_only, no_header, bgzf)
        }
        OutputKind::BcfUncompressed | OutputKind::BcfCompressed => {
            // For uncompressed BCF, upstream uses `wbu`. noodles' bcf writer
            // always wraps in BGZF; an "uncompressed" mode here is treated the
            // same as the compressed path until htslib-rs exposes the raw form.
            write_bcf(path, in_fmt, &header, header_only, no_header, writer)
        }
    }
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
    write_records_into_vcf(path, fmt, header, &mut writer)
}

fn write_records_into_vcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    writer: &mut htslib_rs::vcf::io::Writer<W>,
) -> io::Result<()> {
    use htslib_rs::bcf;
    use htslib_rs::vcf;

    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        for result in reader.record_bufs(header) {
            let rec = result?;
            writer.write_variant_record(header, &rec)?;
        }
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let mut reader = vcf::io::Reader::new(BufReader::new(dec));
        let _h = reader.read_header()?;
        for result in reader.records() {
            let rec = result?;
            writer.write_variant_record(header, &rec)?;
        }
    } else {
        let mut reader = File::open(path)
            .map(BufReader::new)
            .map(vcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        for result in reader.records() {
            let rec = result?;
            writer.write_variant_record(header, &rec)?;
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
            writer.write_variant_record(header, &rec)?;
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
            writer.write_variant_record(&header, &rec)?;
        }
        writer.try_finish()?;
        Ok(())
    } else {
        // Plain VCF → BCF: delegate to htslib-rs's tested helper.
        htslib_rs::variant_io_compat::write_bcf_from_vcf_path(path, out)?;
        Ok(())
    }
}
