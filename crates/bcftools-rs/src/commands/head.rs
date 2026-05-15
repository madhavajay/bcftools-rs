//! Port of `bcftools head` (upstream `vcfhead.c`).
//!
//! Behavior:
//!
//! - `bcftools head FILE` — print the entire header.
//! - `-h N` — print only the first N header lines.
//! - `-n N` — also print the first N variant records.
//! - `-s N` — print N records starting *with* the `#CHROM` header line; this
//!   implies `-h 0` and ensures the `#CHROM` line is present in output even
//!   when `-h` would have truncated before it.
//! - `-v N` — verbosity passthrough to `htslib-rs::log_compat`.
//!
//! Stdin is supported when no file is supplied or when the file is `-`.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read as _, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::vcf::variant::io::Write as _;

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};
use crate::io::apply_verbosity;
use crate::vcf_compat::NormalizeFileformat;

const USAGE: &str = "\n\
About: Displays VCF/BCF headers and optionally the first few variant records\n\
Usage: bcftools head [OPTION]... [FILE]\n\
\n\
Options:\n\
  -h, --headers INT      Display INT header lines [all]\n\
  -n, --records INT      Display INT variant record lines [none]\n\
  -s, --samples INT      Display INT records starting with the #CHROM header line [none]\n\
  -v, --verbosity INT    Verbosity level\n\
\n";

/// Subcommand entry point. `argv[0]` is `"head"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    let long_opts = [
        LongOpt::new("headers", HasArg::Required, b'h' as i32),
        LongOpt::new("records", HasArg::Required, b'n' as i32),
        LongOpt::new("samples", HasArg::Required, b's' as i32),
        LongOpt::new("verbosity", HasArg::Required, b'v' as i32),
    ];

    let mut all_headers = true;
    let mut samples = false;
    let mut nheaders: u64 = 0;
    let mut nrecords: u64 = 0;

    let mut g = Getopt::new("h:n:s:v:", &long_opts, argv);
    loop {
        match g.next() {
            Ok(Some(m)) => match m.code as u8 as char {
                'v' => {
                    if apply_verbosity(m.value.as_deref().unwrap_or("0")).is_err() {
                        eprintln!(
                            "Could not parse argument: --verbosity {}",
                            m.value.as_deref().unwrap_or("")
                        );
                        return ExitCode::from(255);
                    }
                }
                'h' => {
                    all_headers = false;
                    nheaders = parse_u64(m.value.as_deref().unwrap_or("0"));
                }
                'n' => {
                    nrecords = parse_u64(m.value.as_deref().unwrap_or("0"));
                }
                's' => {
                    nrecords = parse_u64(m.value.as_deref().unwrap_or("0"));
                    samples = true;
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

    if samples && all_headers {
        all_headers = false;
    }

    let positional = g.rest();
    let fname = match positional.len() {
        0 => None,
        1 => Some(positional[0].as_os_str().to_string_lossy().into_owned()),
        _ => {
            eprint!("{USAGE}");
            return ExitCode::FAILURE;
        }
    };

    match run_input(fname.as_deref(), all_headers, nheaders, samples, nrecords) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{}", fmt_etag("main_vcfhead", &format!("{e}")));
            ExitCode::FAILURE
        }
    }
}

fn parse_u64(s: &str) -> u64 {
    // Match libc strtoull semantics used in upstream: leading garbage → 0,
    // negative → wrap to large positive (rare in practice), but for our tests
    // we just accept decimal and treat parse failure as 0.
    s.parse::<u64>().unwrap_or(0)
}

fn run_input(
    fname: Option<&str>,
    all_headers: bool,
    nheaders: u64,
    samples: bool,
    nrecords: u64,
) -> io::Result<()> {
    match fname {
        Some(name) if name != "-" => {
            let path = Path::new(name);
            run(path, all_headers, nheaders, samples, nrecords).map_err(|e| {
                io::Error::new(e.kind(), format!("Can't open \"{}\": {e}", path.display()))
            })
        }
        _ => {
            let mut data = Vec::new();
            io::stdin().lock().read_to_end(&mut data)?;
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "No input data",
                ));
            }
            let path = stdin_tmp_path();
            std::fs::write(&path, &data)?;
            let result = run(&path, all_headers, nheaders, samples, nrecords);
            let _ = std::fs::remove_file(&path);
            result
        }
    }
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-head-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn run(
    path: &Path,
    all_headers: bool,
    nheaders: u64,
    samples: bool,
    nrecords: u64,
) -> io::Result<()> {
    let format = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    let header_text = read_header_text(path, format)?;

    let mut stdout = io::stdout().lock();

    if all_headers {
        stdout.write_all(header_text.as_bytes())?;
    } else if nheaders > 0 || samples {
        let lines: Vec<&str> = header_text.lines().collect();
        let take = (nheaders as usize).min(lines.len());
        let mut samples_printed = false;
        for line in &lines[..take] {
            stdout.write_all(line.as_bytes())?;
            stdout.write_all(b"\n")?;
            if line.starts_with("#CHROM\t") {
                samples_printed = true;
            }
        }
        if samples && !samples_printed {
            for line in &lines[take..] {
                if line.starts_with("#CHROM\t") {
                    stdout.write_all(line.as_bytes())?;
                    stdout.write_all(b"\n")?;
                    break;
                }
            }
        }
    }

    if nrecords > 0 {
        write_n_records(path, format, &header_text, nrecords, &mut stdout)?;
    }
    Ok(())
}

/// Read the full header text, serialized to VCF.
fn read_header_text(path: &Path, fmt: format::Format) -> io::Result<String> {
    use htslib_rs::variant_io_compat::read_bcf_header_from_path;

    if fmt.exact != Exact::Bcf {
        // For VCF text, preserve the byte order of raw header lines. Writing a
        // structured header normalizes record categories and fails upstream
        // `test_vcf_head` byte parity.
        if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
            let f = File::open(path)?;
            let dec = flate2::read::MultiGzDecoder::new(f);
            return read_vcf_header_text(BufReader::new(dec));
        }
        return File::open(path)
            .map(BufReader::new)
            .and_then(read_vcf_header_text);
    }

    let header = read_bcf_header_from_path(path)?;
    let mut buf = Vec::new();
    htslib_rs::vcf::io::Writer::new(&mut buf).write_header(&header)?;
    String::from_utf8(buf)
        .map(reorder_bcf_header_text_for_head)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
}

fn reorder_bcf_header_text_for_head(header: String) -> String {
    let mut fileformat = Vec::new();
    let mut filters = Vec::new();
    let mut rest = Vec::new();
    let mut chrom = Vec::new();

    for line in header.split_inclusive('\n') {
        if line.starts_with("##fileformat=") {
            fileformat.push(line);
        } else if line.starts_with("##FILTER=") {
            filters.push(line);
        } else if line.starts_with("#CHROM\t") {
            chrom.push(line);
        } else {
            rest.push(line);
        }
    }

    let mut out = String::with_capacity(header.len());
    for group in [fileformat, filters, rest, chrom] {
        for line in group {
            out.push_str(line);
        }
    }
    out
}

fn read_vcf_header_text<R: BufRead>(mut reader: R) -> io::Result<String> {
    let mut out = String::new();
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 || !line.starts_with('#') {
            break;
        }
        out.push_str(&line);
    }
    Ok(out)
}

fn write_n_records<W: Write>(
    path: &Path,
    fmt: format::Format,
    _header_text: &str,
    nrecords: u64,
    out: &mut W,
) -> io::Result<()> {
    use htslib_rs::bcf;
    use htslib_rs::vcf;

    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let header = reader.read_header()?;
        let mut writer = vcf::io::Writer::new(out);
        for (i, result) in reader.record_bufs(&header).enumerate() {
            if i as u64 >= nrecords {
                break;
            }
            let rec = result?;
            writer.write_variant_record(&header, &rec)?;
        }
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let normalized = NormalizeFileformat::new(BufReader::new(dec))?;
        let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
        let header = reader.read_header()?;
        let mut writer = vcf::io::Writer::new(out);
        for (i, result) in reader.records().enumerate() {
            if i as u64 >= nrecords {
                break;
            }
            let rec = result?;
            writer.write_variant_record(&header, &rec)?;
        }
    } else {
        let f = File::open(path)?;
        let normalized = NormalizeFileformat::new(BufReader::new(f))?;
        let mut reader = vcf::io::Reader::new(BufReader::new(normalized));
        let header = reader.read_header()?;
        let mut writer = vcf::io::Writer::new(out);
        for (i, result) in reader.records().enumerate() {
            if i as u64 >= nrecords {
                break;
            }
            let rec = result?;
            writer.write_variant_record(&header, &rec)?;
        }
    }
    Ok(())
}
