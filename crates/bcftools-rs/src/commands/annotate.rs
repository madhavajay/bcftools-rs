//! Focused `bcftools annotate` implementation (upstream `vcfannotate.c`).
//!
//! This first local slice supports chromosome renaming with `--rename-chrs`.
//! Full annotation transfer/removal and overlap logic remain tracked in
//! `TODO.md`.

use std::collections::HashMap;
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::diagnostics::fmt_etag;
use crate::vcf_compat::normalize_vcf_text;

const USAGE: &str = "\n\
About:   Annotate and edit VCF/BCF files.\n\
Usage:   bcftools annotate [OPTIONS] <in.vcf.gz>\n\
\n\
Options:\n\
        --rename-chrs FILE           Rename chromosomes according to a two-column map\n\
    -o, --output FILE                Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9]   u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
        --no-version                 Accepted for command-shape compatibility\n\
\n";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    rename_chrs: PathBuf,
    output: Option<PathBuf>,
    output_kind: OutputKind,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
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
                eprintln!("{}", fmt_etag("main_vcfannotate", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfannotate", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut input = None;
    let mut rename_chrs = None;
    let mut output = None;
    let mut output_kind = OutputKind::VcfText;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "--rename-chrs" => {
                rename_chrs = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--no-version" => {}
            _ if raw.starts_with("--rename-chrs=") => {
                rename_chrs = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!(
                    "unrecognized option '{raw}' in this local annotate slice"
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
    let rename_chrs =
        rename_chrs.ok_or_else(|| ParseOutcome::Error("expected --rename-chrs FILE".into()))?;
    Ok(Args {
        input,
        rename_chrs,
        output,
        output_kind,
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

fn run(args: &Args) -> io::Result<()> {
    let map = read_rename_map(&args.rename_chrs)?;
    let mut text = read_vcf_text(&args.input)?;
    rename_vcf_text(&mut text, &map);
    write_output(text.as_bytes(), args)
}

fn read_rename_map(path: &Path) -> io::Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for (i, line) in fs::read_to_string(path)?.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid rename-chrs line {}: {line}", i + 1),
            ));
        }
        map.insert(fields[0].to_owned(), fields[1].to_owned());
    }
    Ok(map)
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
        ".bcftools-rs-annotate-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn rename_vcf_text(text: &mut String, map: &HashMap<String, String>) {
    let mut out = String::with_capacity(text.len());
    for line in text.lines() {
        if line.starts_with("##contig=<") {
            out.push_str(&rename_contig_header(line, map));
        } else if line.starts_with('#') {
            out.push_str(line);
        } else {
            out.push_str(&rename_record_chrom(line, map));
        }
        out.push('\n');
    }
    *text = out;
}

fn rename_contig_header(line: &str, map: &HashMap<String, String>) -> String {
    let Some(id_start) = line.find("ID=") else {
        return line.to_owned();
    };
    let value_start = id_start + 3;
    let value_end = line[value_start..]
        .find([',', '>'])
        .map(|idx| value_start + idx)
        .unwrap_or(line.len());
    let old = &line[value_start..value_end];
    let Some(new) = map.get(old) else {
        return line.to_owned();
    };
    let mut renamed = String::with_capacity(line.len() - old.len() + new.len());
    renamed.push_str(&line[..value_start]);
    renamed.push_str(new);
    renamed.push_str(&line[value_end..]);
    renamed
}

fn rename_record_chrom(line: &str, map: &HashMap<String, String>) -> String {
    let Some(tab) = line.find('\t') else {
        return line.to_owned();
    };
    let old = &line[..tab];
    let Some(new) = map.get(old) else {
        return line.to_owned();
    };
    let mut renamed = String::with_capacity(line.len() - old.len() + new.len());
    renamed.push_str(new);
    renamed.push_str(&line[tab..]);
    renamed
}

fn write_output(bytes: &[u8], args: &Args) -> io::Result<()> {
    match &args.output {
        Some(path) if path != Path::new("-") => {
            let file = File::create(path)?;
            write_to(bytes, args.output_kind, file)
        }
        _ => {
            let stdout = io::stdout();
            let out = stdout.lock();
            write_to(bytes, args.output_kind, out)
        }
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
    fn renames_contig_headers_and_records() {
        let mut map = HashMap::new();
        map.insert("1".to_owned(), "chr1".to_owned());
        let mut text = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t2\t.\tA\tC\t.\tPASS\t.\n"
            .to_owned();

        rename_vcf_text(&mut text, &map);
        assert!(text.contains("##contig=<ID=chr1,length=10>"));
        assert!(text.contains("chr1\t2\t.\tA\tC\t.\tPASS\t."));
    }
}
