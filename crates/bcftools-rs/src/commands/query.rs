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

use htslib_rs::format::{self, Compression, Exact};

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};

const USAGE: &str = "\n\
About:   Extract fields from VCF/BCF files and print sample lists.\n\
Usage:   bcftools query [OPTIONS] <in.vcf.gz>|<in.bcf>\n\
\n\
Options:\n\
    -f, --format STR                 format string\n\
    -l, --list-samples               print sample names and exit\n\
\n";

#[derive(Debug, Clone, PartialEq, Eq)]
struct Args {
    list_samples: bool,
    format: Option<String>,
    input: String,
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
        LongOpt::new("list-samples", HasArg::None, b'l' as i32),
    ];

    let mut list_samples = false;
    let mut format = None;

    let mut g = Getopt::new("f:l", &long_opts, argv);
    loop {
        match g.next() {
            Ok(Some(m)) => match m.code {
                v if v == b'l' as i32 => list_samples = true,
                v if v == b'f' as i32 => format = m.value,
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
        input,
    })
}

fn run<W: Write>(args: &Args, mut out: W) -> io::Result<()> {
    let input = materialize_input(&args.input)?;
    if args.list_samples {
        for sample in sample_names_from_path(&input)? {
            writeln!(out, "{sample}")?;
        }
    }
    if let Some(format) = &args.format {
        query_format_from_path(&input, format, &mut out)?;
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

fn sample_names_from_path<P>(path: P) -> io::Result<Vec<String>>
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

fn query_format_from_path<W: Write>(path: &Path, format: &str, out: &mut W) -> io::Result<()> {
    let text = vcf_text_from_path(path)?;
    query_format_text(&text, format, out)
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

fn query_format_text<W: Write>(text: &str, format: &str, out: &mut W) -> io::Result<()> {
    let mut samples = Vec::new();
    for line in text.lines() {
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with("#CHROM\t") {
            samples = line.split('\t').skip(9).map(ToOwned::to_owned).collect();
            continue;
        }
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let record = TextRecord::parse(line, &samples);
        let rendered = render_format(format, &record);
        out.write_all(rendered.as_bytes())?;
    }
    Ok(())
}

#[derive(Debug)]
struct TextRecord<'a> {
    fields: Vec<&'a str>,
    samples: &'a [String],
}

impl<'a> TextRecord<'a> {
    fn parse(line: &'a str, samples: &'a [String]) -> Self {
        Self {
            fields: line.split('\t').collect(),
            samples,
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
            "LINE" => self.fields.join("\t"),
            _ => ".".to_string(),
        }
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

    fn format_value(&self, sample_index: usize, key: &str) -> String {
        let Some(format) = self.fields.get(8) else {
            return ".".into();
        };
        let Some(sample) = self.fields.get(9 + sample_index) else {
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
        for sample_index in 0..record.samples.len() {
            render_segment(block, record, Some(sample_index), &mut out);
        }
        rest = &after_start[end + 1..];
    }
    render_segment(rest, record, None, &mut out);
    out
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
    let token = token.strip_prefix('/').unwrap_or(token);
    if token == "SAMPLE" {
        return sample_index
            .and_then(|i| record.samples.get(i))
            .cloned()
            .unwrap_or_else(|| ".".into());
    }
    if let Some(key) = token.strip_prefix("INFO/") {
        return record.info(key);
    }
    if let Some(key) = token
        .strip_prefix("FMT/")
        .or_else(|| token.strip_prefix("FORMAT/"))
    {
        return sample_index
            .map(|i| record.format_value(i, key))
            .unwrap_or_else(|| ".".into());
    }
    if let Some(i) = sample_index {
        let value = record.format_value(i, token);
        if value != "." {
            return value;
        }
    }
    match token {
        "CHROM" | "POS" | "ID" | "REF" | "ALT" | "QUAL" | "FILTER" | "INFO" | "LINE" => {
            record.core(token)
        }
        _ => record.info(token),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        query_format_text(text, "%CHROM\\t%POS\\t%DP[\\t%SAMPLE=%GT:%DP]\\n", &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "1\t2\t7\tA=0/1:3\tB=0/0:4\n"
        );
    }
}
