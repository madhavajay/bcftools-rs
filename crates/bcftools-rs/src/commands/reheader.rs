//! Partial port of `bcftools reheader` (upstream `reheader.c`).
//!
//! This implements sample renaming, header replacement, and FAI-driven contig
//! updates for VCF text/BGZF input. BCF in-place reheadering and threaded BCF
//! output remain explicitly unsupported.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read as _, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::str::FromStr as _;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::{bcf, vcf};

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};
use crate::io::apply_verbosity;

const USAGE: &str = "\n\
About:   Modify header of VCF/BCF files, change sample names.\n\
Usage:   bcftools reheader [OPTIONS] <in.vcf>\n\
\n\
Options:\n\
    -f, --fai FILE             Update sequences and their lengths from the .fai file\n\
    -h, --header FILE          New header\n\
    -o, --output FILE          Write output to a file [standard output]\n\
    -n, --samples-list LIST    New sample names given as a comma-separated list\n\
    -N, --samples-file FILE    New sample names in a file, see the man page for details\n\
    -s, --samples FILE         Alias for --samples-file\n\
        --in-place             Reheader BCF in place (NOT IMPLEMENTED)\n\
    -T, --temp-prefix PATH     Ignored; was template for temporary file name\n\
        --threads INT          Use multithreading with INT worker threads (BCF only) [0]\n\
    -v, --verbosity INT        Verbosity level\n\
\n";

#[derive(Debug, Clone, PartialEq, Eq)]
struct Args {
    samples: Option<SampleSource>,
    header: Option<PathBuf>,
    fai: Option<PathBuf>,
    output: Option<PathBuf>,
    in_place: bool,
    input: PathBuf,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum SampleSource {
    File(PathBuf),
    List(String),
}

/// Subcommand entry point. `argv[0]` is `"reheader"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_reheader", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_reheader", &message));
            ExitCode::FAILURE
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ParseOutcome {
    Usage,
    Error(String),
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let long_opts = [
        LongOpt::new("temp-prefix", HasArg::Required, b'T' as i32),
        LongOpt::new("fai", HasArg::Required, b'f' as i32),
        LongOpt::new("header", HasArg::Required, b'h' as i32),
        LongOpt::new("output", HasArg::Required, b'o' as i32),
        LongOpt::new("samples", HasArg::Required, b's' as i32),
        LongOpt::new("samples-file", HasArg::Required, b'N' as i32),
        LongOpt::new("samples-list", HasArg::Required, b'n' as i32),
        LongOpt::new("in-place", HasArg::None, 2),
        LongOpt::new("threads", HasArg::Required, 1),
        LongOpt::new("verbosity", HasArg::Required, b'v' as i32),
    ];

    let mut samples: Option<SampleSource> = None;
    let mut header = None;
    let mut fai = None;
    let mut output = None;
    let mut in_place = false;

    let mut g = Getopt::new("s:h:o:f:T:v:N:n:", &long_opts, argv);
    loop {
        match g.next() {
            Ok(Some(m)) => match m.code {
                v if v == b'v' as i32 => {
                    if apply_verbosity(m.value.as_deref().unwrap_or("0")).is_err() {
                        return Err(ParseOutcome::Error(format!(
                            "Could not parse argument: --verbosity {}",
                            m.value.as_deref().unwrap_or("")
                        )));
                    }
                }
                v if v == b'T' as i32 => {}
                v if v == b'o' as i32 => {
                    output = Some(PathBuf::from(m.value.as_deref().unwrap_or_default()));
                }
                v if v == b's' as i32 || v == b'N' as i32 => {
                    samples = Some(SampleSource::File(PathBuf::from(
                        m.value.as_deref().unwrap_or_default(),
                    )));
                }
                v if v == b'n' as i32 => {
                    samples = Some(SampleSource::List(
                        m.value.as_deref().unwrap_or_default().to_string(),
                    ));
                }
                v if v == b'h' as i32 => {
                    header = Some(PathBuf::from(m.value.as_deref().unwrap_or_default()));
                }
                v if v == b'f' as i32 => {
                    fai = Some(PathBuf::from(m.value.as_deref().unwrap_or_default()));
                }
                1 => {
                    return Err(ParseOutcome::Error(
                        "reheader --threads is not yet implemented for BCF output".into(),
                    ));
                }
                2 => in_place = true,
                _ => return Err(ParseOutcome::Usage),
            },
            Ok(None) => break,
            Err(_) => return Err(ParseOutcome::Usage),
        }
    }

    if samples.is_none() && header.is_none() && fai.is_none() {
        return Err(ParseOutcome::Usage);
    }
    let rest = g.rest();
    let input = match rest.len() {
        0 => PathBuf::from("-"),
        1 => PathBuf::from(&rest[0]),
        _ => return Err(ParseOutcome::Usage),
    };
    if in_place && output.is_some() {
        return Err(ParseOutcome::Error(
            "reheader --in-place cannot be combined with --output".into(),
        ));
    }

    Ok(Args {
        samples,
        header,
        fai,
        output,
        in_place,
        input,
    })
}

fn run(args: &Args) -> io::Result<()> {
    if args.input == Path::new("-") {
        if args.in_place {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "reheader --in-place cannot read from standard input",
            ));
        }
        let tmp = stdin_tmp_path();
        let mut data = Vec::new();
        io::stdin().lock().read_to_end(&mut data)?;
        fs::write(&tmp, data)?;
        let mut tmp_args = args.clone();
        tmp_args.input = tmp.clone();
        let result = run(&tmp_args);
        let _ = fs::remove_file(&tmp);
        return result;
    }

    let sample_specs = args.samples.as_ref().map(read_sample_specs).transpose()?;
    let replacement_header = args.header.as_ref().map(read_header_file).transpose()?;
    let fai = args.fai.as_ref().map(read_fai).transpose()?;
    let fmt = format::detect_path(&args.input).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return reheader_bcf(
            args,
            sample_specs.as_deref(),
            replacement_header.as_deref(),
            fai.as_deref(),
        );
    }
    if fmt.compression == Compression::Gzip {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "reheader gzip-compressed VCF input is not supported; use BGZF-compressed VCF",
        ));
    }
    if args.in_place {
        return Err(io::Error::new(
            io::ErrorKind::Unsupported,
            "reheader --in-place is only supported for BCF input",
        ));
    }

    if fmt.compression == Compression::Bgzf {
        let input = File::open(&args.input)?;
        let decoder = flate2::read::MultiGzDecoder::new(input);
        match &args.output {
            Some(path) => {
                let output = File::create(path)?;
                let bgzf = htslib_rs::bgzf::io::Writer::new(output);
                reheader_vcf(
                    BufReader::new(decoder),
                    bgzf,
                    sample_specs.as_deref(),
                    replacement_header.as_deref(),
                    fai.as_deref(),
                )
            }
            None => {
                let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
                reheader_vcf(
                    BufReader::new(decoder),
                    bgzf,
                    sample_specs.as_deref(),
                    replacement_header.as_deref(),
                    fai.as_deref(),
                )
            }
        }
    } else {
        let input = File::open(&args.input)?;
        match &args.output {
            Some(path) => {
                let output = File::create(path)?;
                reheader_vcf(
                    BufReader::new(input),
                    output,
                    sample_specs.as_deref(),
                    replacement_header.as_deref(),
                    fai.as_deref(),
                )
            }
            None => reheader_vcf(
                BufReader::new(input),
                io::stdout().lock(),
                sample_specs.as_deref(),
                replacement_header.as_deref(),
                fai.as_deref(),
            ),
        }
    }
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-reheader-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn reheader_bcf(
    args: &Args,
    sample_specs: Option<&[String]>,
    replacement_header: Option<&str>,
    fai: Option<&[FaiContig]>,
) -> io::Result<()> {
    let mut reader = File::open(&args.input).map(bcf::io::Reader::new)?;
    let source_header = reader.read_header()?;
    let mut source_header_text = Vec::new();
    vcf::io::Writer::new(&mut source_header_text).write_header(&source_header)?;
    let source_header_text = String::from_utf8(source_header_text)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    let header_text = replacement_header.unwrap_or(&source_header_text);
    let rewritten_header_text = rewrite_header_text(header_text, sample_specs, fai)?;
    let rewritten_header = vcf::Header::from_str(&rewritten_header_text)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    if rewritten_header.sample_names().len() != source_header.sample_names().len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "BCF reheader requires the same number of samples",
        ));
    }

    if args.in_place {
        let tmp = in_place_temp_path(&args.input);
        let result = (|| {
            let output = File::create(&tmp)?;
            write_bcf_records(reader, &source_header, &rewritten_header, output)?;
            fs::rename(&tmp, &args.input)
        })();
        if result.is_err() {
            let _ = fs::remove_file(&tmp);
        }
        return result;
    }

    match &args.output {
        Some(path) => {
            let output = File::create(path)?;
            write_bcf_records(reader, &source_header, &rewritten_header, output)
        }
        None => write_bcf_records(
            reader,
            &source_header,
            &rewritten_header,
            io::stdout().lock(),
        ),
    }
}

fn in_place_temp_path(input: &Path) -> PathBuf {
    let parent = input.parent().unwrap_or_else(|| std::path::Path::new("."));
    let name = input
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("input.bcf");
    parent.join(format!(
        ".{name}.bcftools-rs-reheader.{}.tmp",
        std::process::id()
    ))
}

fn write_bcf_records<R, W>(
    mut reader: bcf::io::Reader<R>,
    source_header: &vcf::Header,
    rewritten_header: &vcf::Header,
    output: W,
) -> io::Result<()>
where
    R: io::Read,
    W: Write,
{
    use htslib_rs::vcf::variant::io::Write as _;

    let mut writer = bcf::io::Writer::new(output);
    writer.write_variant_header(rewritten_header)?;
    for result in reader.record_bufs(source_header) {
        let record = result?;
        writer.write_variant_record(rewritten_header, &record)?;
    }
    writer.try_finish()
}

fn read_header_file(path: &PathBuf) -> io::Result<String> {
    let mut data = String::new();
    File::open(path)?.read_to_string(&mut data)?;
    while data.ends_with(char::is_whitespace) {
        data.pop();
    }
    data.push('\n');
    Ok(data)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FaiContig {
    name: String,
    len: String,
}

fn read_fai(path: &PathBuf) -> io::Result<Vec<FaiContig>> {
    let file = File::open(path)?;
    let mut contigs = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let mut fields = line.split('\t');
        let Some(name) = fields.next() else {
            continue;
        };
        let Some(len) = fields.next() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Could not parse {}", path.display()),
            ));
        };
        if !name.is_empty() {
            contigs.push(FaiContig {
                name: name.to_string(),
                len: len.to_string(),
            });
        }
    }
    if contigs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Could not parse {}", path.display()),
        ));
    }
    Ok(contigs)
}

fn read_sample_specs(source: &SampleSource) -> io::Result<Vec<String>> {
    let specs: Vec<String> = match source {
        SampleSource::File(path) => {
            let mut data = String::new();
            File::open(path)?.read_to_string(&mut data)?;
            data.lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(ToOwned::to_owned)
                .collect()
        }
        SampleSource::List(list) => list
            .split(',')
            .map(str::trim)
            .filter(|item| !item.is_empty())
            .map(ToOwned::to_owned)
            .collect(),
    };
    if specs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no sample names were provided",
        ));
    }
    Ok(specs)
}

fn reheader_vcf<R, W>(
    mut reader: R,
    mut writer: W,
    sample_specs: Option<&[String]>,
    replacement_header: Option<&str>,
    fai: Option<&[FaiContig]>,
) -> io::Result<()>
where
    R: BufRead,
    W: Write,
{
    let mut line = String::new();
    let mut header = String::new();
    let mut pending_record = None;

    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }

        if !line.starts_with('#') {
            pending_record = Some(line.clone());
            break;
        }

        header.push_str(&line);
        if line.starts_with("#CHROM\t") {
            break;
        }
    }

    let header = replacement_header.unwrap_or(&header);
    if !header.lines().any(|line| line.starts_with("#CHROM\t")) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Could not parse the header: missing #CHROM line",
        ));
    }
    writer.write_all(rewrite_header_text(header, sample_specs, fai)?.as_bytes())?;

    if let Some(record) = pending_record {
        writer.write_all(record.as_bytes())?;
    }
    io::copy(&mut reader, &mut writer)?;
    Ok(())
}

fn rewrite_header_text(
    header: &str,
    sample_specs: Option<&[String]>,
    fai: Option<&[FaiContig]>,
) -> io::Result<String> {
    let header = ensure_pass_filter(header);
    let header = if let Some(fai) = fai {
        update_header_contigs(&header, fai)?
    } else {
        header
    };

    if let Some(sample_specs) = sample_specs {
        let mut out = String::new();
        let mut saw_chrom = false;
        for line in header.lines() {
            if line.starts_with("#CHROM\t") {
                saw_chrom = true;
                out.push_str(&rename_chrom_line(line, sample_specs)?);
            } else {
                out.push_str(line);
            }
            out.push('\n');
        }
        if !saw_chrom {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "Could not parse the header: missing #CHROM line",
            ));
        }
        Ok(out)
    } else {
        Ok(header)
    }
}

fn ensure_pass_filter(header: &str) -> String {
    if header
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"))
    {
        return header.to_string();
    }

    let mut out = String::new();
    let mut inserted = false;
    for line in header.lines() {
        out.push_str(line);
        out.push('\n');
        if !inserted && line.starts_with("##fileformat=") {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            inserted = true;
        }
    }
    if !inserted {
        format!("##FILTER=<ID=PASS,Description=\"All filters passed\">\n{header}")
    } else {
        out
    }
}

fn update_header_contigs(header: &str, fai: &[FaiContig]) -> io::Result<String> {
    let mut out = String::new();
    let mut seen = Vec::new();
    let mut inserted_missing = false;

    for line in header.lines() {
        if line.starts_with("##contig=<") {
            if let Some(id) = extract_contig_id(line)
                && let Some(contig) = fai.iter().find(|contig| contig.name == id)
            {
                seen.push(contig.name.clone());
                out.push_str(&update_contig_line(line, contig));
                out.push('\n');
            }
            continue;
        }

        if line.starts_with("#CHROM\t") && !inserted_missing {
            for contig in fai {
                if !seen.iter().any(|name| name == &contig.name) {
                    out.push_str(&format_contig_line(contig));
                    out.push('\n');
                }
            }
            inserted_missing = true;
        }

        out.push_str(line);
        out.push('\n');
    }

    if !inserted_missing {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Failed to parse the header, #CHROM not found",
        ));
    }

    Ok(out)
}

fn extract_contig_id(line: &str) -> Option<String> {
    let inner = line.strip_prefix("##contig=<")?.strip_suffix('>')?;
    for field in split_structured_fields(inner) {
        let (key, value) = field.split_once('=')?;
        if key.trim() == "ID" {
            return Some(value.trim().trim_matches('"').to_string());
        }
    }
    None
}

fn format_contig_line(contig: &FaiContig) -> String {
    format!("##contig=<ID={},length={}>", contig.name, contig.len)
}

fn update_contig_line(line: &str, contig: &FaiContig) -> String {
    let Some(inner) = line
        .strip_prefix("##contig=<")
        .and_then(|line| line.strip_suffix('>'))
    else {
        return format_contig_line(contig);
    };

    let mut fields = Vec::new();
    for field in split_structured_fields(inner) {
        let Some((key, value)) = field.split_once('=') else {
            continue;
        };
        let key = key.trim();
        if key == "ID" || key == "length" {
            continue;
        }
        fields.push(format!("{key}={}", value.trim()));
    }

    let mut out = String::with_capacity(line.len() + contig.len.len() + 8);
    out.push_str("##contig=<ID=");
    out.push_str(&contig.name);
    out.push_str(",length=");
    out.push_str(&contig.len);
    for field in fields {
        out.push(',');
        out.push_str(&field);
    }
    out.push('>');
    out
}

fn split_structured_fields(s: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut start = 0;
    let mut quote = false;
    let mut escaped = false;
    let mut angle_depth = 0usize;

    for (idx, ch) in s.char_indices() {
        if quote {
            if escaped {
                escaped = false;
            } else if ch == '\\' {
                escaped = true;
            } else if ch == '"' {
                quote = false;
            }
            continue;
        }

        match ch {
            '"' => quote = true,
            '<' => angle_depth += 1,
            '>' => angle_depth = angle_depth.saturating_sub(1),
            ',' if angle_depth == 0 => {
                fields.push(s[start..idx].trim());
                start = idx + ch.len_utf8();
            }
            _ => {}
        }
    }

    fields.push(s[start..].trim());
    fields
}

fn rename_chrom_line(line: &str, sample_specs: &[String]) -> io::Result<String> {
    let mut columns: Vec<String> = line.split('\t').map(ToOwned::to_owned).collect();
    if columns.len() < 8
        || columns[..8]
            != [
                "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            ]
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Could not parse the header",
        ));
    }
    if columns.len() == 8 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Error: missing FORMAT fields, cowardly refusing to add samples",
        ));
    }
    if columns.get(8).map(String::as_str) != Some("FORMAT") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "Could not parse the header",
        ));
    }

    if let Some(pairs) = parse_sample_pairs(sample_specs) {
        for sample in &mut columns[9..] {
            if let Some((_, new)) = pairs.iter().find(|(old, _)| old == sample) {
                *sample = new.clone();
            }
        }
    } else {
        let current = columns.len().saturating_sub(9);
        if current != sample_specs.len() {
            eprintln!(
                "Warning: different number of samples: {} vs {}",
                sample_specs.len(),
                current
            );
        }
        columns.truncate(9);
        columns.extend(sample_specs.iter().cloned());
    }

    Ok(columns.join("\t"))
}

fn parse_sample_pairs(sample_specs: &[String]) -> Option<Vec<(String, String)>> {
    let mut pairs = Vec::with_capacity(sample_specs.len());
    for spec in sample_specs {
        let fields = split_escaped_whitespace(spec);
        if fields.len() != 2 {
            return None;
        }
        pairs.push((fields[0].clone(), fields[1].clone()));
    }
    Some(pairs)
}

fn split_escaped_whitespace(s: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut escaped = false;
    for ch in s.chars() {
        if escaped {
            current.push(ch);
            escaped = false;
            continue;
        }
        if ch == '\\' {
            escaped = true;
            continue;
        }
        if ch.is_whitespace() {
            if !current.is_empty() {
                fields.push(std::mem::take(&mut current));
            }
            continue;
        }
        current.push(ch);
    }
    if escaped {
        current.push('\\');
    }
    if !current.is_empty() {
        fields.push(current);
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n";

    #[test]
    fn replaces_sample_columns() {
        let mut out = Vec::new();
        reheader_vcf(
            VCF.as_bytes(),
            &mut out,
            Some(&["X".into(), "Y".into()]),
            None,
            None,
        )
        .unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
        assert!(out.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
    }

    #[test]
    fn renames_sample_pairs() {
        let mut out = Vec::new();
        reheader_vcf(VCF.as_bytes(), &mut out, Some(&["B Z".into()]), None, None).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tZ\n"));
    }

    #[test]
    fn replaces_header_text() {
        let header = "##fileformat=VCFv4.3\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n";
        let mut out = Vec::new();
        reheader_vcf(VCF.as_bytes(), &mut out, None, Some(header), None).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.starts_with(
            "##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n#CHROM\t"
        ));
        assert!(out.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
    }

    #[test]
    fn updates_contigs_from_fai() {
        let input = "##fileformat=VCFv4.2\n\
##contig=<ID=old,length=1>\n\
##contig=<assembly=B36,ID=1,species=\"Homo sapiens\",length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n";
        let fai = [
            FaiContig {
                name: "1".into(),
                len: "100".into(),
            },
            FaiContig {
                name: "2".into(),
                len: "200".into(),
            },
        ];
        let mut out = Vec::new();
        reheader_vcf(input.as_bytes(), &mut out, None, None, Some(&fai)).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(!out.contains("ID=old"));
        assert!(out.contains("##contig=<ID=1,length=100,assembly=B36,species=\"Homo sapiens\">\n"));
        assert!(out.contains("##contig=<ID=2,length=200>\n#CHROM\t"));
    }

    #[test]
    fn parses_escaped_pair_fields() {
        assert_eq!(
            split_escaped_whitespace(r"old\ sample new\ sample"),
            ["old sample".to_string(), "new sample".to_string()]
        );
    }
}
