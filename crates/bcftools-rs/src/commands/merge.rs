//! Focused `bcftools merge` implementation (upstream `vcfmerge.c`).
//!
//! This first local slice merges records that are present in the same order in
//! every input and have identical site fields. Full synced-reader merging,
//! allele unification, INFO rules, gVCF mode, and missing-sample synthesis
//! remain tracked in `TODO.md`.

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
use crate::vcf_compat::normalize_vcf_text;

const USAGE: &str = "\n\
About:   Merge VCF/BCF files from non-overlapping sample sets.\n\
Usage:   bcftools merge [OPTIONS] <A.vcf.gz> <B.vcf.gz> [...]\n\
\n\
Options:\n\
    -l, --file-list FILE            Read input file names from FILE\n\
    -m, --merge TYPE                Accepted for command-shape compatibility in this same-site slice\n\
    -o, --output FILE               Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9]  u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
        --force-samples             Allow duplicate sample names by prefixing later inputs with the input index\n\
        --no-version                Accepted for command-shape compatibility\n\
\n";

#[derive(Debug)]
struct Args {
    inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    force_samples: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

#[derive(Debug)]
struct VcfInput {
    meta: Vec<String>,
    fixed_header: Vec<String>,
    samples: Vec<String>,
    records: Vec<RecordLine>,
}

#[derive(Debug)]
struct RecordLine {
    fixed: Vec<String>,
    samples: Vec<String>,
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
                eprintln!("{}", fmt_etag("main_vcfmerge", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfmerge", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut inputs = Vec::new();
    let mut file_list = None;
    let mut output = None;
    let mut output_kind = OutputKind::VcfText;
    let mut force_samples = false;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-l" | "--file-list" => {
                file_list = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-m" | "--merge" => {
                let _ = next_string(&mut iter, raw.as_ref())?;
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--force-samples" => force_samples = true,
            "--no-version" => {}
            _ if raw.starts_with("--file-list=") => {
                file_list = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--merge=") => {}
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-l") && raw.len() > 2 => {
                file_list = Some(PathBuf::from(&raw[2..]))
            }
            _ if raw.starts_with("-m") && raw.len() > 2 => {}
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!(
                    "unrecognized option '{raw}' in this local merge slice"
                )));
            }
            _ => inputs.push(PathBuf::from(raw.as_ref())),
        }
    }

    if let Some(path) = file_list {
        inputs.extend(read_file_list(&path)?);
    }
    if inputs.len() < 2 {
        return Err(ParseOutcome::Error(
            "expected at least two input VCF/BCF paths".into(),
        ));
    }

    Ok(Args {
        inputs,
        output,
        output_kind,
        force_samples,
    })
}

fn read_file_list(path: &Path) -> Result<Vec<PathBuf>, ParseOutcome> {
    let text = fs::read_to_string(path)
        .map_err(|e| ParseOutcome::Error(format!("failed to read {}: {e}", path.display())))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(PathBuf::from)
        .collect())
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
    let mut inputs = Vec::new();
    for path in &args.inputs {
        inputs.push(parse_vcf_text(&read_vcf_text(path)?)?);
    }
    let merged = merge_inputs(&inputs, args.force_samples)?;
    write_output(merged.as_bytes(), args)
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
        ".bcftools-rs-merge-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn parse_vcf_text(text: &str) -> io::Result<VcfInput> {
    let mut meta = Vec::new();
    let mut fixed_header = Vec::new();
    let mut samples = Vec::new();
    let mut records = Vec::new();

    for line in text.lines() {
        if line.starts_with("##") {
            meta.push(line.to_owned());
        } else if line.starts_with("#CHROM") {
            let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
            if fields.len() < 8 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "invalid #CHROM header",
                ));
            }
            fixed_header = fields[..fields.len().min(9)].to_vec();
            if fields.len() > 9 {
                samples = fields[9..].to_vec();
            }
        } else if !line.trim().is_empty() {
            let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
            if fields.len() < fixed_header.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("record has fewer columns than header: {line}"),
                ));
            }
            let fixed_len = fixed_header.len();
            records.push(RecordLine {
                fixed: fields[..fixed_len].to_vec(),
                samples: fields[fixed_len..].to_vec(),
            });
        }
    }

    if fixed_header.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "missing #CHROM header",
        ));
    }
    Ok(VcfInput {
        meta,
        fixed_header,
        samples,
        records,
    })
}

fn merge_inputs(inputs: &[VcfInput], force_samples: bool) -> io::Result<String> {
    let first = inputs
        .first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "no inputs"))?;
    let mut sample_names = Vec::new();
    let mut seen_samples = HashSet::new();
    for (input_idx, input) in inputs.iter().enumerate() {
        for sample in &input.samples {
            let mut name = sample.clone();
            if !seen_samples.insert(name.clone()) {
                if !force_samples {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("duplicate sample name '{sample}'"),
                    ));
                }
                name = format!("{}:{sample}", input_idx + 1);
                seen_samples.insert(name.clone());
            }
            sample_names.push(name);
        }
    }

    let mut out = String::new();
    for line in &first.meta {
        out.push_str(line);
        out.push('\n');
    }
    out.push_str(&first.fixed_header.join("\t"));
    if !sample_names.is_empty() {
        out.push('\t');
        out.push_str(&sample_names.join("\t"));
    }
    out.push('\n');

    let nrecords = first.records.len();
    for input in &inputs[1..] {
        if input.records.len() != nrecords {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "all inputs must contain the same number of records in this local merge slice",
            ));
        }
        if input.fixed_header[..input.fixed_header.len().min(9)] != first.fixed_header {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "inputs must have compatible fixed VCF columns",
            ));
        }
    }

    for idx in 0..nrecords {
        let base = &first.records[idx];
        let mut samples = base.samples.clone();
        for input in &inputs[1..] {
            let record = &input.records[idx];
            if record.fixed != base.fixed {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "record mismatch at input record {}: expected {}:{} {}>{}, found {}:{} {}>{}",
                        idx + 1,
                        base.fixed[0],
                        base.fixed[1],
                        base.fixed[3],
                        base.fixed[4],
                        record.fixed[0],
                        record.fixed[1],
                        record.fixed[3],
                        record.fixed[4],
                    ),
                ));
            }
            samples.extend(record.samples.iter().cloned());
        }
        out.push_str(&base.fixed.join("\t"));
        if !samples.is_empty() {
            out.push('\t');
            out.push_str(&samples.join("\t"));
        }
        out.push('\n');
    }

    Ok(out)
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
    fn merges_same_site_sample_columns() {
        let a = parse_vcf_text(
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        let b = parse_vcf_text(
            "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\n\
1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();
        let merged = merge_inputs(&[a, b], false).unwrap();
        assert!(merged.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"));
        assert!(merged.contains("1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1"));
    }
}
