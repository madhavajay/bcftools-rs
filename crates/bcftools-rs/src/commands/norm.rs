//! Focused `bcftools norm` implementation (upstream `vcfnorm.c`).
//!
//! This first local slice supports duplicate-record removal with
//! `-d/--rm-dup`. Left alignment, reference checks, atomization, split/join
//! multiallelics, and INFO/FORMAT projection remain tracked in `TODO.md`.

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
About:   Left-align and normalize indels; this local slice supports duplicate removal.\n\
Usage:   bcftools norm [OPTIONS] <in.vcf.gz>\n\
\n\
Options:\n\
    -d, --rm-dup TYPE              Remove duplicate records: snps|indels|both|all|exact|none|any\n\
    -f, --fasta-ref FILE           Accepted by the duplicate-removal slice for command compatibility\n\
    -o, --output FILE              Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9] u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
        --no-version               Accepted for command-shape compatibility\n\
\n";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    rm_dup: DupMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DupMode {
    None,
    Exact,
    Snps,
    Indels,
    Both,
    All,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum VariantKind {
    Snp,
    Indel,
    Other,
}

#[derive(Debug, Hash, PartialEq, Eq)]
struct DupKey {
    chrom: String,
    pos: String,
    rest: String,
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
                eprintln!("{}", fmt_etag("main_vcfnorm", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfnorm", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut input = None;
    let mut output = None;
    let mut output_kind = OutputKind::VcfText;
    let mut rm_dup = DupMode::None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-d" | "--rm-dup" | "--rm-dups" => {
                rm_dup = parse_dup_mode(&next_string(&mut iter, raw.as_ref())?)?
            }
            "-f" | "--fasta-ref" => {
                let _ = next_string(&mut iter, raw.as_ref())?;
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--no-version" => {}
            _ if raw.starts_with("--rm-dup=") || raw.starts_with("--rm-dups=") => {
                rm_dup = parse_dup_mode(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--fasta-ref=") => {
                let _ = value_after_equals(&raw);
            }
            _ if raw.starts_with("-f") && raw.len() > 2 => {
                let _ = &raw[2..];
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-d") && raw.len() > 2 => rm_dup = parse_dup_mode(&raw[2..])?,
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!(
                    "unrecognized option '{raw}' in this local norm slice"
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
    Ok(Args {
        input,
        output,
        output_kind,
        rm_dup,
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

fn parse_dup_mode(raw: &str) -> Result<DupMode, ParseOutcome> {
    match raw {
        "none" | "exact" => Ok(DupMode::Exact),
        "snps" => Ok(DupMode::Snps),
        "indels" => Ok(DupMode::Indels),
        "both" | "any" => Ok(DupMode::Both),
        "all" => Ok(DupMode::All),
        _ => Err(ParseOutcome::Error(format!(
            "unknown duplicate mode '{raw}'"
        ))),
    }
}

fn run(args: &Args) -> io::Result<()> {
    let input = read_vcf_text(&args.input)?;
    let output = remove_duplicates(&input, args.rm_dup)?;
    write_output(output.as_bytes(), args)
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
        ".bcftools-rs-norm-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn remove_duplicates(input: &str, mode: DupMode) -> io::Result<String> {
    let mut out = String::with_capacity(input.len());
    let mut seen = HashSet::new();
    let filter_header_seen = input
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    let mut filter_header_inserted = false;

    for line in input.lines() {
        if line.starts_with("##fileformat=") {
            out.push_str(line);
            out.push('\n');
            if !filter_header_seen && !filter_header_inserted {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
                filter_header_inserted = true;
            }
            continue;
        }
        if line.starts_with("#CHROM") {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid VCF record with fewer than 8 columns: {line}"),
            ));
        }
        let keys = duplicate_keys(&fields, mode);
        let is_dup = keys.iter().any(|key| seen.contains(key));
        if !is_dup {
            seen.extend(keys);
            out.push_str(line);
            out.push('\n');
        }
    }

    Ok(out)
}

fn duplicate_keys(fields: &[&str], mode: DupMode) -> Vec<DupKey> {
    let mut keys = Vec::new();
    match mode {
        DupMode::All => keys.push(dup_key(fields, "")),
        DupMode::Exact | DupMode::None => keys.push(dup_key(fields, &exact_rest(fields))),
        DupMode::Snps => {
            if allele_class(fields).has_snp {
                keys.push(dup_key(fields, "Snp"));
            }
        }
        DupMode::Indels => {
            if allele_class(fields).has_indel {
                keys.push(dup_key(fields, "Indel"));
            }
        }
        DupMode::Both => {
            let class = allele_class(fields);
            if class.has_snp {
                keys.push(dup_key(fields, "Snp"));
            }
            if class.has_indel {
                keys.push(dup_key(fields, "Indel"));
            }
        }
    }
    keys
}

fn dup_key(fields: &[&str], rest: &str) -> DupKey {
    DupKey {
        chrom: fields[0].to_owned(),
        pos: fields[1].to_owned(),
        rest: rest.to_owned(),
    }
}

fn exact_rest(fields: &[&str]) -> String {
    let mut alts: Vec<&str> = fields[4].split(',').collect();
    alts.sort_unstable();
    format!("{}:{}", fields[3], alts.join(","))
}

#[derive(Clone, Copy)]
struct AlleleClass {
    has_snp: bool,
    has_indel: bool,
}

fn allele_class(fields: &[&str]) -> AlleleClass {
    let alts: Vec<&str> = fields[4].split(',').collect();
    let kind = classify_kind(fields[3], &alts);
    AlleleClass {
        has_snp: matches!(kind, VariantKind::Snp | VariantKind::Other)
            && alts
                .iter()
                .any(|alt| alt.len() == 1 && fields[3].len() == 1 && *alt != "." && *alt != "*"),
        has_indel: matches!(kind, VariantKind::Indel | VariantKind::Other)
            && alts
                .iter()
                .any(|alt| alt.len() != fields[3].len() && *alt != "." && *alt != "*"),
    }
}

fn classify_kind(reference: &str, alts: &[&str]) -> VariantKind {
    let mut has_snp = false;
    let mut has_indel = false;
    for alt in alts {
        if *alt == "." || *alt == "*" {
            continue;
        }
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
    fn classify_kind_distinguishes_snp_and_indel() {
        assert_eq!(classify_kind("A", &["C"]), VariantKind::Snp);
        assert_eq!(classify_kind("A", &["AT"]), VariantKind::Indel);
        assert_eq!(classify_kind("AT", &["A"]), VariantKind::Indel);
        assert_eq!(classify_kind("A", &["C", "AT"]), VariantKind::Other);
    }

    #[test]
    fn rmdup_all_keeps_first_record_per_position() {
        let input = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t2\t.\tA\tC\t.\t.\t.\n\
1\t2\t.\tA\tG\t.\t.\t.\n";
        let out = remove_duplicates(input, DupMode::All).unwrap();
        assert!(out.contains("1\t2\t.\tA\tC\t.\t.\t."));
        assert!(!out.contains("1\t2\t.\tA\tG\t.\t.\t."));
    }
}
