//! Focused `bcftools consensus` implementation (upstream `consensus.c`).
//!
//! This first local slice applies simple VCF REF -> ALT replacements to FASTA
//! records and preserves FASTA headers/order. Advanced upstream modes such as
//! masks, chains, sample-aware haplotype selection, marks, absent filling, and
//! overlap policy remain tracked in `TODO.md`.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::diagnostics::fmt_etag;
use crate::vcf_compat::NormalizeFileformat;

const USAGE: &str = "\n\
About:   Create consensus sequence by applying VCF variants to a FASTA reference.\n\
Usage:   bcftools consensus [OPTIONS] <in.vcf.gz>\n\
\n\
Options:\n\
   -f, --fasta-ref FILE          Reference sequence in FASTA format\n\
   -H, --haplotype WHICH         Choose ALT allele by 1-based ALT index for this local slice [1]\n\
   -s, --sample NAME             Accepted for command-shape compatibility; sample-aware selection is deferred\n\
   -o, --output FILE             Write output to a file [standard output]\n\
\n";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    fasta: PathBuf,
    output: Option<PathBuf>,
    haplotype: usize,
}

#[derive(Debug)]
enum ParseOutcome {
    Usage,
    Error(String),
}

#[derive(Debug, Clone)]
struct FastaRecord {
    header: String,
    chrom: String,
    start: i64,
    seq: Vec<u8>,
}

#[derive(Debug)]
struct Variant {
    chrom: String,
    pos: i64,
    reference: String,
    alts: Vec<String>,
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_consensus", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_consensus", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut input = None;
    let mut fasta = None;
    let mut output = None;
    let mut haplotype = 1usize;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-f" | "--fasta-ref" => {
                fasta = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-H" | "--haplotype" => {
                haplotype = parse_haplotype(&next_string(&mut iter, raw.as_ref())?)?;
            }
            "-s" | "--sample" => {
                let _ = next_string(&mut iter, raw.as_ref())?;
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            _ if raw.starts_with("--fasta-ref=") => {
                fasta = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--haplotype=") => {
                haplotype = parse_haplotype(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--sample=") => {}
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("-f") && raw.len() > 2 => fasta = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-H") && raw.len() > 2 => {
                haplotype = parse_haplotype(&raw[2..])?;
            }
            _ if raw.starts_with("-s") && raw.len() > 2 => {}
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!("unrecognized option '{raw}'")));
            }
            _ if input.is_none() => input = Some(PathBuf::from(raw.as_ref())),
            _ => {
                return Err(ParseOutcome::Error(format!(
                    "unexpected extra input '{}'",
                    raw
                )));
            }
        }
    }

    let input =
        input.ok_or_else(|| ParseOutcome::Error("expected one input VCF/BCF path".into()))?;
    let fasta = fasta.ok_or_else(|| ParseOutcome::Error("expected -f/--fasta-ref".into()))?;

    Ok(Args {
        input,
        fasta,
        output,
        haplotype,
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

fn parse_haplotype(raw: &str) -> Result<usize, ParseOutcome> {
    let digits: String = raw.chars().take_while(|c| c.is_ascii_digit()).collect();
    if digits.is_empty() {
        return Ok(1);
    }
    let value = digits
        .parse::<usize>()
        .map_err(|_| ParseOutcome::Error(format!("invalid haplotype '{raw}'")))?;
    if value == 0 {
        return Err(ParseOutcome::Error(format!("invalid haplotype '{raw}'")));
    }
    Ok(value)
}

fn run(args: &Args) -> io::Result<()> {
    let mut fasta_records = read_fasta(&args.fasta)?;
    let variants = read_variants(&args.input)?;
    apply_variants(&mut fasta_records, &variants, args.haplotype)?;

    match &args.output {
        Some(path) => {
            let mut out = File::create(path)?;
            write_fasta(&mut out, &fasta_records)
        }
        None => {
            let stdout = io::stdout();
            let mut out = stdout.lock();
            write_fasta(&mut out, &fasta_records)
        }
    }
}

fn read_fasta(path: &Path) -> io::Result<Vec<FastaRecord>> {
    let text = std::fs::read_to_string(path)?;
    let mut records = Vec::new();
    let mut current_header: Option<String> = None;
    let mut current_seq = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(header) = current_header.replace(rest.to_owned()) {
                records.push(fasta_record(header, std::mem::take(&mut current_seq)));
            }
        } else if !line.trim().is_empty() {
            current_seq.extend(line.trim().as_bytes().iter().map(u8::to_ascii_uppercase));
        }
    }
    if let Some(header) = current_header {
        records.push(fasta_record(header, current_seq));
    }
    Ok(records)
}

fn fasta_record(header: String, seq: Vec<u8>) -> FastaRecord {
    let (chrom, start) = parse_fasta_region_header(&header);
    FastaRecord {
        header,
        chrom,
        start,
        seq,
    }
}

fn parse_fasta_region_header(header: &str) -> (String, i64) {
    let name = header.split_whitespace().next().unwrap_or(header);
    if let Some((chrom, range)) = name.rsplit_once(':')
        && let Some((start, _end)) = range.split_once('-')
        && let Ok(start) = start.replace(',', "").parse::<i64>()
    {
        return (chrom.to_owned(), start);
    }
    (name.to_owned(), 1)
}

fn read_variants(path: &Path) -> io::Result<Vec<Variant>> {
    let text = read_vcf_text(path)?;
    let mut variants = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid VCF record with fewer than 8 columns: {line}"),
            ));
        }
        variants.push(Variant {
            chrom: fields[0].to_owned(),
            pos: fields[1].parse::<i64>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid VCF position '{}'", fields[1]),
                )
            })?,
            reference: fields[3].to_ascii_uppercase(),
            alts: fields[4].split(',').map(str::to_owned).collect(),
        });
    }
    Ok(variants)
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }

    let file = File::open(path)?;
    let mut text = String::new();
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let dec = MultiGzDecoder::new(file);
        let mut normalized = NormalizeFileformat::new(BufReader::new(dec))?;
        normalized.read_to_string(&mut text)?;
    } else {
        let mut normalized = NormalizeFileformat::new(BufReader::new(file))?;
        normalized.read_to_string(&mut text)?;
    }
    Ok(text)
}

fn apply_variants(
    fasta_records: &mut [FastaRecord],
    variants: &[Variant],
    haplotype: usize,
) -> io::Result<()> {
    for record in fasta_records {
        let mut local: Vec<&Variant> = variants
            .iter()
            .filter(|variant| variant.chrom == record.chrom)
            .collect();
        local.sort_by_key(|variant| std::cmp::Reverse(variant.pos));

        for variant in local {
            let Some(replacement) = replacement_allele(variant, haplotype) else {
                continue;
            };
            let Some(start) = variant.pos.checked_sub(record.start) else {
                continue;
            };
            if start < 0 {
                continue;
            }
            let start = start as usize;
            if start >= record.seq.len() {
                continue;
            }
            let reference_len = variant_reference_len(variant);
            let end = start.saturating_add(reference_len);
            if end > record.seq.len() {
                continue;
            }
            let observed = String::from_utf8_lossy(&record.seq[start..end]).to_ascii_uppercase();
            if observed != variant.reference {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "reference mismatch at {}:{}: expected {}, found {}",
                        variant.chrom, variant.pos, variant.reference, observed
                    ),
                ));
            }
            record.seq.splice(
                start..end,
                replacement.as_bytes().iter().map(u8::to_ascii_uppercase),
            );
        }
    }
    Ok(())
}

fn replacement_allele(variant: &Variant, haplotype: usize) -> Option<&str> {
    let alt = variant.alts.get(haplotype.saturating_sub(1))?;
    if alt == "." || alt == "*" {
        return None;
    }
    if alt.starts_with('<') && alt.ends_with('>') {
        return None;
    }
    Some(alt.as_str())
}

fn variant_reference_len(variant: &Variant) -> usize {
    variant.reference.len()
}

fn write_fasta<W: Write>(out: &mut W, records: &[FastaRecord]) -> io::Result<()> {
    for record in records {
        writeln!(out, ">{}", record.header)?;
        for chunk in record.seq.chunks(60) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_region_headers() {
        assert_eq!(
            parse_fasta_region_header("1:2-501 description"),
            ("1".to_owned(), 2)
        );
        assert_eq!(parse_fasta_region_header("chr1"), ("chr1".to_owned(), 1));
    }

    #[test]
    fn applies_variants_from_right_to_left() {
        let mut records = vec![FastaRecord {
            header: "chr1".into(),
            chrom: "chr1".into(),
            start: 1,
            seq: b"ACGTACGT".to_vec(),
        }];
        let variants = vec![
            Variant {
                chrom: "chr1".into(),
                pos: 2,
                reference: "C".into(),
                alts: vec!["TT".into()],
            },
            Variant {
                chrom: "chr1".into(),
                pos: 7,
                reference: "G".into(),
                alts: vec!["A".into()],
            },
        ];

        apply_variants(&mut records, &variants, 1).unwrap();
        assert_eq!(records[0].seq, b"ATTGTACAT");
    }
}
