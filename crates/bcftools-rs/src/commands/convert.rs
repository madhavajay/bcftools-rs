//! Focused `bcftools convert` implementation.
//!
//! Implemented: `--tsv2vcf` for explicit VCF-shaped TSV columns such as
//! `CHROM,POS,ID,REF,ALT`, FASTA-backed `AA` genotype columns for the common
//! 23andMe-style input shape, and a text-backed `--gvcf2vcf` block expander.
//! Deferred conversion modes remain tracked in TODO.md.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Read, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    build_bcf_csi_with_min_shift, build_vcf_csi_from_path_with_min_shift, build_vcf_tbi_from_path,
    write_csi, write_tbi,
};

use crate::diagnostics::fmt_etag;
use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::header_version::{build_lines, command_time};
use crate::io::{VariantOutputFormat, apply_verbosity, init_index2, write_index_parse};
use crate::reference::FastaReference;
use crate::tsv2vcf::{Tsv, TsvRecord};

const USAGE: &str = "\n\
About: Convert VCF/BCF files to different formats and back.\n\
Usage: bcftools convert [options]\n\
\n\
Options:\n\
       --tsv2vcf FILE             Convert TSV with CHROM/POS/REF/ALT columns to VCF\n\
       --gvcf2vcf                 Expand gVCF reference blocks to VCF records\n\
   -g, --gensample PREFIX         Convert VCF/BCF to PREFIX.gen.gz and PREFIX.samples\n\
   -G, --gensample2vcf            Convert GEN/SAMPLE input back to VCF\n\
       --hapsample PREFIX         Convert VCF/BCF to PREFIX.hap.gz and PREFIX.samples\n\
   -h, --haplegendsample PREFIX   Convert VCF/BCF to PREFIX.hap.gz, PREFIX.legend.gz, PREFIX.samples\n\
       --hapsample2vcf            Convert HAP/SAMPLE input back to VCF\n\
   -H, --haplegendsample2vcf      Convert HAP/LEGEND/SAMPLE input back to VCF\n\
       --3N6                      Use 3*N+6 GEN columns by adding CHROM first\n\
       --haploid2diploid          Convert haploid genotypes to diploid homozygotes for HAP/SAMPLE\n\
       --sex FILE                 Output sex column in HAP/SAMPLE sample file (Sample\\t[MF])\n\
       --tag STRING               Tag to take values for GEN output: GT,PL,GL,GP [GT]\n\
   -c, --columns STRING           Columns of the input TSV [ID,CHROM,POS,AA]\n\
   -f, --fasta-ref FILE           Reference sequence in FASTA format\n\
   -i, --include EXPR             Include only records for which EXPR is true (gVCF input)\n\
   -e, --exclude EXPR             Exclude records for which EXPR is true (gVCF input)\n\
       --no-version               Do not append version and command line to the header\n\
   -o, --output FILE              Write output to a file [standard output]\n\
   -O, --output-type u|b|v|z[0-9] u/b: BCF, v/z: un/compressed VCF [v]\n\
   -s, --samples LIST             Comma-separated sample names for trailing GT columns\n\
   -S, --samples-file FILE        File with sample names for trailing GT columns\n\
       --threads INT              Use multithreaded BGZF compression for compressed output\n\
   -v, --verbosity INT            Verbosity level\n\
       --vcf-ids                  Use VCF IDs in HAP/SAMPLE marker column\n\
       --keep-duplicates          Keep duplicate positions in GEN output\n\
   -W, --write-index[=FMT]        Automatically index VCF.gz output [off]\n\
\n";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

impl OutputKind {
    fn parse(raw: &str) -> Option<Self> {
        match raw.chars().next()? {
            'v' => Some(Self::VcfText),
            'z' | '0'..='9' => Some(Self::VcfGz),
            'b' | 'u' => Some(Self::Bcf),
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

#[derive(Debug)]
struct Args {
    mode: ConvertMode,
    columns: String,
    fasta_ref: Option<PathBuf>,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    samples: Vec<String>,
    include_expr: Option<String>,
    exclude_expr: Option<String>,
    vcf_ids: bool,
    gen_3n6: bool,
    keep_duplicates: bool,
    tag: String,
    haploid2diploid: bool,
    sex_file: Option<PathBuf>,
    no_version: bool,
    write_index: Option<i32>,
    thread_count: Option<NonZero<usize>>,
}

#[derive(Debug)]
enum ConvertMode {
    Tsv2Vcf(PathBuf),
    Gvcf2Vcf(PathBuf),
    VcfToGenSample { input: PathBuf, output: String },
    GenSampleToVcf(PathBuf),
    VcfToHapSample { input: PathBuf, output: String },
    HapSampleToVcf(PathBuf),
    VcfToHapLegendSample { input: PathBuf, output: String },
    HapLegendSampleToVcf(PathBuf),
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args, argv) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfconvert", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfconvert", &message));
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
    let mut tsv2vcf = None;
    let mut gvcf2vcf = false;
    let mut gensample = None;
    let mut gensample2vcf = false;
    let mut hapsample = None;
    let mut hapsample2vcf = false;
    let mut haplegendsample = None;
    let mut haplegendsample2vcf = false;
    let mut input = None;
    let mut columns = "ID,CHROM,POS,AA".to_owned();
    let mut fasta_ref = None;
    let mut output = None;
    let mut explicit_kind = None;
    let mut samples = Vec::new();
    let mut include_expr = None;
    let mut exclude_expr = None;
    let mut vcf_ids = false;
    let mut gen_3n6 = false;
    let mut keep_duplicates = false;
    let mut tag = "GT".to_owned();
    let mut haploid2diploid = false;
    let mut sex_file = None;
    let mut no_version = false;
    let mut write_index = None;
    let mut thread_count = None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "--help" | "-?" => return Err(ParseOutcome::Usage),
            "--no-version" => no_version = true,
            "--tsv2vcf" => tsv2vcf = Some(PathBuf::from(next_string(&mut iter, "--tsv2vcf")?)),
            "--gvcf2vcf" => gvcf2vcf = true,
            "-g" | "--gensample" => gensample = Some(next_string(&mut iter, "--gensample")?),
            "-G" | "--gensample2vcf" => gensample2vcf = true,
            "--hapsample" => hapsample = Some(next_string(&mut iter, "--hapsample")?),
            "--hapsample2vcf" => hapsample2vcf = true,
            "-h" | "--haplegendsample" => {
                haplegendsample = Some(next_string(&mut iter, "--haplegendsample")?)
            }
            "-H" | "--haplegendsample2vcf" => haplegendsample2vcf = true,
            "-c" | "--columns" => columns = next_string(&mut iter, "--columns")?,
            "-f" | "--fasta-ref" => {
                fasta_ref = Some(PathBuf::from(next_string(&mut iter, "--fasta-ref")?));
            }
            "-o" | "--output" => output = Some(PathBuf::from(next_string(&mut iter, "--output")?)),
            "-i" | "--include" => include_expr = Some(next_string(&mut iter, "--include")?),
            "-e" | "--exclude" => exclude_expr = Some(next_string(&mut iter, "--exclude")?),
            "-O" | "--output-type" => {
                explicit_kind = Some(parse_output_kind(&next_string(
                    &mut iter,
                    "--output-type",
                )?)?);
            }
            "-s" | "--samples" => {
                samples = parse_sample_list(&next_string(&mut iter, "--samples")?)
            }
            "-S" | "--samples-file" => {
                samples = read_sample_file(&next_string(&mut iter, "--samples-file")?)?
            }
            "--vcf-ids" => vcf_ids = true,
            "--3N6" => gen_3n6 = true,
            "--chrom" => {
                return Err(ParseOutcome::Error(
                    "The --chrom option has been deprecated, please use --3N6 instead".into(),
                ));
            }
            "--keep-duplicates" => keep_duplicates = true,
            "--tag" => tag = next_string(&mut iter, "--tag")?,
            "--haploid2diploid" => haploid2diploid = true,
            "--sex" => sex_file = Some(PathBuf::from(next_string(&mut iter, "--sex")?)),
            "-W" | "--write-index" => write_index = parse_write_index(None)?,
            "--threads" => thread_count = parse_threads(&next_string(&mut iter, "--threads")?)?,
            "-v" | "--verbosity" => {
                let value = next_string(&mut iter, "--verbosity")?;
                if apply_verbosity(&value).is_err() {
                    return Err(ParseOutcome::Error(format!(
                        "Could not parse argument: --verbosity {value}"
                    )));
                }
            }
            _ if raw.starts_with("--tsv2vcf=") => {
                tsv2vcf = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--gvcf2vcf=") => {
                gvcf2vcf = true;
                if input.is_some() {
                    return Err(ParseOutcome::Error(
                        "Unexpected positional argument: --gvcf2vcf input".into(),
                    ));
                }
                input = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--gensample=") => {
                gensample = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--hapsample=") => {
                hapsample = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--haplegendsample=") => {
                haplegendsample = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--sex=") => {
                sex_file = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--columns=") => columns = value_after_equals(&raw).to_owned(),
            _ if raw.starts_with("--fasta-ref=") => {
                fasta_ref = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--include=") => {
                include_expr = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--exclude=") => {
                exclude_expr = Some(value_after_equals(&raw).to_owned())
            }
            _ if raw.starts_with("--output-type=") => {
                explicit_kind = Some(parse_output_kind(value_after_equals(&raw))?)
            }
            _ if raw.starts_with("--samples=") => {
                samples = parse_sample_list(value_after_equals(&raw))
            }
            _ if raw.starts_with("--samples-file=") => {
                samples = read_sample_file(value_after_equals(&raw))?
            }
            _ if raw.starts_with("--tag=") => tag = value_after_equals(&raw).to_owned(),
            _ if raw.starts_with("--write-index=") => {
                write_index = parse_write_index(Some(value_after_equals(&raw)))?
            }
            _ if raw.starts_with("--threads=") => {
                thread_count = parse_threads(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                explicit_kind = Some(parse_output_kind(&raw[2..])?)
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-i") && raw.len() > 2 => include_expr = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-e") && raw.len() > 2 => exclude_expr = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-g") && raw.len() > 2 => gensample = Some(raw[2..].to_owned()),
            _ if raw.starts_with("-h") && raw.len() > 2 => {
                haplegendsample = Some(raw[2..].to_owned())
            }
            _ if raw.starts_with("-s") && raw.len() > 2 => samples = parse_sample_list(&raw[2..]),
            _ if raw.starts_with("-S") && raw.len() > 2 => samples = read_sample_file(&raw[2..])?,
            _ if raw.starts_with("-W=") => write_index = parse_write_index(Some(&raw[3..]))?,
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!("Unrecognized option: {raw}")));
            }
            _ => {
                if input.is_some() {
                    return Err(ParseOutcome::Error(format!(
                        "Unexpected positional argument: {raw}"
                    )));
                }
                input = Some(PathBuf::from(raw.as_ref()));
            }
        }
    }

    let mode_count = usize::from(tsv2vcf.is_some())
        + usize::from(gvcf2vcf)
        + usize::from(gensample.is_some())
        + usize::from(gensample2vcf)
        + usize::from(hapsample.is_some())
        + usize::from(hapsample2vcf)
        + usize::from(haplegendsample.is_some())
        + usize::from(haplegendsample2vcf);
    if mode_count == 0 {
        return Err(ParseOutcome::Usage);
    }
    if mode_count > 1 {
        return Err(ParseOutcome::Error(
            "Only one conversion mode can be specified".into(),
        ));
    }
    let mode = match (
        tsv2vcf,
        gvcf2vcf,
        gensample,
        gensample2vcf,
        hapsample,
        hapsample2vcf,
        haplegendsample,
        haplegendsample2vcf,
        input,
    ) {
        (Some(path), false, None, false, None, false, None, false, None) => {
            ConvertMode::Tsv2Vcf(path)
        }
        (None, true, None, false, None, false, None, false, Some(input)) => {
            ConvertMode::Gvcf2Vcf(input)
        }
        (None, false, Some(output), false, None, false, None, false, Some(input)) => {
            ConvertMode::VcfToGenSample { input, output }
        }
        (None, false, None, true, None, false, None, false, Some(input)) => {
            ConvertMode::GenSampleToVcf(input)
        }
        (None, false, None, false, Some(output), false, None, false, Some(input)) => {
            ConvertMode::VcfToHapSample { input, output }
        }
        (None, false, None, false, None, true, None, false, Some(input)) => {
            ConvertMode::HapSampleToVcf(input)
        }
        (None, false, None, false, None, false, Some(output), false, Some(input)) => {
            ConvertMode::VcfToHapLegendSample { input, output }
        }
        (None, false, None, false, None, false, None, true, Some(input)) => {
            ConvertMode::HapLegendSampleToVcf(input)
        }
        (None, true, None, false, None, false, None, false, None) => {
            ConvertMode::Gvcf2Vcf("-".into())
        }
        (None, false, Some(output), false, None, false, None, false, None) => {
            ConvertMode::VcfToGenSample {
                input: "-".into(),
                output,
            }
        }
        (None, false, None, true, None, false, None, false, None) => {
            return Err(ParseOutcome::Error(
                "--gensample2vcf requires GEN/SAMPLE input".into(),
            ));
        }
        (None, false, None, false, Some(output), false, None, false, None) => {
            ConvertMode::VcfToHapSample {
                input: "-".into(),
                output,
            }
        }
        (None, false, None, false, None, true, None, false, None) => {
            return Err(ParseOutcome::Error(
                "--hapsample2vcf requires HAP/SAMPLE input".into(),
            ));
        }
        (None, false, None, false, None, false, Some(output), false, None) => {
            ConvertMode::VcfToHapLegendSample {
                input: "-".into(),
                output,
            }
        }
        (None, false, None, false, None, false, None, true, None) => {
            return Err(ParseOutcome::Error(
                "--haplegendsample2vcf requires HAP/LEGEND/SAMPLE input".into(),
            ));
        }
        (Some(_), false, None, false, None, false, None, false, Some(_)) => {
            return Err(ParseOutcome::Error(
                "Unexpected positional argument for selected conversion mode".into(),
            ));
        }
        _ => unreachable!("conversion mode count checked above"),
    };
    let output_kind = explicit_kind.unwrap_or_else(|| {
        output
            .as_deref()
            .map(OutputKind::for_path)
            .unwrap_or(OutputKind::VcfText)
    });

    Ok(Args {
        mode,
        columns,
        fasta_ref,
        output,
        output_kind,
        samples,
        include_expr,
        exclude_expr,
        vcf_ids,
        gen_3n6,
        keep_duplicates,
        tag,
        haploid2diploid,
        sex_file,
        no_version,
        write_index,
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
    if args.write_index.is_some()
        && !matches!(args.output_kind, OutputKind::VcfGz | OutputKind::Bcf)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W requires compressed VCF or BCF output",
        ));
    }

    let mut buffer = Vec::new();
    match &args.mode {
        ConvertMode::Tsv2Vcf(path) => {
            if args.include_expr.is_some() || args.exclude_expr.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "-i/-e are only supported for VCF/gVCF input",
                ));
            }
            let tsv = Tsv::new(&args.columns);
            let has_aa = has_column(&tsv, "AA");
            validate_columns(&tsv, has_aa)?;
            let reference = if let Some(path) = &args.fasta_ref {
                Some(FastaReference::open(path)?)
            } else {
                None
            };
            if has_aa && reference.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--tsv2vcf requires the --fasta-ref option when AA is used",
                ));
            }
            let (records, stats) = read_tsv_records(path, &tsv, &args.samples, reference.as_ref())?;
            write_vcf_text(&mut buffer, args, argv, reference.as_ref(), &records)?;
            write_output(args, &buffer)?;
            if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
                write_index(path, args.output_kind, index_format)?;
            }
            print_stats(&stats, has_aa && !args.samples.is_empty());
            Ok(())
        }
        ConvertMode::Gvcf2Vcf(path) => {
            let reference_path = args.fasta_ref.as_ref().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--gvcf2vcf requires the --fasta-ref option",
                )
            })?;
            let reference = FastaReference::open(reference_path)?;
            expand_gvcf_to_vcf(path, args, argv, &reference, &mut buffer)?;
            write_output(args, &buffer)?;
            if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
                write_index(path, args.output_kind, index_format)?;
            }
            Ok(())
        }
        ConvertMode::VcfToGenSample { input, output } => {
            write_vcf_to_gensample(input, output, args)
        }
        ConvertMode::GenSampleToVcf(input) => {
            write_gensample_to_vcf(input, args, argv, &mut buffer)?;
            write_output(args, &buffer)?;
            if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
                write_index(path, args.output_kind, index_format)?;
            }
            Ok(())
        }
        ConvertMode::VcfToHapSample { input, output } => {
            write_vcf_to_hapsample(input, output, args)
        }
        ConvertMode::HapSampleToVcf(input) => {
            write_hapsample_to_vcf(input, args, argv, &mut buffer)?;
            write_output(args, &buffer)?;
            if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
                write_index(path, args.output_kind, index_format)?;
            }
            Ok(())
        }
        ConvertMode::VcfToHapLegendSample { input, output } => {
            write_vcf_to_haplegendsample(input, output, args)
        }
        ConvertMode::HapLegendSampleToVcf(input) => {
            write_haplegendsample_to_vcf(input, args, argv, &mut buffer)?;
            write_output(args, &buffer)?;
            if let (Some(index_format), Some(path)) = (args.write_index, args.output.as_deref()) {
                write_index(path, args.output_kind, index_format)?;
            }
            Ok(())
        }
    }
}

fn has_column(tsv: &Tsv, name: &str) -> bool {
    tsv.columns()
        .iter()
        .flatten()
        .any(|column| column.eq_ignore_ascii_case(name))
}

fn expand_gvcf_to_vcf<W: Write>(
    path: &Path,
    args: &Args,
    argv: &[OsString],
    reference: &FastaReference,
    out: &mut W,
) -> io::Result<()> {
    let text = read_vcf_text(path)?;
    let mut wrote_version = false;
    for line in text.lines() {
        if line.starts_with("##") {
            writeln!(out, "{line}")?;
            continue;
        }
        if line.starts_with("#CHROM") {
            if !args.no_version {
                let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
                prog_argv.extend(argv.iter().cloned());
                let lines = build_lines("bcftools_convert", &prog_argv, command_time());
                writeln!(out, "{}", lines.version_line)?;
                writeln!(out, "{}", lines.command_line)?;
            }
            writeln!(out, "{line}")?;
            wrote_version = true;
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        expand_gvcf_record(line, args, reference, out)?;
    }
    if !wrote_version && !args.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_convert", &prog_argv, command_time());
        writeln!(out, "{}", lines.version_line)?;
        writeln!(out, "{}", lines.command_line)?;
    }
    Ok(())
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
    let mut text = if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        read_gzip_text(path)?
    } else {
        fs::read_to_string(path)?
    };
    crate::vcf_compat::normalize_vcf_text(&mut text);
    Ok(text)
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-convert-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn read_text_auto_gzip(path: &Path) -> io::Result<String> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        read_gzip_text(path)
    } else {
        fs::read_to_string(path)
    }
}

fn read_gzip_text(path: &Path) -> io::Result<String> {
    let mut text = String::new();
    {
        let f = File::open(path)?;
        let mut dec = flate2::read::MultiGzDecoder::new(f);
        dec.read_to_string(&mut text)?;
    }
    Ok(text)
}

fn expand_gvcf_record<W: Write>(
    line: &str,
    args: &Args,
    reference: &FastaReference,
    out: &mut W,
) -> io::Result<()> {
    let mut fields: Vec<String> = line.split('\t').map(ToOwned::to_owned).collect();
    if fields.len() < 8 {
        writeln!(out, "{line}")?;
        return Ok(());
    }
    if !expression_pass(&fields, args)? {
        writeln!(out, "{line}")?;
        return Ok(());
    }
    if !is_gvcf_compatible_alt(&fields[4]) {
        writeln!(out, "{line}")?;
        return Ok(());
    }
    let Some(end) = info_end(&fields[7]) else {
        writeln!(out, "{line}")?;
        return Ok(());
    };
    let start = fields[1].parse::<i64>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid VCF POS '{}': {e}", fields[1]),
        )
    })?;
    fields[7] = remove_info_tag(&fields[7], "END");
    for pos in start..=end {
        fields[1] = pos.to_string();
        fields[3] = reference_base(reference, &fields[0], pos)?;
        writeln!(out, "{}", fields.join("\t"))?;
    }
    Ok(())
}

fn expression_pass(fields: &[String], args: &Args) -> io::Result<bool> {
    if let Some(expr) = &args.include_expr
        && !evaluate_expression(expr, fields)?.truthy()
    {
        return Ok(false);
    }
    if let Some(expr) = &args.exclude_expr
        && evaluate_expression(expr, fields)?.truthy()
    {
        return Ok(false);
    }
    Ok(true)
}

fn evaluate_expression(expr: &str, fields: &[String]) -> io::Result<FilterValue> {
    let context = EvalContext::new();
    bcffilter::eval_expression_with(expr, &context, |name, sample_index| {
        if sample_index.is_some() {
            return None;
        }
        super::filter::record_lookup(name, fields)
    })
}

fn is_gvcf_compatible_alt(raw: &str) -> bool {
    raw == "."
        || raw
            .split(',')
            .any(|alt| matches!(alt, "<*>" | "<X>" | "<NON_REF>"))
}

fn info_end(info: &str) -> Option<i64> {
    info.split(';').find_map(|item| {
        item.strip_prefix("END=")
            .and_then(|value| value.parse::<i64>().ok())
    })
}

fn remove_info_tag(info: &str, tag: &str) -> String {
    if info == "." || info.is_empty() {
        return ".".to_owned();
    }
    let kept: Vec<&str> = info
        .split(';')
        .filter(|item| *item != tag && !item.starts_with(&format!("{tag}=")))
        .collect();
    if kept.is_empty() {
        ".".to_owned()
    } else {
        kept.join(";")
    }
}

fn reference_base(reference: &FastaReference, chrom: &str, pos: i64) -> io::Result<String> {
    let region = format!("{chrom}:{pos}-{pos}");
    let bases = reference.fetch_region(&region)?;
    let base = bases
        .first()
        .copied()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty reference fetch"))?;
    Ok((base as char).to_ascii_uppercase().to_string())
}

fn write_vcf_to_gensample(input: &Path, output: &str, args: &Args) -> io::Result<()> {
    if args.write_index.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W is not supported with --gensample output",
        ));
    }
    if !matches!(args.tag.as_str(), "GT" | "GP" | "PL" | "GL") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("todo: --tag {}", args.tag),
        ));
    }
    let (gen_path, sample_path) = gensample_output_paths(output)?;
    let text = read_vcf_text(input)?;
    let header_samples = vcf_sample_names(&text)?;
    let selected_samples = selected_sample_indices(&header_samples, &args.samples)?;
    let samples: Vec<String> = selected_samples
        .iter()
        .map(|idx| header_samples[*idx].clone())
        .collect();
    let sex = match args.sex_file.as_deref() {
        Some(path) => Some(read_hapsample_sex_file(path, &samples)?),
        None => None,
    };
    if let Some(path) = sample_path.as_deref() {
        write_sample_output(path, &samples, sex.as_ref())?;
        eprintln!("Sample file: {}", path.display());
    }
    let Some(gen_path) = gen_path else {
        return Ok(());
    };
    let mut gen_text = String::new();
    let mut prev: Option<(String, String)> = None;
    let mut written = 0usize;
    let mut no_alt = 0usize;
    let mut non_biallelic = 0usize;
    let mut filtered = 0usize;
    let mut duplicated = 0usize;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(ToOwned::to_owned).collect();
        if fields.len() < 8 {
            continue;
        }
        if !expression_pass(&fields, args)? {
            filtered += 1;
            continue;
        }
        let raw_alt = fields.get(4).map(String::as_str).unwrap_or(".");
        if raw_alt == "." || raw_alt.is_empty() {
            no_alt += 1;
            continue;
        }
        if raw_alt.contains(',') {
            non_biallelic += 1;
            continue;
        }
        let key = (fields[0].clone(), fields[1].clone());
        if !args.keep_duplicates && prev.as_ref() == Some(&key) {
            duplicated += 1;
            continue;
        }
        prev = Some(key);
        let format_keys: Vec<&str> = fields
            .get(8)
            .map(|s| s.split(':').collect())
            .unwrap_or_default();
        let Some(tag_index) = format_keys.iter().position(|key| *key == args.tag) else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Error parsing {} tag at {}:{}",
                    args.tag, fields[0], fields[1]
                ),
            ));
        };
        let marker = format!("{}:{}_{}_{}", fields[0], fields[1], fields[3], raw_alt);
        if args.gen_3n6 {
            gen_text.push_str(&fields[0]);
            gen_text.push(' ');
        }
        gen_text.push_str(&marker);
        gen_text.push(' ');
        if args.vcf_ids {
            gen_text.push_str(fields.get(2).map(String::as_str).unwrap_or("."));
        } else {
            gen_text.push_str(&marker);
        }
        gen_text.push(' ');
        gen_text.push_str(&fields[1]);
        gen_text.push(' ');
        gen_text.push_str(&fields[3]);
        gen_text.push(' ');
        gen_text.push_str(raw_alt);
        for sample_idx in &selected_samples {
            let Some(sample) = fields.get(9 + *sample_idx) else {
                gen_text.push_str(" 0.33 0.33 0.33");
                continue;
            };
            let sample_values: Vec<&str> = sample.split(':').collect();
            let tag_value = sample_values.get(tag_index).copied().unwrap_or(".");
            let (aa, ab, bb) = tag_to_prob3(&args.tag, tag_value, &fields[0], &fields[1])?;
            gen_text.push(' ');
            gen_text.push_str(&aa);
            gen_text.push(' ');
            gen_text.push_str(&ab);
            gen_text.push(' ');
            gen_text.push_str(&bb);
        }
        gen_text.push('\n');
        written += 1;
    }
    write_bytes_auto_gzip(&gen_path, gen_text.as_bytes())?;
    eprintln!("Gen file: {}", gen_path.display());
    eprintln!(
        "{written} records written, {} skipped: {no_alt}/{non_biallelic}/{filtered}/{duplicated} no-ALT/non-biallelic/filtered/duplicated",
        no_alt + non_biallelic + filtered + duplicated,
    );
    Ok(())
}

fn write_gensample_to_vcf<W: Write>(
    input: &Path,
    args: &Args,
    argv: &[OsString],
    out: &mut W,
) -> io::Result<()> {
    if args.include_expr.is_some() || args.exclude_expr.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-i/-e are not supported with --gensample2vcf input",
        ));
    }
    let (gen_path, sample_path) = gensample_input_paths(input)?;
    let gen_text = read_text_auto_gzip(&gen_path)?;
    let samples = hapsample_sample_names(&sample_path)?;
    let records: Vec<GenSampleRecord> = gen_text
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(|line| parse_gensample_record(line, args.gen_3n6, args.vcf_ids, samples.len()))
        .collect::<io::Result<_>>()?;
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(
        out,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">"
    )?;
    for chrom in records
        .iter()
        .map(|record| record.chrom.as_str())
        .collect::<std::collections::BTreeSet<_>>()
    {
        writeln!(out, "##contig=<ID={chrom},length=2147483647>")?;
    }
    if !args.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_convert", &prog_argv, command_time());
        writeln!(out, "{}", lines.version_line)?;
        writeln!(out, "{}", lines.command_line)?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample in &samples {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;
    for record in &records {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT:GP",
            record.chrom,
            record.pos,
            record.id,
            record.ref_allele,
            record.alt_allele,
            record
                .end
                .map(|end| format!("END={end}"))
                .unwrap_or_else(|| ".".to_owned())
        )?;
        for (gt, gp) in &record.sample_values {
            write!(out, "\t{gt}:{gp}")?;
        }
        writeln!(out)?;
    }
    eprintln!("Number of processed rows: \t{}", records.len());
    Ok(())
}

#[derive(Debug)]
struct GenSampleRecord {
    chrom: String,
    pos: String,
    id: String,
    ref_allele: String,
    alt_allele: String,
    end: Option<i64>,
    sample_values: Vec<(String, String)>,
}

fn parse_gensample_record(
    line: &str,
    gen_3n6: bool,
    vcf_ids: bool,
    sample_count: usize,
) -> io::Result<GenSampleRecord> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    let min_fields = if gen_3n6 { 6 } else { 5 };
    if fields.len() < min_fields {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Error occurred while parsing: {line}"),
        ));
    }
    let (parsed, id, pos, ref_allele, alt_allele, prob_offset) = if gen_3n6 {
        (
            parse_chrom_pos_ref_alt(fields[1])?,
            if vcf_ids { fields[2] } else { "." },
            fields[3],
            fields[4],
            fields[5],
            6,
        )
    } else {
        match parse_chrom_pos_ref_alt(fields[0]) {
            Ok(parsed) => (
                parsed,
                if vcf_ids { fields[1] } else { "." },
                fields[2],
                fields[3],
                fields[4],
                5,
            ),
            Err(first_err) => {
                let parsed = parse_chrom_pos_ref_alt(fields[1]).map_err(|_| first_err)?;
                (
                    parsed,
                    if vcf_ids { fields[0] } else { "." },
                    fields[2],
                    fields[3],
                    fields[4],
                    5,
                )
            }
        }
    };
    let available_probs = fields.len().saturating_sub(prob_offset);
    if available_probs < sample_count * 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "expected {} genotype-probability fields, found {available_probs}",
                sample_count * 3
            ),
        ));
    }
    let sample_values = (0..sample_count)
        .map(|idx| {
            let offset = prob_offset + idx * 3;
            prob3_to_gt_gp(fields[offset], fields[offset + 1], fields[offset + 2])
        })
        .collect::<io::Result<_>>()?;
    Ok(GenSampleRecord {
        chrom: parsed.0,
        pos: pos.to_owned(),
        id: id.to_owned(),
        ref_allele: ref_allele.to_owned(),
        alt_allele: alt_allele.to_owned(),
        end: parsed.4,
        sample_values,
    })
}

fn prob3_to_gt_gp(aa: &str, ab: &str, bb: &str) -> io::Result<(String, String)> {
    let probs = [
        aa.parse::<f64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        ab.parse::<f64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
        bb.parse::<f64>()
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?,
    ];
    let mut max_idx = 0usize;
    let mut max_value = probs[0];
    for (idx, value) in probs.iter().copied().enumerate().skip(1) {
        if value > max_value {
            max_idx = idx;
            max_value = value;
        }
    }
    let gt = match max_idx {
        0 => "0/0",
        1 => "0/1",
        _ => "1/1",
    };
    Ok((
        gt.to_owned(),
        format!(
            "{},{},{}",
            normalize_probability_text(aa)?,
            normalize_probability_text(ab)?,
            normalize_probability_text(bb)?
        ),
    ))
}

fn normalize_probability_text(raw: &str) -> io::Result<String> {
    let value = raw
        .parse::<f64>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    Ok(value.to_string())
}

fn write_vcf_to_hapsample(input: &Path, output: &str, args: &Args) -> io::Result<()> {
    if args.write_index.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W is not supported with --hapsample output",
        ));
    }
    let (hap_path, sample_path) = hapsample_output_paths(output)?;
    let text = read_vcf_text(input)?;
    let header_samples = vcf_sample_names(&text)?;
    let selected_samples = selected_sample_indices(&header_samples, &args.samples)?;
    let samples: Vec<String> = selected_samples
        .iter()
        .map(|idx| header_samples[*idx].clone())
        .collect();
    let sex = match args.sex_file.as_deref() {
        Some(path) => Some(read_hapsample_sex_file(path, &samples)?),
        None => None,
    };
    if let Some(path) = sample_path.as_deref() {
        write_sample_output(path, &samples, sex.as_ref())?;
        eprintln!("Sample file: {}", path.display());
    }
    let Some(hap_path) = hap_path else {
        return Ok(());
    };
    let mut hap = String::new();
    let mut written = 0usize;
    let mut no_alt = 0usize;
    let mut non_biallelic = 0usize;
    let mut filtered = 0usize;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(ToOwned::to_owned).collect();
        if fields.len() < 8 {
            continue;
        }
        if !expression_pass(&fields, args)? {
            filtered += 1;
            continue;
        }
        let raw_alt = fields.get(4).map(String::as_str).unwrap_or(".");
        if raw_alt == "." || raw_alt.is_empty() {
            no_alt += 1;
            continue;
        }
        if raw_alt.contains(',') {
            non_biallelic += 1;
            continue;
        }
        let first_alt = raw_alt;
        let format_keys: Vec<&str> = fields
            .get(8)
            .map(|s| s.split(':').collect())
            .unwrap_or_default();
        let Some(gt_index) = format_keys.iter().position(|key| *key == "GT") else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FORMAT/GT tag not present at {}:{}", fields[0], fields[1]),
            ));
        };
        let marker = format!("{}:{}_{}_{}", fields[0], fields[1], fields[3], first_alt);
        if args.vcf_ids {
            hap.push_str(&format!(
                "{} {} {} {} {}",
                marker,
                fields.get(2).map(String::as_str).unwrap_or("."),
                fields[1],
                fields[3],
                first_alt
            ));
        } else {
            hap.push_str(&format!(
                "{} {} {} {} {}",
                fields[0], marker, fields[1], fields[3], first_alt
            ));
        }
        for sample_idx in &selected_samples {
            let Some(sample) = fields.get(9 + *sample_idx) else {
                continue;
            };
            let sample_values: Vec<&str> = sample.split(':').collect();
            let gt = sample_values.get(gt_index).copied().unwrap_or(".");
            let (a, b) = gt_to_hap_pair(gt, args.haploid2diploid);
            hap.push(' ');
            hap.push_str(&a);
            hap.push(' ');
            hap.push_str(&b);
        }
        hap.push('\n');
        written += 1;
    }
    write_bytes_auto_gzip(&hap_path, hap.as_bytes())?;
    eprintln!("Hap file: {}", hap_path.display());
    eprintln!(
        "{written} records written, {} skipped: {no_alt}/{non_biallelic}/{filtered} no-ALT/non-biallelic/filtered",
        no_alt + non_biallelic + filtered,
    );
    Ok(())
}

fn write_hapsample_to_vcf<W: Write>(
    input: &Path,
    args: &Args,
    argv: &[OsString],
    out: &mut W,
) -> io::Result<()> {
    if args.include_expr.is_some() || args.exclude_expr.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-i/-e are not supported with --hapsample2vcf input",
        ));
    }
    let (hap_path, sample_path) = hapsample_input_paths(input)?;
    let hap_text = read_text_auto_gzip(&hap_path)?;
    let samples = hapsample_sample_names(&sample_path)?;
    let records: Vec<HapSampleRecord> = hap_text
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .map(|line| parse_hapsample_record(line, args.vcf_ids, samples.len()))
        .collect::<io::Result<_>>()?;
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(
        out,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    for chrom in records
        .iter()
        .map(|record| record.chrom.as_str())
        .collect::<std::collections::BTreeSet<_>>()
    {
        writeln!(out, "##contig=<ID={chrom},length=2147483647>")?;
    }
    if !args.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_convert", &prog_argv, command_time());
        writeln!(out, "{}", lines.version_line)?;
        writeln!(out, "{}", lines.command_line)?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample in &samples {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;
    for record in &records {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT",
            record.chrom,
            record.pos,
            record.id,
            record.ref_allele,
            record.alt_allele,
            record
                .end
                .map(|end| format!("END={end}"))
                .unwrap_or_else(|| ".".to_owned())
        )?;
        for gt in &record.genotypes {
            write!(out, "\t{gt}")?;
        }
        writeln!(out)?;
    }
    eprintln!("Number of processed rows: \t{}", records.len());
    Ok(())
}

#[derive(Debug)]
struct HapSampleRecord {
    chrom: String,
    pos: String,
    id: String,
    ref_allele: String,
    alt_allele: String,
    end: Option<i64>,
    genotypes: Vec<String>,
}

fn parse_hapsample_record(
    line: &str,
    vcf_ids: bool,
    sample_count: usize,
) -> io::Result<HapSampleRecord> {
    let fields: Vec<&str> = line.split_whitespace().collect();
    if fields.len() < 5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Error occurred while parsing: {line}"),
        ));
    }
    let (marker, id, pos, ref_allele, alt_allele, hap_offset) = if vcf_ids {
        (
            fields[0],
            fields.get(1).copied().unwrap_or("."),
            fields.get(2).copied().unwrap_or("."),
            fields.get(3).copied().unwrap_or("."),
            fields.get(4).copied().unwrap_or("."),
            5,
        )
    } else {
        (
            fields.get(1).copied().unwrap_or(fields[0]),
            ".",
            fields.get(2).copied().unwrap_or("."),
            fields.get(3).copied().unwrap_or("."),
            fields.get(4).copied().unwrap_or("."),
            5,
        )
    };
    let parsed = parse_chrom_pos_ref_alt(marker)?;
    let chrom = parsed.0;
    let end = parsed.4;
    let available_haps = fields.len().saturating_sub(hap_offset);
    if available_haps < sample_count * 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "expected {} haplotype fields, found {available_haps}",
                sample_count * 2
            ),
        ));
    }
    let genotypes = (0..sample_count)
        .map(|idx| {
            haps_to_gt(
                fields[hap_offset + idx * 2],
                fields[hap_offset + idx * 2 + 1],
            )
        })
        .collect();
    Ok(HapSampleRecord {
        chrom,
        pos: pos.to_owned(),
        id: id.to_owned(),
        ref_allele: ref_allele.to_owned(),
        alt_allele: alt_allele.to_owned(),
        end,
        genotypes,
    })
}

fn parse_chrom_pos_ref_alt(raw: &str) -> io::Result<(String, String, String, String, Option<i64>)> {
    let (chrom, rest) = raw.split_once(':').ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Could not determine CHROM in {raw}"),
        )
    })?;
    let mut parts = rest.splitn(4, '_');
    let pos = parts.next().unwrap_or_default();
    let ref_allele = parts.next().unwrap_or_default();
    let alt_allele = parts.next().unwrap_or_default();
    if pos.is_empty() || ref_allele.is_empty() || alt_allele.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Could not parse the CHROM:POS_REF_ALT[_END] string: {raw}"),
        ));
    }
    let end = parts.next().and_then(|value| value.parse::<i64>().ok());
    Ok((
        chrom.to_owned(),
        pos.to_owned(),
        ref_allele.to_owned(),
        alt_allele.to_owned(),
        end,
    ))
}

fn haps_to_gt(first: &str, second: &str) -> String {
    let first_unphased = first.ends_with('*');
    let second_unphased = second.ends_with('*');
    let first = first.trim_end_matches('*');
    let second = second.trim_end_matches('*');
    if second == "-" {
        return if first == "?" {
            ".".to_owned()
        } else {
            first.to_owned()
        };
    }
    if first == "?" && second == "?" {
        return ".|.".to_owned();
    }
    let left = if first == "?" { "." } else { first };
    let right = if second == "?" { "." } else { second };
    let sep = if first_unphased || second_unphased {
        '/'
    } else {
        '|'
    };
    format!("{left}{sep}{right}")
}

fn write_vcf_to_haplegendsample(input: &Path, output: &str, args: &Args) -> io::Result<()> {
    if args.write_index.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-W is not supported with --haplegendsample output",
        ));
    }
    let (hap_path, legend_path, sample_path) = haplegendsample_output_paths(output)?;
    let text = read_vcf_text(input)?;
    let header_samples = vcf_sample_names(&text)?;
    let selected_samples = selected_sample_indices(&header_samples, &args.samples)?;
    let samples: Vec<String> = selected_samples
        .iter()
        .map(|idx| header_samples[*idx].clone())
        .collect();
    let sex = match args.sex_file.as_deref() {
        Some(path) => Some(read_hapsample_sex_file(path, &samples)?),
        None => None,
    };
    if let Some(path) = sample_path.as_deref() {
        write_haplegendsample_sample_output(path, &samples, sex.as_ref())?;
        eprintln!("Sample file: {}", path.display());
    }
    let mut hap = String::new();
    let mut legend = String::from("id position a0 a1\n");
    let mut written = 0usize;
    let mut no_alt = 0usize;
    let mut non_biallelic = 0usize;
    let mut filtered = 0usize;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(ToOwned::to_owned).collect();
        if fields.len() < 8 {
            continue;
        }
        if !expression_pass(&fields, args)? {
            filtered += 1;
            continue;
        }
        let raw_alt = fields.get(4).map(String::as_str).unwrap_or(".");
        if raw_alt == "." || raw_alt.is_empty() {
            no_alt += 1;
            continue;
        }
        if raw_alt.contains(',') {
            non_biallelic += 1;
            continue;
        }
        let format_keys: Vec<&str> = fields
            .get(8)
            .map(|s| s.split(':').collect())
            .unwrap_or_default();
        let Some(gt_index) = format_keys.iter().position(|key| *key == "GT") else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("FORMAT/GT tag not present at {}:{}", fields[0], fields[1]),
            ));
        };
        let mut first = true;
        for sample_idx in &selected_samples {
            let Some(sample) = fields.get(9 + *sample_idx) else {
                continue;
            };
            let sample_values: Vec<&str> = sample.split(':').collect();
            let gt = sample_values.get(gt_index).copied().unwrap_or(".");
            let (a, b) = gt_to_hap_pair(gt, args.haploid2diploid);
            if !first {
                hap.push(' ');
            }
            first = false;
            hap.push_str(&a);
            hap.push(' ');
            hap.push_str(&b);
        }
        hap.push('\n');
        let legend_id = if args.vcf_ids && fields.get(2).is_some_and(|id| id != ".") {
            fields[2].clone()
        } else {
            format!("{}:{}_{}_{}", fields[0], fields[1], fields[3], raw_alt)
        };
        legend.push_str(&format!(
            "{} {} {} {}\n",
            legend_id, fields[1], fields[3], raw_alt
        ));
        written += 1;
    }
    if let Some(path) = hap_path.as_deref() {
        write_bytes_auto_gzip(path, hap.as_bytes())?;
        eprintln!("Hap file: {}", path.display());
    }
    if let Some(path) = legend_path.as_deref() {
        write_bytes_auto_gzip(path, legend.as_bytes())?;
        eprintln!("Legend file: {}", path.display());
    }
    eprintln!(
        "{written} records written, {} skipped: {no_alt}/{non_biallelic}/{filtered} no-ALT/non-biallelic/filtered",
        no_alt + non_biallelic + filtered,
    );
    Ok(())
}

fn write_haplegendsample_to_vcf<W: Write>(
    input: &Path,
    args: &Args,
    argv: &[OsString],
    out: &mut W,
) -> io::Result<()> {
    if args.vcf_ids {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "The option --haplegendsample2vcf cannot be combined with --vcf-ids",
        ));
    }
    if args.include_expr.is_some() || args.exclude_expr.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "-i/-e are not supported with --haplegendsample2vcf input",
        ));
    }
    let (hap_path, legend_path, sample_path) = haplegendsample_input_paths(input)?;
    let hap_text = read_text_auto_gzip(&hap_path)?;
    let legend_text = read_text_auto_gzip(&legend_path)?;
    let samples = haplegendsample_sample_names(&sample_path)?;
    let hap_lines: Vec<&str> = hap_text
        .lines()
        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
        .collect();
    let legend_lines: Vec<&str> = legend_text
        .lines()
        .skip(1)
        .filter(|line| !line.trim().is_empty())
        .collect();
    if hap_lines.len() != legend_lines.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Different number of records in {} and {}?",
                legend_path.display(),
                hap_path.display()
            ),
        ));
    }
    let records: Vec<HapSampleRecord> = legend_lines
        .iter()
        .zip(hap_lines.iter())
        .map(|(legend, hap)| parse_haplegend_record(legend, hap, samples.len()))
        .collect::<io::Result<_>>()?;
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(
        out,
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">"
    )?;
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    for chrom in records
        .iter()
        .map(|record| record.chrom.as_str())
        .collect::<std::collections::BTreeSet<_>>()
    {
        writeln!(out, "##contig=<ID={chrom},length=2147483647>")?;
    }
    if !args.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_convert", &prog_argv, command_time());
        writeln!(out, "{}", lines.version_line)?;
        writeln!(out, "{}", lines.command_line)?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")?;
    for sample in &samples {
        write!(out, "\t{sample}")?;
    }
    writeln!(out)?;
    for record in &records {
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t.\t.\t{}\tGT",
            record.chrom,
            record.pos,
            record.id,
            record.ref_allele,
            record.alt_allele,
            record
                .end
                .map(|end| format!("END={end}"))
                .unwrap_or_else(|| ".".to_owned())
        )?;
        for gt in &record.genotypes {
            write!(out, "\t{gt}")?;
        }
        writeln!(out)?;
    }
    eprintln!("Number of processed rows: \t{}", records.len());
    Ok(())
}

fn parse_haplegend_record(
    legend_line: &str,
    hap_line: &str,
    sample_count: usize,
) -> io::Result<HapSampleRecord> {
    let legend_fields: Vec<&str> = legend_line.split_whitespace().collect();
    if legend_fields.len() < 4 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Error occurred while parsing legend row: {legend_line}"),
        ));
    }
    let marker = legend_fields[0];
    let parsed = parse_chrom_pos_ref_alt(marker)?;
    if legend_fields[1] != parsed.1 || legend_fields[2] != parsed.2 || legend_fields[3] != parsed.3
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("REF/ALT/POS mismatch in legend row: {legend_line}"),
        ));
    }
    let hap_fields: Vec<&str> = hap_line.split_whitespace().collect();
    if hap_fields.len() < sample_count * 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "expected {} haplotype fields, found {}",
                sample_count * 2,
                hap_fields.len()
            ),
        ));
    }
    let genotypes = (0..sample_count)
        .map(|idx| haps_to_gt(hap_fields[idx * 2], hap_fields[idx * 2 + 1]))
        .collect();
    Ok(HapSampleRecord {
        chrom: parsed.0,
        pos: parsed.1,
        id: ".".to_owned(),
        ref_allele: parsed.2,
        alt_allele: parsed.3,
        end: parsed.4,
        genotypes,
    })
}

fn gensample_output_paths(raw: &str) -> io::Result<(Option<PathBuf>, Option<PathBuf>)> {
    let parts: Vec<&str> = raw.split(',').collect();
    match parts.as_slice() {
        [prefix] => Ok((
            Some(PathBuf::from(format!("{prefix}.gen.gz"))),
            Some(PathBuf::from(format!("{prefix}.samples"))),
        )),
        [gen_out, sample] => Ok((optional_output_path(gen_out), optional_output_path(sample))),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Error parsing --gensample filenames: {raw}"),
        )),
    }
}

fn gensample_input_paths(raw: &Path) -> io::Result<(PathBuf, PathBuf)> {
    let raw = raw.as_os_str().to_string_lossy();
    let parts: Vec<&str> = raw.split(',').collect();
    match parts.as_slice() {
        [prefix] => Ok((
            PathBuf::from(format!("{prefix}.gen.gz")),
            PathBuf::from(format!("{prefix}.samples")),
        )),
        [gen_out, sample] if !gen_out.is_empty() && !sample.is_empty() => {
            Ok((PathBuf::from(gen_out), PathBuf::from(sample)))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Error parsing --gensample2vcf filenames: {raw}"),
        )),
    }
}

fn hapsample_output_paths(raw: &str) -> io::Result<(Option<PathBuf>, Option<PathBuf>)> {
    let parts: Vec<&str> = raw.split(',').collect();
    match parts.as_slice() {
        [prefix] => Ok((
            Some(PathBuf::from(format!("{prefix}.hap.gz"))),
            Some(PathBuf::from(format!("{prefix}.samples"))),
        )),
        [hap, sample] => Ok((optional_output_path(hap), optional_output_path(sample))),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Error parsing --hapsample filenames: {raw}"),
        )),
    }
}

fn hapsample_input_paths(raw: &Path) -> io::Result<(PathBuf, PathBuf)> {
    let raw = raw.as_os_str().to_string_lossy();
    let parts: Vec<&str> = raw.split(',').collect();
    match parts.as_slice() {
        [prefix] => Ok((
            PathBuf::from(format!("{prefix}.hap.gz")),
            PathBuf::from(format!("{prefix}.samples")),
        )),
        [hap, sample] if !hap.is_empty() && !sample.is_empty() => {
            Ok((PathBuf::from(hap), PathBuf::from(sample)))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Error parsing --hapsample2vcf filenames: {raw}"),
        )),
    }
}

fn haplegendsample_output_paths(
    raw: &str,
) -> io::Result<(Option<PathBuf>, Option<PathBuf>, Option<PathBuf>)> {
    let parts: Vec<&str> = raw.split(',').collect();
    match parts.as_slice() {
        [prefix] => Ok((
            Some(PathBuf::from(format!("{prefix}.hap.gz"))),
            Some(PathBuf::from(format!("{prefix}.legend.gz"))),
            Some(PathBuf::from(format!("{prefix}.samples"))),
        )),
        [hap, legend, sample] => Ok((
            optional_output_path(hap),
            optional_output_path(legend),
            optional_output_path(sample),
        )),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Error parsing --haplegendsample filenames: {raw}"),
        )),
    }
}

fn haplegendsample_input_paths(raw: &Path) -> io::Result<(PathBuf, PathBuf, PathBuf)> {
    let raw = raw.as_os_str().to_string_lossy();
    let parts: Vec<&str> = raw.split(',').collect();
    match parts.as_slice() {
        [prefix] => Ok((
            PathBuf::from(format!("{prefix}.hap.gz")),
            PathBuf::from(format!("{prefix}.legend.gz")),
            PathBuf::from(format!("{prefix}.samples")),
        )),
        [hap, legend, sample] if !hap.is_empty() && !legend.is_empty() && !sample.is_empty() => {
            Ok((
                PathBuf::from(hap),
                PathBuf::from(legend),
                PathBuf::from(sample),
            ))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Could not parse hap/legend/sample file names: {raw}"),
        )),
    }
}

fn optional_output_path(raw: &str) -> Option<PathBuf> {
    if raw.is_empty() || raw == "." {
        None
    } else {
        Some(PathBuf::from(raw))
    }
}

fn vcf_sample_names(text: &str) -> io::Result<Vec<String>> {
    let header = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "missing #CHROM header"))?;
    Ok(header.split('\t').skip(9).map(ToOwned::to_owned).collect())
}

fn hapsample_sample_names(path: &Path) -> io::Result<Vec<String>> {
    let text = read_text_auto_gzip(path)?;
    Ok(text
        .lines()
        .skip(2)
        .filter_map(|line| line.split_whitespace().next())
        .map(str::to_owned)
        .collect())
}

fn haplegendsample_sample_names(path: &Path) -> io::Result<Vec<String>> {
    let text = read_text_auto_gzip(path)?;
    Ok(text
        .lines()
        .skip(1)
        .filter_map(|line| line.split_whitespace().next())
        .map(str::to_owned)
        .collect())
}

fn selected_sample_indices(
    header_samples: &[String],
    requested: &[String],
) -> io::Result<Vec<usize>> {
    if requested.is_empty() {
        return Ok((0..header_samples.len()).collect());
    }
    let exclude = requested
        .first()
        .is_some_and(|sample| sample.strip_prefix('^').is_some());
    let requested_names: std::collections::BTreeSet<&str> = requested
        .iter()
        .map(|sample| sample.strip_prefix('^').unwrap_or(sample.as_str()))
        .collect();
    for name in &requested_names {
        if !header_samples.iter().any(|sample| sample == name) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Sample not found: {name}"),
            ));
        }
    }
    let indices = header_samples
        .iter()
        .enumerate()
        .filter_map(|(idx, sample)| {
            let listed = requested_names.contains(sample.as_str());
            if listed ^ exclude { Some(idx) } else { None }
        })
        .collect();
    Ok(indices)
}

fn read_hapsample_sex_file(
    path: &Path,
    samples: &[String],
) -> io::Result<std::collections::BTreeMap<String, char>> {
    let text = fs::read_to_string(path)?;
    let mut sex = std::collections::BTreeMap::new();
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut fields = trimmed.split_whitespace();
        let Some(sample) = fields.next() else {
            continue;
        };
        let Some(raw_sex) = fields.next() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Could not parse {}: {line}", path.display()),
            ));
        };
        let code = match raw_sex {
            "M" => '1',
            "F" => '2',
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Could not parse {}: {line}", path.display()),
                ));
            }
        };
        sex.insert(sample.to_owned(), code);
    }
    for sample in samples {
        if !sex.contains_key(sample) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Missing sex for sample {sample} in {}", path.display()),
            ));
        }
    }
    Ok(sex)
}

fn write_sample_output(
    path: &Path,
    samples: &[String],
    sex: Option<&std::collections::BTreeMap<String, char>>,
) -> io::Result<()> {
    let mut out = if sex.is_some() {
        String::from("ID_1 ID_2 missing sex\n0 0 0 0\n")
    } else {
        String::from("ID_1 ID_2 missing\n0 0 0\n")
    };
    for sample in samples {
        match sex.and_then(|values| values.get(sample)) {
            Some(code) => {
                out.push_str(sample);
                out.push(' ');
                out.push_str(sample);
                out.push_str(" 0 ");
                out.push(*code);
                out.push('\n');
            }
            None => {
                out.push_str(sample);
                out.push(' ');
                out.push_str(sample);
                out.push_str(" 0\n");
            }
        }
    }
    write_bytes_auto_gzip(path, out.as_bytes())
}

fn write_haplegendsample_sample_output(
    path: &Path,
    samples: &[String],
    sex: Option<&std::collections::BTreeMap<String, char>>,
) -> io::Result<()> {
    let mut out = String::from("sample population group sex\n");
    for sample in samples {
        let code = sex
            .and_then(|values| values.get(sample))
            .copied()
            .unwrap_or('2');
        out.push_str(sample);
        out.push(' ');
        out.push_str(sample);
        out.push(' ');
        out.push_str(sample);
        out.push(' ');
        out.push(code);
        out.push('\n');
    }
    write_bytes_auto_gzip(path, out.as_bytes())
}

fn write_bytes_auto_gzip(path: &Path, bytes: &[u8]) -> io::Result<()> {
    if path == Path::new("-") {
        io::stdout().lock().write_all(bytes)?;
        return Ok(());
    }
    if path
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
    {
        let file = File::create(path)?;
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        encoder.write_all(bytes)?;
        encoder.finish()?;
    } else {
        fs::write(path, bytes)?;
    }
    Ok(())
}

fn gt_to_hap_pair(raw: &str, haploid2diploid: bool) -> (String, String) {
    let gt = raw.split(':').next().unwrap_or(raw);
    if gt == "." {
        return (
            "?".to_owned(),
            if haploid2diploid {
                "?".to_owned()
            } else {
                "-".to_owned()
            },
        );
    }
    if gt == "./." || gt == ".|." || gt.contains('.') {
        return ("?".to_owned(), "?".to_owned());
    }
    let phased = gt.contains('|');
    let mut alleles = gt.split(['/', '|']);
    let Some(first) = alleles.next() else {
        return ("?".to_owned(), "?".to_owned());
    };
    let second = alleles.next();
    if second.is_none() {
        return (
            first.to_owned(),
            if haploid2diploid {
                first.to_owned()
            } else {
                "-".to_owned()
            },
        );
    }
    let left = if phased {
        first.to_owned()
    } else {
        format!("{first}*")
    };
    let right = match second {
        Some(value) if phased => value.to_owned(),
        Some(value) => format!("{value}*"),
        None => unreachable!("haploid GT handled above"),
    };
    (left, right)
}

fn tag_to_prob3(
    tag: &str,
    raw: &str,
    chrom: &str,
    pos: &str,
) -> io::Result<(String, String, String)> {
    match tag {
        "GT" => {
            let (aa, ab, bb) = gt_to_prob3(raw);
            Ok((aa.to_owned(), ab.to_owned(), bb.to_owned()))
        }
        "GP" => gp_to_prob3(raw, chrom, pos),
        "PL" => pl_to_prob3(raw),
        "GL" => gl_to_prob3(raw),
        _ => unreachable!("tag validated before conversion"),
    }
}

fn gt_to_prob3(raw: &str) -> (&'static str, &'static str, &'static str) {
    let gt = raw.split(':').next().unwrap_or(raw);
    if gt == "." || gt == "./." || gt == ".|." || gt.contains('.') {
        if gt.contains('/') || gt.contains('|') {
            return ("0.33", "0.33", "0.33");
        }
        return ("0.5", "0.0", "0.5");
    }
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    match alleles.as_slice() {
        [single] => {
            if *single == "1" {
                ("0", "0", "1")
            } else {
                ("1", "0", "0")
            }
        }
        [a, b] if a != b => ("0", "1", "0"),
        [a, _] if *a == "1" => ("0", "0", "1"),
        [_, _] => ("1", "0", "0"),
        _ => ("0.33", "0.33", "0.33"),
    }
}

fn gp_to_prob3(raw: &str, chrom: &str, pos: &str) -> io::Result<(String, String, String)> {
    let values = parse_float_vector(raw)?;
    if values.len() == 2 {
        return Ok((
            format_prob(values[0]),
            format_prob(0.0),
            format_prob(values[1]),
        ));
    }
    let aa = *values.first().unwrap_or(&0.0);
    let ab = *values.get(1).unwrap_or(&0.0);
    let bb = *values.get(2).unwrap_or(&0.0);
    for value in [aa, ab, bb] {
        if !(0.0..=1.0).contains(&value) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("[{chrom}:{pos}:{value}] GP value outside range [0,1]"),
            ));
        }
    }
    Ok((format_prob(aa), format_prob(ab), format_prob(bb)))
}

fn pl_to_prob3(raw: &str) -> io::Result<(String, String, String)> {
    let values = parse_float_vector(raw)?;
    let weights: Vec<f32> = values
        .iter()
        .take(3)
        .map(|value| 10f32.powf(-0.1 * *value as f32))
        .collect();
    normalize_likelihood_weights(&weights)
}

fn gl_to_prob3(raw: &str) -> io::Result<(String, String, String)> {
    let values = parse_float_vector(raw)?;
    let weights: Vec<f32> = values
        .iter()
        .take(3)
        .map(|value| 10f32.powf(*value as f32))
        .collect();
    normalize_likelihood_weights(&weights)
}

fn normalize_likelihood_weights(weights: &[f32]) -> io::Result<(String, String, String)> {
    let sum: f32 = weights.iter().sum();
    if sum == 0.0 || weights.is_empty() {
        return Ok(("0.33".to_owned(), "0.33".to_owned(), "0.33".to_owned()));
    }
    if weights.len() == 2 {
        return Ok((
            format_prob(weights[0] / sum),
            "0".to_owned(),
            format_prob(weights[1] / sum),
        ));
    }
    Ok((
        format_prob(weights.first().copied().unwrap_or(0.0) / sum),
        format_prob(weights.get(1).copied().unwrap_or(0.0) / sum),
        format_prob(weights.get(2).copied().unwrap_or(0.0) / sum),
    ))
}

fn parse_float_vector(raw: &str) -> io::Result<Vec<f64>> {
    if raw == "." || raw.is_empty() {
        return Ok(Vec::new());
    }
    raw.split(',')
        .map(|value| {
            if value == "." {
                Ok(0.0)
            } else {
                value
                    .parse::<f64>()
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            }
        })
        .collect()
}

fn format_prob<T: Into<f64>>(value: T) -> String {
    let value = value.into();
    format!("{value:.6}")
}

fn validate_columns(tsv: &Tsv, has_aa: bool) -> io::Result<()> {
    let mut has_chrom = false;
    let mut has_pos = false;
    let mut has_ref = false;
    let mut has_alt = false;
    for column in tsv.columns().iter().flatten() {
        match column.to_ascii_uppercase().as_str() {
            "CHROM" => has_chrom = true,
            "POS" => has_pos = true,
            "REF" => has_ref = true,
            "ALT" => has_alt = true,
            "AA" => {}
            _ => {}
        }
    }
    if has_chrom && has_pos && (has_aa || (has_ref && has_alt)) {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--tsv2vcf requires CHROM, POS, and either AA or REF/ALT columns",
        ))
    }
}

#[derive(Debug)]
struct ConvertedRecord {
    record: TsvRecord,
    genotypes: Vec<String>,
}

#[derive(Debug, Default)]
struct ConversionStats {
    total: usize,
    skipped: usize,
    written: usize,
    missing: usize,
    hom_rr: usize,
    het_ra: usize,
    hom_aa: usize,
    het_aa: usize,
}

fn read_tsv_records(
    path: &Path,
    tsv: &Tsv,
    samples: &[String],
    reference: Option<&FastaReference>,
) -> io::Result<(Vec<ConvertedRecord>, ConversionStats)> {
    let file = File::open(path)?;
    let mut records = Vec::new();
    let mut stats = ConversionStats::default();
    for (line_index, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let line_number = line_index + 1;
        stats.total += 1;
        let mut record = match tsv.parse_record(trimmed) {
            Ok(record) => record,
            Err(e) => {
                warn_skip_tsv_line(line_number, &e);
                stats.skipped += 1;
                continue;
            }
        };
        let genotypes = if has_column(tsv, "AA") {
            match populate_from_aa(trimmed, tsv, &mut record, samples, reference, &mut stats) {
                Ok(Some(genotypes)) => genotypes,
                Ok(None) => {
                    stats.skipped += 1;
                    continue;
                }
                Err(e) => {
                    warn_skip_tsv_line(line_number, &e);
                    stats.skipped += 1;
                    continue;
                }
            }
        } else {
            if let Err(e) = validate_record_ready(&record) {
                warn_skip_tsv_line(line_number, &e);
                stats.skipped += 1;
                continue;
            }
            match parse_trailing_genotypes(trimmed, tsv, &record, samples) {
                Ok(genotypes) => genotypes,
                Err(e) => {
                    warn_skip_tsv_line(line_number, &e);
                    stats.skipped += 1;
                    continue;
                }
            }
        };
        if let Err(e) = validate_record_ready(&record) {
            warn_skip_tsv_line(line_number, &e);
            stats.skipped += 1;
            continue;
        }
        records.push(ConvertedRecord { record, genotypes });
        stats.written += 1;
    }
    Ok((records, stats))
}

fn warn_skip_tsv_line(line_number: usize, err: &io::Error) {
    eprintln!("Warning: skipping malformed TSV line {line_number}: {err}");
}

fn validate_record_ready(record: &TsvRecord) -> io::Result<()> {
    required(record.chrom.as_deref(), "CHROM")?;
    record.pos.ok_or_else(|| missing("POS"))?;
    required(record.ref_allele.as_deref(), "REF")?;
    Ok(())
}

fn print_stats(stats: &ConversionStats, include_genotype_counts: bool) {
    eprintln!("Rows total: \t{}", stats.total);
    eprintln!("Rows skipped: \t{}", stats.skipped);
    eprintln!("Sites written: \t{}", stats.written);
    if include_genotype_counts {
        eprintln!("Missing GTs: \t{}", stats.missing);
        eprintln!("Hom RR: \t{}", stats.hom_rr);
        eprintln!("Het RA: \t{}", stats.het_ra);
        eprintln!("Hom AA: \t{}", stats.hom_aa);
        eprintln!("Het AA: \t{}", stats.het_aa);
    }
}

fn write_vcf_text<W: Write>(
    out: &mut W,
    args: &Args,
    argv: &[OsString],
    reference: Option<&FastaReference>,
    records: &[ConvertedRecord],
) -> io::Result<()> {
    writeln!(out, "##fileformat=VCFv4.2")?;
    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
    if let Some(reference) = reference {
        write_contig_headers(out, reference)?;
    } else {
        write_record_contig_headers(out, records)?;
    }
    writeln!(
        out,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    if !args.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_convert", &prog_argv, command_time());
        writeln!(out, "{}", lines.version_line)?;
        writeln!(out, "{}", lines.command_line)?;
    }
    write!(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")?;
    if !args.samples.is_empty() {
        write!(out, "\tFORMAT")?;
        for sample in &args.samples {
            write!(out, "\t{sample}")?;
        }
    }
    writeln!(out)?;
    for converted in records {
        let record = &converted.record;
        write!(
            out,
            "{}\t{}\t{}\t{}\t{}\t.\t.\t.",
            required(record.chrom.as_deref(), "CHROM")?,
            record.pos.ok_or_else(|| missing("POS"))? + 1,
            record.id.as_deref().unwrap_or("."),
            required(record.ref_allele.as_deref(), "REF")?,
            alt_field(record)
        )?;
        if !args.samples.is_empty() {
            write!(out, "\tGT")?;
            for gt in &converted.genotypes {
                write!(out, "\t{gt}")?;
            }
        }
        writeln!(out)?;
    }
    Ok(())
}

fn write_contig_headers<W: Write>(out: &mut W, reference: &FastaReference) -> io::Result<()> {
    for record in reference.index().as_ref() {
        let name = String::from_utf8_lossy(record.name().as_ref());
        writeln!(out, "##contig=<ID={name},length={}>", record.length())?;
    }
    Ok(())
}

fn write_record_contig_headers<W: Write>(
    out: &mut W,
    records: &[ConvertedRecord],
) -> io::Result<()> {
    let mut seen = std::collections::BTreeSet::new();
    for converted in records {
        let Some(chrom) = converted.record.chrom.as_deref() else {
            continue;
        };
        if seen.insert(chrom) {
            writeln!(out, "##contig=<ID={chrom}>")?;
        }
    }
    Ok(())
}

fn alt_field(record: &TsvRecord) -> String {
    if record.alt_alleles.is_empty() {
        ".".to_owned()
    } else {
        record.alt_alleles.join(",")
    }
}

fn parse_trailing_genotypes(
    line: &str,
    tsv: &Tsv,
    record: &TsvRecord,
    samples: &[String],
) -> io::Result<Vec<String>> {
    if samples.is_empty() {
        return Ok(Vec::new());
    }
    let fields: Vec<&str> = line.split_whitespace().collect();
    let offset = tsv.columns().len();
    let available = fields.len().saturating_sub(offset);
    if available < samples.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "expected {} sample genotype columns, found {available}",
                samples.len()
            ),
        ));
    }
    samples
        .iter()
        .enumerate()
        .map(|(i, _)| normalize_gt(fields[offset + i], record))
        .collect()
}

fn populate_from_aa(
    line: &str,
    tsv: &Tsv,
    record: &mut TsvRecord,
    samples: &[String],
    reference: Option<&FastaReference>,
    stats: &mut ConversionStats,
) -> io::Result<Option<Vec<String>>> {
    let reference = reference.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "--tsv2vcf requires the --fasta-ref option when AA is used",
        )
    })?;
    let chrom = required(record.chrom.as_deref(), "CHROM")?;
    let pos0 = record.pos.ok_or_else(|| missing("POS"))?;
    let pos1 = pos0 + 1;
    let region = format!("{chrom}:{pos1}-{pos1}");
    let ref_base = reference.fetch_region(&region)?;
    let ref_allele = String::from_utf8_lossy(&ref_base)
        .to_ascii_uppercase()
        .chars()
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "empty reference fetch"))?
        .to_string();

    let fields: Vec<&str> = line.split_whitespace().collect();
    let aa_offset = tsv
        .columns()
        .iter()
        .position(|column| {
            column
                .as_deref()
                .is_some_and(|name| name.eq_ignore_ascii_case("AA"))
        })
        .ok_or_else(|| missing("AA"))?;
    let sample_count = samples.len();
    let raw_values: Vec<&str> = if sample_count == 0 {
        Vec::new()
    } else {
        let available = fields.len().saturating_sub(aa_offset);
        if available < sample_count {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("expected {sample_count} sample genotype columns, found {available}"),
            ));
        }
        fields[aa_offset..aa_offset + sample_count].to_vec()
    };

    let mut alleles = vec![ref_allele.clone()];
    let mut genotypes = Vec::with_capacity(raw_values.len());
    for raw in raw_values {
        match aa_value_to_gt(raw, &ref_allele, &mut alleles)? {
            Some((gt, class)) => {
                genotypes.push(gt);
                if let Some(class) = class {
                    stats.observe_gt(class);
                }
            }
            None => return Ok(None),
        }
    }

    record.ref_allele = Some(ref_allele);
    record.alt_alleles = alleles.into_iter().skip(1).collect();
    Ok(Some(genotypes))
}

fn aa_value_to_gt(
    raw: &str,
    ref_allele: &str,
    alleles: &mut Vec<String>,
) -> io::Result<Option<(String, Option<GtClass>)>> {
    if raw.eq_ignore_ascii_case("I") || raw.eq_ignore_ascii_case("D") {
        return Ok(None);
    }
    if raw == "." || raw == "-" || raw == "--" {
        return Ok(Some(("./.".to_owned(), Some(GtClass::Missing))));
    }
    if raw.len() > 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "expected AA genotype to have one or two characters",
        ));
    }

    let mut indexes = Vec::new();
    let mut allele_codes = Vec::new();
    for allele in raw.chars() {
        let allele = allele.to_ascii_uppercase();
        if allele == '.' || allele == '-' {
            indexes.push(".".to_owned());
            allele_codes.push(None);
            continue;
        }
        if !matches!(allele, 'A' | 'C' | 'G' | 'T' | 'N') {
            return Err(invalid_gt(raw));
        }
        let allele_s = allele.to_string();
        let idx = if allele_s.eq_ignore_ascii_case(ref_allele) {
            0
        } else if let Some(idx) = alleles.iter().position(|known| known == &allele_s) {
            idx
        } else {
            alleles.push(allele_s);
            alleles.len() - 1
        };
        indexes.push(idx.to_string());
        allele_codes.push(Some(idx));
    }
    let class = classify_gt_indexes(&allele_codes);
    let gt = if indexes.len() == 1 {
        indexes[0].clone()
    } else {
        indexes.join("/")
    };
    Ok(Some((gt, class)))
}

#[derive(Debug, Clone, Copy)]
enum GtClass {
    Missing,
    HomRef,
    HetRefAlt,
    HomAlt,
    HetAlt,
}

impl ConversionStats {
    fn observe_gt(&mut self, class: GtClass) {
        match class {
            GtClass::Missing => self.missing += 1,
            GtClass::HomRef => self.hom_rr += 1,
            GtClass::HetRefAlt => self.het_ra += 1,
            GtClass::HomAlt => self.hom_aa += 1,
            GtClass::HetAlt => self.het_aa += 1,
        }
    }
}

fn classify_gt_indexes(indexes: &[Option<usize>]) -> Option<GtClass> {
    if indexes.is_empty() || indexes.iter().any(Option::is_none) {
        return Some(GtClass::Missing);
    }
    let a0 = indexes[0]?;
    let a1 = indexes.get(1).copied().flatten().unwrap_or(a0);
    if a0 == 0 && a1 == 0 {
        Some(GtClass::HomRef)
    } else if a0 == 0 || a1 == 0 {
        Some(GtClass::HetRefAlt)
    } else if a0 == a1 {
        Some(GtClass::HomAlt)
    } else {
        Some(GtClass::HetAlt)
    }
}

fn normalize_gt(raw: &str, record: &TsvRecord) -> io::Result<String> {
    if raw == "." || raw == "./." || raw == ".|." {
        return Ok(raw.to_owned());
    }
    if raw.contains('/') || raw.contains('|') {
        return validate_vcf_gt(raw);
    }
    if raw.chars().all(|c| c.is_ascii_digit()) {
        return match raw.len() {
            1 => validate_vcf_gt(raw),
            2 => Ok(format!("{}/{}", &raw[0..1], &raw[1..2])),
            _ => Err(invalid_gt(raw)),
        };
    }
    allele_string_to_gt(raw, record)
}

fn validate_vcf_gt(raw: &str) -> io::Result<String> {
    let valid = raw
        .split(['/', '|'])
        .all(|allele| allele == "." || allele.parse::<usize>().is_ok());
    if valid {
        Ok(raw.to_owned())
    } else {
        Err(invalid_gt(raw))
    }
}

fn allele_string_to_gt(raw: &str, record: &TsvRecord) -> io::Result<String> {
    let alleles: Vec<char> = raw.chars().collect();
    if alleles.is_empty() || alleles.len() > 2 {
        return Err(invalid_gt(raw));
    }
    let mut out = Vec::with_capacity(alleles.len());
    for allele in alleles {
        if allele == '.' || allele == '-' || allele.eq_ignore_ascii_case(&'N') {
            out.push(".".to_owned());
            continue;
        }
        let allele_s = allele.to_string();
        if record
            .ref_allele
            .as_deref()
            .is_some_and(|reference| reference.eq_ignore_ascii_case(&allele_s))
        {
            out.push("0".to_owned());
            continue;
        }
        let Some((alt_i, _)) = record
            .alt_alleles
            .iter()
            .enumerate()
            .find(|(_, alt)| alt.eq_ignore_ascii_case(&allele_s))
        else {
            return Err(invalid_gt(raw));
        };
        out.push((alt_i + 1).to_string());
    }
    Ok(out.join("/"))
}

fn invalid_gt(raw: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("invalid sample genotype '{raw}'"),
    )
}

fn required<'a>(value: Option<&'a str>, name: &str) -> io::Result<&'a str> {
    value.ok_or_else(|| missing(name))
}

fn missing(name: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("missing {name} column"))
}

fn write_output(args: &Args, buffer: &[u8]) -> io::Result<()> {
    match args.output_kind {
        OutputKind::VcfText => match &args.output {
            Some(path) => fs::write(path, buffer),
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
        OutputKind::Bcf => {
            use htslib_rs::vcf::variant::io::Write as _;

            let mut reader = htslib_rs::vcf::io::Reader::new(BufReader::new(buffer));
            let header = reader.read_header()?;
            if let (Some(path), Some(thread_count)) = (&args.output, args.thread_count) {
                let file = File::create(path)?;
                let bgzf =
                    htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
                let mut writer = htslib_rs::bcf::io::Writer::from(bgzf);
                writer.write_variant_header(&header)?;
                for result in reader.records() {
                    let record = result?;
                    writer.write_variant_record(&header, &record)?;
                }
                let mut bgzf = writer.into_inner();
                let _file = bgzf.finish()?;
                return Ok(());
            }
            let dst: Box<dyn Write> = match &args.output {
                Some(path) => Box::new(File::create(path)?),
                None => Box::new(io::stdout().lock()),
            };
            let mut writer = htslib_rs::bcf::io::Writer::new(dst);
            writer.write_variant_header(&header)?;
            for result in reader.records() {
                let record = result?;
                writer.write_variant_record(&header, &record)?;
            }
            writer.try_finish()
        }
    }
}

fn parse_output_kind(raw: &str) -> Result<OutputKind, ParseOutcome> {
    OutputKind::parse(raw)
        .ok_or_else(|| ParseOutcome::Error(format!("The output type \"{raw}\" not recognised")))
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

fn parse_sample_list(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(str::trim)
        .filter(|sample| !sample.is_empty())
        .map(str::to_owned)
        .collect()
}

fn read_sample_file(path: &str) -> Result<Vec<String>, ParseOutcome> {
    let text = fs::read_to_string(path)
        .map_err(|e| ParseOutcome::Error(format!("failed to read samples file '{path}': {e}")))?;
    Ok(text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .filter_map(|line| line.split_whitespace().next())
        .map(str::to_owned)
        .collect())
}

fn write_index(path: &Path, output_kind: OutputKind, index_format: i32) -> io::Result<()> {
    let path_str = path
        .to_str()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "non-UTF-8 output path"))?;
    let output_format = match output_kind {
        OutputKind::VcfText | OutputKind::VcfGz => VariantOutputFormat::Vcf,
        OutputKind::Bcf => VariantOutputFormat::Bcf,
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
