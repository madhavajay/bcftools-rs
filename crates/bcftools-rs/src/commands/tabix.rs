//! Port of `bcftools tabix` (upstream `tabix.c`).
//!
//! This command is intentionally hidden from top-level help upstream but is
//! kept for test compatibility. It supports the common preset-driven paths:
//! building TBI/CSI indexes for BGZF-compressed BED/GFF/SAM/VCF and querying
//! one or more regions from an existing associated index.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, BufRead as _, Write as _};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::core::Region;
use htslib_rs::index_compat::{read_csi, read_tbi, write_csi, write_tbi};
use htslib_rs::tabix_compat::{
    TextFormat, build_csi_from_bgzf_path_with_min_shift, build_tbi_from_bgzf_path,
    query_csi_records_from_path, query_records_from_path,
};

const USAGE: &str = "\n\
Usage: bcftools tabix [options] <in.gz> [reg1 [...]]\n\
\n\
Options: -p, --preset STR   preset: gff, bed, sam or vcf [gff]\n\
         -s INT    column number for sequence names (suppressed by -p) [1]\n\
         -b INT    column number for region start [4]\n\
         -e INT    column number for region end (if no end, set INT to -b) [5]\n\
         -0        specify coordinates are zero-based\n\
         -S INT    skip first INT lines [0]\n\
         -c CHAR   skip lines starting with CHAR [null]\n\
         -a, --all          print all records\n\
         -f, --force        force to overwrite existing index\n\
         -C, --csi          generate CSI index\n\
         -m, --min-shift INT    set the minimal interval size to 1<<INT; 0 for the old tabix index [0]\n\
\n";

#[derive(Debug)]
struct Args {
    input: PathBuf,
    regions: Vec<String>,
    format: TextFormat,
    min_shift: i32,
    force: bool,
    all: bool,
}

#[derive(Debug)]
enum ParseOutcome {
    Usage,
    Error(String),
}

/// Subcommand entry point. `argv[0]` is `"tabix"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{e}");
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{message}");
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut all = false;
    let mut force = false;
    let mut min_shift = -1;
    let mut format = TextFormat::Gff;
    let mut detect = true;
    let mut input: Option<PathBuf> = None;
    let mut regions = Vec::new();

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-0" => {}
            "-a" | "--all" => all = true,
            "-f" | "--force" => force = true,
            "-C" | "--csi" => min_shift = 14,
            "-t" | "--tbi" => min_shift = 0,
            "-m" => {
                min_shift = next_string(&mut iter, "-m")?
                    .parse()
                    .map_err(|_| ParseOutcome::Error("Could not parse argument: -m".into()))?;
            }
            "--min-shift" => {
                min_shift = next_string(&mut iter, "--min-shift")?
                    .parse()
                    .map_err(|_| {
                        ParseOutcome::Error("Could not parse argument: --min-shift".into())
                    })?;
            }
            "-p" => {
                format = parse_preset(&next_string(&mut iter, "-p")?)?;
                detect = false;
            }
            "--preset" => {
                format = parse_preset(&next_string(&mut iter, "--preset")?)?;
                detect = false;
            }
            "-s" | "-b" | "-e" | "-S" | "-c" | "--sequence" | "--begin" | "--end"
            | "--skip-lines" | "--comment" => {
                let _ = next_string(&mut iter, raw.as_ref())?;
            }
            _ if raw.starts_with("--preset=") => {
                format = parse_preset(value_after_equals(&raw))?;
                detect = false;
            }
            _ if raw.starts_with("--min-shift=") => {
                min_shift = value_after_equals(&raw).parse().map_err(|_| {
                    ParseOutcome::Error("Could not parse argument: --min-shift".into())
                })?;
            }
            _ if raw.starts_with("-m") && raw.len() > 2 => {
                min_shift = raw[2..]
                    .parse()
                    .map_err(|_| ParseOutcome::Error("Could not parse argument: -m".into()))?;
            }
            _ if raw.starts_with("-p") && raw.len() > 2 => {
                format = parse_preset(&raw[2..])?;
                detect = false;
            }
            _ if raw.starts_with('-') => return Err(ParseOutcome::Usage),
            _ => {
                if input.is_none() {
                    input = Some(PathBuf::from(arg));
                } else {
                    regions.push(raw.into_owned());
                }
            }
        }
    }

    let input = input.ok_or(ParseOutcome::Usage)?;
    if detect {
        format = detect_format(&input).unwrap_or(format);
    }

    Ok(Args {
        input,
        regions,
        format,
        min_shift,
        force,
        all,
    })
}

fn run(args: &Args) -> io::Result<()> {
    if args.all {
        return print_all(&args.input);
    }

    if args.regions.is_empty() {
        return build_index(args);
    }

    query_regions(args)
}

fn build_index(args: &Args) -> io::Result<()> {
    let index_path = index_path(&args.input, args.min_shift);
    if !args.force && index_path.exists() {
        return Err(io::Error::new(
            io::ErrorKind::AlreadyExists,
            "[E::main_tabix] the index file exists; use option '-f' to overwrite",
        ));
    }

    if args.min_shift <= 0 {
        let index = build_tbi_from_bgzf_path(&args.input, args.format).map_err(build_error)?;
        write_tbi(index_path, &index)
    } else {
        let min_shift = u8::try_from(args.min_shift)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "invalid min_shift"))?;
        let index = build_csi_from_bgzf_path_with_min_shift(&args.input, args.format, min_shift)
            .map_err(build_error)?;
        write_csi(index_path, &index)
    }
}

fn query_regions(args: &Args) -> io::Result<()> {
    let mut out = io::stdout().lock();
    let tbi_path = associated_path(&args.input, "tbi");
    if tbi_path.exists() {
        let index = read_tbi(&tbi_path)?;
        for raw in &args.regions {
            let region = parse_region(raw)?;
            for record in query_records_from_path(&args.input, index.clone(), &region)? {
                writeln!(out, "{record}")?;
            }
        }
        return Ok(());
    }

    let csi_path = associated_path(&args.input, "csi");
    let index = read_csi(&csi_path)?;
    for raw in &args.regions {
        let region = parse_region(raw)?;
        for record in query_csi_records_from_path(&args.input, index.clone(), &region)? {
            writeln!(out, "{record}")?;
        }
    }
    Ok(())
}

fn print_all(path: &Path) -> io::Result<()> {
    let file = File::open(path)?;
    let mut reader = io::BufReader::new(htslib_rs::bgzf::io::Reader::new(file));
    let mut out = io::stdout().lock();
    let mut line = String::new();
    while reader.read_line(&mut line)? != 0 {
        out.write_all(line.as_bytes())?;
        line.clear();
    }
    Ok(())
}

fn build_error(_: io::Error) -> io::Error {
    io::Error::other(
        "tbx_index_build failed: Is the file bgzip-compressed? Was wrong -p [type] option used?",
    )
}

fn parse_region(raw: &str) -> io::Result<Region> {
    raw.parse()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))
}

fn index_path(input: &Path, min_shift: i32) -> PathBuf {
    associated_path(input, if min_shift <= 0 { "tbi" } else { "csi" })
}

fn associated_path(input: &Path, ext: &str) -> PathBuf {
    let mut path = input.as_os_str().to_owned();
    path.push(".");
    path.push(ext);
    PathBuf::from(path)
}

fn detect_format(path: &Path) -> Option<TextFormat> {
    let name = path.file_name()?.to_string_lossy();
    if ends_with_ci(&name, ".gff.gz") {
        Some(TextFormat::Gff)
    } else if ends_with_ci(&name, ".bed.gz") {
        Some(TextFormat::Bed)
    } else if ends_with_ci(&name, ".sam.gz") {
        Some(TextFormat::Sam)
    } else if ends_with_ci(&name, ".vcf.gz") {
        Some(TextFormat::Vcf)
    } else {
        None
    }
}

fn ends_with_ci(s: &str, suffix: &str) -> bool {
    let n = suffix.len();
    s.len() >= n && s[s.len() - n..].eq_ignore_ascii_case(suffix)
}

fn value_after_equals(raw: &str) -> &str {
    raw.split_once('=').map(|(_, value)| value).unwrap_or("")
}

fn parse_preset(raw: &str) -> Result<TextFormat, ParseOutcome> {
    match raw {
        "gff" => Ok(TextFormat::Gff),
        "bed" => Ok(TextFormat::Bed),
        "sam" => Ok(TextFormat::Sam),
        "vcf" => Ok(TextFormat::Vcf),
        _ => Err(ParseOutcome::Error(format!(
            "The type '{raw}' not recognised"
        ))),
    }
}

fn next_string<'a, I>(iter: &mut I, name: &str) -> Result<String, ParseOutcome>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .map(|s| s.to_string_lossy().into_owned())
        .ok_or_else(|| ParseOutcome::Error(format!("Missing argument for {name}")))
}
