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
//! Filtering is NOT yet implemented and yields an explicit error if requested.
//! Positional region arguments support the common `CHROM` and `CHROM:START-END`
//! forms by streaming and filtering records.

use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufRead as _, BufReader, Cursor, Read as _, Write};
use std::num::NonZero;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{SystemTime, UNIX_EPOCH};

use htslib_rs::core::Position;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};
use htslib_rs::variant_io_compat::RegionOverlap;
use htslib_rs::vcf::variant::io::Write as _;

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};
use crate::header_version::{build_lines, command_time};
use crate::io::parse_overlap_option;

const USAGE: &str = "\n\
About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.\n\
Usage:   bcftools view [OPTIONS] <in.vcf.gz>|<in.bcf> [REGION...]\n\
\n\
Output options:\n\
    -G, --drop-genotypes              drop individual genotype information\n\
    -f, --apply-filters LIST          require at least one listed FILTER string\n\
    -g, --genotype [^]hom|het|miss    require or exclude genotype class\n\
    -h, --header-only                 print only the header in VCF output\n\
    -H, --no-header                   suppress the header in VCF output\n\
    -l, --compression-level INT       compression level: 0 uncompressed, 1 best speed, 9 best compression [-1]\n\
    -k, --known                       select known sites only (ID is not '.')\n\
    -m, --min-alleles INT             minimum number of alleles listed in REF and ALT\n\
    -M, --max-alleles INT             maximum number of alleles listed in REF and ALT\n\
    -n, --novel                       select novel sites only (ID is '.')\n\
        --no-version                  do not append version and command line to the header\n\
    -p, --phased                      select sites where all samples are phased\n\
    -P, --exclude-phased              exclude sites where all samples are phased\n\
    -o, --output FILE                 output file name [stdout]\n\
    -O, --output-type u|b|v|z[0-9]    u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n\
    -r, --regions REG                 restrict to comma-separated regions\n\
    -R, --regions-file FILE           restrict to regions listed in a file\n\
        --regions-overlap 0|1|2       region overlap mode: 0=POS, 1=record, 2=variant [0]\n\
    -s, --samples LIST                comma-separated sample list, optionally prefixed with ^\n\
    -S, --samples-file FILE           file of samples, optionally prefixed with ^\n\
    -t, --targets REG                 restrict to comma-separated targets\n\
    -T, --targets-file FILE           restrict to targets listed in a file\n\
        --targets-overlap 0|1|2       target overlap mode: 0=POS, 1=record, 2=variant [0]\n\
    -u, --uncalled                    select sites without a called genotype\n\
    -U, --exclude-uncalled            exclude sites without a called genotype\n\
    -v, --types LIST                  select variant types: snps,indels,mnps,other,bnd,overlap,ref\n\
    -V, --exclude-types LIST          exclude variant types\n\
\n";

const OPT_NO_VERSION: i32 = 200;
const OPT_THREADS: i32 = 9;
const OPT_TARGETS_OVERLAP: i32 = 201;
const OPT_REGIONS_OVERLAP: i32 = 202;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    BcfUncompressed,
    BcfCompressed,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Region {
    contig: String,
    start: Option<usize>,
    end: Option<usize>,
}

struct RunOptions<'a> {
    output_kind: OutputKind,
    output_file: Option<&'a str>,
    header_only: bool,
    no_header: bool,
    no_version: bool,
    regions: &'a [Region],
    regions_overlap: RegionOverlap,
    targets: &'a [Region],
    targets_exclude: bool,
    targets_overlap: RegionOverlap,
    apply_filters: Option<Vec<String>>,
    type_filter: Option<TypeFilter>,
    type_filter_exclude: bool,
    min_alleles: Option<usize>,
    max_alleles: Option<usize>,
    phased_filter: Option<bool>,
    known_filter: Option<bool>,
    uncalled_filter: Option<bool>,
    genotype_filter: Option<GenotypeFilter>,
    thread_count: Option<NonZero<usize>>,
    sample_list: Option<&'a str>,
    sample_list_is_file: bool,
    drop_genotypes: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TypeFilter {
    mask: VariantType,
    include_ref: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct GenotypeFilter {
    class: GenotypeClass,
    exclude: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GenotypeClass {
    Hom,
    Het,
    Missing,
}

#[derive(Clone, Copy)]
struct RecordFilters<'a> {
    regions: &'a [Region],
    targets: &'a [Region],
    targets_exclude: bool,
}

impl Region {
    fn parse(raw: &str) -> io::Result<Self> {
        let (contig, interval) = raw.split_once(':').unwrap_or((raw, ""));
        if contig.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Could not parse region \"{raw}\""),
            ));
        }

        if interval.is_empty() {
            return Ok(Self {
                contig: contig.to_string(),
                start: None,
                end: None,
            });
        }

        let (start, end) = interval
            .split_once('-')
            .map(|(start, end)| (start, Some(end)))
            .unwrap_or((interval, None));
        let start = if start.is_empty() {
            None
        } else {
            Some(parse_region_pos(start, raw)?)
        };
        let end = match end {
            Some("") => None,
            Some(end) => Some(parse_region_pos(end, raw)?),
            None => start,
        };

        Ok(Self {
            contig: contig.to_string(),
            start,
            end,
        })
    }

    fn contains(&self, contig: &str, pos: Position) -> bool {
        if self.contig != contig {
            return false;
        }
        let pos = usize::from(pos);
        self.start.is_none_or(|start| pos >= start) && self.end.is_none_or(|end| pos <= end)
    }
}

fn parse_region_pos(s: &str, raw_region: &str) -> io::Result<usize> {
    s.replace(',', "").parse::<usize>().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Could not parse region \"{raw_region}\""),
        )
    })
}

fn extend_regions_from_list(regions: &mut Vec<Region>, raw: &str) -> io::Result<()> {
    for item in raw.split(',') {
        let item = item.trim();
        if !item.is_empty() {
            regions.push(Region::parse(item)?);
        }
    }
    Ok(())
}

fn strip_exclusion_prefix<'a>(raw: &'a str, label: &str) -> io::Result<(bool, &'a str)> {
    let trimmed = raw.trim();
    if let Some(body) = trimmed.strip_prefix('^') {
        if body.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("{label} exclusion is missing a value"),
            ));
        }
        Ok((true, body))
    } else {
        Ok((false, trimmed))
    }
}

fn set_target_exclusion_mode(mode: &mut Option<bool>, exclude: bool) -> io::Result<()> {
    if mode.is_some_and(|existing| existing != exclude) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "cannot mix target inclusion and exclusion lists",
        ));
    }
    *mode = Some(exclude);
    Ok(())
}

fn extend_regions_from_file(regions: &mut Vec<Region>, path: &Path) -> io::Result<()> {
    use crate::regidx::{Parser, parse_bed, parse_tab, parser_for_path};

    let parser = parser_for_path(path);
    let reader: Box<dyn io::BufRead> = if path
        .as_os_str()
        .to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".gz")
    {
        Box::new(io::BufReader::new(flate2::read::MultiGzDecoder::new(
            File::open(path)?,
        )))
    } else {
        Box::new(io::BufReader::new(File::open(path)?))
    };

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let record = match parser {
            Parser::Bed => parse_bed(line),
            _ => parse_tab(line),
        }
        .map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("failed to parse region line: {e:?}"),
            )
        })?;
        if let Some(record) = record {
            regions.push(Region {
                contig: record.seq,
                start: Some((record.start + 1).try_into().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "region start is out of range")
                })?),
                end: Some((record.end + 1).try_into().map_err(|_| {
                    io::Error::new(io::ErrorKind::InvalidData, "region end is out of range")
                })?),
            });
        }
    }

    Ok(())
}

fn extend_targets_from_file(targets: &mut Vec<Region>, path: &Path) -> io::Result<()> {
    use crate::regidx::{Parser, parser_for_path};

    if parser_for_path(path) == Parser::Bed {
        return extend_regions_from_file(targets, path);
    }

    let reader: Box<dyn io::BufRead> = if path
        .as_os_str()
        .to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(".gz")
    {
        Box::new(io::BufReader::new(flate2::read::MultiGzDecoder::new(
            File::open(path)?,
        )))
    } else {
        Box::new(io::BufReader::new(File::open(path)?))
    };

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split_whitespace();
        let Some(contig) = fields.next() else {
            continue;
        };
        let Some(pos) = fields.next() else {
            continue;
        };
        let pos = parse_region_pos(pos, line)?;
        targets.push(Region {
            contig: contig.to_string(),
            start: Some(pos),
            end: Some(pos),
        });
    }

    Ok(())
}

fn parse_threads(raw: &str) -> Result<Option<NonZero<usize>>, std::num::ParseIntError> {
    let n = raw.parse::<usize>()?;
    Ok(NonZero::new(n))
}

fn parse_positive_usize(raw: &str, option: &str) -> io::Result<usize> {
    raw.parse::<usize>().ok().filter(|n| *n > 0).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Could not parse argument: {option} {raw}"),
        )
    })
}

fn parse_variant_type_filter(raw: &str) -> io::Result<TypeFilter> {
    let mut mask = VariantType::REF;
    let mut include_ref = false;
    for item in raw.split(',') {
        let item = item.trim().to_ascii_lowercase();
        let variant_type = match item.as_str() {
            "" => continue,
            "snp" | "snps" => VariantType::SNP,
            "mnp" | "mnps" => VariantType::MNP,
            "indel" | "indels" => VariantType::INDEL,
            "other" => VariantType::OTHER,
            "bnd" => VariantType::BND,
            "overlap" => VariantType::OVERLAP,
            "ins" | "insertion" | "insertions" => VariantType::INS,
            "del" | "deletion" | "deletions" => VariantType::DEL,
            "ref" => {
                include_ref = true;
                continue;
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("The type \"{item}\" not recognised"),
                ));
            }
        };
        mask |= variant_type;
    }

    if mask.bits() == 0 && !include_ref {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing variant type",
        ));
    }

    Ok(TypeFilter { mask, include_ref })
}

fn parse_genotype_filter(raw: &str) -> io::Result<GenotypeFilter> {
    let (exclude, body) = raw
        .strip_prefix('^')
        .map_or((false, raw), |body| (true, body));
    let class = match body.to_ascii_lowercase().as_str() {
        "hom" => GenotypeClass::Hom,
        "het" => GenotypeClass::Het,
        "miss" => GenotypeClass::Missing,
        _ => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"{raw}\"."
                ),
            ));
        }
    };
    Ok(GenotypeFilter { class, exclude })
}

fn parse_apply_filters(raw: &str) -> io::Result<Vec<String>> {
    let filters = raw
        .split(',')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    if filters.is_empty() {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "missing FILTER value",
        ))
    } else {
        Ok(filters)
    }
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
        LongOpt::new("apply-filters", HasArg::Required, b'f' as i32),
        LongOpt::new("genotype", HasArg::Required, b'g' as i32),
        LongOpt::new("known", HasArg::None, b'k' as i32),
        LongOpt::new("min-alleles", HasArg::Required, b'm' as i32),
        LongOpt::new("max-alleles", HasArg::Required, b'M' as i32),
        LongOpt::new("novel", HasArg::None, b'n' as i32),
        LongOpt::new("phased", HasArg::None, b'p' as i32),
        LongOpt::new("exclude-phased", HasArg::None, b'P' as i32),
        LongOpt::new("header-only", HasArg::None, b'h' as i32),
        LongOpt::new("no-header", HasArg::None, b'H' as i32),
        LongOpt::new("regions", HasArg::Required, b'r' as i32),
        LongOpt::new("regions-file", HasArg::Required, b'R' as i32),
        LongOpt::new("regions-overlap", HasArg::Required, OPT_REGIONS_OVERLAP),
        LongOpt::new("samples", HasArg::Required, b's' as i32),
        LongOpt::new("samples-file", HasArg::Required, b'S' as i32),
        LongOpt::new("targets", HasArg::Required, b't' as i32),
        LongOpt::new("targets-file", HasArg::Required, b'T' as i32),
        LongOpt::new("uncalled", HasArg::None, b'u' as i32),
        LongOpt::new("exclude-uncalled", HasArg::None, b'U' as i32),
        LongOpt::new("drop-genotypes", HasArg::None, b'G' as i32),
        LongOpt::new("types", HasArg::Required, b'v' as i32),
        LongOpt::new("exclude-types", HasArg::Required, b'V' as i32),
        LongOpt::new("no-version", HasArg::None, OPT_NO_VERSION),
        LongOpt::new("threads", HasArg::Required, OPT_THREADS),
        LongOpt::new("targets-overlap", HasArg::Required, OPT_TARGETS_OVERLAP),
    ];

    let mut output_kind = OutputKind::VcfText;
    let mut compression_level: Option<u32> = None;
    let mut output_file: Option<String> = None;
    let mut header_only = false;
    let mut no_header = false;
    let mut no_version = false;
    let mut thread_count = None;
    let mut regions_overlap = RegionOverlap::Pos;
    let mut targets_overlap = RegionOverlap::Pos;
    let mut sample_list: Option<String> = None;
    let mut sample_list_is_file = false;
    let mut drop_genotypes = false;
    let mut apply_filters = None;
    let mut type_filter = None;
    let mut type_filter_exclude = false;
    let mut min_alleles = None;
    let mut max_alleles = None;
    let mut phased_filter = None;
    let mut known_filter = None;
    let mut uncalled_filter = None;
    let mut genotype_filter = None;
    let mut region_specs = Vec::new();
    let mut region_files = Vec::new();
    let mut target_specs = Vec::new();
    let mut target_files = Vec::new();

    let mut g = Getopt::new("o:O:l:f:g:km:M:npPhHr:R:s:S:t:T:uUGv:V:", &long_opts, argv);
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
                v if v == b'f' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_apply_filters(raw) {
                        Ok(filters) => apply_filters = Some(filters),
                        Err(e) => {
                            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == b'm' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_positive_usize(raw, "--min-alleles") {
                        Ok(value) => min_alleles = Some(value),
                        Err(e) => {
                            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == b'M' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_positive_usize(raw, "--max-alleles") {
                        Ok(value) => max_alleles = Some(value),
                        Err(e) => {
                            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == b'g' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_genotype_filter(raw) {
                        Ok(filter) => genotype_filter = Some(filter),
                        Err(e) => {
                            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == b'k' as i32 || v == b'n' as i32 => {
                    let include_known = v == b'k' as i32;
                    if known_filter.is_some_and(|existing| existing != include_known) {
                        eprintln!(
                            "{}",
                            fmt_etag("main_vcfview", "Only one of -k or -n can be given.")
                        );
                        return ExitCode::FAILURE;
                    }
                    known_filter = Some(include_known);
                }
                v if v == b'p' as i32 || v == b'P' as i32 => {
                    let include = v == b'p' as i32;
                    if phased_filter.is_some_and(|existing| existing != include) {
                        eprintln!(
                            "{}",
                            fmt_etag("main_vcfview", "Only one of -p or -P can be given.")
                        );
                        return ExitCode::FAILURE;
                    }
                    phased_filter = Some(include);
                }
                v if v == b'u' as i32 || v == b'U' as i32 => {
                    let include = v == b'u' as i32;
                    if uncalled_filter.is_some_and(|existing| existing != include) {
                        eprintln!(
                            "{}",
                            fmt_etag("main_vcfview", "Only one of -u or -U can be given.")
                        );
                        return ExitCode::FAILURE;
                    }
                    uncalled_filter = Some(include);
                }
                v if v == b'h' as i32 => header_only = true,
                v if v == b'H' as i32 => no_header = true,
                v if v == b'r' as i32 => {
                    if let Some(value) = m.value {
                        region_specs.push(value);
                    }
                }
                v if v == b'R' as i32 => {
                    if let Some(value) = m.value {
                        region_files.push(value);
                    }
                }
                v if v == b's' as i32 => {
                    sample_list = m.value;
                    sample_list_is_file = false;
                }
                v if v == b'S' as i32 => {
                    sample_list = m.value;
                    sample_list_is_file = true;
                }
                v if v == b't' as i32 => {
                    if let Some(value) = m.value {
                        target_specs.push(value);
                    }
                }
                v if v == b'T' as i32 => {
                    if let Some(value) = m.value {
                        target_files.push(value);
                    }
                }
                v if v == b'G' as i32 => drop_genotypes = true,
                v if v == b'v' as i32 || v == b'V' as i32 => {
                    let exclude = v == b'V' as i32;
                    if type_filter.is_some() && type_filter_exclude != exclude {
                        eprintln!(
                            "{}",
                            fmt_etag(
                                "main_vcfview",
                                "cannot mix type inclusion and exclusion lists"
                            )
                        );
                        return ExitCode::FAILURE;
                    }
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_variant_type_filter(raw) {
                        Ok(filter) => {
                            type_filter = Some(filter);
                            type_filter_exclude = exclude;
                        }
                        Err(e) => {
                            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == OPT_NO_VERSION => no_version = true,
                v if v == OPT_THREADS => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_threads(raw) {
                        Ok(n) => thread_count = n,
                        Err(_) => {
                            eprintln!(
                                "{}",
                                fmt_etag(
                                    "main_vcfview",
                                    &format!("Could not parse argument: --threads {raw}")
                                )
                            );
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == OPT_REGIONS_OVERLAP => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_overlap_option(raw)
                        .and_then(|mode| RegionOverlap::try_from(mode).ok())
                    {
                        Some(mode) => regions_overlap = mode,
                        None => {
                            eprintln!(
                                "{}",
                                fmt_etag(
                                    "main_vcfview",
                                    &format!("The --regions-overlap mode \"{raw}\" not recognised")
                                )
                            );
                            return ExitCode::FAILURE;
                        }
                    }
                }
                v if v == OPT_TARGETS_OVERLAP => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match parse_overlap_option(raw)
                        .and_then(|mode| RegionOverlap::try_from(mode).ok())
                    {
                        Some(mode) => targets_overlap = mode,
                        None => {
                            eprintln!(
                                "{}",
                                fmt_etag(
                                    "main_vcfview",
                                    &format!("The --targets-overlap mode \"{raw}\" not recognised")
                                )
                            );
                            return ExitCode::FAILURE;
                        }
                    }
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

    let positional = g.rest();
    let fname = positional
        .first()
        .cloned()
        .unwrap_or_else(|| OsString::from("-"));
    let mut regions = Vec::new();
    let parsed_regions = region_specs
        .iter()
        .try_for_each(|raw| extend_regions_from_list(&mut regions, raw));
    if let Err(e) = parsed_regions {
        eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
        return ExitCode::FAILURE;
    }
    for file in &region_files {
        match extend_regions_from_file(&mut regions, Path::new(file)) {
            Ok(()) => {}
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                return ExitCode::FAILURE;
            }
        }
    }
    let parsed_positionals = positional
        .iter()
        .skip(1)
        .try_for_each(|raw| extend_regions_from_list(&mut regions, &raw.to_string_lossy()));
    if let Err(e) = parsed_positionals {
        eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
        return ExitCode::FAILURE;
    }
    let mut targets = Vec::new();
    let mut target_exclusion_mode = None;
    for raw in &target_specs {
        let (exclude, body) = match strip_exclusion_prefix(raw, "target") {
            Ok(value) => value,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                return ExitCode::FAILURE;
            }
        };
        if let Err(e) = set_target_exclusion_mode(&mut target_exclusion_mode, exclude)
            .and_then(|()| extend_regions_from_list(&mut targets, body))
        {
            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
            return ExitCode::FAILURE;
        }
    }
    for file in &target_files {
        let (exclude, path) = match strip_exclusion_prefix(file, "targets-file") {
            Ok(value) => value,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
                return ExitCode::FAILURE;
            }
        };
        if let Err(e) = set_target_exclusion_mode(&mut target_exclusion_mode, exclude)
            .and_then(|()| extend_targets_from_file(&mut targets, Path::new(path)))
        {
            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
            return ExitCode::FAILURE;
        }
    }
    let targets_exclude = target_exclusion_mode.unwrap_or(false);

    let path = Path::new(&fname);
    let _ = compression_level; // consumed by future writers

    let options = RunOptions {
        output_kind,
        output_file: output_file.as_deref(),
        header_only,
        no_header,
        no_version,
        regions: &regions,
        regions_overlap,
        targets: &targets,
        targets_exclude,
        targets_overlap,
        apply_filters,
        type_filter,
        type_filter_exclude,
        min_alleles,
        max_alleles,
        phased_filter,
        known_filter,
        uncalled_filter,
        genotype_filter,
        thread_count,
        sample_list: sample_list.as_deref(),
        sample_list_is_file,
        drop_genotypes,
    };

    match run(path, &options, argv) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{}", fmt_etag("main_vcfview", &format!("{e}")));
            ExitCode::FAILURE
        }
    }
}

fn run(path: &Path, options: &RunOptions<'_>, argv: &[OsString]) -> io::Result<()> {
    if path == Path::new("-") {
        let tmp = stdin_tmp_path();
        let mut data = Vec::new();
        io::stdin().lock().read_to_end(&mut data)?;
        fs::write(&tmp, data)?;
        let result = run(&tmp, options, argv);
        let _ = fs::remove_file(&tmp);
        return result;
    }

    let in_fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    let has_line_filters = options.type_filter.is_some()
        || options.apply_filters.is_some()
        || options.min_alleles.is_some()
        || options.max_alleles.is_some()
        || options.phased_filter.is_some()
        || options.known_filter.is_some()
        || options.uncalled_filter.is_some()
        || options.genotype_filter.is_some();
    let has_text_filters =
        !options.regions.is_empty() || !options.targets.is_empty() || has_line_filters;
    if has_line_filters
        && !(options.sample_list.is_some()
            || options.drop_genotypes
            || (options.output_kind == OutputKind::VcfText
                && options.no_version
                && in_fmt.exact != Exact::Bcf))
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "allele/type filtering currently requires text VCF output with --no-version or sample projection",
        ));
    }
    let filters = RecordFilters {
        regions: options.regions,
        targets: options.targets,
        targets_exclude: options.targets_exclude,
    };

    if options.sample_list.is_some() || options.drop_genotypes {
        return run_sample_subset(path, in_fmt, options, argv);
    }

    let mut header = read_header(path, in_fmt)?;

    if !options.no_version {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        let lines = build_lines("bcftools_view", &prog_argv, command_time());
        // Strip the "##" prefix and the "key=" delimiter from each rendered
        // line, then route both into the header via htslib-rs's typed-wrapper
        // helper. Mirrors upstream `bcf_hdr_append_version` which appends
        // "##bcftools_<cmd>Version" and "##bcftools_<cmd>Command" lines.
        for line in [&lines.version_line, &lines.command_line] {
            htslib_rs::header_compat::append_line(&mut header, line)?;
        }
    }

    if options.output_kind == OutputKind::VcfText
        && options.no_version
        && options.regions.is_empty()
        && options.targets.is_empty()
        && in_fmt.exact == Exact::Bcf
    {
        return match options.output_file {
            Some("-") | None => write_bcf_vcf_text_no_version(
                path,
                options.header_only,
                options.no_header,
                io::stdout().lock(),
            ),
            Some(p) => write_bcf_vcf_text_no_version(
                path,
                options.header_only,
                options.no_header,
                File::create(p)?,
            ),
        };
    }

    match options.output_kind {
        OutputKind::VcfText => match options.output_file {
            Some("-") | None
                if options.no_version && has_text_filters && in_fmt.exact != Exact::Bcf =>
            {
                write_vcf_text_filtered_passthrough(path, in_fmt, options, io::stdout().lock())
            }
            Some("-") | None
                if options.no_version
                    && options.regions.is_empty()
                    && options.targets.is_empty()
                    && in_fmt.exact != Exact::Bcf =>
            {
                write_vcf_text_passthrough(
                    path,
                    in_fmt,
                    options.header_only,
                    options.no_header,
                    io::stdout().lock(),
                )
            }
            Some("-") | None => write_vcf(
                path,
                in_fmt,
                &header,
                options.header_only,
                options.no_header,
                filters,
                io::stdout().lock(),
            ),
            Some(p) if options.no_version && has_text_filters && in_fmt.exact != Exact::Bcf => {
                write_vcf_text_filtered_passthrough(path, in_fmt, options, File::create(p)?)
            }
            Some(p)
                if options.no_version
                    && options.regions.is_empty()
                    && options.targets.is_empty()
                    && in_fmt.exact != Exact::Bcf =>
            {
                write_vcf_text_passthrough(
                    path,
                    in_fmt,
                    options.header_only,
                    options.no_header,
                    File::create(p)?,
                )
            }
            Some(p) => write_vcf(
                path,
                in_fmt,
                &header,
                options.header_only,
                options.no_header,
                filters,
                File::create(p)?,
            ),
        },
        OutputKind::VcfGz
            if options.no_version
                && options.regions.is_empty()
                && options.targets.is_empty()
                && options.type_filter.is_none()
                && in_fmt.exact != Exact::Bcf =>
        {
            match options.output_file {
                Some(p) if p != "-" => {
                    let bgzf = htslib_rs::bgzf::io::Writer::new(File::create(p)?);
                    write_vcf_text_passthrough(
                        path,
                        in_fmt,
                        options.header_only,
                        options.no_header,
                        bgzf,
                    )
                }
                _ => {
                    let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
                    write_vcf_text_passthrough(
                        path,
                        in_fmt,
                        options.header_only,
                        options.no_header,
                        bgzf,
                    )
                }
            }
        }
        OutputKind::VcfGz => write_vcf_gz(path, in_fmt, &header, options),
        OutputKind::BcfUncompressed | OutputKind::BcfCompressed => {
            // For uncompressed BCF, upstream uses `wbu`. noodles' bcf writer
            // always wraps in BGZF; an "uncompressed" mode here is treated the
            // same as the compressed path until htslib-rs exposes the raw form.
            match options.output_file {
                Some("-") | None => write_bcf(
                    path,
                    in_fmt,
                    &header,
                    options.header_only,
                    options.no_header,
                    filters,
                    io::stdout().lock(),
                ),
                Some(p) => write_bcf(
                    path,
                    in_fmt,
                    &header,
                    options.header_only,
                    options.no_header,
                    filters,
                    File::create(p)?,
                ),
            }
        }
    }
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-view-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn run_sample_subset(
    path: &Path,
    in_fmt: format::Format,
    options: &RunOptions<'_>,
    argv: &[OsString],
) -> io::Result<()> {
    let version_lines = if options.no_version {
        None
    } else {
        let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
        prog_argv.extend(argv.iter().cloned());
        Some(build_lines("bcftools_view", &prog_argv, command_time()))
    };

    match options.output_kind {
        OutputKind::VcfText => match options.output_file {
            Some("-") | None => write_sample_subset_vcf(
                path,
                in_fmt,
                options,
                version_lines.as_ref(),
                io::stdout().lock(),
            ),
            Some(p) => write_sample_subset_vcf(
                path,
                in_fmt,
                options,
                version_lines.as_ref(),
                File::create(p)?,
            ),
        },
        OutputKind::VcfGz => match (options.output_file, options.thread_count) {
            (Some(p), Some(thread_count)) if p != "-" => {
                let bgzf = htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(
                    thread_count,
                    File::create(p)?,
                );
                write_sample_subset_vcf(path, in_fmt, options, version_lines.as_ref(), bgzf)
            }
            (Some(p), _) if p != "-" => {
                let bgzf = htslib_rs::bgzf::io::Writer::new(File::create(p)?);
                write_sample_subset_vcf(path, in_fmt, options, version_lines.as_ref(), bgzf)
            }
            _ => {
                let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
                write_sample_subset_vcf(path, in_fmt, options, version_lines.as_ref(), bgzf)
            }
        },
        OutputKind::BcfUncompressed | OutputKind::BcfCompressed => match options.output_file {
            Some("-") | None => {
                write_sample_subset_bcf(path, in_fmt, options, version_lines.as_ref(), io::stdout())
            }
            Some(p) => write_sample_subset_bcf(
                path,
                in_fmt,
                options,
                version_lines.as_ref(),
                File::create(p)?,
            ),
        },
    }
}

fn write_sample_subset_vcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    options: &RunOptions<'_>,
    version_lines: Option<&crate::header_version::HeaderVersionLines>,
    mut out: W,
) -> io::Result<()> {
    let text = vcf_text_from_path(path, fmt)?;
    write_sample_subset_vcf_text(&text, options, version_lines, &mut out)
}

fn write_sample_subset_bcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    options: &RunOptions<'_>,
    version_lines: Option<&crate::header_version::HeaderVersionLines>,
    out: W,
) -> io::Result<()> {
    let text = vcf_text_from_path(path, fmt)?;
    let mut projected = Vec::new();
    let bcf_options = RunOptions {
        output_kind: options.output_kind,
        output_file: options.output_file,
        header_only: options.header_only,
        no_header: false,
        no_version: options.no_version,
        regions: options.regions,
        regions_overlap: options.regions_overlap,
        targets: options.targets,
        targets_exclude: options.targets_exclude,
        targets_overlap: options.targets_overlap,
        apply_filters: options.apply_filters.clone(),
        type_filter: options.type_filter,
        type_filter_exclude: options.type_filter_exclude,
        min_alleles: options.min_alleles,
        max_alleles: options.max_alleles,
        phased_filter: options.phased_filter,
        known_filter: options.known_filter,
        uncalled_filter: options.uncalled_filter,
        genotype_filter: options.genotype_filter,
        thread_count: options.thread_count,
        sample_list: options.sample_list,
        sample_list_is_file: options.sample_list_is_file,
        drop_genotypes: options.drop_genotypes,
    };
    write_sample_subset_vcf_text(&text, &bcf_options, version_lines, &mut projected)?;
    write_bcf_from_vcf_text(&projected, out)
}

fn write_bcf_from_vcf_text<W: Write>(text: &[u8], out: W) -> io::Result<()> {
    use htslib_rs::{bcf, vcf};

    let mut reader = vcf::io::Reader::new(BufReader::new(Cursor::new(text)));
    let header = reader.read_header()?;
    let mut writer = bcf::io::Writer::new(out);
    writer.write_variant_header(&header)?;
    for result in reader.records() {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }
    writer.try_finish()
}

fn vcf_text_from_path(path: &Path, fmt: format::Format) -> io::Result<String> {
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

fn write_sample_subset_vcf_text<W: Write>(
    text: &str,
    options: &RunOptions<'_>,
    version_lines: Option<&crate::header_version::HeaderVersionLines>,
    out: &mut W,
) -> io::Result<()> {
    let mut selected_samples: Option<Vec<usize>> = None;
    let mut inserted_version = false;
    let needs_pass_filter = options.drop_genotypes
        && !text
            .lines()
            .any(|line| line.starts_with("##FILTER=<ID=PASS,"));

    for line in text.split_inclusive('\n') {
        if line.starts_with("##") {
            if options.drop_genotypes && line.starts_with("##FORMAT=") {
                continue;
            }
            if !options.no_header {
                out.write_all(line.as_bytes())?;
                if needs_pass_filter && line.starts_with("##fileformat=") {
                    writeln!(out, "##FILTER=<ID=PASS,Description=\"All filters passed\">")?;
                }
            }
            continue;
        }

        if line.starts_with("#CHROM\t") {
            if !options.no_header {
                if let Some(lines) = version_lines
                    && !inserted_version
                {
                    writeln!(out, "{}", lines.version_line)?;
                    writeln!(out, "{}", lines.command_line)?;
                    inserted_version = true;
                }
                let fields = line_fields(line);
                let selected = selected_sample_indices(&fields, options)?;
                write_projected_vcf_line(&fields, &selected, !options.drop_genotypes, out)?;
                selected_samples = Some(selected);
            } else {
                let fields = line_fields(line);
                selected_samples = Some(selected_sample_indices(&fields, options)?);
            }
            continue;
        }

        if line.starts_with('#') {
            if !options.no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }

        if options.header_only {
            break;
        }
        let fields = line_fields(line);
        if !record_line_matches_filters(&fields, options) {
            continue;
        }
        let selected = selected_samples.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "VCF header is missing the #CHROM sample line",
            )
        })?;
        write_projected_vcf_line(&fields, selected, !options.drop_genotypes, out)?;
    }

    Ok(())
}

fn selected_sample_indices(fields: &[&str], options: &RunOptions<'_>) -> io::Result<Vec<usize>> {
    if options.drop_genotypes {
        return Ok(Vec::new());
    }
    let sample_names = fields[9..]
        .iter()
        .map(|s| s.to_string())
        .collect::<Vec<_>>();
    Ok(crate::smpl_ilist::init(
        &sample_names,
        options.sample_list,
        options.sample_list_is_file,
        crate::smpl_ilist::SMPL_STRICT,
    )?
    .idx)
}

fn line_fields(line: &str) -> Vec<&str> {
    line.trim_end_matches('\n')
        .trim_end_matches('\r')
        .split('\t')
        .collect()
}

fn write_projected_vcf_line<W: Write>(
    fields: &[&str],
    selected_samples: &[usize],
    keep_format_column: bool,
    out: &mut W,
) -> io::Result<()> {
    let fixed_end = fields.len().min(if keep_format_column { 9 } else { 8 });
    let mut projected = fields[..fixed_end].to_vec();
    for &sample_idx in selected_samples {
        if let Some(value) = fields.get(9 + sample_idx) {
            projected.push(value);
        }
    }
    writeln!(out, "{}", projected.join("\t"))
}

fn record_line_matches_filters(fields: &[&str], options: &RunOptions<'_>) -> bool {
    let Some(contig) = fields.first() else {
        return false;
    };
    let Some(pos) = fields.get(1).and_then(|pos| pos.parse::<usize>().ok()) else {
        return false;
    };
    region_line_matches(
        options.regions,
        options.regions_overlap,
        contig,
        pos,
        fields,
    ) && target_line_matches(
        options.targets,
        options.targets_exclude,
        options.targets_overlap,
        contig,
        pos,
        fields,
    ) && filter_line_matches(fields, options.apply_filters.as_deref())
        && allele_count_line_matches(fields, options.min_alleles, options.max_alleles)
        && known_line_matches(fields, options.known_filter)
        && phased_line_matches(fields, options.phased_filter)
        && uncalled_line_matches(fields, options.uncalled_filter)
        && genotype_line_matches(fields, options.genotype_filter)
        && variant_type_line_matches(fields, options.type_filter, options.type_filter_exclude)
}

fn filter_line_matches(fields: &[&str], apply_filters: Option<&[String]>) -> bool {
    let Some(filters) = apply_filters else {
        return true;
    };
    let Some(value) = fields.get(6) else {
        return false;
    };
    filters
        .iter()
        .any(|filter| filter_field_contains(value, filter))
}

fn filter_field_contains(value: &str, wanted: &str) -> bool {
    if value == wanted {
        return true;
    }
    if value == "." || value == "PASS" {
        return false;
    }
    value.split(';').any(|item| item == wanted)
}

fn genotype_line_matches(fields: &[&str], genotype_filter: Option<GenotypeFilter>) -> bool {
    let Some(filter) = genotype_filter else {
        return true;
    };
    let Some(format) = fields.get(8) else {
        return true;
    };
    let Some(gt_index) = format.split(':').position(|key| key == "GT") else {
        return true;
    };
    let has_class = fields[9..].iter().any(|sample| {
        sample
            .split(':')
            .nth(gt_index)
            .and_then(classify_sample_genotype)
            .is_some_and(|class| class == filter.class)
    });
    if filter.exclude {
        !has_class
    } else {
        has_class
    }
}

fn classify_sample_genotype(gt: &str) -> Option<GenotypeClass> {
    if gt.is_empty() {
        return None;
    }
    let alleles = gt
        .split(['/', '|'])
        .filter(|allele| !allele.is_empty())
        .collect::<Vec<_>>();
    if alleles.is_empty() || alleles.iter().all(|allele| *allele == ".") {
        return Some(GenotypeClass::Missing);
    }
    let called = alleles
        .iter()
        .copied()
        .filter(|allele| *allele != ".")
        .collect::<Vec<_>>();
    let first = called.first()?;
    if called.iter().all(|allele| allele == first) {
        Some(GenotypeClass::Hom)
    } else {
        Some(GenotypeClass::Het)
    }
}

fn uncalled_line_matches(fields: &[&str], uncalled_filter: Option<bool>) -> bool {
    let Some(include_uncalled) = uncalled_filter else {
        return true;
    };
    let uncalled = !record_has_called_genotype(fields);
    if include_uncalled {
        uncalled
    } else {
        !uncalled
    }
}

fn record_has_called_genotype(fields: &[&str]) -> bool {
    let Some(format) = fields.get(8) else {
        return false;
    };
    let Some(gt_index) = format.split(':').position(|key| key == "GT") else {
        return false;
    };
    fields[9..].iter().any(|sample| {
        sample
            .split(':')
            .nth(gt_index)
            .is_some_and(genotype_has_called_allele)
    })
}

fn genotype_has_called_allele(gt: &str) -> bool {
    gt.split(['/', '|'])
        .any(|allele| !allele.is_empty() && allele != ".")
}

fn known_line_matches(fields: &[&str], known_filter: Option<bool>) -> bool {
    let Some(include_known) = known_filter else {
        return true;
    };
    let known = fields.get(2).is_some_and(|id| *id != ".");
    if include_known { known } else { !known }
}

fn phased_line_matches(fields: &[&str], phased_filter: Option<bool>) -> bool {
    let Some(include_phased) = phased_filter else {
        return true;
    };
    let all_phased = all_samples_phased(fields);
    if include_phased {
        all_phased
    } else {
        !all_phased
    }
}

fn all_samples_phased(fields: &[&str]) -> bool {
    let Some(format) = fields.get(8) else {
        return true;
    };
    let Some(gt_index) = format.split(':').position(|key| key == "GT") else {
        return true;
    };
    fields[9..].iter().all(|sample| {
        sample
            .split(':')
            .nth(gt_index)
            .is_some_and(sample_genotype_phased)
    })
}

fn sample_genotype_phased(gt: &str) -> bool {
    if gt.is_empty() {
        return false;
    }
    if gt.contains('|') {
        return true;
    }
    !gt.contains('/')
}

fn allele_count_line_matches(
    fields: &[&str],
    min_alleles: Option<usize>,
    max_alleles: Option<usize>,
) -> bool {
    let n_alleles = allele_count_from_fields(fields);
    min_alleles.is_none_or(|min| n_alleles >= min) && max_alleles.is_none_or(|max| n_alleles <= max)
}

fn allele_count_from_fields(fields: &[&str]) -> usize {
    let Some(alts) = fields.get(4) else {
        return 1;
    };
    if *alts == "." || alts.is_empty() {
        1
    } else {
        1 + alts.split(',').count()
    }
}

fn variant_type_line_matches(
    fields: &[&str],
    type_filter: Option<TypeFilter>,
    exclude: bool,
) -> bool {
    let Some(filter) = type_filter else {
        return true;
    };
    let variant_type = variant_type_from_fields(fields);
    let matches = if variant_type.bits() == VariantType::REF.bits() {
        filter.include_ref
    } else {
        (variant_type & filter.mask).bits() != 0
    };
    if exclude { !matches } else { matches }
}

fn variant_type_from_fields(fields: &[&str]) -> VariantType {
    let Some(ref_allele) = fields.get(3) else {
        return VariantType::REF;
    };
    let Some(alts) = fields.get(4) else {
        return VariantType::REF;
    };
    let mut variant_type = VariantType::REF;
    for alt in alts.split(',') {
        variant_type |= classify_variant(ref_allele, alt).variant_type;
    }
    variant_type
}

fn region_line_matches(
    regions: &[Region],
    overlap: RegionOverlap,
    contig: &str,
    pos: usize,
    fields: &[&str],
) -> bool {
    regions.is_empty()
        || regions.iter().any(|region| {
            region.contig == contig && record_overlaps_target(pos, fields, region, overlap)
        })
}

fn target_line_matches(
    targets: &[Region],
    exclude: bool,
    overlap: RegionOverlap,
    contig: &str,
    pos: usize,
    fields: &[&str],
) -> bool {
    if targets.is_empty() {
        return true;
    }
    let matches = targets.iter().any(|target| {
        target.contig == contig && record_overlaps_target(pos, fields, target, overlap)
    });
    if exclude { !matches } else { matches }
}

fn record_overlaps_target(
    pos: usize,
    fields: &[&str],
    target: &Region,
    overlap: RegionOverlap,
) -> bool {
    let start = target.start.unwrap_or(1);
    let end = target.end.unwrap_or(usize::MAX);
    match overlap {
        RegionOverlap::Pos => pos >= start && pos <= end,
        RegionOverlap::Record => {
            let record_end = record_end_from_fields(pos, fields);
            pos <= end && start <= record_end
        }
        RegionOverlap::Variant => variant_overlaps_target(pos, fields, start, end),
    }
}

fn record_end_from_fields(pos: usize, fields: &[&str]) -> usize {
    let ref_len = fields.get(3).map(|s| s.len().max(1)).unwrap_or(1);
    fields
        .get(7)
        .and_then(|info| info.split(';').find_map(|field| field.strip_prefix("END=")))
        .and_then(|value| value.parse::<usize>().ok())
        .map(|end| end.max(pos))
        .unwrap_or(pos + ref_len - 1)
}

fn variant_overlaps_target(
    pos: usize,
    fields: &[&str],
    target_start: usize,
    target_end: usize,
) -> bool {
    if let Some(end) = fields
        .get(7)
        .and_then(|info| info.split(';').find_map(|field| field.strip_prefix("END=")))
        .and_then(|value| value.parse::<usize>().ok())
    {
        return pos <= target_end && target_start <= end;
    }

    let ref_len = fields.get(3).map(|s| s.len().max(1)).unwrap_or(1);
    let Some(alts) = fields.get(4) else {
        return false;
    };
    alts.split(',').any(|alt| {
        let alt_len = alt.len().max(1);
        if ref_len == alt_len {
            pos >= target_start && pos <= target_end
        } else if ref_len > alt_len {
            let start = pos + alt_len;
            let end = pos + ref_len - 1;
            start <= target_end && target_start <= end
        } else {
            pos >= target_start && pos <= target_end && pos < target_end
        }
    })
}

fn write_vcf_gz(
    path: &Path,
    in_fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    options: &RunOptions<'_>,
) -> io::Result<()> {
    match (options.output_file, options.thread_count) {
        (Some(p), Some(thread_count)) if p != "-" => {
            let file = File::create(p)?;
            let bgzf =
                htslib_rs::bgzf::io::MultithreadedWriter::with_worker_count(thread_count, file);
            write_vcf(
                path,
                in_fmt,
                header,
                options.header_only,
                options.no_header,
                RecordFilters {
                    regions: options.regions,
                    targets: options.targets,
                    targets_exclude: options.targets_exclude,
                },
                bgzf,
            )
        }
        (Some(p), _) if p != "-" => {
            let bgzf = htslib_rs::bgzf::io::Writer::new(File::create(p)?);
            write_vcf(
                path,
                in_fmt,
                header,
                options.header_only,
                options.no_header,
                RecordFilters {
                    regions: options.regions,
                    targets: options.targets,
                    targets_exclude: options.targets_exclude,
                },
                bgzf,
            )
        }
        _ => {
            let bgzf = htslib_rs::bgzf::io::Writer::new(io::stdout().lock());
            write_vcf(
                path,
                in_fmt,
                header,
                options.header_only,
                options.no_header,
                RecordFilters {
                    regions: options.regions,
                    targets: options.targets,
                    targets_exclude: options.targets_exclude,
                },
                bgzf,
            )
        }
    }
}

fn write_vcf_text_passthrough<W: Write>(
    path: &Path,
    fmt: format::Format,
    header_only: bool,
    no_header: bool,
    out: W,
) -> io::Result<()> {
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        return write_vcf_text_passthrough_reader(BufReader::new(dec), header_only, no_header, out);
    }
    let reader = File::open(path).map(BufReader::new)?;
    write_vcf_text_passthrough_reader(reader, header_only, no_header, out)
}

fn write_vcf_text_filtered_passthrough<W: Write>(
    path: &Path,
    fmt: format::Format,
    options: &RunOptions<'_>,
    out: W,
) -> io::Result<()> {
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        return write_vcf_text_filtered_passthrough_reader(BufReader::new(dec), options, out);
    }
    let reader = File::open(path).map(BufReader::new)?;
    write_vcf_text_filtered_passthrough_reader(reader, options, out)
}

fn write_vcf_text_filtered_passthrough_reader<R, W>(
    mut reader: R,
    options: &RunOptions<'_>,
    mut out: W,
) -> io::Result<()>
where
    R: io::BufRead,
    W: Write,
{
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with('#') {
            if !options.no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }
        if options.header_only {
            break;
        }
        let fields = line_fields(&line);
        if record_line_matches_filters(&fields, options) {
            out.write_all(line.as_bytes())?;
        }
    }
    Ok(())
}

fn write_vcf_text_passthrough_reader<R, W>(
    mut reader: R,
    header_only: bool,
    no_header: bool,
    mut out: W,
) -> io::Result<()>
where
    R: io::BufRead,
    W: Write,
{
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line)?;
        if n == 0 {
            break;
        }
        if line.starts_with('#') {
            if !no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }
        if header_only {
            break;
        }
        out.write_all(line.as_bytes())?;
        io::copy(&mut reader, &mut out)?;
        break;
    }
    Ok(())
}

fn write_bcf_vcf_text_no_version<W: Write>(
    path: &Path,
    header_only: bool,
    no_header: bool,
    mut out: W,
) -> io::Result<()> {
    let text = htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None)?;
    write_vcf_text_from_string(&text, header_only, no_header, &mut out)
}

fn write_vcf_text_from_string<W: Write>(
    text: &str,
    header_only: bool,
    no_header: bool,
    out: &mut W,
) -> io::Result<()> {
    for line in text.split_inclusive('\n') {
        if line.starts_with('#') {
            if !no_header {
                out.write_all(line.as_bytes())?;
            }
            continue;
        }
        if header_only {
            break;
        }
        out.write_all(line.as_bytes())?;
    }
    Ok(())
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
    filters: RecordFilters<'_>,
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
    write_records_into_vcf(path, fmt, header, filters, &mut writer)
}

fn write_records_into_vcf<W: Write>(
    path: &Path,
    fmt: format::Format,
    header: &htslib_rs::vcf::Header,
    filters: RecordFilters<'_>,
    writer: &mut htslib_rs::vcf::io::Writer<W>,
) -> io::Result<()> {
    use htslib_rs::bcf;
    use htslib_rs::vcf;

    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        for result in reader.record_bufs(header) {
            let rec = result?;
            if record_matches(filters, rec.reference_sequence_name(), rec.variant_start()) {
                writer.write_variant_record(header, &rec)?;
            }
        }
    } else if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let mut reader = vcf::io::Reader::new(BufReader::new(dec));
        let _h = reader.read_header()?;
        for result in reader.records() {
            let rec = result?;
            if record_matches_result(filters, rec.reference_sequence_name(), rec.variant_start())? {
                writer.write_variant_record(header, &rec)?;
            }
        }
    } else {
        let mut reader = File::open(path)
            .map(BufReader::new)
            .map(vcf::io::Reader::new)?;
        let _h = reader.read_header()?;
        for result in reader.records() {
            let rec = result?;
            if record_matches_result(filters, rec.reference_sequence_name(), rec.variant_start())? {
                writer.write_variant_record(header, &rec)?;
            }
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
    filters: RecordFilters<'_>,
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
            if record_matches(filters, rec.reference_sequence_name(), rec.variant_start()) {
                writer.write_variant_record(header, &rec)?;
            }
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
            if record_matches_result(filters, rec.reference_sequence_name(), rec.variant_start())? {
                writer.write_variant_record(&header, &rec)?;
            }
        }
        writer.try_finish()?;
        Ok(())
    } else {
        if filters.regions.is_empty() && filters.targets.is_empty() {
            // Plain VCF → BCF: delegate to htslib-rs's tested helper.
            htslib_rs::variant_io_compat::write_bcf_from_vcf_path(path, out)?;
        } else {
            use htslib_rs::vcf;
            let mut reader = File::open(path)
                .map(BufReader::new)
                .map(vcf::io::Reader::new)?;
            let header = reader.read_header()?;
            let mut writer = bcf::io::Writer::new(out);
            writer.write_variant_header(&header)?;
            for result in reader.records() {
                let rec = result?;
                if record_matches_result(
                    filters,
                    rec.reference_sequence_name(),
                    rec.variant_start(),
                )? {
                    writer.write_variant_record(&header, &rec)?;
                }
            }
            writer.try_finish()?;
        }
        Ok(())
    }
}

fn record_matches(filters: RecordFilters<'_>, contig: &str, pos: Option<Position>) -> bool {
    pos.map(|pos| {
        region_matches(filters.regions, contig, pos)
            && target_matches(filters.targets, filters.targets_exclude, contig, pos)
    })
    .unwrap_or(false)
}

fn region_matches(regions: &[Region], contig: &str, pos: Position) -> bool {
    regions.is_empty() || regions.iter().any(|region| region.contains(contig, pos))
}

fn target_matches(targets: &[Region], exclude: bool, contig: &str, pos: Position) -> bool {
    if targets.is_empty() {
        return true;
    }
    let matches = region_matches(targets, contig, pos);
    if exclude { !matches } else { matches }
}

fn record_matches_result(
    filters: RecordFilters<'_>,
    contig: &str,
    pos: Option<io::Result<Position>>,
) -> io::Result<bool> {
    match pos {
        Some(Ok(pos)) => Ok(record_matches(filters, contig, Some(pos))),
        Some(Err(e)) => Err(e),
        None => Ok(record_matches(filters, contig, None)),
    }
}
