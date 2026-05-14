//! Focused `bcftools isec` implementation (upstream `vcfisec.c`).
//!
//! This text-backed slice covers the common stdout modes used by the upstream
//! `test_vcf_isec` fixtures: set bitmaps, `-w` VCF record output, simple
//! collapse modes, `-n` cardinality filters, `-C`, `-i`/`-e`, and POS-based
//! `-r`/`-R`/`-t`/`-T`. Full multi-reader synced iteration, indexed overlap
//! semantics, prefix directory output, and BCF output remain tracked in
//! `TODO.md` because they depend on broader synced-reader parity.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::ffi::OsString;
use std::fs::{self, File};
use std::io::{self, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    build_bcf_csi_with_min_shift, build_vcf_tbi_from_path, write_csi, write_tbi,
};

use crate::diagnostics::fmt_etag;
use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::header_version::{build_lines, command_time};
use crate::vcf_compat::NormalizeFileformat;

const USAGE: &str = "\n\
About:   Create intersections, unions and complements of VCF files.\n\
Usage:   bcftools isec [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n\
\n\
Options:\n\
   -c, --collapse STRING       Records are compatible by none|snps|indels|both|all|some|id|any [none]\n\
   -C, --complement            Output records private to the first file\n\
   -e, --exclude EXPR          Exclude sites for which the expression is true\n\
   -i, --include EXPR          Include only sites for which the expression is true\n\
   -n, --nfiles [+-=]INT       Output positions present in this many files\n\
       --no-version            Do not append version and command line to VCF headers\n\
   -O, --output-type v|z       Output VCF text or compressed VCF when writing records [v]\n\
   -o, --output FILE           Write site bitmap output to a file [standard output]\n\
   -p, --prefix DIR            Write numbered VCFs plus sites.txt and README.txt to DIR\n\
   -r, --regions REGION        Restrict to comma-separated list of regions\n\
   -R, --regions-file FILE     Restrict to regions listed in a file\n\
   -t, --targets REGION        Similar to -r but streams\n\
   -T, --targets-file FILE     Similar to -R but streams\n\
   -w, --write LIST            Write records from 1-based input indexes as VCF\n\
\n";

#[derive(Debug)]
struct Args {
    inputs: Vec<PathBuf>,
    nfiles: Option<NFiles>,
    complement: bool,
    collapse: CollapseMode,
    write_inputs: Vec<usize>,
    output: Option<PathBuf>,
    prefix: Option<PathBuf>,
    include_expr: Option<String>,
    exclude_expr: Option<String>,
    regions: Vec<RegionSpec>,
    targets: Vec<RegionSpec>,
    output_kind: OutputKind,
    no_version: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NOp {
    Eq,
    Ge,
    Le,
    Exact,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct NFiles {
    op: NOp,
    count: usize,
    exact: Option<Vec<bool>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CollapseMode {
    None,
    Any,
    Both,
    Id,
}

#[derive(Debug, Clone)]
struct RegionSpec {
    contig: String,
    start: Option<i64>,
    end: Option<i64>,
}

#[derive(Debug, Clone)]
struct InputVcf {
    header: String,
    records: Vec<Record>,
}

#[derive(Debug, Clone)]
struct Record {
    line: String,
    key: Key,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct Key {
    chrom: String,
    pos: i64,
    id: String,
    ref_allele: String,
    alt: String,
    class: VariantClass,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum VariantClass {
    Snp,
    Indel,
    Other,
}

#[derive(Debug)]
struct Group {
    id: String,
    present: Vec<bool>,
    representatives: Vec<Option<usize>>,
    first_input: usize,
    first_record: usize,
}

#[derive(Debug)]
enum ParseOutcome {
    Usage,
    Error(String),
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match parse_args(argv) {
        Ok(args) => match run(&args, argv) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfisec", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfisec", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut inputs = Vec::new();
    let mut nfiles = None;
    let mut complement = false;
    let mut collapse = CollapseMode::None;
    let mut write_inputs = Vec::new();
    let mut output = None;
    let mut prefix = None;
    let mut include_expr = None;
    let mut exclude_expr = None;
    let mut regions = Vec::new();
    let mut targets = Vec::new();
    let mut output_kind = OutputKind::VcfText;
    let mut no_version = false;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "--no-version" => no_version = true,
            "-C" | "--complement" => complement = true,
            "-c" | "--collapse" => {
                collapse = parse_collapse(&next_string(&mut iter, raw.as_ref())?)?;
            }
            "-e" | "--exclude" => exclude_expr = Some(next_string(&mut iter, raw.as_ref())?),
            "-i" | "--include" => include_expr = Some(next_string(&mut iter, raw.as_ref())?),
            "-n" | "--nfiles" => {
                nfiles = Some(parse_nfiles(&next_string(&mut iter, raw.as_ref())?)?);
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?;
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?));
            }
            "-p" | "--prefix" => {
                prefix = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?));
            }
            "-r" | "--regions" => {
                parse_region_list(&mut regions, &next_string(&mut iter, raw.as_ref())?)?;
            }
            "-R" | "--regions-file" => {
                load_region_file(&mut regions, &next_string(&mut iter, raw.as_ref())?)?;
            }
            "-t" | "--targets" => {
                parse_region_list(&mut targets, &next_string(&mut iter, raw.as_ref())?)?;
            }
            "-T" | "--targets-file" => {
                load_region_file(&mut targets, &next_string(&mut iter, raw.as_ref())?)?;
            }
            "-w" | "--write" => {
                write_inputs = parse_write_inputs(&next_string(&mut iter, raw.as_ref())?)?;
            }
            _ if raw.starts_with("--collapse=") => {
                collapse = parse_collapse(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--exclude=") => {
                exclude_expr = Some(value_after_equals(&raw).to_owned());
            }
            _ if raw.starts_with("--include=") => {
                include_expr = Some(value_after_equals(&raw).to_owned());
            }
            _ if raw.starts_with("--nfiles=") => {
                nfiles = Some(parse_nfiles(value_after_equals(&raw))?);
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--prefix=") => {
                prefix = Some(PathBuf::from(value_after_equals(&raw)));
            }
            _ if raw.starts_with("--regions=") => {
                parse_region_list(&mut regions, value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--regions-file=") => {
                load_region_file(&mut regions, value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--targets=") => {
                parse_region_list(&mut targets, value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--targets-file=") => {
                load_region_file(&mut targets, value_after_equals(&raw))?;
            }
            _ if raw.starts_with("--write=") => {
                write_inputs = parse_write_inputs(value_after_equals(&raw))?;
            }
            _ if raw.starts_with("-c") && raw.len() > 2 => {
                collapse = parse_collapse(&raw[2..])?;
            }
            _ if raw.starts_with("-e") && raw.len() > 2 => {
                exclude_expr = Some(raw[2..].to_owned());
            }
            _ if raw.starts_with("-i") && raw.len() > 2 => {
                include_expr = Some(raw[2..].to_owned());
            }
            _ if raw.starts_with("-n") && raw.len() > 2 => {
                nfiles = Some(parse_nfiles(&raw[2..])?);
            }
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?;
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => {
                output = Some(PathBuf::from(&raw[2..]));
            }
            _ if raw.starts_with("-p") && raw.len() > 2 => {
                prefix = Some(PathBuf::from(&raw[2..]));
            }
            _ if raw.starts_with("-r") && raw.len() > 2 => {
                parse_region_list(&mut regions, &raw[2..])?;
            }
            _ if raw.starts_with("-t") && raw.len() > 2 => {
                parse_region_list(&mut targets, &raw[2..])?;
            }
            _ if raw.starts_with("-w") && raw.len() > 2 => {
                write_inputs = parse_write_inputs(&raw[2..])?;
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!("unrecognized option '{raw}'")));
            }
            _ => inputs.push(PathBuf::from(raw.as_ref())),
        }
    }

    if inputs.is_empty() {
        return Err(ParseOutcome::Error(
            "expected at least one input VCF".into(),
        ));
    }
    if let Some(mask) = nfiles.as_ref().and_then(|n| n.exact.as_ref())
        && mask.len() != inputs.len()
    {
        return Err(ParseOutcome::Error(format!(
            "the number of files does not match the bitmask: {} vs {}",
            inputs.len(),
            mask.iter()
                .map(|present| if *present { '1' } else { '0' })
                .collect::<String>()
        )));
    }
    if matches!(output_kind, OutputKind::VcfGz | OutputKind::Bcf)
        && write_inputs.is_empty()
        && inputs.len() != 1
        && prefix.is_none()
    {
        return Err(ParseOutcome::Error(
            "compressed VCF/BCF output is only supported for record-writing modes".into(),
        ));
    }
    if write_inputs.len() > 1 && prefix.is_none() {
        return Err(ParseOutcome::Error(
            "expected -p when multiple output files are requested with -w".into(),
        ));
    }
    for &idx in &write_inputs {
        if idx == 0 || idx > inputs.len() {
            return Err(ParseOutcome::Error(format!(
                "-w index {idx} is outside the input range"
            )));
        }
    }

    Ok(Args {
        inputs,
        nfiles,
        complement,
        collapse,
        write_inputs,
        output,
        prefix,
        include_expr,
        exclude_expr,
        regions,
        targets,
        output_kind,
        no_version,
    })
}

fn run(args: &Args, argv: &[OsString]) -> io::Result<()> {
    let inputs = args
        .inputs
        .iter()
        .map(|path| read_input(path, args))
        .collect::<io::Result<Vec<_>>>()?;
    let groups = build_groups(&inputs, args.collapse);
    let selected = if is_prefix_venn(args, inputs.len()) {
        (0..groups.len()).collect()
    } else {
        selected_groups(&groups, args)
    };

    if let Some(prefix) = &args.prefix {
        write_prefix_outputs(prefix, &inputs, &groups, &selected, args, argv)
    } else if !args.write_inputs.is_empty() {
        let mut out = Vec::new();
        write_record_outputs(&mut out, &inputs, &groups, &selected, args, argv)?;
        write_output(&out, args.output_kind)
    } else if inputs.len() == 1 {
        let mut out = Vec::new();
        write_single_input_vcf(&mut out, &inputs[0], &groups, &selected, args, argv)?;
        write_output(&out, args.output_kind)
    } else {
        write_bitmap_destination(&inputs, &groups, &selected, args.output.as_deref())
    }
}

fn read_input(path: &Path, args: &Args) -> io::Result<InputVcf> {
    let text = read_vcf_text(path)?;
    let mut header = String::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            header.push_str(line);
            header.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if fields.len() < 8 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid VCF record with fewer than 8 columns: {line}"),
            ));
        }
        if !record_in_regions(&fields, &args.regions, &args.targets) {
            continue;
        }
        if !evaluate(&fields, args)? {
            continue;
        }
        records.push(Record {
            line: line.to_owned(),
            key: record_key(&fields),
        });
    }
    Ok(InputVcf { header, records })
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    let mut text = String::new();
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }
    let file = File::open(path)?;
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

fn record_key(fields: &[String]) -> Key {
    let ref_allele = fields[3].clone();
    let alt = fields[4].clone();
    Key {
        chrom: fields[0].clone(),
        pos: fields[1].parse().unwrap_or(0),
        id: fields[2].clone(),
        class: variant_class(&ref_allele, &alt),
        ref_allele,
        alt,
    }
}

fn variant_class(reference: &str, alt: &str) -> VariantClass {
    let mut saw_snp = false;
    let mut saw_indel = false;
    for allele in alt.split(',') {
        if reference.len() == 1 && allele.len() == 1 {
            saw_snp = true;
        } else if reference.len() != allele.len() {
            saw_indel = true;
        }
    }
    match (saw_snp, saw_indel) {
        (true, false) => VariantClass::Snp,
        (false, true) => VariantClass::Indel,
        _ => VariantClass::Other,
    }
}

fn collapse_key(key: &Key, mode: CollapseMode) -> String {
    match mode {
        CollapseMode::None => format!(
            "{}\t{}\t{}\t{}",
            key.chrom, key.pos, key.ref_allele, key.alt
        ),
        CollapseMode::Id if key.id != "." => format!("{}\t{}\tID\t{}", key.chrom, key.pos, key.id),
        CollapseMode::Id => format!(
            "{}\t{}\t{}\t{}",
            key.chrom, key.pos, key.ref_allele, key.alt
        ),
        CollapseMode::Any => format!("{}\t{}", key.chrom, key.pos),
        CollapseMode::Both => format!("{}\t{}", key.chrom, key.pos),
    }
}

fn build_groups(inputs: &[InputVcf], collapse: CollapseMode) -> Vec<Group> {
    let mut groups = Vec::new();
    let mut by_key: BTreeMap<String, usize> = BTreeMap::new();
    for (input_i, input) in inputs.iter().enumerate() {
        let mut seen_in_input = HashSet::new();
        for (record_i, record) in input.records.iter().enumerate() {
            let key = collapse_key(&record.key, collapse);
            let group_i = *by_key.entry(key.clone()).or_insert_with(|| {
                let i = groups.len();
                groups.push(Group {
                    id: key.clone(),
                    present: vec![false; inputs.len()],
                    representatives: vec![None; inputs.len()],
                    first_input: input_i,
                    first_record: record_i,
                });
                i
            });
            if seen_in_input.insert(key) {
                groups[group_i].present[input_i] = true;
                groups[group_i].representatives[input_i] = Some(record_i);
            }
        }
    }
    groups
}

fn selected_groups(groups: &[Group], args: &Args) -> HashSet<usize> {
    groups
        .iter()
        .enumerate()
        .filter_map(|(i, group)| {
            let count = group.present.iter().filter(|&&present| present).count();
            let keep = if args.complement {
                group.present.first().copied().unwrap_or(false) && count == 1
            } else if let Some(nfiles) = &args.nfiles {
                match nfiles.op {
                    NOp::Eq => count == nfiles.count,
                    NOp::Ge => count >= nfiles.count,
                    NOp::Le => count <= nfiles.count,
                    NOp::Exact => nfiles
                        .exact
                        .as_ref()
                        .is_some_and(|mask| mask == &group.present),
                }
            } else {
                count == args.inputs.len()
            };
            keep.then_some(i)
        })
        .collect()
}

fn write_bitmap_output<W: Write>(
    out: &mut W,
    inputs: &[InputVcf],
    groups: &[Group],
    selected: &HashSet<usize>,
) -> io::Result<()> {
    if inputs
        .first()
        .is_some_and(|input| !input.records.is_empty())
    {
        let lookup = record_to_group_lookup(inputs, groups);
        for (record_i, record) in inputs[0].records.iter().enumerate() {
            let Some(group_i) = lookup.get(&(0, record_i)).copied() else {
                continue;
            };
            if !selected.contains(&group_i) {
                continue;
            }
            let group = &groups[group_i];
            let bitmap: String = group
                .present
                .iter()
                .map(|present| if *present { '1' } else { '0' })
                .collect();
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                record.key.chrom, record.key.pos, record.key.ref_allele, record.key.alt, bitmap
            )?;
        }
        return Ok(());
    }

    for (group_i, group) in groups.iter().enumerate() {
        if !selected.contains(&group_i) {
            continue;
        }
        let input_i = group.first_input;
        let record_i = group.first_record;
        let record = &inputs[input_i].records[record_i];
        let bitmap: String = group
            .present
            .iter()
            .map(|present| if *present { '1' } else { '0' })
            .collect();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            record.key.chrom, record.key.pos, record.key.ref_allele, record.key.alt, bitmap
        )?;
    }
    Ok(())
}

fn write_bitmap_destination(
    inputs: &[InputVcf],
    groups: &[Group],
    selected: &HashSet<usize>,
    output: Option<&Path>,
) -> io::Result<()> {
    match output {
        Some(path) => {
            let mut file = File::create(path)?;
            write_bitmap_output(&mut file, inputs, groups, selected)
        }
        None => {
            let mut stdout = io::stdout().lock();
            write_bitmap_output(&mut stdout, inputs, groups, selected)
        }
    }
}

fn write_prefix_outputs(
    prefix: &Path,
    inputs: &[InputVcf],
    groups: &[Group],
    selected: &HashSet<usize>,
    args: &Args,
    argv: &[OsString],
) -> io::Result<()> {
    fs::create_dir_all(prefix)?;

    if is_prefix_venn(args, inputs.len()) {
        return write_prefix_venn_outputs(prefix, inputs, groups, selected, args, argv);
    }

    let mut readme = String::new();
    readme.push_str("This file was produced by vcfisec.\n");
    readme.push_str("The command line was:\tbcftools isec");
    for arg in argv.iter().skip(1) {
        readme.push(' ');
        readme.push_str(&arg.to_string_lossy());
    }
    readme.push_str("\n\nUsing the following file names:\n");

    let site_path = prefix.join("sites.txt");
    {
        let mut sites = File::create(&site_path)?;
        write_bitmap_output(&mut sites, inputs, groups, selected)?;
    }

    let group_lookup = record_to_group_lookup(inputs, groups);
    let suffix = match args.output_kind {
        OutputKind::VcfText => "vcf",
        OutputKind::VcfGz => "vcf.gz",
        OutputKind::Bcf => "bcf",
    };
    for input_i in prefix_output_inputs(inputs.len(), args) {
        let input = &inputs[input_i];
        let filename = format!("{input_i:04}.{suffix}");
        let path = prefix.join(filename);
        let mut bytes = Vec::new();
        write_header(&mut bytes, &input.header, args.no_version, argv)?;
        let mut record_indices: Vec<usize> = input
            .records
            .iter()
            .enumerate()
            .filter_map(|(record_i, _)| {
                let group_i = group_lookup.get(&(input_i, record_i)).copied()?;
                selected.contains(&group_i).then_some(record_i)
            })
            .collect();
        if args.collapse == CollapseMode::Id {
            record_indices.sort_by(|&a, &b| {
                input.records[a]
                    .key
                    .id
                    .cmp(&input.records[b].key.id)
                    .then_with(|| a.cmp(&b))
            });
        }
        for record_i in record_indices {
            writeln!(&mut bytes, "{}", input.records[record_i].line)?;
        }
        write_bytes_to_path(&path, &bytes, args.output_kind)?;
        maybe_index_prefix_vcf_gz(&path, args.output_kind)?;
        readme.push_str(&format!(
            "{}\tfor stripped\t{}\n",
            path.display(),
            args.inputs[input_i].display()
        ));
    }

    fs::write(prefix.join("README.txt"), readme)
}

fn write_prefix_venn_outputs(
    prefix: &Path,
    inputs: &[InputVcf],
    groups: &[Group],
    selected: &HashSet<usize>,
    args: &Args,
    argv: &[OsString],
) -> io::Result<()> {
    let mut readme = String::new();
    readme.push_str("This file was produced by vcfisec.\n");
    readme.push_str("The command line was:\tbcftools isec");
    for arg in argv.iter().skip(1) {
        readme.push(' ');
        readme.push_str(&arg.to_string_lossy());
    }
    readme.push_str("\n\nUsing the following file names:\n");

    {
        let mut sites = File::create(prefix.join("sites.txt"))?;
        write_bitmap_all_groups(&mut sites, inputs, groups, selected)?;
    }

    let suffix = match args.output_kind {
        OutputKind::VcfText => "vcf",
        OutputKind::VcfGz => "vcf.gz",
        OutputKind::Bcf => "bcf",
    };
    let outputs = [
        (0usize, 0usize, "private"),
        (1usize, 1usize, "private"),
        (2usize, 0usize, "shared"),
        (3usize, 1usize, "shared"),
    ];
    for (file_i, input_i, kind) in outputs {
        if !prefix_venn_file_requested(file_i, input_i, args) {
            continue;
        }
        let input = &inputs[input_i];
        let filename = format!("{file_i:04}.{suffix}");
        let path = prefix.join(filename);
        let mut bytes = Vec::new();
        write_header(&mut bytes, &input.header, args.no_version, argv)?;
        for (record_i, record) in input.records.iter().enumerate() {
            let Some(group_i) = group_for_record(inputs, groups, input_i, record_i) else {
                continue;
            };
            if !selected.contains(&group_i) {
                continue;
            }
            let group = &groups[group_i];
            let keep = match kind {
                "private" => group.present[input_i] && !group.present[1 - input_i],
                "shared" => group.present[0] && group.present[1],
                _ => false,
            };
            if keep {
                writeln!(&mut bytes, "{}", record.line)?;
            }
        }
        write_bytes_to_path(&path, &bytes, args.output_kind)?;
        maybe_index_prefix_vcf_gz(&path, args.output_kind)?;
        let description = match kind {
            "private" => format!("for records private to\t{}", args.inputs[input_i].display()),
            "shared" => format!(
                "for records from {} shared by both\t{} {}",
                args.inputs[input_i].display(),
                args.inputs[0].display(),
                args.inputs[1].display()
            ),
            _ => String::new(),
        };
        readme.push_str(&format!("{}\t{description}\n", path.display()));
    }

    fs::write(prefix.join("README.txt"), readme)
}

fn write_bitmap_all_groups<W: Write>(
    out: &mut W,
    inputs: &[InputVcf],
    groups: &[Group],
    selected: &HashSet<usize>,
) -> io::Result<()> {
    for (group_i, group) in groups.iter().enumerate() {
        if !selected.contains(&group_i) {
            continue;
        }
        let input_i = group.first_input;
        let record_i = group.first_record;
        let record = &inputs[input_i].records[record_i];
        let bitmap: String = group
            .present
            .iter()
            .map(|present| if *present { '1' } else { '0' })
            .collect();
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}",
            record.key.chrom, record.key.pos, record.key.ref_allele, record.key.alt, bitmap
        )?;
    }
    Ok(())
}

fn is_prefix_venn(args: &Args, n_inputs: usize) -> bool {
    args.prefix.is_some() && n_inputs == 2 && args.nfiles.is_none() && !args.complement
}

fn prefix_venn_file_requested(file_i: usize, input_i: usize, args: &Args) -> bool {
    if args.write_inputs.is_empty() {
        return true;
    }
    let requested = args.write_inputs.iter().any(|&idx| idx == input_i + 1);
    match file_i {
        0 | 2 if input_i == 0 => requested,
        1 | 3 if input_i == 1 => requested,
        _ => false,
    }
}

fn prefix_output_inputs(n_inputs: usize, args: &Args) -> Vec<usize> {
    if args.complement {
        return vec![0];
    }
    if args.write_inputs.is_empty() {
        return (0..n_inputs).collect();
    }
    args.write_inputs.iter().map(|idx| idx - 1).collect()
}

fn write_record_outputs<W: Write>(
    out: &mut W,
    inputs: &[InputVcf],
    groups: &[Group],
    selected: &HashSet<usize>,
    args: &Args,
    argv: &[OsString],
) -> io::Result<()> {
    let group_lookup = record_to_group_lookup(inputs, groups);
    for &one_based in &args.write_inputs {
        let input_i = one_based - 1;
        let input = &inputs[input_i];
        write_header(out, &input.header, args.no_version, argv)?;
        let mut record_indices: Vec<usize> = input
            .records
            .iter()
            .enumerate()
            .filter_map(|(record_i, _)| {
                let group_i = group_lookup.get(&(input_i, record_i)).copied()?;
                selected.contains(&group_i).then_some(record_i)
            })
            .collect();
        if args.collapse == CollapseMode::Id {
            record_indices.sort_by(|&a, &b| {
                input.records[a]
                    .key
                    .id
                    .cmp(&input.records[b].key.id)
                    .then_with(|| a.cmp(&b))
            });
        }
        for record_i in record_indices {
            writeln!(out, "{}", input.records[record_i].line)?;
        }
    }
    Ok(())
}

fn write_single_input_vcf<W: Write>(
    out: &mut W,
    input: &InputVcf,
    groups: &[Group],
    selected: &HashSet<usize>,
    args: &Args,
    argv: &[OsString],
) -> io::Result<()> {
    write_header(out, &input.header, args.no_version, argv)?;
    let group_lookup = record_to_group_lookup(std::slice::from_ref(input), groups);
    for (record_i, record) in input.records.iter().enumerate() {
        let Some(group_i) = group_lookup.get(&(0, record_i)).copied() else {
            continue;
        };
        if selected.contains(&group_i) {
            writeln!(out, "{}", record.line)?;
        }
    }
    Ok(())
}

fn record_to_group_lookup(inputs: &[InputVcf], groups: &[Group]) -> HashMap<(usize, usize), usize> {
    let mut lookup = HashMap::new();
    for (group_i, group) in groups.iter().enumerate() {
        for (input_i, input) in inputs.iter().enumerate() {
            for (record_i, record) in input.records.iter().enumerate() {
                if group.id == collapse_key(&record.key, infer_group_mode(&group.id, &record.key)) {
                    lookup.insert((input_i, record_i), group_i);
                }
            }
        }
    }
    lookup
}

fn group_for_record(
    inputs: &[InputVcf],
    groups: &[Group],
    input_i: usize,
    record_i: usize,
) -> Option<usize> {
    let record = inputs.get(input_i)?.records.get(record_i)?;
    groups.iter().position(|group| {
        group.id == collapse_key(&record.key, infer_group_mode(&group.id, &record.key))
    })
}

fn infer_group_mode(group_id: &str, key: &Key) -> CollapseMode {
    if group_id == format!("{}\t{}", key.chrom, key.pos) {
        CollapseMode::Any
    } else if key.id != "." && group_id == format!("{}\t{}\tID\t{}", key.chrom, key.pos, key.id) {
        CollapseMode::Id
    } else {
        CollapseMode::None
    }
}

fn write_header<W: Write>(
    out: &mut W,
    header: &str,
    no_version: bool,
    argv: &[OsString],
) -> io::Result<()> {
    let header = ensure_pass_filter(header);
    for line in header.lines() {
        if !no_version && line.starts_with("#CHROM\t") {
            let mut prog_argv: Vec<OsString> = vec!["bcftools".into()];
            prog_argv.extend(argv.iter().cloned());
            let lines = build_lines("bcftools_isec", &prog_argv, command_time());
            writeln!(out, "{}", lines.version_line)?;
            writeln!(out, "{}", lines.command_line)?;
        }
        writeln!(out, "{line}")?;
    }
    Ok(())
}

fn ensure_pass_filter(header: &str) -> String {
    if header
        .lines()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"))
    {
        return header.to_owned();
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
    out
}

fn write_output(bytes: &[u8], kind: OutputKind) -> io::Result<()> {
    let stdout = io::stdout();
    let mut out = stdout.lock();
    match kind {
        OutputKind::VcfText => out.write_all(bytes),
        OutputKind::VcfGz => {
            let mut bgzf = htslib_rs::bgzf::io::Writer::new(out);
            bgzf.write_all(bytes)?;
            bgzf.finish().map(|_| ())
        }
        OutputKind::Bcf => write_bcf_from_vcf_text(bytes, out),
    }
}

fn write_bytes_to_path(path: &Path, bytes: &[u8], kind: OutputKind) -> io::Result<()> {
    let file = File::create(path)?;
    match kind {
        OutputKind::VcfText => {
            let mut out = io::BufWriter::new(file);
            out.write_all(bytes)
        }
        OutputKind::VcfGz => {
            let mut bgzf = htslib_rs::bgzf::io::Writer::new(file);
            bgzf.write_all(bytes)?;
            bgzf.finish().map(|_| ())
        }
        OutputKind::Bcf => write_bcf_from_vcf_text(bytes, file),
    }
}

fn maybe_index_prefix_vcf_gz(path: &Path, kind: OutputKind) -> io::Result<()> {
    match kind {
        OutputKind::VcfText => Ok(()),
        OutputKind::VcfGz => {
            let index = build_vcf_tbi_from_path(path)?;
            write_tbi(path_with_added_extension(path, "tbi"), &index)
        }
        OutputKind::Bcf => {
            let index = build_bcf_csi_with_min_shift(path, 14)?;
            write_csi(path_with_added_extension(path, "csi"), &index)
        }
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

fn path_with_added_extension(path: &Path, extension: &str) -> PathBuf {
    let mut out = path.as_os_str().to_os_string();
    out.push(".");
    out.push(extension);
    PathBuf::from(out)
}

fn evaluate(fields: &[String], args: &Args) -> io::Result<bool> {
    if let Some(expr) = &args.include_expr {
        let value = evaluate_expression(expr, fields)?;
        if !value.truthy() {
            return Ok(false);
        }
    }
    if let Some(expr) = &args.exclude_expr {
        let value = evaluate_expression(expr, fields)?;
        if value.truthy() {
            return Ok(false);
        }
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

fn record_in_regions(fields: &[String], regions: &[RegionSpec], targets: &[RegionSpec]) -> bool {
    let pos = fields[1].parse::<i64>().unwrap_or(0);
    (regions.is_empty() || matches_any(regions, &fields[0], pos))
        && (targets.is_empty() || matches_any(targets, &fields[0], pos))
}

fn matches_any(specs: &[RegionSpec], chrom: &str, pos: i64) -> bool {
    specs.iter().any(|spec| {
        spec.contig == chrom
            && spec.start.is_none_or(|start| pos >= start)
            && spec.end.is_none_or(|end| pos <= end)
    })
}

fn parse_region_list(out: &mut Vec<RegionSpec>, raw: &str) -> Result<(), ParseOutcome> {
    for token in raw.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        out.push(parse_region(token)?);
    }
    Ok(())
}

fn parse_region(token: &str) -> Result<RegionSpec, ParseOutcome> {
    let (contig, range) = match token.rsplit_once(':') {
        Some((chrom, rest)) if !chrom.is_empty() => (chrom.to_owned(), Some(rest)),
        _ => (token.to_owned(), None),
    };
    let (start, end) = match range {
        None => (None, None),
        Some(rest) => {
            let clean = rest.replace(',', "");
            if let Some((s, e)) = clean.split_once('-') {
                let start = s
                    .parse::<i64>()
                    .map_err(|_| ParseOutcome::Error(format!("invalid region '{token}'")))?;
                let end = e
                    .parse::<i64>()
                    .map_err(|_| ParseOutcome::Error(format!("invalid region '{token}'")))?;
                (Some(start), Some(end))
            } else {
                let pos = clean
                    .parse::<i64>()
                    .map_err(|_| ParseOutcome::Error(format!("invalid region '{token}'")))?;
                (Some(pos), Some(pos))
            }
        }
    };
    Ok(RegionSpec { contig, start, end })
}

fn load_region_file(out: &mut Vec<RegionSpec>, path: &str) -> Result<(), ParseOutcome> {
    let text = read_text_path(path)
        .map_err(|e| ParseOutcome::Error(format!("failed to read regions file '{path}': {e}")))?;
    for line in text.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = trimmed.split('\t').collect();
        if cols.len() >= 3 && cols[1].parse::<i64>().is_ok() && cols[2].parse::<i64>().is_ok() {
            let start = cols[1].parse::<i64>().unwrap();
            let end = cols[2].parse::<i64>().unwrap();
            out.push(RegionSpec {
                contig: cols[0].to_owned(),
                start: Some(start),
                end: Some(end),
            });
        } else {
            out.push(parse_region(trimmed)?);
        }
    }
    Ok(())
}

fn read_text_path(path: &str) -> io::Result<String> {
    let file = File::open(path)?;
    let mut text = String::new();
    if path.ends_with(".gz") || path.ends_with(".bgz") {
        MultiGzDecoder::new(file).read_to_string(&mut text)?;
    } else {
        BufReader::new(file).read_to_string(&mut text)?;
    }
    Ok(text)
}

fn parse_collapse(raw: &str) -> Result<CollapseMode, ParseOutcome> {
    match raw {
        "none" | "exact" => Ok(CollapseMode::None),
        "any" | "all" | "some" => Ok(CollapseMode::Any),
        "both" | "snps" | "indels" => Ok(CollapseMode::Both),
        "id" => Ok(CollapseMode::Id),
        _ => Err(ParseOutcome::Error(format!(
            "invalid collapse mode '{raw}'"
        ))),
    }
}

fn parse_nfiles(raw: &str) -> Result<NFiles, ParseOutcome> {
    let (op, digits) = match raw.as_bytes().first().copied() {
        Some(b'=') => (NOp::Eq, &raw[1..]),
        Some(b'+') => (NOp::Ge, &raw[1..]),
        Some(b'-') => (NOp::Le, &raw[1..]),
        Some(b'~') => {
            let mask = raw[1..]
                .chars()
                .map(|c| match c {
                    '0' => Ok(false),
                    '1' => Ok(true),
                    _ => Err(ParseOutcome::Error(format!("invalid -n bitmask '{raw}'"))),
                })
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(NFiles {
                op: NOp::Exact,
                count: mask.iter().filter(|&&present| present).count(),
                exact: Some(mask),
            });
        }
        _ => (NOp::Eq, raw),
    };
    let count = digits
        .parse::<usize>()
        .map_err(|_| ParseOutcome::Error(format!("invalid -n value '{raw}'")))?;
    Ok(NFiles {
        op,
        count,
        exact: None,
    })
}

fn parse_write_inputs(raw: &str) -> Result<Vec<usize>, ParseOutcome> {
    raw.split(',')
        .map(|part| {
            part.trim()
                .parse::<usize>()
                .map_err(|_| ParseOutcome::Error(format!("invalid -w index '{part}'")))
        })
        .collect()
}

fn parse_output_kind(raw: &str) -> Result<OutputKind, ParseOutcome> {
    match raw.chars().next() {
        Some('v') => Ok(OutputKind::VcfText),
        Some('z') | Some('0'..='9') => Ok(OutputKind::VcfGz),
        Some('b') | Some('u') => Ok(OutputKind::Bcf),
        _ => Err(ParseOutcome::Error(format!("invalid output type '{raw}'"))),
    }
}

fn next_string<'a, I>(
    iter: &mut std::iter::Peekable<I>,
    option: &str,
) -> Result<String, ParseOutcome>
where
    I: Iterator<Item = &'a OsString>,
{
    iter.next()
        .map(|s| s.to_string_lossy().into_owned())
        .ok_or_else(|| ParseOutcome::Error(format!("expected argument after {option}")))
}

fn value_after_equals(raw: &str) -> &str {
    raw.split_once('=').map(|(_, v)| v).unwrap_or("")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nfiles_parses_upstream_prefixes() {
        assert_eq!(parse_nfiles("=2").unwrap().op, NOp::Eq);
        assert_eq!(parse_nfiles("+1").unwrap().op, NOp::Ge);
        assert_eq!(parse_nfiles("-1").unwrap().op, NOp::Le);
        assert_eq!(
            parse_nfiles("~101").unwrap().exact,
            Some(vec![true, false, true])
        );
    }

    #[test]
    fn parse_region_handles_single_positions() {
        let region = parse_region("20:140").unwrap();
        assert_eq!(region.contig, "20");
        assert_eq!(region.start, Some(140));
        assert_eq!(region.end, Some(140));
    }
}
