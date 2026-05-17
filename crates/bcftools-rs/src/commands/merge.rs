//! Focused `bcftools merge` implementation (upstream `vcfmerge.c`).
//!
//! This local slice merges records that are present in every input or are
//! absent from some inputs and have identical site fields, plus a narrow
//! same-position allele-union slice. Full synced-reader merging, full allele
//! unification, full INFO rules, and gVCF mode remain tracked in `TODO.md`.

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
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
    -i, --info-rules TAG:METHOD,..  Apply AC:sum/AN:sum in the current text ALT-union slice\n\
    -m, --merge TYPE                Support narrow none/both/snp-ins-del slices; other modes accepted for command-shape compatibility\n\
    -o, --output FILE               Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9]  u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
        --force-single              Allow a single input for command-shape compatibility\n\
        --force-samples             Allow duplicate sample names by prefixing later inputs with the input index\n\
        --no-index                  Accepted for command-shape compatibility in this text slice\n\
        --no-version                Accepted for command-shape compatibility\n\
\n";

#[derive(Debug)]
struct Args {
    inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    force_samples: bool,
    info_rules: InfoRules,
    merge_mode: MergeMode,
}

#[derive(Debug, Clone, Copy, Default)]
struct InfoRules {
    sum_ac: bool,
    sum_an: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MergeMode {
    Default,
    None,
    Both,
    SnpInsDel,
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
struct MergedSite {
    fixed: Vec<String>,
    samples_by_input: Vec<Option<Vec<String>>>,
    order: usize,
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
    let mut force_single = false;
    let mut force_samples = false;
    let mut info_rules = InfoRules::default();
    let mut merge_mode = MergeMode::Default;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "-l" | "--file-list" => {
                file_list = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-i" | "--info-rules" => {
                info_rules = parse_info_rules(&next_string(&mut iter, raw.as_ref())?);
            }
            "-m" | "--merge" => {
                merge_mode = parse_merge_mode(&next_string(&mut iter, raw.as_ref())?);
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--force-single" => force_single = true,
            "--force-samples" => force_samples = true,
            "--no-index" => {}
            "--no-version" => {}
            _ if raw.starts_with("--file-list=") => {
                file_list = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--info-rules=") => {
                info_rules = parse_info_rules(value_after_equals(&raw))
            }
            _ if raw.starts_with("--merge=") => {
                merge_mode = parse_merge_mode(value_after_equals(&raw))
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-l") && raw.len() > 2 => {
                file_list = Some(PathBuf::from(&raw[2..]))
            }
            _ if raw.starts_with("-i") && raw.len() > 2 => info_rules = parse_info_rules(&raw[2..]),
            _ if raw.starts_with("-m") && raw.len() > 2 => merge_mode = parse_merge_mode(&raw[2..]),
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
    if inputs.len() < 2 && !force_single {
        return Err(ParseOutcome::Error(
            "expected at least two input VCF/BCF paths".into(),
        ));
    }
    if inputs.is_empty() {
        return Err(ParseOutcome::Error("expected an input VCF/BCF path".into()));
    }

    Ok(Args {
        inputs,
        output,
        output_kind,
        force_samples,
        info_rules,
        merge_mode,
    })
}

fn parse_merge_mode(raw: &str) -> MergeMode {
    match raw {
        "none" => MergeMode::None,
        "both" => MergeMode::Both,
        "snp-ins-del" => MergeMode::SnpInsDel,
        _ => MergeMode::Default,
    }
}

fn parse_info_rules(raw: &str) -> InfoRules {
    let mut rules = InfoRules::default();
    for rule in raw.split(',') {
        let Some((tag, method)) = rule.split_once(':') else {
            continue;
        };
        if method != "sum" {
            continue;
        }
        match tag {
            "AC" => rules.sum_ac = true,
            "AN" => rules.sum_an = true,
            _ => {}
        }
    }
    rules
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
    let merged = merge_inputs(
        &inputs,
        args.force_samples,
        args.info_rules,
        args.merge_mode,
    )?;
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

fn merge_inputs(
    inputs: &[VcfInput],
    force_samples: bool,
    info_rules: InfoRules,
    merge_mode: MergeMode,
) -> io::Result<String> {
    let first = inputs
        .first()
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "no inputs"))?;
    let mut sample_names = Vec::new();
    let mut seen_samples = HashSet::new();
    for (input_idx, input) in inputs.iter().enumerate() {
        for sample in &input.samples {
            let mut name = sample.clone();
            if seen_samples.contains(&name) {
                if !force_samples {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("duplicate sample name '{sample}'"),
                    ));
                }
                while seen_samples.contains(&name) {
                    name = format!("{}:{name}", input_idx + 1);
                }
            }
            seen_samples.insert(name.clone());
            sample_names.push(name);
        }
    }

    let mut out = render_meta_with_pass_filter(&first.meta);
    out.push_str(&first.fixed_header.join("\t"));
    if !sample_names.is_empty() {
        out.push('\t');
        out.push_str(&sample_names.join("\t"));
    }
    out.push('\n');

    for input in &inputs[1..] {
        if !fixed_headers_compatible(&first.fixed_header, &input.fixed_header) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "inputs must have compatible fixed VCF columns",
            ));
        }
    }

    let mut sites = collect_sites(inputs, info_rules, merge_mode)?;
    let contigs = contig_order(&first.meta);
    sites.sort_by(|a, b| compare_sites(a, b, &contigs, merge_mode));

    for site in sites {
        let mut samples = Vec::new();
        for (input_idx, input) in inputs.iter().enumerate() {
            match &site.samples_by_input[input_idx] {
                Some(values) => samples.extend(values.iter().cloned()),
                None => {
                    let missing = missing_sample_value(&site.fixed);
                    samples.extend(std::iter::repeat_n(missing, input.samples.len()));
                }
            }
        }
        out.push_str(&site.fixed.join("\t"));
        if !samples.is_empty() {
            out.push('\t');
            out.push_str(&samples.join("\t"));
        }
        out.push('\n');
    }

    Ok(out)
}

fn fixed_headers_compatible(first: &[String], other: &[String]) -> bool {
    if other == first {
        return true;
    }
    first.len() == 9 && other.len() == 8 && first[8] == "FORMAT" && first[..8] == other[..8]
}

fn render_meta_with_pass_filter(meta: &[String]) -> String {
    let has_pass = meta
        .iter()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"));
    let mut out = String::new();
    let mut inserted = false;
    for line in meta {
        out.push_str(line);
        out.push('\n');
        if !inserted && !has_pass && line.starts_with("##fileformat=") {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            inserted = true;
        }
    }
    if !inserted && !has_pass {
        out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    }
    out
}

fn collect_sites(
    inputs: &[VcfInput],
    info_rules: InfoRules,
    merge_mode: MergeMode,
) -> io::Result<Vec<MergedSite>> {
    let mut sites: Vec<MergedSite> = Vec::new();
    let mut by_key = HashMap::new();
    let mut by_locus: HashMap<String, Vec<usize>> = HashMap::new();
    let mut by_position: HashMap<String, Vec<usize>> = HashMap::new();

    for (input_idx, input) in inputs.iter().enumerate() {
        for record in &input.records {
            let key = site_key(record);
            if let Some(site_idx) = by_key.get(&key).copied() {
                let site: &mut MergedSite = &mut sites[site_idx];
                if record.fixed != site.fixed {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "conflicting records at {}:{} require full merge semantics",
                            record.fixed[0], record.fixed[1]
                        ),
                    ));
                }
                if site.samples_by_input[input_idx].is_some() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!(
                            "duplicate record at {}:{}",
                            record.fixed[0], record.fixed[1]
                        ),
                    ));
                }
                site.samples_by_input[input_idx] = Some(record.samples.clone());
            } else {
                let mut merged = false;
                let mut same_locus_conflict = None;
                if merge_mode != MergeMode::None {
                    if let Some(site_indices) = by_locus.get(&locus_key(record)).cloned() {
                        for site_idx in site_indices {
                            let site: &mut MergedSite = &mut sites[site_idx];
                            if can_merge_same_locus_alt_union(site, record) {
                                merge_sites_only_alt_union(site, record, info_rules);
                                site.samples_by_input[input_idx] = Some(record.samples.clone());
                                merged = true;
                                break;
                            }
                            same_locus_conflict = Some(site_idx);
                        }
                    }
                    if !merged
                        && supports_sampled_same_position_union(merge_mode)
                        && let Some(site_indices) = by_position.get(&position_key(record)).cloned()
                    {
                        for site_idx in site_indices {
                            let site: &mut MergedSite = &mut sites[site_idx];
                            if can_merge_sampled_same_position(site, record, merge_mode) {
                                merge_sampled_same_position(site, record, input_idx)?;
                                merged = true;
                                break;
                            }
                        }
                    }
                    if !merged
                        && merge_mode == MergeMode::Default
                        && let Some(site_idx) = same_locus_conflict
                    {
                        let existing = &sites[site_idx].fixed;
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidInput,
                            format!(
                                "conflicting records at {}:{} require full merge semantics",
                                existing[0], existing[1]
                            ),
                        ));
                    }
                }
                if !merged {
                    let site_idx = sites.len();
                    let mut samples_by_input = vec![None; inputs.len()];
                    samples_by_input[input_idx] = Some(record.samples.clone());
                    by_key.insert(key, site_idx);
                    by_locus
                        .entry(locus_key(record))
                        .or_default()
                        .push(site_idx);
                    by_position
                        .entry(position_key(record))
                        .or_default()
                        .push(site_idx);
                    sites.push(MergedSite {
                        fixed: record.fixed.clone(),
                        samples_by_input,
                        order: site_idx,
                    });
                }
            }
        }
    }

    Ok(sites)
}

fn site_key(record: &RecordLine) -> String {
    record
        .fixed
        .iter()
        .take(5)
        .cloned()
        .collect::<Vec<_>>()
        .join("\t")
}

fn locus_key(record: &RecordLine) -> String {
    record
        .fixed
        .iter()
        .take(4)
        .cloned()
        .collect::<Vec<_>>()
        .join("\t")
}

fn position_key(record: &RecordLine) -> String {
    record
        .fixed
        .iter()
        .take(3)
        .cloned()
        .collect::<Vec<_>>()
        .join("\t")
}

fn can_merge_same_locus_alt_union(site: &MergedSite, record: &RecordLine) -> bool {
    (site.fixed.len() == 8 || site.fixed.len() == 9)
        && record.fixed.len() == 8
        && record.samples.is_empty()
        && site.fixed[..4] == record.fixed[..4]
}

fn supports_sampled_same_position_union(merge_mode: MergeMode) -> bool {
    matches!(merge_mode, MergeMode::Both | MergeMode::SnpInsDel)
}

fn can_merge_sampled_same_position(
    site: &MergedSite,
    record: &RecordLine,
    merge_mode: MergeMode,
) -> bool {
    if site.fixed.len() < 9
        || record.fixed.len() < 9
        || site.fixed[..3] != record.fixed[..3]
        || site.fixed[8] != record.fixed[8]
        || merged_ref(&site.fixed[3], &record.fixed[3]).is_none()
    {
        return false;
    }

    match merge_mode {
        MergeMode::Both => {
            let site_class = coarse_variant_class(&site.fixed[3], &site.fixed[4]);
            site_class != CoarseVariantClass::Other
                && site_class == coarse_variant_class(&record.fixed[3], &record.fixed[4])
        }
        MergeMode::SnpInsDel => {
            let site_class = precise_variant_class(&site.fixed[3], &site.fixed[4]);
            site_class != PreciseVariantClass::Other
                && site_class == precise_variant_class(&record.fixed[3], &record.fixed[4])
        }
        _ => false,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CoarseVariantClass {
    Snp,
    Indel,
    Other,
}

fn coarse_variant_class(reference: &str, alt: &str) -> CoarseVariantClass {
    let classes = split_alt(alt)
        .into_iter()
        .map(|alt| precise_allele_class(reference, &alt))
        .collect::<Vec<_>>();
    if classes.is_empty() {
        return CoarseVariantClass::Other;
    }
    if classes
        .iter()
        .all(|class| *class == PreciseVariantClass::Snp)
    {
        CoarseVariantClass::Snp
    } else if classes.iter().all(|class| {
        matches!(
            class,
            PreciseVariantClass::Insertion | PreciseVariantClass::Deletion
        )
    }) {
        CoarseVariantClass::Indel
    } else {
        CoarseVariantClass::Other
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PreciseVariantClass {
    Snp,
    Insertion,
    Deletion,
    Other,
}

fn precise_variant_class(reference: &str, alt: &str) -> PreciseVariantClass {
    let classes = split_alt(alt)
        .into_iter()
        .map(|alt| precise_allele_class(reference, &alt))
        .collect::<Vec<_>>();
    let Some(first) = classes.first().copied() else {
        return PreciseVariantClass::Other;
    };
    if classes.iter().all(|class| *class == first) {
        first
    } else {
        PreciseVariantClass::Other
    }
}

fn precise_allele_class(reference: &str, alt: &str) -> PreciseVariantClass {
    if reference.len() == alt.len() && reference.len() == 1 {
        PreciseVariantClass::Snp
    } else if alt.len() > reference.len() && alt.starts_with(reference) {
        PreciseVariantClass::Insertion
    } else if reference.len() > alt.len() && reference.starts_with(alt) {
        PreciseVariantClass::Deletion
    } else {
        PreciseVariantClass::Other
    }
}

fn merge_sampled_same_position(
    site: &mut MergedSite,
    record: &RecordLine,
    input_idx: usize,
) -> io::Result<()> {
    if site.samples_by_input[input_idx].is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "duplicate record at {}:{}",
                record.fixed[0], record.fixed[1]
            ),
        ));
    }

    let old_ref = site.fixed[3].clone();
    let old_alts = normalize_alts(&old_ref, &site.fixed[4], &old_ref);
    let new_ref = merged_ref(&old_ref, &record.fixed[3]).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "conflicting records at {}:{} require full merge semantics",
                record.fixed[0], record.fixed[1]
            ),
        )
    })?;
    let normalized_site_alts = normalize_alts(&old_ref, &site.fixed[4], &new_ref);
    let normalized_record_alts = normalize_alts(&record.fixed[3], &record.fixed[4], &new_ref);

    let mut merged_alts = normalized_site_alts.clone();
    for alt in &normalized_record_alts {
        if !merged_alts.contains(alt) {
            merged_alts.push(alt.clone());
        }
    }

    let old_site_map = allele_map(&normalized_site_alts, &merged_alts);
    let record_map = allele_map(&normalized_record_alts, &merged_alts);
    if new_ref != old_ref || merged_alts != old_alts {
        let format = site.fixed[8].clone();
        for values in site.samples_by_input.iter_mut().flatten() {
            *values = remap_sample_values(&format, values, &old_site_map, merged_alts.len());
        }
    }

    site.fixed[3] = new_ref;
    site.fixed[4] = if merged_alts.is_empty() {
        ".".to_owned()
    } else {
        merged_alts.join(",")
    };
    site.samples_by_input[input_idx] = Some(remap_sample_values(
        &record.fixed[8],
        &record.samples,
        &record_map,
        merged_alts.len(),
    ));
    Ok(())
}

fn merged_ref(a: &str, b: &str) -> Option<String> {
    if a.len() >= b.len() && a.starts_with(b) {
        Some(a.to_owned())
    } else if b.starts_with(a) {
        Some(b.to_owned())
    } else {
        None
    }
}

fn normalize_alts(reference: &str, alt: &str, merged_ref: &str) -> Vec<String> {
    let suffix = merged_ref.strip_prefix(reference).unwrap_or("");
    split_alt(alt)
        .into_iter()
        .map(|alt| format!("{alt}{suffix}"))
        .collect()
}

fn allele_map(old_alts: &[String], merged_alts: &[String]) -> Vec<Option<usize>> {
    let mut map = Vec::with_capacity(old_alts.len() + 1);
    map.push(Some(0));
    for alt in old_alts {
        map.push(
            merged_alts
                .iter()
                .position(|merged| merged == alt)
                .map(|idx| idx + 1),
        );
    }
    map
}

fn remap_sample_values(
    format: &str,
    samples: &[String],
    allele_map: &[Option<usize>],
    alt_count: usize,
) -> Vec<String> {
    let keys = format.split(':').collect::<Vec<_>>();
    samples
        .iter()
        .map(|sample| {
            let mut values = sample.split(':').map(str::to_owned).collect::<Vec<_>>();
            for (idx, key) in keys.iter().enumerate() {
                let Some(value) = values.get_mut(idx) else {
                    continue;
                };
                match *key {
                    "GT" => *value = remap_gt(value, allele_map),
                    "AD" => *value = remap_ad(value, allele_map, alt_count),
                    _ => {}
                }
            }
            values.join(":")
        })
        .collect()
}

fn remap_gt(raw: &str, allele_map: &[Option<usize>]) -> String {
    let mut out = String::new();
    let mut chars = raw.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch.is_ascii_digit() {
            let mut allele = ch.to_string();
            while let Some(next) = chars.peek().copied() {
                if next.is_ascii_digit() {
                    allele.push(next);
                    chars.next();
                } else {
                    break;
                }
            }
            let mapped = allele
                .parse::<usize>()
                .ok()
                .and_then(|idx| allele_map.get(idx).copied().flatten());
            match mapped {
                Some(idx) => out.push_str(&idx.to_string()),
                None => out.push('.'),
            }
        } else {
            out.push(ch);
        }
    }
    out
}

fn remap_ad(raw: &str, allele_map: &[Option<usize>], alt_count: usize) -> String {
    if raw == "." {
        return raw.to_owned();
    }
    let old_values = raw.split(',').collect::<Vec<_>>();
    let mut values = vec![".".to_owned(); alt_count + 1];
    for (old_idx, value) in old_values.iter().enumerate() {
        if let Some(new_idx) = allele_map.get(old_idx).copied().flatten()
            && let Some(slot) = values.get_mut(new_idx)
        {
            *slot = (*value).to_owned();
        }
    }
    values.join(",")
}

fn merge_sites_only_alt_union(site: &mut MergedSite, record: &RecordLine, info_rules: InfoRules) {
    let mut alts = split_alt(&site.fixed[4]);
    let next_alts = split_alt(&record.fixed[4]);
    for alt in &next_alts {
        if !alts.contains(alt) {
            alts.push(alt.clone());
        }
    }

    let site_has_no_sample_values = site_has_no_sample_values(site);
    let mut ac_by_alt = HashMap::new();
    add_ac_by_alt(&mut ac_by_alt, &site.fixed[4], &site.fixed[7]);
    if site_has_no_sample_values || info_rules.sum_ac {
        add_ac_by_alt(&mut ac_by_alt, &record.fixed[4], &record.fixed[7]);
    }
    let an = info_i64(&site.fixed[7], "AN").unwrap_or(0)
        + if site_has_no_sample_values || info_rules.sum_an {
            info_i64(&record.fixed[7], "AN").unwrap_or(0)
        } else {
            0
        };
    let ac_values = alts
        .iter()
        .map(|alt| ac_by_alt.get(alt).copied().unwrap_or(0).to_string())
        .collect::<Vec<_>>()
        .join(",");

    site.fixed[4] = alts.join(",");
    site.fixed[7] = if site_has_no_sample_values || info_rules.sum_ac || info_rules.sum_an {
        format!("AC={ac_values};AN={an}")
    } else {
        format!("AN={an};AC={ac_values}")
    };
}

fn site_has_no_sample_values(site: &MergedSite) -> bool {
    site.samples_by_input
        .iter()
        .all(|samples| samples.as_ref().is_none_or(|values| values.is_empty()))
}

fn split_alt(raw: &str) -> Vec<String> {
    if raw == "." || raw.is_empty() {
        Vec::new()
    } else {
        raw.split(',').map(str::to_owned).collect()
    }
}

fn add_ac_by_alt(out: &mut HashMap<String, i64>, alt: &str, info: &str) {
    let alts = split_alt(alt);
    let Some(ac_raw) = info_value(info, "AC") else {
        return;
    };
    for (alt, value) in alts.iter().zip(ac_raw.split(',')) {
        let Ok(value) = value.parse::<i64>() else {
            continue;
        };
        *out.entry(alt.clone()).or_insert(0) += value;
    }
}

fn info_i64(info: &str, key: &str) -> Option<i64> {
    info_value(info, key)?.parse().ok()
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|field| {
        let (name, value) = field.split_once('=')?;
        (name == key).then_some(value)
    })
}

fn contig_order(meta: &[String]) -> HashMap<String, usize> {
    let mut order = HashMap::new();
    for line in meta {
        if let Some(rest) = line.strip_prefix("##contig=<ID=") {
            let id = rest.split([',', '>']).next().unwrap_or("").to_owned();
            if !id.is_empty() && !order.contains_key(&id) {
                order.insert(id, order.len());
            }
        }
    }
    order
}

fn compare_sites(
    a: &MergedSite,
    b: &MergedSite,
    contigs: &HashMap<String, usize>,
    merge_mode: MergeMode,
) -> Ordering {
    let a_chrom = a.fixed.first().map(String::as_str).unwrap_or("");
    let b_chrom = b.fixed.first().map(String::as_str).unwrap_or("");
    match (contigs.get(a_chrom).copied(), contigs.get(b_chrom).copied()) {
        (Some(a_idx), Some(b_idx)) => a_idx.cmp(&b_idx),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => a_chrom.cmp(b_chrom),
    }
    .then_with(|| {
        let a_pos = a
            .fixed
            .get(1)
            .and_then(|pos| pos.parse::<u64>().ok())
            .unwrap_or(0);
        let b_pos = b
            .fixed
            .get(1)
            .and_then(|pos| pos.parse::<u64>().ok())
            .unwrap_or(0);
        a_pos.cmp(&b_pos)
    })
    .then_with(|| {
        if merge_mode != MergeMode::Default {
            a.order.cmp(&b.order)
        } else {
            Ordering::Equal
        }
    })
    .then_with(|| a.fixed.get(2).cmp(&b.fixed.get(2)))
    .then_with(|| a.fixed.get(3).cmp(&b.fixed.get(3)))
    .then_with(|| a.fixed.get(4).cmp(&b.fixed.get(4)))
    .then_with(|| a.order.cmp(&b.order))
}

fn missing_sample_value(fixed: &[String]) -> String {
    let Some(format) = fixed.get(8) else {
        return ".".to_owned();
    };
    if format == "." || format.is_empty() {
        return ".".to_owned();
    }
    format
        .split(':')
        .map(|key| if key == "GT" { "./." } else { "." })
        .collect::<Vec<_>>()
        .join(":")
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
        let merged =
            merge_inputs(&[a, b], false, InfoRules::default(), MergeMode::Default).unwrap();
        assert!(merged.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"));
        assert!(merged.contains("1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1"));
    }
}
