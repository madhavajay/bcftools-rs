//! Focused `bcftools annotate` implementation (upstream `vcfannotate.c`).
//!
//! Local slice supports chromosome renaming via `--rename-chrs` and tag
//! removal via `-x`/`--remove` for `ID`, `QUAL`, `FILTER`, `FILTER/<ID>`,
//! `INFO`, `INFO/<ID>`, `FORMAT`, and `FORMAT/<ID>` targets, including the
//! upstream `^` keep-only form for FILTER/INFO/FORMAT IDs. Full annotation
//! transfer and header injection remain tracked in `TODO.md`.

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
About:   Annotate and edit VCF/BCF files.\n\
Usage:   bcftools annotate [OPTIONS] <in.vcf.gz>\n\
\n\
Options:\n\
        --rename-chrs FILE           Rename chromosomes according to a two-column map\n\
    -x, --remove LIST                Remove annotations (LIST = comma-separated ID, QUAL, FILTER[/<ID>], INFO[/<ID>], FORMAT[/<ID>])\n\
    -o, --output FILE                Write output to a file [standard output]\n\
    -O, --output-type u|b|v|z[0-9]   u/b: BCF, v/z: VCF/BGZF VCF [v]\n\
        --force                      Accepted for command-shape compatibility\n\
        --no-version                 Accepted for command-shape compatibility\n\
\n";

#[derive(Debug, Default)]
struct RemoveSet {
    id: bool,
    qual: bool,
    all_filter: bool,
    all_info: bool,
    all_format: bool,
    filter_ids: HashSet<String>,
    info_ids: HashSet<String>,
    format_ids: HashSet<String>,
    keep_filter_ids: Option<HashSet<String>>,
    keep_info_ids: Option<HashSet<String>>,
    keep_format_ids: Option<HashSet<String>>,
}

impl RemoveSet {
    fn is_empty(&self) -> bool {
        !self.id
            && !self.qual
            && !self.all_filter
            && !self.all_info
            && !self.all_format
            && self.filter_ids.is_empty()
            && self.info_ids.is_empty()
            && self.format_ids.is_empty()
            && self.keep_filter_ids.is_none()
            && self.keep_info_ids.is_none()
            && self.keep_format_ids.is_none()
    }
}

#[derive(Debug)]
struct Args {
    input: PathBuf,
    rename_chrs: Option<PathBuf>,
    remove: RemoveSet,
    output: Option<PathBuf>,
    output_kind: OutputKind,
    no_version: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OutputKind {
    VcfText,
    VcfGz,
    Bcf,
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
                eprintln!("{}", fmt_etag("main_vcfannotate", &format!("{e}")));
                ExitCode::FAILURE
            }
        },
        Err(ParseOutcome::Usage) => {
            eprint!("{USAGE}");
            ExitCode::FAILURE
        }
        Err(ParseOutcome::Error(message)) => {
            eprintln!("{}", fmt_etag("main_vcfannotate", &message));
            ExitCode::FAILURE
        }
    }
}

fn parse_args(argv: &[OsString]) -> Result<Args, ParseOutcome> {
    let mut input = None;
    let mut rename_chrs = None;
    let mut remove = RemoveSet::default();
    let mut output = None;
    let mut output_kind = OutputKind::VcfText;
    let mut no_version = false;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-h" | "--help" | "-?" => return Err(ParseOutcome::Usage),
            "--rename-chrs" => {
                rename_chrs = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-x" | "--remove" => {
                parse_remove_list(&next_string(&mut iter, raw.as_ref())?, &mut remove)?
            }
            "-o" | "--output" => {
                output = Some(PathBuf::from(next_string(&mut iter, raw.as_ref())?))
            }
            "-O" | "--output-type" => {
                output_kind = parse_output_kind(&next_string(&mut iter, raw.as_ref())?)?
            }
            "--force" => {}
            "--no-version" => no_version = true,
            _ if raw.starts_with("--rename-chrs=") => {
                rename_chrs = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--remove=") => {
                parse_remove_list(value_after_equals(&raw), &mut remove)?
            }
            _ if raw.starts_with("--output=") => {
                output = Some(PathBuf::from(value_after_equals(&raw)))
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_output_kind(value_after_equals(&raw))?
            }
            _ if raw.starts_with("-x") && raw.len() > 2 => {
                parse_remove_list(&raw[2..], &mut remove)?
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => output = Some(PathBuf::from(&raw[2..])),
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_output_kind(&raw[2..])?
            }
            _ if raw.starts_with('-') => {
                return Err(ParseOutcome::Error(format!(
                    "unrecognized option '{raw}' in this local annotate slice"
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
    if rename_chrs.is_none() && remove.is_empty() {
        return Err(ParseOutcome::Error(
            "expected --rename-chrs FILE and/or -x LIST".into(),
        ));
    }
    Ok(Args {
        input,
        rename_chrs,
        remove,
        output,
        output_kind,
        no_version,
    })
}

fn parse_remove_list(list: &str, remove: &mut RemoveSet) -> Result<(), ParseOutcome> {
    for entry in list.split(',') {
        let entry = entry.trim();
        if entry.is_empty() {
            continue;
        }
        let keep_only = entry.starts_with('^');
        let entry = entry.strip_prefix('^').unwrap_or(entry);
        match entry {
            "ID" => remove.id = true,
            "QUAL" => remove.qual = true,
            "FILTER" => remove.all_filter = true,
            "INFO" => remove.all_info = true,
            "FORMAT" | "FMT" => remove.all_format = true,
            _ => {
                if let Some(id) = entry.strip_prefix("FILTER/") {
                    if id.is_empty() {
                        return Err(ParseOutcome::Error(format!("invalid -x entry '{entry}'")));
                    }
                    if keep_only || remove.keep_filter_ids.is_some() {
                        remove
                            .keep_filter_ids
                            .get_or_insert_with(HashSet::new)
                            .insert(id.to_owned());
                    } else {
                        remove.filter_ids.insert(id.to_owned());
                    }
                } else if let Some(id) = entry.strip_prefix("INFO/") {
                    if id.is_empty() {
                        return Err(ParseOutcome::Error(format!("invalid -x entry '{entry}'")));
                    }
                    if keep_only || remove.keep_info_ids.is_some() {
                        remove
                            .keep_info_ids
                            .get_or_insert_with(HashSet::new)
                            .insert(id.to_owned());
                    } else {
                        remove.info_ids.insert(id.to_owned());
                    }
                } else if let Some(id) = format_id_entry(entry) {
                    if id.is_empty() {
                        return Err(ParseOutcome::Error(format!("invalid -x entry '{entry}'")));
                    }
                    if keep_only || remove.keep_format_ids.is_some() {
                        remove
                            .keep_format_ids
                            .get_or_insert_with(HashSet::new)
                            .insert(id.to_owned());
                    } else {
                        remove.format_ids.insert(id.to_owned());
                    }
                } else {
                    return Err(ParseOutcome::Error(format!(
                        "unrecognized -x entry '{entry}'"
                    )));
                }
            }
        }
    }
    Ok(())
}

fn format_id_entry(entry: &str) -> Option<&str> {
    entry
        .strip_prefix("FORMAT/")
        .or_else(|| entry.strip_prefix("FMT/"))
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
    let rename = match &args.rename_chrs {
        Some(path) => Some(read_rename_map(path)?),
        None => None,
    };
    let mut text = read_vcf_text(&args.input)?;
    transform_vcf_text(&mut text, rename.as_ref(), &args.remove, args.no_version);
    write_output(text.as_bytes(), args)
}

fn read_rename_map(path: &Path) -> io::Result<HashMap<String, String>> {
    let mut map = HashMap::new();
    for (i, line) in fs::read_to_string(path)?.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() < 2 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid rename-chrs line {}: {line}", i + 1),
            ));
        }
        map.insert(fields[0].to_owned(), fields[1].to_owned());
    }
    Ok(map)
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
        ".bcftools-rs-annotate-{}-{nanos}.tmp",
        std::process::id()
    ))
}

fn transform_vcf_text(
    text: &mut String,
    rename: Option<&HashMap<String, String>>,
    remove: &RemoveSet,
    no_version: bool,
) {
    let empty = HashMap::new();
    let rename = rename.unwrap_or(&empty);
    let inject_pass_header = remove.all_filter
        && !text
            .lines()
            .any(|line| header_id_for_kind(line, "##FILTER=<") == Some("PASS"));
    let mut out = String::with_capacity(text.len());
    for line in text.lines() {
        if no_version && line.starts_with("##bcftools_") {
            continue;
        } else if line.starts_with("##fileformat=") {
            out.push_str(line);
            if inject_pass_header {
                out.push_str("\n##FILTER=<ID=PASS,Description=\"All filters passed\">");
            }
        } else if line.starts_with("##contig=<") {
            if !rename.is_empty() {
                out.push_str(&rename_contig_header(line, rename));
            } else {
                out.push_str(line);
            }
        } else if let Some(id) = header_id_for_kind(line, "##INFO=<") {
            if remove.all_info
                || remove.info_ids.contains(id)
                || remove
                    .keep_info_ids
                    .as_ref()
                    .is_some_and(|keep| !keep.contains(id))
            {
                continue;
            }
            out.push_str(line);
        } else if let Some(id) = header_id_for_kind(line, "##FILTER=<") {
            if id != "PASS"
                && (remove.all_filter
                    || remove.filter_ids.contains(id)
                    || remove
                        .keep_filter_ids
                        .as_ref()
                        .is_some_and(|keep| !keep.contains(id)))
            {
                continue;
            }
            out.push_str(line);
        } else if let Some(id) = header_id_for_kind(line, "##FORMAT=<") {
            if remove.format_ids.contains(id)
                || remove
                    .keep_format_ids
                    .as_ref()
                    .is_some_and(|keep| !keep.contains(id))
                || (remove.all_format && id != "GT")
            {
                continue;
            }
            out.push_str(line);
        } else if line.starts_with('#') {
            out.push_str(line);
        } else {
            let renamed;
            let mut record: &str = line;
            if !rename.is_empty() {
                renamed = rename_record_chrom(line, rename);
                record = &renamed;
            }
            if remove.is_empty() {
                out.push_str(record);
            } else {
                out.push_str(&apply_record_removals(record, remove));
            }
        }
        out.push('\n');
    }
    *text = out;
}

fn rename_contig_header(line: &str, map: &HashMap<String, String>) -> String {
    let Some(id_start) = line.find("ID=") else {
        return line.to_owned();
    };
    let value_start = id_start + 3;
    let value_end = line[value_start..]
        .find([',', '>'])
        .map(|idx| value_start + idx)
        .unwrap_or(line.len());
    let old = &line[value_start..value_end];
    let Some(new) = map.get(old) else {
        return line.to_owned();
    };
    let mut renamed = String::with_capacity(line.len() - old.len() + new.len());
    renamed.push_str(&line[..value_start]);
    renamed.push_str(new);
    renamed.push_str(&line[value_end..]);
    renamed
}

fn rename_record_chrom(line: &str, map: &HashMap<String, String>) -> String {
    let Some(tab) = line.find('\t') else {
        return line.to_owned();
    };
    let old = &line[..tab];
    let Some(new) = map.get(old) else {
        return line.to_owned();
    };
    let mut renamed = String::with_capacity(line.len() - old.len() + new.len());
    renamed.push_str(new);
    renamed.push_str(&line[tab..]);
    renamed
}

fn header_id_for_kind<'a>(line: &'a str, prefix: &str) -> Option<&'a str> {
    let rest = line.strip_prefix(prefix)?;
    let id_start = rest.find("ID=")? + 3;
    let id_rest = &rest[id_start..];
    let id_end = id_rest.find([',', '>']).unwrap_or(id_rest.len());
    Some(&id_rest[..id_end])
}

fn apply_record_removals(line: &str, remove: &RemoveSet) -> String {
    let mut fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if fields.len() < 8 {
        return line.to_owned();
    }
    if remove.id {
        fields[2] = ".".to_owned();
    }
    if remove.qual {
        fields[5] = ".".to_owned();
    }
    if remove.all_filter {
        fields[6] = ".".to_owned();
    } else if let Some(keep) = &remove.keep_filter_ids {
        fields[6] = filter_after_keep(&fields[6], keep);
    } else if !remove.filter_ids.is_empty() {
        fields[6] = filter_after_removal(&fields[6], &remove.filter_ids);
    }
    if remove.all_info {
        fields[7] = ".".to_owned();
    } else if let Some(keep) = &remove.keep_info_ids {
        fields[7] = info_after_keep(&fields[7], keep);
    } else if !remove.info_ids.is_empty() {
        fields[7] = info_after_removal(&fields[7], &remove.info_ids);
    }
    if fields.len() > 8 {
        apply_format_removals(&mut fields, remove);
    }
    fields.join("\t")
}

fn filter_after_removal(current: &str, drop: &HashSet<String>) -> String {
    if current == "." || current == "PASS" {
        return current.to_owned();
    }
    let kept: Vec<&str> = current
        .split(';')
        .filter(|tag| !drop.contains(*tag))
        .collect();
    if kept.is_empty() {
        "PASS".to_owned()
    } else {
        kept.join(";")
    }
}

fn filter_after_keep(current: &str, keep: &HashSet<String>) -> String {
    if current == "." || current == "PASS" {
        return ".".to_owned();
    }
    let kept: Vec<&str> = current
        .split(';')
        .filter(|tag| keep.contains(*tag))
        .collect();
    if kept.is_empty() {
        ".".to_owned()
    } else {
        kept.join(";")
    }
}

fn info_after_removal(current: &str, drop: &HashSet<String>) -> String {
    if current == "." {
        return current.to_owned();
    }
    let kept: Vec<&str> = current
        .split(';')
        .filter(|entry| {
            let key = entry.split_once('=').map(|(k, _)| k).unwrap_or(entry);
            !drop.contains(key)
        })
        .collect();
    if kept.is_empty() {
        ".".to_owned()
    } else {
        kept.join(";")
    }
}

fn info_after_keep(current: &str, keep: &HashSet<String>) -> String {
    if current == "." {
        return current.to_owned();
    }
    let kept: Vec<&str> = current
        .split(';')
        .filter(|entry| {
            let key = entry.split_once('=').map(|(k, _)| k).unwrap_or(entry);
            keep.contains(key)
        })
        .collect();
    if kept.is_empty() {
        ".".to_owned()
    } else {
        kept.join(";")
    }
}

fn apply_format_removals(fields: &mut [String], remove: &RemoveSet) {
    if !remove.all_format && remove.format_ids.is_empty() && remove.keep_format_ids.is_none() {
        return;
    }
    if fields[8] == "." {
        return;
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let kept_indexes: Vec<usize> = format_keys
        .iter()
        .enumerate()
        .filter_map(|(idx, key)| {
            let keep = if let Some(keep) = &remove.keep_format_ids {
                keep.contains(*key)
            } else if remove.all_format {
                *key == "GT"
            } else {
                !remove.format_ids.contains(*key)
            };
            keep.then_some(idx)
        })
        .collect();

    if kept_indexes.is_empty() {
        fields[8] = ".".to_owned();
        for sample in &mut fields[9..] {
            *sample = ".".to_owned();
        }
        return;
    }

    fields[8] = kept_indexes
        .iter()
        .map(|idx| format_keys[*idx])
        .collect::<Vec<_>>()
        .join(":");
    for sample in &mut fields[9..] {
        if *sample == "." {
            continue;
        }
        let values: Vec<&str> = sample.split(':').collect();
        *sample = kept_indexes
            .iter()
            .map(|idx| values.get(*idx).copied().unwrap_or("."))
            .collect::<Vec<_>>()
            .join(":");
    }
}

fn write_output(bytes: &[u8], args: &Args) -> io::Result<()> {
    match &args.output {
        Some(path) if path != Path::new("-") => {
            let file = File::create(path)?;
            write_to(bytes, args.output_kind, file)
        }
        _ => {
            let stdout = io::stdout();
            let out = stdout.lock();
            write_to(bytes, args.output_kind, out)
        }
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

    fn sample_vcf() -> String {
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=10>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n\
##FILTER=<ID=LowQual,Description=\"\">\n\
##FILTER=<ID=q10,Description=\"\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t2\trs1\tA\tC\t99\tLowQual;q10\tAC=1;AN=2;DP=12\n"
            .to_owned()
    }

    #[test]
    fn renames_contig_headers_and_records() {
        let mut map = HashMap::new();
        map.insert("1".to_owned(), "chr1".to_owned());
        let mut text = sample_vcf();
        transform_vcf_text(&mut text, Some(&map), &RemoveSet::default(), false);
        assert!(text.contains("##contig=<ID=chr1,length=10>"));
        assert!(text.contains("chr1\t2\trs1\tA\tC\t99\tLowQual;q10\tAC=1;AN=2;DP=12"));
    }

    #[test]
    fn removes_id_and_qual_columns() {
        let remove = RemoveSet {
            id: true,
            qual: true,
            ..Default::default()
        };
        let mut text = sample_vcf();
        transform_vcf_text(&mut text, None, &remove, false);
        assert!(text.contains("1\t2\t.\tA\tC\t.\tLowQual;q10\tAC=1;AN=2;DP=12"));
    }

    #[test]
    fn removes_specific_info_tags_and_their_headers() {
        let mut remove = RemoveSet::default();
        remove.info_ids.insert("AC".into());
        remove.info_ids.insert("DP".into());
        let mut text = sample_vcf();
        transform_vcf_text(&mut text, None, &remove, false);
        assert!(!text.contains("##INFO=<ID=AC,"));
        assert!(!text.contains("##INFO=<ID=DP,"));
        assert!(text.contains("##INFO=<ID=AN,"));
        assert!(text.contains("AN=2\n"));
        assert!(!text.contains("AC=1"));
        assert!(!text.contains("DP=12"));
    }

    #[test]
    fn removes_specific_filter_and_substitutes_pass() {
        let mut remove = RemoveSet::default();
        remove.filter_ids.insert("LowQual".into());
        remove.filter_ids.insert("q10".into());
        let mut text = sample_vcf();
        transform_vcf_text(&mut text, None, &remove, false);
        assert!(!text.contains("##FILTER=<ID=LowQual,"));
        assert!(!text.contains("##FILTER=<ID=q10,"));
        assert!(text.contains("PASS\tAC=1;AN=2;DP=12"));
    }

    #[test]
    fn removes_all_info() {
        let remove = RemoveSet {
            all_info: true,
            ..Default::default()
        };
        let mut text = sample_vcf();
        transform_vcf_text(&mut text, None, &remove, false);
        assert!(!text.contains("##INFO="));
        assert!(text.contains("LowQual;q10\t.\n"));
    }
}
