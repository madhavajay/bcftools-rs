//! `bcftools +fill-from-fasta` (upstream `bcftools/plugins/fill-from-fasta.c`).
//!
//! Fills the REF allele or an INFO tag from a FASTA reference. Plugin
//! `-i`/`-e` filters route through the shared text filter engine and leave
//! non-matching records unchanged.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::vcf_compat::normalize_vcf_text;

/// What `-c` selects.
enum Column {
    Ref,
    Info(String),
}

#[derive(Clone, Copy)]
pub enum FilterMode {
    Include,
    Exclude,
}

#[derive(Clone, Copy)]
pub struct FilterSpec<'a> {
    pub mode: FilterMode,
    pub expr: &'a str,
}

/// Parses a FASTA into name→sequence (name = first token after `>`).
fn parse_fasta(text: &str) -> HashMap<String, Vec<u8>> {
    let mut map = HashMap::new();
    let mut cur: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = cur.take() {
                map.insert(name, std::mem::take(&mut seq));
            }
            let name = rest.split_whitespace().next().unwrap_or("").to_string();
            cur = Some(name);
        } else {
            seq.extend(line.trim_end().bytes());
        }
    }
    if let Some(name) = cur.take() {
        map.insert(name, seq);
    }
    map
}

/// Reads inputs and returns the fill-from-fasta VCF text.
pub fn run(
    input: &Path,
    fasta: &Path,
    column: &str,
    header_lines: Option<&Path>,
    replace_non_acgtn: bool,
    filter: Option<FilterSpec<'_>>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let fa_text = fs::read_to_string(fasta)?;
    let hdr_lines = match header_lines {
        Some(p) => Some(fs::read_to_string(p)?),
        None => None,
    };
    compute(
        &text,
        &fa_text,
        column,
        hdr_lines.as_deref(),
        replace_non_acgtn,
        filter,
    )
    .map_err(io::Error::other)
}

fn compute(
    text: &str,
    fa_text: &str,
    column: &str,
    header_lines: Option<&str>,
    replace_non_acgtn: bool,
    filter: Option<FilterSpec<'_>>,
) -> Result<String, String> {
    let col = if column.eq_ignore_ascii_case("REF") {
        Column::Ref
    } else {
        let c = column.strip_prefix("INFO/").unwrap_or(column);
        Column::Info(c.to_string())
    };
    let fasta = parse_fasta(fa_text);

    let mut out = String::new();
    for line in text.lines() {
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                if let Some(h) = header_lines {
                    for hl in h.lines() {
                        if !hl.trim().is_empty() {
                            out.push_str(hl);
                            out.push('\n');
                        }
                    }
                }
                out.push_str(line);
                out.push('\n');
                continue;
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let mut f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if f.len() < 8 {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if let Some(filter) = filter
            && !record_should_annotate(&f, filter)?
        {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let pos: usize = f[1].parse().map_err(|_| format!("bad POS: {}", f[1]))?;
        let ref_len = f[3].len();
        let seq = fasta
            .get(&f[0])
            .ok_or_else(|| format!("faidx_fetch_seq failed at {}:{}", f[0], pos))?;
        let beg = pos - 1;
        let end = (beg + ref_len).min(seq.len());
        let mut fa: Vec<u8> = seq[beg..end].to_vec();
        for b in fa.iter_mut() {
            if *b > 96 {
                *b -= 32;
            }
            if replace_non_acgtn
                && *b != b'A'
                && *b != b'C'
                && *b != b'G'
                && *b != b'T'
                && *b != b'N'
            {
                *b = b'N';
            }
        }
        let fa = String::from_utf8_lossy(&fa).into_owned();

        match &col {
            Column::Ref => {
                f[3] = fa;
            }
            Column::Info(tag) => {
                // bcf_update_info_string: set/replace the tag.
                let info = std::mem::take(&mut f[7]);
                let entry = format!("{tag}={fa}");
                f[7] = if info == "." || info.is_empty() {
                    entry
                } else {
                    let mut kept: Vec<&str> = info
                        .split(';')
                        .filter(|kv| kv.split('=').next() != Some(tag.as_str()) && *kv != tag)
                        .collect();
                    kept.push(&entry);
                    kept.join(";")
                };
            }
        }
        out.push_str(&f.join("\t"));
        out.push('\n');
    }
    Ok(out)
}

fn record_should_annotate(fields: &[String], filter: FilterSpec<'_>) -> Result<bool, String> {
    let matched =
        bcffilter::eval_expression_with(filter.expr, &EvalContext::new(), |name, sample_index| {
            if sample_index.is_some() {
                return None;
            }
            fill_from_fasta_lookup(name, fields)
        })
        .map_err(|e| e.to_string())?
        .truthy();
    Ok(match filter.mode {
        FilterMode::Include => matched,
        FilterMode::Exclude => !matched,
    })
}

fn fill_from_fasta_lookup(token: &str, fields: &[String]) -> Option<FilterValue> {
    if token.eq_ignore_ascii_case("TYPE") {
        return Some(variant_type_value(fields));
    }
    crate::commands::filter::record_lookup(token, fields)
}

fn variant_type_value(fields: &[String]) -> FilterValue {
    let reference = fields.get(3).map(String::as_str).unwrap_or("");
    let values: Vec<FilterValue> = fields
        .get(4)
        .map(String::as_str)
        .unwrap_or(".")
        .split(',')
        .filter(|alt| *alt != ".")
        .map(|alt| {
            let kind = if reference.len() == alt.len() {
                if reference.len() == 1 { "snp" } else { "mnp" }
            } else {
                "indel"
            };
            FilterValue::String(kind.to_owned())
        })
        .collect();
    match values.as_slice() {
        [] => FilterValue::String("ref".to_owned()),
        [single] => single.clone(),
        _ => FilterValue::List(values),
    }
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
        ".bcftools-rs-fill-from-fasta-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fasta_parse_and_name() {
        let m = parse_fasta(">2 2:1-9\nACGT\nacgt\n>3 x\nTTTT\n");
        assert_eq!(m.get("2").unwrap(), b"ACGTacgt");
        assert_eq!(m.get("3").unwrap(), b"TTTT");
    }

    #[test]
    fn filter_controls_annotation_not_record_output() {
        let vcf = "\
##fileformat=VCFv4.3
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1\t.\tA\tC\t.\tPASS\t.
chr1\t2\t.\tA\tAT\t.\tPASS\t.
";
        let fasta = ">chr1\nACGT\n";
        let header = "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">\n";

        let include = compute(
            vcf,
            fasta,
            "AA",
            Some(header),
            false,
            Some(FilterSpec {
                mode: FilterMode::Include,
                expr: "TYPE=\"snp\"",
            }),
        )
        .unwrap();
        assert!(include.contains("chr1\t1\t.\tA\tC\t.\tPASS\tAA=A\n"));
        assert!(include.contains("chr1\t2\t.\tA\tAT\t.\tPASS\t.\n"));

        let exclude = compute(
            vcf,
            fasta,
            "AA",
            Some(header),
            false,
            Some(FilterSpec {
                mode: FilterMode::Exclude,
                expr: "TYPE=\"snp\"",
            }),
        )
        .unwrap();
        assert!(exclude.contains("chr1\t1\t.\tA\tC\t.\tPASS\t.\n"));
        assert!(exclude.contains("chr1\t2\t.\tA\tAT\t.\tPASS\tAA=C\n"));
    }
}
