//! `bcftools +split-vep` (upstream `bcftools/plugins/split-vep.c`).
//!
//! First slice: the `-c FIELD[,FIELD...] -s TR:CSQ[:PRN]` path that
//! annotates `INFO/<FIELD>` from a VEP/BCSQ `CSQ` (or `-a TAG`) string,
//! validated by piping through `bcftools query`. Implements the upstream
//! severity scale (`default_severity()`), `csq_to_severity` (lowercase,
//! `&`-split, severity = the first scale line whose token is a substring
//! of the term, unknown → `nscale+1`), transcript selection
//! (`all` / `worst` / `primary`=`CANONICAL=YES` / `pick`=`PICK=1` /
//! `mane`=`MANE_SELECT!=""` / `FIELD<OP>VALUE`), the `CSQ` severity
//! threshold (`+`/`-`/`=`), and `PRN` (`all`/`worst`). The CSQ field
//! list comes from the tag's `##INFO=<...,Description="...Format: a|b">`
//! header.
//!
//! Deferred (tracked in TODO.md): the `-f` format-string output
//! (needs the `convert` engine), `-d`/`--duplicate`, per-sample
//! `[%SAMPLE]` blocks, `-t`/`-T` regions, `-g`/`--gene-list`,
//! `-x`/`--drop-sites` interactions beyond the default, `-S` custom
//! severity file, and `--columns-types`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

/// Upstream `default_severity()` — consequence substrings in ascending
/// severity, one line per severity level (tokens on a line share it).
const DEFAULT_SEVERITY: &[&[&str]] = &[
    &["intergenic"],
    &["feature_truncation", "feature_elongation"],
    &["regulatory"],
    &["tf_binding_site", "tfbs"],
    &["downstream", "upstream"],
    &["non_coding_transcript", "non_coding"],
    &["intron", "nmd_transcript"],
    &["non_coding_transcript_exon"],
    &["5_prime_utr", "3_prime_utr"],
    &["coding_sequence", "mature_mirna"],
    &["stop_retained", "start_retained", "synonymous"],
    &["incomplete_terminal_codon"],
    &["splice_region"],
    &["missense", "inframe", "protein_altering"],
    &["transcript_amplification"],
    &["exon_loss"],
    &["disruptive"],
    &["start_lost", "stop_lost", "stop_gained", "frameshift"],
    &["splice_acceptor", "splice_donor"],
    &["transcript_ablation"],
];

/// Severity scale: ordered scale tokens + a consequence→severity cache
/// that mirrors upstream's dynamic extension for unknown terms.
struct Severity {
    /// (token, severity) in scale order; tokens are lowercased.
    scale: Vec<(String, i32)>,
    nscale: i32,
    cache: std::collections::HashMap<String, i32>,
}

impl Severity {
    fn default_scale() -> Self {
        let mut scale = Vec::new();
        let mut cache = std::collections::HashMap::new();
        for (sev, line) in DEFAULT_SEVERITY.iter().enumerate() {
            for &tok in *line {
                scale.push((tok.to_owned(), sev as i32));
                cache.entry(tok.to_owned()).or_insert(sev as i32);
            }
        }
        let nscale = scale.len() as i32;
        Severity {
            scale,
            nscale,
            cache,
        }
    }

    /// Severity of a single (already lowercased) consequence term,
    /// mirroring upstream `csq_to_severity`'s per-term lookup +
    /// dynamic-extension behavior.
    fn term_severity(&mut self, term: &str) -> i32 {
        if let Some(&s) = self.cache.get(term) {
            return s;
        }
        let mut sev = self.nscale + 1;
        for (tok, s) in &self.scale {
            if term.contains(tok.as_str()) {
                sev = *s;
                break;
            }
        }
        self.nscale += 1;
        self.scale.push((term.to_owned(), sev));
        self.cache.insert(term.to_owned(), sev);
        sev
    }

    /// `csq_to_severity` with `exact_match = -1`: min & max severity
    /// over the `&`-split terms of a consequence string.
    fn min_max(&mut self, csq: &str) -> (i32, i32) {
        let lower = csq.to_ascii_lowercase();
        let (mut mn, mut mx) = (i32::MAX, -1);
        for term in lower.split('&') {
            let s = self.term_severity(term);
            mn = mn.min(s);
            mx = mx.max(s);
        }
        (mn, mx)
    }
}

#[derive(Clone)]
enum SelectTr {
    All,
    Worst,
    /// `FIELD <OP> VALUE` (e.g. `CANONICAL=YES`, `MANE_SELECT!=""`).
    Expr {
        field: String,
        ne: bool,
        value: String,
    },
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum PrnCsq {
    All,
    Worst,
}

struct Select {
    tr: SelectTr,
    /// `SELECT_CSQ_ANY` => no threshold.
    any: bool,
    min_severity: i32,
    max_severity: i32,
    prn: PrnCsq,
}

fn parse_select(spec: &str, sev: &Severity) -> Result<Select, String> {
    let cols: Vec<&str> = spec.split(':').collect();
    let sel_tr = cols
        .first()
        .filter(|s| !s.is_empty())
        .copied()
        .unwrap_or("all");
    let sel_csq = cols
        .get(1)
        .filter(|s| !s.is_empty())
        .copied()
        .unwrap_or("any");
    let prn_csq = cols
        .get(2)
        .filter(|s| !s.is_empty())
        .copied()
        .unwrap_or("all");

    let tr = match sel_tr.to_ascii_lowercase().as_str() {
        "all" => SelectTr::All,
        "worst" => SelectTr::Worst,
        "primary" => SelectTr::Expr {
            field: "CANONICAL".into(),
            ne: false,
            value: "YES".into(),
        },
        "pick" => SelectTr::Expr {
            field: "PICK".into(),
            ne: false,
            value: "1".into(),
        },
        "mane" => SelectTr::Expr {
            field: "MANE_SELECT".into(),
            ne: true,
            value: String::new(),
        },
        _ => parse_tr_expr(sel_tr)?,
    };

    let (any, mut min_severity, mut max_severity) = (sel_csq == "any", -1, -1);
    if !any {
        let mut s = sel_csq.to_string();
        let modifier = s.chars().last().filter(|c| *c == '+' || *c == '-');
        if modifier.is_some() {
            s.pop();
        }
        let severity = *sev
            .cache
            .get(s.to_ascii_lowercase().as_str())
            .ok_or_else(|| format!("the consequence \"{s}\" is not recognised"))?;
        match modifier {
            Some('+') => {
                min_severity = severity;
                max_severity = i32::MAX;
            }
            Some('-') => {
                min_severity = 0;
                max_severity = severity;
            }
            _ => {
                min_severity = severity;
                max_severity = severity;
            }
        }
    }
    let prn = match prn_csq.to_ascii_lowercase().as_str() {
        "all" => PrnCsq::All,
        "worst" => PrnCsq::Worst,
        _ => return Err(format!("could not parse \"{prn_csq}\" in -s \"{spec}\"")),
    };
    Ok(Select {
        tr,
        any,
        min_severity,
        max_severity,
        prn,
    })
}

/// Upstream `init_select_tr_expr`: `FIELD`, `=`/`==`/`!=`, `VALUE`
/// (quotes stripped). Only equality/inequality is exercised by the
/// `primary`/`pick`/`mane` aliases and the upstream `-s` fixtures.
fn parse_tr_expr(s: &str) -> Result<SelectTr, String> {
    let (field, ne, value) = if let Some((f, v)) = s.split_once("!=") {
        (f, true, v)
    } else if let Some((f, v)) = s.split_once("==") {
        (f, false, v)
    } else if let Some((f, v)) = s.split_once('=') {
        (f, false, v)
    } else {
        return Err(format!(
            "could not parse the -s transcript expression \"{s}\""
        ));
    };
    Ok(SelectTr::Expr {
        field: field.to_owned(),
        ne,
        value: value.trim_matches('"').to_owned(),
    })
}

/// `csq_rewrite_worst`: reduce an `&`-joined consequence to its single
/// most severe term (keeps the original order/case of that term).
fn rewrite_worst(csq: &str, sev: &mut Severity) -> String {
    let parts: Vec<&str> = csq.split('&').collect();
    if parts.len() <= 1 {
        return csq.to_owned();
    }
    let mut imax = 0;
    let mut smax = -1;
    for (i, p) in parts.iter().enumerate() {
        let s = sev.term_severity(&p.to_ascii_lowercase());
        if smax < s {
            smax = s;
            imax = i;
        }
    }
    parts[imax].to_owned()
}

pub struct Options<'a> {
    /// `-c` field list (names as they appear in the CSQ `Format:`).
    pub columns: &'a str,
    /// `-s TR:CSQ:PRN` (default `all:any`).
    pub select: &'a str,
    /// `-a`/`--annotation` INFO tag (default `CSQ`).
    pub annotation: &'a str,
    /// `-f`/`--format` query-style format string. When set, the plugin
    /// emits rendered text (records with no severity-passing transcript
    /// are dropped) instead of an annotated VCF, mirroring upstream's
    /// `convert_init` text path.
    pub format: Option<&'a str>,
}

/// Reads the input VCF/BCF and returns either the `-c`-annotated VCF
/// text or, when `-f` is given, the rendered text output (upstream's
/// `convert_init` path: non-matching sites are dropped, and the format
/// string is rendered by our own `bcftools query` engine over the
/// transiently annotated VCF).
pub fn run(input: &Path, opts: Options<'_>) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let tag = opts.annotation;
    let fields = parse_format_tokens(&text, tag).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("Expected INFO/{tag} with a \"Format: a|b|c\" description"),
        )
    })?;
    let field_idx = |name: &str| fields.iter().position(|f| f == name);
    let csq_idx = field_idx("Consequence").ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "the \"Consequence\" subfield is required",
        )
    })?;

    let mut sev = Severity::default_scale();
    let sel = parse_select(opts.select, &sev)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("split-vep: {e}")))?;

    // Fields to annotate as transient INFO: the explicit `-c` columns,
    // plus any CSQ subfield referenced by a `%TOKEN` in the `-f` format
    // string (upstream parses the same set out of the format string).
    let mut names: Vec<String> = Vec::new();
    for c in opts.columns.split(',').filter(|c| !c.is_empty()) {
        if field_idx(c).is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("the field \"{c}\" is not present in INFO/{tag}"),
            ));
        }
        if !names.iter().any(|n| n == c) {
            names.push(c.to_owned());
        }
    }
    if let Some(fmt) = opts.format {
        for tok in format_field_tokens(fmt) {
            if fields.iter().any(|f| f == &tok) && !names.iter().any(|n| n == &tok) {
                names.push(tok);
            }
        }
    }
    let annots: Vec<(String, usize)> = names
        .iter()
        .map(|n| (n.clone(), field_idx(n).unwrap()))
        .collect();

    // Drop non-matching sites only in `-f` (text) mode; the `-c` VCF
    // path keeps every record (annotation just absent when unmatched).
    let drop_unmatched = opts.format.is_some();

    let mut out = String::with_capacity(text.len());
    let mut hdr_done = false;
    for line in text.lines() {
        if line.starts_with('#') {
            if line.starts_with("#CHROM") && !hdr_done {
                for (name, _) in &annots {
                    out.push_str(&format!(
                        "##INFO=<ID={name},Number=.,Type=String,Description=\"The {name} field from INFO/{tag}\">\n"
                    ));
                }
                hdr_done = true;
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let (rendered, passed, all_missing) =
            process_record(line, tag, &fields, csq_idx, &sel, &annots, &mut sev);
        // Upstream `-f` default is `--drop-sites`: skip a record that has
        // no severity-passing transcript, or (when CSQ subfields are
        // requested) whose every requested annotation is missing.
        if drop_unmatched && (!passed || (!annots.is_empty() && all_missing)) {
            continue;
        }
        out.push_str(&rendered);
        out.push('\n');
    }

    match opts.format {
        None => Ok(out),
        Some(fmt) => {
            let mut buf: Vec<u8> = Vec::with_capacity(out.len());
            crate::commands::query::query_format_text(
                &out,
                fmt,
                &crate::commands::query::QueryFormatOptions::plain(),
                &mut buf,
            )?;
            String::from_utf8(buf)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))
        }
    }
}

/// Identifiers referenced by `%NAME` / `%INFO/NAME` tokens in a query
/// format string (used to decide which CSQ subfields to annotate).
fn format_field_tokens(fmt: &str) -> Vec<String> {
    let bytes = fmt.as_bytes();
    let mut toks = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%' {
            let mut j = i + 1;
            while j < bytes.len()
                && (bytes[j].is_ascii_alphanumeric() || bytes[j] == b'_' || bytes[j] == b'/')
            {
                j += 1;
            }
            if j > i + 1 {
                let name = &fmt[i + 1..j];
                // `%INFO/X` is an explicit real-INFO reference (resolved
                // by the query engine), never a CSQ-subfield request, so
                // it must not trigger transient CSQ annotation.
                if !name.starts_with("INFO/") {
                    toks.push(name.to_owned());
                }
            }
            i = j;
        } else {
            i += 1;
        }
    }
    toks
}

#[allow(clippy::too_many_arguments)]
fn process_record(
    line: &str,
    tag: &str,
    fields: &[String],
    csq_idx: usize,
    sel: &Select,
    annots: &[(String, usize)],
    sev: &mut Severity,
) -> (String, bool, bool) {
    let mut f: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if f.len() < 8 {
        return (line.to_owned(), false, true);
    }
    let csq_str = f[7]
        .split(';')
        .find_map(|kv| kv.strip_prefix(&format!("{tag}=")));
    let Some(csq_str) = csq_str else {
        return (line.to_owned(), false, true); // no CSQ: emit unchanged
    };

    // transcripts -> field vectors
    let transcripts: Vec<Vec<&str>> = csq_str.split(',').map(|t| t.split('|').collect()).collect();

    // Select transcripts.
    let selected: Vec<usize> = match &sel.tr {
        SelectTr::All => (0..transcripts.len()).collect(),
        SelectTr::Worst => {
            let mut imax = 0;
            let mut smax = -1;
            for (i, tr) in transcripts.iter().enumerate() {
                let csq = tr.get(csq_idx).copied().unwrap_or("");
                let (_, mx) = sev.min_max(csq);
                if smax < mx {
                    smax = mx;
                    imax = i;
                }
            }
            if transcripts.is_empty() {
                vec![]
            } else {
                vec![imax]
            }
        }
        SelectTr::Expr { field, ne, value } => {
            let fi = fields.iter().position(|x| x == field);
            (0..transcripts.len())
                .filter(|&i| {
                    let v = fi
                        .and_then(|fi| transcripts[i].get(fi).copied())
                        .unwrap_or("");
                    if *ne { v != value } else { v == value }
                })
                .collect()
        }
    };

    // Collect per-annot values across selected transcripts that pass
    // the CSQ severity threshold.
    let mut acc: Vec<Vec<String>> = vec![Vec::new(); annots.len()];
    let mut severity_pass = false;
    // Upstream `all_missing`: every requested annotation empty across all
    // severity-passing transcripts (record dropped in the `-f` path).
    let mut all_missing = true;
    for &ti in &selected {
        let tr = &transcripts[ti];
        let csq = tr.get(csq_idx).copied().unwrap_or("");
        if !csq_severity_pass(csq, sel, sev) {
            continue;
        }
        severity_pass = true;
        for (ai, (_, idx)) in annots.iter().enumerate() {
            let raw = tr.get(*idx).copied().unwrap_or("");
            if !raw.is_empty() {
                all_missing = false;
            }
            let val = if *idx == csq_idx && sel.prn == PrnCsq::Worst {
                rewrite_worst(raw, sev)
            } else {
                raw.to_owned()
            };
            acc[ai].push(if val.is_empty() { ".".to_owned() } else { val });
        }
    }

    if severity_pass {
        for (ai, (name, _)) in annots.iter().enumerate() {
            let joined = acc[ai].join(",");
            if !joined.is_empty() {
                set_info(&mut f, name, &joined);
            }
        }
    }
    (f.join("\t"), severity_pass, all_missing)
}

fn csq_severity_pass(csq: &str, sel: &Select, sev: &mut Severity) -> bool {
    if sel.any {
        return true;
    }
    let (mn, mx) = sev.min_max(csq);
    if mx < sel.min_severity {
        return false;
    }
    if mn > sel.max_severity {
        return false;
    }
    true
}

/// `bcf_update_info_string`: set/replace `INFO/<key>`.
fn set_info(f: &mut [String], key: &str, value: &str) {
    let entry = format!("{key}={value}");
    let info = &f[7];
    f[7] = if info == "." || info.is_empty() {
        entry
    } else {
        let mut kept: Vec<&str> = info
            .split(';')
            .filter(|kv| kv.split('=').next() != Some(key) && *kv != key)
            .collect();
        let e = entry.clone();
        kept.push(&e);
        kept.join(";")
    };
}

/// Extracts the `Format: a|b|c` token list from
/// `##INFO=<ID=tag,...,Description="...Format: ...">` (same rule as
/// `+vcf2table`).
fn parse_format_tokens(text: &str, id: &str) -> Option<Vec<String>> {
    let needle = format!("##INFO=<ID={id},");
    let line = text.lines().find(|l| l.starts_with(&needle))?;
    let dstart = line.find("Description=\"")? + "Description=\"".len();
    let drest = &line[dstart..];
    let dend = drest.find('"').unwrap_or(drest.len());
    let desc = &drest[..dend];
    let fstart = desc.find("Format: ")? + "Format: ".len();
    Some(desc[fstart..].split('|').map(|s| s.to_owned()).collect())
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
        MultiGzDecoder::new(file).read_to_string(&mut text)?;
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
        ".bcftools-rs-split-vep-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn severity_scale_basics() {
        let mut s = Severity::default_scale();
        // missense line is rank 13. Upstream `csq_to_severity` breaks at
        // the first scale-order token that is a substring, so
        // `non_coding_transcript_exon_variant` resolves via the earlier
        // `non_coding_transcript` token (rank 5), not `..._exon` (rank 7).
        assert_eq!(s.term_severity("missense_variant"), 13);
        assert_eq!(s.term_severity("non_coding_transcript_exon_variant"), 5);
        assert_eq!(s.term_severity("intergenic_variant"), 0);
        // `&`-join: min/max across terms.
        assert_eq!(
            s.min_max("missense_variant&splice_region_variant"),
            (12, 13)
        );
    }

    #[test]
    fn select_parsing() {
        let sev = Severity::default_scale();
        let s = parse_select("worst:missense+", &sev).unwrap();
        assert!(matches!(s.tr, SelectTr::Worst));
        assert_eq!(s.min_severity, 13);
        assert_eq!(s.max_severity, i32::MAX);
        assert_eq!(s.prn, PrnCsq::All);
        let s2 = parse_select("primary:missense+:worst", &sev).unwrap();
        assert!(matches!(s2.tr, SelectTr::Expr { .. }));
        assert_eq!(s2.prn, PrnCsq::Worst);
        let s3 = parse_select("worst", &sev).unwrap();
        assert!(s3.any);
    }

    #[test]
    fn rewrite_worst_keeps_most_severe_term() {
        let mut s = Severity::default_scale();
        assert_eq!(
            rewrite_worst("splice_region_variant&missense_variant", &mut s),
            "missense_variant"
        );
        assert_eq!(rewrite_worst("intron_variant", &mut s), "intron_variant");
    }
}
