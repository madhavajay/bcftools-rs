//! `bcftools +fill-tags` (upstream `bcftools/plugins/fill-tags.c`).
//!
//! Implemented: the genotype-derived INFO tags
//! `F_MISSING`/`NS`/`AN`/`AF`/`MAF`/`AC`/`AC_Het`/`AC_Hom`/`AC_Hemi`/
//! `HWE`/`ExcHet`, the `-t LIST` selection, the `all` / default
//! (no-`-t`) set, `-S`/`--samples-file` population grouping (per-group
//! `<TAG>_<pop>` plus the always-present global tag), and
//! `-d`/`--drop-missing`. Counting mirrors upstream `process_fmt`'s
//! per-sample hom/het/hemi/half-missing classification; `HWE`/`ExcHet`
//! port `calc_hwe` (Wigginton 2005). Tags are written in the upstream
//! fixed order and floats use C `%g`/6 over the f32-stored value.
//!
//! `FORMAT/VAF`+`VAF1` are computed from `FORMAT/AD` (upstream
//! `process_vaf_vaf1`), independent of `GT`. The `TAG:Num=EXPR`
//! function engine is partially ported:
//! `[int|float](sum|smpl_sum(INFO/X|FMT/X))` and
//! `[int](F_PASS|N_PASS(EXPR))` (per-pop, EXPR evaluated per sample via
//! the shared filter engine), plus the `Added by +fill-tags
//! expression …` header.
//!
//! Deferred (tracked in TODO.md): `END`/`TYPE`, and the remaining
//! function engine (`ssum`/`fisher`/`binom`/`phred`, arithmetic,
//! `[*:i]` subscripts).

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tag {
    FMissing,
    Ns,
    An,
    Af,
    Maf,
    Ac,
    AcHet,
    AcHom,
    AcHemi,
    Hwe,
    ExcHet,
    Vaf,
    Vaf1,
}

/// Fixed upstream `process_fmt` write order: the `F_MISSING` func is
/// emitted first, then the SET_ block, then `HWE`/`ExcHet`.
const WRITE_ORDER: &[Tag] = &[
    Tag::FMissing,
    Tag::Ns,
    Tag::An,
    Tag::Af,
    Tag::Maf,
    Tag::Ac,
    Tag::AcHet,
    Tag::AcHom,
    Tag::AcHemi,
    Tag::Hwe,
    Tag::ExcHet,
];

/// The `all` / default (no `-t`) tag set in this slice. Upstream `all`
/// is `~(SET_END|SET_TYPE)` plus the `F_MISSING` func; `VAF`/`VAF1` are
/// FORMAT tags that only emit with `FORMAT/AD` (deferred — a no-op for
/// the GT-only fixtures this covers).
const ALL_TAGS: &[Tag] = &[
    Tag::FMissing,
    Tag::Ns,
    Tag::An,
    Tag::Af,
    Tag::Maf,
    Tag::Ac,
    Tag::AcHet,
    Tag::AcHom,
    Tag::AcHemi,
    Tag::Hwe,
    Tag::ExcHet,
];

fn parse_tag(name: &str) -> Result<Tag, String> {
    let n = name
        .strip_prefix("INFO/")
        .or_else(|| name.strip_prefix("FORMAT/"))
        .unwrap_or(name);
    match n.to_ascii_uppercase().as_str() {
        "NS" => Ok(Tag::Ns),
        "AN" => Ok(Tag::An),
        "AF" => Ok(Tag::Af),
        "MAF" => Ok(Tag::Maf),
        "AC" => Ok(Tag::Ac),
        "AC_HET" => Ok(Tag::AcHet),
        "AC_HOM" => Ok(Tag::AcHom),
        "AC_HEMI" => Ok(Tag::AcHemi),
        "HWE" => Ok(Tag::Hwe),
        "EXCHET" => Ok(Tag::ExcHet),
        "F_MISSING" => Ok(Tag::FMissing),
        "VAF" => Ok(Tag::Vaf),
        "VAF1" => Ok(Tag::Vaf1),
        _ => Err(format!(
            "the tag \"{name}\" is not yet ported (this slice supports \
             AN,AC,AC_Hom,AC_Het,AC_Hemi,AF,MAF,NS,HWE,ExcHet,F_MISSING,all)"
        )),
    }
}

impl Tag {
    /// Base ID and the `hdr_append` `##INFO` template (with the `{S}`
    /// suffix and `{IN}` (` in <pop>`) placeholders) used when the tag's
    /// ID is not already declared in the header.
    fn header(self) -> (&'static str, &'static str) {
        match self {
            Tag::An => (
                "AN",
                "##INFO=<ID=AN{S},Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes{IN}\">",
            ),
            Tag::Ac => (
                "AC",
                "##INFO=<ID=AC{S},Number=A,Type=Integer,Description=\"Allele count in genotypes{IN}\">",
            ),
            Tag::Ns => (
                "NS",
                "##INFO=<ID=NS{S},Number=1,Type=Integer,Description=\"Number of samples with data{IN}\">",
            ),
            Tag::AcHom => (
                "AC_Hom",
                "##INFO=<ID=AC_Hom{S},Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes{IN}\">",
            ),
            Tag::AcHet => (
                "AC_Het",
                "##INFO=<ID=AC_Het{S},Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes{IN}\">",
            ),
            Tag::AcHemi => (
                "AC_Hemi",
                "##INFO=<ID=AC_Hemi{S},Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes{IN}\">",
            ),
            Tag::Af => (
                "AF",
                "##INFO=<ID=AF{S},Number=A,Type=Float,Description=\"Allele frequency{IN}\">",
            ),
            Tag::Maf => (
                "MAF",
                "##INFO=<ID=MAF{S},Number=1,Type=Float,Description=\"Frequency of the second most common allele{IN}\">",
            ),
            Tag::FMissing => (
                "F_MISSING",
                "##INFO=<ID=F_MISSING{S},Number=1,Type=Float,Description=\"Added by +fill-tags expression F_MISSING:1=F_MISSING\">",
            ),
            Tag::Hwe => (
                "HWE",
                "##INFO=<ID=HWE{S},Number=A,Type=Float,Description=\"HWE test{IN} (PMID:15789306); 1=good, 0=bad\">",
            ),
            Tag::ExcHet => (
                "ExcHet",
                "##INFO=<ID=ExcHet{S},Number=A,Type=Float,Description=\"Test excess heterozygosity{IN}; 1=good, 0=bad\">",
            ),
            // FORMAT tags — handled by the dedicated VAF step, never via
            // the INFO `HDR_ORDER`/`WRITE_ORDER` paths.
            Tag::Vaf => (
                "VAF",
                "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"The fraction of reads with alternate allele (nALT/nSumAll)\">",
            ),
            Tag::Vaf1 => (
                "VAF1",
                "##FORMAT=<ID=VAF1,Number=1,Type=Float,Description=\"The fraction of reads with alternate alleles (nSumALT/nSumAll)\">",
            ),
        }
    }
    fn base_id(self) -> &'static str {
        self.header().0
    }
}

/// `hdr_append` order: the `F_MISSING` func header is appended first,
/// then upstream lines 590-598, then `HWE`, then `ExcHet`.
const HDR_ORDER: &[Tag] = &[
    Tag::FMissing,
    Tag::An,
    Tag::Ac,
    Tag::Ns,
    Tag::AcHom,
    Tag::AcHet,
    Tag::AcHemi,
    Tag::Af,
    Tag::Maf,
    Tag::Hwe,
    Tag::ExcHet,
];

#[derive(Default, Clone)]
struct Counts {
    nhet: i64,
    nhom: i64,
    nhemi: i64,
    nac: i64,
}

struct Pop {
    /// `""` for the global pop; otherwise the population name.
    name: String,
    /// `""` (global) or `_<name>`.
    suffix: String,
}

#[derive(Clone, Copy, PartialEq)]
enum FuncKind {
    /// `sum(X)` — scalar total over (pop) samples and all values.
    Sum,
    /// `smpl_sum(X)` — per-sample total (a FORMAT result).
    SmplSum,
    /// `F_PASS(EXPR)` — fraction of (pop) samples where EXPR is true.
    FPass,
    /// `N_PASS(EXPR)` — number of (pop) samples where EXPR is true.
    NPass,
}

/// A `[INFO/|FORMAT/]DST[:CNT]=[int(|float(]<KIND>(<OPERAND>))` tag,
/// the subset of upstream's `ftf_filter_expr` engine this slice ports.
struct Func {
    dst: String,
    /// `true` → INFO result, `false` → FORMAT result.
    info: bool,
    /// `Some(n)` → `Number=n` (`:CNT`), `None` → `Number=.`.
    cnt: Option<i64>,
    is_int: bool,
    kind: FuncKind,
    /// Operand: `true` = `INFO/<tag>`, `false` = `FORMAT/<tag>`
    /// (`Sum`/`SmplSum` only).
    operand_info: bool,
    operand_tag: String,
    /// The per-sample filter expression (`FPass`/`NPass` only).
    pass_expr: String,
    /// The original `-t` token (for the `Added by +fill-tags
    /// expression …` header description).
    raw: String,
}

/// Split the `-t` list on commas at parenthesis depth 0 so that
/// function arguments (e.g. `binom(FMT/AD[:0],FMT/AD[:1])`) are kept
/// intact.
fn split_tags(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut depth = 0i32;
    for c in s.chars() {
        match c {
            '(' => {
                depth += 1;
                cur.push(c);
            }
            ')' => {
                depth -= 1;
                cur.push(c);
            }
            ',' if depth == 0 => out.push(std::mem::take(&mut cur)),
            _ => cur.push(c),
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

fn parse_func(token: &str) -> Result<Func, String> {
    let (lhs, expr) = token
        .split_once('=')
        .ok_or_else(|| format!("could not parse \"{token}\""))?;
    let (info, name) = if let Some(r) = lhs
        .strip_prefix("INFO/")
        .or_else(|| lhs.strip_prefix("info/"))
    {
        (true, r)
    } else if let Some(r) = lhs
        .strip_prefix("FORMAT/")
        .or_else(|| lhs.strip_prefix("format/"))
        .or_else(|| lhs.strip_prefix("FMT/"))
        .or_else(|| lhs.strip_prefix("fmt/"))
    {
        (false, r)
    } else {
        (true, lhs)
    };
    let (dst, cnt) = match name.split_once(':') {
        Some((d, c)) => (
            d.to_owned(),
            Some(
                c.parse::<i64>()
                    .map_err(|_| format!("could not parse \"{token}\""))?,
            ),
        ),
        None => (name.to_owned(), None),
    };

    // Optional int()/integer()/float() wrapper.
    let mut e = expr.trim();
    let mut is_int = false;
    for (pfx, isi) in [("int(", true), ("integer(", true), ("float(", false)] {
        if e.len() > pfx.len() && e[..pfx.len()].eq_ignore_ascii_case(pfx) && e.ends_with(')') {
            e = e[pfx.len()..e.len() - 1].trim();
            is_int = isi;
            break;
        }
    }

    // `<KIND>(<INNER>)`.
    if !e.ends_with(')') {
        return Err(format!("could not parse the expression \"{e}\""));
    }
    let (kind, inner) = if let Some(r) = e.strip_prefix("sum(") {
        (FuncKind::Sum, &r[..r.len() - 1])
    } else if let Some(r) = e.strip_prefix("smpl_sum(") {
        (FuncKind::SmplSum, &r[..r.len() - 1])
    } else if let Some(r) = e.strip_prefix("F_PASS(") {
        (FuncKind::FPass, &r[..r.len() - 1])
    } else if let Some(r) = e.strip_prefix("N_PASS(") {
        (FuncKind::NPass, &r[..r.len() - 1])
    } else {
        return Err(format!(
            "the expression \"{e}\" is not yet ported (this slice supports \
             int/float(sum|smpl_sum(INFO/TAG|FMT/TAG)) and F_PASS/N_PASS(EXPR))"
        ));
    };
    let inner = inner.trim();

    if matches!(kind, FuncKind::FPass | FuncKind::NPass) {
        return Ok(Func {
            dst,
            info,
            cnt,
            is_int,
            kind,
            operand_info: false,
            operand_tag: String::new(),
            pass_expr: inner.to_owned(),
            raw: token.to_owned(),
        });
    }

    let (operand_info, operand_tag) = if let Some(r) = inner
        .strip_prefix("INFO/")
        .or_else(|| inner.strip_prefix("info/"))
    {
        (true, r)
    } else if let Some(r) = inner
        .strip_prefix("FORMAT/")
        .or_else(|| inner.strip_prefix("format/"))
        .or_else(|| inner.strip_prefix("FMT/"))
        .or_else(|| inner.strip_prefix("fmt/"))
    {
        (false, r)
    } else {
        return Err(format!(
            "the operand \"{inner}\" must be INFO/TAG or FMT/TAG in this slice"
        ));
    };
    Ok(Func {
        dst,
        info,
        cnt,
        is_int,
        kind,
        operand_info,
        operand_tag: operand_tag.to_owned(),
        pass_expr: String::new(),
        raw: token.to_owned(),
    })
}

/// Evaluate a per-sample filter expression for sample `si` (the
/// `F_PASS`/`N_PASS` predicate), mirroring the `+setGT -t q` wiring:
/// a single-sample [`EvalContext`] over that sample's FORMAT values
/// with record-level fallback.
fn sample_passes(cols: &[&str], si: usize, expr: &str) -> bool {
    let fields: Vec<String> = cols.iter().map(|s| (*s).to_owned()).collect();
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let sample = &fields[9 + si];
    let values: Vec<&str> = sample.split(':').collect();
    let context = EvalContext::new().with_sample(
        format_keys
            .iter()
            .enumerate()
            .map(|(i, key)| {
                let raw = values.get(i).copied().unwrap_or(".");
                let value = if key.eq_ignore_ascii_case("GT") {
                    FilterValue::String(raw.to_owned())
                } else if raw == "." || raw.is_empty() {
                    FilterValue::Missing
                } else if let Ok(n) = raw.parse::<f64>() {
                    FilterValue::Number(n)
                } else if raw.contains(',') && raw.split(',').all(|v| v.parse::<f64>().is_ok()) {
                    FilterValue::List(
                        raw.split(',')
                            .map(|v| FilterValue::Number(v.parse().unwrap()))
                            .collect(),
                    )
                } else {
                    FilterValue::String(raw.to_owned())
                };
                ((*key).to_owned(), value)
            })
            .collect::<Vec<_>>(),
    );
    bcffilter::eval_expression_with(expr, &context, |name, sample_index| {
        if sample_index.is_some() {
            return None;
        }
        crate::commands::filter::record_lookup(name, &fields)
    })
    .map(|v| v.truthy())
    .unwrap_or(false)
}

/// Numeric values of an `INFO/<tag>` array (`.`/non-numeric skipped).
fn info_operand_values(cols: &[&str], tag: &str) -> Vec<f64> {
    let info = cols.get(7).copied().unwrap_or(".");
    if info == "." {
        return Vec::new();
    }
    let pfx = format!("{tag}=");
    for kv in info.split(';') {
        if let Some(v) = kv.strip_prefix(&pfx) {
            return v.split(',').filter_map(|x| x.parse::<f64>().ok()).collect();
        }
    }
    Vec::new()
}

/// Numeric values of `FORMAT/<func.operand_tag>` for sample `si`
/// (`.`/non-numeric skipped).
fn sample_operand_values(cols: &[&str], si: usize, func: &Func) -> Vec<f64> {
    let Some(fmt) = cols.get(8) else {
        return Vec::new();
    };
    let Some(idx) = fmt.split(':').position(|k| k == func.operand_tag) else {
        return Vec::new();
    };
    let Some(col) = cols.get(9 + si) else {
        return Vec::new();
    };
    col.split(':')
        .nth(idx)
        .unwrap_or(".")
        .split(',')
        .filter_map(|x| x.parse::<f64>().ok())
        .collect()
}

/// `sum(...)` for an INFO-result func: scalar over the pop's samples
/// (FMT operand) or the INFO array. `None` ⇒ no values ⇒ missing.
fn eval_sum(
    cols: &[&str],
    func: &Func,
    pop_idx: usize,
    sample_to_pops: &[Vec<usize>],
) -> Option<f64> {
    if func.kind != FuncKind::Sum {
        return None;
    }
    if func.operand_info {
        let v = info_operand_values(cols, &func.operand_tag);
        return if v.is_empty() {
            None
        } else {
            Some(v.iter().sum())
        };
    }
    let nsmpl = cols.len().saturating_sub(9);
    let mut acc = 0.0;
    let mut any = false;
    for si in 0..nsmpl {
        if !sample_to_pops
            .get(si)
            .is_some_and(|ps| ps.contains(&pop_idx))
        {
            continue;
        }
        for v in sample_operand_values(cols, si, func) {
            acc += v;
            any = true;
        }
    }
    if any { Some(acc) } else { None }
}

/// `int32_from_double`: bcftools truncates toward zero; a missing value
/// renders as `.`.
fn fmt_num(v: Option<f64>, is_int: bool) -> String {
    match v {
        None => ".".to_owned(),
        Some(x) if is_int => (x.trunc() as i64).to_string(),
        Some(x) => fmt_float(x),
    }
}

pub struct Options<'a> {
    /// `-t` comma-separated tag list; `"all"` (the default when `-t` is
    /// omitted) expands to [`ALL_TAGS`].
    pub tags: &'a str,
    /// `-S`/`--samples-file` path (`sample<ws>pop1,pop2` per line).
    pub samples_file: Option<&'a Path>,
    /// `-d`/`--drop-missing`: count half-missing `./1` genotypes via
    /// `nac` instead of as hemizygous.
    pub drop_missing: bool,
}

pub fn run(input: &Path, opts: Options<'_>) -> io::Result<String> {
    let mut want: Vec<Tag> = Vec::new();
    // Upstream `all` also enables FORMAT/VAF+VAF1: their ##FORMAT header
    // lines are declared even though values only appear with FORMAT/AD
    // (VAF computation itself deferred).
    let mut vaf_hdr = false;
    let mut funcs: Vec<Func> = Vec::new();
    for t in split_tags(opts.tags).iter().filter(|t| !t.is_empty()) {
        let t = t.as_str();
        if t.eq_ignore_ascii_case("all") {
            for &a in ALL_TAGS {
                if !want.contains(&a) {
                    want.push(a);
                }
            }
            vaf_hdr = true;
            continue;
        }
        if t.contains('=') {
            funcs.push(parse_func(t).map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidInput, format!("fill-tags: {e}"))
            })?);
            continue;
        }
        let tag = parse_tag(t)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("fill-tags: {e}")))?;
        if !want.contains(&tag) {
            want.push(tag);
        }
    }
    if want.is_empty() && funcs.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "fill-tags: empty tag list",
        ));
    }
    let want_vaf = want.contains(&Tag::Vaf);
    let want_vaf1 = want.contains(&Tag::Vaf1);
    if want_vaf || want_vaf1 {
        vaf_hdr = true;
    }

    // `-S` population map: name -> comma list; preserve first-seen order.
    let mut pop_order: Vec<String> = Vec::new();
    let mut sample_pops: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(sf) = opts.samples_file {
        let body = fs::read_to_string(sf)?;
        for line in body.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut it = line.split_whitespace();
            let (Some(sample), Some(pops)) = (it.next(), it.next()) else {
                continue;
            };
            let list: Vec<String> = pops
                .split(',')
                .filter(|p| !p.is_empty())
                .map(|p| {
                    let p = p.to_owned();
                    if !pop_order.contains(&p) {
                        pop_order.push(p.clone());
                    }
                    p
                })
                .collect();
            sample_pops.insert(sample.to_owned(), list);
        }
    }

    // Population array: named pops (file order) then the global "" pop.
    let mut pops: Vec<Pop> = pop_order
        .iter()
        .map(|n| Pop {
            name: n.clone(),
            suffix: format!("_{n}"),
        })
        .collect();
    pops.push(Pop {
        name: String::new(),
        suffix: String::new(),
    });
    let global = pops.len() - 1;

    let text = read_vcf_text(input)?;

    let mut out = String::with_capacity(text.len() + 4096);
    let mut sample_to_pops: Vec<Vec<usize>> = Vec::new();
    let mut gt_warned = false;
    let mut hwe_buf: Vec<f64> = Vec::new();
    // bcftools inserts a PASS FILTER header (after ##fileformat) when
    // the input lacks one.
    let has_pass = text.contains("##FILTER=<ID=PASS,");

    for line in text.lines() {
        if line.starts_with("##") {
            out.push_str(line);
            out.push('\n');
            if !has_pass && line.starts_with("##fileformat=") {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            }
            continue;
        }
        if let Some(cols) = line.strip_prefix("#CHROM") {
            // Append the new ##INFO header lines (IDs not already
            // declared), in upstream hdr_append × pop order.
            let declared: Vec<String> = collect_declared_info_ids(&out);
            let mut hdr = String::new();
            for &tag in HDR_ORDER {
                if !want.contains(&tag) {
                    continue;
                }
                for p in &pops {
                    let id = format!("{}{}", tag.base_id(), p.suffix);
                    if declared.iter().any(|d| d == &id) {
                        continue;
                    }
                    let (_, tmpl) = tag.header();
                    let in_part = if p.name.is_empty() {
                        String::new()
                    } else {
                        format!(" in {}", p.name)
                    };
                    hdr.push_str(&tmpl.replace("{S}", &p.suffix).replace("{IN}", &in_part));
                    hdr.push('\n');
                }
            }
            out.push_str(&hdr);
            if vaf_hdr && !out.contains("##FORMAT=<ID=VAF,") {
                out.push_str(
                    "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"The fraction of reads with alternate allele (nALT/nSumAll)\">\n",
                );
                out.push_str(
                    "##FORMAT=<ID=VAF1,Number=1,Type=Float,Description=\"The fraction of reads with alternate alleles (nSumALT/nSumAll)\">\n",
                );
            }
            // `Added by +fill-tags expression …` headers for the
            // function tags, unless an identically-defined field
            // already exists (upstream `update_hdr` check).
            for func in &funcs {
                let kw = if func.info { "INFO" } else { "FORMAT" };
                let number = match func.cnt {
                    Some(n) => n.to_string(),
                    None => ".".to_owned(),
                };
                let ty = if func.is_int { "Integer" } else { "Float" };
                for p in &pops {
                    let id = format!("{}{}", func.dst, p.suffix);
                    let prefix = format!("##{kw}=<ID={id},Number={number},Type={ty},");
                    if out.contains(&prefix) {
                        continue;
                    }
                    let in_part = if p.name.is_empty() {
                        String::new()
                    } else {
                        format!(" in {}", p.name)
                    };
                    let desc = func.raw.replace('"', "\\\"");
                    out.push_str(&format!(
                        "##{kw}=<ID={id},Number={number},Type={ty},Description=\"Added by +fill-tags expression {desc}{in_part}\">\n"
                    ));
                }
            }
            // Map each sample column to its pop indices.
            let samples: Vec<&str> = cols.split('\t').skip(9).collect();
            sample_to_pops = samples
                .iter()
                .map(|s| {
                    let mut v = Vec::new();
                    if let Some(list) = sample_pops.get(*s) {
                        for n in list {
                            if let Some(i) = pops.iter().position(|p| &p.name == n) {
                                v.push(i);
                            }
                        }
                    }
                    v.push(global); // always the global pop
                    v
                })
                .collect();
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        out.push_str(&process_record(
            line,
            &pops,
            &sample_to_pops,
            &want,
            &funcs,
            opts.drop_missing,
            want_vaf,
            want_vaf1,
            &mut hwe_buf,
            &mut gt_warned,
        ));
        out.push('\n');
    }

    Ok(out)
}

fn collect_declared_info_ids(header: &str) -> Vec<String> {
    let mut ids = Vec::new();
    for l in header.lines() {
        if let Some(rest) = l.strip_prefix("##INFO=<ID=") {
            let id: String = rest.chars().take_while(|&c| c != ',' && c != '>').collect();
            ids.push(id);
        }
    }
    ids
}

#[allow(clippy::too_many_arguments)]
fn process_record(
    line: &str,
    pops: &[Pop],
    sample_to_pops: &[Vec<usize>],
    want: &[Tag],
    funcs: &[Func],
    drop_missing: bool,
    want_vaf: bool,
    want_vaf1: bool,
    hwe_buf: &mut Vec<f64>,
    gt_warned: &mut bool,
) -> String {
    let mut f: Vec<&str> = line.split('\t').collect();
    if f.len() < 10 {
        return line.to_owned();
    }
    let n_allele = if f[4] == "." {
        1
    } else {
        1 + f[4].split(',').count()
    };

    // Locate GT in FORMAT. The genotype-derived tags need it, but the
    // FORMAT/VAF step (below) does not, so a GT-less record still flows
    // through (upstream `process_vaf_vaf1` is independent of GT).
    let gt_idx = f[8].split(':').position(|k| k == "GT");
    if gt_idx.is_none() && !*gt_warned {
        *gt_warned = true;
    }

    let mut info = f[7].to_owned();
    // Genotype-derived INFO tags: only when FORMAT/GT is present.
    if let Some(gt_idx) = gt_idx {
        // Per-pop, per-allele counts + ns.
        let mut counts: Vec<Vec<Counts>> = pops
            .iter()
            .map(|_| vec![Counts::default(); n_allele])
            .collect();
        let mut ns: Vec<i64> = vec![0; pops.len()];
        // F_MISSING bookkeeping: samples mapped to a pop, and those whose
        // GT contains any missing allele (`GT="mis"`).
        let mut npop_smpl: Vec<i64> = vec![0; pops.len()];
        let mut nmiss: Vec<i64> = vec![0; pops.len()];

        for (si, sample) in f[9..].iter().enumerate() {
            let gt = sample.split(':').nth(gt_idx).unwrap_or(".");
            let smpl_pops = sample_to_pops.get(si).map(Vec::as_slice).unwrap_or(&[]);
            let mut distinct: Vec<usize> = Vec::new();
            let mut nals = 0usize;
            let mut islots = 0usize;
            let mut any_missing_allele = false;
            for tok in gt.split(['/', '|']) {
                islots += 1;
                if tok == "." || tok.is_empty() {
                    any_missing_allele = true;
                    continue; // missing allele
                }
                let Ok(idx) = tok.parse::<usize>() else {
                    continue;
                };
                if idx >= n_allele {
                    continue;
                }
                nals += 1;
                if !distinct.contains(&idx) {
                    distinct.push(idx);
                }
            }
            for &pi in smpl_pops {
                npop_smpl[pi] += 1;
                if any_missing_allele {
                    nmiss[pi] += 1;
                }
            }
            if nals == 0 {
                continue; // missing genotype
            }
            let is_hom = distinct.len() == 1;
            // Upstream classification (BRANCH_INT): a partially-missing GT is
            // hemizygous, or — under `-d` — counted via `nac` (`is_half`);
            // a single-allele genotype is hemizygous.
            let (is_half, is_hemi) = if nals != islots {
                if drop_missing {
                    (true, false)
                } else {
                    (false, true)
                }
            } else if nals == 1 {
                (false, true)
            } else {
                (false, false)
            };

            for &pi in smpl_pops {
                for &a in &distinct {
                    let c = &mut counts[pi][a];
                    if is_half {
                        c.nac += 1;
                    } else if !is_hom {
                        c.nhet += 1;
                    } else if !is_hemi {
                        c.nhom += 2;
                    } else {
                        c.nhemi += 1;
                    }
                }
                ns[pi] += 1;
            }
        }

        // Build INFO key=value additions in the fixed write order.
        for &tag in WRITE_ORDER {
            if !want.contains(&tag) {
                continue;
            }
            for (pi, p) in pops.iter().enumerate() {
                let c = &counts[pi];
                let total = |a: usize| c[a].nhet + c[a].nhom + c[a].nhemi + c[a].nac;
                match tag {
                    Tag::FMissing => {
                        let denom = npop_smpl[pi];
                        let v = if denom == 0 {
                            0.0
                        } else {
                            nmiss[pi] as f64 / denom as f64
                        };
                        set_info(&mut info, &format!("F_MISSING{}", p.suffix), &fmt_float(v));
                    }
                    Tag::Hwe | Tag::ExcHet => {
                        let key = format!(
                            "{}{}",
                            if tag == Tag::Hwe { "HWE" } else { "ExcHet" },
                            p.suffix
                        );
                        if n_allele <= 1 {
                            del_info(&mut info, &key);
                            continue;
                        }
                        let mut nref_tot = c[0].nhom;
                        for ci in c.iter().take(n_allele) {
                            nref_tot += ci.nhet;
                        }
                        let vals: Vec<String> = (1..n_allele)
                            .map(|j| {
                                let nref = nref_tot - c[j].nhet;
                                let nalt = c[j].nhet + c[j].nhom;
                                let nhet = c[j].nhet;
                                let (hwe, exc) = if nref > 0 && nalt > 0 {
                                    calc_hwe(nref, nalt, nhet, hwe_buf)
                                } else {
                                    (1.0, 1.0)
                                };
                                fmt_float(if tag == Tag::Hwe { hwe } else { exc })
                            })
                            .collect();
                        set_info(&mut info, &key, &vals.join(","));
                    }
                    Tag::Ns => {
                        set_info(&mut info, &format!("NS{}", p.suffix), &ns[pi].to_string());
                    }
                    Tag::An => {
                        let an: i64 = (0..n_allele).map(total).sum();
                        set_info(&mut info, &format!("AN{}", p.suffix), &an.to_string());
                    }
                    // Number=A / Number=1-second tags: when there is no ALT
                    // allele upstream's `bcf_update_info_*` is called with 0
                    // values, which deletes the tag.
                    Tag::Ac => {
                        let key = format!("AC{}", p.suffix);
                        if n_allele > 1 {
                            let v: Vec<String> =
                                (1..n_allele).map(|a| total(a).to_string()).collect();
                            set_info(&mut info, &key, &v.join(","));
                        } else {
                            del_info(&mut info, &key);
                        }
                    }
                    Tag::AcHet => {
                        let key = format!("AC_Het{}", p.suffix);
                        if n_allele > 1 {
                            let v: Vec<String> =
                                (1..n_allele).map(|a| c[a].nhet.to_string()).collect();
                            set_info(&mut info, &key, &v.join(","));
                        } else {
                            del_info(&mut info, &key);
                        }
                    }
                    Tag::AcHom => {
                        let key = format!("AC_Hom{}", p.suffix);
                        if n_allele > 1 {
                            let v: Vec<String> =
                                (1..n_allele).map(|a| c[a].nhom.to_string()).collect();
                            set_info(&mut info, &key, &v.join(","));
                        } else {
                            del_info(&mut info, &key);
                        }
                    }
                    Tag::AcHemi => {
                        let key = format!("AC_Hemi{}", p.suffix);
                        if n_allele > 1 {
                            let v: Vec<String> =
                                (1..n_allele).map(|a| c[a].nhemi.to_string()).collect();
                            set_info(&mut info, &key, &v.join(","));
                        } else {
                            del_info(&mut info, &key);
                        }
                    }
                    Tag::Af | Tag::Maf => {
                        let key =
                            format!("{}{}", if tag == Tag::Af { "AF" } else { "MAF" }, p.suffix);
                        if n_allele <= 1 {
                            del_info(&mut info, &key);
                            continue;
                        }
                        let mut fr: Vec<f64> = (0..n_allele).map(|a| total(a) as f64).collect();
                        let an: f64 = fr.iter().sum();
                        let missing = an == 0.0;
                        if !missing {
                            for x in &mut fr {
                                *x /= an;
                            }
                        }
                        if tag == Tag::Af {
                            let v: Vec<String> = fr[1..]
                                .iter()
                                .map(|&x| if missing { ".".into() } else { fmt_float(x) })
                                .collect();
                            set_info(&mut info, &key, &v.join(","));
                        } else {
                            // MAF: second most common allele frequency.
                            if !missing {
                                fr.sort_by(|a, b| b.partial_cmp(a).unwrap());
                            }
                            let maf = if missing {
                                ".".into()
                            } else {
                                fmt_float(fr[1])
                            };
                            set_info(&mut info, &key, &maf);
                        }
                    }
                    // FORMAT tags — emitted by the dedicated VAF step below,
                    // never through the per-pop INFO loop.
                    Tag::Vaf | Tag::Vaf1 => {}
                }
            }
        }
    } // end `if let Some(gt_idx)` (genotype-derived INFO tags)

    // INFO function tags: `[int|float](sum(INFO/X | FMT/X))` → a scalar
    // written per pop (upstream `ftf_filter_expr`, `ftf->info`).
    for func in funcs.iter().filter(|f| f.info) {
        for (pi, p) in pops.iter().enumerate() {
            let v = match func.kind {
                FuncKind::FPass | FuncKind::NPass => {
                    let mut npass = 0i64;
                    let mut nsmpl = 0i64;
                    for si in 0..f.len().saturating_sub(9) {
                        if !sample_to_pops.get(si).is_some_and(|ps| ps.contains(&pi)) {
                            continue;
                        }
                        nsmpl += 1;
                        if sample_passes(&f, si, &func.pass_expr) {
                            npass += 1;
                        }
                    }
                    if func.kind == FuncKind::NPass {
                        Some(npass as f64)
                    } else if nsmpl > 0 {
                        Some(npass as f64 / nsmpl as f64)
                    } else {
                        Some(0.0)
                    }
                }
                _ => eval_sum(&f, func, pi, sample_to_pops),
            };
            set_info(
                &mut info,
                &format!("{}{}", func.dst, p.suffix),
                &fmt_num(v, func.is_int),
            );
        }
    }

    f[7] = &info;

    // Appended FORMAT columns: `smpl_sum(...)` function tags and
    // `FORMAT/VAF`+`VAF1`, applied in one pass.
    let mut add_cols: Vec<(String, Vec<String>)> = Vec::new();

    for func in funcs.iter().filter(|f| !f.info) {
        // `smpl_sum`: per-sample total (FORMAT result).
        let vals: Vec<String> = (0..f.len().saturating_sub(9))
            .map(|si| {
                let v = sample_operand_values(&f, si, func);
                if v.is_empty() {
                    ".".to_owned()
                } else {
                    fmt_num(Some(v.iter().sum()), func.is_int)
                }
            })
            .collect();
        if vals.iter().any(|v| v != ".") {
            add_cols.push((func.dst.clone(), vals));
        }
    }

    if (want_vaf || want_vaf1)
        && n_allele > 1
        && let Some(ad_idx) = f[8].split(':').position(|k| k == "AD")
    {
        let mut sample_vaf: Vec<(String, String)> = Vec::new();
        let mut any_valid = false;
        for col in &f[9..] {
            let ad = col.split(':').nth(ad_idx).unwrap_or(".");
            let parsed: Option<Vec<i64>> = if ad == "." || ad.is_empty() {
                None
            } else {
                let v: Vec<&str> = ad.split(',').collect();
                if v.len() == n_allele && v.iter().all(|x| x.parse::<i64>().is_ok()) {
                    Some(v.iter().map(|x| x.parse().unwrap()).collect())
                } else {
                    None
                }
            };
            match parsed {
                Some(ad) => {
                    any_valid = true;
                    let sum: i64 = ad.iter().sum();
                    let vaf = (1..n_allele)
                        .map(|j| {
                            if sum != 0 {
                                fmt_float(ad[j] as f64 / sum as f64)
                            } else {
                                "0".to_owned()
                            }
                        })
                        .collect::<Vec<_>>()
                        .join(",");
                    let vaf1 = if sum != 0 {
                        fmt_float((sum - ad[0]) as f64 / sum as f64)
                    } else {
                        "0".to_owned()
                    };
                    sample_vaf.push((vaf, vaf1));
                }
                None => sample_vaf.push((".".to_owned(), ".".to_owned())),
            }
        }
        if any_valid {
            if want_vaf {
                add_cols.push((
                    "VAF".to_owned(),
                    sample_vaf.iter().map(|(v, _)| v.clone()).collect(),
                ));
            }
            if want_vaf1 {
                add_cols.push((
                    "VAF1".to_owned(),
                    sample_vaf.iter().map(|(_, v)| v.clone()).collect(),
                ));
            }
        }
    }

    let new_format: String;
    let new_samples: Vec<String>;
    if !add_cols.is_empty() {
        let mut fmt = f[8].to_owned();
        for (name, _) in &add_cols {
            fmt.push(':');
            fmt.push_str(name);
        }
        new_format = fmt;
        new_samples = (0..f.len() - 9)
            .map(|si| {
                let mut s = f[9 + si].to_owned();
                for (_, vals) in &add_cols {
                    s.push(':');
                    s.push_str(&vals[si]);
                }
                s
            })
            .collect();
        f[8] = &new_format;
        for (col, ns) in f[9..].iter_mut().zip(new_samples.iter()) {
            *col = ns.as_str();
        }
    }

    // bcftools re-serializes the record; a lone `.` sample column
    // expands to one `.` per FORMAT subfield.
    let nkeys = f[8].split(':').count();
    let dots = vec!["."; nkeys].join(":");
    if nkeys > 1 {
        for col in &mut f[9..] {
            if *col == "." {
                *col = dots.as_str();
            }
        }
    }
    f.join("\t")
}

/// Wigginton 2005 (PMID:15789306) HWE exact test, ported from upstream
/// `calc_hwe`. Returns `(hwe, exc_het)`. `nref`/`nalt` are allele counts
/// (`nalt = nhet + nhom`, `nhom` already doubled per upstream `counts`).
fn calc_hwe(nref: i64, nalt: i64, nhet: i64, probs: &mut Vec<f64>) -> (f64, f64) {
    let ngt = (nref + nalt) / 2;
    let nrare = nref.min(nalt);
    // Upstream asserts these; on violation fall back to the neutral 1.0
    // rather than abort.
    if (nrare & 1) ^ (nhet & 1) != 0 || nrare < nhet || (nref + nalt) & 1 != 0 {
        return (1.0, 1.0);
    }
    let nrare_us = nrare as usize;
    probs.clear();
    probs.resize(nrare_us + 1, 0.0);

    let mut mid = ((nrare as f64) * ((nref + nalt - nrare) as f64) / ((nref + nalt) as f64)) as i64;
    if (nrare & 1) ^ (mid & 1) != 0 {
        mid += 1;
    }

    let mut hom_r = (nrare - mid) / 2;
    let mut hom_c = ngt - mid - hom_r;
    let mut sum = 1.0;
    probs[mid as usize] = 1.0;

    let mut het = mid;
    while het > 1 {
        probs[(het - 2) as usize] = probs[het as usize] * het as f64 * (het as f64 - 1.0)
            / (4.0 * (hom_r as f64 + 1.0) * (hom_c as f64 + 1.0));
        sum += probs[(het - 2) as usize];
        hom_r += 1;
        hom_c += 1;
        het -= 2;
    }

    het = mid;
    hom_r = (nrare - mid) / 2;
    hom_c = ngt - mid - hom_r;
    while het <= nrare - 2 {
        probs[(het + 2) as usize] = probs[het as usize] * 4.0 * hom_r as f64 * hom_c as f64
            / ((het as f64 + 2.0) * (het as f64 + 1.0));
        sum += probs[(het + 2) as usize];
        hom_r -= 1;
        hom_c -= 1;
        het += 2;
    }

    for p in probs.iter_mut() {
        *p /= sum;
    }

    let p_nhet = probs[nhet as usize];
    let mut exc = p_nhet;
    for h in (nhet + 1)..=nrare {
        exc += probs[h as usize];
    }

    let mut hwe = 0.0;
    for h in 0..=nrare {
        let ph = probs[h as usize];
        if ph > p_nhet {
            continue;
        }
        hwe += ph;
    }
    if hwe > 1.0 {
        hwe = 1.0;
    }
    (hwe, exc)
}

/// bcftools float printing: C `printf("%g")` with the default
/// precision 6 (6 significant digits; fixed unless exponent `< -4` or
/// `>= 6`, trailing zeros trimmed; min-2-digit signed exponent).
fn fmt_float(x: f64) -> String {
    // bcftools stores INFO/FORMAT floats as 32-bit; the printed text is
    // the f32 value, so round through f32 before formatting.
    let x = x as f32 as f64;
    if x == 0.0 {
        return "0".to_owned();
    }
    if !x.is_finite() {
        return if x.is_nan() {
            "nan".to_owned()
        } else if x < 0.0 {
            "-inf".to_owned()
        } else {
            "inf".to_owned()
        };
    }
    const P: i32 = 6;
    let exp = x.abs().log10().floor() as i32;
    if !(-4..P).contains(&exp) {
        // Scientific: P-1 mantissa decimals, C-style `e[+-]NN`.
        let s = format!("{:.*e}", (P - 1) as usize, x);
        let (mant, e) = s.split_once('e').unwrap();
        let mant = if mant.contains('.') {
            mant.trim_end_matches('0').trim_end_matches('.')
        } else {
            mant
        };
        let ev: i32 = e.parse().unwrap_or(0);
        return format!("{mant}e{}{:02}", if ev < 0 { '-' } else { '+' }, ev.abs());
    }
    let decimals = (P - 1 - exp).max(0) as usize;
    let s = format!("{x:.decimals$}");
    if s.contains('.') {
        let t = s.trim_end_matches('0').trim_end_matches('.');
        t.to_owned()
    } else {
        s
    }
}

/// `bcf_update_info_*`: replace `INFO/<key>` in place if present, else
/// append (preserving the rest of the INFO column).
fn set_info(info: &mut String, key: &str, value: &str) {
    let entry = format!("{key}={value}");
    if *info == "." || info.is_empty() {
        *info = entry;
        return;
    }
    let mut kept: Vec<String> = Vec::new();
    let mut replaced = false;
    for kv in info.split(';') {
        let k = kv.split('=').next().unwrap_or("");
        if k == key {
            kept.push(entry.clone());
            replaced = true;
        } else {
            kept.push(kv.to_owned());
        }
    }
    if !replaced {
        kept.push(entry);
    }
    *info = kept.join(";");
}

/// `bcf_update_info_*` with zero values: remove `INFO/<key>` if present.
fn del_info(info: &mut String, key: &str) {
    if *info == "." || info.is_empty() {
        return;
    }
    let kept: Vec<&str> = info
        .split(';')
        .filter(|kv| kv.split('=').next() != Some(key))
        .collect();
    *info = if kept.is_empty() {
        ".".to_owned()
    } else {
        kept.join(";")
    };
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
        ".bcftools-rs-fill-tags-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_formatting() {
        assert_eq!(fmt_float(0.5), "0.5");
        assert_eq!(fmt_float(1.0 / 3.0), "0.333333");
        assert_eq!(fmt_float(2.0 / 3.0), "0.666667");
        assert_eq!(fmt_float(0.25), "0.25");
        assert_eq!(fmt_float(1.0), "1");
        assert_eq!(fmt_float(0.0), "0");
    }

    #[test]
    fn set_info_replace_and_append() {
        let mut i = "DP=0".to_owned();
        set_info(&mut i, "AN", "4");
        assert_eq!(i, "DP=0;AN=4");
        set_info(&mut i, "DP", "9");
        assert_eq!(i, "DP=9;AN=4");
        let mut j = ".".to_owned();
        set_info(&mut j, "AC", "2");
        assert_eq!(j, "AC=2");
    }

    #[test]
    fn tag_parsing_aliases() {
        assert!(matches!(parse_tag("AC_Hom"), Ok(Tag::AcHom)));
        assert!(matches!(parse_tag("INFO/AF"), Ok(Tag::Af)));
        assert!(matches!(parse_tag("HWE"), Ok(Tag::Hwe)));
        assert!(matches!(parse_tag("ExcHet"), Ok(Tag::ExcHet)));
        assert!(matches!(parse_tag("F_MISSING"), Ok(Tag::FMissing)));
        // The function engine / END / TYPE remain unported.
        assert!(parse_tag("TYPE").is_err());
    }
}
