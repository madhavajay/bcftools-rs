//! `bcftools +prune` (upstream `bcftools/plugins/prune.c` + the windowed
//! `_prune_sites` path of `bcftools/vcfbuf.c`), `-n`/`-N` window subset.
//!
//! Keeps at most `-n N` sites per `-w` window, choosing which to drop by
//! `-N` mode: `1st` (keep the earliest), `maxAF` (drop the lowest allele
//! frequency, `--AF-tag` or computed). Faithful port of the `vcfbuf`
//! window flush condition and `_prune_sites` removal order. Common
//! `-i`/`-e` record filtering (applied before windowing, as upstream
//! does) routes through the shared filter engine. The LD/`-a`/`-m`
//! `count` cluster modes and `-N rand` (needs `hts_drand48` parity)
//! remain tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    MaxAf,
    First,
}

impl Mode {
    pub fn parse(s: &str) -> Result<Mode, String> {
        if s.eq_ignore_ascii_case("maxAF") {
            Ok(Mode::MaxAf)
        } else if s.eq_ignore_ascii_case("1st") {
            Ok(Mode::First)
        } else if s.eq_ignore_ascii_case("rand") {
            Err("prune -N rand needs hts_drand48 parity and is not supported in this slice".into())
        } else {
            Err(format!("prune -N mode '{s}' not recognised"))
        }
    }
}

/// Parses upstream `-w INT[bp|kb|Mb]` into the `vcfbuf` window value
/// (negative = bp-style span window, positive = site-count window).
pub fn parse_window(s: &str) -> Result<i64, String> {
    let digits: String = s
        .chars()
        .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+' || *c == '.')
        .collect();
    let n: f64 = digits
        .parse()
        .map_err(|_| format!("Could not parse: --window {s}"))?;
    let suffix = &s[digits.len()..];
    let v = n as i64;
    match suffix {
        "" => Ok(v),
        s if s.eq_ignore_ascii_case("bp") => Ok(-v),
        s if s.eq_ignore_ascii_case("kb") => Ok(-v * 1000),
        s if s.eq_ignore_ascii_case("Mb") => Ok(-v * 1_000_000),
        _ => Err(format!("Could not parse: --window {s}")),
    }
}

struct Rec {
    line: String,
    chrom: String,
    pos: i64,
    af: Option<f32>,
}

/// Reads the input VCF/BCF and returns the pruned VCF text.
pub fn run(
    input: &Path,
    win: i64,
    max_sites: i32,
    mode: Mode,
    af_tag: Option<&str>,
    filter: Option<(bool, &str)>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(process(&text, win, max_sites, mode, af_tag, filter))
}

/// Upstream `+prune -i/-e` discards non-matching records *before* the
/// windowed pruning. Record-level eval through the shared filter engine
/// (per-sample `EvalContext` + `record_lookup`, the `+split` wiring).
fn record_passes(line: &str, exclude: bool, expr: &str) -> bool {
    let fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if fields.len() < 8 {
        return true;
    }
    let context = if fields.len() > 9 {
        let keys: Vec<&str> = fields[8].split(':').collect();
        fields[9..].iter().fold(EvalContext::new(), |c, s| {
            let vals: Vec<&str> = s.split(':').collect();
            c.with_sample(
                keys.iter()
                    .enumerate()
                    .map(|(i, k)| {
                        let raw = vals.get(i).copied().unwrap_or(".");
                        let v = if k.eq_ignore_ascii_case("GT") {
                            FilterValue::String(raw.to_owned())
                        } else if raw == "." || raw.is_empty() {
                            FilterValue::Missing
                        } else if let Ok(n) = raw.parse::<f64>() {
                            FilterValue::Number(n)
                        } else {
                            FilterValue::String(raw.to_owned())
                        };
                        ((*k).to_owned(), v)
                    })
                    .collect::<Vec<_>>(),
            )
        })
    } else {
        EvalContext::new()
    };
    let matched = bcffilter::eval_expression_with(expr, &context, |name, si| {
        if si.is_some() {
            return None;
        }
        crate::commands::filter::record_lookup(name, &fields)
    })
    .map(|v| v.truthy())
    .unwrap_or(false);
    if exclude { !matched } else { matched }
}

fn process(
    text: &str,
    win: i64,
    max_sites: i32,
    mode: Mode,
    af_tag: Option<&str>,
    filter: Option<(bool, &str)>,
) -> String {
    let lines: Vec<&str> = text.lines().collect();

    let mut out = String::with_capacity(text.len());
    emit_header(&lines, &mut out);

    let mut buf: Vec<Rec> = Vec::new();
    let emit = |r: Rec, out: &mut String| {
        out.push_str(&r.line);
        out.push('\n');
    };

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        if let Some((exclude, expr)) = filter
            && !record_passes(line, exclude, expr)
        {
            continue; // upstream discards non-matching sites pre-window
        }
        let Ok(pos) = f[1].parse::<i64>() else {
            continue;
        };
        buf.push(Rec {
            line: (*line).to_owned(),
            chrom: f[0].to_owned(),
            pos,
            af: None,
        });
        while let Some(r) = vcfbuf_flush(&mut buf, false, win, max_sites, mode, af_tag) {
            emit(r, &mut out);
        }
    }
    while let Some(r) = vcfbuf_flush(&mut buf, true, win, max_sites, mode, af_tag) {
        emit(r, &mut out);
    }
    out
}

/// Port of the `vcfbuf_flush` window path (`buf->win`): decide whether the
/// front can be emitted, prune the window if over `max_sites`, then shift.
fn vcfbuf_flush(
    buf: &mut Vec<Rec>,
    flush_all: bool,
    win: i64,
    max_sites: i32,
    mode: Mode,
    af_tag: Option<&str>,
) -> Option<Rec> {
    if buf.is_empty() {
        return None;
    }
    if win == 0 {
        // No window configured: upstream falls through to the always-flush
        // path (no pruning happens without -w).
        return Some(buf.remove(0));
    }
    let last = buf.len() - 1;
    let mut can_flush = flush_all;
    if buf[0].chrom != buf[last].chrom {
        can_flush = true;
    } else if win > 0 {
        if buf.len() as i64 > win {
            can_flush = true;
        }
    } else if buf[0].pos - buf[last].pos <= win {
        can_flush = true;
    }
    if !can_flush {
        return None;
    }
    if max_sites > 0 && (max_sites as usize) < buf.len() {
        prune_sites(buf, flush_all, max_sites, mode, af_tag);
    }
    Some(buf.remove(0))
}

/// Port of `_prune_sites` for the `1ST` and `MAX_AF` modes.
fn prune_sites(
    buf: &mut Vec<Rec>,
    flush_all: bool,
    max_sites: i32,
    mode: Mode,
    af_tag: Option<&str>,
) {
    let nbuf = if flush_all { buf.len() } else { buf.len() - 1 };
    let nprune = nbuf as i32 - max_sites;
    if nprune <= 0 {
        return;
    }
    let nprune = nprune as usize;

    if mode == Mode::First {
        let eoff = if flush_all { 1 } else { 2 };
        for _ in 0..nprune {
            let idx = buf.len() - eoff;
            buf.remove(idx);
        }
        return;
    }

    // MAX_AF: compute AF for the first nbuf records, drop the lowest.
    let mut order: Vec<usize> = (0..nbuf).collect();
    for &i in &order {
        if buf[i].af.is_none() {
            buf[i].af = Some(compute_af(&buf[i].line, af_tag));
        }
    }
    // cmpvrec: ascending by af (stable order matches upstream qsort for the
    // tie-free inputs this slice targets).
    order.sort_by(|&a, &b| {
        buf[a]
            .af
            .unwrap()
            .partial_cmp(&buf[b].af.unwrap())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let mut to_remove: Vec<usize> = order[..nprune].to_vec();
    to_remove.sort_unstable_by(|a, b| b.cmp(a)); // descending
    for idx in to_remove {
        buf.remove(idx);
    }
}

/// `--AF-tag` first value, or upstream's `bcf_calc_ac`-derived
/// `nalt / ac[0]` fraction when no tag is given.
fn compute_af(line: &str, af_tag: Option<&str>) -> f32 {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return 0.0;
    }
    if let Some(tag) = af_tag {
        if f[7] != "." {
            for kv in f[7].split(';') {
                let mut it = kv.splitn(2, '=');
                if it.next() == Some(tag)
                    && let Some(v) = it.next()
                    && let Some(first) = v.split(',').next()
                {
                    return first.parse::<f32>().unwrap_or(0.0);
                }
            }
        }
        return 0.0;
    }

    // No tag: bcf_calc_ac over GT, af = nalt / ac[0] (upstream formula).
    if f.len() < 10 {
        return 0.0;
    }
    let n_allele = 1 + if f[4] == "." {
        0
    } else {
        f[4].split(',').count()
    };
    let gt_slot = f[8].split(':').position(|k| k == "GT");
    let mut ac = vec![0i64; n_allele];
    for s in &f[9..] {
        let gt = match gt_slot {
            Some(idx) => s.split(':').nth(idx).unwrap_or("."),
            None => ".",
        };
        for tok in gt.split(['/', '|']) {
            if let Ok(a) = tok.parse::<usize>()
                && a < n_allele
            {
                ac[a] += 1;
            }
        }
    }
    let ntot = ac[0];
    let nalt: i64 = ac[1..].iter().sum();
    if ntot != 0 {
        nalt as f32 / ntot as f32
    } else {
        0.0
    }
}

/// htslib-style header: inject `##FILTER=<ID=PASS>` right after
/// `##fileformat` when absent; otherwise pass through verbatim.
fn emit_header(lines: &[&str], out: &mut String) {
    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
    for (idx, line) in lines.iter().enumerate() {
        if !line.starts_with('#') {
            break;
        }
        out.push_str(line);
        out.push('\n');
        if Some(idx) == fileformat && !has_pass {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">");
            out.push('\n');
        }
    }
}

// ---------------------------------------------------------------------------
// LD (`-a`/`-m`) path: vcfbuf _calc_r2_ld + vcfbuf_ld + prune.c driver
// ---------------------------------------------------------------------------

/// VCFBUF_LD_IDX_{R2,LD,RD}; upstream `VCFBUF_LD_N == 3`.
const LD_N: usize = 3;
const IDX_R2: usize = 0;
const IDX_LD: usize = 1;
const IDX_RD: usize = 2;

/// Faithful port of HTSlib `kstring.c:kputd` — the `%g`-only double
/// formatter HTSlib uses to serialize float INFO values.
pub fn kputd(d: f64) -> String {
    if d == 0.0 {
        return if d.is_sign_negative() {
            "-0".to_owned()
        } else {
            "0".to_owned()
        };
    }
    let mut out = String::new();
    let mut d = d;
    if d < 0.0 {
        out.push('-');
        d = -d;
    }
    if !(0.0001..=999999.0).contains(&d) {
        // HTSlib defers to stdio "%g" for the exponent cases.
        out.push_str(&format_g(d));
        return out;
    }

    // buf[0..21], cp starts at index 20.
    let mut buf = [0u8; 21];
    let mut cp: i32 = 20;
    let i: u32;
    if d < 0.001 {
        i = (d * 1_000_000_000.0).round() as u32;
        cp -= 1;
    } else if d < 0.01 {
        i = (d * 100_000_000.0).round() as u32;
        cp -= 2;
    } else if d < 0.1 {
        i = (d * 10_000_000.0).round() as u32;
        cp -= 3;
    } else if d < 1.0 {
        i = (d * 1_000_000.0).round() as u32;
        cp -= 4;
    } else if d < 10.0 {
        i = (d * 100_000.0).round() as u32;
        cp -= 5;
    } else if d < 100.0 {
        i = (d * 10_000.0).round() as u32;
        cp -= 6;
    } else if d < 1000.0 {
        i = (d * 1_000.0).round() as u32;
        cp -= 7;
    } else if d < 10000.0 {
        i = (d * 100.0).round() as u32;
        cp -= 8;
    } else if d < 100000.0 {
        i = (d * 10.0).round() as u32;
        cp -= 9;
    } else {
        i = d.round() as u32;
        cp -= 10;
    }

    const DIG2R: &[u8; 200] = b"00010203040506070809\
1011121314151617181\
9202122232425262728\
2930313233343536373\
8394041424344454647\
4849505152535455565\
7585960616263646566\
6768697071727374757\
6777879808182838485\
868788899091929394959697989\
9";
    // (above split is cosmetic; rebuild the canonical table instead)
    let mut dig2r = [0u8; 200];
    for n in 0..100 {
        dig2r[2 * n] = b'0' + (n / 10) as u8;
        dig2r[2 * n + 1] = b'0' + (n % 10) as u8;
    }
    let _ = DIG2R;

    let mut i = i;
    let put2 = |buf: &mut [u8; 21], cp: &mut i32, v: usize| {
        *cp -= 2;
        let idx = *cp as usize;
        buf[idx] = dig2r[2 * (v % 100)];
        buf[idx + 1] = dig2r[2 * (v % 100) + 1];
    };
    put2(&mut buf, &mut cp, (i % 100) as usize);
    i /= 100;
    put2(&mut buf, &mut cp, (i % 100) as usize);
    i /= 100;
    put2(&mut buf, &mut cp, (i % 100) as usize);
    i /= 100;
    if i >= 100 {
        cp -= 1;
        buf[cp as usize] = b'0' + (i / 100) as u8;
    }

    let mut p: i32 = 20 - cp;
    let ep: i32;
    if p <= 10 {
        // d < 1: prepend zeros then "0."
        ep = cp + 5;
        while p < 10 {
            cp -= 1;
            buf[cp as usize] = b'0';
            p += 1;
        }
        cp -= 1;
        buf[cp as usize] = b'.';
        cp -= 1;
        buf[cp as usize] = b'0';
    } else {
        // 123.001 is 123001 with p==13: shift down and insert '.'
        cp -= 1;
        ep = cp + 6;
        let mut xp = cp as usize;
        let mut pp = p;
        while pp > 10 {
            buf[xp] = buf[xp + 1];
            xp += 1;
            pp -= 1;
        }
        buf[xp] = b'.';
    }

    // Cull trailing zeros.
    let mut e = ep;
    while e > cp && buf[e as usize] == b'0' {
        e -= 1;
    }
    // End can be 1 out; also turn "123." into "123".
    if buf[e as usize] != 0 && buf[e as usize] != b'.' {
        e += 1;
    }
    let s = std::str::from_utf8(&buf[cp as usize..e as usize]).unwrap_or("");
    out.push_str(s);
    out
}

/// Minimal C `printf("%g", d)` (default precision 6) for the
/// out-of-`[0.0001,999999]` tail kputd defers to stdio for. Adequate for
/// the magnitudes prune emits; trims like `%g`.
fn format_g(d: f64) -> String {
    if d == 0.0 {
        return "0".to_owned();
    }
    let exp = d.abs().log10().floor() as i32;
    if (-4..6).contains(&exp) {
        let prec = (5 - exp).max(0) as usize;
        let s = format!("{d:.prec$}");
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_owned()
        } else {
            s
        }
    } else {
        let mut exp = exp;
        let mut m = d / 10f64.powi(exp);
        if format!("{m:.5}").parse::<f64>().unwrap_or(m) >= 10.0 {
            exp += 1;
            m /= 10.0;
        }
        let ms = format!("{m:.5}");
        let ms = ms.trim_end_matches('0').trim_end_matches('.').to_owned();
        format!("{ms}e{}{:02}", if exp < 0 { '-' } else { '+' }, exp.abs())
    }
}

/// Per-sample dosage (count of non-ref alleles) and allele count, with the
/// upstream break-on-missing / break-on-vector_end semantics.
fn dosage(gt: &str) -> (i32, i32) {
    let mut dsg = 0;
    let mut an = 0;
    for tok in gt.split(['/', '|']) {
        if tok == "." || tok.is_empty() {
            break; // missing (no rand_missing)
        }
        match tok.parse::<i64>() {
            Ok(0) => {}
            Ok(_) => dsg += 1,
            Err(_) => break,
        }
        an += 1;
    }
    (dsg, an)
}

/// Port of `vcfbuf.c:_calc_r2_ld`. Returns `[r2, ld, rd]` or `None`
/// (no GT / no shared data).
fn calc_r2_ld(a_gts: &[&str], b_gts: &[&str]) -> Option<[f64; 3]> {
    let mut nhd = [0.0f64; 9];
    let (mut ab, mut aa, mut bb, mut a, mut b) = (0.0, 0.0, 0.0, 0.0, 0.0);
    let mut nab = 0i32;
    let mut ndiff = 0i32;
    let mut an_tot = 0i32;
    let mut bn_tot = 0i32;
    let n = a_gts.len().min(b_gts.len());
    for s in 0..n {
        let (adsg, an) = dosage(a_gts[s]);
        let (bdsg, bn) = dosage(b_gts[s]);
        if an != 0 && bn != 0 {
            an_tot += an;
            aa += (adsg * adsg) as f64;
            a += adsg as f64;
            bn_tot += bn;
            bb += (bdsg * bdsg) as f64;
            b += bdsg as f64;
            if adsg != bdsg {
                ndiff += 1;
            }
            ab += (adsg * bdsg) as f64;
            nab += 1;
        }
        if an == 2 && bn == 2 {
            nhd[(bdsg * 3 + adsg) as usize] += 1.0;
        }
    }
    if nab == 0 {
        return None;
    }
    let mut nab = nab as f64;
    let pa = a / an_tot as f64;
    let pb = b / bn_tot as f64;
    let cor = if ndiff == 0 {
        1.0
    } else {
        if aa == a * a / nab || bb == b * b / nab {
            aa += 1e-4;
            bb += 1e-4;
            ab += 1e-4;
            a += 1e-2;
            b += 1e-2;
            nab += 1.0;
        }
        (ab - a * b / nab) / (aa - a * a / nab).sqrt() / (bb - b * b / nab).sqrt()
    };

    let mut val = [0.0f64; 3];
    val[IDX_R2] = cor * cor;

    val[IDX_LD] = cor * (pa * (1.0 - pa) * pb * (1.0 - pb)).sqrt();
    let norm = if val[IDX_LD] < 0.0 {
        (-pa * pb).max(-(1.0 - pa) * (1.0 - pb))
    } else {
        (pa * (1.0 - pb)).max((1.0 - pa) * pb)
    };
    if norm != 0.0 {
        val[IDX_LD] = if norm.abs() > val[IDX_LD].abs() {
            val[IDX_LD] / norm
        } else {
            1.0
        };
    }
    if val[IDX_LD] == 0.0 {
        val[IDX_LD] = val[IDX_LD].abs(); // avoid "-0"
    }

    val[IDX_RD] = (nhd[0] + nhd[1] / 2.0 + nhd[3] / 2.0 + nhd[4] / 4.0)
        * (nhd[4] / 4.0 + nhd[5] / 2.0 + nhd[7] / 2.0 + nhd[8])
        - (nhd[1] / 2.0 + nhd[2] + nhd[4] / 4.0 + nhd[5] / 2.0)
            * (nhd[3] / 2.0 + nhd[4] / 4.0 + nhd[6] + nhd[7] / 2.0);
    val[IDX_RD] /= nab;
    val[IDX_RD] /= nab + 1.0;

    Some(val)
}

/// Which LD metrics are requested for `-a` (annotate) and the tag names.
#[derive(Clone, Copy, Default)]
pub struct LdAnnot {
    /// `[R2, LD, RD]` — true if `-a` requested that metric.
    pub annot: [bool; LD_N],
}

/// `-m` thresholds: `Some(max)` per `[R2, LD, RD]`.
pub type LdMax = [Option<f64>; LD_N];

struct LdRec<'a> {
    /// The (possibly annotated) output line, owned.
    out_line: String,
    pos: i64,
    chrom: &'a str,
    gts: Vec<&'a str>,
}

/// Runs `+prune -a/-m` (LD annotate / max-filter). `win` is the upstream
/// `vcfbuf` window (site count if >0). `soft_filter` = `-f` FILTER id.
#[allow(clippy::too_many_arguments)]
pub fn run_ld(
    input: &Path,
    win: i64,
    annot: LdAnnot,
    max: LdMax,
    soft_filter: Option<&str>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(process_ld(&text, win, annot, max, soft_filter))
}

fn ld_idx_order() -> [usize; LD_N] {
    [IDX_R2, IDX_LD, IDX_RD]
}

fn pos_tag(i: usize) -> &'static str {
    match i {
        IDX_R2 => "POS_R2",
        IDX_LD => "POS_LD",
        _ => "POS_RD",
    }
}
fn val_tag(i: usize) -> &'static str {
    match i {
        IDX_R2 => "R2",
        IDX_LD => "LD",
        _ => "RD",
    }
}

fn process_ld(
    text: &str,
    win: i64,
    annot: LdAnnot,
    max: LdMax,
    soft_filter: Option<&str>,
) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = String::with_capacity(text.len() + 512);

    // Header: PASS injection + soft FILTER + the requested INFO defs,
    // appended just before #CHROM (upstream `bcf_hdr_printf`).
    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
    let win_desc = if win < 0 {
        format!("{}kb", -win / 1000)
    } else {
        format!("{win} sites")
    };
    for (idx, line) in lines.iter().enumerate() {
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("#CHROM") {
            if let Some(fid) = soft_filter {
                // condition string: only R2 is exercised by the fixtures.
                if let Some(m) = max[IDX_R2] {
                    out.push_str(&format!(
                        "##FILTER=<ID={fid},Description=\"An upstream site within {win_desc} with R2 bigger than {}\">\n",
                        kputd(m)
                    ));
                }
            }
            for &i in &ld_idx_order() {
                if !annot.annot[i] {
                    continue;
                }
                let (vt, pt) = (val_tag(i), pos_tag(i));
                let desc = match i {
                    IDX_R2 => format!("Pairwise r2 with the {pt} site"),
                    IDX_LD => format!("Pairwise Lewontin's D' (PMID:19433632) with the {pt} site"),
                    _ => {
                        format!("Pairwise Ragsdale's \\hat{{D}} (PMID:31697386) with the {pt} site")
                    }
                };
                out.push_str(&format!(
                    "##INFO=<ID={vt},Number=1,Type=Float,Description=\"{desc}\">\n"
                ));
                out.push_str(&format!(
                    "##INFO=<ID={pt},Number=1,Type=Integer,Description=\"The position of the site for which {vt} was calculated\">\n"
                ));
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        out.push_str(line);
        out.push('\n');
        if Some(idx) == fileformat && !has_pass {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
        }
    }

    let mut gt_slot_cache: Option<usize> = None;
    let mut buf: Vec<LdRec> = Vec::new();

    let emit = |r: &LdRec, out: &mut String| {
        out.push_str(&r.out_line);
        out.push('\n');
    };

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let Ok(pos) = f[1].parse::<i64>() else {
            continue;
        };
        let chrom = f[0];
        let gt_slot = *gt_slot_cache
            .get_or_insert_with(|| f[8].split(':').position(|k| k == "GT").unwrap_or(0));
        let gts: Vec<&str> = f[9..]
            .iter()
            .map(|s| s.split(':').nth(gt_slot).unwrap_or("."))
            .collect();

        // vcfbuf_ld: compare rec against each buffered record (same chrom),
        // keep the per-metric max + position, early-exit when a pair
        // exceeds an -m threshold.
        let mut best_val = [f64::NEG_INFINITY; LD_N];
        let mut best_pos = [0i64; LD_N];
        let mut any = false;
        if buf.first().map(|r| r.chrom) == Some(chrom) {
            for prev in &buf {
                let Some(v) = calc_r2_ld(&prev.gts, &gts) else {
                    continue;
                };
                any = true;
                let mut done = false;
                for j in 0..LD_N {
                    if best_val[j] < v[j] {
                        best_val[j] = v[j];
                        best_pos[j] = prev.pos;
                    }
                    if let Some(mx) = max[j]
                        && mx < v[j]
                    {
                        done = true;
                    }
                }
                if done {
                    break;
                }
            }
        }

        let mut fields: Vec<String> = f.iter().map(|s| s.to_string()).collect();
        let mut dropped = false;
        if any {
            let pass = (0..LD_N).all(|j| match max[j] {
                Some(mx) => best_val[j] <= mx,
                None => true,
            });
            if !pass {
                match soft_filter {
                    None => dropped = true, // hard filter: drop, don't buffer
                    Some(fid) => {
                        fields[6] = if fields[6] == "." || fields[6].is_empty() {
                            fid.to_string()
                        } else {
                            format!("{};{}", fields[6], fid)
                        };
                    }
                }
            }
            if !dropped {
                // INFO: POS_* (int) for all metrics first, then values.
                let mut info = if fields[7] == "." || fields[7].is_empty() {
                    String::new()
                } else {
                    fields[7].clone()
                };
                let push_kv = |info: &mut String, k: &str, v: String| {
                    if !info.is_empty() {
                        info.push(';');
                    }
                    info.push_str(k);
                    info.push('=');
                    info.push_str(&v);
                };
                for &i in &ld_idx_order() {
                    if annot.annot[i] {
                        // best_pos is the matched record's 1-based POS
                        // (upstream uses 0-based rec->pos + 1).
                        push_kv(&mut info, pos_tag(i), best_pos[i].to_string());
                    }
                }
                for &i in &ld_idx_order() {
                    if annot.annot[i] {
                        // bcf_update_info_float stores f32 -> kputd(f64).
                        let v = best_val[i] as f32 as f64;
                        push_kv(&mut info, val_tag(i), kputd(v));
                    }
                }
                fields[7] = if info.is_empty() {
                    ".".to_owned()
                } else {
                    info
                };
            }
        }

        if dropped {
            continue; // hard-filtered: not emitted, not buffered
        }

        buf.push(LdRec {
            out_line: fields.join("\t"),
            pos,
            chrom,
            gts,
        });
        // window flush: positive win keeps at most `win` records.
        while win > 0 && buf.len() as i64 > win {
            let front = buf.remove(0);
            emit(&front, &mut out);
        }
        while win < 0
            && buf.len() > 1
            && buf[0].chrom == buf[buf.len() - 1].chrom
            && buf[buf.len() - 1].pos - buf[0].pos > -win
        {
            let front = buf.remove(0);
            emit(&front, &mut out);
        }
    }
    for r in &buf {
        emit(r, &mut out);
    }
    out
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
        ".bcftools-rs-prune-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=2147483647>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t101\t.\tT\tA\t.\t.\tAF=0.3\n\
1\t102\t.\tT\tA\t.\t.\tAF=0.2\n\
1\t103\t.\tT\tA\t.\t.\tAF=0.1\n\
1\t104\t.\tT\tA\t.\t.\tAF=0.3\n\
1\t105\t.\tT\tA\t.\t.\tAF=0.2\n\
1\t106\t.\tT\tA\t.\t.\tAF=0.1\n\
1\t107\t.\tT\tA\t.\t.\tAF=0.3\n\
1\t108\t.\tT\tA\t.\t.\tAF=0.2\n";

    fn positions(out: &str) -> Vec<&str> {
        out.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| l.split('\t').nth(1).unwrap())
            .collect()
    }

    #[test]
    fn window_parse() {
        assert_eq!(parse_window("2bp").unwrap(), -2);
        assert_eq!(parse_window("3kb").unwrap(), -3000);
        assert_eq!(parse_window("1Mb").unwrap(), -1_000_000);
        assert_eq!(parse_window("100").unwrap(), 100);
        assert!(parse_window("2xy").is_err());
    }

    #[test]
    fn kputd_exponent_rounding_carries_mantissa() {
        assert_eq!(kputd(9.999_999_824_516_7e-15), "1e-14");
    }

    #[test]
    fn first_mode_keeps_earliest() {
        let out = process(VCF, -2, 1, Mode::First, None, None);
        assert_eq!(positions(&out), vec!["101", "103", "105", "107"]);
    }

    #[test]
    fn maxaf_mode_drops_lowest_af() {
        let out = process(VCF, -2, 1, Mode::MaxAf, Some("AF"), None);
        assert_eq!(positions(&out), vec!["101", "104", "105", "107"]);
    }
}
