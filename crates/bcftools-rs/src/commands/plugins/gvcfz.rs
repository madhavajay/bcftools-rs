//! `bcftools +gvcfz` (upstream `bcftools/plugins/gvcfz.c`).
//!
//! Compresses a gVCF by resizing reference blocks according to `-g`
//! group expressions. Each record is assigned to the first group whose
//! expression matches (evaluated through the shared filter engine, the
//! same wiring `+split` uses); consecutive records in the same group are
//! merged into one block whose `INFO/END` is the max end, `FORMAT/DP`
//! the min depth, `FORMAT/GQ`|`RGQ` the min genotype quality, and
//! `FORMAT/PL` the element-wise min. The block's representative record is
//! the first record; non-`PASS` groups stamp their name into `FILTER`.
//! Records that are not gVCF blocks after optional `-a` allele trimming
//! are flushed and emitted unchanged.
//!
//! Text port: single-sample gVCF input, BCF/VCF/VCF.gz via the shared
//! text-view path. Validated byte-for-byte through `bcftools query`
//! against `gvcfz.{1,2}.out` and `gvcfz.2.1.out`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::filter::{self as bcffilter, EvalContext, Value as FilterValue};
use crate::vcf_compat::normalize_vcf_text;

/// One `-g` group: a FILTER label plus an optional include expression
/// (`None` for the literal `-`, which always matches).
struct Group {
    /// FILTER label; `is_pass` short-circuits the FILTER stamp.
    name: String,
    is_pass: bool,
    expr: Option<String>,
}

/// Parse `FLT:expr; FLT:expr; ...` exactly like upstream `init_groups`.
fn parse_groups(spec: &str) -> io::Result<Vec<Group>> {
    let mut groups = Vec::new();
    let mut rest = spec;
    loop {
        let beg = rest.trim_start();
        if beg.is_empty() {
            break;
        }
        let colon = beg.find(':').ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Could not parse the expression: \"{spec}\""),
            )
        })?;
        let flt = &beg[..colon];
        let after = &beg[colon + 1..];
        let (raw_expr, tail) = match after.find(';') {
            Some(semi) => (&after[..semi], Some(&after[semi + 1..])),
            None => (after, None),
        };
        let expr = raw_expr.trim();
        groups.push(Group {
            name: flt.to_owned(),
            is_pass: flt == "PASS",
            expr: if expr == "-" {
                None
            } else {
                Some(expr.to_owned())
            },
        });
        match tail {
            Some(t) => rest = t,
            None => break,
        }
    }
    Ok(groups)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum GqKey {
    None,
    Gq,
    Rgq,
}

struct Block {
    /// Representative record fields (first record of the block).
    fields: Vec<String>,
    rid: String,
    /// 1-based stop position.
    end: i64,
    /// 0-based start (`rec->pos`) of the representative.
    pos0: i64,
    min_dp: i64,
    gq_key: GqKey,
    gq: i64,
    pl: [i64; 3],
    /// Index into the group list (`== groups.len()` when nothing matched).
    grp: usize,
}

pub struct Options<'a> {
    pub group_by: &'a str,
    pub trim_alts: bool,
    /// Common site filter `(exclude, expr)`; `exclude` inverts the test.
    pub site_filter: Option<(bool, &'a str)>,
}

/// Reads the input gVCF and returns the block-compressed VCF text.
pub fn run(input: &Path, opts: Options<'_>) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let groups = parse_groups(opts.group_by)?;
    let ngrp = groups.len();

    let mut out = String::with_capacity(text.len());
    let mut block: Option<Block> = None;
    let mut header_done = false;

    for line in text.lines() {
        if line.starts_with('#') {
            if line.starts_with("#CHROM") && !header_done {
                inject_headers(&mut out, &groups);
                header_done = true;
            }
            if line.starts_with("##INFO=<ID=END,") {
                continue; // re-inserted by inject_headers
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        let mut f: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if f.len() < 10 {
            continue;
        }

        // Common -i/-e site filter.
        if let Some((exclude, expr)) = opts.site_filter {
            let mut pass = record_passes(&f, expr)?;
            if exclude {
                pass = !pass;
            }
            if !pass {
                continue;
            }
        }

        if opts.trim_alts {
            trim_alleles(&mut f);
        }

        // "Not a gVCF block": >1 ALT, or a single non-symbolic ALT.
        if !is_gvcf_block(&f[4]) {
            flush_block(&mut block, &groups, ngrp, Some(&f), &mut out);
            out.push_str(&f.join("\t"));
            out.push('\n');
            continue;
        }

        let pos0 = f[1].parse::<i64>().unwrap_or(1) - 1;
        let end = info_end(&f[7]).unwrap_or(pos0 + 1);
        let fmt: Vec<&str> = f[8].split(':').collect();
        let sample: Vec<&str> = f[9].split(':').collect();
        let (gq_key, gq) = read_gq(&fmt, &sample);
        let min_dp = read_min_dp(&fmt, &sample).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Expected one FORMAT/MIN_DP or FORMAT/DP value at {}:{}",
                    f[0],
                    pos0 + 1
                ),
            )
        })?;
        let pl = read_pl(&fmt, &sample);

        // First matching group; ngrp when nothing matches.
        let mut grp = ngrp;
        for (i, g) in groups.iter().enumerate() {
            let matched = match &g.expr {
                None => true,
                Some(e) => record_passes(&f, e)?,
            };
            if matched {
                grp = i;
                break;
            }
        }

        let rid = f[0].clone();
        if block.as_ref().is_some_and(|b| b.grp != grp) {
            flush_block(&mut block, &groups, ngrp, Some(&f), &mut out);
        }
        if let Some(b) = block.as_ref()
            && b.rid != rid
        {
            flush_block(&mut block, &groups, ngrp, None, &mut out);
        }

        if let Some(b) = block.as_mut() {
            if b.end < end {
                b.end = end;
            }
            if b.gq_key != GqKey::None && gq_key != GqKey::None && b.gq > gq {
                b.gq = gq;
            }
            if b.min_dp > min_dp {
                b.min_dp = min_dp;
            }
            for (bp, p) in b.pl.iter_mut().zip(pl.iter()) {
                if *bp > *p {
                    *bp = *p;
                }
            }
            continue;
        }

        block = Some(Block {
            fields: f.clone(),
            rid,
            end,
            pos0,
            min_dp,
            gq_key,
            gq: if gq_key != GqKey::None { gq } else { 0 },
            pl,
            grp,
        });
    }

    flush_block(&mut block, &groups, ngrp, None, &mut out);
    Ok(out)
}

/// Upstream `flush_block`: clamp end to the next record, write
/// `INFO/END`, `FORMAT/DP|GQ|RGQ|PL`, the group FILTER stamp, emit.
fn flush_block(
    block: &mut Option<Block>,
    groups: &[Group],
    ngrp: usize,
    next: Option<&[String]>,
    out: &mut String,
) {
    let Some(mut b) = block.take() else {
        return;
    };
    if let Some(n) = next {
        let next_pos0 = n[1].parse::<i64>().unwrap_or(1) - 1;
        // Upstream `gvcf->end - 1 >= rec->pos` (clamp to next record).
        if b.end > next_pos0 {
            b.end = next_pos0;
        }
    }

    let mut info = b.fields[7].clone();
    if b.pos0 + 1 < b.end {
        info = set_info_end(&info, b.end);
    }
    b.fields[7] = info;

    let fmt: Vec<String> = b.fields[8].split(':').map(str::to_owned).collect();
    let mut sample: Vec<String> = b.fields[9].split(':').map(str::to_owned).collect();
    sample.resize(fmt.len().max(sample.len()), ".".to_owned());
    let mut fmt = fmt;

    set_format(&mut fmt, &mut sample, "DP", &b.min_dp.to_string());
    match b.gq_key {
        GqKey::Gq => set_format(&mut fmt, &mut sample, "GQ", &b.gq.to_string()),
        GqKey::Rgq => set_format(&mut fmt, &mut sample, "RGQ", &b.gq.to_string()),
        GqKey::None => {}
    }
    if b.pl[0] >= 0 {
        let pl = format!("{},{},{}", b.pl[0], b.pl[1], b.pl[2]);
        set_format(&mut fmt, &mut sample, "PL", &pl);
    }
    b.fields[8] = fmt.join(":");
    b.fields[9] = sample.join(":");

    if b.grp < ngrp && !groups[b.grp].is_pass {
        let name = &groups[b.grp].name;
        b.fields[6] = if b.fields[6] == "." || b.fields[6] == "PASS" || b.fields[6].is_empty() {
            name.clone()
        } else {
            format!("{};{name}", b.fields[6])
        };
    }

    out.push_str(&b.fields.join("\t"));
    out.push('\n');
}

fn is_gvcf_block(alt: &str) -> bool {
    if alt == "." || alt.is_empty() {
        return true;
    }
    let alleles: Vec<&str> = alt.split(',').collect();
    if alleles.len() > 1 {
        return false;
    }
    alleles[0] == "<NON_REF>" || alleles[0] == "<*>"
}

/// Upstream `-a`: drop ALT alleles absent from the genotype, renumber
/// GT, then collapse a now ref-only record to a single-base REF.
fn trim_alleles(f: &mut [String]) {
    let alt = f[4].clone();
    if alt == "." || alt.is_empty() {
        return;
    }
    let alts: Vec<&str> = alt.split(',').collect();
    let fmt: Vec<&str> = f[8].split(':').collect();
    let Some(gt_idx) = fmt.iter().position(|k| *k == "GT") else {
        return;
    };

    let mut used = vec![false; alts.len() + 1];
    for sample in &f[9..] {
        if let Some(gt) = sample.split(':').nth(gt_idx) {
            for a in gt.split(['/', '|']) {
                if let Ok(idx) = a.parse::<usize>()
                    && idx < used.len()
                {
                    used[idx] = true;
                }
            }
        }
    }

    // old allele index -> new index (0 = REF, always kept).
    let mut remap = vec![None; alts.len() + 1];
    remap[0] = Some(0usize);
    let mut kept_alts = Vec::new();
    for (i, a) in alts.iter().enumerate() {
        if used[i + 1] {
            kept_alts.push((*a).to_owned());
            remap[i + 1] = Some(kept_alts.len());
        }
    }

    let samples: Vec<String> = f[9..]
        .iter()
        .map(|s| {
            let mut sub: Vec<String> = s.split(':').map(str::to_owned).collect();
            if gt_idx < sub.len() {
                sub[gt_idx] = remap_gt(&sub[gt_idx], &remap);
            }
            sub.join(":")
        })
        .collect();
    for (i, s) in samples.into_iter().enumerate() {
        f[9 + i] = s;
    }

    // Upstream (`gvcfz.c`): after `bcf_trim_alleles`, if the REF still has
    // a second base it is force-truncated to one base and *all* ALTs are
    // dropped (`bcf_update_alleles(...,1)`) — collapsing a multi-base-REF
    // record into a single-base ref block regardless of which ALTs the
    // genotype used. Only a single-base REF keeps its surviving ALTs.
    if f[3].len() > 1 {
        f[3] = f[3][..1].to_owned();
        f[4] = ".".to_owned();
    } else if kept_alts.is_empty() {
        f[4] = ".".to_owned();
    } else {
        f[4] = kept_alts.join(",");
    }
}

fn remap_gt(gt: &str, remap: &[Option<usize>]) -> String {
    let mut out = String::with_capacity(gt.len());
    let mut num = String::new();
    let flush = |num: &mut String, out: &mut String| {
        if num.is_empty() {
            return;
        }
        match num.parse::<usize>() {
            Ok(idx) => out.push_str(
                &remap
                    .get(idx)
                    .copied()
                    .flatten()
                    .map(|n| n.to_string())
                    .unwrap_or_else(|| ".".to_owned()),
            ),
            Err(_) => out.push_str(num),
        }
        num.clear();
    };
    for c in gt.chars() {
        if c == '/' || c == '|' {
            flush(&mut num, &mut out);
            out.push(c);
        } else if c == '.' {
            out.push('.');
        } else {
            num.push(c);
        }
    }
    flush(&mut num, &mut out);
    out
}

fn info_end(info: &str) -> Option<i64> {
    if info == "." {
        return None;
    }
    info.split(';')
        .find_map(|kv| kv.strip_prefix("END="))
        .and_then(|v| v.parse::<i64>().ok())
}

fn set_info_end(info: &str, end: i64) -> String {
    let kv = format!("END={end}");
    if info == "." || info.is_empty() {
        return kv;
    }
    let mut parts: Vec<String> = info
        .split(';')
        .filter(|p| !p.starts_with("END=") && *p != "END")
        .map(str::to_owned)
        .collect();
    // Upstream bcf_update_info_int32 places END first.
    parts.insert(0, kv);
    parts.join(";")
}

fn read_gq(fmt: &[&str], sample: &[&str]) -> (GqKey, i64) {
    for (key, kind) in [("GQ", GqKey::Gq), ("RGQ", GqKey::Rgq)] {
        if let Some(i) = fmt.iter().position(|k| *k == key)
            && let Some(v) = sample.get(i)
            && let Ok(n) = v.parse::<i64>()
        {
            return (kind, n);
        }
    }
    (GqKey::None, 0)
}

fn read_min_dp(fmt: &[&str], sample: &[&str]) -> Option<i64> {
    for key in ["MIN_DP", "DP"] {
        if let Some(i) = fmt.iter().position(|k| *k == key)
            && let Some(v) = sample.get(i)
            && let Ok(n) = v.parse::<i64>()
        {
            return Some(n);
        }
    }
    None
}

fn read_pl(fmt: &[&str], sample: &[&str]) -> [i64; 3] {
    if let Some(i) = fmt.iter().position(|k| *k == "PL")
        && let Some(v) = sample.get(i)
    {
        let vals: Vec<&str> = v.split(',').collect();
        if vals.len() == 3
            && let (Ok(a), Ok(b), Ok(c)) = (
                vals[0].parse::<i64>(),
                vals[1].parse::<i64>(),
                vals[2].parse::<i64>(),
            )
        {
            return [a, b, c];
        }
    }
    [-1, -1, -1]
}

fn set_format(fmt: &mut Vec<String>, sample: &mut Vec<String>, key: &str, value: &str) {
    if let Some(i) = fmt.iter().position(|k| k == key) {
        if i < sample.len() {
            sample[i] = value.to_owned();
        } else {
            sample.resize(i, ".".to_owned());
            sample.push(value.to_owned());
        }
    } else {
        fmt.push(key.to_owned());
        sample.push(value.to_owned());
    }
}

fn record_passes(fields: &[String], expr: &str) -> io::Result<bool> {
    let context = record_context(fields);
    Ok(
        bcffilter::eval_expression_with(expr, &context, |name, sample_index| {
            if sample_index.is_some() {
                return None;
            }
            crate::commands::filter::record_lookup(name, fields)
        })?
        .truthy(),
    )
}

fn record_context(fields: &[String]) -> EvalContext {
    if fields.len() <= 9 {
        return EvalContext::new();
    }
    let format_keys: Vec<&str> = fields[8].split(':').collect();
    fields[9..]
        .iter()
        .fold(EvalContext::new(), |context, sample| {
            let values: Vec<&str> = sample.split(':').collect();
            context.with_sample(
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
                        } else {
                            FilterValue::String(raw.to_owned())
                        };
                        ((*key).to_owned(), value)
                    })
                    .collect::<Vec<_>>(),
            )
        })
}

fn inject_headers(out: &mut String, groups: &[Group]) {
    out.push_str(
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">\n",
    );
    // Upstream uses the raw group-by string (with `"`→`'`) as every
    // non-PASS FILTER description.
    let desc = groups
        .iter()
        .map(|x| {
            format!(
                "{}:{}",
                x.name,
                x.expr.clone().unwrap_or_else(|| "-".to_owned())
            )
        })
        .collect::<Vec<_>>()
        .join("; ")
        .replace('"', "'");
    for g in groups {
        if !g.is_pass {
            out.push_str(&format!(
                "##FILTER=<ID={},Description=\"{desc}\">\n",
                g.name
            ));
        }
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
        ".bcftools-rs-gvcfz-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_groups() {
        let g = parse_groups("PASS:GQ>10; FLT:-").unwrap();
        assert_eq!(g.len(), 2);
        assert!(g[0].is_pass && g[0].expr.as_deref() == Some("GQ>10"));
        assert!(!g[1].is_pass && g[1].expr.is_none());
    }

    #[test]
    fn gvcf_block_classification() {
        assert!(is_gvcf_block("."));
        assert!(is_gvcf_block("<NON_REF>"));
        assert!(is_gvcf_block("<*>"));
        assert!(!is_gvcf_block("A"));
        assert!(!is_gvcf_block("A,C"));
    }

    #[test]
    fn trim_drops_unused_alt_and_collapses_ref() {
        // GT 0/0 -> all ALTs unused -> ref-only, REF truncated to 1 base.
        let mut f: Vec<String> = "chr1\t10\t.\tAAC\tA,C\t.\t.\t.\tGT:DP\t0/0:7"
            .split('\t')
            .map(str::to_owned)
            .collect();
        trim_alleles(&mut f);
        assert_eq!(f[3], "A");
        assert_eq!(f[4], ".");
    }

    #[test]
    fn trim_keeps_used_alt_and_renumbers_gt() {
        let mut f: Vec<String> = "chr1\t10\t.\tA\tC,G,T\t.\t.\t.\tGT:DP\t0/3:7"
            .split('\t')
            .map(str::to_owned)
            .collect();
        trim_alleles(&mut f);
        assert_eq!(f[4], "T");
        assert_eq!(f[9], "0/1:7");
    }

    #[test]
    fn set_info_end_places_end_first() {
        assert_eq!(set_info_end(".", 99), "END=99");
        assert_eq!(set_info_end("AN=2;DP=5", 99), "END=99;AN=2;DP=5");
        assert_eq!(set_info_end("END=1;AN=2", 99), "END=99;AN=2");
    }
}
