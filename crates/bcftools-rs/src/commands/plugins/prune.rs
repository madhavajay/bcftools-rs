//! `bcftools +prune` (upstream `bcftools/plugins/prune.c` + the windowed
//! `_prune_sites` path of `bcftools/vcfbuf.c`), `-n`/`-N` window subset.
//!
//! Keeps at most `-n N` sites per `-w` window, choosing which to drop by
//! `-N` mode: `1st` (keep the earliest), `maxAF` (drop the lowest allele
//! frequency, `--AF-tag` or computed). Faithful port of the `vcfbuf`
//! window flush condition and `_prune_sites` removal order. The LD/`-a`/`-m`
//! annotation modes, `-N rand` (needs `hts_drand48` parity), and `-i`/`-e`
//! filtering need infrastructure tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

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
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(process(&text, win, max_sites, mode, af_tag))
}

fn process(text: &str, win: i64, max_sites: i32, mode: Mode, af_tag: Option<&str>) -> String {
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
    fn first_mode_keeps_earliest() {
        let out = process(VCF, -2, 1, Mode::First, None);
        assert_eq!(positions(&out), vec!["101", "103", "105", "107"]);
    }

    #[test]
    fn maxaf_mode_drops_lowest_af() {
        let out = process(VCF, -2, 1, Mode::MaxAf, Some("AF"));
        assert_eq!(positions(&out), vec!["101", "104", "105", "107"]);
    }
}
