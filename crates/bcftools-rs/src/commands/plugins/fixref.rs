//! `bcftools +fixref` (upstream `bcftools/plugins/fixref.c`).
//!
//! Determines and fixes REF/ALT strand orientation against a FASTA
//! reference. Implements the self-contained conversion modes: `ref-alt`
//! and `swap` (REF/ALT column changes only) and `flip`/`flip-all` (also
//! flip + swap genotypes). Each record is annotated with INFO/FIXREF
//! (`none`/`swap`/`flip`/`flip,swap`/`skip`/`err`, plus `GT`). The `top`
//! (Illumina TOP-strand walking) and `id`/`-i` (dbSNP) modes need
//! infrastructure tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::reference::FastaReference;
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Mode {
    RefAlt,
    Swap,
    Flip,
    FlipAll,
    Stats,
}

impl Mode {
    pub fn parse(s: &str) -> Result<Mode, String> {
        match s.to_ascii_lowercase().as_str() {
            "ref-alt" => Ok(Mode::RefAlt),
            "swap" => Ok(Mode::Swap),
            "flip" => Ok(Mode::Flip),
            "flip-all" => Ok(Mode::FlipAll),
            "stats" => Ok(Mode::Stats),
            "top" => Err("fixref -m top is not supported in this slice".into()),
            "id" => Err("fixref -m id / -i is not supported in this slice".into()),
            other => Err(format!("fixref -m {other} not recognised")),
        }
    }
}

// FIXREF dirty bits, index order matching upstream `info_annots`.
const FIX_ERR: u32 = 1 << 0;
const FIX_SKIP: u32 = 1 << 1;
const FIX_NONE: u32 = 1 << 2;
const FIX_FLIP: u32 = 1 << 3;
const FIX_SWAP: u32 = 1 << 4;
const FIX_GT: u32 = 1 << 5;
const INFO_ANNOTS: [&str; 6] = ["err", "skip", "none", "flip", "swap", "GT"];

fn nt2int(c: u8) -> i32 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => -1,
    }
}
fn int2nt(i: i32) -> char {
    ['A', 'C', 'G', 'T'][i as usize]
}
/// Complement: A<->T, C<->G (`revint` in upstream: "3210"[x]).
fn revint(i: i32) -> i32 {
    3 - i
}

fn dirty_str(dirty: u32) -> String {
    let mut parts = Vec::new();
    for (i, name) in INFO_ANNOTS.iter().enumerate() {
        if dirty & (1 << i) != 0 {
            parts.push(*name);
        }
    }
    parts.join(",")
}

/// Reads inputs and returns the fixref-annotated VCF text.
pub fn run(input: &Path, fasta: &Path, mode: Mode, discard: bool, tag: &str) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let fa = FastaReference::open(fasta)?;
    Ok(process(&text, &fa, mode, discard, tag))
}

/// Flip 0<->1 in a GT string, preserving `/` and `|` separators.
fn swap_gt(gt: &str) -> String {
    let mut out = String::with_capacity(gt.len());
    for ch in gt.chars() {
        match ch {
            '0' => out.push('1'),
            '1' => out.push('0'),
            c => out.push(c),
        }
    }
    out
}

fn process(text: &str, fa: &FastaReference, mode: Mode, discard: bool, tag: &str) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = String::with_capacity(text.len() + 256);

    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
    let info_line = format!(
        "##INFO=<ID={tag},Number=.,Type=String,Description=\"The change made by bcftools/fixref\">"
    );
    for (idx, line) in lines.iter().enumerate() {
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("#CHROM") {
            out.push_str(&info_line);
            out.push('\n');
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

    let mut skip_chrom: Option<String> = None;
    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if f.len() < 8 {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if mode == Mode::Stats {
            continue; // upstream returns NULL for every record in stats mode
        }
        if skip_chrom.as_deref() == Some(f[0].as_str()) {
            continue; // whole sequence ignored (absent from the FASTA)
        }

        match process_record(&mut f, fa, mode, &mut skip_chrom) {
            Outcome::Write(dirty) => emit(&f, dirty, tag, &mut out),
            Outcome::SkipUnlessDiscard(dirty) => {
                if !discard {
                    emit(&f, dirty, tag, &mut out);
                }
            }
            Outcome::Suppress => {}
        }
    }
    out
}

enum Outcome {
    /// Always written, with this FIXREF dirty mask.
    Write(u32),
    /// `discard ? NULL : ret` — written (with FIXREF) only without `-d`.
    SkipUnlessDiscard(u32),
    /// Never written (sequence absent from the FASTA).
    Suppress,
}

fn emit(f: &[String], dirty: u32, tag: &str, out: &mut String) {
    let mut f = f.to_vec();
    if dirty != 0 {
        let val = dirty_str(dirty);
        let info = if f[7] == "." || f[7].is_empty() {
            format!("{tag}={val}")
        } else {
            format!("{};{tag}={val}", f[7])
        };
        f[7] = info;
    }
    out.push_str(&f.join("\t"));
    out.push('\n');
}

/// Classifies/repairs one record. Mutates `f` (REF/ALT/GT) in place.
fn process_record(
    f: &mut [String],
    fa: &FastaReference,
    mode: Mode,
    skip_chrom: &mut Option<String>,
) -> Outcome {
    let chrom = f[0].clone();
    let reference = f[3].clone();
    let alt = f[4].clone();

    // Skip non-SNPs (variant type must be exactly SNP across all ALTs).
    let vt = alt
        .split(',')
        .filter(|a| *a != ".")
        .fold(VariantType::REF, |acc, a| {
            acc | classify_variant(&reference, a).variant_type
        });
    if vt != VariantType::SNP {
        return Outcome::SkipUnlessDiscard(FIX_SKIP);
    }

    // Reference base at this position.
    if !fa.has_sequence(&chrom) {
        *skip_chrom = Some(chrom);
        return Outcome::Suppress;
    }
    let Ok(pos1) = f[1].parse::<i64>() else {
        return Outcome::SkipUnlessDiscard(FIX_SKIP);
    };
    let ir = match fa.fetch_region(&format!("{chrom}:{pos1}-{pos1}")) {
        Ok(b) if !b.is_empty() => nt2int(b[0]),
        _ => -1,
    };
    if ir < 0 {
        return Outcome::SkipUnlessDiscard(FIX_SKIP); // non-ACGT reference base
    }

    let n_allele = 1 + if alt == "." {
        0
    } else {
        alt.split(',').count()
    };
    if n_allele != 2 {
        return Outcome::SkipUnlessDiscard(FIX_SKIP); // non-biallelic
    }
    let ia = nt2int(reference.as_bytes()[0]);
    let ib = nt2int(alt.as_bytes()[0]);
    if ia < 0 || ib < 0 || ia == ib {
        return Outcome::SkipUnlessDiscard(FIX_SKIP);
    }

    let set = |f: &mut [String], r: i32, a: i32, gt: bool| {
        f[3] = int2nt(r).to_string();
        f[4] = int2nt(a).to_string();
        if gt && f.len() > 9 {
            // GT is the FORMAT field; swap genotypes 0<->1.
            if let Some(slot) = f[8].split(':').position(|k| k == "GT") {
                for s in f[9..].iter_mut() {
                    let mut parts: Vec<String> = s.split(':').map(|x| x.to_string()).collect();
                    if slot < parts.len() {
                        parts[slot] = swap_gt(&parts[slot]);
                    }
                    *s = parts.join(":");
                }
            }
        }
    };

    match mode {
        Mode::RefAlt => {
            if ir == ia {
                Outcome::Write(FIX_NONE)
            } else if ir == ib {
                set(f, ib, ia, false);
                Outcome::Write(FIX_SWAP)
            } else if ir == revint(ia) {
                set(f, revint(ia), revint(ib), false);
                Outcome::Write(FIX_FLIP)
            } else if ir == revint(ib) {
                set(f, revint(ib), revint(ia), false);
                Outcome::Write(FIX_FLIP | FIX_SWAP)
            } else {
                Outcome::Write(FIX_ERR)
            }
        }
        Mode::Swap => {
            if ir == ia {
                Outcome::Write(FIX_NONE)
            } else if ir == ib {
                set(f, ib, ia, false);
                Outcome::Write(FIX_SWAP)
            } else {
                Outcome::Write(FIX_ERR)
            }
        }
        Mode::Flip | Mode::FlipAll => {
            let pair = (1 << ia) | (1 << ib);
            if mode == Mode::Flip && (pair == 0x9 || pair == 0x6) {
                // ambiguous A/T or C/G — unresolved
                return Outcome::SkipUnlessDiscard(FIX_SKIP);
            }
            if ir == ia {
                Outcome::Write(FIX_NONE)
            } else if ir == ib {
                set(f, ib, ia, true);
                Outcome::Write(FIX_SWAP | FIX_GT)
            } else if ir == revint(ia) {
                set(f, revint(ia), revint(ib), false);
                Outcome::Write(FIX_FLIP)
            } else if ir == revint(ib) {
                set(f, revint(ib), revint(ia), true);
                Outcome::Write(FIX_FLIP | FIX_SWAP | FIX_GT)
            } else {
                Outcome::Write(FIX_ERR)
            }
        }
        Mode::Stats => Outcome::Write(0),
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
        ".bcftools-rs-fixref-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nt_helpers() {
        assert_eq!(nt2int(b'a'), 0);
        assert_eq!(nt2int(b'T'), 3);
        assert_eq!(nt2int(b'N'), -1);
        assert_eq!(int2nt(revint(0)), 'T'); // A<->T
        assert_eq!(int2nt(revint(1)), 'G'); // C<->G
    }

    #[test]
    fn dirty_formatting() {
        assert_eq!(dirty_str(FIX_NONE), "none");
        assert_eq!(dirty_str(FIX_SWAP), "swap");
        assert_eq!(dirty_str(FIX_SWAP | FIX_GT), "swap,GT");
        assert_eq!(dirty_str(FIX_FLIP | FIX_SWAP | FIX_GT), "flip,swap,GT");
        assert_eq!(dirty_str(FIX_SKIP), "skip");
    }

    #[test]
    fn swap_gt_keeps_phase() {
        assert_eq!(swap_gt("0/1"), "1/0");
        assert_eq!(swap_gt("1|0"), "0|1");
        assert_eq!(swap_gt("./."), "./.");
    }
}
