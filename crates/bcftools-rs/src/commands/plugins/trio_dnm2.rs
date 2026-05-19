//! `bcftools +trio-dnm2` (upstream `bcftools/plugins/trio-dnm2.c`).
//!
//! First slice: the `--use-NAIVE` GT-only de-novo-mutation model
//! (`-p [1X:|2X:]proband,father,mother --use-NAIVE`). Ports the
//! `seq1`/`seq2`/`seq3` genotype encoding, the autosomal/chrX/chrXX
//! Mendelian-transmission de-novo predicates (the `tprob==0` part of
//! upstream `init_tprob_mprob{,_chrX,_chrXX}`), `set_trio_GT`
//! (incl. the >4-allele remap), the default GRCh37 chrX regions, and
//! `process_record_naive` writing `FORMAT/DNM` (flag) + `FORMAT/VA`
//! (the de-novo allele). Validated by piping through our own
//! `bcftools query` (upstream test.pl rows 768-769).
//!
//! Deferred (TODO.md): the ACM (default) and `--use-DNG` likelihood
//! models, `--ppl`, `--force-AD`, `--with-pAD`, `--strictly-novel`,
//! `--dnm-tag` non-flag types, PED-file `-P`, VAF/VA from AD.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

// Upstream seq1/seq2: genotype-index (0..9) → its two allele indices.
const SEQ1: [usize; 10] = [0, 1, 1, 2, 2, 2, 3, 3, 3, 3];
const SEQ2: [usize; 10] = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3];
// Upstream seq3: (1<<ial)|(1<<jal) bitmask (1..12) → genotype index.
const SEQ3: [i8; 13] = [-1, 0, 2, 1, 5, 3, 4, -1, 9, 6, 7, -1, 8];

#[derive(Clone, Copy)]
enum Kind {
    Autosomal,
    ChrX,
    ChrXX,
}

/// `(is_denovo, denovo_allele)` for a genotype-index trio, mirroring
/// the `tprob==0` branch of upstream `init_tprob_mprob{,_chrX,_chrXX}`
/// (NAIVE uses the non-`--strictly-novel` `is_novel`).
fn denovo(kind: Kind, fi: usize, mi: usize, ci: usize) -> (bool, i32) {
    let (fa, fb) = (SEQ1[fi], SEQ2[fi]);
    let (ma, mb) = (SEQ1[mi], SEQ2[mi]);
    let (ca, cb) = (SEQ1[ci], SEQ2[ci]);
    match kind {
        Kind::Autosomal => {
            let allele = if ca != fa && ca != fb && ca != ma && ca != mb {
                ca
            } else {
                cb
            };
            let is_novel = !(((ca == fa || ca == fb) && (cb == ma || cb == mb))
                || ((ca == ma || ca == mb) && (cb == fa || cb == fb)));
            (is_novel, allele as i32)
        }
        Kind::ChrX => {
            let allele = if ca != ma && ca != mb { ca } else { cb };
            let denovo = if ca != cb {
                // male cannot be heterozygous in X (mosaic); tprob==0
                true
            } else if ca == ma || ca == mb {
                false // inherited
            } else {
                true // de novo
            };
            (denovo, allele as i32)
        }
        Kind::ChrXX => {
            let allele = if ca != fa && ca != fb && ca != ma && ca != mb {
                ca
            } else {
                cb
            };
            if fa != fb {
                // father cannot be het in X → fall back to autosomal.
                return denovo(Kind::Autosomal, fi, mi, ci);
            }
            let inherited =
                (ca == fa && (cb == ma || cb == mb)) || (cb == fa && (ca == ma || ca == mb));
            (!inherited, allele as i32)
        }
    }
}

/// `set_trio_GT` (+ `set_trio_GT_many_alts`): GT strings → the
/// `(1<<allele)` bitmask per member. `gts` order is
/// `[father, mother, child]`. `None` ⇒ skip the trio.
fn set_trio_gt(
    gt_f: &str,
    gt_m: &str,
    gt_c: &str,
    n_allele: usize,
    ignore_father: bool,
) -> Option<[usize; 3]> {
    let raw = [gt_f, gt_m, gt_c];
    let mut gts = [0usize; 3];
    let mut alt_idx: Vec<i32> = vec![-1; n_allele.max(1)];
    let mut nused = 0i32;
    for (j, g) in raw.iter().enumerate() {
        // j: 0=father, 1=mother, 2=child.
        for tok in g.split(['/', '|']) {
            if tok.is_empty() {
                continue;
            }
            if tok == "." {
                if j != 0 || !ignore_father {
                    return None;
                }
                continue; // father ignored (male chrX): missing allowed
            }
            let ial: usize = tok.parse().ok()?;
            let bit = if n_allele <= 4 {
                if ial > 3 {
                    return None;
                }
                ial
            } else {
                if ial >= alt_idx.len() {
                    return None;
                }
                if alt_idx[ial] == -1 {
                    alt_idx[ial] = nused;
                    nused += 1;
                    if nused > 4 {
                        return None;
                    }
                }
                alt_idx[ial] as usize
            };
            gts[j] |= 1 << bit;
            if gts[j] == 0 || gts[j] >= 13 {
                return None;
            }
        }
        if gts[j] == 0 && (j != 0 || !ignore_father) {
            return None;
        }
    }
    Some(gts)
}

fn seq3_of(bitmask: usize) -> Option<usize> {
    if bitmask == 0 || bitmask >= 13 {
        return None;
    }
    let v = SEQ3[bitmask];
    if v < 0 { None } else { Some(v as usize) }
}

/// Default GRCh37/GRCh38 chrX pseudo-/non-autosomal regions
/// (`X:`/`chrX:` 1-based, inclusive). The fixtures use the GRCh37
/// default.
fn chrx_ranges(build: &str) -> Vec<(i64, i64)> {
    match build {
        "GRCh38" => vec![(1, 9999), (2781480, 155701381)],
        _ => vec![(1, 60000), (2699521, 154931043)],
    }
}

fn is_chrx(chrom: &str, pos: i64, reflen: i64, ranges: &[(i64, i64)]) -> bool {
    if chrom != "X" && chrom != "chrX" {
        return false;
    }
    let (lo, hi) = (pos, pos + reflen - 1);
    ranges.iter().any(|&(a, b)| lo <= b && hi >= a)
}

pub struct Options<'a> {
    /// `-p`/`--pfm` value: `[1X:|2X:]proband,father,mother`.
    pub pfm: &'a str,
    /// `--chrX-list` build (`GRCh37`/`GRCh38`) or `None` ⇒ GRCh37.
    pub chrx_build: Option<&'a str>,
}

pub fn run(input: &Path, opts: Options<'_>) -> io::Result<String> {
    // Parse `-p [1X:|2X:]P,F,M`.
    let (pfm, is_male) = if let Some(r) = opts.pfm.strip_prefix("1X:") {
        (r, true)
    } else if let Some(r) = opts.pfm.strip_prefix("2X:") {
        (r, false)
    } else {
        (opts.pfm, false)
    };
    let names: Vec<&str> = pfm.split(',').collect();
    if names.len() != 3 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "trio-dnm2: -p expects proband,father,mother",
        ));
    }
    let (proband, father, mother) = (names[0], names[1], names[2]);
    let ranges = chrx_ranges(opts.chrx_build.unwrap_or("GRCh37"));

    let text = read_vcf_text(input)?;
    let mut out = String::with_capacity(text.len() + 1024);
    // Sample-column indices (0-based within sample columns).
    let mut ci_idx = usize::MAX;
    let mut fi_idx = usize::MAX;
    let mut mi_idx = usize::MAX;

    for line in text.lines() {
        if line.starts_with("##") {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if let Some(rest) = line.strip_prefix("#CHROM") {
            out.push_str(
                "##FORMAT=<ID=DNM,Number=1,Type=Integer,Description=\"De-novo mutation score given as 1 for Mendelian-incompatible genotypes\">\n",
            );
            out.push_str(
                "##FORMAT=<ID=VA,Number=1,Type=Integer,Description=\"The de-novo allele\">\n",
            );
            let samples: Vec<&str> = rest.split('\t').skip(9).collect();
            for (i, s) in samples.iter().enumerate() {
                if *s == proband {
                    ci_idx = i;
                } else if *s == father {
                    fi_idx = i;
                } else if *s == mother {
                    mi_idx = i;
                }
            }
            if ci_idx == usize::MAX || fi_idx == usize::MAX || mi_idx == usize::MAX {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "trio-dnm2: a sample from -p is not present",
                ));
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        out.push_str(&process_record(
            line, ci_idx, fi_idx, mi_idx, is_male, &ranges,
        ));
        out.push('\n');
    }
    Ok(out)
}

#[allow(clippy::too_many_arguments)]
fn process_record(
    line: &str,
    ci_idx: usize,
    fi_idx: usize,
    mi_idx: usize,
    is_male: bool,
    ranges: &[(i64, i64)],
) -> String {
    let mut f: Vec<&str> = line.split('\t').collect();
    if f.len() < 10 {
        return line.to_owned();
    }
    let gt_pos = f[8].split(':').position(|k| k == "GT");
    let Some(gt_pos) = gt_pos else {
        return line.to_owned(); // no GT → unchanged
    };
    let n_allele = if f[4] == "." {
        1
    } else {
        1 + f[4].split(',').count()
    };
    let pos: i64 = f[1].parse().unwrap_or(0);
    let reflen = f[3].len() as i64;
    let chrx = is_chrx(f[0], pos, reflen, ranges);

    let kind = if !chrx {
        Kind::Autosomal
    } else if is_male {
        Kind::ChrX
    } else {
        Kind::ChrXX
    };
    let ignore_father = chrx && is_male;

    let gt_of = |si: usize| -> &str {
        f.get(9 + si)
            .and_then(|c| c.split(':').nth(gt_pos))
            .unwrap_or(".")
    };

    let mut dnm: Option<(i32, i32)> = None; // (score, allele) for child
    if let Some(gts) = set_trio_gt(
        gt_of(fi_idx),
        gt_of(mi_idx),
        gt_of(ci_idx),
        n_allele,
        ignore_father,
    ) && let (Some(fi), Some(mi), Some(c)) = (seq3_of(gts[0]), seq3_of(gts[1]), seq3_of(gts[2]))
    {
        let (is_dnm, allele) = denovo(kind, fi, mi, c);
        if is_dnm {
            dnm = Some((1, allele));
        }
    }

    let Some((score, allele)) = dnm else {
        return line.to_owned(); // no DNM at this site → record unchanged
    };

    // Append FORMAT/DNM + FORMAT/VA: child gets the values, every other
    // sample is missing (upstream `bcf_int32_missing` → `.`).
    let nsmpl = f.len() - 9;
    let new_format = format!("{}:DNM:VA", f[8]);
    let new_samples: Vec<String> = (0..nsmpl)
        .map(|si| {
            let base = f[9 + si];
            if si == ci_idx {
                format!("{base}:{score}:{allele}")
            } else {
                format!("{base}:.:.")
            }
        })
        .collect();
    f[8] = &new_format;
    for (si, col) in f[9..].iter_mut().enumerate() {
        *col = new_samples[si].as_str();
    }
    f.join("\t")
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
        ".bcftools-rs-trio-dnm2-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn autosomal_denovo_basics() {
        // 0/1 child from 0/0 + 0/0 parents → de novo.
        let g00 = seq3_of(1).unwrap(); // {0}
        let g01 = seq3_of(0b11).unwrap(); // {0,1}
        assert!(denovo(Kind::Autosomal, g00, g00, g01).0);
        // 0/1 child from 0/0 + 0/1 → inherited, not de novo.
        assert!(!denovo(Kind::Autosomal, g00, g01, g01).0);
    }

    #[test]
    fn seq3_lookup() {
        assert_eq!(seq3_of(1), Some(0)); // 0/0
        assert_eq!(seq3_of(2), Some(2)); // 1/1
        assert_eq!(seq3_of(0b11), Some(1)); // 0/1
        assert_eq!(seq3_of(7), None); // invalid bitmask
        assert_eq!(seq3_of(0), None);
    }

    #[test]
    fn set_trio_gt_skips_missing() {
        assert!(set_trio_gt("0/0", "0/0", "./.", 2, false).is_none());
        // father missing tolerated when ignore_father.
        assert!(set_trio_gt("./.", "0/0", "0/1", 2, true).is_some());
    }
}
