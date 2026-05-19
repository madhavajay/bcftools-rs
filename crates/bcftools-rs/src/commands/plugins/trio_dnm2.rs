//! `bcftools +trio-dnm2` (upstream `bcftools/plugins/trio-dnm2.c`).
//!
//! Implemented:
//! - `--use-NAIVE` GT-only model: `seq1`/`seq2`/`seq3` genotype
//!   encoding, autosomal/chrX/chrXX Mendelian-transmission de-novo
//!   predicates, `set_trio_GT` (incl. >4-allele remap), GRCh37 chrX
//!   regions, `FORMAT/DNM`(flag)+`VA` (test.pl 768-769).
//! - the default **ACM** likelihood model (autosomal, ≤4 alleles):
//!   `init_mf_priors`/`init_tprob_mprob`/`init_priors`, log-space
//!   helpers (`subtract_log`/`sum_log`/`phred2num`/`phred2log`/
//!   `log2phred`), `set_trio_PL` (normalised log-probs),
//!   `set_trio_QS_noisy` (SNV/INDEL pnoise), `process_trio_ACM`, the
//!   `DNM:log` transform + `FORMAT/DNM`(float)+`VA`+`VAF`-from-AD
//!   (test.pl 758/760/762/766 → `trio-dnm.{4.1,4.2,5.1,7.1}.out`).
//!
//! Deferred (TODO.md): `many_alts_trim` for >4 alleles
//! (`trio-dnm.8.*`), chrX ACM priors (`init_mf_priors_chrX/chrXX`),
//! `--use-DNG`, `--ppl`, `--force-AD`, `--with-pAD`,
//! `--strictly-novel`, `DNM:phred`/`prob`, PED-file `-P`. (Some
//! small-exponent `DNM:log` fixtures, e.g. `trio-dnm.6.2`, differ only
//! in our `query`'s float rendering vs C `%g`, not the model.)

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

// --- ACM (default) likelihood model -----------------------------------

fn phred2num(p: f64) -> f64 {
    10f64.powf(-0.1 * p)
}
fn phred2log(p: f64) -> f64 {
    -p / 4.3429
}
fn log2phred(n: f64) -> f64 {
    (4.3429 * n).abs()
}
/// `log(exp(a) - exp(b))`, upstream `subtract_log`.
fn subtract_log(a: f64, b: f64) -> f64 {
    a + (1.0 - (b - a).exp()).ln()
}
/// `log(exp(a) + exp(b))`, upstream `sum_log`.
fn sum_log(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    if a > b {
        (1.0 + (b - a).exp()).ln() + a
    } else {
        (1.0 + (a - b).exp()).ln() + b
    }
}

/// Upstream `count_unique_alleles` over the father/mother genotype
/// indices; `only_alts` skips the reference allele.
fn count_unique_alleles(fi: usize, mi: usize, only_alts: bool) -> usize {
    let mut als = [0u8; 4];
    for gi in [fi, mi] {
        als[SEQ1[gi]] = 1;
        als[SEQ2[gi]] = 1;
    }
    let beg = if only_alts { 1 } else { 0 };
    (beg..4).map(|i| als[i] as usize).sum()
}

/// Upstream `init_mf_priors` (autosomal parent-genotype prior). The
/// distinct `nref_mf` arms intentionally mirror upstream's separate
/// (commented) cases even where the value coincides.
#[allow(clippy::if_same_then_else)]
fn init_mf_priors(fi: usize, mi: usize) -> f64 {
    let (fa, fb) = (SEQ1[fi], SEQ2[fi]);
    let (ma, mb) = (SEQ1[mi], SEQ2[mi]);
    let nalt_mf = count_unique_alleles(fi, mi, true);
    let nref_mf = (fa == 0) as i32 + (fb == 0) as i32 + (ma == 0) as i32 + (mb == 0) as i32;
    let p_homref = 0.998;
    let p_poly = (1.0 - p_homref) * (1.0 - p_homref);
    let p_nonref = 1.0 - p_homref - p_poly;
    if nalt_mf >= 3 {
        1e-26
    } else if nalt_mf >= 2 {
        p_poly / 57.0
    } else if nref_mf == 4 {
        p_homref
    } else if nref_mf == 3 {
        p_nonref * (4.0 / 15.0) * (1.0 / 3.0)
    } else if nref_mf == 2 && ma == mb {
        p_nonref * (2.0 / 15.0) * (1.0 / 3.0)
    } else if nref_mf == 2 {
        p_nonref * (4.0 / 15.0) * (1.0 / 3.0)
    } else if nref_mf == 1 {
        p_nonref * (4.0 / 15.0) * (1.0 / 3.0)
    } else {
        p_nonref * (1.0 / 15.0) * (1.0 / 3.0)
    }
}

/// Upstream `init_tprob_mprob` (autosomal): `(tprob, mprob,
/// denovo_allele)`. NAIVE uses only `tprob==0`; ACM needs all three.
fn init_tprob_mprob(fi: usize, mi: usize, ci: usize, mrate: f64) -> (f64, f64, i32) {
    let (fa, fb) = (SEQ1[fi], SEQ2[fi]);
    let (ma, mb) = (SEQ1[mi], SEQ2[mi]);
    let (ca, cb) = (SEQ1[ci], SEQ2[ci]);
    let allele = if ca != fa && ca != fb && ca != ma && ca != mb {
        ca
    } else {
        cb
    } as i32;
    // Non-`--strictly-novel` is_novel (ACM default).
    let is_novel = !(((ca == fa || ca == fb) && (cb == ma || cb == mb))
        || ((ca == ma || ca == mb) && (cb == fa || cb == fb)));
    if !is_novel {
        let tprob = if fa == fb && ma == mb {
            1.0
        } else if fa == fb || ma == mb {
            0.5
        } else {
            0.25
        };
        (tprob, 1.0 - mrate, allele)
    } else {
        let mprob = if (ca == fa || ca == fb)
            || (ca == ma || ca == mb)
            || (cb == fa || cb == fb)
            || (cb == ma || cb == mb)
        {
            mrate
        } else {
            mrate * mrate
        };
        (0.0, mprob, allele)
    }
}

/// Autosomal priors tables (`init_priors`).
struct Priors {
    pprob: Vec<f64>,   // [fi*100 + mi*10 + ci]
    denovo: Vec<bool>, //  log(gt_prior*mprob*(tprob==0?1:tprob))
    dnv_allele: Vec<i32>,
}

fn init_priors_autosomal(mrate: f64) -> Priors {
    let mut pprob = vec![0.0f64; 1000];
    let mut denovo = vec![false; 1000];
    let mut dnv_allele = vec![0i32; 1000];
    for fi in 0..10 {
        for mi in 0..10 {
            let gt_prior = init_mf_priors(fi, mi);
            for ci in 0..10 {
                let (tprob, mprob, allele) = init_tprob_mprob(fi, mi, ci, mrate);
                let idx = fi * 100 + mi * 10 + ci;
                denovo[idx] = tprob == 0.0;
                dnv_allele[idx] = allele;
                pprob[idx] = (gt_prior * mprob * if tprob == 0.0 { 1.0 } else { tprob }).ln();
            }
        }
    }
    Priors {
        pprob,
        denovo,
        dnv_allele,
    }
}

/// Upstream `process_trio_ACM`: returns the DNM phred-ish score and
/// sets `(al0, al1)`. `pl`/`qs` are the per-member log-prob arrays
/// (`[father, mother, child]`).
#[allow(clippy::needless_range_loop)]
fn process_trio_acm(
    pr: &Priors,
    nals: usize,
    pl: &[Vec<f64>; 3],
    qs: &[Vec<f64>; 3],
) -> (f64, i32, i32) {
    let (mut al0, mut al1) = (0i32, 0i32);
    let mut sum = f64::NEG_INFINITY;
    let mut max = f64::NEG_INFINITY;
    let mut ci = 0usize;
    for ca in 0..nals {
        for cb in 0..=ca {
            let cals = (1usize << ca) | (1usize << cb);
            let cpl = pl[2][ci];
            let mut fi = 0usize;
            for fa in 0..nals {
                for fb in 0..=fa {
                    let fals = (1usize << fa) | (1usize << fb);
                    let mut fpl = 0.0;
                    for i in 0..nals {
                        if fals & (1 << i) != 0 {
                            fpl += subtract_log(0.0, qs[0][i]);
                        } else if cals & (1 << i) != 0 || fa == fb {
                            fpl += qs[0][i];
                        }
                    }
                    let mut mi = 0usize;
                    for ma in 0..nals {
                        for mb in 0..=ma {
                            let mals = (1usize << ma) | (1usize << mb);
                            let mut mpl = 0.0;
                            for i in 0..nals {
                                if mals & (1 << i) != 0 {
                                    mpl += subtract_log(0.0, qs[1][i]);
                                } else if cals & (1 << i) != 0 || ma == mb {
                                    mpl += qs[1][i];
                                }
                            }
                            let idx = fi * 100 + mi * 10 + ci;
                            let val = cpl + fpl + mpl + pr.pprob[idx];
                            sum = sum_log(sum, val);
                            if pr.denovo[idx] && max < val {
                                max = val;
                                if pr.dnv_allele[idx] == ca as i32 {
                                    al0 = cb as i32;
                                    al1 = ca as i32;
                                } else {
                                    al0 = ca as i32;
                                    al1 = cb as i32;
                                }
                            }
                            mi += 1;
                        }
                    }
                    fi += 1;
                }
            }
            ci += 1;
        }
    }
    (log2phred(subtract_log(0.0, max - sum)), al0, al1)
}

pub struct Options<'a> {
    /// `-p`/`--pfm` value: `[1X:|2X:]proband,father,mother`.
    pub pfm: &'a str,
    /// `--chrX-list` build (`GRCh37`/`GRCh38`) or `None` ⇒ GRCh37.
    pub chrx_build: Option<&'a str>,
    /// `true` ⇒ `--use-NAIVE`; `false` ⇒ the default ACM model.
    pub naive: bool,
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
    let has_ad = text.contains("##FORMAT=<ID=AD,");
    let priors = if opts.naive {
        None
    } else {
        Some(init_priors_autosomal(1e-8))
    };
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
            if opts.naive {
                out.push_str(
                    "##FORMAT=<ID=DNM,Number=1,Type=Integer,Description=\"De-novo mutation score given as 1 for Mendelian-incompatible genotypes\">\n",
                );
            } else {
                out.push_str(
                    "##FORMAT=<ID=DNM,Number=1,Type=Float,Description=\"De-novo mutation score given as log scaled value (bigger value = bigger confidence)\">\n",
                );
            }
            out.push_str(
                "##FORMAT=<ID=VA,Number=1,Type=Integer,Description=\"The de-novo allele\">\n",
            );
            if !opts.naive && has_ad {
                out.push_str(
                    "##FORMAT=<ID=VAF,Number=1,Type=Integer,Description=\"The percentage of ALT reads\">\n",
                );
            }
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
        let rec = match &priors {
            None => process_record(line, ci_idx, fi_idx, mi_idx, is_male, &ranges),
            Some(pr) => process_record_acm(line, ci_idx, fi_idx, mi_idx, pr),
        };
        out.push_str(&rec);
        out.push('\n');
    }
    Ok(out)
}

/// bcftools float→text (`%g`/6 over the f32-stored value); `0.0`→`0`.
fn fmt_float(x: f64) -> String {
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
    let exp = x.abs().log10().floor() as i32;
    if !(-4..6).contains(&exp) {
        let s = format!("{:.*e}", 5usize, x);
        let (m, e) = s.split_once('e').unwrap();
        let m = if m.contains('.') {
            m.trim_end_matches('0').trim_end_matches('.')
        } else {
            m
        };
        let ev: i32 = e.parse().unwrap_or(0);
        return format!("{m}e{}{:02}", if ev < 0 { '-' } else { '+' }, ev.abs());
    }
    let dec = (5 - exp).max(0) as usize;
    let s = format!("{x:.dec$}");
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').to_owned()
    } else {
        s
    }
}

/// FORMAT/AD-derived integer values for `name`-less per-allele arrays.
fn parse_int_field(cols: &[&str], fmt_keys: &[&str], si: usize, key: &str) -> Option<Vec<i64>> {
    let k = fmt_keys.iter().position(|x| *x == key)?;
    let v = cols.get(9 + si)?.split(':').nth(k)?;
    if v == "." {
        return None;
    }
    v.split(',').map(|x| x.parse::<i64>().ok()).collect()
}

/// ACM (default) / log-DNM model — upstream `process_record`
/// (non-naive) for the autosomal case.
fn process_record_acm(
    line: &str,
    ci_idx: usize,
    fi_idx: usize,
    mi_idx: usize,
    pr: &Priors,
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
    // Skip mono-allelic / reference-only sites (upstream `skip_site`).
    if n_allele == 1 || n_allele > 4 {
        return line.to_owned();
    }
    let fmt_keys: Vec<&str> = f[8].split(':').collect();
    let npl1 = n_allele * (n_allele + 1) / 2;

    // FORMAT/PL → normalized log-probs (set_trio_PL); members ordered
    // [father, mother, child].
    let members = [fi_idx, mi_idx, ci_idx];
    let mut pl: [Vec<f64>; 3] = [vec![], vec![], vec![]];
    let mut qs: [Vec<f64>; 3] = [vec![], vec![], vec![]];
    let mut ad: [Vec<i64>; 3] = [vec![], vec![], vec![]];
    for (j, &si) in members.iter().enumerate() {
        let Some(pl_i) = parse_int_field(&f, &fmt_keys, si, "PL") else {
            return line.to_owned();
        };
        if pl_i.len() != npl1 {
            return line.to_owned();
        }
        let mut v: Vec<f64> = pl_i.iter().map(|&p| phred2num(p as f64)).collect();
        let sum: f64 = v.iter().sum();
        for x in &mut v {
            *x = (*x / sum).ln();
        }
        pl[j] = v;
        ad[j] = parse_int_field(&f, &fmt_keys, si, "AD").unwrap_or_default();
        qs[j] = parse_int_field(&f, &fmt_keys, si, "QS")
            .map(|q| q.iter().map(|&x| x as f64).collect())
            .unwrap_or_default();
    }
    if qs.iter().any(|q| q.len() != n_allele) {
        return line.to_owned(); // ACM requires FORMAT/QS
    }
    let has_ad = ad.iter().all(|a| a.len() == n_allele);

    // set_trio_QS_noisy (autosomal): SNV pnoise frac=0.005/frac1=0.045,
    // abs=0; INDEL pnoise all-zero. n_ad kept (frac1≠0).
    let is_indel = {
        let r = f[3];
        f[4].split(',')
            .any(|a| !a.starts_with('<') && a != "*" && a != "." && a.len() != r.len())
    };
    let (pn_frac, pn_frac1) = if is_indel { (0.0, 0.0) } else { (0.005, 0.045) };
    let mut pqs: [Vec<f64>; 3] = [vec![], vec![], vec![]];
    for j in 0..3 {
        let (mut pn, mut pns) = (0.0, 0.0);
        if (pn_frac != 0.0 || pn_frac1 != 0.0) && j != 2 {
            let sum_qs: f64 = qs[j].iter().sum();
            pn = sum_qs * pn_frac;
            pns = sum_qs * pn_frac1;
        }
        pqs[j] = (0..n_allele)
            .map(|k| {
                let val = if has_ad && (ad[0][k] == 0 || ad[1][k] == 0) {
                    qs[j][k] - pns
                } else {
                    qs[j][k] - pn
                };
                phred2log(val.max(0.0))
            })
            .collect();
    }

    let (score, _al0, al1) = process_trio_acm(pr, n_allele, &pl, &pqs);
    // DNM:log output transform.
    let dnm = if score == f64::INFINITY {
        0.0
    } else {
        subtract_log(0.0, phred2log(score))
    };

    let nsmpl = f.len() - 9;
    // VAF: round(AD[al1]*100 / sum(AD)) per member, when al1<n_ad.
    let vaf_set = has_ad && (al1 as usize) < n_allele;
    let new_format = if vaf_set {
        format!("{}:DNM:VA:VAF", f[8])
    } else {
        format!("{}:DNM:VA", f[8])
    };
    let new_samples: Vec<String> = (0..nsmpl)
        .map(|si| {
            let base = f[9 + si];
            if si == ci_idx {
                let mut s = format!("{base}:{}:{al1}", fmt_float(dnm));
                if vaf_set {
                    let m = members.iter().position(|&x| x == si);
                    let adv = m.map(|mm| &ad[mm]);
                    let vaf = adv
                        .map(|a| {
                            let tot: i64 = a.iter().take(n_allele).sum();
                            if tot != 0 {
                                ((a[al1 as usize] * 100) as f64 / tot as f64).round() as i64
                            } else {
                                0
                            }
                        })
                        .unwrap_or(0);
                    s.push(':');
                    s.push_str(&vaf.to_string());
                }
                s
            } else if vaf_set {
                let m = members.iter().position(|&x| x == si);
                if let Some(mm) = m {
                    let a = &ad[mm];
                    let tot: i64 = a.iter().take(n_allele).sum();
                    let vaf = if tot != 0 {
                        ((a[al1 as usize] * 100) as f64 / tot as f64).round() as i64
                    } else {
                        0
                    };
                    format!("{base}:.:.:{vaf}")
                } else {
                    format!("{base}:.:.:.")
                }
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
