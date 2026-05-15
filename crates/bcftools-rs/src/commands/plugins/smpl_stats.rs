//! `bcftools +smpl-stats` (upstream `bcftools/plugins/smpl-stats.c`).
//!
//! Per-sample and per-site genotype statistics for the default (no `-i`/`-e`)
//! "all" filter: counts of passing/non-ref/homRR/homAA/het/hemi genotypes,
//! SNVs, indels, singletons, missing, and ts/tv, with a per-site rollup.
//! Allele counts come from `bcf_calc_ac` (INFO/AC+AN when present, otherwise
//! tallied from FORMAT/GT across every sample); ts/tv classification follows
//! the upstream per-base `bcf_acgt2int` walk. The `-i`/`-e` filter-threshold
//! scanning needs the not-yet-ported bcftools filter engine and is tracked
//! in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::vcf_compat::normalize_vcf_text;

const HEADER: &str = "\
# CMD line shows the command line used to generate this output
# DEF lines define expressions for all tested thresholds
# FLT* lines report numbers for every threshold and every sample:
#   1) filter id
#   2) sample
#   3) number of genotypes which pass the filter
#   4) number of non-reference genotypes
#   5) number of homozygous ref genotypes (0/0 or 0)
#   6) number of homozygous alt genotypes (1/1, 2/2, etc)
#   7) number of heterozygous genotypes (0/1, 1/2, etc)
#   8) number of hemizygous genotypes (0, 1, etc)
#   9) number of SNVs
#   10) number of indels
#   11) number of singletons
#   12) number of missing genotypes (./., ., ./0, etc)
#   13) number of transitions (alt het genotypes such as \"1/2\" are counted twice)
#   14) number of transversions (alt het genotypes such as \"1/2\" are counted twice)
#   15) overall ts/tv
# SITE* lines report numbers for every threshold:
#   1) filter id
#   2) number of sites which pass the filter
#   3) number of SNVs
#   4) number of indels
#   5) number of singletons
#   6) number of transitions (counted at most once at multiallelic sites)
#   7) number of transversions (counted at most once at multiallelic sites)
#   8) overall ts/tv
";

#[derive(Clone, Default)]
struct Stats {
    npass: u32,
    nnon_ref: u32,
    nhom_rr: u32,
    nhom_aa: u32,
    nhemi: u32,
    nhet: u32,
    n_snv: u32,
    n_indel: u32,
    nmissing: u32,
    nsingleton: u32,
    nts: u32,
    ntv: u32,
}

fn acgt2int(c: u8) -> i32 {
    match c.to_ascii_uppercase() {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => -1,
    }
}

fn tstv(nts: u32, ntv: u32) -> String {
    if ntv != 0 {
        format!("{:.2}", nts as f32 / ntv as f32)
    } else {
        // C prints `%.2f` of INFINITY.
        "inf".to_owned()
    }
}

/// Reads the input VCF/BCF and returns the smpl-stats report text (including
/// the `CMD` line the harness strips with `grep -v ^CMD`).
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(compute(&text))
}

/// Genotype parse mirroring upstream `parse_genotype`: returns `Err(())` for
/// a missing genotype, `Ok((a, true))` for hemizygous (one allele /
/// vector_end), `Ok((a0,a1,false))`-style otherwise.
enum Gt {
    Missing,
    Hemi(i32),
    Diploid(i32, i32),
}

fn parse_gt(gt: &str, max_ploidy: usize) -> Gt {
    let toks: Vec<&str> = gt.split(['/', '|']).collect();
    let t0 = toks.first().copied().unwrap_or(".");
    if t0 == "." || t0.is_empty() {
        return Gt::Missing;
    }
    let Ok(a0) = t0.parse::<i32>() else {
        return Gt::Missing;
    };
    if max_ploidy == 1 || toks.len() < 2 {
        return Gt::Hemi(a0);
    }
    let t1 = toks[1];
    if t1 == "." || t1.is_empty() {
        return Gt::Missing;
    }
    let Ok(a1) = t1.parse::<i32>() else {
        return Gt::Missing;
    };
    Gt::Diploid(a0, a1)
}

fn compute(text: &str) -> String {
    let mut samples: Vec<String> = Vec::new();
    let mut per_sample: Vec<Stats> = Vec::new();
    let mut site = Stats::default();

    for line in text.lines() {
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 {
                samples = cols[9..].iter().map(|s| s.to_string()).collect();
            }
            per_sample = vec![Stats::default(); samples.len()];
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if samples.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        process_record(&f, &samples, &mut per_sample, &mut site);
    }

    let mut out = String::new();
    out.push_str(HEADER);
    out.push_str("CMD\tbcftools +smpl-stats\n");
    out.push_str("DEF\tFLT0\tall\n");
    for (i, s) in samples.iter().enumerate() {
        let st = &per_sample[i];
        out.push_str(&format!(
            "FLT0\t{s}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            st.npass,
            st.nnon_ref,
            st.nhom_rr,
            st.nhom_aa,
            st.nhet,
            st.nhemi,
            st.n_snv,
            st.n_indel,
            st.nsingleton,
            st.nmissing,
            st.nts,
            st.ntv,
            tstv(st.nts, st.ntv),
        ));
    }
    out.push_str(&format!(
        "SITE0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        site.npass,
        site.n_snv,
        site.n_indel,
        site.nsingleton,
        site.nts,
        site.ntv,
        tstv(site.nts, site.ntv),
    ));
    out
}

fn process_record(f: &[&str], samples: &[String], per_sample: &mut [Stats], site: &mut Stats) {
    let reference = f[3];
    let alts: Vec<&str> = if f[4] == "." {
        Vec::new()
    } else {
        f[4].split(',').collect()
    };
    let n_allele = 1 + alts.len();

    // allele[i]: 0 = REF, i>=1 = ALT i-1
    let allele = |i: usize| -> &str { if i == 0 { reference } else { alts[i - 1] } };

    let fmt = f[8];
    let gt_slot = fmt.split(':').position(|k| k == "GT");
    let sample_cols = &f[9..];

    // Per-sample GT strings and max ploidy (htslib vector_end padding).
    let mut gts: Vec<&str> = Vec::with_capacity(sample_cols.len());
    let mut max_ploidy = 0usize;
    for s in sample_cols {
        let gt = match gt_slot {
            Some(idx) => s.split(':').nth(idx).unwrap_or("."),
            None => ".",
        };
        let p = gt.split(['/', '|']).count();
        if p > max_ploidy {
            max_ploidy = p;
        }
        gts.push(gt);
    }

    // bcf_calc_ac: INFO/AC+AN if present, else tally FORMAT/GT.
    let ac = match calc_ac_from_info(f[7], n_allele) {
        Some(ac) => ac,
        None => {
            let mut ac = vec![0i64; n_allele];
            let mut any = false;
            for gt in &gts {
                for tok in gt.split(['/', '|']) {
                    if tok == "." || tok.is_empty() {
                        continue;
                    }
                    if let Ok(a) = tok.parse::<usize>()
                        && a < n_allele
                    {
                        ac[a] += 1;
                        any = true;
                    }
                }
            }
            if !any {
                return; // bcf_calc_ac returned 0
            }
            ac
        }
    };

    let ref_code = if reference.len() == 1 {
        acgt2int(reference.as_bytes()[0])
    } else {
        -1
    };
    let star_allele: i32 = (1..n_allele)
        .find(|&i| allele(i) == "*")
        .map(|i| i as i32)
        .unwrap_or(-1);

    let mut site_pass = false;
    let mut site_snv = false;
    let mut site_indel = false;
    let mut site_has_ts = false;
    let mut site_has_tv = false;
    let mut site_singleton = false;

    for (i, gt) in gts.iter().enumerate() {
        let st = &mut per_sample[i];
        let _ = &samples[i];

        let (a0, a1, hemi) = match parse_gt(gt, max_ploidy) {
            Gt::Missing => {
                st.nmissing += 1;
                continue;
            }
            Gt::Hemi(a) => (a, a, true),
            Gt::Diploid(a0, a1) => (a0, a1, false),
        };

        if hemi {
            st.nhemi += 1;
        } else if a0 != a1 {
            st.nhet += 1;
        } else if a0 == 0 {
            st.nhom_rr += 1;
        } else {
            st.nhom_aa += 1;
        }

        st.npass += 1;
        site_pass = true;

        let als = [a0, a1];
        let mut has_nonref = false;
        for &a in &als {
            if a == star_allele || a == 0 {
                continue;
            }
            has_nonref = true;
        }
        if !has_nonref {
            continue;
        }
        st.nnon_ref += 1;

        let mut has_ts = false;
        let mut has_tv = false;
        let mut has_snv = false;
        let mut has_indel = false;
        for &a in &als {
            if a == 0 || a == star_allele {
                continue;
            }
            let ai = a as usize;
            if ai >= n_allele {
                continue;
            }
            if ac.get(ai).copied().unwrap_or(0) == 1 {
                st.nsingleton += 1;
                site_singleton = true;
            }
            let vt = classify_variant(reference, allele(ai)).variant_type;
            let is_snv_mnp = vt.contains(VariantType::SNP) || vt.contains(VariantType::MNP);
            if is_snv_mnp {
                let r = reference.as_bytes();
                let alt = allele(ai).as_bytes();
                let mut k = 0;
                while k < r.len() && k < alt.len() {
                    if r[k] == alt[k] {
                        k += 1;
                        continue;
                    }
                    let alt_code = acgt2int(alt[k]);
                    if (ref_code - alt_code).abs() == 2 {
                        has_ts = true;
                    } else {
                        has_tv = true;
                    }
                    has_snv = true;
                    k += 1;
                }
            } else if vt.contains(VariantType::INDEL) {
                has_indel = true;
            }
        }
        if has_ts {
            st.nts += 1;
            site_has_ts = true;
        }
        if has_tv {
            st.ntv += 1;
            site_has_tv = true;
        }
        if has_snv {
            st.n_snv += 1;
            site_snv = true;
        }
        if has_indel {
            st.n_indel += 1;
            site_indel = true;
        }
    }

    site.npass += site_pass as u32;
    site.n_snv += site_snv as u32;
    site.n_indel += site_indel as u32;
    site.nts += site_has_ts as u32;
    site.ntv += site_has_tv as u32;
    site.nsingleton += site_singleton as u32;
}

/// `bcf_calc_ac` INFO path: needs both INFO/AN and INFO/AC. `ac[0] = AN -
/// sum(AC)`, `ac[i] = AC[i-1]`. Returns `None` to fall back to GT tallying.
fn calc_ac_from_info(info: &str, n_allele: usize) -> Option<Vec<i64>> {
    if info == "." {
        return None;
    }
    let mut an: Option<i64> = None;
    let mut acv: Option<Vec<i64>> = None;
    for kv in info.split(';') {
        let mut it = kv.splitn(2, '=');
        match (it.next(), it.next()) {
            (Some("AN"), Some(v)) => an = v.parse().ok(),
            (Some("AC"), Some(v)) => acv = v.split(',').map(|x| x.parse::<i64>().ok()).collect(),
            _ => {}
        }
    }
    let (an, acv) = (an?, acv?);
    if acv.len() != n_allele - 1 {
        return None;
    }
    let mut ac = vec![0i64; n_allele];
    let sum: i64 = acv.iter().sum();
    ac[0] = an - sum;
    for (i, v) in acv.into_iter().enumerate() {
        ac[i + 1] = v;
    }
    Some(ac)
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
        ".bcftools-rs-smpl-stats-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_transitions() {
        // A<->G and C<->T are transitions (abs diff 2).
        assert_eq!((acgt2int(b'A') - acgt2int(b'G')).abs(), 2);
        assert_eq!((acgt2int(b'C') - acgt2int(b'T')).abs(), 2);
        assert_ne!((acgt2int(b'A') - acgt2int(b'C')).abs(), 2);
    }

    #[test]
    fn tstv_formatting() {
        assert_eq!(tstv(0, 1), "0.00");
        assert_eq!(tstv(2, 1), "2.00");
        assert_eq!(tstv(3, 0), "inf");
    }

    #[test]
    fn parse_gt_modes() {
        assert!(matches!(parse_gt("./.", 2), Gt::Missing));
        assert!(matches!(parse_gt(".", 1), Gt::Missing));
        assert!(matches!(parse_gt("1", 2), Gt::Hemi(1)));
        assert!(matches!(parse_gt("0/1", 2), Gt::Diploid(0, 1)));
        assert!(matches!(parse_gt("2/1", 2), Gt::Diploid(2, 1)));
    }

    #[test]
    fn small_singleton_and_indel() {
        // One sample carries the only ALT copy -> singleton; ALT is an
        // insertion -> indel.
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
20\t1\t.\tA\tATT\t.\t.\t.\tGT\t0/1\t0/0\n";
        let out = compute(vcf);
        assert!(
            out.contains("FLT0\ts1\t1\t1\t0\t0\t1\t0\t0\t1\t1\t0\t0\t0\tinf\n"),
            "{out}"
        );
        assert!(
            out.contains("FLT0\ts2\t1\t0\t1\t0\t0\t0\t0\t0\t0\t0\t0\t0\tinf\n"),
            "{out}"
        );
        assert!(out.contains("SITE0\t1\t0\t1\t1\t0\t0\tinf\n"), "{out}");
    }
}
