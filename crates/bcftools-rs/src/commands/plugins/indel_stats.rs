//! `bcftools +indel-stats` (upstream `bcftools/plugins/indel-stats.c`).
//!
//! Indel summary stats for the default (no `-p`/no `-i`/`-e`) case: site
//! counts (SN), the FORMAT/AD variant-allele-frequency distribution (DVAF),
//! the indel-length distribution (DLEN), and the mean minor-allele fraction
//! at HET indel genotypes vs indel length (DFRAC/NFRAC). Records are
//! pre-filtered to those with at least one INDEL ALT, exactly like upstream's
//! `bcf_get_variant_types & VCF_INDEL` gate. The `-p` trio/de-novo mode and
//! `-i`/`-e` filter-threshold scanning need the not-yet-ported PED/filter
//! infrastructure and are tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::vcf_compat::normalize_vcf_text;

const MAX_LEN: i32 = 20;
const NVAF: usize = 20;

// Rendered verbatim from upstream `destroy()` for the default MAX_LEN=20,
// NVAF=20 (the test passes no --max-len/--nvaf). The harness strips `^CMD`.
const HEADER: &str = "\
# CMD line shows the command line used to generate this output
# DEF lines define expressions for all tested thresholds
# SN* summary number for every threshold:
#   1) SN*, filter id
#   2) number of samples (or trios with -p)
#   3) number of indel sites total
#   4) number of indel sites that pass the filter (and, with -p, have a de novo indel)
#   5) number of indel genotypes that pass the filter (and, with -p, are de novo)
#   6) number of insertions (site-wise, not genotype-wise)
#   7) number of deletions (site-wise, not genotype-wise)
#   8) number of frameshifts (site-wise, not genotype-wise)
#   9) number of inframe indels (site-wise, not genotype-wise)
#
# DVAF* lines report indel variant allele frequency (VAF) distribution for every threshold,
#   k-th bin corresponds to the frequency k/(nVAF-1):
#   1) DVAF*, filter id
#   2) nVAF, number of bins which split the [0,1] VAF interval.
#   3-22) counts of indel genotypes in the VAF bin. For non-reference hets, the VAF of the less supported allele is recorded
#
# DLEN* lines report indel length distribution for every threshold. When genotype fields are available,
#   the counts correspond to the number of genotypes, otherwise the number of sites are given.
#   The k-th bin corresponds to the indel size k-MAX_LEN, negative for deletions, positive for insertions.
#   The first/last bin contains also all deletions/insertions larger than MAX_LEN:
#   1) DLEN*, filter id
#   2) maximum indel length
#   3-43) counts of indel lengths (-max,..,0,..,max), all unique alleles in a genotype are recorded (alt hets increase the counters 2x, alt homs 1x)
#
# DFRAC* lines report the mean minor allele fraction at HET indel genotypes as a function of indel size.
#   The format is the same as for DLEN:
#   1) DFRAC*, filter id
#   2) maximum indel length
#   3-43) mean fraction at indel lengths (-max,..,0,..,max)
#
# NFRAC* lines report the number of indels informing the DFRAC distribution.
#   1) NFRAC*, filter id
#   2) maximum indel length
#   3-43) counts at indel lengths (-max,..,0,..,max)
#
";

fn len2bin(len: i32) -> i32 {
    if len < -MAX_LEN {
        0
    } else if len > MAX_LEN {
        2 * MAX_LEN
    } else {
        MAX_LEN + len
    }
}

fn vaf2bin(vaf: f32) -> usize {
    // C: `return vaf*(NVAF-1);` (float -> int truncation toward zero).
    (vaf * (NVAF as f32 - 1.0)) as i32 as usize
}

struct Stats {
    nvaf: Vec<u32>,
    nlen: Vec<u32>,
    nfrac: Vec<u32>,
    dfrac: Vec<f64>,
    npass_gt: u32,
    npass: u32,
    nsites: u32,
    nins: u32,
    ndel: u32,
    nframeshift: u32,
    ninframe: u32,
}

impl Stats {
    fn new() -> Self {
        let nbins = (2 * MAX_LEN + 1) as usize;
        Stats {
            nvaf: vec![0; NVAF],
            nlen: vec![0; nbins],
            nfrac: vec![0; nbins],
            dfrac: vec![0.0; nbins],
            npass_gt: 0,
            npass: 0,
            nsites: 0,
            nins: 0,
            ndel: 0,
            nframeshift: 0,
            ninframe: 0,
        }
    }
}

/// Reads the input VCF/BCF and returns the indel-stats report text (with the
/// `CMD` line the harness strips via `grep -v ^CMD`).
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(compute(&text))
}

enum Gt {
    Missing,
    Als(i32, i32),
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
        return Gt::Als(a0, a0); // hemizygous: als[1]=als[0]
    }
    let t1 = toks[1];
    if t1 == "." || t1.is_empty() {
        return Gt::Missing;
    }
    let Ok(a1) = t1.parse::<i32>() else {
        return Gt::Missing;
    };
    Gt::Als(a0, a1)
}

fn compute(text: &str) -> String {
    let mut samples = 0usize;
    let mut st = Stats::new();

    for line in text.lines() {
        if line.starts_with("#CHROM") {
            let cols = line.split('\t').count();
            samples = cols.saturating_sub(9);
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 8 {
            continue;
        }
        process_record(&f, samples, &mut st);
    }

    let mut out = String::new();
    out.push_str(HEADER);
    out.push_str("CMD\tbcftools +indel-stats\n");
    out.push_str("DEF\tFLT0\tall\n");

    out.push_str(&format!(
        "SN0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        samples, st.nsites, st.npass, st.npass_gt, st.nins, st.ndel, st.nframeshift, st.ninframe
    ));

    out.push_str(&format!("DVAF0\t{NVAF}"));
    for v in &st.nvaf {
        out.push_str(&format!("\t{v}"));
    }
    out.push('\n');

    out.push_str(&format!("DLEN0\t{MAX_LEN}"));
    for v in &st.nlen {
        out.push_str(&format!("\t{v}"));
    }
    out.push('\n');

    out.push_str(&format!("DFRAC0\t{MAX_LEN}"));
    for (j, &nf) in st.nfrac.iter().enumerate() {
        if nf != 0 {
            out.push_str(&format!("\t{:.2}", st.dfrac[j] / nf as f64));
        } else {
            out.push_str("\t.");
        }
    }
    out.push('\n');

    out.push_str(&format!("NFRAC0\t{MAX_LEN}"));
    for v in &st.nfrac {
        out.push_str(&format!("\t{v}"));
    }
    out.push('\n');

    out
}

/// Net indel length of ALT allele `ai` (>=1): `len(ALT) - len(REF)`.
fn var_n(reference: &str, alt: &str) -> i32 {
    alt.len() as i32 - reference.len() as i32
}

fn is_indel(reference: &str, alleles: &[&str], a: i32) -> bool {
    if a <= 0 || (a as usize) >= alleles.len() {
        return false; // allele 0 is REF; OOB treated as non-indel
    }
    classify_variant(reference, alleles[a as usize])
        .variant_type
        .contains(VariantType::INDEL)
}

fn process_record(f: &[&str], samples: usize, st: &mut Stats) {
    let reference = f[3];
    // alleles[0] = REF, alleles[i>=1] = ALT i-1
    let mut alleles: Vec<&str> = vec![reference];
    if f[4] != "." {
        alleles.extend(f[4].split(','));
    }
    let n_allele = alleles.len();

    // Record-level prefilter: at least one INDEL ALT (bcf_get_variant_types).
    let any_indel = (1..n_allele).any(|i| is_indel(reference, &alleles, i as i32));
    if !any_indel {
        return;
    }
    st.nsites += 1;

    let have_samples = samples > 0 && f.len() >= 10;
    let (gt_slot, ad_slot) = if have_samples {
        let fmt = f[8];
        (
            fmt.split(':').position(|k| k == "GT"),
            fmt.split(':').position(|k| k == "AD"),
        )
    } else {
        (None, None)
    };

    let star_allele: i32 = (1..n_allele)
        .find(|&i| alleles[i] == "*")
        .map(|i| i as i32)
        .unwrap_or(-1);
    let _ = star_allele; // only used by the (unported) trio path

    if have_samples && let Some(gt_slot) = gt_slot {
        let sample_cols = &f[9..];

        // Max ploidy across samples (htslib vector_end padding).
        let gts: Vec<&str> = sample_cols
            .iter()
            .map(|s| s.split(':').nth(gt_slot).unwrap_or("."))
            .collect();
        let max_ploidy = gts
            .iter()
            .map(|g| g.split(['/', '|']).count())
            .max()
            .unwrap_or(1);

        for (i, s) in sample_cols.iter().enumerate() {
            let Gt::Als(a0, a1) = parse_gt(gts[i], max_ploidy) else {
                continue; // missing genotype
            };
            if !is_indel(reference, &alleles, a0) && !is_indel(reference, &alleles, a1) {
                continue; // not an indel
            }
            let ad = ad_slot.and_then(|sl| s.split(':').nth(sl)).map(parse_ad);
            update_indel_stats(reference, &alleles, st, ad.as_deref(), a0, a1);
            st.npass_gt += 1;
        }
    }

    // CSQ-based frameshift/inframe (default tag CSQ).
    if let Some(csq) = info_value(f[7], "CSQ") {
        if csq.contains("inframe") {
            st.ninframe += 1;
        }
        if csq.contains("frameshift") {
            st.nframeshift += 1;
        }
    }

    for i in 1..n_allele {
        if !is_indel(reference, &alleles, i as i32) {
            continue;
        }
        let n = var_n(reference, alleles[i]);
        if n < 0 {
            st.ndel += 1;
        } else if n > 0 {
            st.nins += 1;
        }
        if !have_samples {
            let bin = len2bin(n);
            if bin >= 0 {
                st.nlen[bin as usize] += 1;
            }
        }
    }

    st.npass += 1;
}

/// Port of `update_indel_stats`. `ad` is the sample's FORMAT/AD (one entry
/// per allele; `None` = whole AD missing).
fn update_indel_stats(
    reference: &str,
    alleles: &[&str],
    st: &mut Stats,
    ad: Option<&[Option<i64>]>,
    als0: i32,
    als1: i32,
) {
    let Some(ad) = ad else {
        return; // no AD -> upstream would error; skip gracefully in this slice
    };

    let mut ntot: i64 = 0;
    for v in ad {
        match v {
            Some(x) => ntot += *x,
            None => continue,
        }
    }
    if ntot == 0 {
        return;
    }
    let adv = |a: i32| -> i64 { ad.get(a as usize).copied().flatten().unwrap_or(0) };

    let mut al0 = als0;
    let mut al1 = als1;
    if !is_indel(reference, alleles, al0) {
        // al0 not indel; the caller guarantees al1 is.
        al0 = als1;
        al1 = als0;
    } else if is_indel(reference, alleles, al1) && al0 != al1 {
        if adv(al0) < adv(al1) {
            std::mem::swap(&mut al0, &mut al1);
        }
        let bin = len2bin(var_n(reference, alleles[al1 as usize]));
        if bin >= 0 {
            st.nlen[bin as usize] += 1;
        }
    }

    let vaf = adv(al0) as f32 / ntot as f32;
    st.nvaf[vaf2bin(vaf)] += 1;

    let len_bin = len2bin(var_n(reference, alleles[al0 as usize]));
    if len_bin < 0 {
        return;
    }
    st.nlen[len_bin as usize] += 1;

    if al0 != al1 {
        let ntot2 = adv(al0) + adv(al1);
        if ntot2 != 0 {
            st.nfrac[len_bin as usize] += 1;
            st.dfrac[len_bin as usize] += adv(al0) as f64 / ntot2 as f64;
        }
    }
}

fn parse_ad(field: &str) -> Vec<Option<i64>> {
    field
        .split(',')
        .map(|t| {
            if t == "." {
                None
            } else {
                t.parse::<i64>().ok()
            }
        })
        .collect()
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    if info == "." {
        return None;
    }
    for kv in info.split(';') {
        let mut it = kv.splitn(2, '=');
        if it.next() == Some(key) {
            return it.next();
        }
    }
    None
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
        ".bcftools-rs-indel-stats-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bins() {
        assert_eq!(len2bin(0), 20);
        assert_eq!(len2bin(1), 21);
        assert_eq!(len2bin(-2), 18);
        assert_eq!(len2bin(-100), 0);
        assert_eq!(len2bin(100), 40);
        assert_eq!(vaf2bin(1.0), 19);
        assert_eq!(vaf2bin(0.5), 9);
        assert_eq!(vaf2bin(0.0), 0);
    }

    #[test]
    fn snp_only_record_is_skipped() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
20\t1\t.\tT\tA,C\t.\t.\t.\tGT\t1/2\n";
        let out = compute(vcf);
        assert!(out.contains("SN0\t1\t0\t0\t0\t0\t0\t0\t0\n"), "{out}");
    }

    #[test]
    fn het_insertion_records_vaf_len_frac() {
        // s1 0/1:10,10 with ALT = +1 insertion -> VAF bin 9, len bin 21,
        // frac 0.50; one insertion site.
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
20\t1\t.\tT\tTA\t.\t.\t.\tGT:AD\t0/1:10,10\n";
        let out = compute(vcf);
        assert!(out.contains("SN0\t1\t1\t1\t1\t1\t0\t0\t0\n"), "{out}");
        // DVAF bin 9 == 1
        let dvaf: Vec<&str> = out
            .lines()
            .find(|l| l.starts_with("DVAF0"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(dvaf[2 + 9], "1");
        let nfrac: Vec<&str> = out
            .lines()
            .find(|l| l.starts_with("NFRAC0"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(nfrac[2 + 21], "1");
        let dfrac: Vec<&str> = out
            .lines()
            .find(|l| l.starts_with("DFRAC0"))
            .unwrap()
            .split('\t')
            .collect();
        assert_eq!(dfrac[2 + 21], "0.50");
    }
}
