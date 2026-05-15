//! `bcftools +GTisec` (upstream `bcftools/plugins/GTisec.c`).
//!
//! Counts genotype intersections across all possible sample subsets,
//! emitting subset counts in banker's-sequence order (the upstream
//! `compute_bankers` / `choose` ordering is ported verbatim). Options:
//! `-m` (missing-genotype counts), `-v` (annotate with sample lists),
//! `-H` (human-readable, implies `-v`).

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const MISSING: u8 = 1 << 0;
const VERBOSE: u8 = 1 << 1;
const SMPORDER: u8 = 1 << 2;

/// `bcf_gt_allele(bcf_int32_vector_end)` = `(INT32_MIN >> 1) - 1`,
/// used as the second allele of a haploid/short genotype so its hash
/// key never collides with a diploid one.
const VEND_ALLELE: i64 = ((i32::MIN as i64) >> 1) - 1;

/// Parses `-m`/`-v`/`-H` (also accepts the bundled `-Hm`, `-mv`, ...).
pub fn parse_flags(spec: &str) -> u8 {
    let mut flag = 0;
    for ch in spec.chars() {
        match ch {
            'm' => flag |= MISSING,
            'v' => flag |= VERBOSE,
            'H' => flag |= SMPORDER | VERBOSE,
            _ => {}
        }
    }
    flag
}

struct Bankers {
    nsmp: usize,
    nsmpp2: u64,
    bankers: Vec<u32>,
    quick: HashMap<(u32, u32), u64>,
}

impl Bankers {
    fn new(nsmp: usize) -> Bankers {
        let nsmpp2 = 1u64 << nsmp;
        let mut b = Bankers {
            nsmp,
            nsmpp2,
            bankers: vec![0u32; nsmpp2 as usize],
            quick: HashMap::new(),
        };
        for j in 0..nsmpp2 {
            b.bankers[j as usize] = b.compute(j);
        }
        b
    }

    /// `choose(n, k)` — binomial coefficient (upstream `choose`).
    fn choose(&mut self, n: u32, k: u32) -> u64 {
        if n == 0 {
            return 0;
        }
        if n == k || k == 0 {
            return 1;
        }
        let k = if k > n / 2 { n - k } else { k };
        if let Some(v) = self.quick.get(&(n, k)) {
            return *v;
        }
        let v = self.choose(n - 1, k - 1) + self.choose(n - 1, k);
        self.quick.insert((n, k), v);
        v
    }

    /// `compute_bankers(a)` — the banker's number at position `a`.
    fn compute(&mut self, a: u64) -> u32 {
        if a == 0 {
            return 0;
        }
        if self.bankers[a as usize] != 0 {
            return self.bankers[a as usize];
        }
        if a >= self.nsmpp2 / 2 {
            let v = self.compute(self.nsmpp2 - (a + 1)) ^ ((self.nsmpp2 - 1) as u32);
            self.bankers[a as usize] = v;
            return v;
        }
        let mut c: u32 = 0;
        let mut n: u32 = self.nsmp as u32;
        let mut e: u64 = a;
        let mut binom = self.choose(n, c);
        loop {
            e -= binom;
            c += 1;
            binom = self.choose(n, c);
            if binom > e {
                break;
            }
        }
        let mut val: u32 = 0;
        loop {
            if e == 0 || {
                binom = self.choose(n - 1, c - 1);
                binom > e
            } {
                c -= 1;
                val |= 1;
            } else {
                e -= binom;
            }
            n -= 1;
            if n == 0 || c == 0 {
                break;
            }
            val <<= 1;
        }
        val <<= n;
        self.bankers[a as usize] = val;
        val
    }
}

/// `bcf_alleles2gt(a, b)` from htslib.
fn alleles2gt(a: i64, b: i64) -> i64 {
    if a > b {
        a * (a + 1) / 2 + b
    } else {
        b * (b + 1) / 2 + a
    }
}

/// Reads the input and returns the GTisec report.
pub fn run(input: &Path, flag: u8, argv_tail: &str) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(&text, flag, argv_tail).map_err(io::Error::other)
}

fn compute(text: &str, flag: u8, argv_tail: &str) -> Result<String, String> {
    let lines: Vec<&str> = text.lines().collect();
    let samples: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();
    if samples.is_empty() {
        return Err("No samples in input file.".to_string());
    }
    let nsmp = samples.len();
    if nsmp > 32 {
        return Err("Too many samples. A maximum of 32 is supported.".to_string());
    }
    let nsmpp2 = 1usize << nsmp;

    let bk = Bankers::new(nsmp);
    let mut smp_is = vec![0u64; nsmpp2];
    let mut missing_gts = vec![0u64; nsmp];

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let Some(gslot) = f[8].split(':').position(|k| k == "GT") else {
            return Err(format!("GT not present at {}:{}", f[0], f[1]));
        };

        // Parse alleles per sample; gte_smp = max ploidy.
        let mut toks: Vec<Vec<&str>> = Vec::with_capacity(nsmp);
        let mut gte_smp = 1usize;
        for s in 0..nsmp {
            let gt = f[9 + s].split(':').nth(gslot).unwrap_or(".");
            let parts: Vec<&str> = gt.split(['/', '|']).collect();
            gte_smp = gte_smp.max(parts.len());
            toks.push(parts);
        }
        if gte_smp > 2 {
            return Err("gtisec does not support ploidy higher than 2.".to_string());
        }

        let parse_allele = |t: &str| -> Option<i64> {
            if t == "." || t.is_empty() {
                None
            } else {
                t.parse::<i64>().ok()
            }
        };

        let mut gts: HashMap<i64, u32> = HashMap::new();
        for (i, parts) in toks.iter().enumerate() {
            let a0 = parts.first().copied().unwrap_or(".");
            let a = parse_allele(a0);
            let second_missing = gte_smp == 2 && {
                match parts.get(1) {
                    Some(t) => parse_allele(t).is_none() && *t == ".",
                    None => false, // vector_end, not missing
                }
            };
            if a.is_none() || second_missing {
                if flag & MISSING != 0 {
                    missing_gts[i] += 1;
                }
                continue;
            }
            let a = a.unwrap();
            let b = if gte_smp == 2 {
                match parts.get(1) {
                    Some(t) => parse_allele(t).unwrap_or(VEND_ALLELE),
                    None => VEND_ALLELE,
                }
            } else {
                VEND_ALLELE
            };
            let key = alleles2gt(a, b);
            *gts.entry(key).or_insert(0) |= 1u32 << i;
        }
        for (_, mask) in gts {
            smp_is[mask as usize] += 1;
        }
    }

    let mut out = String::new();
    out.push_str(&format!(
        "# This file was produced by bcftools +GTisec ({})\n",
        env!("CARGO_PKG_VERSION")
    ));
    out.push_str(&format!(
        "# The command line was:\tbcftools +GTisec {argv_tail}\n"
    ));
    out.push_str(
        "# This file can be used as input to the subset plotting tools at:\n\
         #   https://github.com/dlaehnemann/bankers2\n\
         # Genotype intersections across samples:\n",
    );
    out.push_str("@SMPS");
    for s in (0..nsmp).rev() {
        out.push(' ');
        out.push_str(samples[s]);
    }
    out.push('\n');

    if flag & MISSING != 0 {
        if flag & SMPORDER != 0 {
            out.push_str(
                "# The first line of each sample contains its count of missing genotypes, with a '-' appended\n\
                 #   to the sample name.\n",
            );
        } else {
            out.push_str(&format!(
                "# The first {nsmp} lines contain the counts for missing values of each sample in the order provided\n\
                 #   in the SMPS-line above. Intersection counts only start afterwards.\n"
            ));
        }
    }
    if flag & SMPORDER != 0 {
        out.push_str(
            "# Human readable output (-H) was requested. Subset intersection counts are therefore sorted by\n\
             #   sample and repeated for each contained sample. For each sample, counts are in banker's \n\
             #   sequence order regarding all other samples.\n",
        );
    } else {
        out.push_str("# Subset intersection counts are in global banker's sequence order.\n");
        if nsmp > 2 {
            out.push_str(&format!(
                "#   After exclusive sample counts in order of the SMPS-line, banker's sequence continues with:\n\
                 #   {},{}   {},{}   ...\n",
                samples[nsmp - 1],
                samples[nsmp - 2],
                samples[nsmp - 1],
                samples[nsmp - 3],
            ));
        }
    }
    if flag & VERBOSE != 0 {
        out.push_str(
            "# [1] Number of shared non-ref genotypes \t[2] Samples sharing non-ref genotype (GT)\n",
        );
    } else {
        out.push_str("# [1] Number of shared non-ref genotypes\n");
    }

    let bankers = &bk.bankers;
    if flag & SMPORDER != 0 {
        for s in (0..nsmp).rev() {
            if flag & MISSING != 0 {
                out.push_str(&format!("{}\t{}-\n", missing_gts[s], samples[s]));
            }
            for &bi in &bankers[1..nsmpp2] {
                if (bi >> s) & 1 == 1 {
                    out.push_str(&format!("{}\t{}", smp_is[bi as usize], samples[s]));
                    for j in (0..nsmp).rev() {
                        if (bi ^ (1u32 << s)) & (1u32 << j) != 0 {
                            out.push_str(&format!(",{}", samples[j]));
                        }
                    }
                    out.push('\n');
                }
            }
        }
    } else if flag & VERBOSE != 0 {
        if flag & MISSING != 0 {
            for s in (0..nsmp).rev() {
                out.push_str(&format!("{}\t{}-\n", missing_gts[s], samples[s]));
            }
        }
        for &bi in &bankers[1..nsmpp2] {
            out.push_str(&format!("{}\t", smp_is[bi as usize]));
            let mut started = false;
            for s in (0..nsmp).rev() {
                if (bi >> s) & 1 == 1 {
                    out.push_str(&format!("{}{}", if started { "," } else { "" }, samples[s]));
                    started = true;
                }
            }
            out.push('\n');
        }
    } else {
        if flag & MISSING != 0 {
            for s in (0..nsmp).rev() {
                out.push_str(&format!("{}\n", missing_gts[s]));
            }
        }
        for &bi in &bankers[1..nsmpp2] {
            out.push_str(&format!("{}\n", smp_is[bi as usize]));
        }
    }

    Ok(out)
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
        ".bcftools-rs-gtisec-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bankers_three_samples() {
        // Expected exclusive-then-pair-then-triple ordering for nsmp=3.
        let bk = Bankers::new(3);
        assert_eq!(&bk.bankers[1..8], &[4, 2, 1, 6, 5, 3, 7]);
    }

    #[test]
    fn flags_parse() {
        assert_eq!(parse_flags("H"), SMPORDER | VERBOSE);
        assert_eq!(parse_flags("mv"), MISSING | VERBOSE);
    }
}
