//! `bcftools +trio-switch-rate` (upstream `bcftools/plugins/trio-switch-rate.c`).
//!
//! Calculates the phase-switch rate in trio children (children must have
//! phased genotypes). Trios come from a PED file (`familyID sampleID
//! paternalID maternalID sex phenotype [population]`); results are reported
//! per trio and, when a 7th PED column is present, averaged per population.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

struct Trio {
    father: usize,
    mother: usize,
    child: usize,
    ipop: usize,
    // running state
    prev: i32,
    err: u32,
    nswitch: u32,
    ntest: u32,
}

#[derive(Default)]
struct Pop {
    name: String,
    ntrio: u32,
    err: u64,
    nswitch: u64,
    ntest: u64,
    pswitch: f64,
}

/// Reads the inputs and returns the trio-switch-rate report (with the two
/// `bcftools`-tagged provenance lines the harness strips via
/// `grep -v bcftools`).
pub fn run(input: &Path, ped: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let ped_raw = fs::read_to_string(ped)?;
    compute(&text, &ped_raw).map_err(io::Error::other)
}

/// Parsed genotype: alleles `a`,`b` (only 0/1 considered) and phase.
struct Gt {
    a: i32,
    b: i32,
    phased: bool,
}

/// Mirrors upstream `parse_genotype`: `None` for missing/vector_end or an
/// allele index > 1.
fn parse_gt(s: &str) -> Option<Gt> {
    // a SEP b ; SEP '|' => phased. Missing allele or haploid -> None.
    let bytes = s.as_bytes();
    let sep = bytes.iter().position(|&c| c == b'/' || c == b'|')?;
    let a_tok = &s[..sep];
    let b_tok = &s[sep + 1..];
    if a_tok == "." || b_tok == "." || b_tok.is_empty() {
        return None;
    }
    let a: i32 = a_tok.parse().ok()?;
    let b: i32 = b_tok.parse().ok()?;
    if a > 1 || b > 1 {
        return None;
    }
    Some(Gt {
        a,
        b,
        phased: bytes[sep] == b'|',
    })
}

fn compute(text: &str, ped_raw: &str) -> Result<String, String> {
    let lines: Vec<&str> = text.lines().collect();
    let samples: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();
    let idx_of = |n: &str| samples.iter().position(|s| *s == n);

    let mut trios: Vec<Trio> = Vec::new();
    let mut pops: Vec<Pop> = Vec::new();
    let mut pop2i: HashMap<String, usize> = HashMap::new();
    let mut npop = 0usize;
    for line in ped_raw.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        if cols.len() < 4 {
            return Err(format!("Could not parse the ped file: {line}"));
        }
        let (Some(father), Some(mother), Some(child)) =
            (idx_of(cols[2]), idx_of(cols[3]), idx_of(cols[1]))
        else {
            continue;
        };
        let mut ipop = 0;
        if cols.len() > 6 {
            let name = cols[6];
            ipop = *pop2i.entry(name.to_string()).or_insert_with(|| {
                pops.push(Pop {
                    name: name.to_string(),
                    ..Pop::default()
                });
                let i = npop;
                npop += 1;
                i
            });
            pops[ipop].ntrio += 1;
        }
        trios.push(Trio {
            father,
            mother,
            child,
            ipop,
            prev: 0,
            err: 0,
            nswitch: 0,
            ntest: 0,
        });
    }

    let mut prev_chrom: Option<&str> = None;
    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let chrom = f[0];
        if prev_chrom != Some(chrom) {
            prev_chrom = Some(chrom);
            for t in &mut trios {
                t.prev = 0;
            }
        }
        let Some(gt_slot) = f[8].split(':').position(|k| k == "GT") else {
            continue;
        };
        let sample_cols = &f[9..];
        let gt = |i: usize| -> Option<Gt> {
            sample_cols
                .get(i)
                .and_then(|c| c.split(':').nth(gt_slot))
                .and_then(parse_gt)
        };

        for t in &mut trios {
            let Some(child) = gt(t.child) else { continue };
            if !child.phased {
                continue;
            }
            if child.a + child.b != 1 {
                continue; // child not a het
            }
            let Some(father) = gt(t.father) else { continue };
            let Some(mother) = gt(t.mother) else { continue };
            if father.a + father.b == 1 && mother.a + mother.b == 1 {
                continue; // both parents het
            }
            if father.a + father.b == mother.a + mother.b {
                t.err += 1; // mendelian error
                continue;
            }
            let mut test_phase = 0;
            if father.a == father.b {
                test_phase = 1 + i32::from(child.a == father.a);
            } else if mother.a == mother.b {
                test_phase = 1 + i32::from(child.b == mother.a);
            }
            if t.prev > 0 && t.prev != test_phase {
                t.nswitch += 1;
            }
            t.ntest += 1;
            t.prev = test_phase;
        }
    }

    let mut out = String::new();
    out.push_str(
        "# This file was produced by: bcftools +trio-switch-rate(bcftools-rs+htslib-rs)\n",
    );
    out.push_str("# The command line was:\tbcftools +trio-switch-rate\n#\n");
    out.push_str(
        "# TRIO\t[2]Father\t[3]Mother\t[4]Child\t[5]nTested\t[6]nMendelian Errors\t\
[7]nSwitch\t[8]nSwitch (%)\n",
    );
    for t in &trios {
        let pct = if t.ntest != 0 {
            t.nswitch as f64 * 100.0 / t.ntest as f64
        } else {
            0.0
        };
        out.push_str(&format!(
            "TRIO\t{}\t{}\t{}\t{}\t{}\t{}\t{pct:.2}\n",
            samples[t.father], samples[t.mother], samples[t.child], t.ntest, t.err, t.nswitch,
        ));
        if !pops.is_empty() {
            let p = &mut pops[t.ipop];
            p.err += t.err as u64;
            p.nswitch += t.nswitch as u64;
            p.ntest += t.ntest as u64;
            p.pswitch += pct;
        }
    }
    out.push_str(
        "# POP\tpopulation or other grouping defined by an optional 7-th column of the PED file\n",
    );
    out.push_str(
        "# POP\t[2]Name\t[3]Number of trios\t[4]avgTested\t[5]avgMendelian Errors\t\
[6]avgSwitch\t[7]avgSwitch (%)\n",
    );
    for p in &pops {
        let n = p.ntrio as f64;
        out.push_str(&format!(
            "POP\t{}\t{}\t{:.0}\t{:.0}\t{:.0}\t{:.2}\n",
            p.name,
            p.ntrio,
            p.ntest as f32 as f64 / n,
            p.err as f32 as f64 / n,
            p.nswitch as f32 as f64 / n,
            p.pswitch / n,
        ));
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
        ".bcftools-rs-trio-switch-rate-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_parse() {
        let g = parse_gt("0|1").unwrap();
        assert_eq!((g.a, g.b, g.phased), (0, 1, true));
        let g = parse_gt("1/0").unwrap();
        assert_eq!((g.a, g.b, g.phased), (1, 0, false));
        assert!(parse_gt("./.").is_none());
        assert!(parse_gt("0|2").is_none()); // allele > 1
        assert!(parse_gt("1").is_none()); // haploid
    }

    #[test]
    fn small_trio_report() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tF\tM\n\
20\t1\t.\tT\tA\t.\t.\t.\tGT\t0|1\t1|1\t0|0\n\
20\t2\t.\tT\tA\t.\t.\t.\tGT\t1|0\t1|1\t0|0\n";
        // child het phased; father hom-alt, mother hom-ref -> informative.
        let ped = "fam C F M 2 0 POP1\n";
        let out = compute(vcf, ped).unwrap();
        let trio = out.lines().find(|l| l.starts_with("TRIO")).unwrap();
        // 2 informative sites, no mendelian error; test_phase flips
        // between site1 (child.a==father.a? 0==1 ->1) and site2 (1==1->2)
        // => 1 switch over 2 tested.
        assert_eq!(trio, "TRIO\tF\tM\tC\t2\t0\t1\t50.00");
        let pop = out.lines().find(|l| l.starts_with("POP")).unwrap();
        assert_eq!(pop, "POP\tPOP1\t1\t2\t0\t1\t50.00");
    }
}
