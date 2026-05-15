//! `bcftools +trio-stats` (upstream `bcftools/plugins/trio-stats.c`).
//!
//! Per-trio transmission / Mendelian-error / de-novo statistics for the
//! default (no `-i`/`-e`) "all" filter. Trios come from a PED file. Emits
//! the `MERR`/`TRANSMITTED` debug lines interleaved per record (when
//! `-d mendel-errors,transmitted`), then the `DEF`/`FLT0` summary.

use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const HEADER: &str = "\
# CMD line shows the command line used to generate this output
# DEF lines define expressions for all tested thresholds
# FLT* lines report numbers for every threshold and every trio:
#   1) filter id
#   2) child
#   3) father
#   4) mother
#   5) number of valid trio genotypes (all trio members pass filters, all non-missing)
#   6) number of non-reference trio GTs (at least one trio member carries an alternate allele)
#   7) number of DNMs/Mendelian errors
#   8) number of novel singleton alleles in the child (counted also as DNM / Mendelian error)
#   9) number of untransmitted trio singletons (one alternate allele present in one parent)
#   10) number of transmitted trio singletons (one alternate allele present in one parent and the child)
#   11) number of transitions, all distinct ALT alleles present in the trio are considered
#   12) number of transversions, all distinct ALT alleles present in the trio are considered
#   13) overall ts/tv, all distinct ALT alleles present in the trio are considered
#   14) number of homozygous DNMs/Mendelian errors (likely genotyping errors)
#   15) number of recurrent DNMs/Mendelian errors (non-inherited alleles present in other samples; counts GTs, not sites)
";

#[derive(Default, Clone)]
struct Stats {
    npass: u32,
    nnon_ref: u32,
    nmendel_err: u32,
    nnovel: u32,
    nsingleton: u32,
    ndoubleton: u32,
    nts: u32,
    ntv: u32,
    ndnm_hom: u32,
    ndnm_recurrent: u32,
}

struct Trio {
    child: usize,
    father: usize,
    mother: usize,
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

/// Reads inputs and returns the trio-stats report (the `CMD` line the
/// harness strips with `grep -v ^CMD` is emitted too).
pub fn run(
    input: &Path,
    ped: &Path,
    max_alt_trios: i32,
    dbg_mendel: bool,
    dbg_transmitted: bool,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    let ped_raw = fs::read_to_string(ped)?;
    Ok(compute(
        &text,
        &ped_raw,
        max_alt_trios,
        dbg_mendel,
        dbg_transmitted,
    ))
}

/// Parse a sample GT into two allele indices, mirroring upstream
/// `parse_genotype`: `None` if missing; haploid treated as homozygous.
fn parse_gt(gt: &str) -> Option<[i32; 2]> {
    let toks: Vec<&str> = gt.split(['/', '|']).collect();
    let t0 = *toks.first()?;
    if t0 == "." || t0.is_empty() {
        return None;
    }
    let a0: i32 = t0.parse().ok()?;
    if toks.len() < 2 || toks[1].is_empty() {
        return Some([a0, a0]); // haploid -> hom
    }
    let t1 = toks[1];
    if t1 == "." {
        return None;
    }
    let a1: i32 = t1.parse().ok()?;
    Some([a0, a1])
}

fn compute(
    text: &str,
    ped_raw: &str,
    max_alt_trios: i32,
    dbg_mendel: bool,
    dbg_transmitted: bool,
) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let samples: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();
    let idx_of = |n: &str| samples.iter().position(|s| *s == n);

    let mut trios: Vec<Trio> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    for line in ped_raw.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split_whitespace().collect();
        if c.len() < 4 {
            continue;
        }
        let (Some(father), Some(mother), Some(child)) = (idx_of(c[2]), idx_of(c[3]), idx_of(c[1]))
        else {
            continue;
        };
        let key = format!("{} {} {}", c[1], c[2], c[3]);
        if !seen.insert(key) {
            continue;
        }
        trios.push(Trio {
            child,
            father,
            mother,
        });
    }

    let mut stats = vec![Stats::default(); trios.len()];
    let mut out = String::new();
    out.push_str(HEADER);
    out.push_str("CMD\tbcftools +trio-stats\n");

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        let chrom = f[0];
        let pos = f[1];
        let reference = f[3];
        let alleles: Vec<&str> = std::iter::once(reference)
            .chain(if f[4] == "." {
                "".split(',').take(0)
            } else {
                f[4].split(',').take(usize::MAX)
            })
            .collect();
        let n_allele = alleles.len();

        let gt_slot = f[8].split(':').position(|k| k == "GT");
        let scol = &f[9..];
        let gt_of = |s: usize| -> Option<[i32; 2]> {
            let g = match gt_slot {
                Some(i) => scol.get(s).and_then(|c| c.split(':').nth(i)).unwrap_or("."),
                None => ".",
            };
            parse_gt(g)
        };

        // bcf_calc_ac over all samples (INFO/AC+AN else GT tally).
        let ac = calc_ac(&f, n_allele, gt_slot);

        let ref_code = if reference.len() == 1 {
            acgt2int(reference.as_bytes()[0])
        } else {
            -1
        };
        let star_allele: i32 = (1..n_allele)
            .find(|&i| alleles[i] == "*")
            .map(|i| i as i32)
            .unwrap_or(-1);

        // -a: per-allele accumulation across trios at this site.
        // alt_acc[allele] = (nalt, Vec<(trio, is_singleton)>)
        let mut alt_acc: Vec<(i32, Vec<(usize, bool)>)> = vec![(0, Vec::new()); n_allele];

        for (ti, trio) in trios.iter().enumerate() {
            let (Some(ch), Some(fa), Some(mo)) =
                (gt_of(trio.child), gt_of(trio.father), gt_of(trio.mother))
            else {
                continue;
            };
            let st = &mut stats[ti];
            let als = [ch[0], ch[1], fa[0], fa[1], mo[0], mo[1]];
            st.npass += 1;

            let mut has_star = false;
            let mut has_nonref = false;
            let mut ac_trio = vec![0i32; n_allele];
            for &a in &als {
                if a == star_allele {
                    has_star = true;
                    continue;
                }
                if a != 0 {
                    has_nonref = true;
                }
                if (a as usize) < n_allele {
                    ac_trio[a as usize] += 1;
                }
            }
            if !has_nonref {
                continue;
            }
            st.nnon_ref += 1;

            if ref_code != -1 {
                let mut has_ts = false;
                let mut has_tv = false;
                for &a in &als {
                    if a == 0 || a == star_allele {
                        continue;
                    }
                    let ai = a as usize;
                    if ai >= n_allele || alleles[ai].len() != 1 {
                        continue;
                    }
                    let alt = acgt2int(alleles[ai].as_bytes()[0]);
                    if (ref_code - alt).abs() == 2 {
                        has_ts = true;
                    } else {
                        has_tv = true;
                    }
                }
                if has_ts {
                    st.nts += 1;
                }
                if has_tv {
                    st.ntv += 1;
                }
            }

            if has_star {
                continue;
            }

            // Mendelian error
            let a0f = ch[0] == fa[0] || ch[0] == fa[1];
            let a1m = ch[1] == mo[0] || ch[1] == mo[1];
            if !a0f || !a1m {
                let a0m = ch[0] == mo[0] || ch[0] == mo[1];
                let a1f = ch[1] == fa[0] || ch[1] == fa[1];
                if !a0m || !a1f {
                    st.nmendel_err += 1;
                    let dnm_hom = ch[0] == ch[1];
                    if dnm_hom {
                        st.ndnm_hom += 1;
                    }
                    let culprit = if !a0f && !a0m {
                        ch[0]
                    } else if !a1f && !a1m {
                        ch[1]
                    } else if ac[ch[0] as usize] < ac[ch[1] as usize] {
                        ch[0]
                    } else {
                        ch[1]
                    };
                    let acc = ac[culprit as usize];
                    let dnm_recurrent = (!dnm_hom && acc > 1) || (dnm_hom && acc > 2);
                    if dnm_recurrent {
                        st.ndnm_recurrent += 1;
                    }
                    if dbg_mendel {
                        out.push_str(&format!(
                            "MERR\t{chrom}\t{pos}\t{}\t{}\t{}\t{}\t{}\n",
                            samples[trio.child],
                            samples[trio.father],
                            samples[trio.mother],
                            if dnm_hom { "HOM" } else { "-" },
                            if dnm_recurrent { "RECURRENT" } else { "-" },
                        ));
                    }
                }
            }

            for j in 0..n_allele {
                if ac_trio[j] == 0 {
                    continue;
                }
                if max_alt_trios != 0 {
                    alt_acc[j].0 += 1;
                }
                if ac_trio[j] == 1 {
                    if ch[0] == j as i32 || ch[1] == j as i32 {
                        st.nnovel += 1;
                    } else if max_alt_trios == 0 {
                        st.nsingleton += 1;
                        if dbg_transmitted {
                            out.push_str(&transmitted(
                                chrom,
                                pos,
                                samples[trio.child],
                                samples[trio.father],
                                samples[trio.mother],
                                false,
                            ));
                        }
                    } else {
                        alt_acc[j].1.push((ti, true));
                    }
                } else if ac_trio[j] == 2 {
                    let cj = ch[0] == j as i32 || ch[1] == j as i32;
                    let chom = ch[0] == j as i32 && ch[1] == j as i32;
                    if !cj || chom {
                        continue;
                    }
                    if (fa[0] == j as i32 && fa[1] == j as i32)
                        || (mo[0] == j as i32 && mo[1] == j as i32)
                    {
                        continue;
                    }
                    if max_alt_trios == 0 {
                        st.ndoubleton += 1;
                        if dbg_transmitted {
                            out.push_str(&transmitted(
                                chrom,
                                pos,
                                samples[trio.child],
                                samples[trio.father],
                                samples[trio.mother],
                                true,
                            ));
                        }
                    } else {
                        alt_acc[j].1.push((ti, false));
                    }
                }
            }
        }

        if max_alt_trios != 0 {
            for (nalt, recs) in alt_acc.into_iter() {
                if recs.is_empty() || nalt > max_alt_trios {
                    continue;
                }
                for (ti, is_singleton) in recs {
                    let trio = &trios[ti];
                    if is_singleton {
                        stats[ti].nsingleton += 1;
                        if dbg_transmitted {
                            out.push_str(&transmitted(
                                chrom,
                                pos,
                                samples[trio.child],
                                samples[trio.father],
                                samples[trio.mother],
                                false,
                            ));
                        }
                    } else {
                        stats[ti].ndoubleton += 1;
                        if dbg_transmitted {
                            out.push_str(&transmitted(
                                chrom,
                                pos,
                                samples[trio.child],
                                samples[trio.father],
                                samples[trio.mother],
                                true,
                            ));
                        }
                    }
                }
            }
        }
    }

    out.push_str("DEF\tFLT0\tall\n");
    for (ti, trio) in trios.iter().enumerate() {
        let s = &stats[ti];
        let tstv = if s.ntv != 0 {
            format!("{:.2}", s.nts as f32 / s.ntv as f32)
        } else {
            "inf".to_owned()
        };
        out.push_str(&format!(
            "FLT0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{tstv}\t{}\t{}\n",
            samples[trio.child],
            samples[trio.father],
            samples[trio.mother],
            s.npass,
            s.nnon_ref,
            s.nmendel_err,
            s.nnovel,
            s.nsingleton,
            s.ndoubleton,
            s.nts,
            s.ntv,
            s.ndnm_hom,
            s.ndnm_recurrent,
        ));
    }
    out
}

fn transmitted(chrom: &str, pos: &str, c: &str, f: &str, m: &str, yes: bool) -> String {
    format!(
        "TRANSMITTED\t{chrom}\t{pos}\t{c}\t{f}\t{m}\t{}\n",
        if yes { "YES" } else { "NO" }
    )
}

/// `bcf_calc_ac` over all samples: INFO/AC+AN when present, otherwise
/// tally FORMAT/GT (count of each allele across every sample).
fn calc_ac(f: &[&str], n_allele: usize, gt_slot: Option<usize>) -> Vec<i32> {
    if f[7] != "." {
        let mut an: Option<i64> = None;
        let mut acv: Option<Vec<i64>> = None;
        for kv in f[7].split(';') {
            let mut it = kv.splitn(2, '=');
            match (it.next(), it.next()) {
                (Some("AN"), Some(v)) => an = v.parse().ok(),
                (Some("AC"), Some(v)) => {
                    acv = v.split(',').map(|x| x.parse::<i64>().ok()).collect()
                }
                _ => {}
            }
        }
        if let (Some(an), Some(acv)) = (an, acv)
            && acv.len() == n_allele - 1
        {
            let mut ac = vec![0i32; n_allele];
            ac[0] = (an - acv.iter().sum::<i64>()) as i32;
            for (i, v) in acv.into_iter().enumerate() {
                ac[i + 1] = v as i32;
            }
            return ac;
        }
    }
    let mut ac = vec![0i32; n_allele];
    for s in &f[9..] {
        let gt = match gt_slot {
            Some(i) => s.split(':').nth(i).unwrap_or("."),
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
    ac
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
        ".bcftools-rs-trio-stats-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_parse() {
        assert_eq!(parse_gt("0/1"), Some([0, 1]));
        assert_eq!(parse_gt("1|1"), Some([1, 1]));
        assert_eq!(parse_gt("1"), Some([1, 1])); // haploid -> hom
        assert_eq!(parse_gt("./."), None);
    }

    #[test]
    fn ac_from_gt() {
        let f: Vec<&str> = "1\t1\t.\tT\tA\t.\t.\t.\tGT\t0/1\t1/1".split('\t').collect();
        assert_eq!(calc_ac(&f, 2, Some(0)), vec![1, 3]);
    }
}
