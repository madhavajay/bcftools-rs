//! `bcftools +mendelian2` (upstream `bcftools/plugins/mendelian2.c`).
//!
//! Mendelian-consistency checking for a single `-p [1X:|2X:]P,F,M` trio.
//! Like upstream, the built-in GRCh37 ruleset is active by default
//! (`init_rules(args, NULL)` → alias `"GRCh37"`), so chrX/Y/MT use the
//! haploid ploidy-1 inheritance pattern; all other regions inherit from
//! both parents (MF, ploidy 2). Modes: `c` count (default), `a` annotate
//! INFO/MERR, `d` delete inconsistent trio GTs, `e`/`g`/`m` list
//! error/good/missing sites (`E`/`M`/`S` drop variants). The explicit
//! `--rules`/`--rules-file` engine and `-i`/`-e` filtering need
//! infrastructure tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::vcf_compat::normalize_vcf_text;

const MODE_ANNOTATE: u32 = 1 << 0;
const MODE_COUNT: u32 = 1 << 1;
const MODE_DELETE: u32 = 1 << 2;
const MODE_LIST_ERR: u32 = 1 << 3;
const MODE_DROP_ERR: u32 = 1 << 4;
const MODE_LIST_GOOD: u32 = 1 << 5;
const MODE_LIST_MISS: u32 = 1 << 6;
const MODE_DROP_MISS: u32 = 1 << 7;
const MODE_LIST_SKIP: u32 = 1 << 8;
const MODE_DROP_SKIP: u32 = 1 << 9;
const LIST_MODES: u32 = MODE_LIST_ERR | MODE_LIST_GOOD | MODE_LIST_MISS | MODE_LIST_SKIP;

const HAS_GOOD: u32 = 1;
const HAS_MERR: u32 = 2;
const HAS_MISS: u32 = 4;

// Sex-id bit indices, mirroring upstream `#define iMOM 0`, `#define iDAD 1`.
const I_MOM: usize = 0;
const I_DAD: usize = 1;
// Number of sexes; upstream builds str2sex_id={"1X":0,"2X":1} from the
// GRCh37 rules string, so `1X`→0 (male chrX pattern), `2X`→1.
const NSEX_ID: usize = 2;

/// Built-in `GRCh37` ruleset (the default `init_rules(args, NULL)` alias).
/// Tuples are `(sex_id, chrom, beg0, end0, inherits, ploidy)` with 0-based
/// inclusive coordinates (upstream `parse_rules` subtracts 1). `inherits`
/// is a bitmask of `1<<I_MOM` (`M`) / `1<<I_DAD` (`F`); `.` clears both.
/// Regions not covered here inherit from both parents (MF, ploidy 2).
const GRCH37_RULES: &[(usize, &str, u64, u64, u32, u32)] = &[
    (0, "X", 0, 59999, 1 << I_MOM, 1),
    (0, "X", 2699520, 154931042, 1 << I_MOM, 1),
    (0, "Y", 0, 59373565, 1 << I_DAD, 1),
    (1, "Y", 0, 59373565, 0, 0),
    (0, "MT", 0, 16568, 1 << I_MOM, 1),
    (1, "MT", 0, 16568, 1 << I_MOM, 1),
    (0, "chrX", 0, 59999, 1 << I_MOM, 1),
    (0, "chrX", 2699520, 154931042, 1 << I_MOM, 1),
    (0, "chrY", 0, 59373565, 1 << I_DAD, 1),
    (1, "chrY", 0, 59373565, 0, 0),
    (0, "chrM", 0, 16568, 1 << I_MOM, 1),
    (1, "chrM", 0, 16568, 1 << I_MOM, 1),
];

/// Resolves the per-`sex_id` `(inherits, ploidy)` rule for one record by
/// overlapping `[pos0, end0]` against [`GRCH37_RULES`]; unlisted regions
/// keep the MF/ploidy-2 default (upstream `collect_stats`).
fn resolve_rules(chrom: &str, pos0: u64, end0: u64) -> [(u32, u32); NSEX_ID] {
    let mut rule = [((1 << I_MOM) | (1 << I_DAD), 2u32); NSEX_ID];
    for &(sid, rchr, rb, re, inh, plo) in GRCH37_RULES {
        if rchr == chrom && pos0 <= re && end0 >= rb {
            rule[sid] = (inh, plo);
        }
    }
    rule
}

/// Parses the `-m` mode string into the bitmask (`c` default).
pub fn parse_mode(s: &str) -> Result<u32, String> {
    let mut m = 0;
    for ch in s.chars() {
        m |= match ch {
            'g' | '+' => MODE_LIST_GOOD,
            'x' | 'e' => MODE_LIST_ERR,
            'a' => MODE_ANNOTATE,
            'd' => MODE_DELETE,
            'c' => MODE_COUNT,
            'E' => MODE_DROP_ERR,
            'u' | 'm' => MODE_LIST_MISS,
            'M' => MODE_DROP_MISS,
            's' => MODE_LIST_SKIP,
            'S' => MODE_DROP_SKIP,
            other => return Err(format!("Unknown -m mode: {other}")),
        };
    }
    if m == 0 {
        m = MODE_COUNT;
    }
    Ok(m)
}

#[derive(Default)]
struct TrioStats {
    ngood: u32,
    ngood_alt: u32,
    nmerr: u32,
    nmiss: u32,
    nfail: u32,
}

#[derive(Default)]
struct SiteStats {
    nref_only: u32,
    nmany_als: u32,
    nfail: u32,
    nno_gt: u32,
    nnot_diploid: u32,
    nmiss: u32,
    nmerr: u32,
    ngood: u32,
}

/// `-p [1X:|2X:]P,F,M` → `(proband, father, mother, sex_id)` sample
/// names. Mirrors upstream `init_data`: no prefix leaves `sex_id` at the
/// `calloc` default `0` (str2sex `1X`, male chrX pattern); `1X:` sets it
/// to `iDAD` (1) and `2X:` to `iMOM` (0).
fn parse_pfm(spec: &str) -> Result<(String, String, String, usize), String> {
    let (s, sex_id) = if let Some(rest) = spec.strip_prefix("1X:") {
        (rest, I_DAD)
    } else if let Some(rest) = spec.strip_prefix("2X:") {
        (rest, I_MOM)
    } else {
        (spec, 0)
    };
    let parts: Vec<&str> = s.split(',').collect();
    if parts.len() != 3 {
        return Err(format!("Could not parse -p {spec}"));
    }
    Ok((
        parts[0].to_string(),
        parts[1].to_string(),
        parts[2].to_string(),
        sex_id,
    ))
}

/// Reads the input and returns the mendelian2 output (VCF for the
/// list/annotate/delete modes, the stats table for count mode).
pub fn run(input: &Path, pfm: &str, mode: u32) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(&text, pfm, mode).map_err(io::Error::other)
}

/// Allele bitmask parse mirroring upstream `parse_gt`: returns
/// `(a_mask, b_mask, nal)`; for the OR-form pass `or = true`.
fn parse_gt(gt: &str) -> (u64, u64, i32) {
    let toks: Vec<&str> = gt.split(['/', '|']).collect();
    let t0 = toks.first().copied().unwrap_or(".");
    if t0 == "." || t0.is_empty() {
        return (0, 0, 0);
    }
    let Ok(a0) = t0.parse::<u32>() else {
        return (0, 0, 0);
    };
    let amask = 1u64 << a0;
    if toks.len() < 2 || toks[1].is_empty() {
        return (amask, 0, 1);
    }
    let t1 = toks[1];
    if t1 == "." {
        return (0, 0, 0);
    }
    let Ok(a1) = t1.parse::<u32>() else {
        return (0, 0, 0);
    };
    (amask, 1u64 << a1, 2)
}

/// OR-form (`parse_gt(...,&x,&x)`): returns `(mask, nal)`.
fn parse_gt_or(gt: &str) -> (u64, i32) {
    let (a, b, n) = parse_gt(gt);
    (a | b, n)
}

struct Trio {
    kid: usize,
    dad: usize,
    mom: usize,
    sex_id: usize,
}

fn compute(text: &str, pfm: &str, mode: u32) -> Result<String, String> {
    let (p, fa, mo, sex_id) = parse_pfm(pfm)?;
    let lines: Vec<&str> = text.lines().collect();
    let samples: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();
    let idx = |n: &str| -> Result<usize, String> {
        samples
            .iter()
            .position(|s| *s == n)
            .ok_or_else(|| format!("No such sample: \"{n}\""))
    };
    let trio = Trio {
        kid: idx(&p)?,
        dad: idx(&fa)?,
        mom: idx(&mo)?,
        sex_id,
    };

    let mut site = SiteStats::default();
    let mut ts = TrioStats::default();

    let mut out = String::new();

    // Header (VCF-output modes only).
    let vcf_out = mode != MODE_COUNT;
    if vcf_out {
        let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
        let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
        for (i, line) in lines.iter().enumerate() {
            if !line.starts_with('#') {
                break;
            }
            if line.starts_with("#CHROM") {
                if mode & MODE_ANNOTATE != 0 {
                    out.push_str("##INFO=<ID=MERR,Number=1,Type=Integer,Description=\"Number of trios with Mendelian errors\">\n");
                }
                out.push_str(line);
                out.push('\n');
                continue;
            }
            out.push_str(line);
            out.push('\n');
            if Some(i) == fileformat && !has_pass {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            }
        }
    }

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if f.len() < 10 {
            continue;
        }
        let reference = f[3].clone();
        let alts: Vec<&str> = if f[4] == "." {
            Vec::new()
        } else {
            f[4].split(',').collect()
        };
        let n_allele = 1 + alts.len();

        // Skip-site checks.
        let vt = alts.iter().fold(VariantType::REF, |a, x| {
            a | classify_variant(&reference, x).variant_type
        });
        let mut skip = false;
        if n_allele == 1 || vt == VariantType::REF {
            site.nref_only += 1;
            skip = true;
        } else if n_allele > 64 {
            site.nmany_als += 1;
            skip = true;
        }
        if skip {
            if mode & MODE_DROP_SKIP != 0 && vcf_out {
                out.push_str(&f.join("\t"));
                out.push('\n');
            }
            continue;
        }

        let gt_slot = f[8].split(':').position(|k| k == "GT");
        let scol = &f[9..];
        let raw_gt = |i: usize| -> String {
            match gt_slot {
                Some(g) => scol
                    .get(i)
                    .and_then(|c| c.split(':').nth(g))
                    .unwrap_or(".")
                    .to_string(),
                None => ".".to_string(),
            }
        };

        // collect_stats for the single trio, applying the per-record
        // GRCh37 (sex_id, inherits, ploidy) rule (upstream default).
        let mut nmerr_rec = 0i32;
        let mut ret = 0u32;
        let mut has_merr = false;
        {
            let pos0: u64 = f[1].parse::<u64>().unwrap_or(1).saturating_sub(1);
            let rlen = f[3].len().max(1) as u64;
            let end0 = pos0 + rlen - 1;
            let rules = resolve_rules(&f[0], pos0, end0);
            let (inherits, ploidy) = rules[trio.sex_id];
            if inherits == 0 {
                // upstream `if ( !rule->inherits ) { nrule++; continue; }`
                // — the trio is skipped for stats (no HAS_* bit, no site
                // counter); `ret` stays 0 so downstream modes behave as
                // for a record with no good/err/miss trio.
            } else {
                let kid_s = raw_gt(trio.kid);
                let (kid1, kid2, nal_kid) = parse_gt(&kid_s);
                if (nal_kid as u32) < ploidy {
                    ret |= HAS_MISS;
                    ts.nmiss += 1;
                } else if nal_kid == 1 {
                    // Haploid kid: compare against the single inheriting
                    // parent (first set bit of `inherits`).
                    let j = (0..NSEX_ID)
                        .find(|&j| inherits & (1 << j) != 0)
                        .unwrap_or(I_MOM);
                    let parent_smpl = if j == I_MOM { trio.mom } else { trio.dad };
                    let (parent, nal_p) = parse_gt_or(&raw_gt(parent_smpl));
                    if nal_p == 0 {
                        ret |= HAS_MISS;
                        ts.nmiss += 1;
                    } else if parent & kid1 != 0 {
                        ret |= HAS_GOOD;
                        ts.ngood += 1;
                        if !(parent == 1 && kid1 == 1) {
                            ts.ngood_alt += 1;
                        }
                    } else {
                        ret |= HAS_MERR;
                        ts.nmerr += 1;
                        has_merr = true;
                        nmerr_rec += 1;
                    }
                } else {
                    let (mom, nal_mom) = parse_gt_or(&raw_gt(trio.mom));
                    let (dad, nal_dad) = parse_gt_or(&raw_gt(trio.dad));
                    if (kid1 & dad != 0 && kid2 & mom != 0) || (kid1 & mom != 0 && kid2 & dad != 0)
                    {
                        ret |= HAS_GOOD;
                        ts.ngood += 1;
                        if dad != 1 || mom != 1 || (kid1 | kid2) != 1 {
                            ts.ngood_alt += 1;
                        }
                    } else {
                        let mom_miss = nal_mom == 0;
                        let dad_miss = nal_dad == 0;
                        if mom_miss || dad_miss {
                            ret |= HAS_MISS;
                            ts.nmiss += 1;
                        }
                        // Not an error if both parents missing, or one
                        // parent missing while the kid is consistent with
                        // the other (upstream's three `continue` guards).
                        let consistent = (mom_miss && (dad_miss || (kid1 | kid2) & dad != 0))
                            || (dad_miss && (kid1 | kid2) & mom != 0);
                        if !consistent {
                            ret |= HAS_MERR;
                            ts.nmerr += 1;
                            has_merr = true;
                            nmerr_rec += 1;
                        }
                    }
                }
            }
        }

        if ret & HAS_MERR != 0 {
            site.nmerr += 1;
        }
        if ret & HAS_MISS != 0 {
            site.nmiss += 1;
        }
        if ret & HAS_GOOD != 0 {
            site.ngood += 1;
        }

        if mode & MODE_COUNT != 0 {
            continue;
        }
        if mode & MODE_DROP_ERR != 0 && ret & HAS_MERR != 0 {
            continue;
        }
        if mode & MODE_DROP_MISS != 0 && ret & HAS_MISS != 0 {
            continue;
        }

        if mode & MODE_DELETE != 0
            && ret & HAS_MERR != 0
            && has_merr
            && let Some(g) = gt_slot
        {
            for s in [trio.kid, trio.dad, trio.mom] {
                let fi = 9 + s;
                if fi >= f.len() {
                    continue;
                }
                let mut parts: Vec<String> = f[fi].split(':').map(|x| x.to_string()).collect();
                if g < parts.len() {
                    parts[g] = delete_gt(&parts[g]);
                }
                f[fi] = parts.join(":");
            }
        }

        if mode & MODE_ANNOTATE != 0 {
            let info = if f[7] == "." || f[7].is_empty() {
                format!("MERR={nmerr_rec}")
            } else {
                format!("{};MERR={nmerr_rec}", f[7])
            };
            f[7] = info;
        }

        if mode & LIST_MODES != 0 {
            let keep = (mode & MODE_LIST_ERR != 0 && ret & HAS_MERR != 0)
                || (mode & MODE_LIST_MISS != 0 && ret & HAS_MISS != 0)
                || (mode & MODE_LIST_GOOD != 0 && ret & HAS_GOOD != 0);
            if keep {
                out.push_str(&f.join("\t"));
                out.push('\n');
            }
            continue;
        }
        out.push_str(&f.join("\t"));
        out.push('\n');
    }

    if mode & MODE_COUNT != 0 {
        // count mode == exactly MODE_COUNT -> stdout
        out.push_str("# Summary stats\n");
        out.push_str(&format!(
            "sites_ref_only\t{}\t# sites skipped because there was no ALT allele\n",
            site.nref_only
        ));
        out.push_str(&format!(
            "sites_many_als\t{}\t# skipped because of too many ALT alleles\n",
            site.nmany_als
        ));
        out.push_str(&format!(
            "sites_fail\t{}\t# skipped because of failed -i/-e filter\n",
            site.nfail
        ));
        out.push_str(&format!(
            "sites_no_GT\t{}\t# skipped because of absent FORMAT/GT field\n",
            site.nno_gt
        ));
        out.push_str(&format!(
            "sites_not_diploid\t{}\t# skipped because FORMAT/GT not formatted diploid\n",
            site.nnot_diploid
        ));
        out.push_str(&format!(
            "sites_missing\t{}\t# number of sites with at least one trio GT missing\n",
            site.nmiss
        ));
        out.push_str(&format!(
            "sites_merr\t{}\t# number of sites with at least one Mendelian error\n",
            site.nmerr
        ));
        out.push_str(&format!(
            "sites_good\t{}\t# number of sites with at least one good trio\n",
            site.ngood
        ));
        out.push_str(
            "# Per-trio stats, each column corresponds to one trio. List of trios is below.\n",
        );
        out.push_str(
            "# The meaning of per-trio stats is the same as described above, ngood_alt is\n",
        );
        out.push_str(
            "# the number of good genotypes with at least one non-reference allele, and is\n",
        );
        out.push_str("# included in the ngood counter\n");
        out.push_str(&format!("ngood\t{}\n", ts.ngood));
        out.push_str(&format!("ngood_alt\t{}\n", ts.ngood_alt));
        out.push_str(&format!("nmerr\t{}\n", ts.nmerr));
        out.push_str(&format!("nmissing\t{}\n", ts.nmiss));
        out.push_str(&format!("nfail\t{}\n", ts.nfail));
        out.push_str("# List of trios. Their ids are in the same order as the values listed in the stats lines above. For\n");
        out.push_str("# example, the values for the first trio (id=1) and the third trio (id=3) are in the 2nd and the 4th\n");
        out.push_str("# column and their stats can be obtained with the unix command\n");
        out.push_str("#     cat stats.txt | grep ^n | cut -f1,2,4\n");
        out.push_str("# TRIO\t[2]id\t[3]child\t[4]father\t[5]mother\n");
        out.push_str(&format!(
            "TRIO\t1\t{}\t{}\t{}\n",
            samples[trio.kid], samples[trio.dad], samples[trio.mom]
        ));
    }

    Ok(out)
}

/// Set every allele in a GT to missing, preserving separators (e.g.
/// `0/1`→`./.`, `0|1`→`.|.`, `.`→`.`).
fn delete_gt(gt: &str) -> String {
    gt.split_inclusive(['/', '|'])
        .map(|seg| match seg.chars().last() {
            Some(c @ ('/' | '|')) => format!(".{c}"),
            _ => ".".to_string(),
        })
        .collect()
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
        ".bcftools-rs-mendelian2-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pfm_parse() {
        assert_eq!(
            parse_pfm("child1,dad1,mom1").unwrap(),
            ("child1".into(), "dad1".into(), "mom1".into(), 0)
        );
        assert_eq!(
            parse_pfm("1X:c,f,m").unwrap(),
            ("c".into(), "f".into(), "m".into(), I_DAD)
        );
        assert_eq!(
            parse_pfm("2X:c,f,m").unwrap(),
            ("c".into(), "f".into(), "m".into(), I_MOM)
        );
    }

    #[test]
    fn gt_bitmask() {
        assert_eq!(parse_gt("0/1"), (1, 2, 2));
        assert_eq!(parse_gt("1/1"), (2, 2, 2));
        assert_eq!(parse_gt("./."), (0, 0, 0));
        assert_eq!(parse_gt_or("0/0"), (1, 2));
        assert_eq!(parse_gt_or("."), (0, 0));
    }

    #[test]
    fn delete_keeps_sep() {
        assert_eq!(delete_gt("0/0"), "./.");
        assert_eq!(delete_gt("0|1"), ".|.");
        assert_eq!(delete_gt("."), ".");
        assert_eq!(delete_gt("10/1"), "./.");
    }

    #[test]
    fn mode_default_count() {
        assert_eq!(parse_mode("").unwrap(), MODE_COUNT);
        assert_eq!(parse_mode("d").unwrap(), MODE_DELETE);
        assert_eq!(parse_mode("e").unwrap(), MODE_LIST_ERR);
    }
}
