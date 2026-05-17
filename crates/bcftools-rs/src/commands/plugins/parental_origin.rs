//! `bcftools +parental-origin` (upstream `bcftools/plugins/parental-origin.c`).
//!
//! Determines the parental origin of a CNV (deletion or duplication)
//! region in a trio from FORMAT/PL, FORMAT/AD and FORMAT/GT, accumulating
//! a log-likelihood of paternal vs maternal origin over the filtered
//! informative SNP sites inside `-r REGION`. Common `-i`/`-e` filters route
//! through the shared text filter engine before likelihood accumulation.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::math::kf_lgamma;
use htslib_rs::variant::{VariantType, classify_variant};

use crate::filter::{self as bcffilter, EvalContext};
use crate::vcf_compat::normalize_vcf_text;

const KF_GAMMA_EPS: f64 = 1e-14;
const KF_TINY: f64 = 1e-290;

/// CNV type selected with `-t del|dup`.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum CnvType {
    Del,
    Dup,
}

impl CnvType {
    pub fn parse(s: &str) -> Option<CnvType> {
        if s.eq_ignore_ascii_case("dup") {
            Some(CnvType::Dup)
        } else if s.eq_ignore_ascii_case("del") {
            Some(CnvType::Del)
        } else {
            None
        }
    }
    fn label(self) -> &'static str {
        match self {
            CnvType::Del => "del",
            CnvType::Dup => "dup",
        }
    }
}

#[derive(Clone, Copy)]
pub enum FilterMode {
    Include,
    Exclude,
}

#[derive(Clone, Copy)]
pub struct FilterSpec<'a> {
    pub mode: FilterMode,
    pub expr: &'a str,
}

/// Port of HTSlib `kfunc.c` `kf_betai_aux` (modified Lentz continued
/// fraction for the regularized incomplete beta function).
fn kf_betai_aux(a: f64, b: f64, x: f64) -> f64 {
    if x == 0.0 {
        return 0.0;
    }
    if x == 1.0 {
        return 1.0;
    }
    let mut f = 1.0;
    let mut c = f;
    let mut d = 0.0;
    for j in 1..200 {
        let m = (j >> 1) as f64;
        let aa = if j & 1 == 1 {
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        } else {
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        };
        d = 1.0 + aa * d;
        if d < KF_TINY {
            d = KF_TINY;
        }
        c = 1.0 + aa / c;
        if c < KF_TINY {
            c = KF_TINY;
        }
        d = 1.0 / d;
        let dd = c * d;
        f *= dd;
        if (dd - 1.0).abs() < KF_GAMMA_EPS {
            break;
        }
    }
    (kf_lgamma(a + b) - kf_lgamma(a) - kf_lgamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp() / a / f
}

/// Port of HTSlib `kf_betai` (regularized incomplete beta `I_x(a,b)`).
fn kf_betai(a: f64, b: f64, x: f64) -> f64 {
    if x < (a + 1.0) / (a + b + 2.0) {
        kf_betai_aux(a, b, x)
    } else {
        1.0 - kf_betai_aux(b, a, 1.0 - x)
    }
}

/// `bcftools.h` `calc_binom_two_sided`.
fn calc_binom_two_sided(na: i32, nb: i32, aprob: f64) -> f64 {
    if na == 0 && nb == 0 {
        return -1.0;
    }
    if na == nb {
        return 1.0;
    }
    let prob = if na > nb {
        2.0 * kf_betai(na as f64, (nb + 1) as f64, aprob)
    } else {
        2.0 * kf_betai(nb as f64, (na + 1) as f64, aprob)
    };
    prob.min(1.0)
}

/// `bcftools.h` `calc_binom_one_sided`.
fn calc_binom_one_sided(na: i32, nb: i32, aprob: f64, ge: bool) -> f64 {
    if ge {
        kf_betai(na as f64, (nb + 1) as f64, aprob)
    } else {
        kf_betai(nb as f64, (na + 1) as f64, 1.0 - aprob)
    }
}

struct Region {
    chrom: String,
    beg: u64,
    end: u64,
}

/// Parses `chrom`, `chrom:pos`, or `chrom:beg-end` (1-based inclusive).
fn parse_region(spec: &str) -> Region {
    match spec.split_once(':') {
        None => Region {
            chrom: spec.to_string(),
            beg: 1,
            end: u64::MAX,
        },
        Some((chrom, rng)) => {
            let strip = |s: &str| s.replace(',', "");
            let (beg, end) = match rng.split_once('-') {
                None => {
                    let b = strip(rng).parse().unwrap_or(1);
                    (b, b)
                }
                Some((b, e)) => {
                    let b = strip(b).parse().unwrap_or(1);
                    let e = if e.is_empty() {
                        u64::MAX
                    } else {
                        strip(e).parse().unwrap_or(u64::MAX)
                    };
                    (b, e)
                }
            };
            Region {
                chrom: chrom.to_string(),
                beg,
                end,
            }
        }
    }
}

/// Reads the input and returns the parental-origin report.
#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &Path,
    region: &str,
    pfm: &str,
    cnv_type: CnvType,
    greedy: bool,
    min_pbinom: f64,
    argv_tail: &str,
    filter: Option<FilterSpec<'_>>,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(
        &text, region, pfm, cnv_type, greedy, min_pbinom, argv_tail, filter,
    )
    .map_err(io::Error::other)
}

#[allow(clippy::too_many_arguments)]
fn compute(
    text: &str,
    region: &str,
    pfm: &str,
    cnv_type: CnvType,
    greedy: bool,
    min_pbinom: f64,
    argv_tail: &str,
    filter: Option<FilterSpec<'_>>,
) -> Result<String, String> {
    let reg = parse_region(region);
    let parts: Vec<&str> = pfm.split(',').collect();
    if parts.len() != 3 {
        return Err(format!("Expected three sample names with -p: {pfm}"));
    }

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
            .ok_or_else(|| format!("The sample is not present: {n}"))
    };
    // iCHILD=0 (proband), iFATHER=1, iMOTHER=2.
    let trio = [idx(parts[0])?, idx(parts[1])?, idx(parts[2])?];

    let mut ppat = 0.0f64; // accumulated log-prob of paternal origin
    let mut pmat = 0.0f64;
    let mut ntest = 0i32;

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        if let Some(filter) = filter
            && !record_passes_filter(&f, filter)?
        {
            continue;
        }
        if f[0] != reg.chrom {
            continue;
        }
        let pos: u64 = match f[1].parse() {
            Ok(p) => p,
            Err(_) => continue,
        };
        if pos < reg.beg || pos > reg.end {
            continue;
        }
        let reference = f[3];
        let alts: Vec<&str> = if f[4] == "." {
            Vec::new()
        } else {
            f[4].split(',').collect()
        };
        if alts.len() != 1 {
            continue; // rec->n_allele != 2
        }
        if classify_variant(reference, alts[0]).variant_type != VariantType::SNP {
            continue;
        }

        let keys: Vec<&str> = f[8].split(':').collect();
        let kpos = |k: &str| keys.iter().position(|x| *x == k);
        let (Some(ad_k), Some(pl_k), Some(gt_k)) = (kpos("AD"), kpos("PL"), kpos("GT")) else {
            continue;
        };

        // Per-trio gl (normalized GL), dsg (ALT-allele count), ad.
        let mut gl = [0.0f64; 9];
        let mut dsg = [0i32; 3];
        let mut ad = [0i32; 6];
        let mut skip = false;
        for (i, &sidx) in trio.iter().enumerate() {
            let col = f[9 + sidx];
            let sub: Vec<&str> = col.split(':').collect();

            // PL: three values for a biallelic diploid site.
            let pl_field = sub.get(pl_k).copied().unwrap_or(".");
            let pl_v: Vec<&str> = pl_field.split(',').collect();
            if pl_v.len() < 3 {
                skip = true;
                break;
            }
            let mut isum = 0i64;
            let mut sum = 0.0f64;
            for j in 0..3 {
                let Ok(p) = pl_v[j].parse::<i32>() else {
                    skip = true;
                    break;
                };
                let g = 10f64.powf(-0.1 * p as f64);
                gl[3 * i + j] = g;
                sum += g;
                isum += p as i64;
            }
            if skip {
                break;
            }
            if isum == 0 {
                skip = true;
                break;
            }
            for j in 0..3 {
                gl[3 * i + j] /= sum;
            }

            // GT: require a fully-called diploid genotype.
            let gt_field = sub.get(gt_k).copied().unwrap_or(".");
            let alleles: Vec<&str> = gt_field.split(['/', '|']).collect();
            if alleles.len() != 2 {
                skip = true;
                break;
            }
            for a in &alleles {
                match a.parse::<i32>() {
                    Ok(v) => {
                        if v != 0 {
                            dsg[i] += 1;
                        }
                    }
                    Err(_) => {
                        skip = true;
                        break;
                    }
                }
            }
            if skip {
                break;
            }

            // AD: ref/alt depths.
            let ad_field = sub.get(ad_k).copied().unwrap_or(".");
            let ad_v: Vec<&str> = ad_field.split(',').collect();
            ad[2 * i] = ad_v.first().and_then(|s| s.parse().ok()).unwrap_or(0);
            ad[2 * i + 1] = ad_v.get(1).and_then(|s| s.parse().ok()).unwrap_or(0);
        }
        if skip {
            continue;
        }

        let (gl_p, gl_f, gl_m) = (&gl[0..3], &gl[3..6], &gl[6..9]);
        let (dsg_p, dsg_f, dsg_m) = (dsg[0], dsg[1], dsg[2]);

        match cnv_type {
            CnvType::Del => {
                if dsg_p != 0 && dsg_p != 2 {
                    continue; // proband not a hom
                }
                if dsg_f == dsg_m {
                    continue; // cannot distinguish between parents
                }
                if !greedy {
                    if dsg_f == 1 && dsg_p == dsg_m {
                        continue;
                    }
                    if dsg_m == 1 && dsg_p == dsg_f {
                        continue;
                    }
                }
                let p_mat = gl_p[0]
                    * (0.5 * gl_m[0] * gl_f[0]
                        + 2.0 / 3.0 * gl_m[0] * gl_f[1]
                        + gl_m[0] * gl_f[2]
                        + 1.0 / 3.0 * gl_m[1] * gl_f[0]
                        + 0.5 * gl_m[1] * gl_f[1]
                        + gl_m[1] * gl_f[2])
                    + gl_p[2]
                        * (0.5 * gl_m[2] * gl_f[2]
                            + 2.0 / 3.0 * gl_m[2] * gl_f[1]
                            + gl_m[2] * gl_f[0]
                            + 1.0 / 3.0 * gl_m[1] * gl_f[2]
                            + 0.5 * gl_m[1] * gl_f[1]
                            + gl_m[1] * gl_f[0]);
                let p_pat = gl_p[0]
                    * (0.5 * gl_m[0] * gl_f[0]
                        + 2.0 / 3.0 * gl_m[1] * gl_f[0]
                        + gl_m[2] * gl_f[0]
                        + 1.0 / 3.0 * gl_m[0] * gl_f[1]
                        + 0.5 * gl_m[1] * gl_f[1]
                        + gl_m[2] * gl_f[1])
                    + gl_p[2]
                        * (0.5 * gl_m[2] * gl_f[2]
                            + 2.0 / 3.0 * gl_m[1] * gl_f[2]
                            + gl_m[0] * gl_f[2]
                            + 1.0 / 3.0 * gl_m[2] * gl_f[1]
                            + 0.5 * gl_m[1] * gl_f[1]
                            + gl_m[0] * gl_f[1]);
                // NB: deliberate swap (upstream comment): the formulas give
                // the origin of the *observed* allele, the accumulators
                // track the origin of the *deleted* allele.
                pmat += p_pat.ln();
                ppat += p_mat.ln();
                ntest += 1;
            }
            CnvType::Dup => {
                let ad_p = &ad[0..2];
                let ad_f = &ad[2..4];
                let ad_m = &ad[4..6];
                if ad_p[0] == 0 || ad_p[1] == 0 {
                    continue;
                }
                if ad_p[0] == ad_p[1] {
                    continue;
                }
                if dsg_p != 1 {
                    continue;
                }
                if dsg_f == dsg_m {
                    continue;
                }
                if min_pbinom != 0.0 {
                    if dsg_f == 1
                        && ad_f[0] != 0
                        && ad_f[1] != 0
                        && calc_binom_two_sided(ad_f[0], ad_f[1], 0.5) < min_pbinom
                    {
                        continue;
                    }
                    if dsg_m == 1
                        && ad_m[0] != 0
                        && ad_m[1] != 0
                        && calc_binom_two_sided(ad_m[0], ad_m[1], 0.5) < min_pbinom
                    {
                        continue;
                    }
                }
                let prra = gl_p[1] * calc_binom_one_sided(ad_p[1], ad_p[0], 1.0 / 3.0, true);
                let praa = gl_p[1] * calc_binom_one_sided(ad_p[1], ad_p[0], 2.0 / 3.0, false);
                let p_pat = prra
                    * (gl_m[1] * gl_f[0]
                        + gl_m[2] * gl_f[0]
                        + 0.5 * gl_m[1] * gl_f[1]
                        + gl_m[2] * gl_f[1])
                    + praa
                        * (gl_m[1] * gl_f[2]
                            + gl_m[0] * gl_f[2]
                            + 0.5 * gl_m[1] * gl_f[1]
                            + gl_m[0] * gl_f[1]);
                let p_mat = prra
                    * (gl_m[0] * gl_f[1]
                        + gl_m[0] * gl_f[2]
                        + 0.5 * gl_m[1] * gl_f[1]
                        + gl_m[1] * gl_f[2])
                    + praa
                        * (gl_m[2] * gl_f[1]
                            + gl_m[2] * gl_f[0]
                            + 0.5 * gl_m[1] * gl_f[1]
                            + gl_m[1] * gl_f[0]);
                pmat += p_mat.ln();
                ppat += p_pat.ln();
                ntest += 1;
            }
        }
    }

    let qual = 4.3429 * (ppat - pmat).abs();
    let origin = if ppat > pmat {
        "paternal"
    } else if ppat < pmat {
        "maternal"
    } else {
        "uncertain"
    };

    let mut out = String::new();
    out.push_str(&format!("# bcftools +parental-origin{argv_tail}\n"));
    out.push_str("# [1]type\t[2]predicted_origin\t[3]quality\t[4]nmarkers\n");
    out.push_str(&format!(
        "{}\t{}\t{:.6}\t{}\n",
        cnv_type.label(),
        origin,
        qual,
        ntest
    ));
    Ok(out)
}

fn record_passes_filter(fields: &[&str], filter: FilterSpec<'_>) -> Result<bool, String> {
    let fields: Vec<String> = fields.iter().map(|field| (*field).to_owned()).collect();
    let matched =
        bcffilter::eval_expression_with(filter.expr, &EvalContext::new(), |name, sample_index| {
            if sample_index.is_some() {
                return None;
            }
            crate::commands::filter::record_lookup(name, &fields)
        })
        .map_err(|e| e.to_string())?
        .truthy();
    Ok(match filter.mode {
        FilterMode::Include => matched,
        FilterMode::Exclude => !matched,
    })
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
        ".bcftools-rs-parental-origin-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_forms() {
        let r = parse_region("20:100");
        assert_eq!((r.chrom.as_str(), r.beg, r.end), ("20", 100, 100));
        let r = parse_region("20:100-200");
        assert_eq!((r.chrom.as_str(), r.beg, r.end), ("20", 100, 200));
        let r = parse_region("X");
        assert_eq!((r.chrom.as_str(), r.beg), ("X", 1));
    }

    #[test]
    fn betai_symmetry() {
        // I_x(a,b) + I_{1-x}(b,a) == 1
        let v = kf_betai(2.0, 3.0, 0.4) + kf_betai(3.0, 2.0, 0.6);
        assert!((v - 1.0).abs() < 1e-12, "got {v}");
    }

    #[test]
    fn binom_edges() {
        assert_eq!(calc_binom_two_sided(0, 0, 0.5), -1.0);
        assert_eq!(calc_binom_two_sided(5, 5, 0.5), 1.0);
    }

    #[test]
    fn include_exclude_filters_records_before_likelihoods() {
        let vcf = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=20,length=81195210>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype Likelihoods\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tproband\tfather\tmother\n\
20\t100\t.\tA\tC\t.\t.\t.\tGT:PL:AD\t0/0:0,30,50:10,0\t0/0:0,30,50:10,0\t1/1:50,30,0:0,10\n\
20\t101\t.\tA\tC\t.\t.\t.\tGT:PL:AD\t0/0:0,30,50:10,0\t0/1:30,0,50:10,0\t1/1:50,30,0:0,10\n\
20\t102\t.\tA\tC\t.\t.\t.\tGT:PL:AD\t1/1:50,30,0:10,0\t0/1:30,0,50:10,0\t0/0:0,30,50:0,10\n";
        let pfm = "proband,father,mother";
        let nmarkers = |report: &str| {
            report
                .lines()
                .find(|l| !l.starts_with('#'))
                .and_then(|l| l.split('\t').nth(3))
                .unwrap()
                .to_owned()
        };

        let include = compute(
            vcf,
            "20:100-102",
            pfm,
            CnvType::Del,
            false,
            1e-2,
            "",
            Some(FilterSpec {
                mode: FilterMode::Include,
                expr: "POS=100",
            }),
        )
        .unwrap();
        assert_eq!(nmarkers(&include), "1");

        let exclude = compute(
            vcf,
            "20:100-102",
            pfm,
            CnvType::Del,
            false,
            1e-2,
            "",
            Some(FilterSpec {
                mode: FilterMode::Exclude,
                expr: "POS=100",
            }),
        )
        .unwrap();
        assert_eq!(nmarkers(&exclude), "2");
    }
}
