//! `bcftools +contrast` (upstream `bcftools/plugins/contrast.c`).
//!
//! Basic association test + novel allele/genotype detection between a
//! control (`-0`) and a case (`-1`) sample group. Adds INFO annotations:
//! `PASSOC` (Fisher's exact test, REF vs non-REF), `FASSOC` (non-REF
//! proportion in controls/cases), `NASSOC` (the 4 allele counts),
//! `NOVELAL`/`NOVELGT` (case samples with an allele/genotype not seen in
//! controls). Fisher's exact test routes through
//! `htslib_rs::math::kt_fisher_exact`; float INFO uses the shared HTSlib
//! `kputd` formatter. `-i`/`-e` filtering and `-f` rare-allele enrichment
//! need infrastructure tracked in `TODO.md`.

use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::math::kt_fisher_exact;

use super::prune::kputd;
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, Default)]
pub struct Annots {
    pub passoc: bool,
    pub fassoc: bool,
    pub nassoc: bool,
    pub novelal: bool,
    pub novelgt: bool,
}

impl Annots {
    pub fn parse(list: &str) -> Result<Annots, String> {
        let mut a = Annots::default();
        for t in list.split(',') {
            match t.to_ascii_uppercase().as_str() {
                "PASSOC" => a.passoc = true,
                "FASSOC" => a.fassoc = true,
                "NASSOC" => a.nassoc = true,
                "NOVELAL" => a.novelal = true,
                "NOVELGT" => a.novelgt = true,
                other => return Err(format!("Unknown annotation: {other}")),
            }
        }
        Ok(a)
    }
}

/// Reads the input and returns the contrast-annotated VCF text.
pub fn run(
    input: &Path,
    annots: Annots,
    control: &str,
    case: &str,
    force_samples: bool,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(&text, annots, control, case, force_samples).map_err(io::Error::other)
}

/// Resolve a `-0`/`-1` argument: a comma list of sample names, or (only if
/// not all tokens are VCF samples) a file with one sample per line.
/// `--force-samples` drops names absent from the VCF instead of erroring.
fn resolve_samples(arg: &str, sample_idx: &[&str], force: bool) -> Result<Vec<usize>, String> {
    let idx_of = |n: &str| sample_idx.iter().position(|s| *s == n);
    let list_toks: Vec<&str> = arg.split(',').collect();
    let all_samples = list_toks.iter().all(|t| idx_of(t).is_some());
    let toks: Vec<String> = if all_samples {
        list_toks.iter().map(|s| s.to_string()).collect()
    } else if Path::new(arg).is_file() {
        fs::read_to_string(arg)
            .map_err(|e| e.to_string())?
            .lines()
            .map(|l| l.trim().to_string())
            .filter(|l| !l.is_empty())
            .collect()
    } else {
        list_toks.iter().map(|s| s.to_string()).collect()
    };
    let mut out = Vec::new();
    for t in &toks {
        match idx_of(t) {
            Some(i) => out.push(i),
            None => {
                if !force {
                    return Err(format!("The sample \"{t}\" is not present in the VCF"));
                }
            }
        }
    }
    Ok(out)
}

fn parse_alleles(gt: &str) -> Vec<Option<i32>> {
    // None = missing ('.'); the upstream loop *continues* on missing and
    // *breaks* on vector_end (a shorter token list).
    gt.split(['/', '|'])
        .map(|t| {
            if t == "." || t.is_empty() {
                None
            } else {
                t.parse::<i32>().ok()
            }
        })
        .collect()
}

fn compute(
    text: &str,
    annots: Annots,
    control: &str,
    case: &str,
    force_samples: bool,
) -> Result<String, String> {
    let lines: Vec<&str> = text.lines().collect();
    let sample_idx: Vec<&str> = lines
        .iter()
        .find(|l| l.starts_with("#CHROM"))
        .map(|l| l.split('\t').skip(9).collect())
        .unwrap_or_default();

    let control_smpl = resolve_samples(control, &sample_idx, force_samples)?;
    let case_smpl = resolve_samples(case, &sample_idx, force_samples)?;

    let mut out = String::with_capacity(text.len() + 512);
    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
    for (idx, line) in lines.iter().enumerate() {
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("#CHROM") {
            if annots.passoc {
                out.push_str("##INFO=<ID=PASSOC,Number=1,Type=Float,Description=\"Fisher's exact test probability of genotypic association (REF vs non-REF allele)\">\n");
            }
            if annots.fassoc {
                out.push_str("##INFO=<ID=FASSOC,Number=2,Type=Float,Description=\"Proportion of non-REF allele in controls and cases\">\n");
            }
            if annots.nassoc {
                out.push_str("##INFO=<ID=NASSOC,Number=4,Type=Integer,Description=\"Number of control-ref, control-alt, case-ref and case-alt alleles\">\n");
            }
            if annots.novelal {
                out.push_str("##INFO=<ID=NOVELAL,Number=.,Type=String,Description=\"List of samples with novel alleles. Note that samples listed here are not listed in NOVELGT again.\">\n");
            }
            if annots.novelgt {
                out.push_str("##INFO=<ID=NOVELGT,Number=.,Type=String,Description=\"List of samples with novel genotypes\">\n");
            }
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

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        let gt_slot = f[8].split(':').position(|k| k == "GT");
        let sample_cols = &f[9..];
        let gt_of = |s: usize| -> Vec<Option<i32>> {
            let raw = match gt_slot {
                Some(g) => sample_cols
                    .get(s)
                    .and_then(|c| c.split(':').nth(g))
                    .unwrap_or("."),
                None => ".",
            };
            parse_alleles(raw)
        };

        let mut control_als: u32 = 0;
        let mut nals = [0i32; 4]; // ctrl-ref, ctrl-alt, case-ref, case-alt
        let mut control_gts: BTreeSet<u32> = BTreeSet::new();
        for &s in &control_smpl {
            let mut gt: u32 = 0;
            for al in gt_of(s) {
                let Some(ial) = al else {
                    continue; // missing
                };
                let ial = ial as u32;
                control_als |= 1 << ial;
                gt |= 1 << ial;
                if ial != 0 {
                    nals[1] += 1;
                } else {
                    nals[0] += 1;
                }
            }
            if annots.novelgt {
                control_gts.insert(gt);
            }
        }

        let mut skipped = false;
        if control_als == 0 && !control_smpl.is_empty() {
            skipped = true;
        }

        let mut novelal: Vec<&str> = Vec::new();
        let mut novelgt: Vec<&str> = Vec::new();
        let mut has_gt = false;
        if !skipped {
            for &s in &case_smpl {
                let mut case_al = false;
                let mut gt: u32 = 0;
                for al in gt_of(s) {
                    let Some(ial) = al else {
                        continue;
                    };
                    let ial = ial as u32;
                    if control_als & (1 << ial) == 0 {
                        case_al = true;
                    }
                    gt |= 1 << ial;
                    if ial != 0 {
                        nals[3] += 1;
                    } else {
                        nals[2] += 1;
                    }
                }
                if gt == 0 {
                    continue;
                }
                has_gt = true;
                let name = sample_idx[s];
                if case_al && annots.novelal {
                    novelal.push(name);
                } else if annots.novelgt && !control_gts.contains(&gt) {
                    novelgt.push(name);
                }
            }
            if !has_gt && !case_smpl.is_empty() {
                skipped = true;
            }
        }

        // Build the annotated INFO. Skipped records are written verbatim.
        if !skipped {
            let mut info = if f[7] == "." || f[7].is_empty() {
                String::new()
            } else {
                f[7].to_string()
            };
            let push = |info: &mut String, k: &str, v: &str| {
                if !info.is_empty() {
                    info.push(';');
                }
                info.push_str(k);
                info.push('=');
                info.push_str(v);
            };
            if annots.passoc && !control_smpl.is_empty() && !case_smpl.is_empty() {
                let p = kt_fisher_exact(nals[0], nals[1], nals[2], nals[3]).two_tail;
                push(&mut info, "PASSOC", &kputd(p as f32 as f64));
            }
            if annots.fassoc && !control_smpl.is_empty() && !case_smpl.is_empty() {
                let v0 = if nals[0] + nals[1] != 0 {
                    kputd((nals[1] as f32 / (nals[0] + nals[1]) as f32) as f64)
                } else {
                    ".".to_owned()
                };
                let v1 = if nals[2] + nals[3] != 0 {
                    kputd((nals[3] as f32 / (nals[2] + nals[3]) as f32) as f64)
                } else {
                    ".".to_owned()
                };
                push(&mut info, "FASSOC", &format!("{v0},{v1}"));
            }
            if annots.nassoc {
                push(
                    &mut info,
                    "NASSOC",
                    &format!("{},{},{},{}", nals[0], nals[1], nals[2], nals[3]),
                );
            }
            if !novelal.is_empty() {
                push(&mut info, "NOVELAL", &novelal.join(","));
            }
            if !novelgt.is_empty() {
                push(&mut info, "NOVELGT", &novelgt.join(","));
            }
            let owned = if info.is_empty() {
                ".".to_owned()
            } else {
                info
            };
            f[7] = "";
            let mut joined: Vec<String> = f.iter().map(|s| s.to_string()).collect();
            joined[7] = owned;
            out.push_str(&joined.join("\t"));
            out.push('\n');
        } else {
            out.push_str(line);
            out.push('\n');
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
        ".bcftools-rs-contrast-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n\
1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/0\n\
1\t101\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n\
1\t102\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\n\
1\t103\t.\tA\tG\t.\t.\t.\tGT\t1/1\t1/1\t0/1\n";

    fn body(out: &str) -> Vec<&str> {
        out.lines().filter(|l| !l.starts_with('#')).collect()
    }

    #[test]
    fn matches_upstream_contrast_out() {
        let a = Annots::parse("PASSOC,FASSOC,NOVELAL,NOVELGT").unwrap();
        let out = compute(VCF, a, "a,b", "c", false).unwrap();
        let d = body(&out);
        assert_eq!(
            d[0],
            "1\t100\t.\tA\tG\t.\t.\tPASSOC=1;FASSOC=0,0\tGT\t0/0\t0/0\t0/0"
        );
        assert_eq!(
            d[1],
            "1\t101\t.\tA\tG\t.\t.\tPASSOC=0.333333;FASSOC=0,0.5;NOVELAL=c\tGT\t0/0\t0/0\t0/1"
        );
        assert_eq!(
            d[2],
            "1\t102\t.\tA\tG\t.\t.\tPASSOC=0.4;FASSOC=0.25,1;NOVELGT=c\tGT\t0/0\t0/1\t1/1"
        );
        assert_eq!(
            d[3],
            "1\t103\t.\tA\tG\t.\t.\tPASSOC=0.333333;FASSOC=1,0.5;NOVELAL=c\tGT\t1/1\t1/1\t0/1"
        );
    }

    #[test]
    fn nassoc_with_force_samples_missing_case() {
        let a = Annots::parse("NASSOC").unwrap();
        let out = compute(VCF, a, "a,b,c", "d", true).unwrap();
        let d = body(&out);
        assert_eq!(
            d[0],
            "1\t100\t.\tA\tG\t.\t.\tNASSOC=6,0,0,0\tGT\t0/0\t0/0\t0/0"
        );
        assert_eq!(
            d[3],
            "1\t103\t.\tA\tG\t.\t.\tNASSOC=1,5,0,0\tGT\t1/1\t1/1\t0/1"
        );
    }
}
