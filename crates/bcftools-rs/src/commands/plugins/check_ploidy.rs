//! `bcftools +check-ploidy` (upstream `bcftools/plugins/check-ploidy.c`).
//!
//! Reports, per sample, contiguous regions of constant genotype ploidy as
//! `Sample  Chrom  Start  End  Ploidy`. By default genotypes containing any
//! missing allele are ignored (they neither establish nor extend a region);
//! `-m`/`--use-missing` counts missing-allele slots toward ploidy instead.
//!
//! Faithful to upstream's flush timing: a chromosome change flushes open
//! regions under the *previous* chromosome name; a ploidy change within a
//! chromosome flushes under the current one.

use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const HEADER: &str = "# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n";

#[derive(Clone)]
struct Region {
    ploidy: usize, // 0 = no open region
    beg: i64,
    end: i64,
}

/// Reads the input VCF/BCF and returns the upstream-shaped ploidy report.
/// `use_missing` corresponds to upstream `-m`/`--use-missing`.
pub fn run(input: &Path, use_missing: bool) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(report(&text, use_missing))
}

fn report(text: &str, use_missing: bool) -> String {
    let ignore_missing = !use_missing;
    let mut out = String::from(HEADER);

    let mut samples: Vec<&str> = Vec::new();
    let mut state: Vec<Region> = Vec::new();
    let mut cur_chrom: Option<String> = None;

    for line in text.lines() {
        if line.starts_with("#CHROM") {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() > 9 {
                samples = cols[9..].to_vec();
            }
            state = samples
                .iter()
                .map(|_| Region {
                    ploidy: 0,
                    beg: 0,
                    end: 0,
                })
                .collect();
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if samples.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }
        let Some(gt_idx) = fields[8].split(':').position(|k| k == "GT") else {
            continue; // no GT tag -> record skipped entirely
        };
        let chrom = fields[0];
        let Ok(pos) = fields[1].parse::<i64>() else {
            continue;
        };

        if let Some(prev) = &cur_chrom
            && prev != chrom
        {
            for (i, st) in state.iter_mut().enumerate() {
                if st.ploidy != 0 {
                    let _ = writeln!(
                        out,
                        "{}\t{}\t{}\t{}\t{}",
                        samples[i], prev, st.beg, st.end, st.ploidy
                    );
                }
                st.ploidy = 0;
            }
        }
        cur_chrom = Some(chrom.to_owned());

        for (i, sample) in fields[9..].iter().enumerate() {
            if i >= state.len() {
                break;
            }
            let Some(gt) = sample.split(':').nth(gt_idx) else {
                continue;
            };
            let mut nal = 0usize;
            let mut missing = false;
            for tok in gt.split(['/', '|']) {
                if tok == "." && ignore_missing {
                    missing = true;
                    break;
                }
                nal += 1;
            }
            if nal == 0 || missing {
                continue;
            }
            let st = &mut state[i];
            if st.ploidy == nal {
                st.end = pos;
                continue;
            }
            if st.ploidy != 0 {
                let _ = writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    samples[i], chrom, st.beg, st.end, st.ploidy
                );
            }
            st.ploidy = nal;
            st.beg = pos;
            st.end = pos;
        }
    }

    if let Some(chrom) = &cur_chrom {
        for (i, st) in state.iter().enumerate() {
            if st.ploidy != 0 {
                let _ = writeln!(
                    out,
                    "{}\t{}\t{}\t{}\t{}",
                    samples[i], chrom, st.beg, st.end, st.ploidy
                );
            }
        }
    }

    out
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
        ".bcftools-rs-check-ploidy-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIX1: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
22\t1\t1\tC\tG\t.\tPASS\t.\tGT\t1/1\n\
X\t1\t2\tC\tG\t.\tPASS\t.\tGT\t1\n";

    const FIX2: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
20\t1\t1\tC\tG\t.\tPASS\t.\tGT\t0/0\n\
21\t1\t2\tC\tG\t.\tPASS\t.\tGT\t./0\n\
22\t1\t3\tC\tG\t.\tPASS\t.\tGT\t./.\n\
X\t1\t4\tC\tG\t.\tPASS\t.\tGT\t.\n";

    #[test]
    fn default_reports_fully_called_runs() {
        assert_eq!(
            report(FIX1, false),
            "# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n\
S1\t22\t1\t1\t2\n\
S1\tX\t1\t1\t1\n"
        );
    }

    #[test]
    fn default_skips_any_missing_allele() {
        assert_eq!(
            report(FIX2, false),
            "# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n\
S1\t20\t1\t1\t2\n"
        );
    }

    #[test]
    fn use_missing_counts_missing_slots() {
        assert_eq!(
            report(FIX2, true),
            "# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n\
S1\t20\t1\t1\t2\n\
S1\t21\t1\t1\t2\n\
S1\t22\t1\t1\t2\n\
S1\tX\t1\t1\t1\n"
        );
    }

    #[test]
    fn extends_region_across_constant_ploidy() {
        let vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
1\t10\t.\tC\tG\t.\t.\t.\tGT\t0/1\n\
1\t20\t.\tC\tG\t.\t.\t.\tGT\t1/1\n\
1\t30\t.\tC\tG\t.\t.\t.\tGT\t1\n";
        assert_eq!(
            report(vcf, false),
            "# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n\
S1\t1\t10\t20\t2\n\
S1\t1\t30\t30\t1\n"
        );
    }
}
