//! `bcftools +impute-info` (upstream `bcftools/plugins/impute-info.c`).
//!
//! Adds the IMPUTE2 `INFO/INFO` info-score metric, computed from
//! `FORMAT/GP` (biallelic-diploid genotype probabilities, 3 values per
//! sample), to each record. Sites without `FORMAT/GP`, or whose GP width
//! is not 3 (not biallelic diploid), are emitted unchanged. The
//! `##INFO=<ID=INFO,...>` header line is inserted after the last existing
//! `##INFO` line, mirroring HTSlib `bcf_hdr_append` grouping. A
//! total/added/skipped summary is written to stderr by `destroy()`, like
//! upstream.
//!
//! No upstream `test.pl` row exercises this plugin, so parity is held by
//! the synthetic integration suite in `tests/plugin_impute_info.rs` plus
//! the unit tests below; the info score itself is computed in `f64` and
//! stored as `f32` (matching upstream `bcf_update_info_float`), then
//! serialized via the shared HTSlib `kputd` formatter.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use super::prune::kputd;
use crate::vcf_compat::normalize_vcf_text;

const INFO_HEADER: &str = "##INFO=<ID=INFO,Number=1,Type=Float,Description=\"IMPUTE2 info score\">";

/// Result of an `impute-info` run: the rewritten VCF plus the stderr
/// lines upstream emits (first-occurrence warnings, then the `destroy()`
/// summary). The caller writes `vcf` to the output and `stderr` to stderr.
pub struct Output {
    pub vcf: String,
    pub stderr: String,
}

/// Reads the input VCF/BCF and returns the INFO-score-annotated VCF text
/// together with the upstream stderr summary/warnings.
pub fn run(input: &Path) -> io::Result<Output> {
    let text = read_vcf_text(input)?;
    Ok(annotate(&text))
}

fn annotate(text: &str) -> Output {
    let lines: Vec<&str> = text.lines().collect();
    let last_info = lines
        .iter()
        .rposition(|l| l.starts_with("##INFO="))
        .or_else(|| lines.iter().position(|l| l.starts_with("#CHROM")));

    let mut out = String::with_capacity(text.len() + 128);
    let mut nrec = 0u64; // info-added
    let mut nskip_gp = 0u64; // unchanged, no GP tag
    let mut nskip_dip = 0u64; // unchanged, not biallelic diploid

    for (idx, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            out.push_str(line);
            out.push('\n');
            if Some(idx) == last_info {
                if line.starts_with("##INFO=") {
                    out.push_str(INFO_HEADER);
                    out.push('\n');
                } else {
                    // No ##INFO lines: insert just before #CHROM.
                    out.truncate(out.len() - line.len() - 1);
                    out.push_str(INFO_HEADER);
                    out.push('\n');
                    out.push_str(line);
                    out.push('\n');
                }
            }
            continue;
        }
        if line.trim().is_empty() {
            out.push('\n');
            continue;
        }
        match annotate_record(line) {
            RecordOutcome::Added(rewritten) => {
                out.push_str(&rewritten);
                nrec += 1;
            }
            RecordOutcome::NoGp => {
                out.push_str(line);
                nskip_gp += 1;
            }
            RecordOutcome::NotBiallelicDiploid => {
                out.push_str(line);
                nskip_dip += 1;
            }
        }
        out.push('\n');
    }

    // Upstream prints each warning once, on first occurrence, then a final
    // tally in destroy(). Order relative to stdout is irrelevant (separate
    // streams); we emit warnings (if any) followed by the summary.
    let mut stderr = String::new();
    if nskip_gp > 0 {
        stderr.push_str("[impute-info.c] Warning: info tag not added to sites without GP tag\n");
    }
    if nskip_dip > 0 {
        stderr.push_str(
            "[impute-info.c] Warning: info tag not added to sites that are not biallelic diploid\n",
        );
    }
    stderr.push_str(&format!(
        "Lines total/info-added/unchanged-no-tag/unchanged-not-biallelic-diploid:\t{}/{}/{}/{}\n",
        nrec + nskip_gp + nskip_dip,
        nrec,
        nskip_gp,
        nskip_dip
    ));

    Output { vcf: out, stderr }
}

enum RecordOutcome {
    Added(String),
    NoGp,
    NotBiallelicDiploid,
}

fn annotate_record(line: &str) -> RecordOutcome {
    let mut fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 10 {
        // No FORMAT/sample columns => no GP available.
        return RecordOutcome::NoGp;
    }
    let Some(gp_idx) = fields[8].split(':').position(|k| k == "GP") else {
        return RecordOutcome::NoGp;
    };

    let n_sample = fields.len() - 9;
    // Per-sample GP value tokens; HTSlib pads to the max width across
    // samples with vector-end, so the effective width is that maximum.
    let mut per_sample: Vec<Vec<&str>> = Vec::with_capacity(n_sample);
    let mut width = 0usize;
    for sample in &fields[9..] {
        let toks: Vec<&str> = if *sample == "." {
            Vec::new()
        } else {
            match sample.split(':').nth(gp_idx) {
                Some(".") | None => Vec::new(),
                Some(gp) => gp.split(',').collect(),
            }
        };
        width = width.max(toks.len());
        per_sample.push(toks);
    }

    // Upstream: nret = total_values / n_sample, must equal 3.
    if width != 3 {
        return RecordOutcome::NotBiallelicDiploid;
    }

    let mut esum = 0.0f64;
    let mut e2sum = 0.0f64;
    let mut fsum = 0.0f64;
    for toks in &per_sample {
        let mut vals = [0.0f64; 3];
        for (j, v) in vals.iter_mut().enumerate() {
            // Missing or vector-end (short sample) breaks the inner loop,
            // leaving the remaining entries at 0 (upstream BRANCH macro).
            match toks.get(j) {
                None | Some(&".") => break,
                Some(&t) => match t.parse::<f32>() {
                    Ok(f) if !f.is_nan() => *v = f as f64,
                    _ => break,
                },
            }
        }
        let norm = vals[0] + vals[1] + vals[2];
        if norm != 0.0 {
            for v in &mut vals {
                *v /= norm;
            }
        }
        let e = vals[1] + 2.0 * vals[2];
        esum += e;
        e2sum += e * e;
        fsum += vals[1] + 4.0 * vals[2];
    }

    let nval = n_sample as f64;
    let theta = esum / (2.0 * nval);
    let info: f32 = if theta > 0.0 && theta < 1.0 {
        (1.0 - (fsum - e2sum) / (2.0 * nval * theta * (1.0 - theta))) as f32
    } else {
        1.0
    };

    let info_str = strip_info_key(fields[7]);
    let kv = format!("INFO={}", kputd(info as f64));
    let new_info = if info_str.is_empty() || info_str == "." {
        kv
    } else {
        format!("{info_str};{kv}")
    };
    fields[7] = new_info.as_str();
    RecordOutcome::Added(fields.join("\t"))
}

/// Removes any existing `INFO=...` key from the INFO column so the freshly
/// computed value replaces it (matching `bcf_update_info_float`).
fn strip_info_key(info: &str) -> String {
    if info == "." || info.is_empty() {
        return String::new();
    }
    info.split(';')
        .filter(|kv| {
            let key = kv.split_once('=').map(|(k, _)| k).unwrap_or(kv);
            key != "INFO"
        })
        .collect::<Vec<_>>()
        .join(";")
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
        ".bcftools-rs-impute-info-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn added(line: &str) -> String {
        match annotate_record(line) {
            RecordOutcome::Added(s) => s,
            RecordOutcome::NoGp => panic!("expected Added, got NoGp"),
            RecordOutcome::NotBiallelicDiploid => {
                panic!("expected Added, got NotBiallelicDiploid")
            }
        }
    }

    #[test]
    fn no_gp_format_key_is_unchanged() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:DP\t0/1:9";
        assert!(matches!(annotate_record(line), RecordOutcome::NoGp));
    }

    #[test]
    fn non_triplet_gp_is_not_biallelic_diploid() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/1:0.1,0.9";
        assert!(matches!(
            annotate_record(line),
            RecordOutcome::NotBiallelicDiploid
        ));
    }

    #[test]
    fn certain_homref_gives_info_one() {
        // All samples certain 0/0 => esum=0 => theta=0 => info=1.
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/0:1,0,0\t0/0:1,0,0";
        assert_eq!(
            added(line),
            "1\t10\t.\tC\tT\t.\tPASS\tINFO=1\tGT:GP\t0/0:1,0,0\t0/0:1,0,0"
        );
    }

    #[test]
    fn perfectly_called_dosages_give_info_one() {
        // Two samples, certain 0/0 and 1/1 => theta=0.5, fsum-e2sum=0 => I=1.
        let line = "1\t10\t.\tC\tT\t.\tPASS\tDP=5\tGT:GP\t0/0:1,0,0\t1/1:0,0,1";
        assert_eq!(
            added(line),
            "1\t10\t.\tC\tT\t.\tPASS\tDP=5;INFO=1\tGT:GP\t0/0:1,0,0\t1/1:0,0,1"
        );
    }

    #[test]
    fn unnormalized_gp_is_normalized() {
        // GP need not sum to 1; values are renormalized per sample.
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/1:0,2,0\t0/1:0,2,0";
        // Each sample: norm 2 -> [0,1,0]; e=1, theta=0.5, fsum=1,e2sum=1.
        // I = 1 - (2 - 2) / (2*2*0.25) = 1.
        assert_eq!(
            added(line),
            "1\t10\t.\tC\tT\t.\tPASS\tINFO=1\tGT:GP\t0/1:0,2,0\t0/1:0,2,0"
        );
    }

    #[test]
    fn existing_info_key_is_replaced() {
        let line = "1\t10\t.\tC\tT\t.\tPASS\tDP=3;INFO=9;AC=1\tGT:GP\t0/0:1,0,0";
        assert_eq!(
            added(line),
            "1\t10\t.\tC\tT\t.\tPASS\tDP=3;AC=1;INFO=1\tGT:GP\t0/0:1,0,0"
        );
    }

    #[test]
    fn uncertain_calls_reduce_info_below_one() {
        // Two samples GP 0.3,0.4,0.3: per sample e=1, e2=1, f=1.6.
        // theta = 2/(2*2) = 0.5 (in (0,1)); denom = 2*2*0.5*0.5 = 1.0.
        // I = 1 - (3.2 - 2.0)/1.0 = -0.2. The IMPUTE2 score is unbounded
        // below 1 and can be negative for poorly-imputed sites.
        let line = "1\t10\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/1:0.3,0.4,0.3\t0/1:0.3,0.4,0.3";
        let s = added(line);
        let info = s
            .split('\t')
            .nth(7)
            .and_then(|i| i.strip_prefix("INFO="))
            .and_then(|v| v.parse::<f64>().ok())
            .expect("INFO float");
        assert!(info < 1.0, "info={info} should be < 1");
        assert!(
            (info - (-0.2f32 as f64)).abs() < 1e-6,
            "info={info} should be ~-0.2"
        );
    }

    #[test]
    fn header_inserted_after_last_info_and_summary_counts() {
        let vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##contig=<ID=1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t1\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/0:1,0,0\n\
1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n\
1\t3\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/1:0.1,0.9\n";
        let o = annotate(vcf);
        let dp = o.vcf.find("##INFO=<ID=DP").unwrap();
        let inf = o.vcf.find("##INFO=<ID=INFO").unwrap();
        let chrom = o.vcf.find("#CHROM").unwrap();
        assert!(dp < inf && inf < chrom);
        assert!(o.vcf.contains("\tINFO=1\tGT:GP\t0/0:1,0,0\n"));
        // One added, one no-GP, one not-biallelic-diploid.
        assert!(o.stderr.contains(
            "Lines total/info-added/unchanged-no-tag/unchanged-not-biallelic-diploid:\t3/1/1/1\n"
        ));
        assert!(o.stderr.contains("without GP tag"));
        assert!(o.stderr.contains("not biallelic diploid"));
    }
}
