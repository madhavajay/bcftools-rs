//! `bcftools +guess-ploidy` (upstream `bcftools/plugins/guess-ploidy.c`).
//!
//! Determines sample sex from genotype likelihoods (PL/GL) or genotypes (GT)
//! in a haploid region (typically non-PAR chrX). For each SNP site it derives
//! per-genotype probabilities, an observed alternate-allele frequency, then
//! accumulates `log P(haploid)` / `log P(diploid)` per sample and reports the
//! per-sample mean log-likelihoods + predicted sex. All probability math is
//! in `f64` to match upstream's `double` precision. The default `-t PL`
//! auto-switches to `GL` then `GT` when the tag is absent from the header,
//! exactly like upstream `run()`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Tag {
    Gt,
    Pl,
    Gl,
}

impl Tag {
    pub fn parse(s: &str) -> Result<Tag, String> {
        if s.eq_ignore_ascii_case("GT") {
            Ok(Tag::Gt)
        } else if s.eq_ignore_ascii_case("PL") {
            Ok(Tag::Pl)
        } else if s.eq_ignore_ascii_case("GL") {
            Ok(Tag::Gl)
        } else {
            Err(format!(
                "The argument not recognised, expected --tag GT, PL or GL: {s}"
            ))
        }
    }
}

#[derive(Default, Clone)]
struct Counts {
    phap: f64,
    pdip: f64,
    ncount: i64,
}

#[derive(Debug, PartialEq, Eq)]
struct RegionSpec<'a> {
    chrom: &'a str,
    start: Option<u64>,
    end: Option<u64>,
}

impl<'a> RegionSpec<'a> {
    fn parse(raw: &'a str) -> RegionSpec<'a> {
        let Some((chrom, interval)) = raw.split_once(':') else {
            return RegionSpec {
                chrom: raw,
                start: None,
                end: None,
            };
        };
        let (start, end) = interval.split_once('-').unwrap_or((interval, interval));
        RegionSpec {
            chrom,
            start: start.parse::<u64>().ok(),
            end: end.parse::<u64>().ok(),
        }
    }

    fn matches(&self, chrom: &str, pos: &str) -> bool {
        if chrom != self.chrom {
            return false;
        }
        if self.start.is_none() && self.end.is_none() {
            return true;
        }
        let Ok(pos) = pos.parse::<u64>() else {
            return false;
        };
        self.start.is_none_or(|start| pos >= start) && self.end.is_none_or(|end| pos <= end)
    }
}

#[derive(Clone, Copy)]
pub struct Options<'a> {
    pub tag: Tag,
    pub region: Option<&'a str>,
    pub af_tag: Option<&'a str>,
    pub gt_err_prob: f64,
    pub af_dflt: f64,
    pub include_indels: bool,
    pub verbose: u32,
}

/// Reads the input and returns the guess-ploidy report (with the two
/// `bcftools`-tagged provenance lines the harness strips via
/// `grep -v bcftools`).
pub fn run(input: &Path, opts: Options<'_>) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    compute(&text, opts).map_err(io::Error::other)
}

fn compute(text: &str, opts: Options<'_>) -> Result<String, String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut tag = opts.tag;
    let region = opts.region.map(RegionSpec::parse);

    let mut has_pl = false;
    let mut has_gl = false;
    let mut has_gt = false;
    let mut has_af_tag = opts.af_tag.is_none();
    let mut samples: Vec<&str> = Vec::new();
    for l in &lines {
        if l.starts_with("##FORMAT=<ID=PL,") {
            has_pl = true;
        } else if l.starts_with("##FORMAT=<ID=GL,") {
            has_gl = true;
        } else if l.starts_with("##FORMAT=<ID=GT,") {
            has_gt = true;
        } else if let Some(tag) = opts.af_tag
            && info_header_has_id(l, tag)
        {
            has_af_tag = true;
        } else if l.starts_with("#CHROM") {
            samples = l.split('\t').skip(9).collect();
        }
    }
    if let Some(tag) = opts.af_tag
        && !has_af_tag
    {
        return Err(format!("No such INFO tag: {tag}"));
    }
    // PL -> GL -> GT auto-switch (upstream run()).
    if tag == Tag::Pl && !has_pl {
        tag = Tag::Gl;
    }
    if tag == Tag::Gl && !has_gl {
        tag = Tag::Gt;
    }
    let _ = has_gt;

    // pl2p[i] = 10^(-i/10), i in 0..256.
    let pl2p: Vec<f64> = (0..256).map(|i| 10f64.powf(-(i as f64) / 10.0)).collect();
    let pl2p_at = |v: i64| -> f64 {
        if !(0..256).contains(&v) {
            pl2p[255]
        } else {
            pl2p[v as usize]
        }
    };

    let nsample = samples.len();
    let mut counts = vec![Counts::default(); nsample];

    for l in &lines {
        if l.starts_with('#') || l.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = l.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        if let Some(r) = &region
            && !r.matches(f[0], f[1])
        {
            continue;
        }
        let reference = f[3];
        let alts: Vec<&str> = if f[4] == "." {
            Vec::new()
        } else {
            f[4].split(',').collect()
        };
        let n_allele = 1 + alts.len();
        if n_allele == 1 {
            continue;
        }
        if !opts.include_indels {
            let is_snp = alts.iter().any(|a| {
                classify_variant(reference, a)
                    .variant_type
                    .contains(VariantType::SNP)
            });
            if !is_snp {
                continue;
            }
        }

        let fmt = f[8];
        let sample_cols = &f[9..];
        // tmp[sample] = Some([pRR,pRA,pAA]) or None (non-informative/missing).
        let mut tmp: Vec<Option<[f64; 3]>> = vec![None; nsample];
        let mut freq = [0.0f64, 0.0f64];
        let ndip_gt = n_allele * (n_allele + 1) / 2;

        match tag {
            Tag::Gt => {
                let Some(slot) = fmt.split(':').position(|k| k == "GT") else {
                    continue;
                };
                for (i, s) in sample_cols.iter().enumerate() {
                    let gt = s.split(':').nth(slot).unwrap_or(".");
                    let toks: Vec<&str> = gt.split(['/', '|']).collect();
                    let a0 = toks.first().copied().unwrap_or(".");
                    if a0 == "." || a0.is_empty() {
                        continue; // missing -> tmp stays None
                    }
                    let p = opts.gt_err_prob;
                    let t = if toks.len() < 2 || toks[1] == "." {
                        // haploid
                        if a0 == "0" {
                            [1.0 - 2.0 * p, p, p]
                        } else {
                            [p, p, 1.0 - 2.0 * p]
                        }
                    } else {
                        let a1 = toks[1];
                        if a0 == "0" && a1 == "0" {
                            [1.0 - 2.0 * p, p, p]
                        } else if a0 == a1 {
                            [p, p, 1.0 - 2.0 * p]
                        } else {
                            [p, 1.0 - 2.0 * p, p]
                        }
                    };
                    tmp[i] = Some(t);
                    freq[0] += 2.0 * t[0] + t[1];
                    freq[1] += t[1] + 2.0 * t[2];
                }
            }
            Tag::Pl | Tag::Gl => {
                let key = if tag == Tag::Pl { "PL" } else { "GL" };
                let Some(slot) = fmt.split(':').position(|k| k == key) else {
                    continue;
                };
                // Per-sample value vectors; n = max length (htslib padding).
                let raw: Vec<Vec<&str>> = sample_cols
                    .iter()
                    .map(|s| {
                        s.split(':')
                            .nth(slot)
                            .unwrap_or(".")
                            .split(',')
                            .collect::<Vec<_>>()
                    })
                    .collect();
                let n = raw.iter().map(|v| v.len()).max().unwrap_or(0);
                let diploid = n == ndip_gt;
                let haploid = n == n_allele;
                if !diploid && !haploid {
                    continue;
                }
                for (i, vals) in raw.iter().enumerate() {
                    let parse = |idx: usize| -> Option<f64> {
                        vals.get(idx).and_then(|t| {
                            if *t == "." {
                                None
                            } else {
                                t.parse::<f64>().ok()
                            }
                        })
                    };
                    // Restrict to the first ALT: indices 0,1,2 of the
                    // diploid GL order; haploid uses 0,1.
                    let v0 = parse(0);
                    let v1 = parse(1);
                    let mut t: [f64; 3];
                    let mut vec_end = false;
                    if diploid {
                        let v2 = parse(2);
                        if v0.is_none() || v1.is_none() {
                            continue;
                        }
                        // ptr[2]==vector_end: this sample is haploid.
                        if vals.len() < 3 {
                            vec_end = true;
                        } else if v2.is_none() {
                            continue; // ptr[2] missing
                        }
                        let (a, b, c) = (v0.unwrap(), v1.unwrap(), v2.unwrap_or(0.0));
                        if !vec_end && a == b && a == c {
                            continue; // non-informative
                        }
                        if tag == Tag::Pl {
                            if vec_end {
                                t = [pl2p_at(a as i64), pl2p[255], pl2p_at(b as i64)];
                            } else {
                                t = [pl2p_at(a as i64), pl2p_at(b as i64), pl2p_at(c as i64)];
                            }
                        } else if vec_end {
                            t = [10f64.powf(a), 1e-26, 10f64.powf(b)];
                        } else {
                            t = [10f64.powf(a), 10f64.powf(b), 10f64.powf(c)];
                        }
                    } else {
                        // all-haploid record
                        if v0.is_none() || v1.is_none() {
                            continue;
                        }
                        let (a, b) = (v0.unwrap(), v1.unwrap());
                        if tag == Tag::Pl {
                            t = [pl2p_at(a as i64), pl2p[255], pl2p_at(b as i64)];
                        } else {
                            t = [10f64.powf(a), 1e-26, 10f64.powf(b)];
                        }
                        vec_end = true;
                    }
                    let sum = t[0] + t[1] + t[2];
                    if sum != 0.0 {
                        for x in t.iter_mut() {
                            *x /= sum;
                        }
                    }
                    tmp[i] = Some(t);
                    if vec_end {
                        freq[0] += t[0];
                        freq[1] += t[2];
                    } else {
                        freq[0] += 2.0 * t[0] + t[1];
                        freq[1] += t[1] + 2.0 * t[2];
                    }
                }
            }
        }

        if freq[0] == 0.0 && freq[1] == 0.0 {
            freq[0] = 1.0 - opts.af_dflt;
            freq[1] = opts.af_dflt;
        }
        if let Some(tag) = opts.af_tag
            && let Some(af) = parse_info_float_first(f[7], tag)
        {
            freq[0] = 1.0 - af;
            freq[1] = af;
        }
        let sum = freq[0] + freq[1];
        freq[0] /= sum;
        freq[1] /= sum;

        for (i, c) in counts.iter_mut().enumerate() {
            let Some(t) = tmp[i] else {
                continue;
            };
            let phap = freq[0] * t[0] + freq[1] * t[2];
            let pdip = freq[0] * freq[0] * t[0]
                + 2.0 * freq[0] * freq[1] * t[1]
                + freq[1] * freq[1] * t[2];
            c.phap += phap.ln();
            c.pdip += pdip.ln();
            c.ncount += 1;
        }
    }

    let mut out = String::new();
    if opts.verbose > 0 {
        out.push_str(
            "# This file was produced by: bcftools +guess-ploidy(bcftools-rs+htslib-rs)\n",
        );
        out.push_str("# The command line was:\tbcftools +guess-ploidy\n");
        out.push_str(
            "# [1]SEX\t[2]Sample\t[3]Predicted sex\t[4]log P(Haploid)/nSites\t\
[5]log P(Diploid)/nSites\t[6]nSites\t[7]Score: F < 0 < M ($4-$5)\n",
        );
    }
    for (i, s) in samples.iter().enumerate() {
        let c = &counts[i];
        let phap = if c.ncount != 0 {
            c.phap / c.ncount as f64
        } else {
            0.5
        };
        let pdip = if c.ncount != 0 {
            c.pdip / c.ncount as f64
        } else {
            0.5
        };
        let sex = if phap > pdip {
            'M'
        } else if phap < pdip {
            'F'
        } else {
            'U'
        };
        if opts.verbose > 0 {
            out.push_str(&format!(
                "SEX\t{s}\t{sex}\t{phap:.6}\t{pdip:.6}\t{}\t{:.6}\n",
                c.ncount,
                phap - pdip
            ));
        } else {
            out.push_str(&format!("{s}\t{sex}\n"));
        }
    }
    Ok(out)
}

fn info_header_has_id(line: &str, tag: &str) -> bool {
    line.strip_prefix("##INFO=<ID=")
        .and_then(|rest| rest.split([',', '>']).next())
        == Some(tag)
}

fn parse_info_float_first(info: &str, tag: &str) -> Option<f64> {
    if info == "." {
        return None;
    }
    for kv in info.split(';') {
        let Some((key, val)) = kv.split_once('=') else {
            continue;
        };
        if key != tag {
            continue;
        }
        let first = val.split(',').next()?;
        if first == "." {
            return None;
        }
        return first.parse::<f64>().ok();
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
        ".bcftools-rs-guess-ploidy-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts<'a>(region: Option<&'a str>, af_tag: Option<&'a str>) -> Options<'a> {
        Options {
            tag: Tag::Pl,
            region,
            af_tag,
            gt_err_prob: 1e-3,
            af_dflt: 0.5,
            include_indels: false,
            verbose: 1,
        }
    }

    #[test]
    fn tag_parse() {
        assert!(matches!(Tag::parse("PL"), Ok(Tag::Pl)));
        assert!(matches!(Tag::parse("gl"), Ok(Tag::Gl)));
        assert!(matches!(Tag::parse("GT"), Ok(Tag::Gt)));
        assert!(Tag::parse("XX").is_err());
    }

    #[test]
    fn region_spec_matches_chrom_and_interval() {
        assert!(RegionSpec::parse("X").matches("X", "1"));
        assert!(RegionSpec::parse("X:10-20").matches("X", "10"));
        assert!(RegionSpec::parse("X:10-20").matches("X", "20"));
        assert!(!RegionSpec::parse("X:10-20").matches("X", "9"));
        assert!(!RegionSpec::parse("X:10-20").matches("X", "21"));
        assert!(!RegionSpec::parse("X:10-20").matches("Y", "15"));
    }

    #[test]
    fn header_and_region_filter() {
        let vcf = "##fileformat=VCFv4.2\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
11\t100\t.\tA\tC\t.\t.\t.\tPL\t0,5,10\n\
X\t200\t.\tA\tC\t.\t.\t.\tPL\t0,5,10\n";
        let out = compute(vcf, opts(Some("X"), None)).unwrap();
        // Only the X site counts -> nSites == 1 for S1.
        let sex_line = out.lines().find(|l| l.starts_with("SEX\tS1")).unwrap();
        assert_eq!(sex_line.split('\t').nth(5).unwrap(), "1");
        assert!(out.contains("# [1]SEX\t[2]Sample"));
    }

    #[test]
    fn header_and_interval_region_filter() {
        let vcf = "##fileformat=VCFv4.2\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
X\t9\t.\tA\tC\t.\t.\t.\tPL\t0,5,10\n\
X\t10\t.\tA\tC\t.\t.\t.\tPL\t0,5,10\n\
X\t21\t.\tA\tC\t.\t.\t.\tPL\t0,5,10\n";
        let out = compute(vcf, opts(Some("X:10-20"), None)).unwrap();
        let sex_line = out.lines().find(|l| l.starts_with("SEX\tS1")).unwrap();
        assert_eq!(sex_line.split('\t').nth(5).unwrap(), "1");
    }

    #[test]
    fn af_tag_overrides_observed_frequency() {
        let vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n\
##INFO=<ID=OTHER,Number=A,Type=Float,Description=\"Other AF\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
X\t200\t.\tA\tC\t.\t.\tAF=0.01;OTHER=0.9\tPL\t0,10,100\n";
        let default = compute(vcf, opts(Some("X"), None)).unwrap();
        let tagged = compute(vcf, opts(Some("X"), Some("OTHER"))).unwrap();
        let default_score = default
            .lines()
            .find(|l| l.starts_with("SEX\tS1"))
            .unwrap()
            .split('\t')
            .nth(6)
            .unwrap()
            .to_owned();
        let tagged_score = tagged
            .lines()
            .find(|l| l.starts_with("SEX\tS1"))
            .unwrap()
            .split('\t')
            .nth(6)
            .unwrap()
            .to_owned();
        assert_ne!(default_score, tagged_score);
    }

    #[test]
    fn af_tag_must_exist_in_info_header() {
        let vcf = "##fileformat=VCFv4.2\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
X\t200\t.\tA\tC\t.\t.\tAF=0.01\tPL\t0,10,100\n";
        let err = compute(vcf, opts(Some("X"), Some("AF"))).unwrap_err();
        assert_eq!(err, "No such INFO tag: AF");
    }

    #[test]
    fn parse_info_float_first_value_or_none() {
        assert_eq!(
            parse_info_float_first("AN=2;AF=0.25,0.75", "AF"),
            Some(0.25)
        );
        assert_eq!(parse_info_float_first("AN=2;AF=.", "AF"), None);
        assert_eq!(parse_info_float_first("AN=2", "AF"), None);
        assert_eq!(parse_info_float_first(".", "AF"), None);
    }
}
