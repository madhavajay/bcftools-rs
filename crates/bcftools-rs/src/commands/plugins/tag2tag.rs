//! `bcftools +tag2tag` (upstream `bcftools/plugins/tag2tag.c`).
//!
//! Local slice: core GL/PL/GP/GT conversions
//! - `--gl-to-pl`: `PL = lround(-10 * GL)`, missing preserved.
//! - `--gl-to-gp`: `GP = 10^GL / sum(10^GL)`, missing preserved, with
//!   upstream `float` arithmetic and `%g`-style formatting.
//! - `--gp-to-gt`: hard-call `GT` from normalized `GP` with `-t`/`--threshold`
//!   (call iff max posterior >= 1 - threshold).
//! - `--LXX-to-XX`: expand localized `LAD`/`LPL` vectors to `AD`/`PL` through
//!   per-sample `LAA` maps.
//!
//! `-r`/`--replace` drops the source FORMAT tag (and its header line) and
//! adds the destination tag's header line as the last `##` line, mirroring
//! HTSlib `bcf_hdr_remove` + `bcf_hdr_append`. The localized `--LXX-to-XX`
//! family remains tracked in `TODO.md`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use super::prune::kputd;
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Conversion {
    GlToPl,
    GlToGp,
    GpToGt,
    LxxToXx,
}

const PL_HEADER: &str =
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">";
const GP_HEADER: &str =
    "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype probabilities\">";
const GT_HEADER: &str = "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">";
const AD_HEADER: &str = "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">";
const PASS_HEADER: &str = "##FILTER=<ID=PASS,Description=\"All filters passed\">";

#[derive(Clone, Copy, Debug, Default)]
pub struct LocalizedOptions {
    default_ad: Option<i32>,
    default_pl: Option<i32>,
    skip_nalt: usize,
}

impl LocalizedOptions {
    pub fn from_defaults(defaults: Option<&str>, skip_nalt: usize) -> io::Result<Self> {
        let mut opts = Self {
            default_ad: None,
            default_pl: None,
            skip_nalt,
        };
        let Some(defaults) = defaults else {
            return Ok(opts);
        };
        for part in defaults.split(',').filter(|p| !p.is_empty()) {
            let Some((tag, value)) = part.split_once(':') else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("Could not parse: --defaults {defaults}"),
                ));
            };
            let parsed = if value == "." {
                None
            } else {
                Some(value.parse::<i32>().map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("Could not parse: --defaults {defaults}"),
                    )
                })?)
            };
            match tag {
                "AD" => opts.default_ad = parsed,
                "PL" => opts.default_pl = parsed,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("Could not parse: --defaults {defaults}"),
                    ));
                }
            }
        }
        Ok(opts)
    }
}

struct Plan {
    src: &'static str,
    dst: &'static str,
    dst_header: &'static str,
}

fn plan(conv: Conversion) -> Plan {
    match conv {
        Conversion::GlToPl => Plan {
            src: "GL",
            dst: "PL",
            dst_header: PL_HEADER,
        },
        Conversion::GlToGp => Plan {
            src: "GL",
            dst: "GP",
            dst_header: GP_HEADER,
        },
        Conversion::GpToGt => Plan {
            src: "GP",
            dst: "GT",
            dst_header: GT_HEADER,
        },
        Conversion::LxxToXx => Plan {
            src: "LXX",
            dst: "XX",
            dst_header: "",
        },
    }
}

/// Reads the input VCF/BCF and returns the converted VCF text.
pub fn run(
    input: &Path,
    conv: Conversion,
    replace: bool,
    threshold: f64,
    localized: LocalizedOptions,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(convert_with_options(
        &text, conv, replace, threshold, localized,
    ))
}

#[cfg(test)]
fn convert(text: &str, conv: Conversion, replace: bool, threshold: f64) -> String {
    convert_with_options(text, conv, replace, threshold, LocalizedOptions::default())
}

fn convert_with_options(
    text: &str,
    conv: Conversion,
    replace: bool,
    threshold: f64,
    localized: LocalizedOptions,
) -> String {
    if conv == Conversion::LxxToXx {
        return convert_lxx_to_xx(text, replace, localized);
    }

    let p = plan(conv);
    let lines: Vec<&str> = text.lines().collect();
    let src_hdr_prefix = format!("##FORMAT=<ID={},", p.src);
    let chrom_idx = lines.iter().position(|l| l.starts_with("#CHROM"));

    let mut out = String::with_capacity(text.len() + 64);
    for (idx, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            if replace && line.starts_with(src_hdr_prefix.as_str()) {
                continue; // drop the source FORMAT header line
            }
            if Some(idx) == chrom_idx {
                out.push_str(p.dst_header);
                out.push('\n');
                out.push_str(line);
                out.push('\n');
                continue;
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            out.push('\n');
            continue;
        }
        out.push_str(&convert_record(line, conv, replace, threshold));
        out.push('\n');
    }
    out
}

fn convert_lxx_to_xx(text: &str, replace: bool, opts: LocalizedOptions) -> String {
    let lines: Vec<&str> = text.lines().collect();
    let chrom_idx = lines.iter().position(|l| l.starts_with("#CHROM"));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
    let remove_src_headers = replace && opts.skip_nalt == 0;

    let mut out = String::with_capacity(text.len() + 128);
    for (idx, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            if remove_src_headers
                && (line.starts_with("##FORMAT=<ID=LAD,")
                    || line.starts_with("##FORMAT=<ID=LPL,")
                    || line.starts_with("##FORMAT=<ID=LAA,"))
            {
                continue;
            }
            if Some(idx) == chrom_idx {
                out.push_str(PL_HEADER);
                out.push('\n');
                out.push_str(AD_HEADER);
                out.push('\n');
                out.push_str(line);
                out.push('\n');
                continue;
            }
            out.push_str(line);
            out.push('\n');
            if !has_pass && line.starts_with("##fileformat=") {
                out.push_str(PASS_HEADER);
                out.push('\n');
            }
            continue;
        }
        if line.trim().is_empty() {
            out.push('\n');
            continue;
        }
        out.push_str(&convert_lxx_record(line, replace, opts));
        out.push('\n');
    }
    out
}

fn convert_record(line: &str, conv: Conversion, replace: bool, threshold: f64) -> String {
    let p = plan(conv);
    let mut fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if fields.len() < 10 {
        return line.to_owned();
    }
    let fmt_keys: Vec<&str> = fields[8].split(':').collect();
    let Some(src_idx) = fmt_keys.iter().position(|k| *k == p.src) else {
        return line.to_owned();
    };

    let new_samples: Vec<String> = fields[9..]
        .iter()
        .map(|s| {
            let parts: Vec<&str> = s.split(':').collect();
            let src_val = parts.get(src_idx).copied().unwrap_or(".");
            let converted = match conv {
                Conversion::GlToPl => gl_to_pl(src_val),
                Conversion::GlToGp => gl_to_gp(src_val),
                Conversion::GpToGt => gp_to_gt(src_val, threshold),
                Conversion::LxxToXx => unreachable!("handled by convert_with_options"),
            };
            if replace {
                let mut out: Vec<String> = Vec::with_capacity(parts.len());
                for (i, val) in parts.iter().enumerate() {
                    if i == src_idx {
                        out.push(converted.clone());
                    } else {
                        out.push((*val).to_owned());
                    }
                }
                out.join(":")
            } else {
                let mut out: Vec<String> = parts.iter().map(|v| (*v).to_owned()).collect();
                out.push(converted);
                out.join(":")
            }
        })
        .collect();

    let new_format: String = if replace {
        fmt_keys
            .iter()
            .enumerate()
            .map(|(i, k)| if i == src_idx { p.dst } else { *k })
            .collect::<Vec<_>>()
            .join(":")
    } else {
        format!("{}:{}", fields[8], p.dst)
    };

    fields[8] = new_format;
    for (i, s) in new_samples.into_iter().enumerate() {
        fields[9 + i] = s;
    }
    fields.join("\t")
}

fn convert_lxx_record(line: &str, replace: bool, opts: LocalizedOptions) -> String {
    let mut fields: Vec<String> = line.split('\t').map(str::to_owned).collect();
    if fields.len() < 10 {
        return line.to_owned();
    }
    let nals = allele_count(&fields[3], &fields[4]);
    if opts.skip_nalt != 0 && nals > opts.skip_nalt {
        return line.to_owned();
    }

    let fmt_keys: Vec<&str> = fields[8].split(':').collect();
    let Some(laa_idx) = fmt_keys.iter().position(|k| *k == "LAA") else {
        return line.to_owned();
    };
    let lad_idx = fmt_keys.iter().position(|k| *k == "LAD");
    let lpl_idx = fmt_keys.iter().position(|k| *k == "LPL");
    if lad_idx.is_none() && lpl_idx.is_none() {
        return line.to_owned();
    }

    let new_samples: Vec<String> = fields[9..]
        .iter()
        .map(|sample| {
            let parts: Vec<&str> = sample.split(':').collect();
            let laa = parse_int_list(parts.get(laa_idx).copied().unwrap_or("."));
            let ad = lad_idx.and_then(|idx| {
                parts
                    .get(idx)
                    .map(|value| expand_lad(value, &laa, nals, opts.default_ad))
            });
            let pl = lpl_idx.and_then(|idx| {
                parts
                    .get(idx)
                    .map(|value| expand_lpl(value, &laa, nals, opts.default_pl))
            });

            if replace {
                let mut out = Vec::with_capacity(parts.len());
                for (i, val) in parts.iter().enumerate() {
                    if Some(i) == lad_idx || Some(i) == lpl_idx || i == laa_idx {
                        continue;
                    }
                    out.push((*val).to_owned());
                }
                if let Some(ad) = ad {
                    out.push(ad);
                }
                if let Some(pl) = pl {
                    out.push(pl);
                }
                out.join(":")
            } else {
                let mut out: Vec<String> = parts.iter().map(|v| (*v).to_owned()).collect();
                if let Some(ad) = ad {
                    out.push(ad);
                }
                if let Some(pl) = pl {
                    out.push(pl);
                }
                out.join(":")
            }
        })
        .collect();

    let new_format = if replace {
        let mut keys = Vec::with_capacity(fmt_keys.len());
        for (i, key) in fmt_keys.iter().enumerate() {
            if Some(i) == lad_idx || Some(i) == lpl_idx || i == laa_idx {
                continue;
            }
            keys.push((*key).to_owned());
        }
        if lad_idx.is_some() {
            keys.push("AD".to_owned());
        }
        if lpl_idx.is_some() {
            keys.push("PL".to_owned());
        }
        keys.join(":")
    } else {
        let mut keys = fields[8].clone();
        if lad_idx.is_some() {
            keys.push_str(":AD");
        }
        if lpl_idx.is_some() {
            keys.push_str(":PL");
        }
        keys
    };

    fields[8] = new_format;
    for (i, sample) in new_samples.into_iter().enumerate() {
        fields[9 + i] = sample;
    }
    fields.join("\t")
}

fn allele_count(_ref_allele: &str, alt: &str) -> usize {
    if alt == "." || alt.is_empty() {
        1
    } else {
        1 + alt.split(',').count()
    }
}

fn parse_int_list(value: &str) -> Vec<i32> {
    if value == "." || value.is_empty() {
        return Vec::new();
    }
    value
        .split(',')
        .filter_map(|v| v.parse::<i32>().ok())
        .collect()
}

fn format_int_list(values: &[Option<i32>]) -> String {
    values
        .iter()
        .map(|v| v.map_or_else(|| ".".to_owned(), |n| n.to_string()))
        .collect::<Vec<_>>()
        .join(",")
}

fn expand_lad(lad: &str, laa: &[i32], nals: usize, default: Option<i32>) -> String {
    let src = parse_int_list(lad);
    if src.is_empty() {
        return format_int_list(&vec![default; nals]);
    }
    let mut dst = vec![default; nals];
    dst[0] = Some(src[0]);
    for (src_idx, allele) in laa.iter().enumerate() {
        let src_pos = src_idx + 1;
        if src_pos >= src.len() {
            break;
        }
        if *allele >= 0 {
            let allele = *allele as usize;
            if allele < nals {
                dst[allele] = Some(src[src_pos]);
            }
        }
    }
    format_int_list(&dst)
}

fn expand_lpl(lpl: &str, laa: &[i32], nals: usize, default: Option<i32>) -> String {
    let src = parse_int_list(lpl);
    let ndst = nals * (nals + 1) / 2;
    let mut dst = vec![default; ndst];
    let mut localized = Vec::with_capacity(laa.len() + 1);
    localized.push(0);
    localized.extend_from_slice(laa);

    let mut src_idx = 0usize;
    for j in 0..localized.len() {
        let aj = localized[j];
        if aj < 0 || aj as usize >= nals {
            break;
        }
        for ak in localized.iter().take(j + 1) {
            if src_idx >= src.len() {
                break;
            }
            if *ak < 0 {
                src_idx += 1;
                continue;
            }
            let idx = aj as usize * (aj as usize + 1) / 2 + *ak as usize;
            if idx < dst.len() {
                dst[idx] = Some(src[src_idx]);
            }
            src_idx += 1;
        }
    }
    format_int_list(&dst)
}

fn gl_to_pl(gl: &str) -> String {
    if gl == "." {
        return ".".to_owned();
    }
    gl.split(',')
        .map(|v| {
            if v == "." {
                ".".to_owned()
            } else {
                match v.parse::<f64>() {
                    Ok(g) => ((-10.0 * g).round() as i64).to_string(),
                    Err(_) => ".".to_owned(),
                }
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn gl_to_gp(gl: &str) -> String {
    if gl == "." {
        return ".".to_owned();
    }

    let mut values: Vec<Option<f32>> = Vec::new();
    let mut sum = 0.0f32;
    for token in gl.split(',') {
        if token == "." {
            values.push(None);
            continue;
        }
        let Ok(g) = token.parse::<f32>() else {
            values.push(None);
            continue;
        };
        let p = 10.0f32.powf(g);
        sum += p;
        values.push(Some(p));
    }

    if sum > 0.0 {
        values
            .into_iter()
            .map(|v| match v {
                Some(p) => kputd((p / sum) as f64),
                None => ".".to_owned(),
            })
            .collect::<Vec<_>>()
            .join(",")
    } else {
        values
            .into_iter()
            .map(|v| match v {
                Some(p) => kputd(p as f64),
                None => ".".to_owned(),
            })
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn gp_to_gt(gp: &str, threshold: f64) -> String {
    if gp == "." {
        return "./.".to_owned();
    }
    let raw: Vec<f64> = gp
        .split(',')
        .map(|v| v.parse::<f64>().unwrap_or(f64::NAN))
        .collect();
    if raw.is_empty() || raw[0].is_nan() {
        return "./.".to_owned();
    }
    let sum: f64 = raw.iter().filter(|x| !x.is_nan()).sum();
    let norm: Vec<f64> = if sum > 0.0 {
        raw.iter().map(|x| x / sum).collect()
    } else {
        raw.clone()
    };
    let mut jmax = 0usize;
    for (j, v) in norm.iter().enumerate() {
        if v.is_nan() {
            break;
        }
        if *v > norm[jmax] {
            jmax = j;
        }
    }
    let n = norm.iter().take_while(|x| !x.is_nan()).count();
    let nals = n_alleles_for_diploid_count(n);
    let called = norm[jmax] >= 1.0 - threshold;
    if let Some(nals) = nals {
        // diploid
        if !called {
            return "./.".to_owned();
        }
        let (a, b) = gt2alleles(jmax);
        let _ = nals;
        format!("{a}/{b}")
    } else if n >= 1 {
        // treat as haploid: n == number of alleles
        if !called {
            ".".to_owned()
        } else {
            jmax.to_string()
        }
    } else {
        "./.".to_owned()
    }
}

/// If `count` equals `nals*(nals+1)/2` for some `nals >= 2`, returns `nals`
/// (diploid). `count == 1`/`2`/... that is also a valid triangular number for
/// small nals is resolved as diploid when `count >= 3` or `count == 1`.
fn n_alleles_for_diploid_count(count: usize) -> Option<usize> {
    // diploid genotype value count = nals*(nals+1)/2
    let mut nals = 1usize;
    loop {
        let tri = nals * (nals + 1) / 2;
        if tri == count && nals >= 2 {
            return Some(nals);
        }
        if tri > count {
            return None;
        }
        nals += 1;
        if nals > 64 {
            return None;
        }
    }
}

/// Genotype index -> (a, b) with a <= b, matching HTSlib `bcf_gt2alleles`.
fn gt2alleles(idx: usize) -> (usize, usize) {
    let mut k = 0usize;
    let mut b = 0usize;
    loop {
        for a in 0..=b {
            if k == idx {
                return (a, b);
            }
            k += 1;
        }
        b += 1;
    }
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
        ".bcftools-rs-tag2tag-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gl_to_pl_values() {
        assert_eq!(gl_to_pl("0"), "0");
        assert_eq!(gl_to_pl("0,0,0"), "0,0,0");
        assert_eq!(gl_to_pl("-25.5,0,-25.5"), "255,0,255");
        assert_eq!(gl_to_pl("."), ".");
        assert_eq!(gl_to_pl("-1,.,-2"), "10,.,20");
    }

    #[test]
    fn gl_to_gp_values() {
        assert_eq!(gl_to_gp("0"), "1");
        assert_eq!(gl_to_gp("0,0,0"), "0.333333,0.333333,0.333333");
        assert_eq!(gl_to_gp("-25.5,0,-25.5"), "3.16228e-26,1,3.16228e-26");
        assert_eq!(gl_to_gp(".,.,."), ".,.,.");
        assert_eq!(gl_to_gp("."), ".");
    }

    #[test]
    fn gp_to_gt_threshold() {
        assert_eq!(gp_to_gt("0.962,0.038,0", 0.2), "0/0");
        assert_eq!(gp_to_gt("0,1,0", 0.2), "0/1");
        assert_eq!(gp_to_gt("1,0,0", 0.2), "0/0");
        assert_eq!(gp_to_gt("0,0.443,0.557", 0.2), "./."); // max 0.557 < 0.8
        assert_eq!(gp_to_gt(".", 0.2), "./.");
    }

    #[test]
    fn localized_expansion_uses_laa_mapping() {
        assert_eq!(expand_lad("13,12", &[2], 3, None), "13,.,12");
        assert_eq!(expand_lpl("82,0,55", &[2], 3, None), "82,.,.,0,.,55");
        assert_eq!(
            expand_lpl("305,339,82,220,0,61", &[5, 2], 8, None),
            "305,.,.,220,.,61,.,.,0,.,.,.,.,.,.,339,.,.,.,.,82,.,.,.,.,.,.,.,.,.,.,.,.,.,.,."
        );
    }

    #[test]
    fn localized_defaults_fill_missing_slots() {
        assert_eq!(expand_lad("31", &[], 3, Some(0)), "31,0,0");
        assert_eq!(expand_lpl("0", &[], 3, Some(255)), "0,255,255,255,255,255");
    }

    #[test]
    fn gt2alleles_layout() {
        assert_eq!(gt2alleles(0), (0, 0));
        assert_eq!(gt2alleles(1), (0, 1));
        assert_eq!(gt2alleles(2), (1, 1));
        assert_eq!(gt2alleles(3), (0, 2));
        assert_eq!(gt2alleles(4), (1, 2));
        assert_eq!(gt2alleles(5), (2, 2));
    }

    #[test]
    fn replace_rewrites_format_and_header() {
        let vcf = "##fileformat=VCFv4.1\n\
##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"x\">\n\
##FILTER=<ID=q,Description=\"y\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
11\t5\t.\tC\tT\t.\tPASS\t.\tGL\t-25.5,0,-25.5\n";
        let out = convert(vcf, Conversion::GlToPl, true, 0.1);
        assert!(!out.contains("##FORMAT=<ID=GL,"));
        // dst header is the last ## line, just before #CHROM
        let pl_pos = out.find(PL_HEADER).unwrap();
        let chrom_pos = out.find("#CHROM").unwrap();
        let filter_pos = out.find("##FILTER=<ID=q,").unwrap();
        assert!(filter_pos < pl_pos && pl_pos < chrom_pos);
        assert!(out.contains("\tPL\t255,0,255\n"));
    }
}
