//! `bcftools +fill-tags` (upstream `bcftools/plugins/fill-tags.c`).
//!
//! First slice: the genotype-derived INFO count tags
//! `AN`/`AC`/`AC_Hom`/`AC_Het`/`AC_Hemi`/`AF`/`MAF`/`NS`, the `-t LIST`
//! tag selection, and `-S`/`--samples-file` population grouping
//! (per-group `<TAG>_<pop>` plus the always-present global tag). The
//! counting mirrors upstream `process_fmt`'s per-sample classification
//! (hom/het/hemi/half-missing), and tags are written in the upstream
//! fixed order (NS, AN, AF, MAF, AC, AC_Het, AC_Hom, AC_Hemi),
//! replacing an existing INFO key in place or appending otherwise.
//!
//! Deferred (tracked in TODO.md): `HWE`/`ExcHet` (kf functions),
//! `END`/`TYPE`/`F_MISSING`/`VAF`/`VAF1`, the `TAG:Num=EXPR` function
//! engine (`sum`/`ssum`/`fisher`/`binom`/`F_PASS`/`N_PASS`/…), the
//! `all` tag set, and `-d`/`--drop-missing`.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Tag {
    Ns,
    An,
    Af,
    Maf,
    Ac,
    AcHet,
    AcHom,
    AcHemi,
}

/// `-t` tags in the fixed upstream `process_fmt` write order.
const WRITE_ORDER: &[Tag] = &[
    Tag::Ns,
    Tag::An,
    Tag::Af,
    Tag::Maf,
    Tag::Ac,
    Tag::AcHet,
    Tag::AcHom,
    Tag::AcHemi,
];

fn parse_tag(name: &str) -> Result<Tag, String> {
    let n = name.strip_prefix("INFO/").unwrap_or(name);
    match n.to_ascii_uppercase().as_str() {
        "NS" => Ok(Tag::Ns),
        "AN" => Ok(Tag::An),
        "AF" => Ok(Tag::Af),
        "MAF" => Ok(Tag::Maf),
        "AC" => Ok(Tag::Ac),
        "AC_HET" => Ok(Tag::AcHet),
        "AC_HOM" => Ok(Tag::AcHom),
        "AC_HEMI" => Ok(Tag::AcHemi),
        _ => Err(format!(
            "the tag \"{name}\" is not yet ported (this slice supports \
             AN,AC,AC_Hom,AC_Het,AC_Hemi,AF,MAF,NS)"
        )),
    }
}

impl Tag {
    /// Base ID and the `hdr_append` `##INFO` template (with the `{S}`
    /// suffix and `{IN}` (` in <pop>`) placeholders) used when the tag's
    /// ID is not already declared in the header.
    fn header(self) -> (&'static str, &'static str) {
        match self {
            Tag::An => (
                "AN",
                "##INFO=<ID=AN{S},Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes{IN}\">",
            ),
            Tag::Ac => (
                "AC",
                "##INFO=<ID=AC{S},Number=A,Type=Integer,Description=\"Allele count in genotypes{IN}\">",
            ),
            Tag::Ns => (
                "NS",
                "##INFO=<ID=NS{S},Number=1,Type=Integer,Description=\"Number of samples with data{IN}\">",
            ),
            Tag::AcHom => (
                "AC_Hom",
                "##INFO=<ID=AC_Hom{S},Number=A,Type=Integer,Description=\"Allele counts in homozygous genotypes{IN}\">",
            ),
            Tag::AcHet => (
                "AC_Het",
                "##INFO=<ID=AC_Het{S},Number=A,Type=Integer,Description=\"Allele counts in heterozygous genotypes{IN}\">",
            ),
            Tag::AcHemi => (
                "AC_Hemi",
                "##INFO=<ID=AC_Hemi{S},Number=A,Type=Integer,Description=\"Allele counts in hemizygous genotypes{IN}\">",
            ),
            Tag::Af => (
                "AF",
                "##INFO=<ID=AF{S},Number=A,Type=Float,Description=\"Allele frequency{IN}\">",
            ),
            Tag::Maf => (
                "MAF",
                "##INFO=<ID=MAF{S},Number=1,Type=Float,Description=\"Frequency of the second most common allele{IN}\">",
            ),
        }
    }
    fn base_id(self) -> &'static str {
        self.header().0
    }
}

/// `hdr_append` order (upstream lines 590-598).
const HDR_ORDER: &[Tag] = &[
    Tag::An,
    Tag::Ac,
    Tag::Ns,
    Tag::AcHom,
    Tag::AcHet,
    Tag::AcHemi,
    Tag::Af,
    Tag::Maf,
];

#[derive(Default, Clone)]
struct Counts {
    nhet: i64,
    nhom: i64,
    nhemi: i64,
    nac: i64,
}

struct Pop {
    /// `""` for the global pop; otherwise the population name.
    name: String,
    /// `""` (global) or `_<name>`.
    suffix: String,
}

pub struct Options<'a> {
    /// `-t` comma-separated tag list (required in this slice).
    pub tags: &'a str,
    /// `-S`/`--samples-file` path (`sample<ws>pop1,pop2` per line).
    pub samples_file: Option<&'a Path>,
}

pub fn run(input: &Path, opts: Options<'_>) -> io::Result<String> {
    let mut want: Vec<Tag> = Vec::new();
    for t in opts.tags.split(',').filter(|t| !t.is_empty()) {
        let tag = parse_tag(t)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, format!("fill-tags: {e}")))?;
        if !want.contains(&tag) {
            want.push(tag);
        }
    }
    if want.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "fill-tags: this slice requires -t with a supported tag list",
        ));
    }

    // `-S` population map: name -> comma list; preserve first-seen order.
    let mut pop_order: Vec<String> = Vec::new();
    let mut sample_pops: HashMap<String, Vec<String>> = HashMap::new();
    if let Some(sf) = opts.samples_file {
        let body = fs::read_to_string(sf)?;
        for line in body.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut it = line.split_whitespace();
            let (Some(sample), Some(pops)) = (it.next(), it.next()) else {
                continue;
            };
            let list: Vec<String> = pops
                .split(',')
                .filter(|p| !p.is_empty())
                .map(|p| {
                    let p = p.to_owned();
                    if !pop_order.contains(&p) {
                        pop_order.push(p.clone());
                    }
                    p
                })
                .collect();
            sample_pops.insert(sample.to_owned(), list);
        }
    }

    // Population array: named pops (file order) then the global "" pop.
    let mut pops: Vec<Pop> = pop_order
        .iter()
        .map(|n| Pop {
            name: n.clone(),
            suffix: format!("_{n}"),
        })
        .collect();
    pops.push(Pop {
        name: String::new(),
        suffix: String::new(),
    });
    let global = pops.len() - 1;

    let text = read_vcf_text(input)?;

    let mut out = String::with_capacity(text.len() + 4096);
    let mut sample_to_pops: Vec<Vec<usize>> = Vec::new();
    let mut gt_warned = false;
    // bcftools inserts a PASS FILTER header (after ##fileformat) when
    // the input lacks one.
    let has_pass = text.contains("##FILTER=<ID=PASS,");

    for line in text.lines() {
        if line.starts_with("##") {
            out.push_str(line);
            out.push('\n');
            if !has_pass && line.starts_with("##fileformat=") {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            }
            continue;
        }
        if let Some(cols) = line.strip_prefix("#CHROM") {
            // Append the new ##INFO header lines (IDs not already
            // declared), in upstream hdr_append × pop order.
            let declared: Vec<String> = collect_declared_info_ids(&out);
            let mut hdr = String::new();
            for &tag in HDR_ORDER {
                if !want.contains(&tag) {
                    continue;
                }
                for p in &pops {
                    let id = format!("{}{}", tag.base_id(), p.suffix);
                    if declared.iter().any(|d| d == &id) {
                        continue;
                    }
                    let (_, tmpl) = tag.header();
                    let in_part = if p.name.is_empty() {
                        String::new()
                    } else {
                        format!(" in {}", p.name)
                    };
                    hdr.push_str(&tmpl.replace("{S}", &p.suffix).replace("{IN}", &in_part));
                    hdr.push('\n');
                }
            }
            out.push_str(&hdr);
            // Map each sample column to its pop indices.
            let samples: Vec<&str> = cols.split('\t').skip(9).collect();
            sample_to_pops = samples
                .iter()
                .map(|s| {
                    let mut v = Vec::new();
                    if let Some(list) = sample_pops.get(*s) {
                        for n in list {
                            if let Some(i) = pops.iter().position(|p| &p.name == n) {
                                v.push(i);
                            }
                        }
                    }
                    v.push(global); // always the global pop
                    v
                })
                .collect();
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        out.push_str(&process_record(
            line,
            &pops,
            &sample_to_pops,
            &want,
            &mut gt_warned,
        ));
        out.push('\n');
    }

    Ok(out)
}

fn collect_declared_info_ids(header: &str) -> Vec<String> {
    let mut ids = Vec::new();
    for l in header.lines() {
        if let Some(rest) = l.strip_prefix("##INFO=<ID=") {
            let id: String = rest.chars().take_while(|&c| c != ',' && c != '>').collect();
            ids.push(id);
        }
    }
    ids
}

fn process_record(
    line: &str,
    pops: &[Pop],
    sample_to_pops: &[Vec<usize>],
    want: &[Tag],
    gt_warned: &mut bool,
) -> String {
    let mut f: Vec<&str> = line.split('\t').collect();
    if f.len() < 10 {
        return line.to_owned();
    }
    let n_allele = if f[4] == "." {
        1
    } else {
        1 + f[4].split(',').count()
    };

    // Locate GT in FORMAT.
    let gt_idx = f[8].split(':').position(|k| k == "GT");
    let Some(gt_idx) = gt_idx else {
        if !*gt_warned {
            *gt_warned = true;
        }
        return line.to_owned();
    };

    // Per-pop, per-allele counts + ns.
    let mut counts: Vec<Vec<Counts>> = pops
        .iter()
        .map(|_| vec![Counts::default(); n_allele])
        .collect();
    let mut ns: Vec<i64> = vec![0; pops.len()];

    for (si, sample) in f[9..].iter().enumerate() {
        let gt = sample.split(':').nth(gt_idx).unwrap_or(".");
        let mut distinct: Vec<usize> = Vec::new();
        let mut nals = 0usize;
        let mut islots = 0usize;
        for tok in gt.split(['/', '|']) {
            islots += 1;
            if tok == "." || tok.is_empty() {
                continue; // missing allele
            }
            let Ok(idx) = tok.parse::<usize>() else {
                continue;
            };
            if idx >= n_allele {
                continue;
            }
            nals += 1;
            if !distinct.contains(&idx) {
                distinct.push(idx);
            }
        }
        if nals == 0 {
            continue; // missing genotype
        }
        let is_hom = distinct.len() == 1;
        // Half-missing (no `-d`) or a single-allele genotype is
        // hemizygous; `is_half` is only ever set under `-d` (deferred).
        let is_half = false;
        let is_hemi = nals != islots || nals == 1;

        let smpl_pops = sample_to_pops.get(si).map(Vec::as_slice).unwrap_or(&[]);
        for &pi in smpl_pops {
            for &a in &distinct {
                let c = &mut counts[pi][a];
                if is_half {
                    c.nac += 1;
                } else if !is_hom {
                    c.nhet += 1;
                } else if !is_hemi {
                    c.nhom += 2;
                } else {
                    c.nhemi += 1;
                }
            }
            ns[pi] += 1;
        }
    }

    // Build INFO key=value additions in the fixed write order.
    let mut info = f[7].to_owned();
    for &tag in WRITE_ORDER {
        if !want.contains(&tag) {
            continue;
        }
        for (pi, p) in pops.iter().enumerate() {
            let c = &counts[pi];
            let total = |a: usize| c[a].nhet + c[a].nhom + c[a].nhemi + c[a].nac;
            match tag {
                Tag::Ns => {
                    set_info(&mut info, &format!("NS{}", p.suffix), &ns[pi].to_string());
                }
                Tag::An => {
                    let an: i64 = (0..n_allele).map(total).sum();
                    set_info(&mut info, &format!("AN{}", p.suffix), &an.to_string());
                }
                // Number=A / Number=1-second tags: when there is no ALT
                // allele upstream's `bcf_update_info_*` is called with 0
                // values, which deletes the tag.
                Tag::Ac => {
                    let key = format!("AC{}", p.suffix);
                    if n_allele > 1 {
                        let v: Vec<String> = (1..n_allele).map(|a| total(a).to_string()).collect();
                        set_info(&mut info, &key, &v.join(","));
                    } else {
                        del_info(&mut info, &key);
                    }
                }
                Tag::AcHet => {
                    let key = format!("AC_Het{}", p.suffix);
                    if n_allele > 1 {
                        let v: Vec<String> = (1..n_allele).map(|a| c[a].nhet.to_string()).collect();
                        set_info(&mut info, &key, &v.join(","));
                    } else {
                        del_info(&mut info, &key);
                    }
                }
                Tag::AcHom => {
                    let key = format!("AC_Hom{}", p.suffix);
                    if n_allele > 1 {
                        let v: Vec<String> = (1..n_allele).map(|a| c[a].nhom.to_string()).collect();
                        set_info(&mut info, &key, &v.join(","));
                    } else {
                        del_info(&mut info, &key);
                    }
                }
                Tag::AcHemi => {
                    let key = format!("AC_Hemi{}", p.suffix);
                    if n_allele > 1 {
                        let v: Vec<String> =
                            (1..n_allele).map(|a| c[a].nhemi.to_string()).collect();
                        set_info(&mut info, &key, &v.join(","));
                    } else {
                        del_info(&mut info, &key);
                    }
                }
                Tag::Af | Tag::Maf => {
                    let key = format!("{}{}", if tag == Tag::Af { "AF" } else { "MAF" }, p.suffix);
                    if n_allele <= 1 {
                        del_info(&mut info, &key);
                        continue;
                    }
                    let mut fr: Vec<f64> = (0..n_allele).map(|a| total(a) as f64).collect();
                    let an: f64 = fr.iter().sum();
                    let missing = an == 0.0;
                    if !missing {
                        for x in &mut fr {
                            *x /= an;
                        }
                    }
                    if tag == Tag::Af {
                        let v: Vec<String> = fr[1..]
                            .iter()
                            .map(|&x| if missing { ".".into() } else { fmt_float(x) })
                            .collect();
                        set_info(&mut info, &key, &v.join(","));
                    } else {
                        // MAF: second most common allele frequency.
                        if !missing {
                            fr.sort_by(|a, b| b.partial_cmp(a).unwrap());
                        }
                        let maf = if missing {
                            ".".into()
                        } else {
                            fmt_float(fr[1])
                        };
                        set_info(&mut info, &key, &maf);
                    }
                }
            }
        }
    }

    f[7] = &info;

    // bcftools re-serializes the record; a lone `.` sample column
    // expands to one `.` per FORMAT subfield.
    let nkeys = f[8].split(':').count();
    let dots = vec!["."; nkeys].join(":");
    if nkeys > 1 {
        for col in &mut f[9..] {
            if *col == "." {
                *col = dots.as_str();
            }
        }
    }
    f.join("\t")
}

/// bcftools float printing: `%.6f` with trailing zeros (and a lone
/// trailing `.`) trimmed.
fn fmt_float(x: f64) -> String {
    let s = format!("{x:.6}");
    let t = s.trim_end_matches('0').trim_end_matches('.');
    if t.is_empty() {
        "0".to_owned()
    } else {
        t.to_owned()
    }
}

/// `bcf_update_info_*`: replace `INFO/<key>` in place if present, else
/// append (preserving the rest of the INFO column).
fn set_info(info: &mut String, key: &str, value: &str) {
    let entry = format!("{key}={value}");
    if *info == "." || info.is_empty() {
        *info = entry;
        return;
    }
    let mut kept: Vec<String> = Vec::new();
    let mut replaced = false;
    for kv in info.split(';') {
        let k = kv.split('=').next().unwrap_or("");
        if k == key {
            kept.push(entry.clone());
            replaced = true;
        } else {
            kept.push(kv.to_owned());
        }
    }
    if !replaced {
        kept.push(entry);
    }
    *info = kept.join(";");
}

/// `bcf_update_info_*` with zero values: remove `INFO/<key>` if present.
fn del_info(info: &mut String, key: &str) {
    if *info == "." || info.is_empty() {
        return;
    }
    let kept: Vec<&str> = info
        .split(';')
        .filter(|kv| kv.split('=').next() != Some(key))
        .collect();
    *info = if kept.is_empty() {
        ".".to_owned()
    } else {
        kept.join(";")
    };
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
        MultiGzDecoder::new(file).read_to_string(&mut text)?;
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
        ".bcftools-rs-fill-tags-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_formatting() {
        assert_eq!(fmt_float(0.5), "0.5");
        assert_eq!(fmt_float(1.0 / 3.0), "0.333333");
        assert_eq!(fmt_float(2.0 / 3.0), "0.666667");
        assert_eq!(fmt_float(0.25), "0.25");
        assert_eq!(fmt_float(1.0), "1");
        assert_eq!(fmt_float(0.0), "0");
    }

    #[test]
    fn set_info_replace_and_append() {
        let mut i = "DP=0".to_owned();
        set_info(&mut i, "AN", "4");
        assert_eq!(i, "DP=0;AN=4");
        set_info(&mut i, "DP", "9");
        assert_eq!(i, "DP=9;AN=4");
        let mut j = ".".to_owned();
        set_info(&mut j, "AC", "2");
        assert_eq!(j, "AC=2");
    }

    #[test]
    fn tag_parsing_aliases() {
        assert!(matches!(parse_tag("AC_Hom"), Ok(Tag::AcHom)));
        assert!(matches!(parse_tag("INFO/AF"), Ok(Tag::Af)));
        assert!(parse_tag("HWE").is_err());
    }
}
