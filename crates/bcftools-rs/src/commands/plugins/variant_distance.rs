//! `bcftools +variant-distance` (upstream `bcftools/plugins/variant-distance.c`).
//!
//! Annotates each site with the distance to the nearest variant via an
//! `INFO/DIST` tag (configurable name). Directionality:
//! - `nearest` (default, Number=1): min of the previous/next distinct-position
//!   distances on the same chromosome.
//! - `fwd` (Number=1): distance to the next distinct position, missing if none.
//! - `rev` (Number=1): distance to the previous distinct position, missing if
//!   none.
//! - `both` (Number=2): `<rev>,<fwd>` with `0` substituted for a missing side.
//!
//! Records sharing a POS are duplicates and all receive the same distance.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Nearest,
    Fwd,
    Rev,
    Both,
}

impl Direction {
    pub fn parse(s: &str) -> Option<Self> {
        match s {
            "nearest" => Some(Self::Nearest),
            "fwd" => Some(Self::Fwd),
            "rev" => Some(Self::Rev),
            "both" => Some(Self::Both),
            _ => None,
        }
    }

    fn number(self) -> u8 {
        if self == Self::Both { 2 } else { 1 }
    }

    fn desc(self) -> &'static str {
        match self {
            Self::Nearest => "nearest",
            Self::Fwd => "next",
            Self::Rev => "previous",
            Self::Both => "previous and next",
        }
    }
}

/// Reads the input VCF/BCF and returns the distance-annotated VCF text.
pub fn run(input: &Path, dir: Direction, tag: &str) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(annotate(&text, dir, tag))
}

fn annotate(text: &str, dir: Direction, tag: &str) -> String {
    let lines: Vec<&str> = text.lines().collect();

    // Per-chromosome ordered distinct positions (input assumed coord-sorted).
    let mut distinct: HashMap<&str, Vec<i64>> = HashMap::new();
    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let (Some(chrom), Some(pos)) = (f.next(), f.next()) else {
            continue;
        };
        let Ok(p) = pos.parse::<i64>() else { continue };
        let v = distinct.entry(chrom).or_default();
        if v.last() != Some(&p) {
            v.push(p);
        }
    }

    let info_header = format!(
        "##INFO=<ID={tag},Number={},Type=Integer,Description=\"Distance to the {} variant\">",
        dir.number(),
        dir.desc()
    );
    let last_info = lines
        .iter()
        .rposition(|l| l.starts_with("##INFO="))
        .or_else(|| lines.iter().position(|l| l.starts_with("#CHROM")));
    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));

    let mut out = String::with_capacity(text.len() + 256);
    for (idx, line) in lines.iter().enumerate() {
        if line.starts_with('#') {
            if line.starts_with("#CHROM") && Some(idx) == last_info && !line.starts_with("##") {
                out.push_str(&info_header);
                out.push('\n');
            }
            out.push_str(line);
            out.push('\n');
            if Some(idx) == fileformat && !has_pass {
                out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">");
                out.push('\n');
            }
            if line.starts_with("##INFO=") && Some(idx) == last_info {
                out.push_str(&info_header);
                out.push('\n');
            }
            continue;
        }
        if line.trim().is_empty() {
            out.push('\n');
            continue;
        }
        out.push_str(&annotate_record(line, dir, tag, &distinct));
        out.push('\n');
    }
    out
}

fn annotate_record(
    line: &str,
    dir: Direction,
    tag: &str,
    distinct: &HashMap<&str, Vec<i64>>,
) -> String {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return line.to_owned();
    }
    let chrom = fields[0];
    let Ok(pos) = fields[1].parse::<i64>() else {
        return line.to_owned();
    };
    let Some(positions) = distinct.get(chrom) else {
        return line.to_owned();
    };
    let Ok(k) = positions.binary_search(&pos) else {
        return line.to_owned();
    };
    let rev = if k > 0 {
        Some(pos - positions[k - 1])
    } else {
        None
    };
    let fwd = if k + 1 < positions.len() {
        Some(positions[k + 1] - pos)
    } else {
        None
    };

    let value: Option<String> = match dir {
        Direction::Nearest => match (rev, fwd) {
            (Some(r), Some(f)) => Some(r.min(f).to_string()),
            (Some(r), None) => Some(r.to_string()),
            (None, Some(f)) => Some(f.to_string()),
            (None, None) => None,
        },
        Direction::Fwd => fwd.map(|f| f.to_string()),
        Direction::Rev => rev.map(|r| r.to_string()),
        Direction::Both => Some(format!("{},{}", rev.unwrap_or(0), fwd.unwrap_or(0))),
    };

    let stripped = strip_tag(fields[7], tag);
    let new_info = match value {
        Some(v) => {
            let kv = format!("{tag}={v}");
            if stripped == "." || stripped.is_empty() {
                kv
            } else {
                format!("{stripped};{kv}")
            }
        }
        None => {
            if stripped.is_empty() {
                ".".to_owned()
            } else {
                stripped
            }
        }
    };

    let mut out = fields.clone();
    out[7] = new_info.as_str();
    out.join("\t")
}

fn strip_tag(info: &str, tag: &str) -> String {
    if info == "." {
        return info.to_owned();
    }
    info.split(';')
        .filter(|kv| {
            let key = kv.split_once('=').map(|(k, _)| k).unwrap_or(kv);
            key != tag
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
        ".bcftools-rs-variant-distance-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t101\t.\tT\tA\t.\t.\t.\n\
1\t101\t.\tT\tC\t.\t.\t.\n\
1\t104\t.\tT\tA\t.\t.\t.\n\
1\t115\t.\tT\tA\t.\t.\t.\n";

    #[test]
    fn nearest_uses_min_of_both_sides() {
        let out = annotate(VCF, Direction::Nearest, "DIST");
        let data: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert!(data[0].ends_with("\tDIST=3"), "{}", data[0]); // 101 -> 104
        assert!(data[1].ends_with("\tDIST=3"), "{}", data[1]);
        assert!(data[2].ends_with("\tDIST=3"), "{}", data[2]); // 104 -> 101
        assert!(data[3].ends_with("\tDIST=11"), "{}", data[3]); // 115 -> 104
    }

    #[test]
    fn fwd_missing_at_chrom_end() {
        let out = annotate(VCF, Direction::Fwd, "DIST");
        let data: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert!(data[0].ends_with("\tDIST=3"), "{}", data[0]);
        assert!(data[3].ends_with("\t."), "115 has no fwd: {}", data[3]);
    }

    #[test]
    fn rev_missing_at_chrom_start() {
        let out = annotate(VCF, Direction::Rev, "DIST");
        let data: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert!(data[0].ends_with("\t."), "101 has no rev: {}", data[0]);
        assert!(data[3].ends_with("\tDIST=11"), "{}", data[3]);
    }

    #[test]
    fn both_substitutes_zero_for_missing_side() {
        let out = annotate(VCF, Direction::Both, "DIST");
        let data: Vec<&str> = out.lines().filter(|l| !l.starts_with('#')).collect();
        assert!(data[0].ends_with("\tDIST=0,3"), "{}", data[0]);
        assert!(data[3].ends_with("\tDIST=11,0"), "{}", data[3]);
        assert!(out.contains(
            "##INFO=<ID=DIST,Number=2,Type=Integer,Description=\"Distance to the previous and next variant\">"
        ));
    }

    #[test]
    fn injects_pass_filter_header_when_absent() {
        let out = annotate(VCF, Direction::Nearest, "DIST");
        let ff = out.find("##fileformat=").unwrap();
        let pass = out
            .find("##FILTER=<ID=PASS,Description=\"All filters passed\">")
            .unwrap();
        let contig = out.find("##contig=").unwrap();
        assert!(ff < pass && pass < contig, "PASS must follow fileformat");
    }
}
