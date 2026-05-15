//! `bcftools +scatter` (upstream `bcftools/plugins/scatter.c`).
//!
//! Splits a VCF into multiple VCFs, either by fixed-size chunks
//! (`-n N`) or by a comma-separated region list (`-s`, with `-x` for an
//! "extra" file holding records that match no region). Each output
//! file gets the full input header. `-i`/`-e` filtering needs the
//! not-yet-ported filter engine (the fixtures do not use it).

use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

struct Region {
    chrom: String,
    beg0: u64,
    end0: u64,
}

fn parse_region(tok: &str) -> Region {
    match tok.split_once(':') {
        None => Region {
            chrom: tok.to_string(),
            beg0: 0,
            end0: u64::MAX,
        },
        Some((chrom, rng)) => {
            let (b, e) = match rng.split_once('-') {
                None => {
                    let b: u64 = rng.replace(',', "").parse().unwrap_or(1);
                    (b, b)
                }
                Some((b, e)) => {
                    let b: u64 = b.replace(',', "").parse().unwrap_or(1);
                    let e = if e.is_empty() {
                        u64::MAX
                    } else {
                        e.replace(',', "").parse().unwrap_or(u64::MAX)
                    };
                    (b, e)
                }
            };
            Region {
                chrom: chrom.to_string(),
                beg0: b.saturating_sub(1),
                end0: if e == u64::MAX { u64::MAX } else { e - 1 },
            }
        }
    }
}

/// Reads the input and writes the scattered VCF files into `output_dir`.
#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &Path,
    output_dir: &Path,
    nsites: Option<usize>,
    scatter: Option<&str>,
    scatter_file: Option<&Path>,
    extra: Option<&str>,
    prefix: Option<&str>,
) -> io::Result<()> {
    let text = read_vcf_text(input)?;
    let scatter_file_text = match scatter_file {
        Some(p) => Some(fs::read_to_string(p)?),
        None => None,
    };
    let files = compute(&text, nsites, scatter, scatter_file_text.as_deref(), extra)
        .map_err(io::Error::other)?;

    fs::create_dir_all(output_dir)?;
    let pre = prefix.unwrap_or("");
    for (fname, content) in files {
        let sanitized: String = fname
            .chars()
            .map(|c| if c.is_whitespace() { '_' } else { c })
            .collect();
        let path = output_dir.join(format!("{pre}{sanitized}.vcf"));
        let mut fh = File::create(&path)?;
        fh.write_all(content.as_bytes())?;
    }
    Ok(())
}

/// Returns the ordered list of `(fname, file_content)` to write.
fn compute(
    text: &str,
    nsites: Option<usize>,
    scatter: Option<&str>,
    scatter_file: Option<&str>,
    extra: Option<&str>,
) -> Result<Vec<(String, String)>, String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut header = String::new();
    for l in &lines {
        if l.starts_with('#') {
            header.push_str(l);
            header.push('\n');
        }
    }
    let records: Vec<&&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .collect();

    let mut out: Vec<(String, String)> = Vec::new();

    if let Some(n) = nsites {
        let n = n.max(1);
        let mut chunk = 0usize;
        let mut cnt = 0usize;
        let mut buf = String::new();
        for rec in &records {
            if cnt == 0 {
                buf.clear();
                buf.push_str(&header);
            }
            buf.push_str(rec);
            buf.push('\n');
            cnt += 1;
            if cnt == n {
                out.push((chunk.to_string(), std::mem::take(&mut buf)));
                chunk += 1;
                cnt = 0;
            }
        }
        if cnt > 0 {
            out.push((chunk.to_string(), buf));
        }
        return Ok(out);
    }

    // Region-scatter mode. Build the ordered set list.
    let mut set_names: Vec<String> = Vec::new();
    let mut regions: Vec<(Region, usize)> = Vec::new(); // (region, set index)
    let mut name_to_idx: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();

    // (region spec, set name) pairs in encounter order.
    let mut specs: Vec<(String, String)> = Vec::new();
    if let Some(sf) = scatter_file {
        for line in sf.lines() {
            let s = line.trim();
            if s.is_empty() || s.starts_with('#') {
                continue;
            }
            let mut it = s.split_whitespace();
            let reg = it.next().unwrap_or("");
            let name = it.next().unwrap_or(reg);
            specs.push((reg.to_string(), name.to_string()));
        }
    } else if let Some(s) = scatter {
        for tok in s.split(',') {
            specs.push((tok.to_string(), tok.to_string()));
        }
    } else {
        return Err("Missing either the -n or one of the -s or -S options".to_string());
    }
    for (spec, name) in &specs {
        let idx = *name_to_idx.entry(name.clone()).or_insert_with(|| {
            set_names.push(name.clone());
            set_names.len() - 1
        });
        regions.push((parse_region(spec), idx));
    }

    let extra_idx = extra.map(|x| {
        set_names.push(x.to_string());
        set_names.len() - 1
    });

    let mut buffers: Vec<String> = set_names.iter().map(|_| header.clone()).collect();

    for rec in &records {
        let f: Vec<&str> = rec.split('\t').collect();
        if f.len() < 2 {
            continue;
        }
        let chrom = f[0];
        let pos0 = f[1].parse::<u64>().unwrap_or(1).saturating_sub(1);
        let mut matched = false;
        for (reg, idx) in &regions {
            if reg.chrom == chrom && pos0 >= reg.beg0 && pos0 <= reg.end0 {
                matched = true;
                buffers[*idx].push_str(rec);
                buffers[*idx].push('\n');
            }
        }
        if !matched && let Some(ei) = extra_idx {
            buffers[ei].push_str(rec);
            buffers[ei].push('\n');
        }
    }

    for (name, buf) in set_names.into_iter().zip(buffers) {
        out.push((name, buf));
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
        ".bcftools-rs-scatter-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chunk_by_n() {
        let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t1\n1\t2\n1\t3\n1\t4\n";
        let out = compute(vcf, Some(3), None, None, None).unwrap();
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].0, "0");
        assert_eq!(out[1].0, "1");
        assert!(out[1].1.ends_with("1\t4\n"));
    }

    #[test]
    fn scatter_regions_and_extra() {
        let vcf = "#CHROM\tPOS\n21\t10\n22\t10\nX\t10\n";
        let out = compute(vcf, None, Some("21,22"), None, Some("X")).unwrap();
        let names: Vec<&str> = out.iter().map(|(n, _)| n.as_str()).collect();
        assert_eq!(names, vec!["21", "22", "X"]);
        assert!(out[0].1.contains("21\t10"));
        assert!(out[2].1.contains("X\t10"));
    }
}
