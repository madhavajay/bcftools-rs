//! `bcftools +variantkey-hex` (upstream `bcftools/plugins/variantkey-hex.c`).
//!
//! Generates the three unsorted VariantKey lookup-table files in `dir`
//! (`vkrs.unsorted.hex`, `rsvk.unsorted.hex`, `nrvk.unsorted.tsv`) and
//! returns the upstream summary printed by `destroy()`. The VCF/BCF header
//! and records are suppressed (upstream `init` returns 1, `process` returns
//! NULL), so only the summary reaches stdout.

use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use super::variantkey::{rsid_u32, variantkey};
use crate::vcf_compat::normalize_vcf_text;

const FILE_VKRS: &str = "vkrs.unsorted.hex";
const FILE_RSVK: &str = "rsvk.unsorted.hex";
const FILE_NRVK: &str = "nrvk.unsorted.tsv";

/// Reads the input VCF/BCF, writes the lookup files under `dir` (joined by
/// raw string concatenation like upstream's `strcat`, so a trailing `/` in
/// `dir` is required for a directory), and returns the summary text.
pub fn run(input: &Path, dir: &str) -> io::Result<String> {
    let text = read_vcf_text(input)?;

    let mut fp_vkrs = File::create(format!("{dir}{FILE_VKRS}"))?;
    let mut fp_rsvk = File::create(format!("{dir}{FILE_RSVK}"))?;
    let mut fp_nrvk = File::create(format!("{dir}{FILE_NRVK}"))?;

    let mut numvar: u64 = 0;
    let mut nrv: u64 = 0;

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 5 {
            continue;
        }
        let chrom = fields[0];
        let Ok(pos1) = fields[1].parse::<u32>() else {
            continue;
        };
        let pos0 = pos1.wrapping_sub(1);
        let id = fields[2];
        let reference = fields[3];
        let alt = fields[4].split(',').next().unwrap_or(fields[4]);

        let vk = variantkey(chrom.as_bytes(), pos0, reference.as_bytes(), alt.as_bytes());
        let rs = rsid_u32(id);

        writeln!(fp_vkrs, "{vk:016x}\t{rs:08x}")?;
        writeln!(fp_rsvk, "{rs:08x}\t{vk:016x}")?;
        if vk & 1 == 1 {
            writeln!(fp_nrvk, "{vk:016x}\t{reference}\t{alt}")?;
            nrv += 1;
        }
        numvar += 1;
    }

    fp_vkrs.flush()?;
    fp_rsvk.flush()?;
    fp_nrvk.flush()?;

    Ok(format!(
        "VariantKeys: {numvar}\nNon-reversible VariantKeys: {nrv}\n"
    ))
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
        ".bcftools-rs-variantkey-hex-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_dir() -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let d = std::env::temp_dir().join(format!("vkhex-test-{}-{nanos}", std::process::id()));
        fs::create_dir_all(&d).unwrap();
        d
    }

    #[test]
    fn writes_lookup_files_and_summary() {
        let dir = tmp_dir();
        let vcf = dir.join("in.vcf");
        fs::write(
            &vcf,
            "##fileformat=VCFv4.0\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10019\trs775809821\tTA\tT\t.\t.\t.\n\
1\t10228\trs200462216\tTAACCCCTAACCCTAACCCTAAACCCTA\tT\t.\t.\t.\n",
        )
        .unwrap();
        let prefix = format!("{}/", dir.display());
        let summary = run(&vcf, &prefix).unwrap();
        assert_eq!(summary, "VariantKeys: 2\nNon-reversible VariantKeys: 1\n");

        let vkrs = fs::read_to_string(dir.join(FILE_VKRS)).unwrap();
        assert_eq!(
            vkrs,
            "0800139110e60000\t2e3deb1d\n080013f9a00e1d03\t0bf2cf88\n"
        );
        let rsvk = fs::read_to_string(dir.join(FILE_RSVK)).unwrap();
        assert_eq!(
            rsvk,
            "2e3deb1d\t0800139110e60000\n0bf2cf88\t080013f9a00e1d03\n"
        );
        let nrvk = fs::read_to_string(dir.join(FILE_NRVK)).unwrap();
        assert_eq!(nrvk, "080013f9a00e1d03\tTAACCCCTAACCCTAACCCTAAACCCTA\tT\n");

        let _ = fs::remove_dir_all(&dir);
    }
}
