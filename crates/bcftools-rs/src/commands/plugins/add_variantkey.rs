//! `bcftools +add-variantkey` (upstream `bcftools/plugins/add-variantkey.c`).
//!
//! Appends two INFO fields to every record: `VKX` (16-hex-digit 64-bit
//! VariantKey over CHROM, 0-based POS, REF and the first ALT allele) and
//! `RSX` (8-hex-digit form of the numeric part of the `rs` ID). The two
//! `##INFO` definitions are injected immediately before `#CHROM`, matching
//! upstream's `bcf_hdr_append` ordering after the harness strips
//! `##bcftools_*` lines.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use super::variantkey::{rsid_u32, rsx_hex, variantkey, variantkey_hex};
use crate::vcf_compat::normalize_vcf_text;

const VKX_HDR: &str = "##INFO=<ID=VKX,Number=1,Type=String,Description=\"Hexadecimal representation of 64 bit VariantKey\">";
const RSX_HDR: &str = "##INFO=<ID=RSX,Number=1,Type=String,Description=\"Hexadecimal representation of ID minus the 'rs' prefix (32bit)\">";

/// Reads the input VCF/BCF and returns the VKX/RSX-annotated VCF text.
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(rewrite(&text))
}

fn rewrite(text: &str) -> String {
    let mut out = String::with_capacity(text.len() + text.len() / 8);
    for line in text.lines() {
        if line.starts_with("#CHROM") {
            out.push_str(VKX_HDR);
            out.push('\n');
            out.push_str(RSX_HDR);
            out.push('\n');
            out.push_str(line);
            out.push('\n');
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            out.push_str(line);
            out.push('\n');
            continue;
        }
        out.push_str(&rewrite_record(line));
        out.push('\n');
    }
    out
}

fn rewrite_record(line: &str) -> String {
    let mut fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 8 {
        return line.to_owned();
    }
    let chrom = fields[0];
    // POS is 1-based in VCF; the algorithm wants the 0-based position.
    let pos0 = fields[1].parse::<u32>().map(|p| p.wrapping_sub(1)).ok();
    let id = fields[2];
    let reference = fields[3];
    let alt = fields[4].split(',').next().unwrap_or(fields[4]);

    let Some(pos0) = pos0 else {
        return line.to_owned();
    };

    let vk = variantkey(chrom.as_bytes(), pos0, reference.as_bytes(), alt.as_bytes());
    let vkx = variantkey_hex(vk);
    let rsx = rsx_hex(rsid_u32(id));

    let info = fields[7];
    let new_info = if info == "." || info.is_empty() {
        format!("VKX={vkx};RSX={rsx}")
    } else {
        format!("{info};VKX={vkx};RSX={rsx}")
    };
    fields[7] = new_info.as_str();
    fields.join("\t")
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
        ".bcftools-rs-add-variantkey-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn appends_vkx_rsx_to_info() {
        let line = "1\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;ASP";
        let got = rewrite_record(line);
        assert_eq!(
            got,
            "1\t10019\trs775809821\tTA\tT\t.\t.\tRS=775809821;ASP;VKX=0800139110e60000;RSX=2e3deb1d"
        );
    }

    #[test]
    fn empty_info_does_not_get_leading_dot() {
        let line = "1\t10019\trs775809821\tTA\tT\t.\t.\t.";
        let got = rewrite_record(line);
        assert_eq!(
            got,
            "1\t10019\trs775809821\tTA\tT\t.\t.\tVKX=0800139110e60000;RSX=2e3deb1d"
        );
    }

    #[test]
    fn first_alt_allele_used_for_multiallelic() {
        let line = "1\t10019\trs775809821\tTA\tT,TAA\t.\t.\t.";
        assert!(rewrite_record(line).contains("VKX=0800139110e60000;"));
    }

    #[test]
    fn header_injection_before_chrom() {
        let vcf = "##fileformat=VCFv4.0\n\
##bcftools_normVersion=1.6\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10019\trs775809821\tTA\tT\t.\t.\t.\n";
        let out = rewrite(vcf);
        let lines: Vec<&str> = out.lines().collect();
        // VKX/RSX INFO lines land immediately before #CHROM.
        let chrom_idx = lines.iter().position(|l| l.starts_with("#CHROM")).unwrap();
        assert!(lines[chrom_idx - 2].starts_with("##INFO=<ID=VKX"));
        assert!(lines[chrom_idx - 1].starts_with("##INFO=<ID=RSX"));
        // The pre-existing ##bcftools_ line is preserved (the harness greps
        // it out downstream, not this plugin).
        assert!(lines.iter().any(|l| l.starts_with("##bcftools_norm")));
    }
}
