//! End-to-end tests for `bcftools_rs::commands::tabix`.

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

use tempfile::TempDir;

fn bin_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("target");
    p.push("debug");
    p.push("bcftools");
    p
}

fn ensure_binary_built() {
    let status = Command::new(env!("CARGO"))
        .args(["build", "-p", "bcftools-rs-cli"])
        .status()
        .expect("cargo build");
    assert!(status.success(), "failed to build bcftools-rs-cli");
}

fn run(args: &[&str]) -> (String, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

fn write_bgzf(path: &Path, text: &str) {
    let file = std::fs::File::create(path).unwrap();
    let mut writer = htslib_rs::bgzf::io::Writer::new(file);
    writer.write_all(text.as_bytes()).unwrap();
    writer.finish().unwrap();
}

const VCF: &str = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tG\t.\tPASS\t.\n\
1\t20\t.\tC\tT\t.\tPASS\t.\n\
2\t15\t.\tG\tA\t.\tPASS\t.\n";

#[test]
fn tabix_builds_tbi_and_queries_vcf_regions() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.vcf.gz");
    write_bgzf(&input, VCF);

    let (_out, err, code) = run(&["tabix", "-f", "-p", "vcf", input.to_str().unwrap()]);
    assert_eq!(code, 0, "tabix index failed: {err}");
    assert!(dir.path().join("in.vcf.gz.tbi").exists());

    let (out, err, code) = run(&["tabix", input.to_str().unwrap(), "1:11-20"]);
    assert_eq!(code, 0, "tabix query failed: {err}");
    assert_eq!(out, "1\t20\t.\tC\tT\t.\tPASS\t.\n");
}

#[test]
fn tabix_all_streams_all_bgzf_lines() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.vcf.gz");
    write_bgzf(&input, VCF);

    let (out, err, code) = run(&["tabix", "-a", input.to_str().unwrap()]);
    assert_eq!(code, 0, "tabix -a failed: {err}");
    assert_eq!(out, VCF);
}

#[test]
fn tabix_refuses_existing_index_without_force() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.vcf.gz");
    write_bgzf(&input, VCF);
    std::fs::write(dir.path().join("in.vcf.gz.tbi"), b"old").unwrap();

    let (_out, err, code) = run(&["tabix", "-p", "vcf", input.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("the index file exists"));
}

#[test]
fn tabix_builds_sam_preset_csi() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.sam.gz");
    write_bgzf(
        &input,
        "@SQ\tSN:chr1\tLN:100\nr1\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\t*\n",
    );

    let (_out, err, code) = run(&[
        "tabix",
        "-f",
        "-p",
        "sam",
        "-m",
        "14",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "tabix SAM CSI failed: {err}");
    assert!(dir.path().join("in.sam.gz.csi").exists());

    let (out, err, code) = run(&["tabix", input.to_str().unwrap(), "chr1:10-10"]);
    assert_eq!(code, 0, "tabix SAM query failed: {err}");
    assert_eq!(out, "r1\t0\tchr1\t10\t60\t1M\t*\t0\t0\tA\t*\n");
}

#[test]
fn tabix_builds_bed_and_gff_presets() {
    let dir = TempDir::new().expect("tempdir");
    let bed = dir.path().join("in.bed.gz");
    let gff = dir.path().join("in.gff.gz");
    write_bgzf(&bed, "chr1\t10\t20\tbed-hit\nchr1\t30\t40\tbed-miss\n");
    write_bgzf(
        &gff,
        "chr1\tsrc\tgene\t11\t20\t.\t+\t.\tID=gff-hit\nchr1\tsrc\tgene\t30\t40\t.\t+\t.\tID=gff-miss\n",
    );

    let (_out, err, code) = run(&["tabix", "-f", "-p", "bed", bed.to_str().unwrap()]);
    assert_eq!(code, 0, "tabix BED index failed: {err}");
    let (out, err, code) = run(&["tabix", bed.to_str().unwrap(), "chr1:11-20"]);
    assert_eq!(code, 0, "tabix BED query failed: {err}");
    assert_eq!(out, "chr1\t10\t20\tbed-hit\n");

    let (_out, err, code) = run(&["tabix", "-f", "-p", "gff", gff.to_str().unwrap()]);
    assert_eq!(code, 0, "tabix GFF index failed: {err}");
    let (out, err, code) = run(&["tabix", gff.to_str().unwrap(), "chr1:11-20"]);
    assert_eq!(code, 0, "tabix GFF query failed: {err}");
    assert_eq!(out, "chr1\tsrc\tgene\t11\t20\t.\t+\t.\tID=gff-hit\n");
}
