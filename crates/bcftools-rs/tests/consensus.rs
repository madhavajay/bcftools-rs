//! End-to-end tests for `bcftools_rs::commands::consensus`.

use std::path::PathBuf;
use std::process::Command;

use tempfile::TempDir;

fn fixture_path(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("bcftools");
    p.push("test");
    p.push(name);
    p
}

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

#[test]
fn consensus_applies_simple_alt_alleles_to_fasta() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&fasta, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tC\tTT\t.\tPASS\t.\n\
chr1\t7\t.\tG\tA\t.\tPASS\t.\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "consensus",
        "-f",
        fasta.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "consensus failed: {err}");
    assert_eq!(out, ">chr1\nATTGTACAT\n");
}

#[test]
fn consensus_empty_vcf_preserves_upstream_fixture_fasta() {
    let fasta = fixture_path("consensus.fa");
    let vcf = fixture_path("empty.vcf");
    let expected = std::fs::read_to_string(fixture_path("consensus.5.out")).unwrap();

    let (out, err, code) = run(&[
        "consensus",
        "-s",
        "-",
        "-f",
        fasta.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "consensus empty fixture failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn consensus_reads_bcf_input() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let vcf = dir.path().join("in.vcf");
    let bcf = dir.path().join("in.bcf");
    std::fs::write(&fasta, ">chr1\nACGTACGT\n").unwrap();
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=8>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t4\t.\tT\tG\t.\tPASS\t.\n",
    )
    .unwrap();

    let (_view_out, view_err, view_code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(view_code, 0, "view -Ob failed: {view_err}");

    let (out, err, code) = run(&[
        "consensus",
        "-f",
        fasta.to_str().unwrap(),
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "consensus BCF failed: {err}");
    assert_eq!(out, ">chr1\nACGGACGT\n");
}

#[test]
fn consensus_reference_mismatch_errors() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let vcf = dir.path().join("in.vcf");
    std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\t.\tA\tT\t.\tPASS\t.\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "consensus",
        "-f",
        fasta.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_ne!(code, 0);
    assert!(err.contains("reference mismatch"), "{err}");
}
