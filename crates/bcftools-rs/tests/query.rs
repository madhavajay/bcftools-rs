//! End-to-end tests for `bcftools_rs::commands::query`.

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
fn query_list_samples_from_vcf() {
    let path = fixture_path("annotate2.vcf");
    let (out, err, code) = run(&["query", "-l", path.to_str().unwrap()]);
    assert_eq!(code, 0, "query -l failed: {err}");
    assert_eq!(out, "A\nB\nC\n");
}

#[test]
fn query_list_samples_from_bcf() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("samples.vcf");
    let bcf = dir.path().join("annotate2.bcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["query", "--list-samples", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "query --list-samples BCF failed: {err}");
    assert_eq!(out, "A\nB\nC\n");
}

#[test]
fn query_format_core_fields_from_vcf() {
    let path = fixture_path("annotate2.vcf");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM\\t%POS\\t%REF\\t%ALT\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f failed: {err}");
    assert!(out.starts_with("1\t3000001\tC\tT\n"));
    assert!(out.lines().count() > 1);
}
