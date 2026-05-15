//! End-to-end tests for `bcftools_rs::commands::annotate`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

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

fn run(args: &[&str]) -> (Vec<u8>, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (out.stdout, stderr, out.status.code().unwrap_or(-1))
}

fn write_fixture(dir: &TempDir) -> (PathBuf, PathBuf) {
    let vcf = dir.path().join("in.vcf");
    let map = dir.path().join("rename.map");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t2\t.\tA\tC\t.\tPASS\t.\n",
    )
    .unwrap();
    std::fs::write(&map, "1\tchr1\n").unwrap();
    (vcf, map)
}

#[test]
fn annotate_rename_chrs_updates_vcf_text() {
    let dir = TempDir::new().unwrap();
    let (vcf, map) = write_fixture(&dir);

    let (out, err, code) = run(&[
        "annotate",
        "--rename-chrs",
        map.to_str().unwrap(),
        "-Ov",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(out.contains("##contig=<ID=chr1,length=10>"), "{out}");
    assert!(out.contains("chr1\t2\t.\tA\tC\t.\tPASS\t."), "{out}");
}

#[test]
fn annotate_rename_chrs_writes_bcf_that_query_reads() {
    let dir = TempDir::new().unwrap();
    let (vcf, map) = write_fixture(&dir);

    let (bcf, err, code) = run(&[
        "annotate",
        "--rename-chrs",
        map.to_str().unwrap(),
        "-Ob",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate -Ob failed: {err}");

    ensure_binary_built();
    let mut child = Command::new(bin_path())
        .args(["query", "-f%CHROM\\t%POS\\n", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn query");
    child.stdin.as_mut().unwrap().write_all(&bcf).unwrap();
    let out = child.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "query failed: {stderr}");
    assert_eq!(String::from_utf8(out.stdout).unwrap(), "chr1\t2\n");
}

#[test]
fn annotate_rename_chrs_reads_bcf_input() {
    let dir = TempDir::new().unwrap();
    let (vcf, map) = write_fixture(&dir);
    let bcf = dir.path().join("in.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&[
        "annotate",
        "--rename-chrs",
        map.to_str().unwrap(),
        "-Ov",
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate BCF input failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(out.contains("chr1\t2\t.\tA\tC\t.\tPASS\t."), "{out}");
}
