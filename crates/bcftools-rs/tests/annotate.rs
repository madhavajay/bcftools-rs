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

fn write_annotated_fixture(dir: &TempDir) -> PathBuf {
    let vcf = dir.path().join("ann.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=10>\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"x\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"x\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"x\">\n\
##FILTER=<ID=LowQual,Description=\"x\">\n\
##FILTER=<ID=q10,Description=\"x\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t5\trs9\tA\tC\t99\tLowQual;q10\tAC=1;AN=2;DP=12\n",
    )
    .unwrap();
    vcf
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

#[test]
fn annotate_remove_specific_info_tags() {
    let dir = TempDir::new().unwrap();
    let vcf = write_annotated_fixture(&dir);

    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "INFO/AC,INFO/DP",
        "-Ov",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate -x failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(!out.contains("##INFO=<ID=AC,"), "{out}");
    assert!(!out.contains("##INFO=<ID=DP,"), "{out}");
    assert!(out.contains("##INFO=<ID=AN,"), "{out}");
    assert!(out.contains("\tAN=2\n"), "{out}");
}

#[test]
fn annotate_remove_id_and_qual() {
    let dir = TempDir::new().unwrap();
    let vcf = write_annotated_fixture(&dir);

    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "ID,QUAL",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate -x ID,QUAL failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(
        out.contains("1\t5\t.\tA\tC\t.\tLowQual;q10\tAC=1;AN=2;DP=12"),
        "{out}"
    );
}

#[test]
fn annotate_remove_specific_filter_becomes_pass() {
    let dir = TempDir::new().unwrap();
    let vcf = write_annotated_fixture(&dir);

    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "FILTER/LowQual,FILTER/q10",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate -x FILTER/... failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(!out.contains("##FILTER=<ID=LowQual,"), "{out}");
    assert!(out.contains("\tPASS\tAC=1;AN=2;DP=12"), "{out}");
}

#[test]
fn annotate_remove_all_info() {
    let dir = TempDir::new().unwrap();
    let vcf = write_annotated_fixture(&dir);

    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "INFO",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate -x INFO failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(!out.contains("##INFO="), "{out}");
    assert!(out.contains("\tLowQual;q10\t.\n"), "{out}");
}

#[test]
fn annotate_rename_and_remove_combined() {
    let dir = TempDir::new().unwrap();
    let vcf = write_annotated_fixture(&dir);
    let map = dir.path().join("m.map");
    std::fs::write(&map, "1\tchr1\n").unwrap();

    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "--rename-chrs",
        map.to_str().unwrap(),
        "-x",
        "INFO/AC",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate combined failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(out.contains("##contig=<ID=chr1,"), "{out}");
    assert!(!out.contains("##INFO=<ID=AC,"), "{out}");
    assert!(
        out.contains("chr1\t5\trs9\tA\tC\t99\tLowQual;q10\tAN=2;DP=12"),
        "{out}"
    );
}

#[test]
fn annotate_keep_only_remove_fixture_matches_upstream_text_output() {
    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "ID,QUAL,^FILTER/fltA,FILTER/fltB,^INFO/AA,INFO/BB,^FMT/GT,FMT/PL",
        "../../bcftools/test/annotate3.vcf",
    ]);
    assert_eq!(code, 0, "annotate6 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/annotate6.out").unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn annotate_format_remove_fixture_matches_upstream_text_output() {
    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "FORMAT",
        "../../bcftools/test/annotate3.vcf",
    ]);
    assert_eq!(code, 0, "annotate7 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/annotate7.out").unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn annotate_force_remove_unknown_tags_fixture_matches_upstream_text_output() {
    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "FILTER/XX,INFO/XX",
        "--force",
        "../../bcftools/test/annotate14.vcf",
    ]);
    assert_eq!(code, 0, "annotate25 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/annotate25.out").unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn annotate_remove_filter_fixture_matches_upstream_text_output() {
    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "FILTER",
        "../../bcftools/test/annotate16.vcf",
    ]);
    assert_eq!(code, 0, "annotate28 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/annotate28.out").unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn annotate_keep_only_form_on_local_fixture() {
    let dir = TempDir::new().unwrap();
    let vcf = write_annotated_fixture(&dir);

    let (out, err, code) = run(&[
        "annotate",
        "--no-version",
        "-x",
        "^INFO/AC,INFO/DP",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "annotate keep-only failed: {err}");
    let out = String::from_utf8(out).unwrap();
    assert!(out.contains("##INFO=<ID=AC,"), "{out}");
    assert!(out.contains("##INFO=<ID=DP,"), "{out}");
    assert!(!out.contains("##INFO=<ID=AN,"), "{out}");
    assert!(out.contains("\tAC=1;DP=12\n"), "{out}");
}
