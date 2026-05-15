//! End-to-end tests for the `+counts` plugin.

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

fn run(args: &[&str]) -> (String, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stdout = String::from_utf8(out.stdout).unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

fn write_vcf(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("in.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\n\
1\t20\t.\tA\tAT\t.\t.\t.\tGT\t0/1\t0/0\n\
1\t30\t.\tAC\tGT\t.\t.\t.\tGT\t0/1\t0/0\n",
    )
    .unwrap();
    path
}

const EXPECTED: &str = "Number of samples: 2\n\
Number of SNPs:    1\n\
Number of INDELs:  1\n\
Number of MNPs:    1\n\
Number of others:  0\n\
Number of sites:   3\n";

#[test]
fn plugin_counts_via_plus_shorthand() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir);
    let (out, err, code) = run(&["+counts", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "+counts failed: {err}");
    assert_eq!(out, EXPECTED);
}

#[test]
fn plugin_counts_via_plugin_subcommand() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir);
    let (out, err, code) = run(&["plugin", "counts", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "plugin counts failed: {err}");
    assert_eq!(out, EXPECTED);
}

#[test]
fn plugin_counts_reads_bcf_input() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir);
    let bcf = dir.path().join("in.bcf");

    let (_o, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["+counts", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "+counts on bcf failed: {err}");
    assert_eq!(out, EXPECTED);
}

#[test]
fn plugin_counts_reads_stdin() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir);
    let data = std::fs::read(&vcf).unwrap();

    ensure_binary_built();
    let mut child = Command::new(bin_path())
        .args(["+counts", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn +counts");
    child.stdin.as_mut().unwrap().write_all(&data).unwrap();
    let out = child.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "+counts - failed: {stderr}");
    assert_eq!(String::from_utf8(out.stdout).unwrap(), EXPECTED);
}

#[test]
fn plugin_counts_help_still_lists_registry() {
    let (out, _err, code) = run(&["plugin", "-l"]);
    assert_eq!(code, 0);
    assert!(out.contains("counts"), "registry listing missing counts");
}
