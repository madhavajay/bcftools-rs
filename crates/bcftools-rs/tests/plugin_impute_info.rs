//! End-to-end tests for the `+impute-info` plugin.
//!
//! Upstream `test.pl` has no `impute-info` row, so these use synthetic
//! fixtures written to a tempdir and assert the IMPUTE2 `INFO/INFO`
//! annotation, header insertion, passthrough of non-GP / non-triplet
//! sites, and the stderr summary line.

use std::path::PathBuf;
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

fn run_on(vcf: &str, args: &[&str]) -> (String, String, i32) {
    ensure_binary_built();
    let dir = TempDir::new().unwrap();
    let path = dir.path().join("in.vcf");
    std::fs::write(&path, vcf).unwrap();
    let mut full: Vec<&str> = args.to_vec();
    let p = path.to_str().unwrap().to_owned();
    full.push(&p);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
    (
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8_lossy(&out.stderr).into_owned(),
        out.status.code().unwrap_or(-1),
    )
}

const VCF: &str = "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"gp\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t10\t.\tC\tT\t.\tPASS\tDP=5\tGT:GP\t0/0:1,0,0\t1/1:0,0,1\n\
1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n\
1\t30\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/1:0.5,0.5\t0/1:0.5,0.5\n";

#[test]
fn annotates_header_and_records() {
    let (out, err, code) = run_on(VCF, &["+impute-info", "--no-version"]);
    assert_eq!(code, 0, "+impute-info failed: {err}");

    // Header line inserted after the last ##INFO line, before #CHROM.
    let dp = out.find("##INFO=<ID=DP").unwrap();
    let info = out
        .find("##INFO=<ID=INFO,Number=1,Type=Float,Description=\"IMPUTE2 info score\">")
        .expect("INFO header inserted");
    let chrom = out.find("#CHROM").unwrap();
    assert!(dp < info && info < chrom);

    // Perfectly called 0/0 + 1/1 => info score 1.
    assert!(
        out.contains("\tDP=5;INFO=1\tGT:GP\t0/0:1,0,0\t1/1:0,0,1\n"),
        "biallelic-diploid record annotated; got:\n{out}"
    );
    // No FORMAT/GP => unchanged.
    assert!(out.contains("1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t1/1\n"));
    // GP width != 3 => unchanged.
    assert!(out.contains("1\t30\t.\tC\tT\t.\tPASS\t.\tGT:GP\t0/1:0.5,0.5\t0/1:0.5,0.5\n"));
}

#[test]
fn emits_upstream_summary_and_warnings_on_stderr() {
    let (_out, err, code) = run_on(VCF, &["+impute-info", "--no-version"]);
    assert_eq!(code, 0);
    assert!(
        err.contains(
            "Lines total/info-added/unchanged-no-tag/unchanged-not-biallelic-diploid:\t3/1/1/1"
        ),
        "summary line on stderr; got:\n{err}"
    );
    assert!(err.contains("[impute-info.c] Warning: info tag not added to sites without GP tag"));
    assert!(err.contains(
        "[impute-info.c] Warning: info tag not added to sites that are not biallelic diploid"
    ));
}
