//! End-to-end parity tests for `+mendelian2` against the upstream
//! `mendelian.{1,3,4,6,7,8}.out` fixtures (test.pl rows 744-749).

use std::path::PathBuf;
use std::process::Command;

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

/// `grep_hash`: strip lines starting with `#` (test.pl row 749 only).
fn check(args: &[&str], expected_fixture: &str, grep_hash: bool) {
    ensure_binary_built();
    let input = fixture_path("mendelian.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec!["+mendelian2", input.to_str().unwrap()];
    full.extend_from_slice(args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !(l.starts_with("##bcftools_") || grep_hash && l.starts_with('#')))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn mendelian2_delete() {
    check(&["-p", "child1,dad1,mom1", "-md"], "mendelian.1.out", false);
}

#[test]
fn mendelian2_list_good() {
    check(&["-p", "child1,dad1,mom1", "-mg"], "mendelian.6.out", false);
}

#[test]
fn mendelian2_list_err() {
    check(&["-p", "child1,dad1,mom1", "-me"], "mendelian.3.out", false);
}

#[test]
fn mendelian2_annotate() {
    check(&["-p", "child1,dad1,mom1", "-ma"], "mendelian.4.out", false);
}

#[test]
fn mendelian2_list_miss() {
    check(&["-p", "child1,dad1,mom1", "-mm"], "mendelian.7.out", false);
}

#[test]
fn mendelian2_count() {
    check(&["-p", "child1,dad1,mom1"], "mendelian.8.out", true);
}

#[test]
fn mendelian2_include_filter_counts_failed_sites() {
    ensure_binary_built();
    let input = fixture_path("mendelian.vcf");
    let out = Command::new(bin_path())
        .args([
            "+mendelian2",
            input.to_str().unwrap(),
            "-p",
            "child1,dad1,mom1",
            "-i",
            "CHROM=\"1\" && POS=100",
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "include filter failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(stdout.contains("sites_fail\t29\t# skipped because of failed -i/-e filter\n"));
    assert!(stdout.contains("sites_good\t1\t# number of sites with at least one good trio\n"));
    assert!(
        stdout.contains("sites_merr\t0\t# number of sites with at least one Mendelian error\n")
    );
    assert!(stdout.contains("ngood\t1\n"));
    assert!(stdout.contains("nfail\t0\n"));
}
