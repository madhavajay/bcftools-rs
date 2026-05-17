//! End-to-end parity tests for `+parental-origin` against the upstream
//! `parental-origin.{1..5}.out` fixtures (test.pl rows 857-861; the
//! harness pipes `grep -v ^#`).

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

fn check(region: &str, cnv: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("parental-origin.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+parental-origin",
            input.to_str().unwrap(),
            "-r",
            region,
            "-p",
            "proband,father,mother",
            "-t",
            cnv,
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{region} {cnv} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {region} {cnv}");
}

#[test]
fn parental_origin_del_1() {
    check("20:100", "del", "parental-origin.1.out");
}

#[test]
fn parental_origin_del_2() {
    check("20:101", "del", "parental-origin.2.out");
}

#[test]
fn parental_origin_del_3() {
    check("20:102", "del", "parental-origin.3.out");
}

#[test]
fn parental_origin_dup_4() {
    check("20:103", "dup", "parental-origin.4.out");
}

#[test]
fn parental_origin_dup_5() {
    check("20:104", "dup", "parental-origin.5.out");
}

#[test]
fn parental_origin_include_filter_limits_informative_sites() {
    ensure_binary_built();
    let input = fixture_path("parental-origin.vcf");
    let expected = std::fs::read_to_string(fixture_path("parental-origin.1.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+parental-origin",
            input.to_str().unwrap(),
            "-r",
            "20:100-102",
            "-p",
            "proband,father,mother",
            "-t",
            "del",
            "-i",
            "POS=100",
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
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected);
}
