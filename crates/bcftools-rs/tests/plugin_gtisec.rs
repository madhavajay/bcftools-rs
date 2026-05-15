//! End-to-end parity tests for `+GTisec` against the upstream
//! `view.GTisec.*.out` fixtures (test.pl rows 712-719; the harness
//! pipes `grep -v bcftools`).

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

fn check(flag: Option<&str>, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("view.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec!["+GTisec".to_string(), input.to_str().unwrap().to_string()];
    if let Some(fl) = flag {
        full.push("--".to_string());
        full.push(fl.to_string());
    }
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
    // Harness pipes `grep -v bcftools`.
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.contains("bcftools"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn gtisec_plain() {
    check(None, "view.GTisec.out");
}

#[test]
fn gtisec_h() {
    check(Some("-H"), "view.GTisec.H.out");
}

#[test]
fn gtisec_hm() {
    check(Some("-Hm"), "view.GTisec.Hm.out");
}

#[test]
fn gtisec_hmv() {
    check(Some("-Hmv"), "view.GTisec.Hmv.out");
}

#[test]
fn gtisec_hv() {
    check(Some("-Hv"), "view.GTisec.Hv.out");
}

#[test]
fn gtisec_m() {
    check(Some("-m"), "view.GTisec.m.out");
}

#[test]
fn gtisec_mv() {
    check(Some("-mv"), "view.GTisec.mv.out");
}

#[test]
fn gtisec_v() {
    check(Some("-v"), "view.GTisec.v.out");
}
