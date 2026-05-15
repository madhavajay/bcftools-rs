//! End-to-end parity tests for `+dosage` against the upstream
//! `dosage.{1,2,3}.out` fixtures (`-t PL`, `-t GL`, `-t GT`).

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

fn check(tag: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("dosage.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args(["+dosage", input.to_str().unwrap(), "--", "-t", tag])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for -t {tag}");
}

#[test]
fn dosage_pl() {
    check("PL", "dosage.1.out");
}

#[test]
fn dosage_gl() {
    check("GL", "dosage.2.out");
}

#[test]
fn dosage_gt() {
    check("GT", "dosage.3.out");
}
