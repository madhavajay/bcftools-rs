//! End-to-end parity tests for `+GTsubset` against the upstream
//! `view.GTsubset.*.out` fixtures (test.pl rows 741-743).

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

fn check(sample_spec: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("view.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+GTsubset",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-s",
            sample_spec,
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{sample_spec} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for -s {sample_spec}");
}

#[test]
fn gtsubset_na1() {
    check("NA00001", "view.GTsubset.NA1.out");
}

#[test]
fn gtsubset_na1_na2() {
    check("NA00001,NA00002", "view.GTsubset.NA1NA2.out");
}

#[test]
fn gtsubset_na1_na2_na3() {
    check("NA00001,NA00002,NA00003", "view.GTsubset.NA1NA2NA3.out");
}
