//! End-to-end parity tests for `+fixref` conversion modes against the
//! upstream `fixref.{4,5,6,7}.out` fixtures (test.pl rows 734-737).

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

fn check(mode: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("fixref.2a.vcf");
    let fa = fixture_path("norm.fa");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+fixref",
            input.to_str().unwrap(),
            "--",
            "-f",
            fa.to_str().unwrap(),
            "-m",
            mode,
        ])
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
    assert_eq!(filtered, expected, "mismatch for -m {mode}");
}

#[test]
fn fixref_ref_alt() {
    check("ref-alt", "fixref.4.out");
}

#[test]
fn fixref_flip() {
    check("flip", "fixref.5.out");
}

#[test]
fn fixref_flip_all() {
    check("flip-all", "fixref.6.out");
}

#[test]
fn fixref_swap() {
    check("swap", "fixref.7.out");
}
