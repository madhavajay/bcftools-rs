//! End-to-end parity tests for `+guess-ploidy -v -rX` against the upstream
//! `guess-ploidy.{PL,GL}.out` fixtures (harness pipes `grep -v bcftools`).

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

fn check(input_vcf: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input_vcf);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args(["+guess-ploidy", input.to_str().unwrap(), "-v", "-rX"])
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
        .filter(|l| !l.contains("bcftools"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {input_vcf}");
}

#[test]
fn guess_ploidy_pl() {
    check("view.PL.vcf", "guess-ploidy.PL.out");
}

#[test]
fn guess_ploidy_gl() {
    check("view.GL.vcf", "guess-ploidy.GL.out");
}
