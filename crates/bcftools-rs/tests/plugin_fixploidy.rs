//! End-to-end parity test for `+fixploidy` against the upstream
//! `fixploidy.out` fixture (test.pl row 673).

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

#[test]
fn fixploidy_samples_and_ploidy() {
    ensure_binary_built();
    let input = fixture_path("fixploidy.vcf");
    let samples = fixture_path("fixploidy.samples");
    let ploidy = fixture_path("fixploidy.ploidy");
    let expected = std::fs::read_to_string(fixture_path("fixploidy.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+fixploidy",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-s",
            samples.to_str().unwrap(),
            "-p",
            ploidy.to_str().unwrap(),
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "fixploidy failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for fixploidy");
}
