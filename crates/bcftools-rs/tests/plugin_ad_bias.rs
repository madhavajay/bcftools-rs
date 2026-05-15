//! End-to-end parity test for `+ad-bias` (report mode) against the upstream
//! `ad-bias.out` fixture. Per test.pl, both `ad-bias.vcf` and `ad-bias.2.vcf`
//! produce `ad-bias.out` with `-s ad-bias.samples`; the harness pipes
//! `grep -v bcftools`.

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

fn check(input_vcf: &str) {
    ensure_binary_built();
    let input = fixture_path(input_vcf);
    let samples = fixture_path("ad-bias.samples");
    let expected = std::fs::read_to_string(fixture_path("ad-bias.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+ad-bias",
            input.to_str().unwrap(),
            "--",
            "-s",
            samples.to_str().unwrap(),
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
        .filter(|l| !l.contains("bcftools"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {input_vcf}");
}

#[test]
fn ad_bias_matches_fixture_input1() {
    check("ad-bias.vcf");
}

#[test]
fn ad_bias_matches_fixture_input2() {
    check("ad-bias.2.vcf");
}
