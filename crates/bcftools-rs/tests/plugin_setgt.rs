//! End-to-end test for the `+setGT` plugin (first slice).
//!
//! Mirrors the upstream `test_vcf_plugin` row
//! `in=>'plugin1', out=>'missing2ref.out', cmd=>'+setGT --no-version',
//! args=>'-- -t . -n 0'`: compare stdout (after `grep -v ^##bcftools_`)
//! byte-for-byte against the fixture.

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
fn setgt_missing_to_ref_matches_upstream_fixture() {
    ensure_binary_built();
    let input = fixture_path("plugin1.vcf");
    let expected = std::fs::read_to_string(fixture_path("missing2ref.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+setGT",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-t",
            ".",
            "-n",
            "0",
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "+setGT failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected);
}
