//! End-to-end test for the `+add-variantkey` plugin against the upstream
//! `query.add-variantkey.vcf` fixture (66 records, 3 hash/non-reversible).

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
fn add_variantkey_matches_upstream_fixture() {
    ensure_binary_built();
    let input = fixture_path("query.variantkey.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.add-variantkey.vcf")).unwrap();

    let out = Command::new(bin_path())
        .args(["+add-variantkey", input.to_str().unwrap()])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // The upstream harness pipes plugin output through `grep -v ^##bcftools_`
    // before diffing; the checked-in fixture is post-filter.
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();

    assert_eq!(filtered, expected);
}
