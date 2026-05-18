//! End-to-end test for the `+vcf2table` plugin (first slice).
//!
//! Mirrors the upstream `test_vcf_plugin` row
//! `in=>'merge.4.b', out=>'vcf2table.1.out', cmd=>'+vcf2table'`: run the
//! plugin and compare stdout (after `grep -v ^##bcftools_`, a no-op for
//! this table output) byte-for-byte against the checked-in fixture.

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
fn vcf2table_matches_upstream_fixture() {
    ensure_binary_built();
    let input = fixture_path("merge.4.b.vcf");
    let expected = std::fs::read_to_string(fixture_path("vcf2table.1.out")).unwrap();

    let out = Command::new(bin_path())
        .args(["+vcf2table", input.to_str().unwrap()])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "+vcf2table failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Upstream pipes through `grep -v ^##bcftools_`; the table output
    // contains no such lines, so filtering is a no-op kept for fidelity.
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected);
}

#[test]
fn vcf2table_vep_bcsq_matches_upstream_fixture() {
    // Upstream row: in=>'split-vep.2', out=>'vcf2table.2.out',
    // cmd=>'+vcf2table', args=>'-- --hide INFO,URL'.
    ensure_binary_built();
    let input = fixture_path("split-vep.2.vcf");
    let expected = std::fs::read_to_string(fixture_path("vcf2table.2.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+vcf2table",
            input.to_str().unwrap(),
            "--",
            "--hide",
            "INFO,URL",
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "+vcf2table failed: {}",
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
