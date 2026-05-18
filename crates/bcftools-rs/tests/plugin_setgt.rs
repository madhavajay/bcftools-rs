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

#[test]
fn setgt_query_per_sample_filter_matches_upstream_fixture() {
    // Upstream row: in=>'setGT', out=>'setGT.1.out',
    // args=>'-- -t q -n 0 -i \'GT~"." && FMT/DP=30 && GQ=150\''.
    ensure_binary_built();
    let input = fixture_path("setGT.vcf");
    let expected = std::fs::read_to_string(fixture_path("setGT.1.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+setGT",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-t",
            "q",
            "-n",
            "0",
            "-i",
            r#"GT~"." && FMT/DP=30 && GQ=150"#,
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "+setGT -t q failed: {}",
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

fn run_setgt_subset(out_fixture: &str, expr: &str) {
    ensure_binary_built();
    let input = fixture_path("setGT.2.vcf");
    let samples = fixture_path("setGT.samples.txt");
    let expr = expr.replace("{S}", samples.to_str().unwrap());
    let expected = std::fs::read_to_string(fixture_path(out_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+setGT",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-t",
            "q",
            "-n",
            ".",
            "-i",
            &expr,
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "+setGT -t q [@file] failed: {}",
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

#[test]
fn setgt_query_sample_subset_het_matches_upstream_fixture() {
    // in=>'setGT.2', out=>'setGT.2.out',
    // args=>'-- -t q -n . -i \'GT[@.../setGT.samples.txt]="het"\''.
    run_setgt_subset("setGT.2.out", r#"GT[@{S}]="het""#);
}

#[test]
fn setgt_query_sample_subset_het_binom_matches_upstream_fixture() {
    // in=>'setGT.2', out=>'setGT.3.out', adds & binom(AD[@file])<0.1.
    run_setgt_subset("setGT.3.out", r#"GT[@{S}]="het" & binom(AD[@{S}])<0.1"#);
}

#[test]
fn setgt_invert_phase_matches_upstream_fixture() {
    // Upstream row: in=>'setGT.2', out=>'setGT.2.1.out',
    // args=>'-- -t a -n i' (invert allele order, separator preserved).
    ensure_binary_built();
    let input = fixture_path("setGT.2.vcf");
    let expected = std::fs::read_to_string(fixture_path("setGT.2.1.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+setGT",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-t",
            "a",
            "-n",
            "i",
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "+setGT -t a -n i failed: {}",
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
