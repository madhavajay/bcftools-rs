//! End-to-end tests for the `+allele-length` plugin.

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

fn run(args: &[&str]) -> (String, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stdout = String::from_utf8(out.stdout).unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

#[test]
fn allele_length_matches_upstream_fixture() {
    let input = fixture_path("query.nucleotide.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.allele-length.tsv")).unwrap();

    let (out, err, code) = run(&["+allele-length", input.to_str().unwrap()]);
    assert_eq!(code, 0, "+allele-length failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn allele_length_via_plugin_subcommand() {
    let input = fixture_path("query.nucleotide.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.allele-length.tsv")).unwrap();

    let (out, err, code) = run(&["plugin", "allele-length", input.to_str().unwrap()]);
    assert_eq!(code, 0, "plugin allele-length failed: {err}");
    assert_eq!(out, expected);
}
