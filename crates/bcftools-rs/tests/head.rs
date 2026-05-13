//! End-to-end tests for `bcftools_rs::commands::head`.
//!
//! Tests build the `bcftools` binary from `crates/bcftools-rs-cli` and
//! invoke it as a subprocess to validate output. This mirrors how the
//! upstream Perl parity gate exercises the CLI and avoids fighting the
//! cargo test harness over stdout.

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
    // Per `cargo test` convention, sibling binaries live in the same target dir.
    // CARGO_BIN_EXE_<name> is set when the test crate depends (dev or build) on
    // the binary's package; we don't, so derive the path manually from
    // CARGO_MANIFEST_DIR.
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("target");
    p.push("debug");
    p.push("bcftools");
    p
}

fn ensure_binary_built() {
    let p = bin_path();
    if !p.exists() {
        let status = Command::new(env!("CARGO"))
            .args(["build", "-p", "bcftools-rs-cli"])
            .status()
            .expect("cargo build");
        assert!(status.success(), "failed to build bcftools-rs-cli");
        assert!(p.exists(), "binary not at expected path: {}", p.display());
    }
}

fn run(args: &[&str]) -> (String, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

#[test]
fn head_default_prints_full_header_no_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(out.starts_with("##fileformat=VCFv"), "got: {out:?}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(
        record_lines.is_empty(),
        "unexpected records: {record_lines:?}"
    );
}

#[test]
fn head_with_n2_emits_two_records_after_header() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", "-n", "2", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(record_lines.len(), 2);
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn head_with_h_truncates_header_lines() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", "-h", "3", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let header_lines: Vec<_> = out.lines().collect();
    assert_eq!(header_lines.len(), 3);
    assert!(header_lines[0].starts_with("##fileformat="));
}

#[test]
fn head_with_s_emits_chrom_line_then_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", "-s", "1", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let lines: Vec<_> = out.lines().collect();
    assert!(lines.iter().any(|l| l.starts_with("#CHROM\t")));
    let record_lines: Vec<_> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(record_lines.len(), 1);
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn version_flag_prints_block() {
    let (out, _err, code) = run(&["--version"]);
    assert_eq!(code, 0);
    assert!(out.contains("bcftools "));
    assert!(out.contains("htslib "));
}

#[test]
fn version_only_one_line() {
    let (out, _err, code) = run(&["--version-only"]);
    assert_eq!(code, 0);
    assert!(out.contains("+htslib-"));
    assert_eq!(out.lines().count(), 1);
}

#[test]
fn unknown_subcommand_errors() {
    let (_out, err, code) = run(&["bogus"]);
    assert_ne!(code, 0);
    assert!(err.contains("unrecognized command 'bogus'"));
}
