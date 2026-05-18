//! End-to-end parity tests for the `+remove-overlaps` plugin against the
//! upstream `remove-overlaps.1.*` fixtures (overlap/dup modes).

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

fn check(args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("remove-overlaps.1.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec!["+remove-overlaps", input.to_str().unwrap()];
    full.extend_from_slice(args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // The harness compares plugin output after `grep -v ^##bcftools_`.
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

fn check_in(input_fixture: &str, args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input_fixture);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+remove-overlaps", input.to_str().unwrap(), "--"];
    full.extend_from_slice(args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn min_qual_overlap_marking() {
    // Upstream rows: in=>'remove-overlaps.2'→2.1.out, 'remove-overlaps.3'
    // →3.1.out (also with --missing 0), `-m 'min(QUAL)' -M rmme`.
    check_in(
        "remove-overlaps.2.vcf",
        &["-m", "min(QUAL)", "-M", "rmme"],
        "remove-overlaps.2.1.out",
    );
    check_in(
        "remove-overlaps.3.vcf",
        &["-m", "min(QUAL)", "-M", "rmme"],
        "remove-overlaps.3.1.out",
    );
    check_in(
        "remove-overlaps.3.vcf",
        &["-m", "min(QUAL)", "-M", "rmme", "--missing", "0"],
        "remove-overlaps.3.1.out",
    );
}

#[test]
fn overlap_remove() {
    check(&["-m", "overlap"], "remove-overlaps.1.1.out");
}

#[test]
fn overlap_mark() {
    check(
        &["-m", "overlap", "-M", "overlap"],
        "remove-overlaps.1.2.out",
    );
}

#[test]
fn overlap_text_list() {
    check(&["-m", "overlap", "-O", "t"], "remove-overlaps.1.3.out");
}

#[test]
fn overlap_reverse() {
    check(&["-m", "overlap", "--reverse"], "remove-overlaps.1.4.out");
}

#[test]
fn dup_mark() {
    check(&["-m", "dup", "-M", "DUP"], "remove-overlaps.1.5.out");
}

#[test]
fn dup_mark_unique_reverse() {
    check(
        &["-m", "dup", "-M", "unique", "--reverse"],
        "remove-overlaps.1.6.out",
    );
}
