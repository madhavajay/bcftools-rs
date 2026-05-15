//! End-to-end tests for the `+variant-distance` plugin.

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

fn check(args: &[&str], expected_fixture: &str) {
    let input = fixture_path("variant-distance.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+variant-distance"];
    full.extend_from_slice(args);
    full.push(input.to_str().unwrap());
    let (out, err, code) = run(&full);
    assert_eq!(code, 0, "{full:?} failed: {err}");
    assert_eq!(out, expected, "mismatch for {full:?}");
}

#[test]
fn variant_distance_default_nearest() {
    check(&[], "variant-distance.1.out");
}

#[test]
fn variant_distance_explicit_nearest() {
    check(&["-d", "nearest"], "variant-distance.1.out");
}

#[test]
fn variant_distance_fwd() {
    check(&["-d", "fwd"], "variant-distance.2.out");
}

#[test]
fn variant_distance_rev() {
    check(&["-d", "rev"], "variant-distance.3.out");
}

#[test]
fn variant_distance_both() {
    check(&["-d", "both"], "variant-distance.4.out");
}

#[test]
fn variant_distance_via_plugin_subcommand() {
    let input = fixture_path("variant-distance.vcf");
    let expected = std::fs::read_to_string(fixture_path("variant-distance.1.out")).unwrap();
    let (out, err, code) = run(&["plugin", "variant-distance", input.to_str().unwrap()]);
    assert_eq!(code, 0, "plugin variant-distance failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn variant_distance_rejects_bad_direction() {
    let input = fixture_path("variant-distance.vcf");
    let (_o, err, code) = run(&[
        "+variant-distance",
        "-d",
        "sideways",
        input.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "expected failure for bad direction");
    assert!(err.contains("unknown -d direction"), "stderr: {err}");
}
