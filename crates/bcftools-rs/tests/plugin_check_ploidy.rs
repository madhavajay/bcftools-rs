//! End-to-end tests for the `+check-ploidy` plugin.

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
fn check_ploidy_fixture_1() {
    let input = fixture_path("checkploidy.vcf");
    let expected = std::fs::read_to_string(fixture_path("checkploidy.out")).unwrap();
    let (out, err, code) = run(&["+check-ploidy", "--no-version", input.to_str().unwrap()]);
    assert_eq!(code, 0, "+check-ploidy failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn check_ploidy_fixture_2_default() {
    let input = fixture_path("checkploidy.2.vcf");
    let expected = std::fs::read_to_string(fixture_path("checkploidy.2.out")).unwrap();
    let (out, err, code) = run(&["+check-ploidy", "--no-version", input.to_str().unwrap()]);
    assert_eq!(code, 0, "+check-ploidy 2 failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn check_ploidy_fixture_2_use_missing() {
    let input = fixture_path("checkploidy.2.vcf");
    let expected = std::fs::read_to_string(fixture_path("checkploidy.3.out")).unwrap();
    let (out, err, code) = run(&[
        "+check-ploidy",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "-m",
    ]);
    assert_eq!(code, 0, "+check-ploidy -m failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn check_ploidy_via_plugin_subcommand() {
    let input = fixture_path("checkploidy.vcf");
    let expected = std::fs::read_to_string(fixture_path("checkploidy.out")).unwrap();
    let (out, err, code) = run(&[
        "plugin",
        "check-ploidy",
        "--no-version",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "plugin check-ploidy failed: {err}");
    assert_eq!(out, expected);
}
