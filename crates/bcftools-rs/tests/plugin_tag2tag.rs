//! End-to-end tests for the `+tag2tag` plugin (gl-to-pl, gl-to-gp, gp-to-gt).

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
fn tag2tag_gl_to_pl_matches_upstream_fixture() {
    let input = fixture_path("view.GL.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.PL.vcf")).unwrap();
    let (out, err, code) = run(&[
        "+tag2tag",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "-r",
        "--gl-to-pl",
    ]);
    assert_eq!(code, 0, "+tag2tag --gl-to-pl failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn tag2tag_gp_to_gt_matches_upstream_fixture() {
    let input = fixture_path("view.GP.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.GT.vcf")).unwrap();
    let (out, err, code) = run(&[
        "+tag2tag",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "-r",
        "--gp-to-gt",
        "-t",
        "0.2",
    ]);
    assert_eq!(code, 0, "+tag2tag --gp-to-gt failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn tag2tag_gl_to_gp_matches_upstream_fixture() {
    let input = fixture_path("view.GL.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.GL-GP.vcf")).unwrap();
    let (out, err, code) = run(&[
        "+tag2tag",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "--gl-to-gp",
    ]);
    assert_eq!(code, 0, "+tag2tag --gl-to-gp failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn tag2tag_via_plugin_subcommand() {
    let input = fixture_path("view.GL.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.PL.vcf")).unwrap();
    let (out, err, code) = run(&[
        "plugin",
        "tag2tag",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "-r",
        "--gl-to-pl",
    ]);
    assert_eq!(code, 0, "plugin tag2tag failed: {err}");
    assert_eq!(out, expected);
}
