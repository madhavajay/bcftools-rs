//! End-to-end tests for the `+fill-AN-AC` plugin.

use std::path::PathBuf;
use std::process::Command;

use tempfile::TempDir;

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
fn fill_an_ac_matches_upstream_fixture() {
    let input = fixture_path("plugin1.vcf");
    let expected = std::fs::read_to_string(fixture_path("fill-AN-AC.out")).unwrap();

    let (out, err, code) = run(&["+fill-AN-AC", "--no-version", input.to_str().unwrap()]);
    assert_eq!(code, 0, "+fill-AN-AC failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn fill_an_ac_via_plugin_subcommand() {
    let input = fixture_path("plugin1.vcf");
    let expected = std::fs::read_to_string(fixture_path("fill-AN-AC.out")).unwrap();

    let (out, err, code) = run(&[
        "plugin",
        "fill-AN-AC",
        "--no-version",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "plugin fill-AN-AC failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn fill_an_ac_reads_bcf_input() {
    let dir = TempDir::new().unwrap();
    let input = fixture_path("plugin1.vcf");
    let bcf = dir.path().join("in.bcf");

    let (_o, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["+fill-AN-AC", "--no-version", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "+fill-AN-AC on bcf failed: {err}");
    assert!(out.contains(";AN=4;AC=2\t"), "{out}");
    assert!(out.contains("\tAN=0;AC=0\t"), "{out}");
}

#[test]
fn fill_an_ac_writes_bgzf_output() {
    let dir = TempDir::new().unwrap();
    let input = fixture_path("plugin1.vcf");
    let out_path = dir.path().join("out.vcf.gz");

    let (_o, err, code) = run(&[
        "+fill-AN-AC",
        "--no-version",
        "-Oz",
        "-o",
        out_path.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "+fill-AN-AC -Oz failed: {err}");
    let bytes = std::fs::read(&out_path).unwrap();
    assert!(bytes.starts_with(&[0x1f, 0x8b]), "expected gzip magic");
}
