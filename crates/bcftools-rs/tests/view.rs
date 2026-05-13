//! End-to-end tests for `bcftools_rs::commands::view`.

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
    let p = bin_path();
    if !p.exists() {
        let status = Command::new(env!("CARGO"))
            .args(["build", "-p", "bcftools-rs-cli"])
            .status()
            .expect("cargo build");
        assert!(status.success(), "failed to build bcftools-rs-cli");
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
fn view_text_round_trip_emits_all_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    // The fixture has 22 records (lines 11..=32 of the file).
    assert_eq!(record_lines.len(), 21);
    // First record is `1\t105\t.\tTAAACCCTA\t...`
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn view_header_only_drops_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(record_lines.is_empty());
    assert!(out.contains("#CHROM\t"));
}

#[test]
fn view_no_header_drops_header() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-H", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let header_lines: Vec<_> = out.lines().filter(|l| l.starts_with('#')).collect();
    assert!(header_lines.is_empty(), "header leaked: {header_lines:?}");
    let record_lines: Vec<_> = out.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(record_lines.len(), 21);
}

#[test]
fn view_unknown_output_type_errors() {
    let path = fixture_path("aa.vcf");
    let (_out, err, code) = run(&["view", "-Oq", path.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("not recognised"), "stderr: {err}");
}

#[test]
fn view_default_injects_bcftools_version_and_command_lines() {
    // Without --no-version, `bcftools view` must emit
    // ##bcftools_viewVersion=<v> and ##bcftools_viewCommand=<cmdline>; Date=...
    // header lines, mirroring upstream `bcf_hdr_append_version`.
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);

    let header_lines: Vec<_> = out.lines().filter(|l| l.starts_with("##")).collect();
    let version_line = header_lines
        .iter()
        .find(|l| l.starts_with("##bcftools_viewVersion="))
        .unwrap_or_else(|| panic!("missing version line in header:\n{out}"));
    let command_line = header_lines
        .iter()
        .find(|l| l.starts_with("##bcftools_viewCommand="))
        .unwrap_or_else(|| panic!("missing command line in header:\n{out}"));

    // The version line ends with the htslib-rs version we're built against.
    assert!(version_line.contains("+htslib-"), "got: {version_line}");
    // The command line names the subcommand and includes a `Date=` field.
    assert!(command_line.contains("view"), "got: {command_line}");
    assert!(command_line.contains("; Date="), "got: {command_line}");
}

#[test]
fn view_no_version_suppresses_injected_header_lines() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(
        !out.contains("##bcftools_viewVersion="),
        "version line leaked despite --no-version:\n{out}"
    );
    assert!(
        !out.contains("##bcftools_viewCommand="),
        "command line leaked despite --no-version:\n{out}"
    );
}
