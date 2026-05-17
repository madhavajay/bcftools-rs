//! End-to-end tests for `bcftools_rs::commands::norm`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

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

fn run(args: &[&str]) -> (Vec<u8>, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (out.stdout, stderr, out.status.code().unwrap_or(-1))
}

#[test]
fn norm_rmdup_snps_matches_upstream_fixture() {
    let input = fixture_path("norm.rmdup.vcf");
    let expected = std::fs::read_to_string(fixture_path("norm.rmdup.1.out")).unwrap();

    let (out, err, code) = run(&[
        "norm",
        "--no-version",
        "-d",
        "snps",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "norm -d snps failed: {err}");
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn norm_rmdup_exact_matches_upstream_fixture() {
    let input = fixture_path("norm.rmdup.vcf");
    let expected = std::fs::read_to_string(fixture_path("norm.rmdup.5.out")).unwrap();

    let (out, err, code) = run(&[
        "norm",
        "--no-version",
        "-d",
        "exact",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "norm -d exact failed: {err}");
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn norm_rmdup_indels_both_and_all_match_upstream_fixtures() {
    let input = fixture_path("norm.rmdup.vcf");
    for (mode, fixture) in [
        ("indels", "norm.rmdup.2.out"),
        ("both", "norm.rmdup.3.out"),
        ("all", "norm.rmdup.4.out"),
    ] {
        let expected = std::fs::read_to_string(fixture_path(fixture)).unwrap();
        let (out, err, code) = run(&["norm", "--no-version", "-d", mode, input.to_str().unwrap()]);
        assert_eq!(code, 0, "norm -d {mode} failed: {err}");
        assert_eq!(String::from_utf8(out).unwrap(), expected, "mode {mode}");
    }
}

#[test]
fn norm_rmdup2_numbered_fixtures_match_upstream_text_output() {
    let input = fixture_path("norm.rmdup.2.vcf");
    for (mode, fixture) in [
        ("none", "norm.rmdup.2.1.out"),
        ("exact", "norm.rmdup.2.1.out"),
        ("indels", "norm.rmdup.2.1.out"),
        ("any", "norm.rmdup.2.2.out"),
        ("both", "norm.rmdup.2.2.out"),
        ("snps", "norm.rmdup.2.2.out"),
    ] {
        let expected = std::fs::read_to_string(fixture_path(fixture)).unwrap();
        let (out, err, code) = run(&["norm", "--no-version", "-d", mode, input.to_str().unwrap()]);
        assert_eq!(code, 0, "norm.rmdup.2 fixture failed for -d {mode}: {err}");
        assert_eq!(
            String::from_utf8(out).unwrap(),
            expected,
            "norm.rmdup.2 fixture differed for -d {mode}"
        );
    }
}

#[test]
fn norm_rmdup3_left_aligned_duplicate_fixtures_match_upstream_text_output() {
    let input = fixture_path("norm.rmdup.3.vcf");
    let reference = fixture_path("norm.rmdup.3.fa");
    for (mode, fixture) in [
        ("exact", "norm.rmdup.3.1.out"),
        ("all", "norm.rmdup.3.2.out"),
    ] {
        let expected = std::fs::read_to_string(fixture_path(fixture)).unwrap();
        let (out, err, code) = run(&[
            "norm",
            "--no-version",
            "-d",
            mode,
            "-f",
            reference.to_str().unwrap(),
            input.to_str().unwrap(),
        ]);
        assert_eq!(code, 0, "norm.rmdup.3 fixture failed for -d {mode}: {err}");
        assert_eq!(
            String::from_utf8(out).unwrap(),
            expected,
            "norm.rmdup.3 fixture differed for -d {mode}"
        );
    }
}

#[test]
fn norm_filter_duplicate_removal_fixtures_match_upstream_text_output() {
    let input = fixture_path("norm.filter.vcf");
    let filter_file = fixture_path("norm.filter.txt");
    let expected = std::fs::read_to_string(fixture_path("norm.filter.2.out")).unwrap();

    let membership_expr = format!("ID=@{}", filter_file.display());
    for include_expr in [membership_expr.as_str(), "ALT!=\"C\""] {
        let args = [
            "norm",
            "--no-version",
            "-d",
            "both",
            "-i",
            include_expr,
            input.to_str().unwrap(),
        ];
        let (out, err, code) = run(&args);
        assert_eq!(
            code, 0,
            "norm.filter fixture failed for -i {include_expr}: {err}"
        );
        assert_eq!(
            String::from_utf8(out).unwrap(),
            expected,
            "norm.filter fixture differed for -i {include_expr}"
        );
    }
}

#[test]
fn norm_check_ref_swap_matches_upstream_text_output() {
    let input = fixture_path("norm.check-ref.vcf");
    let reference = fixture_path("norm.check-ref.fa");
    let expected = std::fs::read_to_string(fixture_path("norm.check-ref.1.out")).unwrap();

    let (out, err, code) = run(&[
        "norm",
        "--no-version",
        "-c",
        "s",
        "-f",
        reference.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "norm -c s failed: {err}");
    assert_eq!(String::from_utf8(out).unwrap(), expected);
}

#[test]
fn norm_sort_split_multiallelic_fixtures_match_upstream_text_output() {
    let input = fixture_path("norm.sort.vcf");
    for (args, fixture) in [
        (
            vec!["norm", "--no-version", "-m", "-", input.to_str().unwrap()],
            "norm.sort.1.out",
        ),
        (
            vec![
                "norm",
                "--no-version",
                "-m",
                "-",
                "-S",
                "lex",
                input.to_str().unwrap(),
            ],
            "norm.sort.2.out",
        ),
    ] {
        let expected = std::fs::read_to_string(fixture_path(fixture)).unwrap();
        let (out, err, code) = run(&args);
        assert_eq!(code, 0, "norm split fixture {fixture} failed: {err}");
        assert_eq!(
            String::from_utf8(out).unwrap(),
            expected,
            "norm split fixture {fixture} differed"
        );
    }
}

#[test]
fn norm_rmdup_reads_bcf_and_writes_bcf() {
    let dir = TempDir::new().unwrap();
    let input = fixture_path("norm.rmdup.vcf");
    let bcf = dir.path().join("in.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (bcf_out, err, code) = run(&[
        "norm",
        "--no-version",
        "-d",
        "all",
        "-Ou",
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "norm -Ou failed: {err}");

    ensure_binary_built();
    let mut child = Command::new(bin_path())
        .args(["query", "-f%CHROM\\t%POS\\t%REF\\t%ALT\\n", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn query");
    child.stdin.as_mut().unwrap().write_all(&bcf_out).unwrap();
    let out = child.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "query failed: {stderr}");
    let text = String::from_utf8(out.stdout).unwrap();
    assert!(text.contains("1\t789241\tC\tG\n"), "{text}");
    assert!(!text.contains("1\t789245\tC\tT\n"), "{text}");
}
