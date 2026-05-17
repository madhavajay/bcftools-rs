//! End-to-end parity tests for `+fill-from-fasta` against the upstream
//! `aa.out` / `ref.out` / `aa.2.out` fixtures (test.pl rows 738-740).

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

fn check(in_name: &str, args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(&format!("{in_name}.vcf"));
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec![
        "+fill-from-fasta".to_string(),
        input.to_str().unwrap().to_string(),
        "--".to_string(),
    ];
    for a in args {
        if *a == "{FA}" {
            full.push("-f".to_string());
            full.push(fixture_path("placeholder").to_str().unwrap().to_string());
        } else {
            full.push(a.to_string());
        }
    }
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
fn fill_from_fasta_ref() {
    let fa = fixture_path("norm.fa");
    check("ref", &["-f", fa.to_str().unwrap(), "-c", "REF"], "ref.out");
}

#[test]
fn fill_from_fasta_ref_replace_n() {
    let fa = fixture_path("aa.fa");
    check(
        "aa",
        &["-f", fa.to_str().unwrap(), "-c", "REF", "-N"],
        "aa.2.out",
    );
}

#[test]
fn fill_from_fasta_info_filter_matches_upstream_fixture() {
    let fa = fixture_path("aa.fa");
    let hdr = fixture_path("aa.hdr");
    check(
        "aa",
        &[
            "-f",
            fa.to_str().unwrap(),
            "-c",
            "AA",
            "-h",
            hdr.to_str().unwrap(),
            "-i",
            "TYPE=\"snp\"",
        ],
        "aa.out",
    );
}
