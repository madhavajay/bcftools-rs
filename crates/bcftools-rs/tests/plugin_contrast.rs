//! End-to-end parity tests for `+contrast` against the upstream
//! `contrast*.out` fixtures (test.pl rows 750-754).

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

fn check(input: &str, args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let inp = fixture_path(input);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+contrast", inp.to_str().unwrap()];
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
fn contrast_passoc_fassoc_novel_list() {
    check(
        "contrast.vcf",
        &[
            "-a",
            "PASSOC,FASSOC,NOVELAL,NOVELGT",
            "-0",
            "a,b",
            "-1",
            "c",
        ],
        "contrast.out",
    );
}

#[test]
fn contrast_passoc_fassoc_novel_files() {
    let c0 = fixture_path("contrast0.txt");
    let c1 = fixture_path("contrast1.txt");
    check(
        "contrast.vcf",
        &[
            "-a",
            "PASSOC,FASSOC,NOVELAL,NOVELGT",
            "-0",
            c0.to_str().unwrap(),
            "-1",
            c1.to_str().unwrap(),
        ],
        "contrast.out",
    );
}

#[test]
fn contrast_nassoc_force_samples() {
    check(
        "contrast.vcf",
        &["-a", "NASSOC", "-0", "a,b,c", "-1", "d", "--force-samples"],
        "contrast.1.out",
    );
}

#[test]
fn contrast_novelal_novelgt() {
    check(
        "contrast.1.vcf",
        &["-a", "NOVELAL,NOVELGT", "-0", "A", "-1", "B"],
        "contrast.1.1.out",
    );
}

#[test]
fn contrast_novelgt_only() {
    check(
        "contrast.1.vcf",
        &["-a", "NOVELGT", "-0", "A", "-1", "B"],
        "contrast.1.2.out",
    );
}

#[test]
fn contrast_rare_allele_summary_to_stderr() {
    ensure_binary_built();
    let input = fixture_path("contrast.vcf");
    let out = Command::new(bin_path())
        .args([
            "+contrast",
            input.to_str().unwrap(),
            "-a",
            "NASSOC",
            "-0",
            "a,b",
            "-1",
            "c",
            "-f",
            "1",
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert!(
        stderr.contains(
            "max_AC/PASSOC/FASSOC/NASSOC:\t1\t9.803922e-02\t0.000000,0.333333\t12,0,4,2\n"
        )
    );
}
