//! End-to-end parity test for `+af-dist` against the upstream `af-dist.out`
//! fixture (the harness pipes through `grep -v bcftools`).

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

#[test]
fn af_dist_matches_upstream_fixture() {
    ensure_binary_built();
    let input = fixture_path("af-dist.vcf");
    let expected = std::fs::read_to_string(fixture_path("af-dist.out")).unwrap();

    let out = Command::new(bin_path())
        .args(["+af-dist", input.to_str().unwrap()])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Upstream harness: `... +af-dist | grep -v bcftools`.
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.contains("bcftools"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected);
}

#[test]
fn af_dist_list_emits_matching_genotypes_before_histograms() {
    ensure_binary_built();
    let input = fixture_path("af-dist.vcf");

    let out = Command::new(bin_path())
        .args(["+af-dist", input.to_str().unwrap(), "--", "-l", "0.5,0.5"])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        stdout.contains("# GT, genotypes with P(AF) in [0.500000,0.500000];"),
        "{stdout}"
    );
    assert!(
        stdout.contains("GT\t20\t326891\tNA00001\t1\t0.500000\n"),
        "{stdout}"
    );
    assert!(
        stdout.contains("GT\t20\t326891\tNA00002\t1\t0.500000\n"),
        "{stdout}"
    );
    assert!(
        stdout
            .find("GT\t20\t326891\tNA00001\t1\t0.500000\n")
            .unwrap()
            < stdout.find("# PROB_DIST").unwrap(),
        "{stdout}"
    );
}

#[test]
fn af_dist_accepts_bin_definition_files() {
    ensure_binary_built();
    let input = fixture_path("af-dist.vcf");
    let tmp = TempDir::new().expect("tempdir");
    let bins = tmp.path().join("bins.txt");
    std::fs::write(&bins, "0.25\n0.5\n0.75\n").unwrap();

    let out = Command::new(bin_path())
        .args([
            "+af-dist",
            input.to_str().unwrap(),
            "--",
            "-d",
            bins.to_str().unwrap(),
            "-p",
            bins.to_str().unwrap(),
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        stdout.contains("PROB_DIST\t0.000000\t0.250000\t"),
        "{stdout}"
    );
    assert!(
        stdout.contains("PROB_DIST\t0.750000\t1.000000\t"),
        "{stdout}"
    );
    assert!(
        stdout.contains("DEV_DIST\t0.000000\t0.250000\t"),
        "{stdout}"
    );
    assert!(
        stdout.contains("DEV_DIST\t0.750000\t1.000000\t"),
        "{stdout}"
    );
}
