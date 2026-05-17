//! End-to-end parity test for `+fixploidy` against the upstream
//! `fixploidy.out` fixture (test.pl row 673).

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
fn fixploidy_samples_and_ploidy() {
    ensure_binary_built();
    let input = fixture_path("fixploidy.vcf");
    let samples = fixture_path("fixploidy.samples");
    let ploidy = fixture_path("fixploidy.ploidy");
    let expected = std::fs::read_to_string(fixture_path("fixploidy.out")).unwrap();

    let out = Command::new(bin_path())
        .args([
            "+fixploidy",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-s",
            samples.to_str().unwrap(),
            "-p",
            ploidy.to_str().unwrap(),
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "fixploidy failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for fixploidy");
}

#[test]
fn fixploidy_honors_default_ploidy_lines() {
    ensure_binary_built();
    let tmp = TempDir::new().expect("tempdir");
    let input = tmp.path().join("input.vcf");
    let samples = tmp.path().join("samples.txt");
    let ploidy = tmp.path().join("ploidy.txt");

    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tM\tF\tU\n\
1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/1\t0/1\n",
    )
    .unwrap();
    std::fs::write(&samples, "M M\nF F\nU U\n").unwrap();
    std::fs::write(&ploidy, "* * * * 3\n* * * M 1\n* * * F 2\n").unwrap();

    let out = Command::new(bin_path())
        .args([
            "+fixploidy",
            "--no-version",
            input.to_str().unwrap(),
            "--",
            "-d",
            "4",
            "-s",
            samples.to_str().unwrap(),
            "-p",
            ploidy.to_str().unwrap(),
        ])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "fixploidy failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(
        stdout.contains("\n1\t10\t.\tA\tC\t.\t.\t.\tGT\t0\t0/1\t0/1/1\n"),
        "{stdout}"
    );
}
