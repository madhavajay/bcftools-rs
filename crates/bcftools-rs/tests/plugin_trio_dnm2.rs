//! End-to-end parity for the `+trio-dnm2` `--use-NAIVE` slice
//! (upstream test.pl rows 768-769): `+trio-dnm2 -p [1X:|2X:]P,F,M
//! --use-NAIVE` piped through our own `bcftools query`.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

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

/// `+trio-dnm2 <in> -- -p <pfm> --use-NAIVE | bcftools query -f<fmt>`.
fn check(pfm: &str, fmt: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("trio-dnm/trio-dnm.9.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let plugin = Command::new(bin_path())
        .args([
            "+trio-dnm2",
            input.to_str().unwrap(),
            "--",
            "-p",
            pfm,
            "--use-NAIVE",
        ])
        .output()
        .expect("spawn +trio-dnm2");
    assert_eq!(
        plugin.status.code(),
        Some(0),
        "+trio-dnm2 -p {pfm} failed: {}",
        String::from_utf8_lossy(&plugin.stderr)
    );

    let mut q = Command::new(bin_path())
        .args(["query", "-f", fmt, "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("spawn query");
    q.stdin
        .take()
        .unwrap()
        .write_all(&plugin.stdout)
        .expect("pipe to query");
    let out = q.wait_with_output().expect("query output");
    assert_eq!(
        out.status.code(),
        Some(0),
        "query failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        expected,
        "mismatch for -p {pfm}"
    );
}

#[test]
fn naive_male_proband_chrx() {
    check(
        "1X:proband,father,mother",
        "[\t%DNM]\n",
        "trio-dnm/trio-dnm.9.1.out",
    );
}

#[test]
fn naive_female_proband_chrxx() {
    check(
        "2X:proband,father,mother",
        "[\t%DNM]\n",
        "trio-dnm/trio-dnm.9.2.out",
    );
}
