//! End-to-end parity tests for `+trio-stats` against the upstream
//! `trio-stats.out` / `trio-stats.2.out` fixtures (test.pl rows 721-722;
//! harness pipes `grep -v ^CMD`).

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

fn check(args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("trio-stats.vcf");
    let ped = fixture_path("trio-stats.ped");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec!["+trio-stats", input.to_str().unwrap(), "--"];
    let mut sub: Vec<String> = Vec::new();
    for a in args {
        sub.push(if *a == "{PED}" {
            ped.to_str().unwrap().to_string()
        } else {
            a.to_string()
        });
    }
    let subref: Vec<&str> = sub.iter().map(|s| s.as_str()).collect();
    full.extend_from_slice(&subref);

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
        .filter(|l| !l.starts_with("CMD"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn trio_stats_alt1() {
    check(
        &["-a", "1", "-p", "{PED}", "-d", "mendel-errors,transmitted"],
        "trio-stats.out",
    );
}

#[test]
fn trio_stats_no_alt() {
    check(
        &["-p", "{PED}", "-d", "mendel-errors,transmitted"],
        "trio-stats.2.out",
    );
}
