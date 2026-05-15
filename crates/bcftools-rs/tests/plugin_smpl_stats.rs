//! End-to-end parity test for `+smpl-stats` against the upstream
//! `smpl-stats.1.out` fixture (input `indel-stats.vcf`, harness pipes
//! through `grep -v ^CMD`).

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

#[test]
fn smpl_stats_matches_upstream_fixture() {
    ensure_binary_built();
    let input = fixture_path("indel-stats.vcf");
    let expected = std::fs::read_to_string(fixture_path("smpl-stats.1.out")).unwrap();

    let out = Command::new(bin_path())
        .args(["+smpl-stats", input.to_str().unwrap()])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    // Upstream harness: `... +smpl-stats | grep -v ^CMD`.
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("CMD"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected);
}
