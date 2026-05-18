//! End-to-end parity for the `+split-vep` first slice: the
//! `-c FIELD -s TR:CSQ[:PRN]` path that annotates `INFO/<FIELD>` from the
//! VEP `CSQ` string, validated by piping through our own `bcftools query`
//! (upstream test.pl rows 786/787/788/790).

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

/// `+split-vep <in> -c Consequence -s <sel>` piped through
/// `bcftools query -f '%POS\t%Consequence\n' [query_args]`.
fn check(sel: &str, query_args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("split-vep.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let plugin = Command::new(bin_path())
        .args([
            "+split-vep",
            input.to_str().unwrap(),
            "-c",
            "Consequence",
            "-s",
            sel,
        ])
        .output()
        .expect("spawn +split-vep");
    assert_eq!(
        plugin.status.code(),
        Some(0),
        "+split-vep -s {sel} failed: {}",
        String::from_utf8_lossy(&plugin.stderr)
    );

    let mut q = Command::new(bin_path())
        .arg("query")
        .arg("-f")
        .arg("%POS\t%Consequence\n")
        .args(query_args)
        .arg("-")
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

    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for -s {sel} {query_args:?}");
}

#[test]
fn worst_missense_threshold() {
    check("worst:missense+", &[], "split-vep.1.out");
}

#[test]
fn worst_missense_threshold_prn_worst() {
    check("worst:missense+:worst", &[], "split-vep.1.1.out");
}

#[test]
fn worst_missense_threshold_query_filter() {
    check(
        "worst:missense+",
        &["-i", "Consequence!=\".\""],
        "split-vep.2.out",
    );
}

#[test]
fn worst_missense_threshold_prn_worst_query_filter() {
    check(
        "worst:missense+:worst",
        &["-i", "Consequence!=\".\""],
        "split-vep.2.1.out",
    );
}
