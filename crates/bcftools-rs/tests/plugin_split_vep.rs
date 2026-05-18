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

/// `+split-vep <in> [extra args] -f <fmt>` — the text-output path.
/// Upstream `test_vcf_plugin` compares after `grep -v ^##bcftools_`,
/// which is a no-op for `-f` rendered text.
fn check_f(input: &str, extra: &[&str], fmt: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+split-vep", input.to_str().unwrap()];
    full.extend_from_slice(extra);
    let farg = format!("-f{fmt}");
    full.push(&farg);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn +split-vep -f");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        expected,
        "mismatch for {full:?}"
    );
}

#[test]
fn format_worst_missense() {
    check_f(
        "split-vep.vcf",
        &["-s", "worst:missense+"],
        r"%POS\t%Consequence\n",
        "split-vep.2.out",
    );
    check_f(
        "split-vep.vcf",
        &["-s", "worst:missense+:worst"],
        r"%POS\t%Consequence\n",
        "split-vep.2.1.out",
    );
}

#[test]
fn format_primary_missense() {
    check_f(
        "split-vep.vcf",
        &["-s", "primary:missense+"],
        r"%POS\t%Consequence\n",
        "split-vep.3.out",
    );
    check_f(
        "split-vep.vcf",
        &["-s", "primary:missense+:worst"],
        r"%POS\t%Consequence\n",
        "split-vep.3.1.out",
    );
    // No CSQ field referenced: nannot==0, severity-pass gate only.
    check_f(
        "split-vep.vcf",
        &["-s", "primary:missense+"],
        r"%POS\n",
        "split-vep.4.out",
    );
}

#[test]
fn format_csq_subfield_vs_real_info() {
    // `%AF` resolves to the CSQ subfield (rows with empty AF dropped).
    check_f(
        "split-vep.2.vcf",
        &["-s", "worst"],
        r"%POS\t%AF\n",
        "split-vep.5.out",
    );
    // `-a BCSQ` has no `AF` subfield → `%AF` falls back to real INFO/AF.
    check_f(
        "split-vep.2.vcf",
        &["-s", "worst", "-a", "BCSQ"],
        r"%POS\t%AF\n",
        "split-vep.6.out",
    );
    // `%INFO/AF` is always the real INFO tag, never the CSQ subfield.
    check_f(
        "split-vep.2.vcf",
        &["-s", "worst"],
        r"%POS\t%INFO/AF\n",
        "split-vep.6.out",
    );
    check_f(
        "split-vep.2.vcf",
        &["-s", "worst", "-a", "BCSQ"],
        r"%POS\t%INFO/AF\n",
        "split-vep.6.out",
    );
}

#[test]
fn format_worst_unfiltered() {
    check_f(
        "split-vep.3.vcf",
        &["-s", "worst"],
        r"%POS\t%Consequence\n",
        "split-vep.7.out",
    );
    check_f(
        "split-vep.3.vcf",
        &["-s", "worst::worst"],
        r"%POS\t%Consequence\n",
        "split-vep.7.1.out",
    );
}

#[test]
fn format_region_and_duplicate() {
    // `-t 1:14464`: only that site; default select → all transcripts,
    // annotations comma-joined across them.
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464"],
        r"%POS\t%CANONICAL\t%Consequence\n",
        "split-vep.9.out",
    );
    // `-d`: one output row per selected severity-passing transcript.
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464", "-d"],
        r"%POS\t%CANONICAL\t%Consequence\n",
        "split-vep.10.out",
    );
}

#[test]
fn format_raw_tag_expansion_and_headers() {
    // `%CSQ` + `-A tab`: expand to every subfield (comma-joined across
    // transcripts when not `-d`).
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464", "-A", "tab"],
        r"%POS\t%CSQ\n",
        "split-vep.11.out",
    );
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464", "-A", "tab", "-d"],
        r"%POS\t%CSQ\n",
        "split-vep.12.out",
    );
    // `-H` / `-HH` header rows.
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464", "-A", "tab", "-d", "-H"],
        r"%POS\t%CSQ\n",
        "split-vep.12.2.out",
    );
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464", "-A", "tab", "-d", "-HH"],
        r"%POS\t%CSQ\n",
        "split-vep.12.3.out",
    );
    // Custom `-A` delimiter (the POS→block separator stays the format's
    // own TAB; the delimiter joins the expanded subfields).
    check_f(
        "split-vep.vcf",
        &["-t", "1:14464", "-A", "@@@", "-d", "-H"],
        r"%POS\t%CSQ\n",
        "split-vep.12.4.out",
    );
}

#[test]
fn format_bcsq_tag_autodetect() {
    // `-a BCSQ`, and tag auto-detection (no `-a`: CSQ absent → BCSQ).
    check_f(
        "split-vep.4.vcf",
        &["-a", "BCSQ", "-A", "tab", "-d"],
        r"%POS\t%BCSQ\n",
        "split-vep.13.out",
    );
    check_f(
        "split-vep.4.vcf",
        &["-A", "tab", "-d"],
        r"%POS\t%BCSQ\n",
        "split-vep.13.out",
    );
    // `-s ::worst` collapses the `&`-joined Consequence subfield.
    check_f(
        "split-vep.4.vcf",
        &["-A", "tab", "-d", "-s", "::worst"],
        r"%POS\t%BCSQ\n",
        "split-vep.13.1.out",
    );
}

#[test]
fn format_per_sample_block_with_filter() {
    // Per-sample `[%SAMPLE]` block + `-i 'GT="alt"'` — the split-vep
    // filter is handed to the query engine, which also drives per-sample
    // inclusion.
    check_f(
        "split-vep.3.vcf",
        &["-s", "worst", "-i", r#"GT="alt""#],
        r"[%POS\t%SAMPLE\t%GT\t%Consequence\n]",
        "split-vep.8.out",
    );
}
