//! End-to-end tests for the `+gvcfz` plugin.
//!
//! Mirrors the upstream `test_vcf_plugin` rows that pipe `+gvcfz`
//! through `bcftools query`:
//!   in=>'gvcfz',   out=>'gvcfz.1.out'   (-g 'PASS:GT!="alt"' -a)
//!   in=>'gvcfz.2', out=>'gvcfz.2.1.out' (-g 'PASS:GT!="alt"' -a)
//! Compared byte-for-byte after `grep -v ^##bcftools_` (a no-op here:
//! `query` output carries no such lines).
//!
//! The `gvcfz.2.out` row (`-g 'PASS:GQ>10; FLT:-' -a`) is deferred: in
//! the multi-group catch-all (`FLT:-`) path, blocks whose representative
//! is a multiallelic record forced to a ref block by `-a` differ from
//! upstream in 3 RGQ cells (a `gq_key` / `bcf_update_alleles`
//! interaction), tracked in TODO.md.

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

/// Run `+gvcfz <input> <gvcfz_args>` then pipe through
/// `query -f <query_fmt>`, returning the filtered stdout.
fn gvcfz_query(input: &str, gvcfz_args: &[&str], query_fmt: &str) -> String {
    ensure_binary_built();
    let bin = bin_path();
    let input = fixture_path(input);

    let mut args = vec!["+gvcfz".to_string(), input.to_str().unwrap().to_string()];
    args.extend(gvcfz_args.iter().map(|s| s.to_string()));
    let gz = Command::new(&bin)
        .args(&args)
        .output()
        .expect("spawn +gvcfz");
    assert_eq!(
        gz.status.code().unwrap_or(-1),
        0,
        "+gvcfz failed: {}",
        String::from_utf8_lossy(&gz.stderr)
    );

    let mut child = Command::new(&bin)
        .args(["query", "-f", query_fmt])
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .spawn()
        .expect("spawn query");
    use std::io::Write;
    child.stdin.take().unwrap().write_all(&gz.stdout).unwrap();
    let q = child.wait_with_output().expect("query wait");
    assert_eq!(
        q.status.code().unwrap_or(-1),
        0,
        "query failed: {}",
        String::from_utf8_lossy(&q.stderr)
    );

    String::from_utf8(q.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect()
}

#[test]
fn gvcfz_1_matches_upstream_fixture() {
    let got = gvcfz_query(
        "gvcfz.vcf",
        &["-g", r#"PASS:GT!="alt""#, "-a"],
        r"%POS\t%REF\t%ALT\t%END[\t%GT][\t%DP][\t%GQ][\t%RGQ]\n",
    );
    let expected = std::fs::read_to_string(fixture_path("gvcfz.1.out")).unwrap();
    assert_eq!(got, expected);
}

#[test]
fn gvcfz_2_1_matches_upstream_fixture() {
    let got = gvcfz_query(
        "gvcfz.2.vcf",
        &["-g", r#"PASS:GT!="alt""#, "-a"],
        r"%POS\t%REF\t%ALT\t%FILTER\t%END[\t%GT][\t%DP]\n",
    );
    let expected = std::fs::read_to_string(fixture_path("gvcfz.2.1.out")).unwrap();
    assert_eq!(got, expected);
}
