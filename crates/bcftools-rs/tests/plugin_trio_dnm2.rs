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

/// `+trio-dnm2 <in> -- <plugin_args> | bcftools query -f<fmt>`.
fn check(input_fixture: &str, plugin_args: &[&str], fmt: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input_fixture);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut args: Vec<&str> = vec!["+trio-dnm2", input.to_str().unwrap(), "--"];
    args.extend_from_slice(plugin_args);
    let plugin = Command::new(bin_path())
        .args(&args)
        .output()
        .expect("spawn +trio-dnm2");
    assert_eq!(
        plugin.status.code(),
        Some(0),
        "{args:?} failed: {}",
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
        "mismatch for {args:?}"
    );
}

#[test]
fn naive_male_proband_chrx() {
    check(
        "trio-dnm/trio-dnm.9.vcf",
        &["-p", "1X:proband,father,mother", "--use-NAIVE"],
        "[\t%DNM]\n",
        "trio-dnm/trio-dnm.9.1.out",
    );
}

#[test]
fn naive_female_proband_chrxx() {
    check(
        "trio-dnm/trio-dnm.9.vcf",
        &["-p", "2X:proband,father,mother", "--use-NAIVE"],
        "[\t%DNM]\n",
        "trio-dnm/trio-dnm.9.2.out",
    );
}

#[test]
fn acm_default_model() {
    // The default ACM likelihood model (DNM:log), validated through
    // our own `bcftools query` (test.pl rows 758, 760, 762, 766).
    check(
        "trio-dnm/trio-dnm.4.vcf",
        &["-p", "proband,father,mother"],
        "[\t%DNM]\t[\t%VAF]\n",
        "trio-dnm/trio-dnm.4.1.out",
    );
    check(
        "trio-dnm/trio-dnm.4.vcf",
        &["-p", "proband,father,mother", "--dnm-tag", "DNM:log"],
        "[\t%DNM]\t[\t%VAF]\n",
        "trio-dnm/trio-dnm.4.2.out",
    );
    check(
        "trio-dnm/trio-dnm.5.vcf",
        &["-p", "proband,father,mother", "--dnm-tag", "DNM:log"],
        "[\t%DNM]\t[\t%VAF]\n",
        "trio-dnm/trio-dnm.5.1.out",
    );
    check(
        "trio-dnm/trio-dnm.7.vcf",
        &["-p", "proband,father,mother", "--dnm-tag", "DNM:log"],
        "[\t%DNM]\t[\t%VAF]\n",
        "trio-dnm/trio-dnm.7.1.out",
    );
}
