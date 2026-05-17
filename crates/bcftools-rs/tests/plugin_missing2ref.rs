//! End-to-end tests for the `+missing2ref` plugin.

use std::io::Write;
use std::path::PathBuf;
use std::process::{Command, Stdio};

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

fn run(args: &[&str]) -> (String, String, i32) {
    ensure_binary_built();
    let out = Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools");
    let stdout = String::from_utf8(out.stdout).unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

#[test]
fn missing2ref_matches_upstream_fixture() {
    let input = fixture_path("plugin1.vcf");
    let expected = std::fs::read_to_string(fixture_path("missing2ref.out")).unwrap();

    let (out, err, code) = run(&["+missing2ref", "--no-version", input.to_str().unwrap()]);
    assert_eq!(code, 0, "+missing2ref failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn missing2ref_via_plugin_subcommand() {
    let input = fixture_path("plugin1.vcf");
    let expected = std::fs::read_to_string(fixture_path("missing2ref.out")).unwrap();

    let (out, err, code) = run(&[
        "plugin",
        "missing2ref",
        "--no-version",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "plugin missing2ref failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn missing2ref_reads_bcf_and_round_trips() {
    let dir = TempDir::new().unwrap();
    let input = fixture_path("plugin1.vcf");
    let bcf = dir.path().join("in.bcf");

    let (_o, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["+missing2ref", "--no-version", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "+missing2ref on bcf failed: {err}");
    assert!(out.contains("\tGT:GQ\t0/0:245\t0/0:245"), "{out}");
    assert!(
        !out.contains("./."),
        "missing genotypes should be gone: {out}"
    );
}

#[test]
fn missing2ref_writes_bgzf_output() {
    let dir = TempDir::new().unwrap();
    let input = fixture_path("plugin1.vcf");
    let out_path = dir.path().join("out.vcf.gz");

    let (_o, err, code) = run(&[
        "+missing2ref",
        "--no-version",
        "-Oz",
        "-o",
        out_path.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "+missing2ref -Oz failed: {err}");
    let bytes = std::fs::read(&out_path).unwrap();
    assert!(bytes.starts_with(&[0x1f, 0x8b]), "expected gzip magic");
}

#[test]
fn missing2ref_reads_stdin() {
    let input = fixture_path("plugin1.vcf");
    let expected = std::fs::read_to_string(fixture_path("missing2ref.out")).unwrap();
    let data = std::fs::read(&input).unwrap();

    ensure_binary_built();
    let mut child = Command::new(bin_path())
        .args(["+missing2ref", "--no-version", "-"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn +missing2ref");
    child.stdin.as_mut().unwrap().write_all(&data).unwrap();
    let out = child.wait_with_output().unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "+missing2ref - failed: {stderr}");
    assert_eq!(String::from_utf8(out.stdout).unwrap(), expected);
}

#[test]
fn missing2ref_common_include_exclude_filters_records() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("filter.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
1\t10\t.\tC\tT\t.\tPASS\tDP=5\tGT\t./.\n\
1\t20\t.\tC\tG\t.\tPASS\tDP=8\tGT\t./.\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "+missing2ref",
        "--no-version",
        "-i",
        "DP=5",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "+missing2ref -i failed: {err}");
    assert!(out.contains("1\t10\t.\tC\tT\t.\tPASS\tDP=5\tGT\t0/0\n"));
    assert!(!out.contains("1\t20\t"));

    let (out, err, code) = run(&[
        "+missing2ref",
        "--no-version",
        "--exclude=DP=5",
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "+missing2ref -e failed: {err}");
    assert!(!out.contains("1\t10\t"));
    assert!(out.contains("1\t20\t.\tC\tG\t.\tPASS\tDP=8\tGT\t0/0\n"));
}

#[test]
fn missing2ref_phased_and_major_modes() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("major.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\tD\n\
1\t10\t.\tC\tT,G\t.\tPASS\t.\tGT\t1/1\t0/1\t./.\t2/.\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "+missing2ref",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "-m",
    ]);
    assert_eq!(code, 0, "+missing2ref -m failed: {err}");
    assert!(
        out.contains("\tGT\t1/1\t0/1\t1/1\t2/1\n"),
        "major allele should fill missing allele tokens: {out}"
    );

    let (out, err, code) = run(&[
        "+missing2ref",
        "--no-version",
        input.to_str().unwrap(),
        "--",
        "-p",
        "-m",
    ]);
    assert_eq!(code, 0, "+missing2ref -p -m failed: {err}");
    assert!(
        out.contains("\tGT\t1/1\t0/1\t1|1\t2|1\n"),
        "phased major allele should phase only replaced alleles: {out}"
    );
}
