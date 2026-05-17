//! End-to-end parity tests for `+guess-ploidy -v -rX` against the upstream
//! `guess-ploidy.{PL,GL}.out` fixtures (harness pipes `grep -v bcftools`).

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

fn check(input_vcf: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input_vcf);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let out = Command::new(bin_path())
        .args(["+guess-ploidy", input.to_str().unwrap(), "-v", "-rX"])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.contains("bcftools"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {input_vcf}");
}

#[test]
fn guess_ploidy_pl() {
    check("view.PL.vcf", "guess-ploidy.PL.out");
}

#[test]
fn guess_ploidy_gl() {
    check("view.GL.vcf", "guess-ploidy.GL.out");
}

#[test]
fn guess_ploidy_accepts_af_tag() {
    ensure_binary_built();
    let tmp = TempDir::new().expect("tempdir");
    let input = tmp.path().join("guess-ploidy-af-tag.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=CUSTOM_AF,Number=A,Type=Float,Description=\"Custom AF\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
X\t200\t.\tA\tC\t.\t.\tCUSTOM_AF=0.9\tPL\t0,10,100\n",
    )
    .unwrap();

    let out = Command::new(bin_path())
        .args([
            "+guess-ploidy",
            input.to_str().unwrap(),
            "-v",
            "-rX",
            "--AF-tag",
            "CUSTOM_AF",
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
    let sex_line = stdout
        .lines()
        .find(|line| line.starts_with("SEX\tS1"))
        .expect("SEX line");
    assert_eq!(sex_line.split('\t').nth(5), Some("1"));
}

#[test]
fn guess_ploidy_genome_shortcut_filters_interval() {
    ensure_binary_built();
    let tmp = TempDir::new().expect("tempdir");
    let input = tmp.path().join("guess-ploidy-genome.vcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
X\t2699520\t.\tA\tC\t.\t.\t.\tPL\t0,10,100\n\
X\t2699521\t.\tA\tC\t.\t.\t.\tPL\t0,10,100\n",
    )
    .unwrap();

    let out = Command::new(bin_path())
        .args(["+guess-ploidy", input.to_str().unwrap(), "-v", "-g", "b37"])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let sex_line = stdout
        .lines()
        .find(|line| line.starts_with("SEX\tS1"))
        .expect("SEX line");
    assert_eq!(sex_line.split('\t').nth(5), Some("1"));
}
