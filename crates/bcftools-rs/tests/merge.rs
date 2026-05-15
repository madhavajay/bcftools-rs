//! End-to-end tests for `bcftools_rs::commands::merge`.

use std::path::PathBuf;
use std::process::Command;

use tempfile::TempDir;

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

fn write_vcf(dir: &TempDir, name: &str, sample: &str, gt: &str) -> PathBuf {
    let path = dir.path().join(name);
    let body = format!(
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n\
1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t{gt}\n"
    );
    std::fs::write(&path, body).unwrap();
    path
}

#[test]
fn merge_combines_same_site_samples() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b = write_vcf(&dir, "b.vcf", "SAMPLE_B", "1/1");

    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge failed: {err}");
    assert!(
        out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B"),
        "merged header missing: {out}"
    );
    assert!(
        out.contains("1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1"),
        "merged record missing: {out}"
    );
}

#[test]
fn merge_rejects_duplicate_sample_without_force() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "DUP", "0/1");
    let b = write_vcf(&dir, "b.vcf", "DUP", "1/1");

    let (_out, err, code) = run(&[
        "merge",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "expected duplicate-sample failure, got success");
    assert!(
        err.contains("duplicate sample name") && err.contains("DUP"),
        "stderr should mention duplicate sample 'DUP': {err}"
    );
}

#[test]
fn merge_force_samples_prefixes_duplicates() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "DUP", "0/1");
    let b = write_vcf(&dir, "b.vcf", "DUP", "1/1");

    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "--force-samples",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge --force-samples failed: {err}");
    assert!(
        out.contains("FORMAT\tDUP\t2:DUP"),
        "expected DUP and 2:DUP columns: {out}"
    );
}

#[test]
fn merge_rejects_record_set_mismatch() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b_path = dir.path().join("b.vcf");
    std::fs::write(
        &b_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_B\n\
1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "merge",
        "--no-version",
        a.to_str().unwrap(),
        b_path.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "expected mismatch failure, got success");
    assert!(
        err.contains("record mismatch") || err.contains("compatible"),
        "stderr should mention record mismatch: {err}"
    );
}

#[test]
fn merge_writes_bgzf_vcf_output() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b = write_vcf(&dir, "b.vcf", "SAMPLE_B", "1/1");
    let out_path = dir.path().join("merged.vcf.gz");

    let (_stdout, err, code) = run(&[
        "merge",
        "--no-version",
        "-Oz",
        "-o",
        out_path.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge -Oz failed: {err}");

    let bytes = std::fs::read(&out_path).unwrap();
    assert!(
        bytes.starts_with(&[0x1f, 0x8b]),
        "output should start with gzip magic: {:?}",
        &bytes[..bytes.len().min(4)]
    );
}

#[test]
fn merge_reads_file_list() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b = write_vcf(&dir, "b.vcf", "SAMPLE_B", "1/1");
    let list = dir.path().join("inputs.txt");
    std::fs::write(&list, format!("{}\n{}\n", a.display(), b.display())).unwrap();

    let (out, err, code) = run(&["merge", "--no-version", "-l", list.to_str().unwrap()]);
    assert_eq!(code, 0, "merge -l failed: {err}");
    assert!(out.contains("SAMPLE_A\tSAMPLE_B"), "{out}");
}

#[test]
fn merge_rejects_single_input() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");

    let (_out, err, code) = run(&["merge", "--no-version", a.to_str().unwrap()]);
    assert_ne!(code, 0, "expected single-input rejection");
    assert!(
        err.contains("at least two") || err.contains("expected at least"),
        "stderr should request at least two inputs: {err}"
    );
}
