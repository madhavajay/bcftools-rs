//! Synthetic integration coverage for `+isecGT`.
//!
//! Upstream has no dedicated `test.pl` fixture for this plugin. These tests
//! cover the core behavior from `plugins/isecGT.c`: compare matching records
//! across two files, map samples by name, and set non-identical genotypes in
//! the first file to missing.

use std::path::{Path, PathBuf};
use std::process::Command;

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

fn run_bcftools(args: &[String]) -> std::process::Output {
    Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools")
}

fn write_inputs(root: &Path) -> (PathBuf, PathBuf) {
    let a = root.join("a.vcf");
    let b = root.join("b.vcf");
    std::fs::write(
        &a,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:8\t0/0:9
1\t2\t.\tG\tT\t.\t.\t.\tGT:DP\t1|0:7\t0/0:6
1\t3\t.\tC\tG\t.\t.\t.\tGT:DP\t0/1:5\t0/1:4
1\t4\t.\tT\tC\t.\t.\t.\tGT:DP\t0/0:3\t0/1:2
",
    )
    .unwrap();
    std::fs::write(
        &b,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tB\tA
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/1
1\t2\t.\tG\tT\t.\t.\t.\tGT\t0/1\t0|1
1\t4\t.\tT\tC\t.\t.\t.\tGT\t0/0\t0/0
",
    )
    .unwrap();
    (a, b)
}

#[test]
fn isecgt_sets_non_identical_genotypes_to_missing() {
    ensure_binary_built();
    let root = std::env::temp_dir().join(format!("bcftools-rs-isecgt-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let (a, b) = write_inputs(&root);

    let out = run_bcftools(&[
        "+isecGT".into(),
        a.display().to_string(),
        b.display().to_string(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+isecGT failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let actual = String::from_utf8(out.stdout).unwrap();
    let expected = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
1\t1\t.\tA\tC\t.\t.\t.\tGT:DP\t0/1:8\t0/0:9
1\t2\t.\tG\tT\t.\t.\t.\tGT:DP\t./.:7\t./.:6
1\t3\t.\tC\tG\t.\t.\t.\tGT:DP\t0/1:5\t0/1:4
1\t4\t.\tT\tC\t.\t.\t.\tGT:DP\t0/0:3\t./.:2
";
    assert_eq!(actual, expected);
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn isecgt_writes_bgzf_vcf_output() {
    ensure_binary_built();
    let root = std::env::temp_dir().join(format!("bcftools-rs-isecgt-bgzf-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let (a, b) = write_inputs(&root);
    let output = root.join("out.vcf.gz");

    let out = run_bcftools(&[
        "+isecGT".into(),
        a.display().to_string(),
        b.display().to_string(),
        "-Oz".into(),
        "-o".into(),
        output.display().to_string(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+isecGT -Oz failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let viewed = run_bcftools(&["view".into(), "-H".into(), output.display().to_string()]);
    assert_eq!(
        viewed.status.code(),
        Some(0),
        "view -H failed: {}",
        String::from_utf8_lossy(&viewed.stderr)
    );
    let records = String::from_utf8(viewed.stdout).unwrap();
    assert!(records.contains("1\t2\t.\tG\tT\t.\t.\t.\tGT:DP\t./.:7\t./.:6\n"));
    assert!(records.contains("1\t4\t.\tT\tC\t.\t.\t.\tGT:DP\t0/0:3\t./.:2\n"));
    let _ = std::fs::remove_dir_all(&root);
}
