//! Synthetic integration coverage for `+frameshifts`.
//!
//! Upstream has no direct fixture row for this plugin in `test.pl`, so the
//! tests exercise the core OOF annotation behavior and compressed output.

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
    let input = root.join("input.vcf");
    let exons = root.join("exons.bed");
    std::fs::write(
        &input,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t10\t.\tA\tAT,ATGC,C\t.\t.\t.
1\t20\t.\tATGC\tA\t.\t.\tDP=1
1\t40\t.\tA\tAT\t.\t.\t.
",
    )
    .unwrap();
    std::fs::write(&exons, "1\t9\t30\n").unwrap();
    (input, exons)
}

#[test]
fn frameshifts_annotates_oof_info() {
    ensure_binary_built();
    let root = std::env::temp_dir().join(format!("bcftools-rs-frameshifts-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let (input, exons) = write_inputs(&root);

    let out = run_bcftools(&[
        "+frameshifts".into(),
        input.display().to_string(),
        "--".into(),
        "-e".into(),
        exons.display().to_string(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+frameshifts failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    assert!(stdout.contains("##INFO=<ID=OOF,"));
    assert!(stdout.contains("1\t10\t.\tA\tAT,ATGC,C\t.\t.\tOOF=1,0,-1\n"));
    assert!(stdout.contains("1\t20\t.\tATGC\tA\t.\t.\tDP=1;OOF=0\n"));
    assert!(stdout.contains("1\t40\t.\tA\tAT\t.\t.\t.\n"));
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn frameshifts_writes_bgzf_vcf_output() {
    ensure_binary_built();
    let root = std::env::temp_dir().join(format!(
        "bcftools-rs-frameshifts-bgzf-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let (input, exons) = write_inputs(&root);
    let output = root.join("out.vcf.gz");

    let out = run_bcftools(&[
        "+frameshifts".into(),
        input.display().to_string(),
        "-Oz".into(),
        "-o".into(),
        output.display().to_string(),
        "--".into(),
        "-e".into(),
        exons.display().to_string(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+frameshifts -Oz failed: {}",
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
    assert!(records.contains("1\t10\t.\tA\tAT,ATGC,C\t.\t.\tOOF=1,0,-1\n"));
    let _ = std::fs::remove_dir_all(&root);
}
