//! Synthetic integration coverage for `+check-sparsity`.
//!
//! Upstream does not exercise this plugin directly in `test.pl`; these tests
//! cover per-contig reports, `-n`, and region filtering.

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

fn write_input(root: &Path) -> PathBuf {
    let input = root.join("input.vcf");
    std::fs::write(
        &input,
        "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t./.\t.
1\t2\t.\tA\tG\t.\t.\t.\tGT:DP\t0/1:5\t0/.:7\t.:9
2\t1\t.\tG\tT\t.\t.\t.\tGT\t./.\t0/0\t0/1
",
    )
    .unwrap();
    input
}

fn run_bcftools(args: &[String]) -> std::process::Output {
    Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools")
}

#[test]
fn check_sparsity_reports_samples_by_chromosome() {
    ensure_binary_built();
    let root =
        std::env::temp_dir().join(format!("bcftools-rs-check-sparsity-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let input = write_input(&root);

    let out = run_bcftools(&["+check-sparsity".into(), input.display().to_string()]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+check-sparsity failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(String::from_utf8(out.stdout).unwrap(), "1\tC\n2\tA\n");
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn check_sparsity_supports_threshold_and_regions_file() {
    ensure_binary_built();
    let root = std::env::temp_dir().join(format!(
        "bcftools-rs-check-sparsity-regions-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    let input = write_input(&root);
    let regions = root.join("regions.txt");
    std::fs::write(&regions, "1:1-2\n2:1\n").unwrap();

    let out = run_bcftools(&[
        "+check-sparsity".into(),
        input.display().to_string(),
        "--".into(),
        "-n".into(),
        "2".into(),
        "-R".into(),
        regions.display().to_string(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+check-sparsity -R failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        "1:1-2\tB\n1:1-2\tC\n2:1\tA\n2:1\tB\n2:1\tC\n"
    );
    let _ = std::fs::remove_dir_all(&root);
}
