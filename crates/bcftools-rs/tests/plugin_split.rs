//! End-to-end parity tests for the filter-free `+split` upstream fixtures.
//!
//! Mirrors `test_plugin_split`: sort output files, then for each file append
//! the file name, `bcftools query -l`, and `bcftools view -H`.

use std::path::{Path, PathBuf};
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

fn run_bcftools(args: &[String]) -> std::process::Output {
    Command::new(bin_path())
        .args(args)
        .output()
        .expect("spawn bcftools")
}

fn query_list(path: &Path) -> String {
    let out = run_bcftools(&["query".into(), "-l".into(), path.display().to_string()]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "query -l failed for {}: {}",
        path.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn view_records(path: &Path) -> String {
    let out = run_bcftools(&["view".into(), "-H".into(), path.display().to_string()]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "view -H failed for {}: {}",
        path.display(),
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn check(input_fixture: &str, args: &[&str], tmp_name: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input_fixture);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let dir = std::env::temp_dir().join(format!(
        "bcftools-rs-split-test-{}-{tmp_name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);

    let mut full = vec![
        "+split".to_string(),
        input.to_str().unwrap().to_string(),
        "-o".to_string(),
        dir.to_str().unwrap().to_string(),
    ];
    full.extend(args.iter().map(|s| s.to_string()));
    let out = run_bcftools(&full);
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| !n.starts_with('.'))
        .collect();
    files.sort();

    let mut actual = String::new();
    for file in &files {
        let path = dir.join(file);
        actual.push_str(file);
        actual.push('\n');
        actual.push_str(&query_list(&path));
        actual.push_str(&view_records(&path));
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(actual, expected, "mismatch for {full:?}");
}

#[test]
fn split_default_single_sample_outputs() {
    check("split.1.vcf", &[], "default", "split.1.1.out");
}

#[test]
fn split_samples_file_single_sample_renames() {
    check(
        "split.1.vcf",
        &["-S", fixture_path("split.smpl.1.2.txt").to_str().unwrap()],
        "samples12",
        "split.1.2.out",
    );
}

#[test]
fn split_samples_file_multi_sample_renames() {
    check(
        "split.1.vcf",
        &["-S", fixture_path("split.smpl.1.3.txt").to_str().unwrap()],
        "samples13",
        "split.1.3.out",
    );
}

#[test]
fn split_groups_file() {
    check(
        "split.1.vcf",
        &["-G", fixture_path("split.grp.1.1.txt").to_str().unwrap()],
        "groups",
        "split.1.7.out",
    );
}

#[test]
fn split_sanitizes_duplicate_file_names() {
    check("split.2.vcf", &[], "sanitize", "split.2.1.out");
}

#[test]
fn split_keep_tags_projects_info_and_format() {
    ensure_binary_built();
    let root = std::env::temp_dir().join(format!(
        "bcftools-rs-split-keep-tags-{}",
        std::process::id()
    ));
    let out_dir = root.join("out");
    let input = root.join("input.vcf");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    std::fs::write(
        &input,
        "\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
##INFO=<ID=AA,Number=1,Type=String,Description=\"Allele\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
1\t1\t.\tA\tC\t.\t.\tDP=5;AA=T\tGT:AD:GQ\t0/1:2,3:9\t0/0:4,0:8
",
    )
    .unwrap();

    let out = run_bcftools(&[
        "+split".into(),
        input.display().to_string(),
        "-o".into(),
        out_dir.display().to_string(),
        "-k".into(),
        "INFO/DP,FMT/GT,AD".into(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+split -k failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let a = std::fs::read_to_string(out_dir.join("A.vcf")).unwrap();
    assert!(a.contains("##INFO=<ID=DP,"));
    assert!(!a.contains("##INFO=<ID=AA,"));
    assert!(a.contains("##FORMAT=<ID=GT,"));
    assert!(a.contains("##FORMAT=<ID=AD,"));
    assert!(!a.contains("##FORMAT=<ID=GQ,"));
    assert!(a.contains("1\t1\t.\tA\tC\t.\t.\tDP=5\tGT:AD\t0/1:2,3\n"));

    let b = std::fs::read_to_string(out_dir.join("B.vcf")).unwrap();
    assert!(b.contains("1\t1\t.\tA\tC\t.\t.\tDP=5\tGT:AD\t0/0:4,0\n"));
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn split_writes_bgzf_vcf_outputs() {
    ensure_binary_built();
    let input = fixture_path("split.1.vcf");
    let dir = std::env::temp_dir().join(format!("bcftools-rs-split-bgzf-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    let out = run_bcftools(&[
        "+split".into(),
        input.display().to_string(),
        "-o".into(),
        dir.display().to_string(),
        "-Oz".into(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+split -Oz failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let a1 = dir.join("A1.vcf.gz");
    assert!(a1.exists(), "missing {}", a1.display());
    assert_eq!(query_list(&a1), "A1\n");
    assert!(view_records(&a1).contains("22\t10\t.\tC\tA\t.\t.\t.\tGT\t./.\n"));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn split_include_filter_applies_to_each_projected_output() {
    ensure_binary_built();
    let root =
        std::env::temp_dir().join(format!("bcftools-rs-split-filter-{}", std::process::id()));
    let out_dir = root.join("out");
    let input = root.join("input.vcf");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).unwrap();
    std::fs::write(
        &input,
        "\
##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB
1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/1\t0/0
1\t2\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1
",
    )
    .unwrap();

    let out = run_bcftools(&[
        "+split".into(),
        input.display().to_string(),
        "-o".into(),
        out_dir.display().to_string(),
        "-i".into(),
        r#"GT="alt""#.into(),
    ]);
    assert_eq!(
        out.status.code(),
        Some(0),
        "+split -i failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let a_records = view_records(&out_dir.join("A.vcf"));
    assert!(a_records.contains("1\t1\t.\tA\tC"));
    assert!(!a_records.contains("1\t2\t.\tA\tG"));

    let b_records = view_records(&out_dir.join("B.vcf"));
    assert!(!b_records.contains("1\t1\t.\tA\tC"));
    assert!(b_records.contains("1\t2\t.\tA\tG"));
    let _ = std::fs::remove_dir_all(&root);
}
