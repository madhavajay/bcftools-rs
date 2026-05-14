//! End-to-end tests for `bcftools_rs::commands::isec`.

use std::io::Read as _;
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
    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

const VCF_A: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tC\t.\t.\tDP=5\n\
1\t20\tid2\tG\tGA\t.\t.\tDP=8\n\
1\t30\t.\tT\tG\t.\t.\tDP=2\n";

const VCF_B: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tC\t.\t.\tDP=7\n\
1\t20\tid2\tG\tGT\t.\t.\tDP=9\n\
1\t40\t.\tC\tA\t.\t.\tDP=4\n";

const VCF_C: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tC\t.\t.\tDP=6\n\
1\t30\t.\tT\tG\t.\t.\tDP=3\n";

fn write_temp(dir: &TempDir, name: &str, body: &str) -> PathBuf {
    let p = dir.path().join(name);
    std::fs::write(&p, body).unwrap();
    p
}

#[test]
fn isec_default_prints_exact_intersection_bitmap() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&["isec", "-n", "=2", a.to_str().unwrap(), b.to_str().unwrap()]);
    assert_eq!(code, 0, "isec failed: {err}");
    assert_eq!(out, "1\t10\tA\tC\t11\n");
}

#[test]
fn isec_collapse_any_matches_same_position_records() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-c",
        "any",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -c any failed: {err}");
    let records: Vec<_> = out.lines().collect();
    assert_eq!(records, ["1\t10\tA\tC\t11", "1\t20\tG\tGA\t11"]);
}

#[test]
fn isec_complement_reports_first_input_private_records() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&[
        "isec",
        "-C",
        "-c",
        "any",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -C failed: {err}");
    assert_eq!(out, "1\t30\tT\tG\t10\n");
}

#[test]
fn isec_include_expression_filters_before_intersection() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-iDP>6",
        "-c",
        "any",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -i failed: {err}");
    assert_eq!(out, "1\t20\tG\tGA\t11\n");
}

#[test]
fn isec_write_input_vcf_preserves_requested_input_records() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-w",
        "2",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -w failed: {err}");
    assert!(out.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">"));
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    assert!(out.contains("1\t10\t.\tA\tC\t.\t.\tDP=7"));
    assert!(!out.contains("1\t40\t.\tC\tA"));
    assert!(!out.contains("##bcftools_isecVersion="));
}

#[test]
fn isec_single_input_targets_emit_filtered_vcf() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let targets = write_temp(&dir, "targets.tab", "1\t10\t20\n");
    let (out, err, code) = run(&[
        "isec",
        "--no-version",
        "-T",
        targets.to_str().unwrap(),
        a.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -T failed: {err}");
    let records: Vec<_> = out.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(
        records,
        [
            "1\t10\t.\tA\tC\t.\t.\tDP=5",
            "1\t20\tid2\tG\tGA\t.\t.\tDP=8"
        ]
    );
}

#[test]
fn isec_output_writes_site_bitmap_to_file() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("sites.txt");
    let (out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-o",
        out_path.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -o failed: {err}");
    assert!(out.is_empty());
    assert_eq!(
        std::fs::read_to_string(out_path).unwrap(),
        "1\t10\tA\tC\t11\n"
    );
}

#[test]
fn isec_prefix_writes_sites_readme_and_numbered_vcfs() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let prefix = dir.path().join("isec-out");
    let (_out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-c",
        "any",
        "--no-version",
        "-p",
        prefix.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -p failed: {err}");
    assert_eq!(
        std::fs::read_to_string(prefix.join("sites.txt")).unwrap(),
        "1\t10\tA\tC\t11\n1\t20\tG\tGA\t11\n"
    );
    let readme = std::fs::read_to_string(prefix.join("README.txt")).unwrap();
    assert!(readme.contains("This file was produced by vcfisec."));
    assert!(readme.contains("0000.vcf"));
    assert!(readme.contains("0001.vcf"));
    let first = std::fs::read_to_string(prefix.join("0000.vcf")).unwrap();
    let second = std::fs::read_to_string(prefix.join("0001.vcf")).unwrap();
    assert!(first.contains("1\t10\t.\tA\tC\t.\t.\tDP=5"));
    assert!(first.contains("1\t20\tid2\tG\tGA\t.\t.\tDP=8"));
    assert!(second.contains("1\t10\t.\tA\tC\t.\t.\tDP=7"));
    assert!(second.contains("1\t20\tid2\tG\tGT\t.\t.\tDP=9"));
    assert!(!first.contains("##bcftools_isecVersion="));
}

#[test]
fn isec_prefix_output_type_z_writes_numbered_gzip_vcfs() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let prefix = dir.path().join("isec-gz");
    let (_out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-w",
        "1",
        "--no-version",
        "-O",
        "z",
        "-p",
        prefix.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -p -O z failed: {err}");
    assert!(prefix.join("0000.vcf.gz").exists());
    assert!(prefix.join("0000.vcf.gz.tbi").exists());
    assert!(!prefix.join("0001.vcf.gz").exists());
    assert!(!prefix.join("0001.vcf.gz.tbi").exists());
    let mut decoded = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(prefix.join("0000.vcf.gz")).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    assert!(decoded.contains("1\t10\t.\tA\tC\t.\t.\tDP=5"));
}

#[test]
fn isec_reads_bcf_inputs_and_writes_indexed_prefix_bcf() {
    let dir = TempDir::new().unwrap();
    let a_vcf = write_temp(&dir, "a.vcf", VCF_A);
    let b_vcf = write_temp(&dir, "b.vcf", VCF_B);
    let a_bcf = dir.path().join("a.bcf");
    let b_bcf = dir.path().join("b.bcf");
    let (_a_out, a_err, a_code) = run(&[
        "view",
        "--no-version",
        "-O",
        "b",
        "-o",
        a_bcf.to_str().unwrap(),
        a_vcf.to_str().unwrap(),
    ]);
    assert_eq!(a_code, 0, "view -Ob a failed: {a_err}");
    let (_b_out, b_err, b_code) = run(&[
        "view",
        "--no-version",
        "-O",
        "b",
        "-o",
        b_bcf.to_str().unwrap(),
        b_vcf.to_str().unwrap(),
    ]);
    assert_eq!(b_code, 0, "view -Ob b failed: {b_err}");

    let prefix = dir.path().join("isec-bcf");
    let (_out, err, code) = run(&[
        "isec",
        "-n",
        "=2",
        "-w",
        "1",
        "--no-version",
        "-O",
        "b",
        "-p",
        prefix.to_str().unwrap(),
        a_bcf.to_str().unwrap(),
        b_bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec BCF failed: {err}");
    let out_bcf = prefix.join("0000.bcf");
    assert!(out_bcf.exists());
    assert!(prefix.join("0000.bcf.csi").exists());

    let (view_out, view_err, view_code) = run(&["view", "--no-version", out_bcf.to_str().unwrap()]);
    assert_eq!(view_code, 0, "view isec BCF failed: {view_err}");
    assert!(view_out.contains("1\t10\t.\tA\tC\t.\t.\tDP=5"));
    assert!(!view_out.contains("1\t20\tid2\tG\tGA"));
}

#[test]
fn isec_prefix_default_two_input_venn_writes_private_and_shared_files() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let prefix = dir.path().join("venn");
    let (_out, err, code) = run(&[
        "isec",
        "--no-version",
        "-p",
        prefix.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec default -p failed: {err}");
    assert_eq!(
        std::fs::read_to_string(prefix.join("sites.txt")).unwrap(),
        "1\t10\tA\tC\t11\n1\t20\tG\tGA\t10\n1\t30\tT\tG\t10\n1\t20\tG\tGT\t01\n1\t40\tC\tA\t01\n"
    );
    let private_a = std::fs::read_to_string(prefix.join("0000.vcf")).unwrap();
    let private_b = std::fs::read_to_string(prefix.join("0001.vcf")).unwrap();
    let shared_a = std::fs::read_to_string(prefix.join("0002.vcf")).unwrap();
    let shared_b = std::fs::read_to_string(prefix.join("0003.vcf")).unwrap();
    assert!(private_a.contains("1\t20\tid2\tG\tGA\t.\t.\tDP=8"));
    assert!(private_a.contains("1\t30\t.\tT\tG\t.\t.\tDP=2"));
    assert!(private_b.contains("1\t20\tid2\tG\tGT\t.\t.\tDP=9"));
    assert!(private_b.contains("1\t40\t.\tC\tA\t.\t.\tDP=4"));
    assert!(shared_a.contains("1\t10\t.\tA\tC\t.\t.\tDP=5"));
    assert!(shared_b.contains("1\t10\t.\tA\tC\t.\t.\tDP=7"));
}

#[test]
fn isec_nfiles_exact_bitmask_selects_requested_presence_pattern() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let c = write_temp(&dir, "c.vcf", VCF_C);
    let (out, err, code) = run(&[
        "isec",
        "-n~101",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        c.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "isec -n~ failed: {err}");
    assert_eq!(out, "1\t30\tT\tG\t101\n");
}
