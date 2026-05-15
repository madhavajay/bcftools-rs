//! End-to-end tests for `bcftools_rs::commands::reheader`.

use std::io::Read as _;
use std::path::PathBuf;
use std::process::{Command, Stdio};

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

fn run_with_stdin(args: &[&str], input: &[u8]) -> (String, String, i32) {
    ensure_binary_built();
    let mut child = Command::new(bin_path())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn bcftools");
    {
        use std::io::Write as _;
        child
            .stdin
            .as_mut()
            .expect("stdin")
            .write_all(input)
            .expect("write stdin");
    }
    let out = child.wait_with_output().expect("wait bcftools");
    let stdout = String::from_utf8_lossy(&out.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

fn write_vcf(dir: &TempDir) -> PathBuf {
    let path = dir.path().join("input.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n",
    )
    .unwrap();
    path
}

#[test]
fn reheader_samples_file_replaces_all_sample_names() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "X\nY\n").unwrap();

    let (out, err, code) = run(&[
        "reheader",
        "-s",
        samples.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader -s failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
    assert!(out.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
}

#[test]
fn reheader_samples_list_renames_matching_pairs() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);

    let (out, err, code) = run(&["reheader", "-n", "B Z", input.to_str().unwrap()]);
    assert_eq!(code, 0, "reheader -n failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tZ\n"));
}

#[test]
fn reheader_samples_list_long_option_replaces_all_sample_names() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);

    let (out, err, code) = run(&["reheader", "--samples-list", "X,Y", input.to_str().unwrap()]);
    assert_eq!(code, 0, "reheader --samples-list failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
    assert!(out.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
}

#[test]
fn reheader_output_writes_file() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let output = dir.path().join("renamed.vcf");

    let (stdout, err, code) = run(&[
        "reheader",
        "-n",
        "X,Y",
        "-o",
        output.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader -o failed: {err}");
    assert!(stdout.is_empty());
    let out = std::fs::read_to_string(output).unwrap();
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
}

#[test]
fn reheader_reads_plain_vcf_from_stdin() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "X\nY\n").unwrap();
    let data = std::fs::read(input).unwrap();

    let (out, err, code) = run_with_stdin(&["reheader", "-s", samples.to_str().unwrap()], &data);
    assert_eq!(code, 0, "reheader stdin failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
}

#[test]
fn reheader_bgzf_vcf_writes_bgzf_output() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let compressed = dir.path().join("input.vcf.gz");
    let output = dir.path().join("renamed.vcf.gz");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Oz",
        "-o",
        compressed.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Oz failed: {err}");

    let (_out, err, code) = run(&[
        "reheader",
        "-n",
        "X,Y",
        "-o",
        output.to_str().unwrap(),
        compressed.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader .vcf.gz failed: {err}");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    assert!(decoded.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
    assert!(decoded.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
}

#[test]
fn reheader_threads_writes_bgzf_vcf_output() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let compressed = dir.path().join("input.vcf.gz");
    let output = dir.path().join("renamed.vcf.gz");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Oz",
        "-o",
        compressed.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Oz failed: {err}");

    let (_out, err, code) = run(&[
        "reheader",
        "--threads",
        "2",
        "-n",
        "X,Y",
        "-o",
        output.to_str().unwrap(),
        compressed.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader --threads .vcf.gz failed: {err}");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    assert!(decoded.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tX\tY\n"));
    assert!(decoded.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
}

#[test]
fn reheader_bcf_samples_file_rewrites_bcf_header() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let bcf = dir.path().join("input.bcf");
    let output = dir.path().join("renamed.bcf");
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "X\nY\n").unwrap();

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (_out, err, code) = run(&[
        "reheader",
        "-s",
        samples.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader BCF failed: {err}");

    let (out, err, code) = run(&["query", "-l", output.to_str().unwrap()]);
    assert_eq!(code, 0, "query -l renamed BCF failed: {err}");
    assert_eq!(out, "X\nY\n");
}

#[test]
fn reheader_threads_writes_bcf_output() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let bcf = dir.path().join("input.bcf");
    let output = dir.path().join("renamed.bcf");
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "X\nY\n").unwrap();

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (_out, err, code) = run(&[
        "reheader",
        "--threads",
        "2",
        "-s",
        samples.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader --threads BCF failed: {err}");

    let (out, err, code) = run(&["query", "-l", output.to_str().unwrap()]);
    assert_eq!(code, 0, "query -l threaded renamed BCF failed: {err}");
    assert_eq!(out, "X\nY\n");
}

#[test]
fn reheader_bcf_in_place_rewrites_same_file() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let bcf = dir.path().join("input.bcf");
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "X\nY\n").unwrap();

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (_out, err, code) = run(&[
        "reheader",
        "--in-place",
        "-s",
        samples.to_str().unwrap(),
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader --in-place failed: {err}");

    let (out, err, code) = run(&["query", "-l", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "query -l in-place BCF failed: {err}");
    assert_eq!(out, "X\nY\n");
}

#[test]
fn reheader_threads_rejects_non_integer_argument() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);

    let (_out, err, code) = run(&[
        "reheader",
        "--threads",
        "abc",
        "-n",
        "X,Y",
        input.to_str().unwrap(),
    ]);
    assert_ne!(code, 0);
    assert!(err.contains("Could not parse argument: --threads abc"));
}

#[test]
fn reheader_header_replacement_uses_new_header() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let header = dir.path().join("header.txt");
    std::fs::write(
        &header,
        "##fileformat=VCFv4.3\n\
##source=reheader-test\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "reheader",
        "-h",
        header.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader -h failed: {err}");
    assert!(out.starts_with(
        "##fileformat=VCFv4.3\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##source=reheader-test\n#CHROM\t"
    ));
    assert!(out.ends_with("1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\n"));
}

#[test]
fn reheader_fai_updates_contig_lines() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let fai = dir.path().join("ref.fa.fai");
    std::fs::write(&fai, "1\t1000\t0\t80\t81\n2\t2000\t1009\t80\t81\n").unwrap();

    let (out, err, code) = run(&[
        "reheader",
        "-f",
        fai.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "reheader -f failed: {err}");
    assert!(out.contains("##contig=<ID=1,length=1000>\n"));
    assert!(out.contains("##contig=<ID=2,length=2000>\n#CHROM\t"));
}

#[test]
fn reheader_in_place_rejects_text_vcf() {
    let dir = TempDir::new().expect("tempdir");
    let input = write_vcf(&dir);
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "X\nY\n").unwrap();

    let (_out, err, code) = run(&[
        "reheader",
        "--in-place",
        "-s",
        samples.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_ne!(code, 0);
    assert!(err.contains("--in-place is only supported for BCF input"));
}
