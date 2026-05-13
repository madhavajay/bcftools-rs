//! End-to-end tests for `bcftools_rs::commands::view`.

use std::io::Read as _;
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

#[test]
fn view_text_round_trip_emits_all_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    // The fixture has 22 records (lines 11..=32 of the file).
    assert_eq!(record_lines.len(), 21);
    // First record is `1\t105\t.\tTAAACCCTA\t...`
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn view_reads_plain_vcf_from_stdin_without_filename() {
    let input = std::fs::read(fixture_path("aa.vcf")).unwrap();
    let (out, err, code) = run_with_stdin(&["view", "--no-version"], &input);
    assert_eq!(code, 0, "view stdin failed: {err}");
    let record_lines: Vec<_> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(record_lines.len(), 21);
}

#[test]
fn view_header_only_drops_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(record_lines.is_empty());
    assert!(out.contains("#CHROM\t"));
}

#[test]
fn view_no_header_drops_header() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-H", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let header_lines: Vec<_> = out.lines().filter(|l| l.starts_with('#')).collect();
    assert!(header_lines.is_empty(), "header leaked: {header_lines:?}");
    let record_lines: Vec<_> = out.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(record_lines.len(), 21);
}

#[test]
fn view_region_filters_by_chrom_and_position_interval() {
    let path = fixture_path("aa.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        path.to_str().unwrap(),
        "20:80-95",
    ]);
    assert_eq!(code, 0, "view region failed: {err}");
    let records: Vec<_> = out.lines().filter(|line| !line.is_empty()).collect();
    assert_eq!(
        records,
        [
            "20\t81\t.\tA\tC\t999\tPASS\t.",
            "20\t84\t.\tG\tT\t999\tPASS\t.",
            "20\t95\t.\tT\tA\t999\tPASS\t.",
            "20\t95\t.\tTCACCG\tAAAAAA\t999\tPASS\t.",
        ]
    );
}

#[test]
fn view_region_filters_bcf_input() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("aa.vcf");
    let bcf = tmp.path().join("aa.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["view", "--no-version", "-H", bcf.to_str().unwrap(), "20"]);
    assert_eq!(code, 0, "view BCF region failed: {err}");
    let records: Vec<_> = out.lines().filter(|line| !line.is_empty()).collect();
    assert_eq!(records.len(), 12);
    assert!(records.iter().all(|line| line.starts_with("20\t")));
}

#[test]
fn view_samples_list_subsets_vcf_columns() {
    let path = fixture_path("query.smpl.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-s", "11", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -s failed: {err}");
    assert_eq!(
        out,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t11\n\
chr1\t10000\t.\tA\tC\t.\t.\t.\tGT\t1/1\n"
    );
}

#[test]
fn view_samples_file_exclusion_subsets_bcf_input() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("query.smpl.vcf");
    let bcf = tmp.path().join("query.smpl.bcf");
    let samples = fixture_path("query.smpl.11.txt");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let excluded = format!("^{}", samples.display());
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-S",
        &excluded,
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -S ^file BCF failed: {err}");
    let chrom = out
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("#CHROM line");
    assert_eq!(
        chrom,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t00"
    );
    assert!(out.contains("chr1\t10000\t.\tA\tC\t.\t.\t.\tGT\t0/0\n"));
}

#[test]
fn view_drop_genotypes_matches_upstream_fixture() {
    let path = fixture_path("view.omitgenotypes.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.dropgenotypes.out")).unwrap();
    let (out, err, code) = run(&["view", "--no-version", "-G", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -G failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_drop_genotypes_no_header_matches_upstream_fixture() {
    let path = fixture_path("view.omitgenotypes.vcf");
    let expected =
        std::fs::read_to_string(fixture_path("view.dropgenotypes.noheader.out")).unwrap();
    let (out, err, code) = run(&["view", "--no-version", "-HG", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -HG failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_threads_writes_bgzf_vcf_output() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("aa.vcf");
    let output = tmp.path().join("aa.vcf.gz");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "--threads",
        "2",
        "-Oz",
        "-o",
        output.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view --threads -Oz failed: {err}");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    let records = decoded
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .count();
    assert_eq!(records, 21);
}

#[test]
fn view_threads_rejects_non_integer_argument() {
    let path = fixture_path("aa.vcf");
    let (_out, err, code) = run(&["view", "--threads", "abc", path.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("Could not parse argument: --threads abc"));
}

#[test]
fn view_unknown_output_type_errors() {
    let path = fixture_path("aa.vcf");
    let (_out, err, code) = run(&["view", "-Oq", path.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("not recognised"), "stderr: {err}");
}

#[test]
fn view_default_injects_bcftools_version_and_command_lines() {
    // Without --no-version, `bcftools view` must emit
    // ##bcftools_viewVersion=<v> and ##bcftools_viewCommand=<cmdline>; Date=...
    // header lines, mirroring upstream `bcf_hdr_append_version`.
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);

    let header_lines: Vec<_> = out.lines().filter(|l| l.starts_with("##")).collect();
    let version_line = header_lines
        .iter()
        .find(|l| l.starts_with("##bcftools_viewVersion="))
        .unwrap_or_else(|| panic!("missing version line in header:\n{out}"));
    let command_line = header_lines
        .iter()
        .find(|l| l.starts_with("##bcftools_viewCommand="))
        .unwrap_or_else(|| panic!("missing command line in header:\n{out}"));

    // The version line ends with the htslib-rs version we're built against.
    assert!(version_line.contains("+htslib-"), "got: {version_line}");
    // The command line names the subcommand and includes a `Date=` field.
    assert!(command_line.contains("view"), "got: {command_line}");
    assert!(command_line.contains("; Date="), "got: {command_line}");
}

#[test]
fn view_no_version_suppresses_injected_header_lines() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(
        !out.contains("##bcftools_viewVersion="),
        "version line leaked despite --no-version:\n{out}"
    );
    assert!(
        !out.contains("##bcftools_viewCommand="),
        "command line leaked despite --no-version:\n{out}"
    );
}
