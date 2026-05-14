//! End-to-end tests for `bcftools_rs::commands::sort`.

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

const UNSORTED_VCF: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##contig=<ID=2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
2\t25\t.\tA\tT\t100\tPASS\t.\n\
1\t20\t.\tC\tT\t100\tPASS\t.\n\
1\t10\t.\tA\tG\t100\tPASS\t.\n\
1\t10\t.\tA\tC\t100\tPASS\t.\n\
2\t15\t.\tT\tA\t100\tPASS\t.\n";

const KESTREL_HEADER_VCF: &str = "##fileformat=VCF4.2\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\tC\t.\tPASS\t.\n";

#[test]
fn sort_writes_records_in_contig_position_ref_alt_order() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.vcf");
    std::fs::write(&input, UNSORTED_VCF).unwrap();

    let (out, err, code) = run(&["sort", input.to_str().unwrap()]);
    assert_eq!(code, 0, "sort failed: {err}");

    let records: Vec<_> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t10\t.\tA\tC\t100\tPASS\t.",
            "1\t10\t.\tA\tG\t100\tPASS\t.",
            "1\t20\t.\tC\tT\t100\tPASS\t.",
            "2\t15\t.\tT\tA\t100\tPASS\t.",
            "2\t25\t.\tA\tT\t100\tPASS\t.",
        ]
    );
}

#[test]
fn sort_supports_vntyper_compressed_write_index_command_shape() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("output_indel.vcf");
    let output = dir.path().join("output_indel.vcf.gz");
    std::fs::write(&input, UNSORTED_VCF).unwrap();

    let (_out, err, code) = run(&[
        "sort",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-W",
        "-O",
        "z",
    ]);
    assert_eq!(code, 0, "sort -W -O z failed: {err}");
    assert!(output.exists(), "compressed VCF output not produced");

    let csi = dir.path().join("output_indel.vcf.gz.csi");
    assert!(csi.exists(), "CSI index not produced for -W");
    assert!(!std::fs::read(&csi).unwrap().is_empty(), "CSI is empty");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(&output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    let records: Vec<_> = decoded
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(records[0], "1\t10\t.\tA\tC\t100\tPASS\t.");
    assert_eq!(
        records.last().copied(),
        Some("2\t25\t.\tA\tT\t100\tPASS\t.")
    );
}

#[test]
fn sort_threads_writes_bgzf_vcf_output() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.vcf");
    let output = dir.path().join("threaded.vcf.gz");
    std::fs::write(&input, UNSORTED_VCF).unwrap();

    let (_out, err, code) = run(&[
        "sort",
        "--threads",
        "2",
        "-O",
        "z",
        "-o",
        output.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "sort --threads -O z failed: {err}");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(&output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    let records: Vec<_> = decoded
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(records[0], "1\t10\t.\tA\tC\t100\tPASS\t.");
    assert_eq!(
        records.last().copied(),
        Some("2\t25\t.\tA\tT\t100\tPASS\t.")
    );
}

#[test]
fn sort_threads_rejects_non_integer_argument() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.vcf");
    std::fs::write(&input, UNSORTED_VCF).unwrap();

    let (_out, err, code) = run(&["sort", "--threads", "abc", input.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("Could not parse argument: --threads abc"));
}

#[test]
fn sort_accepts_kestrel_non_canonical_fileformat_header() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("kestrel-style.vcf");
    let output = dir.path().join("kestrel-style.sorted.vcf");
    std::fs::write(&input, KESTREL_HEADER_VCF).unwrap();

    let (_out, err, code) = run(&[
        "sort",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-O",
        "v",
    ]);
    assert_eq!(code, 0, "sort rejected Kestrel header: {err}");
    assert!(
        err.contains("[W::bcf_get_version] Couldn't get VCF version, considering as 4.2"),
        "expected upstream-style warning, got: {err}"
    );

    let sorted = std::fs::read_to_string(&output).unwrap();
    assert!(
        sorted.contains("##fileformat=VCFv4.2"),
        "fileformat header not normalized in output: {sorted}"
    );
    let records: Vec<_> = sorted
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(records, ["chr1\t1\t.\tA\tC\t.\tPASS\t."]);
}

#[test]
fn sort_accepts_kestrel_header_with_compressed_write_index() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("kestrel-style.vcf");
    let output = dir.path().join("kestrel-style.sorted.vcf.gz");
    std::fs::write(&input, KESTREL_HEADER_VCF).unwrap();

    let (_out, err, code) = run(&[
        "sort",
        input.to_str().unwrap(),
        "-o",
        output.to_str().unwrap(),
        "-W",
        "-O",
        "z",
    ]);
    assert_eq!(code, 0, "sort -W -O z rejected Kestrel header: {err}");
    assert!(output.exists(), "compressed VCF output not produced");
    let csi = dir.path().join("kestrel-style.sorted.vcf.gz.csi");
    assert!(csi.exists(), "CSI index not produced for -W");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(&output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    assert!(
        decoded.contains("##fileformat=VCFv4.2"),
        "fileformat header not normalized in output: {decoded}"
    );
}

#[test]
fn sort_does_not_warn_for_canonical_fileformat_header() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("canonical.vcf");
    std::fs::write(&input, UNSORTED_VCF).unwrap();

    let (_out, err, code) = run(&["sort", input.to_str().unwrap()]);
    assert_eq!(code, 0, "sort failed: {err}");
    assert!(
        !err.contains("Couldn't get VCF version"),
        "unexpected warning for canonical header: {err}"
    );
}

#[test]
fn sort_spills_to_temp_dir_when_max_mem_is_small() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("unsorted.vcf");
    let temp = dir.path().join("sort-tmp");
    std::fs::write(&input, UNSORTED_VCF).unwrap();

    let (out, err, code) = run(&[
        "sort",
        "-m",
        "1",
        "-T",
        temp.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "sort external path failed: {err}");

    let records: Vec<_> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t10\t.\tA\tC\t100\tPASS\t.",
            "1\t10\t.\tA\tG\t100\tPASS\t.",
            "1\t20\t.\tC\tT\t100\tPASS\t.",
            "2\t15\t.\tT\tA\t100\tPASS\t.",
            "2\t25\t.\tA\tT\t100\tPASS\t.",
        ]
    );
    assert!(temp.exists(), "temp directory was not created");
    let leftovers = std::fs::read_dir(&temp).unwrap().count();
    assert_eq!(leftovers, 0, "temporary chunk files were not cleaned up");
}
