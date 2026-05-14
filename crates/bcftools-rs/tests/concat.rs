//! End-to-end tests for `bcftools_rs::commands::concat`.

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
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##contig=<ID=2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tC\t100\tPASS\t.\n\
1\t200\t.\tG\tT\t100\tPASS\t.\n";

const VCF_B: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##contig=<ID=2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
2\t10\t.\tT\tA\t100\tPASS\t.\n\
2\t20\t.\tC\tG\t100\tPASS\t.\n";

const VCF_OVERLAP_B: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##contig=<ID=2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t150\t.\tT\tA\t100\tPASS\t.\n\
2\t20\t.\tC\tG\t100\tPASS\t.\n";

const VCF_DIFF_SAMPLES: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\n\
1\t100\t.\tA\tC\t100\tPASS\t.\tGT\t0/1\n";

const VCF_DIFF_SAMPLES_B: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA2\n\
1\t100\t.\tA\tC\t100\tPASS\t.\tGT\t0/1\n";

fn write_temp(dir: &TempDir, name: &str, body: &str) -> PathBuf {
    let p = dir.path().join(name);
    std::fs::write(&p, body).unwrap();
    p
}

#[test]
fn concat_two_text_vcfs_emits_records_in_input_order() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&["concat", a.to_str().unwrap(), b.to_str().unwrap()]);
    assert_eq!(code, 0, "concat failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t100\t.\tA\tC\t100\tPASS\t.",
            "1\t200\t.\tG\tT\t100\tPASS\t.",
            "2\t10\t.\tT\tA\t100\tPASS\t.",
            "2\t20\t.\tC\tG\t100\tPASS\t.",
        ]
    );
}

#[test]
fn concat_rejects_overlapping_inputs_unless_allow_overlaps() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_OVERLAP_B);
    let (_out, err, code) = run(&[
        "concat",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_ne!(code, 0);
    assert!(
        err.contains("Input files overlap at 1:150"),
        "expected overlap error, got: {err}"
    );

    let (out, err, code) = run(&[
        "concat",
        "--no-version",
        "--allow-overlaps",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "concat --allow-overlaps failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t100\t.\tA\tC\t100\tPASS\t.",
            "1\t200\t.\tG\tT\t100\tPASS\t.",
            "1\t150\t.\tT\tA\t100\tPASS\t.",
            "2\t20\t.\tC\tG\t100\tPASS\t.",
        ]
    );
}

#[test]
fn concat_default_injects_bcftools_version_and_command_lines() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let (out, err, code) = run(&["concat", a.to_str().unwrap()]);
    assert_eq!(code, 0, "concat failed: {err}");
    let version_line = out
        .lines()
        .find(|l| l.starts_with("##bcftools_concatVersion="))
        .unwrap_or_else(|| panic!("missing concat version line:\n{out}"));
    let command_line = out
        .lines()
        .find(|l| l.starts_with("##bcftools_concatCommand="))
        .unwrap_or_else(|| panic!("missing concat command line:\n{out}"));
    assert!(version_line.contains("+htslib-"), "got: {version_line}");
    assert!(command_line.contains("concat"), "got: {command_line}");
    assert!(command_line.contains("; Date="), "got: {command_line}");
}

#[test]
fn concat_no_version_suppresses_bcftools_header_lines() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let (out, err, code) = run(&["concat", "--no-version", a.to_str().unwrap()]);
    assert_eq!(code, 0, "concat --no-version failed: {err}");
    assert!(
        !out.contains("##bcftools_concatVersion="),
        "version line leaked despite --no-version:\n{out}"
    );
    assert!(
        !out.contains("##bcftools_concatCommand="),
        "command line leaked despite --no-version:\n{out}"
    );
}

#[test]
fn concat_writes_bgzf_with_output_type_z() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.vcf.gz");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
    ]);
    assert_eq!(code, 0, "concat -O z failed: {err}");

    let mut decoded = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&out_path).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    let records: Vec<&str> = decoded
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(records.len(), 4);
}

#[test]
fn concat_threads_writes_bgzf_vcf_output() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.vcf.gz");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "--threads",
        "2",
    ]);
    assert_eq!(code, 0, "concat --threads -O z failed: {err}");
    let mut decoded = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&out_path).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    assert_eq!(
        decoded
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty())
            .count(),
        4
    );
}

#[test]
fn concat_write_index_creates_csi_for_bgzf_vcf() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.vcf.gz");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "-W",
    ]);
    assert_eq!(code, 0, "concat -W -O z failed: {err}");
    assert!(
        dir.path().join("merged.vcf.gz.csi").exists(),
        "CSI index not created"
    );
}

#[test]
fn concat_write_index_tbi_creates_tbi_for_bgzf_vcf() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.vcf.gz");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "--write-index=tbi",
    ]);
    assert_eq!(code, 0, "concat --write-index=tbi failed: {err}");
    assert!(
        dir.path().join("merged.vcf.gz.tbi").exists(),
        "TBI index not created"
    );
}

#[test]
fn concat_writes_bcf_with_output_type_b() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.bcf");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "b",
    ]);
    assert_eq!(code, 0, "concat -O b failed: {err}");
    let bytes = std::fs::read(&out_path).unwrap();
    // BCF starts with the BGZF magic 1f 8b 08 04.
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04]);
}

#[test]
fn concat_threads_writes_bcf_output() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.bcf");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "b",
        "--threads=2",
    ]);
    assert_eq!(code, 0, "concat --threads -O b failed: {err}");
    let bytes = std::fs::read(&out_path).unwrap();
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04]);
}

#[test]
fn concat_write_index_creates_csi_for_bcf() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let out_path = dir.path().join("merged.bcf");
    let (_out, err, code) = run(&[
        "concat",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "b",
        "-W",
    ]);
    assert_eq!(code, 0, "concat -W -O b failed: {err}");
    assert!(
        dir.path().join("merged.bcf.csi").exists(),
        "BCF CSI index not created"
    );
}

#[test]
fn concat_file_list_reads_inputs_from_file() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let list = dir.path().join("files.list");
    std::fs::write(&list, format!("{}\n{}\n", a.display(), b.display())).unwrap();
    let (out, err, code) = run(&["concat", "-f", list.to_str().unwrap()]);
    assert_eq!(code, 0, "concat -f failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(records.len(), 4);
}

#[test]
fn concat_naive_preserves_first_header_and_appends_bodies() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&[
        "concat",
        "--naive",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "concat --naive failed: {err}");
    assert_eq!(
        out.lines().filter(|l| l.starts_with("#CHROM\t")).count(),
        1,
        "naive output should contain exactly one header line:\n{out}"
    );
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(records.len(), 4);
    assert_eq!(records[0], "1\t100\t.\tA\tC\t100\tPASS\t.");
    assert_eq!(records[3], "2\t20\t.\tC\tG\t100\tPASS\t.");
}

#[test]
fn concat_naive_rejects_different_headers_without_force() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_DIFF_SAMPLES);
    let (_out, err, code) = run(&["concat", "-n", a.to_str().unwrap(), b.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("Different headers"),
        "expected header mismatch error, got: {err}"
    );
}

#[test]
fn concat_naive_force_allows_different_headers() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_DIFF_SAMPLES);
    let (out, err, code) = run(&[
        "concat",
        "--naive-force",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "concat --naive-force failed: {err}");
    assert!(
        out.lines()
            .any(|l| l.starts_with("1\t100\t.\tA\tC\t100\tPASS\t.\tGT\t0/1")),
        "missing forced second record:\n{out}"
    );
}

#[test]
fn concat_regions_restricts_records_by_position() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let (out, err, code) = run(&[
        "concat",
        "-r",
        "1:150-250,2:1-15",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "concat -r failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t200\t.\tG\tT\t100\tPASS\t.",
            "2\t10\t.\tT\tA\t100\tPASS\t.",
        ]
    );
}

#[test]
fn concat_regions_file_supports_bed_coordinates() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_B);
    let bed = dir.path().join("regions.bed");
    std::fs::write(&bed, "1\t199\t200\n2\t9\t10\n").unwrap();
    let (out, err, code) = run(&[
        "concat",
        "-R",
        bed.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "concat -R BED failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t200\t.\tG\tT\t100\tPASS\t.",
            "2\t10\t.\tT\tA\t100\tPASS\t.",
        ]
    );
}

#[test]
fn concat_regions_overlap_record_includes_spanning_deletion() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tACGTACGTAC\tA\t100\tPASS\t.\n\
1\t200\t.\tG\tT\t100\tPASS\t.\n";
    let a = write_temp(&dir, "a.vcf", body);
    let (out, err, code) = run(&[
        "concat",
        "-r",
        "1:105-105",
        "--regions-overlap",
        "record",
        a.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "concat --regions-overlap record failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(records, ["1\t100\t.\tACGTACGTAC\tA\t100\tPASS\t."]);
}

#[test]
fn concat_regions_overlap_pos_excludes_spanning_deletion_by_default() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tACGTACGTAC\tA\t100\tPASS\t.\n";
    let a = write_temp(&dir, "a.vcf", body);
    let (out, err, code) = run(&["concat", "-r", "1:105-105", a.to_str().unwrap()]);
    assert_eq!(code, 0, "concat default regions-overlap failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(records.is_empty(), "unexpected records: {records:?}");
}

#[test]
fn concat_drop_genotypes_strips_format_and_samples() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_DIFF_SAMPLES);
    let (out, err, code) = run(&["concat", "-G", a.to_str().unwrap()]);
    assert_eq!(code, 0, "concat -G failed: {err}");
    let header_line = out.lines().find(|l| l.starts_with("#CHROM\t")).unwrap();
    assert!(
        !header_line.contains("\tFORMAT\t"),
        "FORMAT column not dropped: {header_line}"
    );
    assert!(
        !header_line.contains("\tNA1"),
        "sample column not dropped: {header_line}"
    );
    let record = out
        .lines()
        .find(|l| !l.starts_with('#') && !l.is_empty())
        .unwrap();
    let cols: Vec<&str> = record.split('\t').collect();
    assert_eq!(cols.len(), 8);
}

#[test]
fn concat_rejects_inputs_with_different_sample_columns() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_DIFF_SAMPLES);
    let b = write_temp(&dir, "b.vcf", VCF_DIFF_SAMPLES_B);
    let (_out, err, code) = run(&["concat", a.to_str().unwrap(), b.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("Different sample columns"),
        "expected sample-mismatch error, got: {err}"
    );
}

#[test]
fn concat_rm_dups_exact_drops_byte_identical_duplicates() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let b = write_temp(&dir, "b.vcf", VCF_A);
    let (out, err, code) = run(&["concat", "-D", a.to_str().unwrap(), b.to_str().unwrap()]);
    assert_eq!(code, 0, "concat -D failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "1\t100\t.\tA\tC\t100\tPASS\t.",
            "1\t200\t.\tG\tT\t100\tPASS\t."
        ]
    );
}

#[test]
fn concat_no_args_errors_with_usage() {
    let (_out, err, code) = run(&["concat"]);
    assert_ne!(code, 0);
    assert!(err.contains("Usage:"), "usage missing in stderr: {err}");
}

#[test]
fn concat_write_index_requires_output_file() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let (_out, err, code) = run(&["concat", "-W", a.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("-W requires an output file"),
        "expected output-file error, got: {err}"
    );
}

#[test]
fn concat_unknown_output_type_errors() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let (_out, err, code) = run(&["concat", a.to_str().unwrap(), "-O", "Q"]);
    assert_ne!(code, 0);
    assert!(err.contains("not recognised"), "got: {err}");
}

#[test]
fn concat_threads_rejects_non_integer_argument() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let (_out, err, code) = run(&["concat", "--threads", "abc", a.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("Could not parse argument: --threads abc"),
        "got: {err}"
    );
}

#[test]
fn concat_regions_overlap_rejects_unknown_mode() {
    let dir = TempDir::new().unwrap();
    let a = write_temp(&dir, "a.vcf", VCF_A);
    let (_out, err, code) = run(&["concat", "--regions-overlap", "bad", a.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("Could not parse --regions-overlap bad"),
        "got: {err}"
    );
}
