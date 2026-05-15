//! End-to-end tests for `bcftools_rs::commands::index`.

use std::path::PathBuf;
use std::process::{Command, Stdio};

use tempfile::TempDir;

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

/// A small fully sorted VCF — sorted both by contig declaration order and by
/// position within each contig. Suitable for CSI indexing.
const SORTED_VCF: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##contig=<ID=2,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tA\tG\t100\tPASS\t.\n\
1\t20\t.\tC\tT\t100\tPASS\t.\n\
1\t30\t.\tG\tA\t100\tPASS\t.\n\
2\t15\t.\tT\tA\t100\tPASS\t.\n\
2\t25\t.\tA\tT\t100\tPASS\t.\n";

const IDX_STYLE_VCF: &str = "##fileformat=VCFv4.1\n\
##contig=<ID=1,length=249250621>\n\
##contig=<ID=11,length=135006516>\n\
##contig=<ID=20,length=63025520>\n\
##contig=<ID=X,length=155270560>\n\
##contig=<ID=Y,length=59373566>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
11\t2343543\t.\tA\t.\t999\tPASS\t.\n\
11\t5464562\t.\tC\tT\t999\tPASS\t.\n\
20\t76962\t.\tT\tC\t999\tPASS\t.\n\
20\t126310\t.\tACC\tA\t999\tPASS\t.\n\
20\t138125\t.\tG\tT\t999\tPASS\t.\n\
20\t138148\t.\tC\tT\t999\tPASS\t.\n\
20\t271225\t.\tT\tTTTA,TA\t999\tPASS\t.\n\
20\t304568\t.\tC\tT\t999\tPASS\t.\n\
20\t326891\t.\tA\tAC\t999\tPASS\t.\n\
X\t2928329\t.\tC\tT\t999\tPASS\t.\n\
X\t2933066\t.\tG\tC\t999\tPASS\t.\n\
X\t2942109\t.\tT\tC\t999\tPASS\t.\n\
X\t3048719\t.\tT\tC\t999\tPASS\t.\n\
Y\t8657215\t.\tC\tA\t999\tPASS\t.\n\
Y\t10011673\t.\tG\tA\t999\tPASS\t.\n";

#[test]
fn index_bcf_writes_csi() {
    let dir = TempDir::new().expect("tempdir");
    let vcf_path = dir.path().join("sorted.vcf");
    let bcf_path = dir.path().join("sorted.bcf");
    std::fs::write(&vcf_path, SORTED_VCF).unwrap();

    // Convert the inline sorted VCF to BCF using our own `view`.
    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf_path.to_str().unwrap(),
        vcf_path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");
    assert!(bcf_path.exists());

    // Now build a CSI index for it.
    let (_out, err, code) = run(&["index", "-f", bcf_path.to_str().unwrap()]);
    assert_eq!(code, 0, "index failed: {err}");

    let csi_path = dir.path().join("sorted.bcf.csi");
    assert!(csi_path.exists(), "csi not produced");

    let bytes = std::fs::read(&csi_path).expect("read csi");
    assert!(!bytes.is_empty(), "csi file is empty");
}

#[test]
fn index_vcf_gz_writes_csi_by_default() {
    // Build a real bgzf-encoded VCF and index it via `bcftools index`.
    // Mirrors `bcftools index -c in.vcf.gz` upstream.
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = dir.path().join("sorted.vcf.gz");
    let bytes = bgzf_encode(SORTED_VCF.as_bytes());
    std::fs::write(&vcf_gz, &bytes).unwrap();

    let (_out, err, code) = run(&["index", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index failed: {err}");

    let csi_path = dir.path().join("sorted.vcf.gz.csi");
    assert!(csi_path.exists(), "csi index not produced");
    let csi_bytes = std::fs::read(&csi_path).expect("read csi");
    assert!(!csi_bytes.is_empty());
}

#[test]
fn index_vcf_gz_with_t_flag_writes_tbi() {
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = dir.path().join("sorted.vcf.gz");
    let bytes = bgzf_encode(SORTED_VCF.as_bytes());
    std::fs::write(&vcf_gz, &bytes).unwrap();

    let (_out, err, code) = run(&["index", "-t", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -t failed: {err}");

    let tbi_path = dir.path().join("sorted.vcf.gz.tbi");
    assert!(tbi_path.exists(), "tbi index not produced");
    let tbi_bytes = std::fs::read(&tbi_path).expect("read tbi");
    assert!(!tbi_bytes.is_empty());
    let mut decoder = flate2::read::MultiGzDecoder::new(&tbi_bytes[..]);
    let mut head = [0u8; 4];
    use std::io::Read as _;
    decoder.read_exact(&mut head).expect("decode tbi head");
    assert_eq!(&head, b"TBI\x01", "wrong tbi magic: {head:?}");
}

#[test]
fn index_vcf_gz_does_not_rewrite_input() {
    // The new htslib-rs build_vcf_*_from_path helpers do not re-encode the
    // input. Verify the input bytes are byte-for-byte identical after `index`.
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = dir.path().join("sorted.vcf.gz");
    let bytes = bgzf_encode(SORTED_VCF.as_bytes());
    std::fs::write(&vcf_gz, &bytes).unwrap();

    let before = std::fs::read(&vcf_gz).expect("read before");
    let (_out, err, code) = run(&["index", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index failed: {err}");
    let after = std::fs::read(&vcf_gz).expect("read after");

    assert_eq!(before, after, "index unexpectedly rewrote input");
}

#[test]
fn index_stats_reports_per_contig_counts_for_vcf_gz_indexes() {
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = dir.path().join("sorted.vcf.gz");
    let bytes = bgzf_encode(SORTED_VCF.as_bytes());
    std::fs::write(&vcf_gz, &bytes).unwrap();

    let (_out, err, code) = run(&["index", "-t", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -t failed: {err}");

    let (out, err, code) = run(&["index", "-s", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -s failed: {err}");
    assert_eq!(out, "1\t1000\t3\n2\t1000\t2\n");
}

#[test]
fn index_nrecords_reports_total_for_data_or_index_path() {
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = dir.path().join("sorted.vcf.gz");
    let bytes = bgzf_encode(SORTED_VCF.as_bytes());
    std::fs::write(&vcf_gz, &bytes).unwrap();

    let (_out, err, code) = run(&["index", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index failed: {err}");

    let (out, err, code) = run(&["index", "-n", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -n failed: {err}");
    assert_eq!(out, "5\n");

    let csi_path = dir.path().join("sorted.vcf.gz.csi");
    let (out, err, code) = run(&["index", "-n", csi_path.to_str().unwrap()]);
    assert_eq!(code, 0, "index -n index-path failed: {err}");
    assert_eq!(out, "5\n");
}

#[test]
fn index_can_build_from_stdin_when_output_path_is_given() {
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = bgzf_encode(SORTED_VCF.as_bytes());
    let csi = dir.path().join("streamed.vcf.gz.csi");

    let (_out, err, code) = run_with_stdin(&["index", "-f", "-o", csi.to_str().unwrap()], &vcf_gz);
    assert_eq!(code, 0, "index stdin failed: {err}");
    assert!(csi.exists(), "CSI not produced from stdin");

    let (out, err, code) = run(&["index", "-n", csi.to_str().unwrap()]);
    assert_eq!(code, 0, "index -n on stdin-built CSI failed: {err}");
    assert_eq!(out, "5\n");
}

#[test]
fn index_stats_match_upstream_idx_fixture_when_header_has_empty_leading_contig() {
    let dir = TempDir::new().expect("tempdir");
    let vcf_gz = dir.path().join("idx.vcf.gz");
    let bytes = bgzf_encode(IDX_STYLE_VCF.as_bytes());
    std::fs::write(&vcf_gz, &bytes).unwrap();

    let (_out, err, code) = run(&["index", "-t", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -t failed: {err}");

    let (out, err, code) = run(&["index", "-s", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -s failed: {err}");
    assert_eq!(
        out,
        "11\t135006516\t2\n20\t63025520\t7\nX\t155270560\t4\nY\t59373566\t2\n"
    );

    let (out, err, code) = run(&["index", "-n", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index -n failed: {err}");
    assert_eq!(out, "15\n");
}

/// Encode `data` as a BGZF stream using htslib-rs's bgzf writer.
fn bgzf_encode(data: &[u8]) -> Vec<u8> {
    use std::io::Write as _;
    let mut writer = htslib_rs::bgzf::io::Writer::new(Vec::new());
    writer.write_all(data).unwrap();
    writer.finish().unwrap()
}

#[test]
fn index_min_shift_out_of_range_errors() {
    // No need to produce a real BCF — the option-validation path errors before
    // any I/O happens.
    let (_out, err, code) = run(&["index", "-f", "-m", "99", "/tmp/does-not-matter.bcf"]);
    assert_ne!(code, 0);
    assert!(err.contains("expected min_shift in range"), "stderr: {err}");
}

#[test]
fn index_threads_rejects_non_integer_argument() {
    let (_out, err, code) = run(&["index", "--threads", "abc", "/tmp/does-not-matter.bcf"]);
    assert_eq!(code, 255);
    assert!(err.contains("Could not parse argument: --threads abc"));
}

#[test]
fn index_rejects_extra_input_paths() {
    let (_out, err, code) = run(&[
        "index",
        "/tmp/first-does-not-matter.vcf.gz",
        "/tmp/second-does-not-matter.vcf.gz",
    ]);
    assert_ne!(code, 0);
    assert!(
        err.contains("multiple input files are not supported"),
        "stderr: {err}"
    );
}

#[test]
fn index_and_view_large_chromosome_fixture_matches_upstream_output() {
    let dir = TempDir::new().expect("tempdir");
    let input = fixture_path("large_chrom_csi_limit.vcf");
    let expected =
        std::fs::read_to_string(fixture_path("large_chrom_csi_limit.20.1.2147483647.out")).unwrap();

    let vcf_gz = dir.path().join("large_chrom_csi_limit.vcf.gz");
    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Oz",
        "-o",
        vcf_gz.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Oz failed: {err}");

    let (_out, err, code) = run(&["index", "-f", vcf_gz.to_str().unwrap()]);
    assert_eq!(code, 0, "index VCF.gz failed: {err}");

    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        vcf_gz.to_str().unwrap(),
        "chr20:1-2147483647",
    ]);
    assert_eq!(code, 0, "view VCF.gz region failed: {err}");
    assert_eq!(out, expected);

    let bcf = dir.path().join("large_chrom_csi_limit.bcf");
    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (_out, err, code) = run(&["index", "-f", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "index BCF failed: {err}");

    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        bcf.to_str().unwrap(),
        "chr20:1-2147483647",
    ]);
    assert_eq!(code, 0, "view BCF region failed: {err}");
    assert_eq!(out, expected);
}
