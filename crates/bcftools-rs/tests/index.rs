//! End-to-end tests for `bcftools_rs::commands::index`.

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
    let p = bin_path();
    if !p.exists() {
        let status = Command::new(env!("CARGO"))
            .args(["build", "-p", "bcftools-rs-cli"])
            .status()
            .expect("cargo build");
        assert!(status.success(), "failed to build bcftools-rs-cli");
    }
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
