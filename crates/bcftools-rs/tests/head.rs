//! End-to-end tests for `bcftools_rs::commands::head`.
//!
//! Tests build the `bcftools` binary from `crates/bcftools-rs-cli` and
//! invoke it as a subprocess to validate output. This mirrors how the
//! upstream Perl parity gate exercises the CLI and avoids fighting the
//! cargo test harness over stdout.

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
    // Per `cargo test` convention, sibling binaries live in the same target dir.
    // CARGO_BIN_EXE_<name> is set when the test crate depends (dev or build) on
    // the binary's package; we don't, so derive the path manually from
    // CARGO_MANIFEST_DIR.
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
    assert!(
        bin_path().exists(),
        "binary not at expected path: {}",
        bin_path().display()
    );
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
fn head_default_prints_full_header_no_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(out.starts_with("##fileformat=VCFv"), "got: {out:?}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"));
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(
        record_lines.is_empty(),
        "unexpected records: {record_lines:?}"
    );
}

#[test]
fn head_with_n2_emits_two_records_after_header() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", "-n", "2", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(record_lines.len(), 2);
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn head_with_h_truncates_header_lines() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", "-h", "3", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let header_lines: Vec<_> = out.lines().collect();
    assert_eq!(header_lines.len(), 3);
    assert!(header_lines[0].starts_with("##fileformat="));
}

#[test]
fn head_with_s_emits_chrom_line_then_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["head", "-s", "1", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let lines: Vec<_> = out.lines().collect();
    assert!(lines.iter().any(|l| l.starts_with("#CHROM\t")));
    let record_lines: Vec<_> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(record_lines.len(), 1);
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn head_reads_plain_vcf_from_stdin_without_filename() {
    let input = std::fs::read(fixture_path("mpileup.2.vcf")).unwrap();
    let (out, err, code) = run_with_stdin(&["head", "-s2", "-h2"], &input);
    assert_eq!(code, 0, "head stdin failed: {err}");
    assert!(out.starts_with("##fileformat=VCFv4.2\n"));
    assert!(
        out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n")
    );
    let records: Vec<_> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(records.len(), 2);
    assert!(records[0].starts_with("chr1\t212740\t"));
}

#[test]
fn head_reads_bcf_from_stdin_pipe() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("mpileup.2.vcf");
    let bcf = tmp.path().join("mpileup.2.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let data = std::fs::read(&bcf).unwrap();
    let (out, err, code) = run_with_stdin(&["head", "-s1"], &data);
    assert_eq!(code, 0, "head BCF stdin failed: {err}");
    let expected = std::fs::read_to_string(fixture_path("head.2.out")).unwrap();
    assert_eq!(out, expected);
}

#[test]
fn head_reads_bgzf_vcf_file_with_records() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("mpileup.2.vcf");
    let compressed = tmp.path().join("mpileup.2.vcf.gz");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Oz",
        "-o",
        compressed.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Oz failed: {err}");

    let (out, err, code) = run(&["head", "-n", "2", compressed.to_str().unwrap()]);
    assert_eq!(code, 0, "head -n on BGZF VCF failed: {err}");
    assert!(out.starts_with("##fileformat=VCFv4.2\n"));
    let records: Vec<_> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(records.len(), 2);
    assert!(records[0].starts_with("chr1\t212740\t"));
}

#[test]
fn head_with_s_accepts_kestrel_non_canonical_fileformat_header() {
    let dir = tempfile::TempDir::new().expect("tempdir");
    let path = dir.path().join("kestrel.vcf");
    let kestrel = "##fileformat=VCF4.2\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\tC\t.\tPASS\t.\n";
    std::fs::write(&path, kestrel).unwrap();

    let (out, err, code) = run(&["head", "-s", "1", path.to_str().unwrap()]);
    assert_eq!(code, 0, "head -s rejected Kestrel header: {err}");
    assert!(
        err.contains("[W::bcf_get_version] Couldn't get VCF version, considering as 4.2"),
        "missing upstream-style warning: {err}"
    );
    let lines: Vec<_> = out.lines().collect();
    assert!(lines.iter().any(|l| l.starts_with("#CHROM\t")));
    let records: Vec<&str> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .copied()
        .collect();
    assert_eq!(records, ["chr1\t1\t.\tA\tC\t.\tPASS\t."]);
}

#[test]
fn version_flag_prints_block() {
    let (out, _err, code) = run(&["--version"]);
    assert_eq!(code, 0);
    assert!(out.contains("bcftools "));
    assert!(out.contains("htslib "));
}

#[test]
fn version_only_one_line() {
    let (out, _err, code) = run(&["--version-only"]);
    assert_eq!(code, 0);
    assert!(out.contains("+htslib-"));
    assert_eq!(out.lines().count(), 1);
}

#[test]
fn help_lists_upstream_dispatch_sections() {
    let (out, _err, code) = run(&["--help"]);
    assert_eq!(code, 0);
    assert!(out.contains(" -- Indexing"));
    assert!(out.contains(" -- VCF/BCF manipulation"));
    assert!(out.contains(" -- VCF/BCF analysis"));
    assert!(out.contains(" -- Plugins"));
    assert!(out.contains("    annotate"));
    assert!(out.contains("    mpileup"));
    assert!(out.contains("    plugin"));
    assert!(out.contains("41 plugins available"));
    assert!(!out.contains("    tabix"));
    assert!(!out.contains("    som"));
}

#[test]
fn plugin_lists_static_registry() {
    let (out, err, code) = run(&["plugin", "-l"]);
    assert_eq!(code, 0, "plugin -l failed: {err}");
    let names: Vec<_> = out.lines().collect();
    assert_eq!(names.len(), 41);
    assert!(names.contains(&"fill-tags"));
    assert!(names.contains(&"missing2ref"));
    assert!(names.contains(&"trio-dnm2"));
}

#[test]
fn plugin_verbose_list_includes_descriptions() {
    let (out, err, code) = run(&["plugin", "-lv"]);
    assert_eq!(code, 0, "plugin -lv failed: {err}");
    assert!(out.contains("-- fill-tags --"));
    assert!(out.contains("Fill INFO tags"));
    assert!(out.contains("-- split-vep --"));
}

#[test]
fn plugin_shortcut_help_uses_registry() {
    let (_out, err, code) = run(&["+fill-tags", "--help"]);
    assert_eq!(code, 0, "+fill-tags --help failed: {err}");
    assert!(err.contains("About:   Fill INFO tags"));
    assert!(err.contains("registered but its record-processing implementation is not yet ported"));
}

#[test]
fn plugin_known_name_without_implementation_errors() {
    let (_out, err, code) = run(&["+fill-tags"]);
    assert_ne!(code, 0);
    assert!(err.contains("plugin 'fill-tags' is registered but not yet implemented"));
}

#[test]
fn unknown_subcommand_errors() {
    let (_out, err, code) = run(&["bogus"]);
    assert_ne!(code, 0);
    assert!(err.contains("unrecognized command 'bogus'"));
}
