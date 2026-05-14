//! End-to-end tests for `bcftools_rs::commands::convert`.

use std::io::{Read as _, Write as _};
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

fn fixture_path(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("bcftools");
    p.push("test");
    p.push(name);
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

fn without_meta_headers(text: &str) -> String {
    text.lines()
        .filter(|line| !line.starts_with("##"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

#[test]
fn convert_tsv2vcf_explicit_columns_writes_vcf_records() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    std::fs::write(&tsv, "#comment\nchr1\t10\trs1\tA\tC,G\nchr2\t5\t.\tT\tA\n").unwrap();
    let (out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,ID,REF,ALT",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --tsv2vcf failed: {err}");
    assert!(out.starts_with("##fileformat=VCFv4.2\n"));
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t10\trs1\tA\tC,G\t.\t.\t.",
            "chr2\t5\t.\tT\tA\t.\t.\t."
        ]
    );
}

#[test]
fn convert_tsv2vcf_writes_bgzf_and_index() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    let out_path = dir.path().join("out.vcf.gz");
    std::fs::write(&tsv, "chr1\t10\trs1\tA\tC\n").unwrap();
    let (_out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,ID,REF,ALT",
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "-W",
        "--threads",
        "2",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert -Oz -W failed: {err}");
    assert!(dir.path().join("out.vcf.gz.csi").exists());
    let mut decoded = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(out_path).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    assert!(decoded.contains("chr1\t10\trs1\tA\tC\t.\t.\t."));
}

#[test]
fn convert_tsv2vcf_writes_bcf_and_index() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    let out_path = dir.path().join("out.bcf");
    std::fs::write(&tsv, "chr1\t10\trs1\tA\tC\n").unwrap();
    let (_out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,ID,REF,ALT",
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "b",
        "-W",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert -Ob -W failed: {err}");
    assert!(out_path.exists());
    assert!(dir.path().join("out.bcf.csi").exists());

    let (view_out, view_err, view_code) =
        run(&["view", "--no-version", out_path.to_str().unwrap()]);
    assert_eq!(view_code, 0, "view of BCF failed: {view_err}");
    assert!(view_out.contains("chr1\t10\trs1\tA\tC\t.\t.\t."));
}

#[test]
fn convert_tsv2vcf_samples_emit_gt_columns_from_vcf_style_values() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    std::fs::write(&tsv, "chr1\t10\trs1\tA\tC\t0/1\t1/1\n").unwrap();
    let (out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,ID,REF,ALT",
        "-s",
        "S1,S2",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --tsv2vcf -s failed: {err}");
    assert!(out.contains("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"));
    assert!(out.contains("chr1\t10\trs1\tA\tC\t.\t.\t.\tGT\t0/1\t1/1"));
}

#[test]
fn convert_tsv2vcf_samples_file_converts_allele_letter_pairs_to_gt() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    let samples = dir.path().join("samples.txt");
    std::fs::write(&tsv, "chr1\t10\trs1\tA\tC,G\tAC\tGG\n").unwrap();
    std::fs::write(&samples, "S1\nS2\n").unwrap();
    let (out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,ID,REF,ALT",
        "-S",
        samples.to_str().unwrap(),
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --tsv2vcf -S failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"));
    assert!(out.contains("chr1\t10\trs1\tA\tC,G\t.\t.\t.\tGT\t0/1\t2/2"));
}

#[test]
fn convert_tsv2vcf_ignored_columns_drop_ref_matching_alt() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let tsv = dir.path().join("in.tsv");
    std::fs::write(&fasta, ">1\nAAAAAAAAAAAA\n").unwrap();
    std::fs::write(&tsv, "rs001\t1\t2\tA   A\nrs002\t1\t10\tA   G\n").unwrap();

    let (out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "-,CHROM,POS,REF,ALT",
        "-f",
        fasta.to_str().unwrap(),
        "--no-version",
    ]);
    assert_eq!(
        code, 0,
        "convert --tsv2vcf ignored-column input failed: {err}"
    );
    assert!(out.contains("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"));
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert!(out.contains("1\t2\t.\tA\t.\t.\t.\t."));
    assert!(out.contains("1\t10\t.\tA\tG\t.\t.\t."));
}

#[test]
fn convert_tsv2vcf_skips_malformed_rows_and_keeps_writing() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    std::fs::write(
        &tsv,
        "chr1\t10\trs1\tA\tC\nchr1\tbad\trs_bad\tA\tC\nchr1\t20\trs2\tG\tT\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,ID,REF,ALT",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert should skip malformed rows: {err}");
    assert!(out.contains("chr1\t10\trs1\tA\tC\t.\t.\t."));
    assert!(out.contains("chr1\t20\trs2\tG\tT\t.\t.\t."));
    assert!(!out.contains("rs_bad"));
    assert!(
        err.contains("Warning: skipping malformed TSV line 2"),
        "got: {err}"
    );
    assert!(err.contains("Rows total: \t3"), "got: {err}");
    assert!(err.contains("Rows skipped: \t1"), "got: {err}");
    assert!(err.contains("Sites written: \t2"), "got: {err}");
}

#[test]
fn convert_tsv2vcf_aa_with_reference_derives_ref_alt_and_gt() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let tsv = dir.path().join("in.tsv");
    std::fs::write(&fasta, ">chr1\nAGTAC\n>chr2\nCCCC\n").unwrap();
    std::fs::write(
        &tsv,
        "rs1\tchr1\t2\tAG\nrs2\tchr1\t3\tTT\nrs3\tchr2\t1\t--\nrs4\tchr1\t4\tI\n",
    )
    .unwrap();
    let (out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-f",
        fasta.to_str().unwrap(),
        "-c",
        "ID,CHROM,POS,AA",
        "-s",
        "SAMPLE1",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --tsv2vcf AA failed: {err}");
    assert!(out.contains("##contig=<ID=chr1,length=5>"));
    assert!(out.contains("##contig=<ID=chr2,length=4>"));
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"));
    assert!(out.contains("chr1\t2\trs1\tG\tA\t.\t.\t.\tGT\t1/0"));
    assert!(out.contains("chr1\t3\trs2\tT\t.\t.\t.\t.\tGT\t0/0"));
    assert!(out.contains("chr2\t1\trs3\tC\t.\t.\t.\t.\tGT\t./."));
    assert!(!out.contains("rs4"));
    assert!(err.contains("Rows total: \t4"), "got: {err}");
    assert!(err.contains("Rows skipped: \t1"), "got: {err}");
    assert!(err.contains("Sites written: \t3"), "got: {err}");
    assert!(err.contains("Missing GTs: \t1"), "got: {err}");
    assert!(err.contains("Hom RR: \t1"), "got: {err}");
    assert!(err.contains("Het RA: \t1"), "got: {err}");
    assert!(err.contains("Hom AA: \t0"), "got: {err}");
    assert!(err.contains("Het AA: \t0"), "got: {err}");
}

#[test]
fn convert_tsv2vcf_aa_requires_reference() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    std::fs::write(&tsv, "rs1\tchr1\t10\tA\n").unwrap();
    let (_out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "ID,CHROM,POS,AA",
    ]);
    assert_ne!(code, 0);
    assert!(
        err.contains("--tsv2vcf requires the --fasta-ref option when AA is used"),
        "got: {err}"
    );
}

#[test]
fn convert_gvcf2vcf_expands_reference_blocks_with_fasta_ref() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let gvcf = dir.path().join("in.g.vcf");
    std::fs::write(&fasta, ">chr1\nACGTAC\n").unwrap();
    std::fs::write(
        &gvcf,
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##contig=<ID=chr1,length=6>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t2\t.\tC\t<NON_REF>\t.\tPASS\tEND=4;BLOCK\tGT:DP\t0/0:5\n\
chr1\t5\tvar\tA\tT\t50\tPASS\t.\tGT:DP\t0/1:8\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "convert",
        "--gvcf2vcf",
        "-f",
        fasta.to_str().unwrap(),
        "--no-version",
        gvcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --gvcf2vcf failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t2\t.\tC\t<NON_REF>\t.\tPASS\tBLOCK\tGT:DP\t0/0:5",
            "chr1\t3\t.\tG\t<NON_REF>\t.\tPASS\tBLOCK\tGT:DP\t0/0:5",
            "chr1\t4\t.\tT\t<NON_REF>\t.\tPASS\tBLOCK\tGT:DP\t0/0:5",
            "chr1\t5\tvar\tA\tT\t50\tPASS\t.\tGT:DP\t0/1:8",
        ]
    );
    assert!(
        !out.contains("END=4"),
        "expanded rows should drop INFO/END:\n{out}"
    );
}

#[test]
fn convert_gvcf2vcf_reads_gzip_input() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let gvcf = dir.path().join("in.g.vcf.gz");
    std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();
    let input = "##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\t.\t.\tPASS\tEND=2\n";
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    encoder.write_all(input.as_bytes()).unwrap();
    std::fs::write(&gvcf, encoder.finish().unwrap()).unwrap();

    let (out, err, code) = run(&[
        "convert",
        "--gvcf2vcf",
        gvcf.to_str().unwrap(),
        "-f",
        fasta.to_str().unwrap(),
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --gvcf2vcf .gz failed: {err}");
    assert!(out.contains("chr1\t1\t.\tA\t.\t.\tPASS\t."));
    assert!(out.contains("chr1\t2\t.\tC\t.\t.\tPASS\t."));
}

#[test]
fn convert_gvcf2vcf_reads_bcf_input() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let gvcf = dir.path().join("in.g.vcf");
    let bcf = dir.path().join("in.g.bcf");
    std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();
    std::fs::write(
        &gvcf,
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=2\n",
    )
    .unwrap();
    let (_view_out, view_err, view_code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        gvcf.to_str().unwrap(),
    ]);
    assert_eq!(view_code, 0, "view -Ob failed: {view_err}");

    let (out, err, code) = run(&[
        "convert",
        "--gvcf2vcf",
        bcf.to_str().unwrap(),
        "-f",
        fasta.to_str().unwrap(),
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --gvcf2vcf BCF failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\t.",
            "chr1\t2\t.\tC\t<NON_REF>\t.\tPASS\t.",
        ]
    );
}

#[test]
fn convert_gvcf2vcf_reads_bcf_from_stdin() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let gvcf = dir.path().join("in.g.vcf");
    let bcf = dir.path().join("in.g.bcf");
    std::fs::write(&fasta, ">chr1\nACGT\n").unwrap();
    std::fs::write(
        &gvcf,
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##contig=<ID=chr1,length=4>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\tEND=2\n",
    )
    .unwrap();
    let (_view_out, view_err, view_code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        gvcf.to_str().unwrap(),
    ]);
    assert_eq!(view_code, 0, "view -Ob failed: {view_err}");

    ensure_binary_built();
    let mut child = Command::new(bin_path())
        .args([
            "convert",
            "--gvcf2vcf",
            "-f",
            fasta.to_str().unwrap(),
            "--no-version",
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn bcftools convert");
    child
        .stdin
        .as_mut()
        .unwrap()
        .write_all(&std::fs::read(&bcf).unwrap())
        .unwrap();
    let out = child.wait_with_output().unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "convert stdin BCF failed: {stderr}");
    assert!(stdout.contains("chr1\t1\t.\tA\t<NON_REF>\t.\tPASS\t."));
    assert!(stdout.contains("chr1\t2\t.\tC\t<NON_REF>\t.\tPASS\t."));
}

#[test]
fn convert_gvcf2vcf_filters_before_expansion_and_keeps_failing_records() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    let gvcf = dir.path().join("in.g.vcf");
    std::fs::write(&fasta, ">chr1\nACGTAC\n").unwrap();
    std::fs::write(
        &gvcf,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##contig=<ID=chr1,length=6>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t1\t.\tA\t<NON_REF>\t.\tLowGQX\tEND=2;DP=9\n\
chr1\t3\t.\tG\t<NON_REF>\t.\tPASS\tEND=4;DP=9\n\
chr1\t5\tvar\tA\tT\t50\tPASS\tDP=3\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "convert",
        "--gvcf2vcf",
        gvcf.to_str().unwrap(),
        "-f",
        fasta.to_str().unwrap(),
        "-i",
        "FILTER=\"PASS\"",
        "--no-version",
    ]);
    assert_eq!(code, 0, "convert --gvcf2vcf filtered failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t1\t.\tA\t<NON_REF>\t.\tLowGQX\tEND=2;DP=9",
            "chr1\t3\t.\tG\t<NON_REF>\t.\tPASS\tDP=9",
            "chr1\t4\t.\tT\t<NON_REF>\t.\tPASS\tDP=9",
            "chr1\t5\tvar\tA\tT\t50\tPASS\tDP=3",
        ]
    );
}

#[test]
fn convert_gvcf2vcf_matches_upstream_filtered_fixture() {
    let gvcf = fixture_path("convert.gvcf.vcf");
    let fasta = fixture_path("gvcf.fa");
    let expected = std::fs::read_to_string(fixture_path("convert.gvcf.out")).unwrap();

    let (out, err, code) = run(&[
        "convert",
        "--gvcf2vcf",
        "-i",
        "FILTER=\"PASS\"",
        "-f",
        fasta.to_str().unwrap(),
        "--no-version",
        gvcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "fixture --gvcf2vcf failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn convert_hapsample_writes_hap_gz_and_samples_from_vcf() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let prefix = dir.path().join("shapeit");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t2\trs1\tC\tT\t.\tPASS\tDP=8\tGT\t0|1\t1/1\n\
chr1\t3\trs2\tG\tA,C\t.\tPASS\tDP=8\tGT\t0/1\t0/2\n\
chr1\t4\trs3\tT\t.\t.\tPASS\tDP=8\tGT\t0/0\t0/0\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        prefix.to_str().unwrap(),
        "-i",
        "DP>=8",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --hapsample failed: {err}");
    assert!(err.contains("1 records written, 2 skipped: 1/1/0 no-ALT/non-biallelic/filtered"));

    let samples = std::fs::read_to_string(prefix.with_extension("samples")).unwrap();
    assert_eq!(samples, "ID_1 ID_2 missing\n0 0 0\nS1 S1 0\nS2 S2 0\n");

    let mut hap = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("shapeit.hap.gz")).unwrap(),
    )
    .read_to_string(&mut hap)
    .unwrap();
    assert_eq!(hap, "chr1 chr1:2_C_T 2 C T 0 1 1* 1*\n");
}

#[test]
fn convert_hapsample_reads_bcf_and_supports_vcf_ids_and_explicit_outputs() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let bcf = dir.path().join("in.bcf");
    let hap = dir.path().join("out.hap");
    let samples = dir.path().join("out.samples");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t2\trs1\tC\tT\t.\tPASS\tDP=8\tGT\t0/1\n\
chr1\t5\trs2\tA\tG\t.\tPASS\tDP=2\tGT\t0|1\n",
    )
    .unwrap();
    let (_view_out, view_err, view_code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(view_code, 0, "view -Ob failed: {view_err}");

    let outputs = format!("{},{}", hap.display(), samples.display());
    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        &outputs,
        "--vcf-ids",
        "-e",
        "DP<5",
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --hapsample BCF failed: {err}");
    assert_eq!(
        std::fs::read_to_string(&samples).unwrap(),
        "ID_1 ID_2 missing\n0 0 0\nS1 S1 0\n"
    );
    assert_eq!(
        std::fs::read_to_string(&hap).unwrap(),
        "chr1:2_C_T rs1 2 C T 0* 1*\n"
    );
}

#[test]
fn convert_hapsample_matches_upstream_stdout_fixtures() {
    let vcf = fixture_path("convert.vcf");
    let expected_hap = std::fs::read_to_string(fixture_path("convert.hs.hap")).unwrap();
    let expected_ids_hap = std::fs::read_to_string(fixture_path("convert.hs.ids.hap")).unwrap();
    let expected_sample = std::fs::read_to_string(fixture_path("convert.hs.sample")).unwrap();

    let (hap_out, hap_err, hap_code) =
        run(&["convert", "--hapsample", "-,.", vcf.to_str().unwrap()]);
    assert_eq!(hap_code, 0, "fixture --hapsample -,. failed: {hap_err}");
    assert_eq!(hap_out, expected_hap);

    let (ids_hap_out, ids_hap_err, ids_hap_code) = run(&[
        "convert",
        "--hapsample",
        "-,.",
        "--vcf-ids",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(
        ids_hap_code, 0,
        "fixture --hapsample -,. --vcf-ids failed: {ids_hap_err}"
    );
    assert_eq!(ids_hap_out, expected_ids_hap);

    let (sample_out, sample_err, sample_code) =
        run(&["convert", "--hapsample", ".,-", vcf.to_str().unwrap()]);
    assert_eq!(
        sample_code, 0,
        "fixture --hapsample .,- failed: {sample_err}"
    );
    assert_eq!(sample_out, expected_sample);
}

#[test]
fn convert_hapsample2vcf_writes_vcf_from_hap_and_samples() {
    let dir = TempDir::new().unwrap();
    let hap = dir.path().join("in.hap");
    let samples = dir.path().join("in.samples");
    std::fs::write(
        &hap,
        "chr1 chr1:2_C_T 2 C T 0 1 1* 0*\nchr1 chr1:5_A_G_9 5 A G ? ? 1 -\n",
    )
    .unwrap();
    std::fs::write(&samples, "ID_1 ID_2 missing\n0 0 0\nA A 0\nB B 0\n").unwrap();
    let input = format!("{},{}", hap.display(), samples.display());
    let (out, err, code) = run(&["convert", "--hapsample2vcf", "--no-version", &input]);
    assert_eq!(code, 0, "convert --hapsample2vcf failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"));
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t2\t.\tC\tT\t.\t.\t.\tGT\t0|1\t1/0",
            "chr1\t5\t.\tA\tG\t.\t.\tEND=9\tGT\t.|.\t1"
        ]
    );
    assert!(err.contains("Number of processed rows: \t2"), "got: {err}");
}

#[test]
fn convert_hapsample2vcf_supports_vcf_ids_bcf_and_index() {
    let dir = TempDir::new().unwrap();
    let hap = dir.path().join("in.hap");
    let samples = dir.path().join("in.samples");
    let bcf = dir.path().join("out.bcf");
    std::fs::write(&hap, "chr1:2_C_T rs1 2 C T 0 1\n").unwrap();
    std::fs::write(&samples, "ID_1 ID_2 missing\n0 0 0\nA A 0\n").unwrap();
    let input = format!("{},{}", hap.display(), samples.display());
    let (_out, err, code) = run(&[
        "convert",
        "--hapsample2vcf",
        "--vcf-ids",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        "-W",
        &input,
    ]);
    assert_eq!(code, 0, "convert --hapsample2vcf -Ob failed: {err}");
    assert!(dir.path().join("out.bcf.csi").exists());
    let (view_out, view_err, view_code) = run(&["view", "--no-version", bcf.to_str().unwrap()]);
    assert_eq!(view_code, 0, "view of converted BCF failed: {view_err}");
    assert!(view_out.contains("chr1\t2\trs1\tC\tT\t.\t.\t.\tGT\t0|1"));
}

#[test]
fn convert_hapsample2vcf_matches_upstream_hap_sample_fixtures() {
    let hap = fixture_path("convert.hs.gt.hap");
    let ids_hap = fixture_path("convert.hs.gt.ids.hap");
    let samples = fixture_path("convert.hs.gt.samples");
    let expected = std::fs::read_to_string(fixture_path("convert.gt.noHead.vcf")).unwrap();
    let expected_ids = std::fs::read_to_string(fixture_path("convert.gt.noHead.ids.vcf")).unwrap();

    let input = format!("{},{}", hap.display(), samples.display());
    let (out, err, code) = run(&["convert", "--hapsample2vcf", "--no-version", &input]);
    assert_eq!(code, 0, "fixture --hapsample2vcf failed: {err}");
    assert_eq!(without_meta_headers(&out), expected);

    let ids_input = format!("{},{}", ids_hap.display(), samples.display());
    let (ids_out, ids_err, ids_code) = run(&[
        "convert",
        "--vcf-ids",
        "--hapsample2vcf",
        "--no-version",
        &ids_input,
    ]);
    assert_eq!(
        ids_code, 0,
        "fixture --vcf-ids --hapsample2vcf failed: {ids_err}"
    );
    assert_eq!(without_meta_headers(&ids_out), expected_ids);
}

#[test]
fn convert_hapsample_respects_sample_selection() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let samples_file = dir.path().join("samples.txt");
    let prefix = dir.path().join("subset");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0/0\n",
    )
    .unwrap();
    std::fs::write(&samples_file, "C\nA\n").unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        prefix.to_str().unwrap(),
        "-S",
        samples_file.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --hapsample -S failed: {err}");
    assert_eq!(
        std::fs::read_to_string(prefix.with_extension("samples")).unwrap(),
        "ID_1 ID_2 missing\n0 0 0\nA A 0\nC C 0\n"
    );

    let mut hap = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("subset.hap.gz")).unwrap(),
    )
    .read_to_string(&mut hap)
    .unwrap();
    assert_eq!(hap, "chr1 chr1:2_C_T 2 C T 0 1 0* 0*\n");
}

#[test]
fn convert_hapsample_supports_sample_exclusion() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let prefix = dir.path().join("exclude");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        prefix.to_str().unwrap(),
        "-s",
        "^B",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --hapsample -s ^ failed: {err}");
    assert_eq!(
        std::fs::read_to_string(prefix.with_extension("samples")).unwrap(),
        "ID_1 ID_2 missing\n0 0 0\nA A 0\nC C 0\n"
    );
    let mut hap = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("exclude.hap.gz")).unwrap(),
    )
    .read_to_string(&mut hap)
    .unwrap();
    assert_eq!(hap, "chr1 chr1:2_C_T 2 C T 0 1 0 0\n");
}

#[test]
fn convert_hapsample_haploid2diploid_duplicates_haploid_genotypes() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let default_prefix = dir.path().join("default");
    let diploid_prefix = dir.path().join("diploid");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tH0\tH1\tHM\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0\t1\t.\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        default_prefix.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --hapsample haploid failed: {err}");
    let mut default_hap = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("default.hap.gz")).unwrap(),
    )
    .read_to_string(&mut default_hap)
    .unwrap();
    assert_eq!(default_hap, "chr1 chr1:2_C_T 2 C T 0 - 1 - ? -\n");

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        diploid_prefix.to_str().unwrap(),
        "--haploid2diploid",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(
        code, 0,
        "convert --hapsample --haploid2diploid failed: {err}"
    );
    let mut diploid_hap = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("diploid.hap.gz")).unwrap(),
    )
    .read_to_string(&mut diploid_hap)
    .unwrap();
    assert_eq!(diploid_hap, "chr1 chr1:2_C_T 2 C T 0 0 1 1 ? ?\n");
}

#[test]
fn convert_hapsample_sex_file_adds_sample_sex_column() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let sex = dir.path().join("sex.txt");
    let prefix = dir.path().join("sexed");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1\t0|0\n",
    )
    .unwrap();
    std::fs::write(&sex, "A M\nC F\n").unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        prefix.to_str().unwrap(),
        "-s",
        "A,C",
        "--sex",
        sex.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --hapsample --sex failed: {err}");
    assert_eq!(
        std::fs::read_to_string(prefix.with_extension("samples")).unwrap(),
        "ID_1 ID_2 missing sex\n0 0 0 0\nA A 0 1\nC C 0 2\n"
    );
}

#[test]
fn convert_hapsample_sex_file_requires_selected_samples() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let sex = dir.path().join("sex.txt");
    let prefix = dir.path().join("missing-sex");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0|1\t1|1\n",
    )
    .unwrap();
    std::fs::write(&sex, "A M\n").unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--hapsample",
        prefix.to_str().unwrap(),
        "--sex",
        sex.to_str().unwrap(),
        vcf.to_str().unwrap(),
    ]);
    assert_ne!(code, 0);
    assert!(err.contains("Missing sex for sample B"), "got: {err}");
}

#[test]
fn convert_haplegendsample_writes_hap_legend_and_samples() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let sex = dir.path().join("sex.txt");
    let prefix = dir.path().join("oxford");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t2\trs1\tC\tT\t.\tPASS\tDP=8\tGT\t0|1\t1/1\n\
chr1\t3\trs2\tG\tA,C\t.\tPASS\tDP=8\tGT\t0/1\t0/2\n",
    )
    .unwrap();
    std::fs::write(&sex, "A M\nB F\n").unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--haplegendsample",
        prefix.to_str().unwrap(),
        "--sex",
        sex.to_str().unwrap(),
        "-i",
        "DP>=8",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --haplegendsample failed: {err}");
    assert!(err.contains("1 records written, 1 skipped: 0/1/0 no-ALT/non-biallelic/filtered"));
    assert_eq!(
        std::fs::read_to_string(prefix.with_extension("samples")).unwrap(),
        "sample population group sex\nA A A 1\nB B B 2\n"
    );

    let mut hap = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("oxford.hap.gz")).unwrap(),
    )
    .read_to_string(&mut hap)
    .unwrap();
    assert_eq!(hap, "0 1 1* 1*\n");

    let mut legend = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("oxford.legend.gz")).unwrap(),
    )
    .read_to_string(&mut legend)
    .unwrap();
    assert_eq!(legend, "id position a0 a1\nchr1:2_C_T 2 C T\n");
}

#[test]
fn convert_haplegendsample_supports_explicit_outputs_and_vcf_ids() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let hap = dir.path().join("out.hap");
    let legend = dir.path().join("out.legend");
    let samples = dir.path().join("out.samples");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0/1\t1|1\t0|0\n",
    )
    .unwrap();

    let outputs = format!(
        "{},{},{}",
        hap.display(),
        legend.display(),
        samples.display()
    );
    let (_out, err, code) = run(&[
        "convert",
        "--haplegendsample",
        &outputs,
        "--vcf-ids",
        "-s",
        "^B",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --haplegendsample explicit failed: {err}");
    assert_eq!(
        std::fs::read_to_string(&samples).unwrap(),
        "sample population group sex\nA A A 2\nC C C 2\n"
    );
    assert_eq!(std::fs::read_to_string(&hap).unwrap(), "0* 1* 0 0\n");
    assert_eq!(
        std::fs::read_to_string(&legend).unwrap(),
        "id position a0 a1\nrs1 2 C T\n"
    );
}

#[test]
fn convert_haplegendsample_matches_upstream_stdout_fixtures() {
    let vcf = fixture_path("convert.vcf");
    let hap_missing_vcf = fixture_path("convert.hap-missing.vcf");
    let expected_haps = std::fs::read_to_string(fixture_path("convert.hls.haps")).unwrap();
    let expected_legend = std::fs::read_to_string(fixture_path("convert.hls.legend")).unwrap();
    let expected_ids_legend =
        std::fs::read_to_string(fixture_path("convert.hls.ids.legend")).unwrap();
    let expected_samples = std::fs::read_to_string(fixture_path("convert.hls.samples")).unwrap();
    let expected_hap_missing =
        std::fs::read_to_string(fixture_path("convert.hap-missing.haps")).unwrap();

    let (haps_out, haps_err, haps_code) = run(&["convert", "-h", "-,.,.", vcf.to_str().unwrap()]);
    assert_eq!(haps_code, 0, "fixture -h -,.,. failed: {haps_err}");
    assert_eq!(haps_out, expected_haps);

    let (legend_out, legend_err, legend_code) =
        run(&["convert", "-h", ".,-,.", vcf.to_str().unwrap()]);
    assert_eq!(legend_code, 0, "fixture -h .,-,. failed: {legend_err}");
    assert_eq!(legend_out, expected_legend);

    let (ids_legend_out, ids_legend_err, ids_legend_code) =
        run(&["convert", "-h", ".,-,.", "--vcf-ids", vcf.to_str().unwrap()]);
    assert_eq!(
        ids_legend_code, 0,
        "fixture -h .,-,. --vcf-ids failed: {ids_legend_err}"
    );
    assert_eq!(ids_legend_out, expected_ids_legend);

    let (samples_out, samples_err, samples_code) =
        run(&["convert", "-h", ".,.,-", vcf.to_str().unwrap()]);
    assert_eq!(samples_code, 0, "fixture -h .,.,- failed: {samples_err}");
    assert_eq!(samples_out, expected_samples);

    let (hap_missing_out, hap_missing_err, hap_missing_code) = run(&[
        "convert",
        "--haplegendsample",
        "-,.,.",
        hap_missing_vcf.to_str().unwrap(),
    ]);
    assert_eq!(
        hap_missing_code, 0,
        "fixture --haplegendsample -,.,. missing GT failed: {hap_missing_err}"
    );
    assert_eq!(hap_missing_out, expected_hap_missing);
}

#[test]
fn convert_haplegendsample2vcf_writes_vcf_from_hap_legend_and_samples() {
    let dir = TempDir::new().unwrap();
    let hap = dir.path().join("in.hap");
    let legend = dir.path().join("in.legend");
    let samples = dir.path().join("in.samples");
    std::fs::write(&hap, "0 1 1* 0*\n? ? 1 -\n").unwrap();
    std::fs::write(
        &legend,
        "id position a0 a1\nchr1:2_C_T 2 C T\nchr1:5_A_G_9 5 A G\n",
    )
    .unwrap();
    std::fs::write(&samples, "sample population group sex\nA A A 2\nB B B 1\n").unwrap();
    let input = format!(
        "{},{},{}",
        hap.display(),
        legend.display(),
        samples.display()
    );
    let (out, err, code) = run(&["convert", "-H", "--no-version", &input]);
    assert_eq!(code, 0, "convert --haplegendsample2vcf failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"));
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t2\t.\tC\tT\t.\t.\t.\tGT\t0|1\t1/0",
            "chr1\t5\t.\tA\tG\t.\t.\tEND=9\tGT\t.|.\t1"
        ]
    );
    assert!(err.contains("Number of processed rows: \t2"), "got: {err}");
}

#[test]
fn convert_haplegendsample2vcf_supports_bcf_and_rejects_vcf_ids() {
    let dir = TempDir::new().unwrap();
    let hap = dir.path().join("in.hap");
    let legend = dir.path().join("in.legend");
    let samples = dir.path().join("in.samples");
    let bcf = dir.path().join("out.bcf");
    std::fs::write(&hap, "0 1\n").unwrap();
    std::fs::write(&legend, "id position a0 a1\nchr1:2_C_T 2 C T\n").unwrap();
    std::fs::write(&samples, "sample population group sex\nA A A 2\n").unwrap();
    let input = format!(
        "{},{},{}",
        hap.display(),
        legend.display(),
        samples.display()
    );
    let (_out, err, code) = run(&[
        "convert",
        "--haplegendsample2vcf",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        "-W",
        &input,
    ]);
    assert_eq!(code, 0, "convert --haplegendsample2vcf -Ob failed: {err}");
    assert!(dir.path().join("out.bcf.csi").exists());
    let (view_out, view_err, view_code) = run(&["view", "--no-version", bcf.to_str().unwrap()]);
    assert_eq!(view_code, 0, "view of converted BCF failed: {view_err}");
    assert!(view_out.contains("chr1\t2\t.\tC\tT\t.\t.\t.\tGT\t0|1"));

    let (_out, err, code) = run(&["convert", "-H", "--vcf-ids", &input]);
    assert_ne!(code, 0);
    assert!(
        err.contains("The option --haplegendsample2vcf cannot be combined with --vcf-ids"),
        "got: {err}"
    );
}

#[test]
fn convert_haplegendsample2vcf_matches_upstream_hap_legend_sample_fixture() {
    let hap = fixture_path("convert.hls.gt.hap");
    let legend = fixture_path("convert.hls.gt.legend");
    let samples = fixture_path("convert.hls.gt.samples");
    let expected = std::fs::read_to_string(fixture_path("convert.gt.noHead.vcf")).unwrap();

    let input = format!(
        "{},{},{}",
        hap.display(),
        legend.display(),
        samples.display()
    );
    let (out, err, code) = run(&["convert", "-H", "--no-version", &input]);
    assert_eq!(code, 0, "fixture -H failed: {err}");
    assert_eq!(without_meta_headers(&out), expected);
}

#[test]
fn convert_gensample_writes_gen_gz_and_samples_from_gt() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let sex = dir.path().join("sex.txt");
    let prefix = dir.path().join("impute");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
chr1\t2\trs1\tC\tT\t.\tPASS\tDP=8\tGT\t0/0\t0/1\t1/1\n\
chr1\t3\trs2\tG\tA,C\t.\tPASS\tDP=8\tGT\t0/1\t0/2\t./.\n",
    )
    .unwrap();
    std::fs::write(&sex, "A M\nC F\n").unwrap();

    let (_out, err, code) = run(&[
        "convert",
        "--gensample",
        prefix.to_str().unwrap(),
        "-s",
        "A,C",
        "--sex",
        sex.to_str().unwrap(),
        "-i",
        "DP>=8",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --gensample failed: {err}");
    assert!(err.contains(
        "1 records written, 1 skipped: 0/1/0/0 no-ALT/non-biallelic/filtered/duplicated"
    ));
    assert_eq!(
        std::fs::read_to_string(prefix.with_extension("samples")).unwrap(),
        "ID_1 ID_2 missing sex\n0 0 0 0\nA A 0 1\nC C 0 2\n"
    );
    let mut gen_text = String::new();
    flate2::read::MultiGzDecoder::new(
        std::fs::File::open(dir.path().join("impute.gen.gz")).unwrap(),
    )
    .read_to_string(&mut gen_text)
    .unwrap();
    assert_eq!(gen_text, "chr1:2_C_T chr1:2_C_T 2 C T 1 0 0 0 0 1\n");
}

#[test]
fn convert_gensample_supports_3n6_vcf_ids_duplicates_and_explicit_outputs() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let gen_path = dir.path().join("out.gen");
    let samples = dir.path().join("out.samples");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0|1\t1\n\
chr1\t2\trsdup\tC\tG\t.\tPASS\t.\tGT\t.\t0\n",
    )
    .unwrap();

    let outputs = format!("{},{}", gen_path.display(), samples.display());
    let (_out, err, code) = run(&[
        "convert",
        "--gensample",
        &outputs,
        "--3N6",
        "--vcf-ids",
        "--keep-duplicates",
        "-s",
        "^B",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --gensample explicit failed: {err}");
    assert_eq!(
        std::fs::read_to_string(&samples).unwrap(),
        "ID_1 ID_2 missing\n0 0 0\nA A 0\n"
    );
    assert_eq!(
        std::fs::read_to_string(&gen_path).unwrap(),
        "chr1 chr1:2_C_T rs1 2 C T 0 1 0\n\
chr1 chr1:2_C_G rsdup 2 C G 0.5 0.0 0.5\n"
    );
}

#[test]
fn convert_gensample_matches_upstream_stdout_fixtures() {
    let vcf = fixture_path("convert.vcf");
    let expected_gen = std::fs::read_to_string(fixture_path("convert.gs.gt.gen")).unwrap();
    let expected_ids = std::fs::read_to_string(fixture_path("convert.gs.gt.ids.gen")).unwrap();
    let expected_ids_3n6 = std::fs::read_to_string(fixture_path("convert.gs.gt.ids.gen6")).unwrap();
    let expected_samples = std::fs::read_to_string(fixture_path("convert.gs.gt.samples")).unwrap();
    let expected_pl_gen = std::fs::read_to_string(fixture_path("convert.gs.pl.gen")).unwrap();
    let expected_pl_samples =
        std::fs::read_to_string(fixture_path("convert.gs.pl.samples")).unwrap();

    let (gen_out, gen_err, gen_code) = run(&["convert", "-g", "-,.", vcf.to_str().unwrap()]);
    assert_eq!(gen_code, 0, "fixture -g -,. failed: {gen_err}");
    assert_eq!(gen_out, expected_gen);

    let (ids_out, ids_err, ids_code) =
        run(&["convert", "-g", "-,.", "--vcf-ids", vcf.to_str().unwrap()]);
    assert_eq!(ids_code, 0, "fixture -g -,. --vcf-ids failed: {ids_err}");
    assert_eq!(ids_out, expected_ids);

    let (ids_3n6_out, ids_3n6_err, ids_3n6_code) = run(&[
        "convert",
        "-g",
        "-,.",
        "--vcf-ids",
        "--3N6",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(
        ids_3n6_code, 0,
        "fixture -g -,. --vcf-ids --3N6 failed: {ids_3n6_err}"
    );
    assert_eq!(ids_3n6_out, expected_ids_3n6);

    let (samples_out, samples_err, samples_code) =
        run(&["convert", "-g", ".,-", vcf.to_str().unwrap()]);
    assert_eq!(samples_code, 0, "fixture -g .,- failed: {samples_err}");
    assert_eq!(samples_out, expected_samples);

    let (pl_gen_out, pl_gen_err, pl_gen_code) =
        run(&["convert", "-g", "-,.", "--tag", "PL", vcf.to_str().unwrap()]);
    assert_eq!(
        pl_gen_code, 0,
        "fixture -g -,. --tag PL failed: {pl_gen_err}"
    );
    assert_eq!(pl_gen_out, expected_pl_gen);

    let (pl_samples_out, pl_samples_err, pl_samples_code) =
        run(&["convert", "-g", ".,-", "--tag", "PL", vcf.to_str().unwrap()]);
    assert_eq!(
        pl_samples_code, 0,
        "fixture -g .,- --tag PL failed: {pl_samples_err}"
    );
    assert_eq!(pl_samples_out, expected_pl_samples);
}

#[test]
fn convert_gensample_matches_upstream_check_fixtures() {
    let vcf = fixture_path("check.vcf");
    let vcf = vcf.to_str().unwrap();
    let cases = [
        (
            vec!["convert", "-g", "-,.", "--vcf-ids", vcf],
            "check.gs.vcfids.gen",
        ),
        (
            vec!["convert", "-g", ".,-", "--vcf-ids", vcf],
            "check.gs.vcfids.samples",
        ),
        (
            vec!["convert", "-g", "-,.", "--3N6", vcf],
            "check.gs.chrom.gen",
        ),
        (
            vec!["convert", "-g", ".,-", "--3N6", vcf],
            "check.gs.chrom.samples",
        ),
        (
            vec!["convert", "-g", "-,.", "--3N6", "--vcf-ids", vcf],
            "check.gs.vcfids_chrom.gen",
        ),
        (
            vec!["convert", "-g", ".,-", "--3N6", "--vcf-ids", vcf],
            "check.gs.vcfids_chrom.samples",
        ),
    ];

    for (args, fixture) in cases {
        let expected = std::fs::read_to_string(fixture_path(fixture)).unwrap();
        let (out, err, code) = run(&args);
        assert_eq!(code, 0, "fixture {fixture} failed: {err}");
        assert_eq!(out, expected, "fixture {fixture} mismatch");
    }
}

#[test]
fn convert_gensample2vcf_writes_gt_gp_from_gen_and_samples() {
    let dir = TempDir::new().unwrap();
    let gen_path = dir.path().join("in.gen");
    let samples = dir.path().join("in.samples");
    std::fs::write(
        &gen_path,
        "chr1:2_C_T rs1 2 C T 0.9 0.1 0 0.1 0.8 0.1\n\
chr1:5_A_G_8 rs2 5 A G 0 0.2 0.8 0.34 0.33 0.33\n",
    )
    .unwrap();
    std::fs::write(&samples, "ID_1 ID_2 missing\n0 0 0\nA A 0\nB B 0\n").unwrap();
    let input = format!("{},{}", gen_path.display(), samples.display());
    let (out, err, code) = run(&["convert", "-G", "--vcf-ids", "--no-version", &input]);
    assert_eq!(code, 0, "convert --gensample2vcf failed: {err}");
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB"));
    let records: Vec<&str> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(
        records,
        [
            "chr1\t2\trs1\tC\tT\t.\t.\t.\tGT:GP\t0/0:0.9,0.1,0\t0/1:0.1,0.8,0.1",
            "chr1\t5\trs2\tA\tG\t.\t.\tEND=8\tGT:GP\t1/1:0,0.2,0.8\t0/0:0.34,0.33,0.33"
        ]
    );
    assert!(err.contains("Number of processed rows: \t2"), "got: {err}");
}

#[test]
fn convert_gensample2vcf_supports_3n6_bcf_and_index() {
    let dir = TempDir::new().unwrap();
    let gen_path = dir.path().join("in.gen");
    let samples = dir.path().join("in.samples");
    let bcf = dir.path().join("out.bcf");
    std::fs::write(&gen_path, "chr1 chr1:2_C_T rs1 2 C T 0 1 0\n").unwrap();
    std::fs::write(&samples, "ID_1 ID_2 missing\n0 0 0\nA A 0\n").unwrap();
    let input = format!("{},{}", gen_path.display(), samples.display());
    let (_out, err, code) = run(&[
        "convert",
        "--gensample2vcf",
        "--3N6",
        "--vcf-ids",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        "-W",
        &input,
    ]);
    assert_eq!(code, 0, "convert --gensample2vcf -Ob failed: {err}");
    assert!(dir.path().join("out.bcf.csi").exists());
    let (view_out, view_err, view_code) = run(&["view", "--no-version", bcf.to_str().unwrap()]);
    assert_eq!(view_code, 0, "view of converted BCF failed: {view_err}");
    assert!(view_out.contains("chr1\t2\trs1\tC\tT\t.\t.\t.\tGT:GP\t0/1:0,1,0"));
}

#[test]
fn convert_gensample2vcf_matches_upstream_gen_sample_fixtures() {
    let samples = fixture_path("convert.gs.gt.samples");
    let gen_path = fixture_path("convert.gs.gt.ids.gen");
    let gen_3n6 = fixture_path("convert.gs.gt.ids.3N6.gen");
    let gen_rev = fixture_path("convert.gs.gt.ids.gen.rev");
    let expected_ids = std::fs::read_to_string(fixture_path("convert.gs.vcf")).unwrap();
    let expected_noids = std::fs::read_to_string(fixture_path("convert.gs.noids.vcf")).unwrap();

    let input = format!("{},{}", gen_path.display(), samples.display());
    let (ids_out, ids_err, ids_code) = run(&["convert", "--vcf-ids", "-G", "--no-version", &input]);
    assert_eq!(ids_code, 0, "fixture --vcf-ids -G failed: {ids_err}");
    assert_eq!(without_meta_headers(&ids_out), expected_ids);

    let (noids_out, noids_err, noids_code) = run(&["convert", "-G", "--no-version", &input]);
    assert_eq!(noids_code, 0, "fixture -G failed: {noids_err}");
    assert_eq!(without_meta_headers(&noids_out), expected_noids);

    let input_3n6 = format!("{},{}", gen_3n6.display(), samples.display());
    let (out_3n6, err_3n6, code_3n6) = run(&["convert", "--3N6", "-G", "--no-version", &input_3n6]);
    assert_eq!(code_3n6, 0, "fixture --3N6 -G failed: {err_3n6}");
    assert_eq!(without_meta_headers(&out_3n6), expected_noids);

    let input_rev = format!("{},{}", gen_rev.display(), samples.display());
    let (rev_ids_out, rev_ids_err, rev_ids_code) =
        run(&["convert", "--vcf-ids", "-G", "--no-version", &input_rev]);
    assert_eq!(
        rev_ids_code, 0,
        "fixture reversed --vcf-ids -G failed: {rev_ids_err}"
    );
    assert_eq!(without_meta_headers(&rev_ids_out), expected_ids);

    let (rev_out, rev_err, rev_code) = run(&["convert", "-G", "--no-version", &input_rev]);
    assert_eq!(rev_code, 0, "fixture reversed -G failed: {rev_err}");
    assert_eq!(without_meta_headers(&rev_out), expected_noids);
}

#[test]
fn convert_gensample_skips_duplicate_positions_by_default() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let gen_path = dir.path().join("out.gen");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGT\t0/1\n\
chr1\t2\trsdup\tC\tG\t.\tPASS\t.\tGT\t1/1\n",
    )
    .unwrap();
    let outputs = format!("{},.", gen_path.display());
    let (_out, err, code) = run(&["convert", "--gensample", &outputs, vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "convert --gensample duplicate failed: {err}");
    assert!(err.contains(
        "1 records written, 1 skipped: 0/0/0/1 no-ALT/non-biallelic/filtered/duplicated"
    ));
    assert_eq!(
        std::fs::read_to_string(&gen_path).unwrap(),
        "chr1:2_C_T chr1:2_C_T 2 C T 0 1 0\n"
    );
}

#[test]
fn convert_gensample_supports_gp_tag_probabilities() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let gen_path = dir.path().join("out.gen");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype posterior probabilities\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGP\t0.8,0.1,0.1\t0.2,0.3,0.5\n",
    )
    .unwrap();
    let outputs = format!("{},.", gen_path.display());
    let (_out, err, code) = run(&[
        "convert",
        "--gensample",
        &outputs,
        "--tag",
        "GP",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --gensample --tag GP failed: {err}");
    assert_eq!(
        std::fs::read_to_string(&gen_path).unwrap(),
        "chr1:2_C_T chr1:2_C_T 2 C T 0.800000 0.100000 0.100000 0.200000 0.300000 0.500000\n"
    );
}

#[test]
fn convert_gensample_supports_pl_tag_likelihoods() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let gen_path = dir.path().join("out.gen");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tPL\t0,20,40\t40,20,0\n",
    )
    .unwrap();
    let outputs = format!("{},.", gen_path.display());
    let (_out, err, code) = run(&[
        "convert",
        "--gensample",
        &outputs,
        "--tag=PL",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --gensample --tag PL failed: {err}");
    assert_eq!(
        std::fs::read_to_string(&gen_path).unwrap(),
        "chr1:2_C_T chr1:2_C_T 2 C T 0.990001 0.009900 0.000099 0.000099 0.009900 0.990001\n"
    );
}

#[test]
fn convert_gensample_supports_gl_tag_likelihoods() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    let gen_path = dir.path().join("out.gen");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods\">\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\tGL\t0,-2,-4\t-4,-2,0\n",
    )
    .unwrap();
    let outputs = format!("{},.", gen_path.display());
    let (_out, err, code) = run(&[
        "convert",
        "--gensample",
        &outputs,
        "--tag",
        "GL",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "convert --gensample --tag GL failed: {err}");
    assert_eq!(
        std::fs::read_to_string(&gen_path).unwrap(),
        "chr1:2_C_T chr1:2_C_T 2 C T 0.990001 0.009900 0.000099 0.000099 0.009900 0.990001\n"
    );
}

#[test]
fn convert_gensample_rejects_deprecated_chrom_option() {
    let dir = TempDir::new().unwrap();
    let vcf = dir.path().join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=10>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t2\trs1\tC\tT\t.\tPASS\t.\n",
    )
    .unwrap();
    let outputs = format!("{},.", dir.path().join("out.gen").display());
    let (_out, err, code) = run(&[
        "convert",
        "--gensample",
        &outputs,
        "--chrom",
        vcf.to_str().unwrap(),
    ]);
    assert_ne!(code, 0);
    assert!(
        err.contains("The --chrom option has been deprecated, please use --3N6 instead"),
        "got: {err}"
    );
}

#[test]
fn convert_tsv2vcf_rejects_include_exclude_filters() {
    let dir = TempDir::new().unwrap();
    let tsv = dir.path().join("in.tsv");
    std::fs::write(&tsv, "chr1\t1\tA\tC\n").unwrap();
    let (_out, err, code) = run(&[
        "convert",
        "--tsv2vcf",
        tsv.to_str().unwrap(),
        "-c",
        "CHROM,POS,REF,ALT",
        "-i",
        "POS=1",
    ]);
    assert_ne!(code, 0);
    assert!(
        err.contains("-i/-e are only supported for VCF/gVCF input"),
        "got: {err}"
    );
}

#[test]
fn convert_no_args_prints_usage() {
    let (_out, err, code) = run(&["convert"]);
    assert_ne!(code, 0);
    assert!(err.contains("Usage:"), "got: {err}");
}
