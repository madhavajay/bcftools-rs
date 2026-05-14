//! End-to-end tests for `bcftools_rs::commands::stats`.

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

const VCF: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\trs1\tA\tG\t100\tPASS\tAF=0.05\n\
1\t200\t.\tT\tA\t100\tPASS\tAF=0.5\n\
1\t300\t.\tA\tAT\t100\tPASS\tAF=0.2\n\
1\t400\t.\tC\tT,G\t100\tPASS\tAF=0.7\n\
1\t500\t.\tA\t.\t100\tPASS\t.\n";

fn write_vcf(dir: &TempDir, body: &str) -> PathBuf {
    let p = dir.path().join("in.vcf");
    std::fs::write(&p, body).unwrap();
    p
}

fn extract_value<'a>(out: &'a str, prefix: &str) -> Option<&'a str> {
    out.lines()
        .find(|l| l.starts_with(prefix))
        .map(|l| l.rsplit_once('\t').map(|(_, v)| v).unwrap_or(""))
}

#[test]
fn stats_emits_expected_summary_numbers() {
    let dir = TempDir::new().unwrap();
    let v = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["stats", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats failed: {err}");

    assert!(out.contains("# SN, Summary numbers:"));
    assert_eq!(extract_value(&out, "SN\t0\tnumber of records:"), Some("5"));
    assert_eq!(extract_value(&out, "SN\t0\tnumber of no-ALTs:"), Some("1"));
    // 1@100 (A>G), 1@200 (T>A), 1@400 (C>T,C>G - both SNPs)
    assert_eq!(extract_value(&out, "SN\t0\tnumber of SNPs:"), Some("3"));
    assert_eq!(extract_value(&out, "SN\t0\tnumber of indels:"), Some("1"));
    assert_eq!(
        extract_value(&out, "SN\t0\tnumber of multiallelic sites:"),
        Some("1")
    );
    assert_eq!(
        extract_value(&out, "SN\t0\tnumber of multiallelic SNP sites:"),
        Some("1")
    );
}

#[test]
fn stats_tstv_section_lists_ts_and_tv() {
    let dir = TempDir::new().unwrap();
    let v = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["stats", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats failed: {err}");
    assert!(out.contains("# TSTV"));
    let tstv_line = out.lines().find(|l| l.starts_with("TSTV\t0\t")).unwrap();
    let cols: Vec<&str> = tstv_line.split('\t').collect();
    // ts: A>G(100), C>T(400) = 2 transitions; tv: T>A(200), C>G(400) = 2 transversions.
    assert_eq!(cols[2], "2");
    assert_eq!(cols[3], "2");
}

#[test]
fn stats_st_section_lists_zero_count_substitution_rows() {
    let dir = TempDir::new().unwrap();
    let v = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["stats", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats failed: {err}");
    let substitution_rows: Vec<&str> = out
        .lines()
        .filter(|line| line.starts_with("ST\t0\t"))
        .collect();
    assert_eq!(substitution_rows.len(), 12, "unexpected ST rows:\n{out}");
    assert!(
        substitution_rows.contains(&"ST\t0\tA>C\t0"),
        "missing zero-count A>C row:\n{out}"
    );
    assert!(
        substitution_rows.contains(&"ST\t0\tA>G\t1"),
        "missing observed A>G row:\n{out}"
    );
}

#[test]
fn stats_apply_filters_drops_records() {
    let dir = TempDir::new().unwrap();
    let body = VCF.replace(
        "1\t100\trs1\tA\tG\t100\tPASS",
        "1\t100\trs1\tA\tG\t100\tQ10",
    );
    let v = write_vcf(&dir, &body);
    let (out, err, code) = run(&["stats", "-f", "PASS", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -f failed: {err}");
    // The PASS filter should drop the Q10 record at pos 100.
    assert_eq!(extract_value(&out, "SN\t0\tnumber of records:"), Some("4"));
}

#[test]
fn stats_split_by_id_separates_known_from_novel() {
    let dir = TempDir::new().unwrap();
    let v = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["stats", "-I", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -I failed: {err}");
    // id=1 is "known": only the rs1 record at pos 100.
    assert_eq!(extract_value(&out, "SN\t1\tnumber of records:"), Some("1"));
    // id=2 is "novel": the remaining 4 records.
    assert_eq!(extract_value(&out, "SN\t2\tnumber of records:"), Some("4"));
}

#[test]
fn stats_first_allele_only_ignores_extra_alts() {
    let dir = TempDir::new().unwrap();
    let v = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["stats", "-1", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -1 failed: {err}");
    // multiallelic SNP site (C>T,G) is still counted once for the SNP type but ts/tv
    // should now include only the first allele (C>T = transition).
    let tstv_line = out.lines().find(|l| l.starts_with("TSTV\t0\t")).unwrap();
    let cols: Vec<&str> = tstv_line.split('\t').collect();
    // ts: A>G(100), C>T(400) = 2 transitions; tv: T>A(200) = 1 transversion (C>G dropped).
    assert_eq!(cols[2], "2");
    assert_eq!(cols[3], "1");
}

#[test]
fn stats_include_expression_filters_records() {
    let dir = TempDir::new().unwrap();
    let v = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["stats", "-i", "POS<400", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -i failed: {err}");
    assert_eq!(extract_value(&out, "SN\t0\tnumber of records:"), Some("3"));
}

#[test]
fn stats_af_tag_uses_named_info_tag_for_af_bins() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Default allele frequency\">\n\
##INFO=<ID=XAF,Number=A,Type=Float,Description=\"Alternate allele frequency\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tG\t100\tPASS\tAF=0.9;XAF=0.05\n\
1\t200\t.\tC\tT\t100\tPASS\tAF=0.05;XAF=0.9\n";
    let v = write_vcf(&dir, body);
    let (out, err, code) = run(&[
        "stats",
        "--af-tag",
        "XAF",
        "--af-bins",
        "0.1,0.5",
        v.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "stats --af-tag failed: {err}");
    let af_01 = out
        .lines()
        .find(|l| l.starts_with("AF\t0\t0.100000\t"))
        .expect("0.1 AF bin");
    let af_10 = out
        .lines()
        .find(|l| l.starts_with("AF\t0\t1.000000\t"))
        .expect("overflow AF bin");
    assert_eq!(af_01.split('\t').nth(3), Some("1"));
    assert_eq!(af_10.split('\t').nth(3), Some("1"));
}

#[test]
fn stats_depth_distribution_counts_info_dp_and_sample_depths() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tG\t100\tPASS\tDP=7\tGT:DP\t0/1:5\t0/0:12\n\
1\t200\t.\tC\tT\t100\tPASS\tDP=22\tGT:AD\t0/1:3,4\t0/1:1,1\n";
    let v = write_vcf(&dir, body);
    let (out, err, code) = run(&["stats", "--depth", "0,20,10", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats --depth failed: {err}");
    assert!(
        out.contains("# DP, Depth distribution"),
        "missing DP section:\n{out}"
    );
    let bin0 = out
        .lines()
        .find(|l| l.starts_with("DP\t0\t0\t"))
        .expect("0-depth bin");
    let bin1 = out
        .lines()
        .find(|l| l.starts_with("DP\t0\t1\t"))
        .expect("10-depth bin");
    let overflow = out
        .lines()
        .find(|l| l.starts_with("DP\t0\t>20\t"))
        .expect("overflow-depth bin");
    let cols0: Vec<&str> = bin0.split('\t').collect();
    let cols1: Vec<&str> = bin1.split('\t').collect();
    let cols_overflow: Vec<&str> = overflow.split('\t').collect();
    assert_eq!(cols0[3], "3", "genotypes in <=10 bucket:\n{out}");
    assert_eq!(cols0[5], "1", "sites in <=10 bucket:\n{out}");
    assert_eq!(cols1[3], "1", "genotypes in 11-20 bucket:\n{out}");
    assert_eq!(cols_overflow[5], "1", "sites in >20 bucket:\n{out}");
}

#[test]
fn stats_samples_selection_emits_psc_for_selected_sample() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tG\t100\tPASS\tDP=10\tGT:DP\t0/0:4\t0/1:8\n\
1\t200\t.\tC\tT\t100\tPASS\tDP=20\tGT:DP\t1/1:6\t./.:.\n\
1\t300\t.\tA\tAT\t100\tPASS\tDP=30\tGT:DP\t0/1:5\t0/1:7\n\
1\t400\t.\tG\tA\t100\tPASS\tDP=40\tGT:DP\t0/1:6\t0/0:6\n";
    let v = write_vcf(&dir, body);
    let (out, err, code) = run(&["stats", "-s", "S1", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -s failed: {err}");
    assert_eq!(extract_value(&out, "SN\t0\tnumber of samples:"), Some("1"));
    assert!(out.contains("# PSC"), "missing PSC section:\n{out}");
    let psc = out
        .lines()
        .find(|l| l.starts_with("PSC\t0\tS1\t"))
        .expect("S1 PSC row");
    let cols: Vec<&str> = psc.split('\t').collect();
    assert_eq!(cols[3], "1", "ref hom count:\n{out}");
    assert_eq!(cols[4], "1", "non-ref hom count:\n{out}");
    assert_eq!(cols[5], "1", "SNP het count:\n{out}");
    assert_eq!(cols[6], "2", "transition count:\n{out}");
    assert_eq!(cols[8], "1", "indel count:\n{out}");
    assert_eq!(cols[9], "5.2", "average depth:\n{out}");
    assert_eq!(cols[10], "2", "singleton count:\n{out}");
}

#[test]
fn stats_samples_selection_emits_psi_indel_counts() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tAT\t100\tPASS\t.\tGT\t0/1\t1/1\n\
1\t200\t.\tAT\tA\t100\tPASS\t.\tGT\t0/1\t1/1\n\
1\t300\t.\tA\tAT,ATT\t100\tPASS\t.\tGT\t1/2\t0/0\n\
1\t400\t.\tAT\tA,ATT\t100\tPASS\t.\tGT\t1/2\t0/0\n";
    let v = write_vcf(&dir, body);
    let (out, err, code) = run(&["stats", "-s", "S1,S2", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -s failed: {err}");
    assert!(out.contains("# PSI"), "missing PSI section:\n{out}");

    let s1 = out
        .lines()
        .find(|line| line.starts_with("PSI\t0\tS1\t"))
        .expect("S1 PSI row");
    let s1_cols: Vec<&str> = s1.split('\t').collect();
    assert_eq!(s1_cols[7], "3", "S1 insertion hets:\n{out}");
    assert_eq!(s1_cols[8], "2", "S1 deletion hets:\n{out}");
    assert_eq!(s1_cols[9], "0", "S1 insertion homs:\n{out}");
    assert_eq!(s1_cols[10], "0", "S1 deletion homs:\n{out}");

    let s2 = out
        .lines()
        .find(|line| line.starts_with("PSI\t0\tS2\t"))
        .expect("S2 PSI row");
    let s2_cols: Vec<&str> = s2.split('\t').collect();
    assert_eq!(s2_cols[7], "0", "S2 insertion hets:\n{out}");
    assert_eq!(s2_cols[8], "0", "S2 deletion hets:\n{out}");
    assert_eq!(s2_cols[9], "1", "S2 insertion homs:\n{out}");
    assert_eq!(s2_cols[10], "1", "S2 deletion homs:\n{out}");
}

#[test]
fn stats_exons_emit_frameshift_counts_and_sample_psi_frames() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
1\t100\t.\tA\tAT\t100\tPASS\t.\tGT\t0/1\n\
1\t200\t.\tA\tATTT\t100\tPASS\t.\tGT\t1/1\n\
1\t300\t.\tAT\tA\t100\tPASS\t.\tGT\t0/1\n";
    let v = write_vcf(&dir, body);
    let exons = dir.path().join("exons.txt");
    std::fs::write(&exons, "1\t100\t250\n").unwrap();

    let (out, err, code) = run(&[
        "stats",
        "-E",
        exons.to_str().unwrap(),
        "-s",
        "S1",
        v.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "stats -E failed: {err}");

    let fs = out
        .lines()
        .find(|line| line.starts_with("FS\t0\t"))
        .expect("FS row");
    let fs_cols: Vec<&str> = fs.split('\t').collect();
    assert_eq!(fs_cols[2], "1", "in-frame site count:\n{out}");
    assert_eq!(fs_cols[3], "1", "out-frame site count:\n{out}");
    assert_eq!(fs_cols[4], "1", "not-applicable site count:\n{out}");

    let psi = out
        .lines()
        .find(|line| line.starts_with("PSI\t0\tS1\t"))
        .expect("S1 PSI row");
    let psi_cols: Vec<&str> = psi.split('\t').collect();
    assert_eq!(psi_cols[3], "2", "sample in-frame alleles:\n{out}");
    assert_eq!(psi_cols[4], "1", "sample out-frame alleles:\n{out}");
    assert_eq!(psi_cols[5], "1", "sample not-applicable alleles:\n{out}");
}

#[test]
fn stats_fasta_ref_emits_indel_context_sections() {
    let dir = TempDir::new().unwrap();
    let fasta = dir.path().join("ref.fa");
    std::fs::write(
        &fasta,
        ">1\nAATATATATATATATATATATATATATATATATATATATATATATATATAT\n",
    )
    .unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=60>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tAAT\tA\t100\tPASS\t.\n\
1\t1\t.\tA\tATT\t100\tPASS\t.\n\
1\t1\t.\tA\tAT\t100\tPASS\t.\n";
    let v = write_vcf(&dir, body);

    let (out, err, code) = run(&["stats", "-F", fasta.to_str().unwrap(), v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -F failed: {err}");
    assert!(out.contains("# ICS"), "missing ICS section:\n{out}");
    assert!(out.contains("# ICL"), "missing ICL section:\n{out}");

    let ics = out
        .lines()
        .find(|line| line.starts_with("ICS\t0\t"))
        .expect("ICS row");
    let ics_cols: Vec<&str> = ics.split('\t').collect();
    assert_eq!(ics_cols[2], "2", "consistent indels:\n{out}");
    assert_eq!(ics_cols[3], "1", "inconsistent indels:\n{out}");
    assert_eq!(ics_cols[4], "0", "not-applicable indels:\n{out}");

    let icl = out
        .lines()
        .find(|line| line.starts_with("ICL\t0\t2\t"))
        .expect("ICL repeat length 2 row");
    let icl_cols: Vec<&str> = icl.split('\t').collect();
    assert_eq!(icl_cols[3], "1", "consistent deletions:\n{out}");
    assert_eq!(icl_cols[4], "0", "inconsistent deletions:\n{out}");
    assert_eq!(icl_cols[5], "1", "consistent insertions:\n{out}");
    assert_eq!(icl_cols[6], "1", "inconsistent insertions:\n{out}");
}

#[test]
fn stats_samples_all_emit_singleton_section() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/1\t0/0\n\
1\t200\t.\tC\tA\t100\tPASS\t.\tGT\t0/0\t0/1\n\
1\t300\t.\tA\tAT\t100\tPASS\t.\tGT\t0/1\t0/0\n\
1\t400\t.\tG\tT\t100\tPASS\t.\tGT\t0/1\t0/1\n";
    let v = write_vcf(&dir, body);

    let (out, err, code) = run(&["stats", "-s", "-", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -s - failed: {err}");
    assert!(out.contains("# SiS"), "missing SiS section:\n{out}");
    let sis = out
        .lines()
        .find(|line| line.starts_with("SiS\t0\t"))
        .expect("SiS row");
    let cols: Vec<&str> = sis.split('\t').collect();
    assert_eq!(cols[3], "2", "singleton SNP count:\n{out}");
    assert_eq!(cols[4], "1", "singleton transition count:\n{out}");
    assert_eq!(cols[5], "1", "singleton transversion count:\n{out}");
    assert_eq!(cols[6], "1", "singleton indel count:\n{out}");
    assert_eq!(cols[9], "1", "singleton indel context n/a:\n{out}");
}

#[test]
fn stats_samples_emit_vaf_distribution_from_format_ad() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tG\t100\tPASS\t.\tGT:AD\t0/1:8,2\t0/0:10,0\n\
1\t200\t.\tA\tAT\t100\tPASS\t.\tGT:AD\t0/1:4,6\t0/1:9,1\n";
    let v = write_vcf(&dir, body);

    let (out, err, code) = run(&["stats", "-s", "S1,S2", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -s failed: {err}");
    assert!(out.contains("# VAF"), "missing VAF section:\n{out}");

    let s1 = out
        .lines()
        .find(|line| line.starts_with("VAF\t0\tS1\t"))
        .expect("S1 VAF row");
    let s1_cols: Vec<&str> = s1.split('\t').collect();
    let s1_snv: Vec<&str> = s1_cols[3].split(',').collect();
    let s1_indel: Vec<&str> = s1_cols[4].split(',').collect();
    assert_eq!(s1_snv[4], "1", "S1 SNV VAF 0.20 bin:\n{out}");
    assert_eq!(s1_indel[12], "1", "S1 indel VAF 0.60 bin:\n{out}");

    let s2 = out
        .lines()
        .find(|line| line.starts_with("VAF\t0\tS2\t"))
        .expect("S2 VAF row");
    let s2_cols: Vec<&str> = s2.split('\t').collect();
    let s2_indel: Vec<&str> = s2_cols[4].split(',').collect();
    assert_eq!(s2_indel[2], "1", "S2 indel VAF 0.10 bin:\n{out}");
}

#[test]
fn stats_emits_quality_and_indel_distribution_sections() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tG\t12.5\tPASS\t.\tGT:AD\t0/1:8,2\t0/0:10,0\n\
1\t200\t.\tC\tA\t9.0\tPASS\t.\tGT:AD\t0/1:5,5\t0/0:10,0\n\
1\t300\t.\tATTT\tA\t20.0\tPASS\t.\tGT:AD\t0/1:7,3\t1/1:2,8\n\
1\t400\t.\tA\tAT,ATTT\t.\tPASS\t.\tGT:AD\t1/2:4,2,4\t0/0:10,0,0\n";
    let v = write_vcf(&dir, body);

    let (out, err, code) = run(&["stats", "-s", "S1,S2", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -s failed: {err}");
    assert!(
        out.contains("# QUAL, Stats by quality"),
        "missing QUAL:\n{out}"
    );
    assert!(
        out.contains("# IDD, InDel distribution:"),
        "missing IDD:\n{out}"
    );

    let qual_9 = out
        .lines()
        .find(|line| line.starts_with("QUAL\t0\t9.0\t"))
        .expect("quality 9.0 row");
    let qual_9_cols: Vec<&str> = qual_9.split('\t').collect();
    assert_eq!(qual_9_cols[3], "1", "SNPs at quality 9.0:\n{out}");
    assert_eq!(qual_9_cols[4], "0", "transitions at quality 9.0:\n{out}");
    assert_eq!(qual_9_cols[5], "1", "transversions at quality 9.0:\n{out}");

    let qual_20 = out
        .lines()
        .find(|line| line.starts_with("QUAL\t0\t20.0\t"))
        .expect("quality 20.0 row");
    assert_eq!(
        qual_20.split('\t').nth(6),
        Some("1"),
        "indels at quality 20.0:\n{out}"
    );

    let deletion = out
        .lines()
        .find(|line| line.starts_with("IDD\t0\t-3\t"))
        .expect("deletion length -3 row");
    let deletion_cols: Vec<&str> = deletion.split('\t').collect();
    assert_eq!(deletion_cols[3], "1", "deletion sites:\n{out}");
    assert_eq!(deletion_cols[4], "2", "deletion genotypes:\n{out}");
    assert_eq!(deletion_cols[5], "0.55", "deletion mean VAF:\n{out}");

    let insertion_one = out
        .lines()
        .find(|line| line.starts_with("IDD\t0\t1\t"))
        .expect("insertion length 1 row");
    let insertion_one_cols: Vec<&str> = insertion_one.split('\t').collect();
    assert_eq!(
        insertion_one_cols[3], "1",
        "insertion length 1 sites:\n{out}"
    );
    assert_eq!(
        insertion_one_cols[4], "1",
        "insertion length 1 genotypes:\n{out}"
    );
    assert_eq!(
        insertion_one_cols[5], "0.20",
        "insertion length 1 mean VAF:\n{out}"
    );
}

#[test]
fn stats_samples_emit_hwe_quartiles_from_gt_counts() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t0/1\t1/1\n\
1\t200\t.\tC\tT\t100\tPASS\t.\tGT\t0/0\t0/0\t0/1\t0/1\n\
1\t300\t.\tA\tAT\t100\tPASS\t.\tGT\t0/1\t0/1\t0/0\t0/0\n";
    let v = write_vcf(&dir, body);

    let (out, err, code) = run(&["stats", "-s", "-", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats -s - failed: {err}");
    assert!(out.contains("# HWE"), "missing HWE section:\n{out}");
    let hwe_49 = out
        .lines()
        .find(|line| line.starts_with("HWE\t0\t0.490000\t"))
        .expect("HWE 0.49 AF bin");
    let cols: Vec<&str> = hwe_49.split('\t').collect();
    assert_eq!(cols[3], "1", "one SNP in 0.49 AF bin:\n{out}");

    let hwe_25 = out
        .lines()
        .find(|line| line.starts_with("HWE\t0\t0.250000\t"))
        .expect("HWE 0.25 AF bin");
    let cols: Vec<&str> = hwe_25.split('\t').collect();
    assert_eq!(cols[3], "1", "one SNP in 0.25 AF bin:\n{out}");
    assert!(
        !out.lines()
            .any(|line| line.starts_with("HWE\t0\t0.120000\t"))
    );
}

#[test]
fn stats_samples_file_and_exclusion_select_expected_samples() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
1\t100\t.\tA\tG\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1\n";
    let v = write_vcf(&dir, body);
    let samples = dir.path().join("samples.txt");
    std::fs::write(&samples, "S2\n").unwrap();
    let (out, err, code) = run(&[
        "stats",
        "-S",
        samples.to_str().unwrap(),
        v.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "stats -S failed: {err}");
    assert!(out.lines().any(|l| l.starts_with("PSC\t0\tS2\t")));
    assert!(!out.lines().any(|l| l.starts_with("PSC\t0\tS1\t")));

    let (excluded_out, excluded_err, excluded_code) =
        run(&["stats", "-s", "^S2", v.to_str().unwrap()]);
    assert_eq!(excluded_code, 0, "stats -s ^ failed: {excluded_err}");
    assert!(excluded_out.lines().any(|l| l.starts_with("PSC\t0\tS1\t")));
    assert!(excluded_out.lines().any(|l| l.starts_with("PSC\t0\tS3\t")));
    assert!(!excluded_out.lines().any(|l| l.starts_with("PSC\t0\tS2\t")));
}

#[test]
fn stats_user_tstv_bins_first_alt_tstv_by_info_tag() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=QS,Number=1,Type=Float,Description=\"Quality score\">\n\
##INFO=<ID=PV4,Number=4,Type=Float,Description=\"Indexed values\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tG\t100\tPASS\tQS=0.0;PV4=9,0.0,9,9\n\
1\t200\t.\tC\tA\t100\tPASS\tQS=1.0;PV4=9,1.0,9,9\n\
1\t300\t.\tA\tAT\t100\tPASS\tQS=0.5;PV4=9,0.5,9,9\n";
    let v = write_vcf(&dir, body);
    let (out, err, code) = run(&["stats", "--user-tstv", "PV4[1]:0:1:3", v.to_str().unwrap()]);
    assert_eq!(code, 0, "stats --user-tstv failed: {err}");
    assert!(out.contains("# USR:PV4/1"), "missing USR header:\n{out}");
    let ts_bin = out
        .lines()
        .find(|l| l.starts_with("USR:PV4/1\t0\t0.000000\t"))
        .expect("transition user bin");
    let tv_bin = out
        .lines()
        .find(|l| l.starts_with("USR:PV4/1\t0\t1.000000\t"))
        .expect("transversion user bin");
    assert_eq!(ts_bin.split('\t').collect::<Vec<_>>()[3..], ["1", "1", "0"]);
    assert_eq!(tv_bin.split('\t').collect::<Vec<_>>()[3..], ["1", "0", "1"]);
    assert!(
        !out.lines()
            .any(|line| line.starts_with("USR:PV4/1\t0\t0.500000\t")),
        "non-SNP indel should not contribute to user-tstv:\n{out}"
    );
}

#[test]
fn stats_pairwise_reports_left_right_and_shared_sets() {
    let dir = TempDir::new().unwrap();
    let left = dir.path().join("left.vcf");
    let right = dir.path().join("right.vcf");
    let header = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    std::fs::write(
        &left,
        format!(
            "{header}\
1\t100\t.\tA\tG\t100\tPASS\t.\n\
1\t200\t.\tC\tT\t100\tPASS\t.\n"
        ),
    )
    .unwrap();
    std::fs::write(
        &right,
        format!(
            "{header}\
1\t100\t.\tA\tG\t100\tPASS\t.\n\
1\t300\t.\tG\tA\t100\tPASS\t.\n"
        ),
    )
    .unwrap();

    let (out, err, code) = run(&["stats", left.to_str().unwrap(), right.to_str().unwrap()]);
    assert_eq!(code, 0, "pairwise stats failed: {err}");
    assert_eq!(extract_value(&out, "SN\t0\tnumber of records:"), Some("1"));
    assert_eq!(extract_value(&out, "SN\t1\tnumber of records:"), Some("1"));
    assert_eq!(extract_value(&out, "SN\t2\tnumber of records:"), Some("1"));
    assert!(
        out.contains("ID\t2\t"),
        "missing shared set definition:\n{out}"
    );
}

#[test]
fn stats_pairwise_collapse_snps_matches_same_position_snps() {
    let dir = TempDir::new().unwrap();
    let left = dir.path().join("left.vcf");
    let right = dir.path().join("right.vcf");
    let header = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    std::fs::write(&left, format!("{header}1\t100\t.\tA\tG\t100\tPASS\t.\n")).unwrap();
    std::fs::write(&right, format!("{header}1\t100\t.\tA\tT\t100\tPASS\t.\n")).unwrap();

    let (without, without_err, without_code) =
        run(&["stats", left.to_str().unwrap(), right.to_str().unwrap()]);
    assert_eq!(without_code, 0, "pairwise stats failed: {without_err}");
    assert_eq!(
        extract_value(&without, "SN\t2\tnumber of records:"),
        Some("0"),
        "exact matching should not collapse different ALTs:\n{without}"
    );

    let (collapsed, err, code) = run(&[
        "stats",
        "--collapse",
        "snps",
        left.to_str().unwrap(),
        right.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "pairwise collapse stats failed: {err}");
    assert_eq!(
        extract_value(&collapsed, "SN\t2\tnumber of records:"),
        Some("1"),
        "SNP collapse should match same-position SNPs:\n{collapsed}"
    );
}

#[test]
fn stats_no_args_prints_usage() {
    let (_out, err, code) = run(&["stats"]);
    assert_ne!(code, 0);
    assert!(err.contains("Usage:"), "got: {err}");
}
