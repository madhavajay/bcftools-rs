//! End-to-end tests for `bcftools_rs::commands::filter`.

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

const VCF: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tC\t100\tPASS\tDP=10;AF=0.05\n\
1\t200\t.\tG\tT\t10\tPASS\tDP=2;AF=0.5\n\
1\t300\t.\tT\tA\t50\tPASS\tDP=15;AF=0.20\n";

fn write_vcf(dir: &TempDir, body: &str) -> PathBuf {
    let p = dir.path().join("in.vcf");
    std::fs::write(&p, body).unwrap();
    p
}

#[test]
fn filter_include_keeps_only_matching_records() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["filter", "-i", "DP>=10", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter -i failed: {err}");
    let recs: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(recs.len(), 2);
    assert!(recs[0].starts_with("1\t100\t"));
    assert!(recs[1].starts_with("1\t300\t"));
}

#[test]
fn filter_default_injects_bcftools_version_and_command_lines() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["filter", "-i", "DP>=10", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter failed: {err}");
    let version_line = out
        .lines()
        .find(|l| l.starts_with("##bcftools_filterVersion="))
        .unwrap_or_else(|| panic!("missing filter version line:\n{out}"));
    let command_line = out
        .lines()
        .find(|l| l.starts_with("##bcftools_filterCommand="))
        .unwrap_or_else(|| panic!("missing filter command line:\n{out}"));
    assert!(version_line.contains("+htslib-"), "got: {version_line}");
    assert!(command_line.contains("filter"), "got: {command_line}");
    assert!(command_line.contains("; Date="), "got: {command_line}");
}

#[test]
fn filter_no_version_suppresses_bcftools_header_lines() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, err, code) = run(&[
        "filter",
        "--no-version",
        "-i",
        "DP>=10",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --no-version failed: {err}");
    assert!(
        !out.contains("##bcftools_filterVersion="),
        "version line leaked despite --no-version:\n{out}"
    );
    assert!(
        !out.contains("##bcftools_filterCommand="),
        "command line leaked despite --no-version:\n{out}"
    );
}

#[test]
fn filter_exclude_drops_matching_records() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, err, code) = run(&["filter", "-e", "AF<0.1", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter -e failed: {err}");
    let recs: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(recs.len(), 2);
    assert!(recs.iter().all(|r| !r.starts_with("1\t100\t")));
}

#[test]
fn filter_soft_filter_annotates_failed_records() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, err, code) = run(&[
        "filter",
        "-e",
        "QUAL<50",
        "-s",
        "LowQual",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter soft failed: {err}");
    assert!(
        out.contains("##FILTER=<ID=LowQual,"),
        "missing FILTER header: {out}"
    );
    let pos200 = out
        .lines()
        .find(|l| l.starts_with("1\t200\t"))
        .expect("pos200 record");
    let cols: Vec<&str> = pos200.split('\t').collect();
    assert_eq!(cols[6], "LowQual");
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let cols: Vec<&str> = pos100.split('\t').collect();
    assert_eq!(cols[6], "PASS");
}

#[test]
fn filter_mode_plus_appends_to_existing_filter() {
    let dir = TempDir::new().unwrap();
    let body = VCF.replace("1\t200\t.\tG\tT\t10\tPASS", "1\t200\t.\tG\tT\t10\tOldTag");
    let vcf = write_vcf(&dir, &body);
    let (out, err, code) = run(&[
        "filter",
        "-e",
        "QUAL<50",
        "-s",
        "LowQual",
        "-m",
        "+",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter -m+ failed: {err}");
    let pos200 = out.lines().find(|l| l.starts_with("1\t200\t")).unwrap();
    let cols: Vec<&str> = pos200.split('\t').collect();
    assert_eq!(cols[6], "OldTag;LowQual");
}

#[test]
fn filter_mode_x_resets_pass_for_passing_records() {
    let dir = TempDir::new().unwrap();
    // Simulate a record that already has a non-PASS filter and passes the expr.
    let body = VCF.replace("1\t100\t.\tA\tC\t100\tPASS", "1\t100\t.\tA\tC\t100\tStale");
    let vcf = write_vcf(&dir, &body);
    let (out, err, code) = run(&["filter", "-i", "DP>=10", "-m", "x", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter -mx failed: {err}");
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let cols: Vec<&str> = pos100.split('\t').collect();
    assert_eq!(cols[6], "PASS");
}

#[test]
fn filter_mask_soft_filters_overlapping_records() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, err, code) = run(&[
        "filter",
        "--mask",
        "1:150-250",
        "-s",
        "Masked",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --mask failed: {err}");
    assert!(
        out.contains("##FILTER=<ID=Masked,"),
        "missing mask FILTER header:\n{out}"
    );
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let pos200 = out.lines().find(|l| l.starts_with("1\t200\t")).unwrap();
    let pos300 = out.lines().find(|l| l.starts_with("1\t300\t")).unwrap();
    assert_eq!(pos100.split('\t').nth(6), Some("PASS"));
    assert_eq!(pos200.split('\t').nth(6), Some("Masked"));
    assert_eq!(pos300.split('\t').nth(6), Some("PASS"));
}

#[test]
fn filter_mask_file_and_negated_mask_match_upstream_shapes() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let mask = dir.path().join("mask.bed");
    std::fs::write(&mask, "1\t149\t250\n").unwrap();
    let mask_arg = format!("^{}", mask.display());
    let (out, err, code) = run(&[
        "filter",
        "--mask-file",
        &mask_arg,
        "-s",
        "Outside",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --mask-file ^file failed: {err}");
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let pos200 = out.lines().find(|l| l.starts_with("1\t200\t")).unwrap();
    let pos300 = out.lines().find(|l| l.starts_with("1\t300\t")).unwrap();
    assert_eq!(pos100.split('\t').nth(6), Some("Outside"));
    assert_eq!(pos200.split('\t').nth(6), Some("PASS"));
    assert_eq!(pos300.split('\t').nth(6), Some("Outside"));
}

#[test]
fn filter_mask_overlap_record_span_can_include_deletions() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tACGT\tA\t100\tPASS\t.\n";
    let vcf = write_vcf(&dir, body);
    let (pos_out, pos_err, pos_code) = run(&[
        "filter",
        "--mask",
        "1:103-103",
        "--mask-overlap",
        "0",
        "-s",
        "Masked",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(pos_code, 0, "filter --mask-overlap 0 failed: {pos_err}");
    let pos_record = pos_out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    assert_eq!(pos_record.split('\t').nth(6), Some("PASS"));

    let (span_out, span_err, span_code) = run(&[
        "filter",
        "--mask",
        "1:103-103",
        "--mask-overlap",
        "1",
        "-s",
        "Masked",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(span_code, 0, "filter --mask-overlap 1 failed: {span_err}");
    let span_record = span_out
        .lines()
        .find(|l| l.starts_with("1\t100\t"))
        .unwrap();
    assert_eq!(span_record.split('\t').nth(6), Some("Masked"));
}

#[test]
fn filter_mask_requires_soft_filter() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (_out, err, code) = run(&["filter", "--mask", "1:150-250", vcf.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("--soft-filter is required"),
        "expected soft-filter error, got: {err}"
    );
}

#[test]
fn filter_set_gts_missing_rewrites_failed_site_genotypes() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tC\t100\tPASS\tAC=3;AN=4;DP=10\tGT:DP\t0/1:8\t1|1:9\n\
1\t200\t.\tG\tT\t10\tPASS\tAC=1;AN=4;DP=2\tGT:DP\t0/1:3\t0|0:4\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&["filter", "-i", "DP>=10", "-S", ".", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter -S . failed: {err}");
    let records: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(
        records.len(),
        2,
        "set-GTs should retain failed records:\n{out}"
    );
    let pass_cols: Vec<&str> = records[0].split('\t').collect();
    assert_eq!(pass_cols[9], "0/1:8");
    assert_eq!(pass_cols[10], "1|1:9");
    let fail_cols: Vec<&str> = records[1].split('\t').collect();
    assert_eq!(fail_cols[7], "AC=0;AN=0;DP=2");
    assert_eq!(fail_cols[9], "./.:3");
    assert_eq!(fail_cols[10], ".|.:4");
}

#[test]
fn filter_set_gts_ref_rewrites_failed_site_genotypes() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tC\t100\tPASS\tAC=3;AN=4;DP=10\tGT\t0/1\t1|1\n\
1\t200\t.\tG\tT\t10\tPASS\tAC=1;AN=1;DP=2\tGT\t./1\t.|.\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&[
        "filter",
        "-e",
        "DP<10",
        "--set-GTs=0",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --set-GTs=0 failed: {err}");
    let fail = out
        .lines()
        .find(|l| l.starts_with("1\t200\t"))
        .expect("failed site retained");
    let cols: Vec<&str> = fail.split('\t').collect();
    assert_eq!(cols[7], "AC=0;AN=4;DP=2");
    assert_eq!(cols[9], "0/0");
    assert_eq!(cols[10], "0|0");
}

#[test]
fn filter_set_gts_missing_rewrites_only_failed_samples_for_format_expression() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Sample depth\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
1\t100\t.\tA\tC\t100\tPASS\tAC=4;AN=6\tGT:DP\t0/1:12\t1/1:3\t0|1:9\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&[
        "filter",
        "-i",
        "FMT/DP>=10",
        "-S",
        ".",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter FORMAT-scoped -S failed: {err}");
    let record = out
        .lines()
        .find(|line| line.starts_with("1\t100\t"))
        .expect("record retained");
    let cols: Vec<&str> = record.split('\t').collect();
    assert_eq!(cols[7], "AC=1;AN=2");
    assert_eq!(cols[9], "0/1:12", "passing sample should be untouched");
    assert_eq!(cols[10], "./.:3", "failed sample should be rewritten");
    assert_eq!(cols[11], ".|.:9", "failed sample should preserve separator");
}

#[test]
fn filter_set_gts_rejects_unknown_target() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (_out, err, code) = run(&["filter", "-S", "1", vcf.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("argument to -S not recognised"),
        "expected set-GTs parse error, got: {err}"
    );
}

#[test]
fn filter_snp_gap_soft_filters_snps_near_indels() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tC\t100\tPASS\t.\n\
1\t105\t.\tG\tGA\t100\tPASS\t.\n\
1\t115\t.\tT\tC\t100\tPASS\t.\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&["filter", "-g", "5", "-s", "Keep", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter -g failed: {err}");
    assert!(
        out.contains("##FILTER=<ID=SnpGap,"),
        "missing SnpGap header:\n{out}"
    );
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let pos105 = out.lines().find(|l| l.starts_with("1\t105\t")).unwrap();
    let pos115 = out.lines().find(|l| l.starts_with("1\t115\t")).unwrap();
    assert_eq!(pos100.split('\t').nth(6), Some("SnpGap"));
    assert_eq!(pos105.split('\t').nth(6), Some("PASS"));
    assert_eq!(pos115.split('\t').nth(6), Some("PASS"));
}

#[test]
fn filter_snp_gap_hard_drops_tagged_snps() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tC\t100\tPASS\t.\n\
1\t105\t.\tG\tGA\t100\tPASS\t.\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&["filter", "--SnpGap", "5", vcf.to_str().unwrap()]);
    assert_eq!(code, 0, "filter --SnpGap failed: {err}");
    assert!(!out.lines().any(|l| l.starts_with("1\t100\t")));
    assert!(out.lines().any(|l| l.starts_with("1\t105\t")));
}

#[test]
fn filter_indel_gap_keeps_highest_qual_indel_in_cluster() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100\t.\tA\tAT\t10\tPASS\t.\n\
1\t103\t.\tG\tGA\t90\tPASS\t.\n\
1\t120\t.\tC\tCA\t30\tPASS\t.\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&[
        "filter",
        "--IndelGap",
        "5",
        "-s",
        "Keep",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --IndelGap failed: {err}");
    assert!(
        out.contains("##FILTER=<ID=IndelGap,"),
        "missing IndelGap header:\n{out}"
    );
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let pos103 = out.lines().find(|l| l.starts_with("1\t103\t")).unwrap();
    let pos120 = out.lines().find(|l| l.starts_with("1\t120\t")).unwrap();
    assert_eq!(pos100.split('\t').nth(6), Some("IndelGap"));
    assert_eq!(pos103.split('\t').nth(6), Some("PASS"));
    assert_eq!(pos120.split('\t').nth(6), Some("PASS"));
}

#[test]
fn filter_indel_gap_ties_by_first_record_and_falls_back_to_ac() {
    let dir = TempDir::new().unwrap();
    let body = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t100\t.\tA\tAT\t50\tPASS\t.\tGT\t0/1\t0/0\n\
1\t103\t.\tG\tGA\t50\tPASS\t.\tGT\t0/1\t0/1\n\
1\t200\t.\tC\tCA\t.\tPASS\t.\tGT\t0/1\t0/0\n\
1\t203\t.\tT\tTA\t.\tPASS\t.\tGT\t0/1\t1/1\n";
    let vcf = write_vcf(&dir, body);
    let (out, err, code) = run(&[
        "filter",
        "--IndelGap",
        "5",
        "-s",
        "Keep",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --IndelGap failed: {err}");
    let pos100 = out.lines().find(|l| l.starts_with("1\t100\t")).unwrap();
    let pos103 = out.lines().find(|l| l.starts_with("1\t103\t")).unwrap();
    let pos200 = out.lines().find(|l| l.starts_with("1\t200\t")).unwrap();
    let pos203 = out.lines().find(|l| l.starts_with("1\t203\t")).unwrap();

    assert_eq!(
        pos100.split('\t').nth(6),
        Some("PASS"),
        "first equal-QUAL indel should pass:\n{out}"
    );
    assert_eq!(
        pos103.split('\t').nth(6),
        Some("IndelGap"),
        "second equal-QUAL indel should be filtered:\n{out}"
    );
    assert_eq!(
        pos200.split('\t').nth(6),
        Some("IndelGap"),
        "lower-AC missing-QUAL indel should be filtered:\n{out}"
    );
    assert_eq!(
        pos203.split('\t').nth(6),
        Some("PASS"),
        "higher-AC missing-QUAL indel should pass:\n{out}"
    );
}

#[test]
fn filter_gap_rejects_bad_arguments() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (_out, snp_err, snp_code) = run(&["filter", "-g", "5:bogus", vcf.to_str().unwrap()]);
    assert_ne!(snp_code, 0);
    assert!(
        snp_err.contains("Could not parse \"bogus\""),
        "got: {snp_err}"
    );
    let (_out, indel_err, indel_code) = run(&["filter", "-G", "abc", vcf.to_str().unwrap()]);
    assert_ne!(indel_code, 0);
    assert!(
        indel_err.contains("Could not parse argument: --IndelGap abc"),
        "got: {indel_err}"
    );
}

#[test]
fn filter_region_restricts_to_named_chrom() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (out, _err, code) = run(&["filter", "-r", "1:150-250", vcf.to_str().unwrap()]);
    assert_eq!(code, 0);
    let recs: Vec<&str> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(recs.len(), 1);
    assert!(recs[0].starts_with("1\t200\t"));
}

#[test]
fn filter_write_index_creates_csi_for_bgzf_vcf() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let out_path = dir.path().join("filtered.vcf.gz");
    let (_out, err, code) = run(&[
        "filter",
        "-i",
        "DP>=10",
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "-W",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter -W -O z failed: {err}");
    assert!(
        dir.path().join("filtered.vcf.gz.csi").exists(),
        "CSI index not created"
    );
}

#[test]
fn filter_write_index_tbi_creates_tbi_for_bgzf_vcf() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let out_path = dir.path().join("filtered.vcf.gz");
    let (_out, err, code) = run(&[
        "filter",
        "-i",
        "DP>=10",
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "--write-index=tbi",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --write-index=tbi failed: {err}");
    assert!(
        dir.path().join("filtered.vcf.gz.tbi").exists(),
        "TBI index not created"
    );
}

#[test]
fn filter_write_index_requires_output_file() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (_out, err, code) = run(&["filter", "-W", vcf.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("-W requires an output file"),
        "expected output-file error, got: {err}"
    );
}

#[test]
fn filter_threads_writes_bgzf_vcf_output() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let out_path = dir.path().join("filtered.vcf.gz");
    let (_out, err, code) = run(&[
        "filter",
        "-i",
        "DP>=10",
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "z",
        "--threads",
        "2",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --threads -O z failed: {err}");
    let mut decoded = String::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(&out_path).unwrap())
        .read_to_string(&mut decoded)
        .unwrap();
    let records: Vec<&str> = decoded
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert_eq!(records.len(), 2);
}

#[test]
fn filter_threads_writes_bcf_output() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let out_path = dir.path().join("filtered.bcf");
    let (_out, err, code) = run(&[
        "filter",
        "-i",
        "DP>=10",
        "-o",
        out_path.to_str().unwrap(),
        "-O",
        "b",
        "--threads=2",
        vcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "filter --threads -O b failed: {err}");
    let bytes = std::fs::read(&out_path).unwrap();
    assert_eq!(&bytes[..4], &[0x1f, 0x8b, 0x08, 0x04]);
}

#[test]
fn filter_threads_rejects_non_integer_argument() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (_out, err, code) = run(&["filter", "--threads", "abc", vcf.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(
        err.contains("Could not parse argument: --threads abc"),
        "got: {err}"
    );
}

#[test]
fn filter_sample_fraction_matches_sample_count_on_upstream_fixture() {
    for expr in ["N_PASS(DP>32)=1", "F_PASS(DP>32)=0.5"] {
        let (out, err, code) = run(&[
            "filter",
            "--no-version",
            "-i",
            expr,
            "../../bcftools/test/filter.2.vcf",
        ]);
        assert_eq!(code, 0, "filter -i {expr} failed: {err}");
        let records: Vec<Vec<&str>> = out
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(|line| line.split('\t').collect())
            .collect();
        assert_eq!(records.len(), 3, "unexpected records for {expr}:\n{out}");
        assert_eq!(records[0][1], "3062915");
        assert_eq!(records[0][9], "0/1:25:35:-20,-5,-20");
        assert_eq!(records[0][10], "0/1:45:11:-20,-5,-20");
        assert_eq!(records[1][1], "3106154");
        assert_eq!(records[1][9], "0/1:245:32");
        assert_eq!(records[1][10], "0/1:25:300");
        assert_eq!(records[2][1], "3106154");
        assert_eq!(records[2][9], "0/1:25:12");
        assert_eq!(records[2][10], "0/1:245:310");
    }
}

#[test]
fn filter_missing_fraction_matches_upstream_filter_28_fixture() {
    for expr in [
        "F_MISSING>=1/5",
        "F_MISSING>=0.2",
        "F_PASS(GT==\"mis\")>=1/5",
        "F_PASS(GT==\"mis\")>=0.2",
    ] {
        let (out, err, code) = run(&[
            "filter",
            "--no-version",
            "-i",
            expr,
            "../../bcftools/test/filter.6.vcf",
        ]);
        assert_eq!(code, 0, "filter -i {expr} failed: {err}");
        let records: Vec<Vec<&str>> = out
            .lines()
            .filter(|line| !line.starts_with('#') && !line.is_empty())
            .map(|line| line.split('\t').collect())
            .collect();
        assert_eq!(records.len(), 1, "unexpected records for {expr}:\n{out}");
        assert_eq!(records[0][1], "3162007");
    }
}

#[test]
fn filter_single_amp_pipe_and_set_gts_match_upstream_filter_2_fixture() {
    let expected = std::fs::read_to_string("../../bcftools/test/filter.2.out").unwrap();
    let (out, err, code) = run(&[
        "filter",
        "--no-version",
        "-e",
        "QUAL==59.2 || (INDEL=0 & (FMT/GQ=25 | FMT/DP=10))",
        "-sModified",
        "-S.",
        "../../bcftools/test/filter.2.vcf",
    ]);

    assert_eq!(code, 0, "filter.2 fixture command failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn filter_no_args_prints_usage() {
    let (_out, err, code) = run(&["filter"]);
    assert_ne!(code, 0);
    assert!(err.contains("Usage:"), "no usage in stderr: {err}");
}

#[test]
fn filter_unknown_output_type_errors() {
    let dir = TempDir::new().unwrap();
    let vcf = write_vcf(&dir, VCF);
    let (_out, err, code) = run(&["filter", "-O", "Q", vcf.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("not recognised"), "got: {err}");
}
