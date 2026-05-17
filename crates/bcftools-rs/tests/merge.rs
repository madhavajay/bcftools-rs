//! End-to-end tests for `bcftools_rs::commands::merge`.

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
    let stdout = String::from_utf8(out.stdout).unwrap();
    let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
    (stdout, stderr, out.status.code().unwrap_or(-1))
}

fn write_vcf(dir: &TempDir, name: &str, sample: &str, gt: &str) -> PathBuf {
    let path = dir.path().join(name);
    let body = format!(
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n\
1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t{gt}\n"
    );
    std::fs::write(&path, body).unwrap();
    path
}

#[test]
fn merge_combines_same_site_samples() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b = write_vcf(&dir, "b.vcf", "SAMPLE_B", "1/1");

    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge failed: {err}");
    assert!(
        out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\tSAMPLE_B"),
        "merged header missing: {out}"
    );
    assert!(
        out.contains("1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1"),
        "merged record missing: {out}"
    );
}

#[test]
fn merge_rejects_duplicate_sample_without_force() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "DUP", "0/1");
    let b = write_vcf(&dir, "b.vcf", "DUP", "1/1");

    let (_out, err, code) = run(&[
        "merge",
        "--no-version",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_ne!(code, 0, "expected duplicate-sample failure, got success");
    assert!(
        err.contains("duplicate sample name") && err.contains("DUP"),
        "stderr should mention duplicate sample 'DUP': {err}"
    );
}

#[test]
fn merge_force_samples_prefixes_duplicates() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "DUP", "0/1");
    let b = write_vcf(&dir, "b.vcf", "DUP", "1/1");

    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "--force-samples",
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge --force-samples failed: {err}");
    assert!(
        out.contains("FORMAT\tDUP\t2:DUP"),
        "expected DUP and 2:DUP columns: {out}"
    );
}

#[test]
fn merge_force_samples_repeats_prefix_until_unique_like_upstream_fixture() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "--force-samples",
        "../../bcftools/test/merge.7.a.vcf",
        "../../bcftools/test/merge.7.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.7 --force-samples fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.9.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_unions_missing_sites_with_missing_sample_values() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b_path = dir.path().join("b.vcf");
    std::fs::write(
        &b_path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_B\n\
1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t1/1\n",
    )
    .unwrap();

    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        a.to_str().unwrap(),
        b_path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge should union missing sites: {err}");
    assert!(
        out.contains("1\t2\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t./."),
        "missing SAMPLE_B value should be synthesized at first site: {out}"
    );
    assert!(
        out.contains("1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t./.\t1/1"),
        "missing SAMPLE_A value should be synthesized at second site: {out}"
    );
}

#[test]
fn merge_writes_bgzf_vcf_output() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b = write_vcf(&dir, "b.vcf", "SAMPLE_B", "1/1");
    let out_path = dir.path().join("merged.vcf.gz");

    let (_stdout, err, code) = run(&[
        "merge",
        "--no-version",
        "-Oz",
        "-o",
        out_path.to_str().unwrap(),
        a.to_str().unwrap(),
        b.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "merge -Oz failed: {err}");

    let bytes = std::fs::read(&out_path).unwrap();
    assert!(
        bytes.starts_with(&[0x1f, 0x8b]),
        "output should start with gzip magic: {:?}",
        &bytes[..bytes.len().min(4)]
    );
}

#[test]
fn merge_reads_file_list() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");
    let b = write_vcf(&dir, "b.vcf", "SAMPLE_B", "1/1");
    let list = dir.path().join("inputs.txt");
    std::fs::write(&list, format!("{}\n{}\n", a.display(), b.display())).unwrap();

    let (out, err, code) = run(&["merge", "--no-version", "-l", list.to_str().unwrap()]);
    assert_eq!(code, 0, "merge -l failed: {err}");
    assert!(out.contains("SAMPLE_A\tSAMPLE_B"), "{out}");
}

#[test]
fn merge_noidx_fixture_matches_upstream_text_output() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "--no-index",
        "../../bcftools/test/merge.noidx.a.vcf",
        "../../bcftools/test/merge.noidx.b.vcf",
        "../../bcftools/test/merge.noidx.c.vcf",
    ]);
    assert_eq!(code, 0, "merge noidx fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.noidx.abc.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_force_single_fixture_matches_upstream_text_output() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "--force-single",
        "../../bcftools/test/merge.LPL.a.vcf",
    ]);
    assert_eq!(code, 0, "merge --force-single fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.LPL.0.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_sites_only_alt_union_matches_upstream_fixture() {
    for extra_args in [Vec::<&str>::new(), vec!["-i", "AN:sum,AC:sum"]] {
        let mut args = vec![
            "merge",
            "--no-version",
            "../../bcftools/test/merge.8.a.vcf",
            "../../bcftools/test/merge.8.b.vcf",
        ];
        args.splice(2..2, extra_args);
        let (out, err, code) = run(&args);
        assert_eq!(code, 0, "merge.8 fixture failed for {args:?}: {err}");

        let expected = std::fs::read_to_string("../../bcftools/test/merge.8.out").unwrap();
        assert_eq!(out, expected, "arguments {args:?}");
    }
}

#[test]
fn merge_sampled_with_sites_only_alt_union_matches_upstream_fixture() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "../../bcftools/test/merge.9.a.vcf",
        "../../bcftools/test/merge.9.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.9 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.9.1.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_sampled_sites_only_alt_union_info_rules_matches_upstream_fixture() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "-i",
        "AN:sum,AC:sum",
        "../../bcftools/test/merge.9.a.vcf",
        "../../bcftools/test/merge.9.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.9 -i fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.9.2.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_mode_none_keeps_conflicting_same_position_records_as_separate_rows() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "-m",
        "none",
        "../../bcftools/test/merge.10.a.vcf",
        "../../bcftools/test/merge.10.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.10 -m none fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.10.1.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_mode_both_and_snp_ins_del_match_upstream_text_fixtures() {
    for (mode, fixture) in [
        ("both", "../../bcftools/test/merge.10.2.out"),
        ("snp-ins-del", "../../bcftools/test/merge.10.3.out"),
    ] {
        let (out, err, code) = run(&[
            "merge",
            "--no-version",
            "-m",
            mode,
            "../../bcftools/test/merge.10.a.vcf",
            "../../bcftools/test/merge.10.b.vcf",
        ]);
        assert_eq!(code, 0, "merge.10 -m {mode} fixture failed: {err}");

        let expected = std::fs::read_to_string(fixture).unwrap();
        assert_eq!(out, expected, "mode {mode}");
    }
}

#[test]
fn merge_ad_vector_allele_union_matches_upstream_fixture() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "../../bcftools/test/merge.11.a.vcf",
        "../../bcftools/test/merge.11.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.11 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.11.1.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_non_ref_symbolic_allele_union_matches_upstream_fixture() {
    for mode in ["none", "both"] {
        let (out, err, code) = run(&[
            "merge",
            "--no-version",
            "--merge",
            mode,
            "../../bcftools/test/merge.12.a.vcf",
            "../../bcftools/test/merge.12.b.vcf",
        ]);
        assert_eq!(code, 0, "merge.12 --merge {mode} fixture failed: {err}");

        let expected = std::fs::read_to_string("../../bcftools/test/merge.12.1.out").unwrap();
        assert_eq!(out, expected, "mode {mode}");
    }
}

#[test]
fn merge_info_af_join_matches_upstream_fixture() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "-i",
        "AF:join",
        "../../bcftools/test/merge.join.a.vcf",
        "../../bcftools/test/merge.join.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.join -i AF:join fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.join.1.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_symbolic_records_use_highest_input_fileformat() {
    let (out, err, code) = run(&[
        "merge",
        "--no-version",
        "../../bcftools/test/merge.symbolic.1.a.vcf",
        "../../bcftools/test/merge.symbolic.1.b.vcf",
    ]);
    assert_eq!(code, 0, "merge.symbolic.1 fixture failed: {err}");

    let expected = std::fs::read_to_string("../../bcftools/test/merge.symbolic.1.1.out").unwrap();
    assert_eq!(out, expected);
}

#[test]
fn merge_rejects_single_input() {
    let dir = TempDir::new().unwrap();
    let a = write_vcf(&dir, "a.vcf", "SAMPLE_A", "0/1");

    let (_out, err, code) = run(&["merge", "--no-version", a.to_str().unwrap()]);
    assert_ne!(code, 0, "expected single-input rejection");
    assert!(
        err.contains("at least two") || err.contains("expected at least"),
        "stderr should request at least two inputs: {err}"
    );
}
