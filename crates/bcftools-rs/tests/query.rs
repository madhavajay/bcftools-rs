//! End-to-end tests for `bcftools_rs::commands::query`.

use std::path::PathBuf;
use std::process::Command;

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

#[test]
fn query_list_samples_from_vcf() {
    let path = fixture_path("annotate2.vcf");
    let (out, err, code) = run(&["query", "-l", path.to_str().unwrap()]);
    assert_eq!(code, 0, "query -l failed: {err}");
    assert_eq!(out, "A\nB\nC\n");
}

#[test]
fn query_list_samples_from_bcf() {
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("samples.vcf");
    let bcf = dir.path().join("annotate2.bcf");
    std::fs::write(
        &input,
        "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n\
1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t0/0\t1/1\n",
    )
    .unwrap();

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["query", "--list-samples", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "query --list-samples BCF failed: {err}");
    assert_eq!(out, "A\nB\nC\n");
}

#[test]
fn query_format_core_fields_from_vcf() {
    let path = fixture_path("annotate2.vcf");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM\\t%POS\\t%REF\\t%ALT\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f failed: {err}");
    assert!(out.starts_with("1\t3000001\tC\tT\n"));
    assert!(out.lines().count() > 1);
}

#[test]
fn query_samples_file_filters_list_samples_in_header_order() {
    let path = fixture_path("query.smpl.vcf");
    let samples = fixture_path("query.smpl.txt");
    let (out, err, code) = run(&[
        "query",
        "-l",
        "-S",
        samples.to_str().unwrap(),
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -l -S failed: {err}");
    assert_eq!(out, "00\n11\n");
}

#[test]
fn query_samples_file_reorders_format_loops() {
    let path = fixture_path("query.smpl.vcf");
    let samples = fixture_path("query.smpl.txt");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "[%SAMPLE %GT\\n]",
        "-S",
        samples.to_str().unwrap(),
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f -S failed: {err}");
    assert_eq!(out, "11 1/1\n00 0/0\n");
}

#[test]
fn query_samples_file_exclusion_filters_format_loops() {
    let path = fixture_path("query.smpl.vcf");
    let samples = fixture_path("query.smpl.11.txt");
    let excluded = format!("^{}", samples.display());
    let (out, err, code) = run(&[
        "query",
        "-f",
        "[%SAMPLE %GT\\n]",
        "-S",
        &excluded,
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f -S ^ failed: {err}");
    assert_eq!(out, "00 0/0\n");
}

#[test]
fn query_print_header_adds_indexed_column_names() {
    let path = fixture_path("query.header.vcf");
    let (out, err, code) = run(&[
        "query",
        "-H",
        "-f",
        "%CHROM %POS[ %SAMPLE %DP %GT]\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -H failed: {err}");
    assert_eq!(
        out,
        "#[1]CHROM [2]POS [3]C:SAMPLE [4]C:DP [5]C:GT [6]D:SAMPLE [7]D:DP [8]D:GT\n\
4 3258449 C 1 1/1 D 0 0/0\n"
    );
}

#[test]
fn query_print_header_twice_omits_column_indices() {
    let path = fixture_path("query.header.vcf");
    let (out, err, code) = run(&[
        "query",
        "-HH",
        "-f",
        "%CHROM %POS[ %SAMPLE][ %DP][ %GT]",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -HH failed: {err}");
    assert_eq!(
        out,
        "#CHROM POS C:SAMPLE D:SAMPLE C:DP D:DP C:GT D:GT\n\
4 3258449 C D 1 0 1/1 0/0"
    );
}

#[test]
fn query_regions_file_filters_records() {
    let path = fixture_path("regions.vcf");
    let regions = fixture_path("regions.tab");
    let expected = std::fs::read_to_string(fixture_path("regions.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM %POS %REF,%ALT\\n",
        "-R",
        regions.to_str().unwrap(),
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -R failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_inline_regions_filter_records() {
    let path = fixture_path("regions.vcf");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM %POS %REF,%ALT\\n",
        "-r",
        "1:3062915-3106154,2:3199815-3199815",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -r failed: {err}");
    assert_eq!(
        out,
        "1 3062915 GTT,G\n\
1 3062915 G,T\n\
1 3106154 CA,C\n\
1 3106154 C,T,CT\n\
2 3199815 C,T\n"
    );
}

#[test]
fn query_targets_file_filters_records() {
    let path = fixture_path("regions.vcf");
    let targets = fixture_path("regions.tab");
    let expected = std::fs::read_to_string(fixture_path("regions.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM %POS %REF,%ALT\\n",
        "-T",
        targets.to_str().unwrap(),
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -T failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_targets_exclusion_filters_records() {
    let path = fixture_path("regions.vcf");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM %POS %REF,%ALT\\n",
        "-t",
        "^1:3062915-3184885",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -t ^ failed: {err}");
    assert_eq!(
        out,
        "2 3199812 G,T\n\
2 3199815 C,T\n\
3 3212016 C,A\n\
3 3212026 C,A\n\
3 3212036 C,A\n"
    );
}

#[test]
fn query_include_filters_core_and_info_fields() {
    let path = fixture_path("annotate2.vcf");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS %FILTER %IINT\\n",
        "-i",
        "IINT=11 && FILTER=\"PASS\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i failed: {err}");
    assert_eq!(out, "3000001 PASS 11\n");
}

#[test]
fn query_exclude_filters_string_info_fields() {
    let path = fixture_path("query.string.vcf");
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%CLNREVSTAT\\n",
        "-e",
        "CLNREVSTAT=\"criteria_provided,_single_submitter\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -e failed: {err}");
    assert_eq!(
        out,
        "865568\tcriteria_provided,_conflicting_interpretations\n\
865628\tcriteria_provided,_multiple_submitters,_no_conflicts\n"
    );
}

#[test]
fn query_computed_n_alt_filter_matches_upstream_fixture() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.6.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS %REF %ALT\\n",
        "-i",
        "N_ALT=2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i N_ALT failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_computed_n_samples_filter_matches_upstream_fixture() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.7.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS %AN\\n",
        "-i",
        "AN!=2*N_SAMPLES",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i N_SAMPLES failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_format_token_respects_sample_reordering() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.64.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM\\t%POS\\t%INFO\\t%FORMAT\\n",
        "-s",
        "D,C",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query %FORMAT failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_computed_type_filter_matches_upstream_exact_fixture() {
    let path = fixture_path("query.filter-type.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.26.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\n",
        "-i",
        "type=\"snp\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i type=\"snp\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_computed_type_filter_matches_upstream_regex_fixture() {
    let path = fixture_path("query.filter-type.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.27.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\n",
        "-i",
        "type~\"snp\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i type~\"snp\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_info_type_still_prefers_info_namespace() {
    let path = fixture_path("query.filter-type.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.67.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\n",
        "-i",
        "INFO/TYPE=\"xxx\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i INFO/TYPE failed: {err}");
    assert_eq!(out, expected);
}
