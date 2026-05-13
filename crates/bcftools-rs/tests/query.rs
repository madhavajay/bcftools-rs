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
4 3258449 C D 1 0 1/1 0/0\n"
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
fn query_alt_vector_regex_filter_matches_upstream_fixture() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.5.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS %REF %ALT\\n",
        "-i",
        "REF~\"C\" && ALT[*]~\"CT\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i ALT[*] regex failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_vector_missing_predicates_match_upstream_fixtures() {
    let path = fixture_path("query.filter.15.vcf");
    let expected_missing = std::fs::read_to_string(fixture_path("query.filter.15.1.out")).unwrap();
    let expected_present = std::fs::read_to_string(fixture_path("query.filter.15.2.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%TAG",
        "-i",
        "TAG[*]=\".\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i TAG[*]=missing failed: {err}");
    assert_eq!(out, expected_missing);

    let (out, err, code) = run(&[
        "query",
        "-f",
        "%TAG",
        "-i",
        "TAG[*]!=\".\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i TAG[*]!=missing failed: {err}");
    assert_eq!(out, expected_present);

    let (out, err, code) = run(&[
        "query",
        "-f",
        "%TAG",
        "-i",
        "TAG[*]~\"\\.\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i TAG[*] regex missing failed: {err}");
    assert_eq!(out, expected_missing);

    let (out, err, code) = run(&[
        "query",
        "-f",
        "%TAG",
        "-i",
        "TAG[*]!~\"\\.\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i TAG[*] not-regex missing failed: {err}");
    assert_eq!(out, expected_present);
}

#[test]
fn query_alt_scalar_filter_matches_any_alternate_allele() {
    let path = fixture_path("query.filter.4.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.55.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT[\\t%GT]\\n",
        "-e",
        "TYPE!=\"snp\" || ALT=\"*\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -e ALT=\"*\" failed: {err}");
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
fn query_computed_type_filter_matches_upstream_negated_exact_fixture() {
    let path = fixture_path("query.filter-type.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.28.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\n",
        "-i",
        "type!=\"snp\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i type!=\"snp\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_computed_type_filter_matches_upstream_negated_regex_fixture() {
    let path = fixture_path("query.filter-type.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.29.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\n",
        "-i",
        "type!~\"snp\"",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i type!~\"snp\" failed: {err}");
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

#[test]
fn query_percent_ilen_filter_uses_computed_length() {
    let path = fixture_path("query.filter.8.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.69.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\t%ILEN\\n",
        "-i",
        "%ILEN==1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i %ILEN failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_bare_ilen_filter_uses_info_tag() {
    let path = fixture_path("query.filter.8.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.70.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS\\t%REF\\t%ALT\\t%ILEN\\n",
        "-i",
        "ILEN==1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i ILEN failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_filter_exact_match_matches_upstream_fixture() {
    let path = fixture_path("filter.11.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.76.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-i",
        "FILTER=\"A\"",
        "-f",
        "%POS %FILTER\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i FILTER=\"A\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_filter_negated_exact_match_matches_upstream_fixture() {
    let path = fixture_path("filter.11.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.77.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-i",
        "FILTER!=\"A\"",
        "-f",
        "%POS %FILTER\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i FILTER!=\"A\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_filter_id_match_matches_upstream_fixture() {
    let path = fixture_path("filter.11.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.78.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-i",
        "FILTER~\"A\"",
        "-f",
        "%POS %FILTER\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i FILTER~\"A\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_filter_negated_id_match_matches_upstream_fixture() {
    let path = fixture_path("filter.11.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.79.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-i",
        "FILTER!~\"A\"",
        "-f",
        "%POS %FILTER\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i FILTER!~\"A\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_missing_integer_info_matches_upstream_fixtures() {
    let path = fixture_path("missing.vcf");
    let expected_missing = std::fs::read_to_string(fixture_path("query.18.out")).unwrap();
    let expected_present = std::fs::read_to_string(fixture_path("query.19.out")).unwrap();

    let (out, err, code) = run(&[
        "query",
        "-i",
        "IINT=\".\"",
        "-f",
        "%POS %IINT\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i IINT missing failed: {err}");
    assert_eq!(out, expected_missing);

    let (out, err, code) = run(&[
        "query",
        "-i",
        "IINT!=\".\"",
        "-f",
        "%POS %IINT\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i IINT present failed: {err}");
    assert_eq!(out, expected_present);
}

#[test]
fn query_missing_float_info_matches_upstream_fixtures() {
    let path = fixture_path("missing.vcf");
    let expected_missing = std::fs::read_to_string(fixture_path("query.20.out")).unwrap();
    let expected_present = std::fs::read_to_string(fixture_path("query.21.out")).unwrap();

    let (out, err, code) = run(&[
        "query",
        "-i",
        "IFLT=\".\"",
        "-f",
        "%POS %IFLT\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i IFLT missing failed: {err}");
    assert_eq!(out, expected_missing);

    let (out, err, code) = run(&[
        "query",
        "-i",
        "IFLT!=\".\"",
        "-f",
        "%POS %IFLT\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i IFLT present failed: {err}");
    assert_eq!(out, expected_present);
}

#[test]
fn query_missing_string_info_matches_upstream_fixtures() {
    let path = fixture_path("missing.vcf");
    let expected_missing = std::fs::read_to_string(fixture_path("query.22.out")).unwrap();
    let expected_present = std::fs::read_to_string(fixture_path("query.23.out")).unwrap();

    let (out, err, code) = run(&[
        "query",
        "-i",
        "ISTR=\".\"",
        "-f",
        "%POS %ISTR\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i ISTR missing failed: {err}");
    assert_eq!(out, expected_missing);

    let (out, err, code) = run(&[
        "query",
        "-i",
        "ISTR!=\".\"",
        "-f",
        "%POS %ISTR\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i ISTR present failed: {err}");
    assert_eq!(out, expected_present);
}

#[test]
fn query_filter_exact_missing_fixture_matches_upstream() {
    let path = fixture_path("missing.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.24.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-i",
        "FILTER=\"q11\"",
        "-f",
        "%POS %ISTR\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i FILTER=\"q11\" failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_line_token_matches_upstream_fixture() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.25.out")).unwrap();
    let (out, err, code) = run(&["query", "-f", "%LINE", path.to_str().unwrap()]);
    assert_eq!(code, 0, "query -f %LINE failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_info_namespace_tokens_match_upstream_fixture() {
    let path = fixture_path("query.3.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.3.1.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM %POS %ID %REF %ALT %QUAL %FILTER \\t %INFO/CHROM %INFO/POS %INFO/ID %INFO/REF %INFO/ALT %INFO/QUAL %INFO/FILTER",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f INFO/TAG failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_sample_loop_prefers_format_namespace_for_bare_tags() {
    let path = fixture_path("query.3.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.3.2.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "[ %CHROM] \\t [ %POS] \\t [ %ID] \\t [ %REF] \\t [ %ALT] \\t [ %QUAL] \\t [ %FILTER]",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f sample-loop bare TAG failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_slash_prefix_forces_record_namespace_in_sample_loop() {
    let path = fixture_path("query.3.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.3.3.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "[ %/CHROM] \\t [ %/POS] \\t [ %/ID] \\t [ %/REF] \\t [ %/ALT] \\t [ %/QUAL] \\t [ %/FILTER]",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f %/TAG failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_n_pass_formatter_counts_selected_samples_matching_predicate() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.75.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%CHROM:%POS\\t%N_PASS(GT=\"alt\" & GQ>110)\\t[\\t%GT]\\t[\\t%GQ]\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -f %N_PASS failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_n_pass_filter_counts_numeric_format_predicates() {
    let path = fixture_path("query.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.63.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "[%POS\\t%SAMPLE\\t%GQ\\n]",
        "-i",
        "N_PASS(GQ<20)==1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i N_PASS(GQ<20) failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_n_pass_filter_counts_gt_class_predicates() {
    let path = fixture_path("query.filter.11.vcf");
    let include_expected = std::fs::read_to_string(fixture_path("query.80.out")).unwrap();
    let exclude_expected = std::fs::read_to_string(fixture_path("query.81.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "[%POS\\t%SAMPLE\\t%GT\\n]",
        "-i",
        "N_PASS(GT=\"alt\")==1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i N_PASS(GT=alt) failed: {err}");
    assert_eq!(out, include_expected);

    let (out, err, code) = run(&[
        "query",
        "-f",
        "[%POS\\t%SAMPLE\\t%GT\\n]",
        "-e",
        "N_PASS(GT=\"alt\")==1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -e N_PASS(GT=alt) failed: {err}");
    assert_eq!(out, exclude_expected);
}

#[test]
fn query_count_filter_counts_info_vector_values() {
    let path = fixture_path("query.filter.10.vcf");
    let numeric_expected = std::fs::read_to_string(fixture_path("query.73.out")).unwrap();
    let string_expected = std::fs::read_to_string(fixture_path("query.74.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS  %NUM_TAG\\n",
        "-i",
        "COUNT(INFO/NUM_TAG)=2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i COUNT(INFO/NUM_TAG) failed: {err}");
    assert_eq!(out, numeric_expected);

    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS  %STR_TAG\\n",
        "-i",
        "COUNT(INFO/STR_TAG)=2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i COUNT(INFO/STR_TAG) failed: {err}");
    assert_eq!(out, string_expected);
}

#[test]
fn query_count_filter_counts_gt_class_predicates() {
    let path = fixture_path("query.filter.3.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.53.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-f",
        "%POS[\\t%GT]\\n",
        "-i",
        "COUNT(GT=\"het\")=1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i COUNT(GT=het) failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn query_modulo_filter_evaluates_format_values() {
    let path = fixture_path("filter.10.vcf");
    let expected = std::fs::read_to_string(fixture_path("query.91.out")).unwrap();
    let (out, err, code) = run(&[
        "query",
        "-i",
        "DP%10==2",
        "-f",
        "[ %DP]\\n",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "query -i DP%10 failed: {err}");
    assert_eq!(out, expected);
}
