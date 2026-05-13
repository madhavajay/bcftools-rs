//! End-to-end tests for `bcftools_rs::commands::view`.

use std::io::Read as _;
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

#[test]
fn view_text_round_trip_emits_all_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    // The fixture has 22 records (lines 11..=32 of the file).
    assert_eq!(record_lines.len(), 21);
    // First record is `1\t105\t.\tTAAACCCTA\t...`
    assert!(record_lines[0].starts_with("1\t105\t"));
}

#[test]
fn view_reads_plain_vcf_from_stdin_without_filename() {
    let input = std::fs::read(fixture_path("aa.vcf")).unwrap();
    let (out, err, code) = run_with_stdin(&["view", "--no-version"], &input);
    assert_eq!(code, 0, "view stdin failed: {err}");
    let record_lines: Vec<_> = out
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect();
    assert_eq!(record_lines.len(), 21);
}

#[test]
fn view_header_only_drops_records() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let record_lines: Vec<_> = out
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .collect();
    assert!(record_lines.is_empty());
    assert!(out.contains("#CHROM\t"));
}

#[test]
fn view_no_header_drops_header() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-H", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    let header_lines: Vec<_> = out.lines().filter(|l| l.starts_with('#')).collect();
    assert!(header_lines.is_empty(), "header leaked: {header_lines:?}");
    let record_lines: Vec<_> = out.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(record_lines.len(), 21);
}

#[test]
fn view_region_filters_by_chrom_and_position_interval() {
    let path = fixture_path("aa.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        path.to_str().unwrap(),
        "20:80-95",
    ]);
    assert_eq!(code, 0, "view region failed: {err}");
    let records: Vec<_> = out.lines().filter(|line| !line.is_empty()).collect();
    assert_eq!(
        records,
        [
            "20\t81\t.\tA\tC\t999\tPASS\t.",
            "20\t84\t.\tG\tT\t999\tPASS\t.",
            "20\t95\t.\tT\tA\t999\tPASS\t.",
            "20\t95\t.\tTCACCG\tAAAAAA\t999\tPASS\t.",
        ]
    );
}

#[test]
fn view_regions_option_filters_comma_separated_regions() {
    let path = fixture_path("regions.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-r",
        "1:3062915,2:3199815",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -r failed: {err}");
    let records: Vec<_> = out.lines().filter(|line| !line.is_empty()).collect();
    assert_eq!(
        records,
        [
            "1\t3062915\t.\tGTT\tG\t1806\tq10\tDP=35;DP4=1,2,3,4;AN=2;AC=1\tGT:GQ:DP:GL\t0/1:409:35:-20,-5,-20",
            "1\t3062915\t.\tG\tT\t1806\tq10\tDP=35;DP4=1,2,3,4;AN=2;AC=1\tGT:GQ:DP:GL\t0/1:409:35:-20,-5,-20",
            "2\t3199815\t.\tC\tT\t481\tPASS\tDP=26;AN=2;AC=1\tGT:GQ:DP\t1/2:322:26",
        ]
    );
}

#[test]
fn view_regions_file_filters_tab_regions() {
    let path = fixture_path("regions.vcf");
    let regions = fixture_path("regions.tab");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-R",
        regions.to_str().unwrap(),
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -R failed: {err}");
    let projected = out
        .lines()
        .filter(|line| !line.is_empty())
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{} {} {},{}", fields[0], fields[1], fields[3], fields[4])
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    let expected = std::fs::read_to_string(fixture_path("regions.out")).unwrap();
    assert_eq!(projected, expected);
}

#[test]
fn view_regions_support_braced_contig_names_with_colons() {
    let path = fixture_path("weird-chr-names.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-r",
        "{1:1}:1,{1:1}:2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -r braced colon contig failed: {err}");
    let records = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(records, ["1:1:1", "1:1:2"]);
}

#[test]
fn view_regions_support_braced_contig_names_with_intervals() {
    let path = fixture_path("weird-chr-names.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-r",
        "{1:1-1}:1-1",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -r braced interval contig failed: {err}");
    let records = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(records, ["1:1-1:1"]);
}

#[test]
fn view_targets_option_filters_contig_targets() {
    let path = fixture_path("view-t.vcf");
    let expected = std::fs::read_to_string(fixture_path("view-t.1.out")).unwrap();
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-t",
        "2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -t failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_targets_option_excludes_contig_targets() {
    let path = fixture_path("view-t.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-t",
        "^2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -t ^ failed: {err}");
    assert_eq!(
        out,
        "1\t1\t.\tA\tC\t.\t.\t.\n\
3\t2\t.\tA\tC\t.\t.\t.\n"
    );
}

#[test]
fn view_targets_file_filters_site_targets() {
    let path = fixture_path("view.sites.vcf");
    let targets = fixture_path("view.sites.txt");
    let expected = std::fs::read_to_string(fixture_path("view.sites.1.out")).unwrap();
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-T",
        targets.to_str().unwrap(),
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -T failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_targets_file_excludes_site_targets() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let path = fixture_path("view.sites.vcf");
    let targets = tmp.path().join("exclude-sites.txt");
    std::fs::write(&targets, "1\t10002\n").unwrap();
    let excluded = format!("^{}", targets.display());

    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-T",
        &excluded,
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -T ^file failed: {err}");
    assert_eq!(
        out,
        "1\t10001\t.\tG\tC\t40\t.\t.\n\
1\t10003\t.\tA\tG\t60\t.\t.\n\
1\t10004\t.\tA\tT\t70\t.\t.\n"
    );
}

#[test]
fn view_targets_overlap_modes_match_upstream_fixtures() {
    let path = fixture_path("overlap.vcf");
    for (mode, expected_name) in [
        ("0", "overlap.0.out"),
        ("1", "overlap.1.out"),
        ("2", "overlap.2.out"),
    ] {
        let expected = std::fs::read_to_string(fixture_path(expected_name)).unwrap();
        let (out, err, code) = run(&[
            "view",
            "--no-version",
            "-H",
            "-t",
            "chr1:100-200",
            "--targets-overlap",
            mode,
            path.to_str().unwrap(),
        ]);
        assert_eq!(code, 0, "view --targets-overlap {mode} failed: {err}");
        assert_eq!(out, expected, "unexpected overlap mode {mode} output");
    }
}

#[test]
fn view_regions_overlap_modes_match_upstream_fixtures() {
    let path = fixture_path("overlap.vcf");
    for (mode, expected_name) in [
        ("0", "overlap.0.out"),
        ("1", "overlap.1.out"),
        ("2", "overlap.2.out"),
    ] {
        let expected = std::fs::read_to_string(fixture_path(expected_name)).unwrap();
        let (out, err, code) = run(&[
            "view",
            "--no-version",
            "-H",
            "-r",
            "chr1:100-200",
            "--regions-overlap",
            mode,
            path.to_str().unwrap(),
        ]);
        assert_eq!(code, 0, "view --regions-overlap {mode} failed: {err}");
        assert_eq!(
            out, expected,
            "unexpected region overlap mode {mode} output"
        );
    }
}

#[test]
fn view_targets_overlap_exclusion_matches_upstream_fixture() {
    let path = fixture_path("overlap.vcf");
    let expected = std::fs::read_to_string(fixture_path("overlap.neg2.out")).unwrap();
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-t",
        "^chr1:100-200",
        "--targets-overlap",
        "2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view --targets-overlap exclusion failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_region_filters_bcf_input() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("aa.vcf");
    let bcf = tmp.path().join("aa.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let (out, err, code) = run(&["view", "--no-version", "-H", bcf.to_str().unwrap(), "20"]);
    assert_eq!(code, 0, "view BCF region failed: {err}");
    let records: Vec<_> = out.lines().filter(|line| !line.is_empty()).collect();
    assert_eq!(records.len(), 12);
    assert!(records.iter().all(|line| line.starts_with("20\t")));
}

#[test]
fn view_samples_list_subsets_vcf_columns() {
    let path = fixture_path("query.smpl.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-s", "11", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -s failed: {err}");
    assert_eq!(
        out,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t11\n\
chr1\t10000\t.\tA\tC\t.\t.\t.\tGT\t1/1\n"
    );
}

#[test]
fn view_samples_file_exclusion_subsets_bcf_input() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("query.smpl.vcf");
    let bcf = tmp.path().join("query.smpl.bcf");
    let samples = fixture_path("query.smpl.11.txt");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob failed: {err}");

    let excluded = format!("^{}", samples.display());
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-S",
        &excluded,
        bcf.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -S ^file BCF failed: {err}");
    let chrom = out
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("#CHROM line");
    assert_eq!(
        chrom,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t00"
    );
    assert!(out.contains("chr1\t10000\t.\tA\tC\t.\t.\t.\tGT\t0/0\n"));
}

#[test]
fn view_samples_list_subsets_bcf_output() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("query.smpl.vcf");
    let bcf = tmp.path().join("query.smpl.11.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "-Ob",
        "-s",
        "11",
        "-o",
        bcf.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -Ob -s failed: {err}");

    let (out, err, code) = run(&["view", "--no-version", bcf.to_str().unwrap()]);
    assert_eq!(code, 0, "view subset BCF failed: {err}");
    assert!(out.contains("##contig=<ID=chr1>\n"));
    assert!(out.contains("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"));
    assert!(out.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t11\n"));
    assert!(out.contains("chr1\t10000\t.\tA\tC\t.\t.\t.\tGT\t1/1\n"));
}

#[test]
fn view_drop_genotypes_matches_upstream_fixture() {
    let path = fixture_path("view.omitgenotypes.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.dropgenotypes.out")).unwrap();
    let (out, err, code) = run(&["view", "--no-version", "-G", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -G failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_exclude_snps_matches_upstream_fixture() {
    let path = fixture_path("view.vcf");
    let expected = std::fs::read_to_string(fixture_path("view.9.out")).unwrap();
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-G",
        "-V",
        "snps",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -G -V snps failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_min_alleles_filters_multiallelic_records() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-m",
        "3",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -m 3 failed: {err}");
    assert_eq!(
        out,
        "20\t271225\t.\tT\tTTTA,TA\t999\tStrandBias\tDP4=29281,42401,27887,29245;DP=272732;INDEL;IS=95,0.748031;MQ=47;PV4=0,1,0,1;QD=0.0948;AN=6;AC=2,2\tGT:DP:GQ:PL\t0/2:33:49:151,53,203,0,52,159\t0/1:51:99:255,0,213,255,255,255\t1/2:47:99:255,255,255,255,0,241\n"
    );
}

#[test]
fn view_min_max_alleles_filters_biallelic_records() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-m2",
        "-M2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -m2 -M2 failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:5464562",
            "20:76962",
            "20:126310",
            "20:138125",
            "20:138148",
            "20:304568",
            "20:326891",
            "X:2928329",
            "X:2933066",
            "X:2942109",
            "X:3048719",
            "Y:8657215",
            "Y:10011673",
        ]
    );
}

#[test]
fn view_max_ac_filters_nonmajor_allele_counts_from_genotypes() {
    let path = fixture_path("view.minmaxac.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-C5:nonmajor",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -C5:nonmajor failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(positions, ["20:1234567"]);
}

#[test]
fn view_min_ac_filters_nonmajor_allele_counts_from_genotypes() {
    let path = fixture_path("view.minmaxac.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-c6:nonmajor",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -c6:nonmajor failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(positions, ["20:1234568"]);
}

#[test]
fn view_min_af_filters_major_allele_frequency_from_genotypes() {
    let path = fixture_path("view.minmaxac.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-q0.3:major",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -q0.3:major failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(positions, ["20:1234567"]);
}

#[test]
fn view_known_filter_selects_records_with_ids() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-H", "-k", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -k failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "20:76962",
            "20:138125",
            "20:138148",
            "X:2928329",
            "X:2933066",
            "X:2942109",
            "Y:10011673",
        ]
    );
}

#[test]
fn view_apply_filters_selects_pass_records() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-f",
        "PASS",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -f PASS failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:2343543",
            "11:5464562",
            "20:76962",
            "20:138125",
            "20:138148",
            "20:304568",
            "20:326891",
            "X:2928329",
            "X:2933066",
            "X:2942109",
            "X:3048719",
            "Y:8657215",
        ]
    );
}

#[test]
fn view_apply_filters_selects_semicolon_delimited_filter_tags() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-f",
        "StrandBias",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -f StrandBias failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(positions, ["20:126310", "20:271225"]);
}

#[test]
fn view_include_expression_filters_core_and_info_fields() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-i",
        "QUAL==999 && (FS<20 || FS>=41.02) && ICF>-0.1 && HWE*2>1.2",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -i expression failed: {err}");
    let positions = out
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(positions, ["X:2942109", "X:3048719"]);
}

#[test]
fn view_exclude_expression_filters_info_fields() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-e",
        "FS<20",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -e expression failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:2343543",
            "11:5464562",
            "20:76962",
            "20:126310",
            "20:138125",
            "20:138148",
            "20:271225",
            "20:304568",
            "20:326891",
            "X:2942109",
            "X:3048719",
            "Y:10011673",
        ]
    );
}

#[test]
fn view_expression_filters_indexed_info_vectors() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-n",
        "-e",
        "INDEL=1 || PV4[0]<0.006",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -n -e indexed expression failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        ["11:2343543", "11:5464562", "X:3048719", "Y:8657215"]
    );
}

#[test]
fn view_novel_filter_selects_records_without_ids() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-H", "-n", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -n failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:2343543",
            "11:5464562",
            "20:126310",
            "20:271225",
            "20:304568",
            "20:326891",
            "X:3048719",
            "Y:8657215",
        ]
    );
}

#[test]
fn view_uncalled_filter_selects_records_without_called_gt() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-H", "-u", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -u failed: {err}");
    assert_eq!(
        out,
        "11\t5464562\t.\tC\tT\t999\tPASS\tDP=0\tGT:PL:DP:GQ\t./.:0,0,0:.:.\t./.:0,0,0:.:.\t./.:0,0,0:.:.\n"
    );
}

#[test]
fn view_exclude_uncalled_filter_drops_records_without_called_gt() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-H", "-U", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -U failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:2343543",
            "20:76962",
            "20:126310",
            "20:138125",
            "20:138148",
            "20:271225",
            "20:304568",
            "20:326891",
            "X:2928329",
            "X:2933066",
            "X:2942109",
            "X:3048719",
            "Y:8657215",
            "Y:10011673",
        ]
    );
}

#[test]
fn view_genotype_filter_selects_records_with_missing_gt() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-g",
        "miss",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -g miss failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        ["11:5464562", "20:326891", "Y:8657215", "Y:10011673"]
    );
}

#[test]
fn view_genotype_filter_excludes_records_with_het_gt() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&[
        "view",
        "--no-version",
        "-H",
        "-g",
        "^het",
        path.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view -g ^het failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:2343543",
            "11:5464562",
            "X:2942109",
            "Y:8657215",
            "Y:10011673",
        ]
    );
}

#[test]
fn view_phased_filter_selects_all_phased_records() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-H", "-p", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -p failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        ["20:304568", "X:3048719", "Y:8657215", "Y:10011673"]
    );
}

#[test]
fn view_exclude_phased_filter_drops_all_phased_records() {
    let path = fixture_path("view.vcf");
    let (out, err, code) = run(&["view", "--no-version", "-H", "-P", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -P failed: {err}");
    let positions = out
        .lines()
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            format!("{}:{}", fields[0], fields[1])
        })
        .collect::<Vec<_>>();
    assert_eq!(
        positions,
        [
            "11:2343543",
            "11:5464562",
            "20:76962",
            "20:126310",
            "20:138125",
            "20:138148",
            "20:271225",
            "20:326891",
            "X:2928329",
            "X:2933066",
            "X:2942109",
        ]
    );
}

#[test]
fn view_drop_genotypes_no_header_matches_upstream_fixture() {
    let path = fixture_path("view.omitgenotypes.vcf");
    let expected =
        std::fs::read_to_string(fixture_path("view.dropgenotypes.noheader.out")).unwrap();
    let (out, err, code) = run(&["view", "--no-version", "-HG", path.to_str().unwrap()]);
    assert_eq!(code, 0, "view -HG failed: {err}");
    assert_eq!(out, expected);
}

#[test]
fn view_threads_writes_bgzf_vcf_output() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("aa.vcf");
    let output = tmp.path().join("aa.vcf.gz");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "--threads",
        "2",
        "-Oz",
        "-o",
        output.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view --threads -Oz failed: {err}");

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(output).unwrap());
    let mut decoded = String::new();
    decoder.read_to_string(&mut decoded).unwrap();
    let records = decoded
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .count();
    assert_eq!(records, 21);
}

#[test]
fn view_threads_writes_bcf_output() {
    let tmp = tempfile::TempDir::new().expect("tempdir");
    let input = fixture_path("aa.vcf");
    let output = tmp.path().join("aa.bcf");

    let (_out, err, code) = run(&[
        "view",
        "--no-version",
        "--threads",
        "2",
        "-Ob",
        "-o",
        output.to_str().unwrap(),
        input.to_str().unwrap(),
    ]);
    assert_eq!(code, 0, "view --threads -Ob failed: {err}");

    let (out, err, code) = run(&["view", "--no-version", "-H", output.to_str().unwrap()]);
    assert_eq!(code, 0, "view threaded BCF failed: {err}");
    let records = out.lines().filter(|line| !line.is_empty()).count();
    assert_eq!(records, 21);
}

#[test]
fn view_threads_rejects_non_integer_argument() {
    let path = fixture_path("aa.vcf");
    let (_out, err, code) = run(&["view", "--threads", "abc", path.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("Could not parse argument: --threads abc"));
}

#[test]
fn view_unknown_output_type_errors() {
    let path = fixture_path("aa.vcf");
    let (_out, err, code) = run(&["view", "-Oq", path.to_str().unwrap()]);
    assert_ne!(code, 0);
    assert!(err.contains("not recognised"), "stderr: {err}");
}

#[test]
fn view_default_injects_bcftools_version_and_command_lines() {
    // Without --no-version, `bcftools view` must emit
    // ##bcftools_viewVersion=<v> and ##bcftools_viewCommand=<cmdline>; Date=...
    // header lines, mirroring upstream `bcf_hdr_append_version`.
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);

    let header_lines: Vec<_> = out.lines().filter(|l| l.starts_with("##")).collect();
    let version_line = header_lines
        .iter()
        .find(|l| l.starts_with("##bcftools_viewVersion="))
        .unwrap_or_else(|| panic!("missing version line in header:\n{out}"));
    let command_line = header_lines
        .iter()
        .find(|l| l.starts_with("##bcftools_viewCommand="))
        .unwrap_or_else(|| panic!("missing command line in header:\n{out}"));

    // The version line ends with the htslib-rs version we're built against.
    assert!(version_line.contains("+htslib-"), "got: {version_line}");
    // The command line names the subcommand and includes a `Date=` field.
    assert!(command_line.contains("view"), "got: {command_line}");
    assert!(command_line.contains("; Date="), "got: {command_line}");
}

#[test]
fn view_no_version_suppresses_injected_header_lines() {
    let path = fixture_path("aa.vcf");
    let (out, _err, code) = run(&["view", "--no-version", "-h", path.to_str().unwrap()]);
    assert_eq!(code, 0);
    assert!(
        !out.contains("##bcftools_viewVersion="),
        "version line leaked despite --no-version:\n{out}"
    );
    assert!(
        !out.contains("##bcftools_viewCommand="),
        "command line leaked despite --no-version:\n{out}"
    );
}
