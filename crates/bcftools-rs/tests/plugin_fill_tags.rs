//! End-to-end parity for the `+fill-tags` first slice: the
//! genotype-derived count tags `AN/AC/AC_Hom/AC_Het/AC_Hemi/AF/MAF/NS`,
//! `-t` selection, and `-S` population grouping (upstream test.pl rows
//! 695-697, 699). The harness compares full output (no `grep`).

use std::path::PathBuf;
use std::process::Command;

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

fn check(input: &str, plugin_args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec!["+fill-tags", "--no-version", input.to_str().unwrap(), "--"];
    full.extend_from_slice(plugin_args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn +fill-tags");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

/// Like [`check`] but drops all `#`-prefixed lines, mirroring the
/// upstream `| grep -v ^#` harness filter (the `fmissing.*` rows).
fn check_nohdr(input: &str, plugin_args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+fill-tags", "--no-version", input.to_str().unwrap(), "--"];
    full.extend_from_slice(plugin_args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn +fill-tags");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8(out.stdout).unwrap();
    let filtered: String = stdout
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn count_tags_merge_a() {
    check(
        "merge.a.vcf",
        &["-t", "AN,AC,AC_Hom,AC_Het,AC_Hemi"],
        "fill-tags.out",
    );
}

#[test]
fn af_maf_ns_view() {
    check("view.vcf", &["-t", "AC,AN,AF,MAF,NS"], "fill-tags.2.out");
}

#[test]
fn per_population_grouping() {
    let smpl = fixture_path("fill-tags.3.smpl");
    check(
        "view.vcf",
        &["-t", "AC", "-S", smpl.to_str().unwrap()],
        "fill-tags.3.out",
    );
}

#[test]
fn many_alt_alleles() {
    check("many-alts.vcf", &["-t", "AN,AC"], "fill-tags.4.out");
}

#[test]
fn default_all_tag_set_with_hwe_exchet() {
    // No `-t` ⇒ the `all` set: F_MISSING/NS/AN/AF/MAF/AC/AC_Het/AC_Hom/
    // AC_Hemi/HWE/ExcHet (+ the VAF/VAF1 ##FORMAT header lines).
    check("fill-tags-hemi.vcf", &[], "fill-tags-hemi.1.out");
    check("fill-tags-hwe.vcf", &[], "fill-tags-hwe.out");
}

#[test]
fn drop_missing_flag() {
    check("fill-tags-hemi.vcf", &["-d"], "fill-tags-hemi.2.out");
}

#[test]
fn format_vaf_vaf1() {
    // FORMAT/VAF + VAF1 from FORMAT/AD, no GT in the record.
    check(
        "fill-tags-VAF.vcf",
        &["-t", "VAF,VAF1"],
        "fill-tags-VAF.out",
    );
}

#[test]
fn func_sum_smpl_sum() {
    // The TAG:Num=EXPR engine: int(sum(...)) / int(smpl_sum(...)).
    check(
        "fill-tags-AD.vcf",
        &["-t", "INFO/DP:1=int(sum(FMT/AD))"],
        "fill-tags-AD.1.out",
    );
    check(
        "fill-tags-AD.vcf",
        &["-t", "INFO/DP:1=int(sum(INFO/AD))"],
        "fill-tags-AD.2.out",
    );
    check(
        "fill-tags-AD.vcf",
        &["-t", "FORMAT/DP:1=int(smpl_sum(FMT/AD))"],
        "fill-tags-AD.3.out",
    );
}

#[test]
fn func_with_population_grouping() {
    let smpl = fixture_path("fill-tags.3.smpl");
    check(
        "view.vcf",
        &[
            "-t",
            "DP:1=int(sum(FORMAT/DP))",
            "-S",
            smpl.to_str().unwrap(),
        ],
        "fill-tags.5.out",
    );
}

#[test]
fn func_f_pass() {
    // F_PASS(EXPR): fraction of samples where the per-sample filter
    // expression holds (full output, no `grep`).
    check(
        "fill-tags-hwe.vcf",
        &["-t", r#"XX:1=F_PASS(GT="alt")"#],
        "fill-tags-func.out",
    );
}

#[test]
fn func_f_pass_n_pass_missing() {
    let smpl = fixture_path("fmissing.txt");
    let s = smpl.to_str().unwrap();
    // `-t F_MISSING` (builtin), the `F_PASS(GT="mis")` func form, and
    // `int(N_PASS(GT="mis"))` — the harness compares after `grep -v ^#`.
    check_nohdr(
        "fmissing.vcf",
        &["-S", s, "-t", "F_MISSING"],
        "fmissing.1.out",
    );
    check_nohdr(
        "fmissing.vcf",
        &["-S", s, "-t", r#"F_MISSING:1=F_PASS(GT="mis")"#],
        "fmissing.1.out",
    );
    check_nohdr(
        "fmissing.vcf",
        &["-S", s, "-t", r#"N_MISSING:1=int(N_PASS(GT="mis"))"#],
        "fmissing.2.out",
    );
}

#[test]
fn func_n_pass_subscript_and_binom() {
    // Comma-joined func list with subscripts and `binom(...)` (whose
    // own comma must not split the `-t` list).
    check(
        "fill-tags-AD.vcf",
        &["-t", "XX=N_PASS(FMT/AD[:0]<=10),YY=N_PASS(FMT/AD[:0]>10)"],
        "fill-tags-AD.4.out",
    );
    check(
        "fill-tags-AD.vcf",
        &[
            "-t",
            "good=N_PASS(binom(FMT/AD[:0],FMT/AD[:1])>=1e-5),bad=N_PASS(binom(FMT/AD[:0],FMT/AD[:1])<1e-5)",
        ],
        "fill-tags-AD.5.out",
    );
}

#[test]
fn end_type_with_all_set() {
    // `-t all,END,TYPE,F_MISSING` — same output for an input that
    // pre-declares every tag (`fill-tags-rw`) and one that declares
    // only NS/AN/AF/AC (`fill-tags-AN0`).
    check(
        "fill-tags-rw.vcf",
        &["-t", "all,END,TYPE,F_MISSING"],
        "fill-tags-AN0.out",
    );
    check(
        "fill-tags-AN0.vcf",
        &["-t", "all,END,TYPE,F_MISSING"],
        "fill-tags-AN0.out",
    );
}

#[test]
fn func_phred_fisher() {
    // The general-expression path: phred(fisher(...)) over a 4-value
    // INFO/DP4, INFO arrays, and a per-sample FMT/DP4.
    check(
        "fisher.vcf",
        &["-t", "FT:1=phred(fisher(INFO/DP4))"],
        "fisher.1.out",
    );
    check(
        "fisher.vcf",
        &["-t", "FT:1=phred(fisher(INFO/ADF,INFO/ADR))"],
        "fisher.3.out",
    );
    check(
        "fisher.vcf",
        &["-t", "FMT/FT:1=phred(fisher(FMT/DP4))"],
        "fisher.4.out",
    );
    check(
        "fisher.vcf",
        &["-t", "FT:1=phred(fisher(INFO/ADF[0,2],INFO/ADR[0,2]))"],
        "fisher.2.out",
    );
    check(
        "fisher.vcf",
        &["-t", "FMT/FT:1=phred(fisher(FORMAT/ADF,FORMAT/ADR))"],
        "fisher.5.out",
    );
    check(
        "fisher.vcf",
        &[
            "-t",
            "FMT/FT:1=phred(fisher(FORMAT/ADF[:0,1],FORMAT/ADR[:0,1]))",
        ],
        "fisher.6.out",
    );
}
