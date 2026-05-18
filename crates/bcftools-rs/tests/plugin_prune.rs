//! End-to-end parity tests for `+prune` window mode against the upstream
//! `prune.1.{4,6}.out` fixtures (harness strips `^##bcftools_`).

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

fn check(args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("prune.1.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();

    let mut full = vec!["+prune", input.to_str().unwrap()];
    full.extend_from_slice(args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
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

#[test]
fn prune_maxaf_af_tag() {
    check(&["-w", "2bp", "-n", "1", "--AF-tag", "AF"], "prune.1.4.out");
}

#[test]
fn prune_first_mode() {
    check(&["-w", "2bp", "-n", "1", "-N", "1st"], "prune.1.6.out");
}

#[test]
fn prune_include_filter_discards_ref_only_sites() {
    // Upstream row prune.1.5: same as prune.1.4 but `-i 'GT="alt"'`
    // discards REF-only sites before the windowed pruning.
    check(
        &["-w", "2bp", "-n", "1", "--AF-tag", "AF", "-i", "GT=\"alt\""],
        "prune.1.5.out",
    );
}

fn check3(args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("prune.3.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+prune", input.to_str().unwrap()];
    full.extend_from_slice(args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "{full:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let filtered: String = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with("##bcftools_"))
        .map(|l| format!("{l}\n"))
        .collect();
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn prune_cluster_count_modes() {
    // Upstream rows in=>'prune.3': `-a count` (CLUSTER_SIZE annot, with
    // and without `-k`) and `-m count=2` (drop clusters > 2 sites).
    check3(
        &["-w", "3bp", "-a", "count", "-i", "XX!=0", "-k"],
        "prune.3.1.out",
    );
    check3(
        &["-w", "3bp", "-a", "count", "-i", "XX!=0"],
        "prune.3.2.out",
    );
    check3(
        &["-w", "3bp", "-m", "count=2", "-i", "XX!=0"],
        "prune.3.3.out",
    );
}

fn check_in(input: &str, args: &[&str], expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path(input);
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let mut full = vec!["+prune", input.to_str().unwrap()];
    full.extend_from_slice(args);
    let out = Command::new(bin_path())
        .args(&full)
        .output()
        .expect("spawn bcftools");
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

#[test]
fn prune_ld_annotate_r2_ld_hd() {
    check_in(
        "prune.1.vcf",
        &["-w", "1", "-a", "r2,LD,HD"],
        "prune.1.1.out",
    );
}

#[test]
fn prune_ld_max_soft_filter() {
    check_in(
        "prune.1.vcf",
        &["-w", "2", "-a", "r2", "-m", "0.5", "-f", "MaxR2"],
        "prune.1.2.out",
    );
}

#[test]
fn prune_ld_max_hard_filter() {
    check_in(
        "prune.1.vcf",
        &["-w", "2", "-a", "r2", "-m", "0.5"],
        "prune.1.3.out",
    );
}

#[test]
fn prune_ld_annotate_multisample() {
    check_in(
        "prune.2.vcf",
        &["-w", "1", "-a", "r2,LD,HD"],
        "prune.2.1.out",
    );
}
