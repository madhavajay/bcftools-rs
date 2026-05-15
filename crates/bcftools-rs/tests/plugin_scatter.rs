//! End-to-end parity tests for `+scatter` against the upstream
//! `scatter.1.{1,2,3}.out` fixtures (test.pl rows 887-889; the harness
//! sorts the output dir, `cat`s the files and pipes `grep -v ^##`).

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

fn check(args: &[&str], tmp_name: &str, expected_fixture: &str) {
    ensure_binary_built();
    let input = fixture_path("scatter.1.vcf");
    let expected = std::fs::read_to_string(fixture_path(expected_fixture)).unwrap();
    let dir = std::env::temp_dir().join(format!(
        "bcftools-rs-scatter-test-{}-{tmp_name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);

    let mut full = vec![
        "+scatter".to_string(),
        input.to_str().unwrap().to_string(),
        "-o".to_string(),
        dir.to_str().unwrap().to_string(),
    ];
    full.extend(args.iter().map(|s| s.to_string()));
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

    // Sort the dir entries, cat them, strip `##` lines.
    let mut files: Vec<String> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| !n.starts_with('.'))
        .collect();
    files.sort();
    let mut cat = String::new();
    for f in &files {
        let c = std::fs::read_to_string(dir.join(f)).unwrap();
        cat.push_str(&c);
    }
    let filtered: String = cat
        .lines()
        .filter(|l| !l.starts_with("##"))
        .map(|l| format!("{l}\n"))
        .collect();
    let _ = std::fs::remove_dir_all(&dir);
    assert_eq!(filtered, expected, "mismatch for {full:?}");
}

#[test]
fn scatter_nsites_3() {
    check(&["-n", "3"], "n3", "scatter.1.1.out");
}

#[test]
fn scatter_regions() {
    check(&["-s", "21,22"], "s", "scatter.1.2.out");
}

#[test]
fn scatter_regions_extra() {
    check(&["-s", "21,22", "-x", "X"], "sx", "scatter.1.3.out");
}
