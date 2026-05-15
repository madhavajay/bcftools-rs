//! End-to-end test for the `+variantkey-hex` plugin against the upstream
//! `variantkey-hex.out` summary fixture.

use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

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

#[test]
fn variantkey_hex_summary_and_lookup_files() {
    ensure_binary_built();
    let input = fixture_path("query.variantkey.vcf");
    let expected = std::fs::read_to_string(fixture_path("variantkey-hex.out")).unwrap();

    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    let dir = std::env::temp_dir().join(format!("vkhex-it-{}-{nanos}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let dir_arg = format!("{}/", dir.display());

    let out = Command::new(bin_path())
        .args(["+variantkey-hex", input.to_str().unwrap(), dir_arg.as_str()])
        .output()
        .expect("spawn bcftools");
    assert_eq!(
        out.status.code(),
        Some(0),
        "stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let stdout = String::from_utf8(out.stdout).unwrap();
    assert_eq!(stdout, expected);

    // 66 input records -> 66 vkrs/rsvk rows; the 3 hash-mode keys -> 3 nrvk.
    let vkrs = std::fs::read_to_string(dir.join("vkrs.unsorted.hex")).unwrap();
    let rsvk = std::fs::read_to_string(dir.join("rsvk.unsorted.hex")).unwrap();
    let nrvk = std::fs::read_to_string(dir.join("nrvk.unsorted.tsv")).unwrap();
    assert_eq!(vkrs.lines().count(), 66);
    assert_eq!(rsvk.lines().count(), 66);
    assert_eq!(nrvk.lines().count(), 3);

    let _ = std::fs::remove_dir_all(&dir);
}
