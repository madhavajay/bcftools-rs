use std::path::PathBuf;

fn fixture_path(name: &str) -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("..");
    p.push("..");
    p.push("bcftools");
    p.push("test");
    p.push(name);
    p
}

#[test]
fn rust_gate_can_locate_top_level_and_nested_upstream_fixtures() {
    for name in [
        "mpileup.2.vcf",
        "csq/sort-csq",
        "mpileup/mpileup.ref.fa",
        "mpileup/mpileup.1.bam",
    ] {
        let path = fixture_path(name);
        assert!(
            path.exists(),
            "missing upstream fixture {} at {}",
            name,
            path.display()
        );
    }
}
