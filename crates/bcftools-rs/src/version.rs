//! bcftools version constants and helpers.
//!
//! The pinned upstream bcftools tag is recorded here. The string emitted by
//! `bcftools --version-only` is `<bcftools-version>+htslib-<htslib-version>`,
//! reproduced by [`version_only_string`].

/// Pinned upstream bcftools version that this port tracks.
///
/// Sourced from `bcftools/version.sh` in the vendored submodule.
pub const BCFTOOLS_VERSION: &str = "1.23.1";

/// Version of the `htslib-rs` workspace this binary links against.
///
/// Reported by `bcftools --version` after the literal `htslib `.
pub const HTSLIB_RS_VERSION: &str = "0.1.0";

/// String reproduced by upstream's `bcftools --version-only`.
pub fn version_only_string() -> String {
    format!("{}+htslib-{}", BCFTOOLS_VERSION, HTSLIB_RS_VERSION)
}

/// Multi-line block reproduced by upstream's `bcftools version` / `--version`.
pub fn version_block() -> String {
    format!(
        "bcftools {bv}\nUsing htslib {hv}\nCopyright (C) 2025 Genome Research Ltd.\n\
         License Expat: The MIT/Expat license\n\
         This is free software: you are free to change and redistribute it.\n\
         There is NO WARRANTY, to the extent permitted by law.\n",
        bv = BCFTOOLS_VERSION,
        hv = HTSLIB_RS_VERSION,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_only_format() {
        let s = version_only_string();
        assert!(s.starts_with(BCFTOOLS_VERSION));
        assert!(s.contains("+htslib-"));
        assert!(s.ends_with(HTSLIB_RS_VERSION));
    }

    #[test]
    fn version_block_contains_both_versions() {
        let s = version_block();
        assert!(s.contains(BCFTOOLS_VERSION));
        assert!(s.contains(HTSLIB_RS_VERSION));
        assert!(s.contains("Copyright"));
    }
}
