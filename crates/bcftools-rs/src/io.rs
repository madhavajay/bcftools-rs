//! HTSlib write-mode and verbosity helpers.
//!
//! Ports of `version.c::hts_bcf_wmode`, `hts_bcf_wmode2`, `set_wmode`,
//! `apply_verbosity`, `parse_overlap_option`, and `write_index_parse`.
//!
//! These are the mode-string conventions the rest of bcftools speaks:
//!
//! | Mode | Meaning                  |
//! | ---- | ------------------------ |
//! | `w`  | uncompressed VCF text    |
//! | `wz` | bgzf-compressed VCF      |
//! | `wbu`| uncompressed BCF         |
//! | `wb` | compressed BCF (default) |

/// HTSlib file-type bits, mirroring `bcftools.h:44-50`.
pub mod file_type {
    /// Plain TSV/text file (custom non-VCF formats).
    pub const TAB_TEXT: i32 = 0;
    /// Bgzf-compressed flag.
    pub const GZ: i32 = 1;
    /// VCF text format.
    pub const VCF: i32 = 2;
    /// Bgzf-compressed VCF.
    pub const VCF_GZ: i32 = GZ | VCF;
    /// Binary BCF (uncompressed sentinel; combine with `GZ` for compressed BCF).
    pub const BCF: i32 = 1 << 2;
    /// Bgzf-compressed BCF.
    pub const BCF_GZ: i32 = GZ | BCF;
    /// Stdin/stdout sentinel.
    pub const STDIN: i32 = 1 << 3;
}

/// HTSlib's `##idx##` separator between a data path and a paired index path.
pub const HTS_IDX_DELIM: &str = "##idx##";

/// HTSlib `HTS_FMT_TBI` constant (mirrors upstream).
pub const HTS_FMT_TBI: i32 = 2;
/// HTSlib `HTS_FMT_CSI` constant (mirrors upstream; zero by design).
pub const HTS_FMT_CSI: i32 = 0;

/// Bit set on the parsed `--write-index` value to signal "enabled".
///
/// HTSlib's `HTS_FMT_CSI` is 0, so callers cannot use the format value alone
/// as a boolean flag — upstream layers an extra `128` bit on top.
pub const WRITE_INDEX_ENABLED_BIT: i32 = 128;

/// Port of `hts_bcf_wmode`.
pub fn hts_bcf_wmode(file_type: i32) -> &'static str {
    if file_type == file_type::BCF {
        "wbu"
    } else if file_type & file_type::BCF != 0 {
        "wb"
    } else if file_type & file_type::GZ != 0 {
        "wz"
    } else {
        "w"
    }
}

/// Strip the optional `##idx##<index-path>` suffix from an output path,
/// returning just the data path.
fn strip_idx_suffix(fname: &str) -> &str {
    match fname.find(HTS_IDX_DELIM) {
        Some(i) => &fname[..i],
        None => fname,
    }
}

fn ends_with_ci(s: &str, suffix: &str) -> bool {
    let n = suffix.len();
    s.len() >= n && s[s.len() - n..].eq_ignore_ascii_case(suffix)
}

/// Port of `hts_bcf_wmode2`.
pub fn hts_bcf_wmode2(file_type: i32, fname: Option<&str>) -> &'static str {
    let Some(fname) = fname else {
        return hts_bcf_wmode(file_type);
    };
    let head = strip_idx_suffix(fname);
    if ends_with_ci(head, ".bcf") {
        return hts_bcf_wmode(file_type::BCF | file_type::GZ);
    }
    if ends_with_ci(head, ".vcf") {
        return hts_bcf_wmode(file_type::VCF);
    }
    if ends_with_ci(head, ".vcf.gz") || ends_with_ci(head, ".vcf.bgz") {
        return hts_bcf_wmode(file_type::VCF | file_type::GZ);
    }
    hts_bcf_wmode(file_type)
}

/// Port of `set_wmode`.
///
/// Returns the assembled mode string (e.g. `"wb6"`). `clevel` of `-1` means
/// "do not append a compression level". Returns `Err` if a level is supplied
/// for an uncompressed-only mode (matches upstream's `error()` exit there).
pub fn set_wmode(file_type: i32, fname: Option<&str>, clevel: i32) -> Result<String, String> {
    let head = fname.map(strip_idx_suffix);
    let base = match head {
        Some(h) if ends_with_ci(h, ".bcf") => hts_bcf_wmode(if file_type & file_type::BCF != 0 {
            file_type
        } else {
            file_type::BCF | file_type::GZ
        }),
        Some(h) if ends_with_ci(h, ".vcf") => hts_bcf_wmode(file_type::VCF),
        Some(h) if ends_with_ci(h, ".vcf.gz") || ends_with_ci(h, ".vcf.bgz") => {
            hts_bcf_wmode(file_type::VCF | file_type::GZ)
        }
        _ => hts_bcf_wmode(file_type),
    };
    if (0..=9).contains(&clevel) {
        if base.contains('v') || base.contains('u') {
            return Err(format!(
                "Error: compression level ({}) cannot be set on uncompressed streams ({})\n",
                clevel,
                fname.unwrap_or("-"),
            ));
        }
        Ok(format!("{base}{clevel}"))
    } else {
        Ok(base.to_string())
    }
}

/// Error returned by [`apply_verbosity`] when the input is not a non-negative
/// integer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InvalidVerbosity;

impl std::fmt::Display for InvalidVerbosity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("invalid verbosity (expected non-negative integer)")
    }
}

impl std::error::Error for InvalidVerbosity {}

/// Port of `apply_verbosity`.
///
/// Returns `Ok(level)` on a valid non-negative integer string. Levels above 3
/// are treated as a tracing dial in upstream; we round-trip the value so a
/// future logging bridge can apply it.
pub fn apply_verbosity(s: &str) -> Result<u32, InvalidVerbosity> {
    s.parse::<u32>().map_err(|_| InvalidVerbosity)
}

/// Port of `parse_overlap_option`.
pub fn parse_overlap_option(arg: &str) -> Option<u8> {
    match arg {
        a if a.eq_ignore_ascii_case("pos") || a == "0" => Some(0),
        a if a.eq_ignore_ascii_case("record") || a == "1" => Some(1),
        a if a.eq_ignore_ascii_case("variant") || a == "2" => Some(2),
        _ => None,
    }
}

/// Port of `write_index_parse`. `arg = None` matches upstream calling with a
/// null pointer (default to CSI).
///
/// Returns `Some(value)` to set `args->write_index` to, or `None` for an
/// unparseable argument (upstream returns 0).
pub fn write_index_parse(arg: Option<&str>) -> Option<i32> {
    let fmt = match arg {
        None => HTS_FMT_CSI,
        Some(a) if a == "csi" || a == "=csi" => HTS_FMT_CSI,
        Some(a) if a == "tbi" || a == "=tbi" => HTS_FMT_TBI,
        Some(_) => return None,
    };
    Some(WRITE_INDEX_ENABLED_BIT | fmt)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wmode_strings() {
        assert_eq!(hts_bcf_wmode(file_type::BCF), "wbu");
        assert_eq!(hts_bcf_wmode(file_type::BCF_GZ), "wb");
        assert_eq!(hts_bcf_wmode(file_type::VCF_GZ), "wz");
        assert_eq!(hts_bcf_wmode(file_type::VCF), "w");
        assert_eq!(hts_bcf_wmode(0), "w");
    }

    #[test]
    fn wmode2_extension_dispatch() {
        assert_eq!(hts_bcf_wmode2(0, Some("foo.bcf")), "wb");
        assert_eq!(hts_bcf_wmode2(0, Some("FOO.BCF")), "wb");
        assert_eq!(hts_bcf_wmode2(0, Some("foo.vcf")), "w");
        assert_eq!(hts_bcf_wmode2(0, Some("foo.vcf.gz")), "wz");
        assert_eq!(hts_bcf_wmode2(0, Some("foo.vcf.bgz")), "wz");
        assert_eq!(
            hts_bcf_wmode2(file_type::VCF, Some("foo.unknown")),
            "w",
            "unknown extension falls back to file_type"
        );
        assert_eq!(
            hts_bcf_wmode2(file_type::BCF_GZ, Some("foo.bcf##idx##foo.bcf.csi")),
            "wb",
            "##idx## suffix is stripped before dispatch"
        );
        assert_eq!(hts_bcf_wmode2(file_type::BCF, None), "wbu");
    }

    #[test]
    fn set_wmode_appends_clevel() {
        assert_eq!(set_wmode(0, Some("foo.bcf"), 6).unwrap(), "wb6");
        assert_eq!(set_wmode(0, Some("foo.vcf.gz"), 0).unwrap(), "wz0");
        assert_eq!(set_wmode(0, Some("foo.bcf"), -1).unwrap(), "wb");
    }

    #[test]
    fn set_wmode_rejects_clevel_only_when_mode_string_marks_uncompressed() {
        // Quirk of upstream `set_wmode`: it tests `strchr(ret,'v') || strchr(ret,'u')`,
        // but `hts_bcf_wmode` never includes 'v' in its returned strings, so the
        // only mode that triggers the error is uncompressed BCF ("wbu"). For
        // VCF (mode "w"), the C version silently produces "w6". We match that.
        assert_eq!(set_wmode(0, Some("foo.vcf"), 6).unwrap(), "w6");
        // BCF|GZ + `.bcf` extension → mode "wb" → clevel appended cleanly.
        assert_eq!(
            set_wmode(file_type::BCF_GZ, Some("foo.bcf"), 6).unwrap(),
            "wb6"
        );
        // The error path: input file_type is plain BCF (not BCF|GZ) and the
        // `.bcf` extension dispatch keeps it as `file_type` (BCF without GZ)
        // → mode "wbu" → contains 'u' → clevel rejected.
        assert!(set_wmode(file_type::BCF, Some("foo.bcf"), 6).is_err());
        // Path with no recognized extension and plain BCF → mode "wbu" → also
        // contains 'u' → clevel rejected.
        assert!(set_wmode(file_type::BCF, Some("foo.unknown"), 6).is_err());
    }

    #[test]
    fn overlap_option_parses() {
        assert_eq!(parse_overlap_option("pos"), Some(0));
        assert_eq!(parse_overlap_option("0"), Some(0));
        assert_eq!(parse_overlap_option("RECORD"), Some(1));
        assert_eq!(parse_overlap_option("variant"), Some(2));
        assert_eq!(parse_overlap_option("2"), Some(2));
        assert_eq!(parse_overlap_option("nope"), None);
    }

    #[test]
    fn write_index_parse_modes() {
        assert_eq!(
            write_index_parse(None),
            Some(WRITE_INDEX_ENABLED_BIT | HTS_FMT_CSI),
        );
        assert_eq!(
            write_index_parse(Some("csi")),
            Some(WRITE_INDEX_ENABLED_BIT | HTS_FMT_CSI),
        );
        assert_eq!(
            write_index_parse(Some("=tbi")),
            Some(WRITE_INDEX_ENABLED_BIT | HTS_FMT_TBI),
        );
        assert_eq!(write_index_parse(Some("garbage")), None);
    }

    #[test]
    fn apply_verbosity_accepts_nonneg_ints() {
        assert_eq!(apply_verbosity("0").unwrap(), 0);
        assert_eq!(apply_verbosity("4").unwrap(), 4);
        assert!(apply_verbosity("-1").is_err());
        assert!(apply_verbosity("abc").is_err());
    }
}
