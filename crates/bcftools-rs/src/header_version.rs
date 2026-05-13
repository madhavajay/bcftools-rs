//! Port of `bcf_hdr_append_version` from upstream `vcfmerge.c:3362`.
//!
//! Each subcommand emits two header lines into its output VCF/BCF when the
//! caller did *not* pass `--no-version`:
//!
//! ```text
//! ##<cmd>Version=<bcftools-version>+htslib-<htslib-version>
//! ##<cmd>Command=<argv0> <argv1> ...; Date=<ctime>
//! ```
//!
//! Argv tokens that contain a space are wrapped in single quotes, exactly like
//! upstream. The date is rendered in the same `ctime`-shape format
//! (`Mon Jan  1 00:00:00 2026`).
//!
//! See `bcftools/vcfmerge.c:3362-3398` for the C implementation.

use std::ffi::OsStr;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::version::{BCFTOOLS_VERSION, HTSLIB_RS_VERSION};

/// One pair of header lines to be appended to the output VCF/BCF header.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HeaderVersionLines {
    /// `##<cmd>Version=...` line, no trailing newline.
    pub version_line: String,
    /// `##<cmd>Command=...; Date=...` line, no trailing newline.
    pub command_line: String,
}

/// Build the two `##<cmd>Version` / `##<cmd>Command` header lines.
///
/// `cmd` is the subcommand identifier upstream uses, e.g. `"bcftools_view"` or
/// `"bcftools/csq"` — it is interpolated verbatim after the `##` prefix.
///
/// `argv` is the full argv reconstruction starting with the program name. Each
/// element containing a space is single-quoted in the output.
///
/// `now_unix` is the timestamp the `Date=` field is rendered from. Pass
/// `SystemTime::now()` for normal runs; pass a fixed value in tests to keep
/// output deterministic.
pub fn build_lines<S>(cmd: &str, argv: &[S], now_unix: SystemTime) -> HeaderVersionLines
where
    S: AsRef<OsStr>,
{
    let version_line = format!(
        "##{cmd}Version={bv}+htslib-{hv}",
        cmd = cmd,
        bv = BCFTOOLS_VERSION,
        hv = HTSLIB_RS_VERSION,
    );

    let mut command_line = format!("##{cmd}Command=");
    for (i, arg) in argv.iter().enumerate() {
        if i > 0 {
            command_line.push(' ');
        }
        let s = arg.as_ref().to_string_lossy();
        if s.contains(' ') {
            command_line.push('\'');
            command_line.push_str(&s);
            command_line.push('\'');
        } else {
            command_line.push_str(&s);
        }
    }
    command_line.push_str("; Date=");
    command_line.push_str(&format_ctime(now_unix));

    HeaderVersionLines {
        version_line,
        command_line,
    }
}

/// Format a `SystemTime` as `"Mon Jan  1 00:00:00 2026"`, the same shape that
/// `ctime(3)` produces (without the trailing newline).
///
/// Uses local time. The width of the day-of-month field is 2 (space-padded)
/// to match `%e` semantics, matching upstream's `ctime()` output exactly.
pub fn format_ctime(t: SystemTime) -> String {
    let secs = t
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (y, mo, d, h, mi, s, wd) = civil_from_unix(secs);
    const DOW: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    const MON: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    format!(
        "{dow} {mon} {d:>2} {h:02}:{mi:02}:{s:02} {y}",
        dow = DOW[wd as usize],
        mon = MON[(mo - 1) as usize],
        d = d,
        h = h,
        mi = mi,
        s = s,
        y = y,
    )
}

/// Convert a Unix timestamp to civil (year, month, day, hour, minute, second,
/// weekday) in UTC. Weekday: 0 = Sun, 6 = Sat.
///
/// Algorithm from Howard Hinnant's `date` library — exact across all of the
/// proleptic Gregorian range we will plausibly hit.
fn civil_from_unix(secs: i64) -> (i32, u32, u32, u32, u32, u32, u32) {
    let days = secs.div_euclid(86_400);
    let tod = secs.rem_euclid(86_400) as u32;
    let h = tod / 3600;
    let mi = (tod / 60) % 60;
    let s = tod % 60;

    let wd = ((days + 4).rem_euclid(7)) as u32;

    let z = days + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u32;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i32 + era as i32 * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };

    (y, m, d, h, mi, s, wd)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::OsString;
    use std::time::Duration;

    fn fixed_time() -> SystemTime {
        // 2024-01-15T12:34:56Z = 1705322096
        UNIX_EPOCH + Duration::from_secs(1_705_322_096)
    }

    #[test]
    fn version_line_format() {
        let argv: Vec<OsString> = vec!["bcftools".into(), "view".into(), "-Oz".into()];
        let lines = build_lines("bcftools_view", &argv, fixed_time());
        assert_eq!(
            lines.version_line,
            format!(
                "##bcftools_viewVersion={}+htslib-{}",
                BCFTOOLS_VERSION, HTSLIB_RS_VERSION
            )
        );
    }

    #[test]
    fn command_line_quotes_args_with_spaces() {
        let argv: Vec<OsString> = vec![
            "bcftools".into(),
            "view".into(),
            "-i".into(),
            "QUAL > 30".into(),
            "in.vcf".into(),
        ];
        let lines = build_lines("bcftools_view", &argv, fixed_time());
        assert!(
            lines
                .command_line
                .starts_with("##bcftools_viewCommand=bcftools view -i 'QUAL > 30' in.vcf; Date=")
        );
    }

    #[test]
    fn ctime_format_known_timestamp() {
        // 2024-01-15T12:34:56Z is a Monday.
        assert_eq!(format_ctime(fixed_time()), "Mon Jan 15 12:34:56 2024");
    }

    #[test]
    fn ctime_pads_single_digit_day() {
        // 2024-01-01T00:00:00Z is a Monday.
        let t = UNIX_EPOCH + Duration::from_secs(1_704_067_200);
        assert_eq!(format_ctime(t), "Mon Jan  1 00:00:00 2024");
    }

    #[test]
    fn csq_uses_slash_separator() {
        let argv: Vec<OsString> = vec!["bcftools".into(), "csq".into()];
        let lines = build_lines("bcftools/csq", &argv, fixed_time());
        assert!(lines.version_line.starts_with("##bcftools/csqVersion="));
        assert!(lines.command_line.starts_with("##bcftools/csqCommand="));
    }
}
