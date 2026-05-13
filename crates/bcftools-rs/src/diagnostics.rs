//! Error reporting helpers matching upstream `bcftools` stderr conventions.
//!
//! Upstream `version.c` defines two functions:
//!
//! - `error(fmt, ...)` — prints to stderr and `exit(-1)`. Caller controls the
//!   trailing newline.
//! - `error_errno(fmt, ...)` — prints to stderr, then appends `: strerror(errno)`
//!   if `errno != 0`, then a newline, then `exit(-1)`. Caller MUST NOT append a
//!   newline.
//!
//! The Rust analogues here panic-by-process-exit the same way. The
//! `[E::funcname] message` shape used pervasively in upstream subcommands is
//! produced by callers; these helpers do not prepend anything.

use std::io::{self, Write};
use std::process::exit;

/// Write `msg` to stderr verbatim and exit with status 255 (mirroring
/// upstream's `exit(-1)` truncation to `255`).
pub fn error(msg: impl AsRef<str>) -> ! {
    let mut stderr = io::stderr().lock();
    let _ = stderr.write_all(msg.as_ref().as_bytes());
    let _ = stderr.flush();
    exit(255)
}

/// Write `msg` to stderr, then `: <reason>\n` if `reason` is `Some`, else `\n`,
/// then exit with status 255.
///
/// Matches upstream `error_errno`: the caller's message must NOT end in a
/// newline; this helper always appends one.
pub fn error_errno(msg: impl AsRef<str>, reason: Option<&io::Error>) -> ! {
    let mut stderr = io::stderr().lock();
    let _ = stderr.write_all(msg.as_ref().as_bytes());
    if let Some(e) = reason {
        let _ = writeln!(stderr, ": {e}");
    } else {
        let _ = stderr.write_all(b"\n");
    }
    let _ = stderr.flush();
    exit(255)
}

/// Format an upstream-style `[E::funcname] msg` error tag.
pub fn fmt_etag(funcname: &str, msg: &str) -> String {
    format!("[E::{funcname}] {msg}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn etag_shape() {
        assert_eq!(fmt_etag("main_vcfhead", "boom"), "[E::main_vcfhead] boom");
    }
}
