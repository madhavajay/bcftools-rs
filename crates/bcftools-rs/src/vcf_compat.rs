//! VCF text-stream compatibility helpers used by subcommand readers.
//!
//! Upstream HTSlib accepts certain non-canonical VCF headers that strict
//! parsers reject. The most visible case is Kestrel's Java VCF writer, which
//! emits `##fileformat=VCF4.2` instead of the canonical `##fileformat=VCFv4.2`.
//! HTSlib logs `[W::bcf_get_version] Couldn't get VCF version, considering as
//! 4.2` and proceeds. We reproduce the same behavior here without modifying
//! the underlying parser: a small `Read` adapter inspects the first line of
//! the stream, rewrites a non-canonical `##fileformat=VCF<x>.<y>` line in
//! place, and emits the matching warning. All other bytes pass through.
//!
//! The wrapper is safe to apply to any uncompressed VCF text stream — if the
//! file does not begin with `##fileformat=VCF` it leaves the bytes untouched.

use std::io::{self, Read};

const FILEFORMAT_PREFIX: &[u8] = b"##fileformat=VCF";
const FILEFORMAT_SCAN_LIMIT: usize = 4096;

/// Read adapter that normalizes a non-canonical `##fileformat=VCF<x>.<y>`
/// header line on the first read. Subsequent bytes are forwarded as-is.
pub struct NormalizeFileformat<R> {
    inner: R,
    prefix: Vec<u8>,
    cursor: usize,
}

impl<R: Read> NormalizeFileformat<R> {
    /// Wrap `inner` and pre-read its first line for inspection. If the line is
    /// a non-canonical VCF fileformat declaration, rewrite it to the canonical
    /// `##fileformat=VCFv<x>.<y>` form and emit the upstream-compatible
    /// warning to stderr. The buffered bytes are returned to the consumer in
    /// the order they would have appeared.
    pub fn new(mut inner: R) -> io::Result<Self> {
        let mut prefix = Vec::with_capacity(64);
        let mut byte = [0u8; 1];
        while prefix.len() < FILEFORMAT_SCAN_LIMIT {
            let n = inner.read(&mut byte)?;
            if n == 0 {
                break;
            }
            prefix.push(byte[0]);
            if byte[0] == b'\n' {
                break;
            }
        }

        if let Some((normalized, version)) = normalize_fileformat_line(&prefix) {
            emit_version_warning(&version);
            prefix = normalized;
        }

        Ok(Self {
            inner,
            prefix,
            cursor: 0,
        })
    }
}

impl<R: Read> Read for NormalizeFileformat<R> {
    fn read(&mut self, dst: &mut [u8]) -> io::Result<usize> {
        if self.cursor < self.prefix.len() {
            let remaining = &self.prefix[self.cursor..];
            let n = remaining.len().min(dst.len());
            dst[..n].copy_from_slice(&remaining[..n]);
            self.cursor += n;
            return Ok(n);
        }
        self.inner.read(dst)
    }
}

/// If `line` is a non-canonical `##fileformat=VCF<x>.<y>` declaration,
/// return the canonical replacement (with the original line terminator
/// preserved) along with the version string we considered it as.
fn normalize_fileformat_line(line: &[u8]) -> Option<(Vec<u8>, String)> {
    let (content, terminator) = match line.iter().position(|&b| b == b'\n') {
        Some(i) => (&line[..i], &line[i..]),
        None => (line, &b""[..]),
    };
    let (content, cr) = match content.last() {
        Some(b'\r') => (&content[..content.len() - 1], &b"\r"[..]),
        _ => (content, &b""[..]),
    };

    if !content.starts_with(FILEFORMAT_PREFIX) {
        return None;
    }

    let rest = &content[FILEFORMAT_PREFIX.len()..];
    if rest.is_empty() {
        return None;
    }

    if rest[0] == b'v' || rest[0] == b'V' {
        return None;
    }

    let version = std::str::from_utf8(rest).ok()?.to_owned();

    let mut out = Vec::with_capacity(line.len() + 1);
    out.extend_from_slice(FILEFORMAT_PREFIX);
    out.push(b'v');
    out.extend_from_slice(rest);
    out.extend_from_slice(cr);
    out.extend_from_slice(terminator);

    Some((out, version))
}

fn emit_version_warning(version: &str) {
    eprintln!("[W::bcf_get_version] Couldn't get VCF version, considering as {version}");
}

/// In-place normalize a slurped VCF text. If the first line is a
/// non-canonical `##fileformat=VCF<x>.<y>` declaration, rewrite it to the
/// canonical `##fileformat=VCFv<x>.<y>` and emit the upstream-style warning.
/// Otherwise, leave the text untouched.
pub fn normalize_vcf_text(text: &mut String) {
    let first_line_end = text.find('\n').map(|i| i + 1).unwrap_or(text.len());
    let first_line = &text.as_bytes()[..first_line_end];
    if let Some((normalized, version)) = normalize_fileformat_line(first_line) {
        emit_version_warning(&version);
        let rest = text[first_line_end..].to_owned();
        text.clear();
        text.push_str(std::str::from_utf8(&normalized).unwrap_or(""));
        text.push_str(&rest);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Read};

    fn read_to_string<R: Read>(mut r: R) -> String {
        let mut s = String::new();
        r.read_to_string(&mut s).unwrap();
        s
    }

    #[test]
    fn normalizes_kestrel_fileformat_line() {
        let input = "##fileformat=VCF4.2\n##contig=<ID=chr1,length=10>\n";
        let wrapper = NormalizeFileformat::new(Cursor::new(input)).unwrap();
        assert_eq!(
            read_to_string(wrapper),
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10>\n"
        );
    }

    #[test]
    fn passes_through_canonical_fileformat_line() {
        let input = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10>\n";
        let wrapper = NormalizeFileformat::new(Cursor::new(input)).unwrap();
        assert_eq!(read_to_string(wrapper), input);
    }

    #[test]
    fn preserves_crlf_line_terminator() {
        let input = "##fileformat=VCF4.2\r\nrest\r\n";
        let wrapper = NormalizeFileformat::new(Cursor::new(input)).unwrap();
        assert_eq!(read_to_string(wrapper), "##fileformat=VCFv4.2\r\nrest\r\n");
    }

    #[test]
    fn handles_missing_line_terminator() {
        let input = "##fileformat=VCF4.2";
        let wrapper = NormalizeFileformat::new(Cursor::new(input)).unwrap();
        assert_eq!(read_to_string(wrapper), "##fileformat=VCFv4.2");
    }

    #[test]
    fn leaves_non_fileformat_first_line_untouched() {
        let input = "##contig=<ID=chr1>\n##fileformat=VCF4.2\n";
        let wrapper = NormalizeFileformat::new(Cursor::new(input)).unwrap();
        assert_eq!(read_to_string(wrapper), input);
    }

    #[test]
    fn does_not_double_prefix_already_canonical_form() {
        assert!(normalize_fileformat_line(b"##fileformat=VCFv4.2\n").is_none());
        assert!(normalize_fileformat_line(b"##fileformat=VCFV4.2\n").is_none());
    }

    #[test]
    fn normalize_vcf_text_rewrites_first_line_only() {
        let mut text = String::from("##fileformat=VCF4.2\n##fileformat=VCF4.2\nrest\n");
        normalize_vcf_text(&mut text);
        assert_eq!(text, "##fileformat=VCFv4.2\n##fileformat=VCF4.2\nrest\n");
    }

    #[test]
    fn normalize_vcf_text_no_op_for_canonical() {
        let original = "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n";
        let mut text = String::from(original);
        normalize_vcf_text(&mut text);
        assert_eq!(text, original);
    }
}
