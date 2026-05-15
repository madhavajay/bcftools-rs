//! `bcftools +allele-length` (upstream `bcftools/plugins/allele-length.c`).
//!
//! Counts the frequency of REF, first-ALT, and REF+ALT allele lengths into
//! fixed `MAXLEN`-bounded histograms, plus a "non-base" tally for alleles
//! containing characters outside `[ACGTacgt]`. Only the first ALT allele is
//! considered, matching upstream's `rec->d.allele[1]` access. Text report.

use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

const MAXLEN: usize = 512;

struct Hist {
    reflen: Vec<u64>,
    altlen: Vec<u64>,
    refaltlen: Vec<u64>,
    xrefaltlen: Vec<u64>,
    numvar: u64,
    numxvar: u64,
}

impl Hist {
    fn new() -> Self {
        Self {
            reflen: vec![0; MAXLEN],
            altlen: vec![0; MAXLEN],
            refaltlen: vec![0; MAXLEN],
            xrefaltlen: vec![0; MAXLEN],
            numvar: 0,
            numxvar: 0,
        }
    }

    fn report(&self) -> String {
        let mut out = String::with_capacity(MAXLEN * 16);
        out.push_str("LENGTH\tREF\tALT\tREF+ALT\tREF+ALT WITH NON-BASE NUCLEOTIDES\n");
        for i in 0..MAXLEN {
            let _ = writeln!(
                out,
                "{i}\t{}\t{}\t{}\t{}",
                self.reflen[i], self.altlen[i], self.refaltlen[i], self.xrefaltlen[i]
            );
        }
        let _ = writeln!(out, "\t\t\t{}\t{}", self.numvar, self.numxvar);
        out
    }
}

fn contain_non_base(s: &str) -> bool {
    s.bytes()
        .any(|c| !matches!(c, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
}

/// Reads the input VCF/BCF and returns the upstream-shaped histogram report.
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(tally(&text).report())
}

fn tally(text: &str) -> Hist {
    let mut h = Hist::new();
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let (Some(_c), Some(_p), Some(_i), Some(reference), Some(alt)) =
            (f.next(), f.next(), f.next(), f.next(), f.next())
        else {
            continue;
        };
        let first_alt = match alt {
            "." | "" => "",
            _ => alt.split(',').next().unwrap_or(""),
        };
        let rl = reference.len().min(MAXLEN - 1);
        let al = first_alt.len().min(MAXLEN - 1);
        let ral = (reference.len() + first_alt.len()).min(MAXLEN - 1);
        h.reflen[rl] += 1;
        h.altlen[al] += 1;
        h.refaltlen[ral] += 1;
        if contain_non_base(reference) || contain_non_base(first_alt) {
            h.xrefaltlen[ral] += 1;
            h.numxvar += 1;
        }
        h.numvar += 1;
    }
    h
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    if path == Path::new("-") {
        let tmp = stdin_tmp_path();
        let mut data = Vec::new();
        io::stdin().lock().read_to_end(&mut data)?;
        fs::write(&tmp, data)?;
        let result = read_vcf_text(&tmp);
        let _ = fs::remove_file(&tmp);
        return result;
    }

    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }

    let mut text = String::new();
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let file = File::open(path)?;
        let mut dec = MultiGzDecoder::new(file);
        dec.read_to_string(&mut text)?;
    } else {
        text = fs::read_to_string(path)?;
    }
    normalize_vcf_text(&mut text);
    Ok(text)
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-allele-length-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_base_detection() {
        assert!(!contain_non_base("ACGTacgt"));
        assert!(contain_non_base("ACGN"));
        assert!(contain_non_base("<DEL>"));
        assert!(contain_non_base("A.T"));
    }

    #[test]
    fn tallies_lengths_first_alt_only() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tTA\tT\t.\t.\t.\n\
1\t2\t.\tC\tT,G\t.\t.\t.\n\
1\t3\t.\tC\tN\t.\t.\t.\n";
        let h = tally(vcf);
        // Row1: REF len 2, ALT len 1, REF+ALT 3.
        // Row2: REF len 1, first ALT len 1, REF+ALT 2.
        // Row3: REF len 1, ALT 'N' len 1 (non-base), REF+ALT 2.
        assert_eq!(h.reflen[2], 1);
        assert_eq!(h.reflen[1], 2);
        assert_eq!(h.altlen[1], 3);
        assert_eq!(h.refaltlen[3], 1);
        assert_eq!(h.refaltlen[2], 2);
        assert_eq!(h.xrefaltlen[2], 1);
        assert_eq!(h.numvar, 3);
        assert_eq!(h.numxvar, 1);
    }

    #[test]
    fn report_shape() {
        let h = Hist::new();
        let r = h.report();
        let lines: Vec<&str> = r.lines().collect();
        assert_eq!(
            lines[0],
            "LENGTH\tREF\tALT\tREF+ALT\tREF+ALT WITH NON-BASE NUCLEOTIDES"
        );
        assert_eq!(lines[1], "0\t0\t0\t0\t0");
        // header + 512 rows + 1 summary line.
        assert_eq!(lines.len(), 1 + MAXLEN + 1);
        assert_eq!(*lines.last().unwrap(), "\t\t\t0\t0");
    }
}
