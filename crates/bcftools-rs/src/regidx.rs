//! bcftools-shaped region index helpers.
//!
//! This is a thin facade over [`htslib_rs::regidx`] that keeps bcftools-rs
//! call sites from depending directly on the lower-level parser selection and
//! line/list loading details. Coordinates follow HTSlib's regidx convention:
//! stored intervals are 0-based inclusive, while tab/region text inputs are
//! parsed with their upstream 1-based semantics.

use std::{
    fs::File,
    io::{self, BufRead, BufReader},
    path::Path,
};

pub use htslib_rs::regidx::{ParseError, Parser, REGIDX_MAX, RegionRecord};

/// A bcftools-facing in-memory region index.
#[derive(Clone, Debug, Default)]
pub struct RegionIndex {
    inner: htslib_rs::regidx::RegionIndex,
}

impl RegionIndex {
    /// Creates an empty region index.
    pub fn new() -> Self {
        Self::default()
    }

    /// Loads a local region file, selecting BED parsing for `.bed` /
    /// `.bed.gz`-style names and tab parsing otherwise.
    pub fn from_path<P>(path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();
        let parser = parser_for_path(path);
        let mut index = Self::new();
        let file = File::open(path)?;
        for line in BufReader::new(file).lines() {
            index.insert_line(&line?, parser).map_err(parse_io_error)?;
        }
        Ok(index)
    }

    /// Builds an index from a delimiter-separated region list, using
    /// `Parser::Region` for each item.
    pub fn from_region_list(list: &str, delimiter: char) -> io::Result<Self> {
        let mut index = Self::new();
        for raw in list.split(delimiter) {
            let raw = raw.trim();
            if raw.is_empty() {
                continue;
            }
            index
                .insert_line(raw, Parser::Region)
                .map_err(parse_io_error)?;
        }
        Ok(index)
    }

    /// Parses and inserts one line with the selected parser.
    pub fn insert_line(&mut self, line: &str, parser: Parser) -> Result<bool, ParseError> {
        self.inner.insert_line(line, parser)
    }

    /// Inserts a pre-parsed record.
    pub fn push(&mut self, record: RegionRecord) {
        self.inner.push(record);
    }

    /// Returns all records overlapping `seq:start-end`, where coordinates are
    /// 0-based inclusive.
    pub fn overlaps<'a>(&'a self, seq: &str, start: i64, end: i64) -> Vec<&'a RegionRecord> {
        self.inner.overlaps(seq, start, end)
    }

    /// Returns whether any record overlaps `seq:start-end`, where coordinates
    /// are 0-based inclusive.
    pub fn has_overlap(&self, seq: &str, start: i64, end: i64) -> bool {
        self.inner.has_overlap(seq, start, end)
    }

    /// Iterates all indexed records in sequence-name order and coordinate
    /// order.
    pub fn iter(&self) -> impl Iterator<Item = &RegionRecord> {
        self.inner.iter()
    }

    /// Exposes the wrapped htslib-rs index for low-level consumers.
    pub fn as_inner(&self) -> &htslib_rs::regidx::RegionIndex {
        &self.inner
    }
}

/// Parses a BED-style line: `CHROM FROM TO`, 0-based, right-open.
pub fn parse_bed(line: &str) -> Result<Option<RegionRecord>, ParseError> {
    htslib_rs::regidx::parse_bed(line)
}

/// Parses a tabular line: `CHROM POS [TO]`, 1-based, inclusive.
pub fn parse_tab(line: &str) -> Result<Option<RegionRecord>, ParseError> {
    htslib_rs::regidx::parse_tab(line)
}

/// Parses a region line: `CHROM`, `CHROM:POS`, `CHROM:FROM-TO`, or
/// `CHROM:FROM-`.
pub fn parse_region_line(line: &str) -> Result<Option<RegionRecord>, ParseError> {
    htslib_rs::regidx::parse_region_line(line)
}

/// Selects the default parser bcftools uses for local region files.
pub fn parser_for_path(path: &Path) -> Parser {
    let raw = path.as_os_str().to_string_lossy().to_ascii_lowercase();
    if raw.ends_with(".bed") || raw.ends_with(".bed.gz") || raw.ends_with(".bed.bgz") {
        Parser::Bed
    } else {
        Parser::Tab
    }
}

fn parse_io_error(error: ParseError) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("failed to parse region line: {error:?}"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parser_uses_bed_for_bed_extensions() {
        assert_eq!(parser_for_path(Path::new("targets.bed")), Parser::Bed);
        assert_eq!(parser_for_path(Path::new("targets.bed.gz")), Parser::Bed);
        assert_eq!(parser_for_path(Path::new("targets.tsv")), Parser::Tab);
    }

    #[test]
    fn parses_bed_tab_and_region_shapes() {
        let bed = parse_bed("chr1\t10\t20\tpayload").unwrap().unwrap();
        assert_eq!((bed.seq.as_str(), bed.start, bed.end), ("chr1", 10, 19));
        assert_eq!(bed.payload.as_deref(), Some("payload"));

        let tab = parse_tab("chr1\t11\t20\tpayload").unwrap().unwrap();
        assert_eq!((tab.seq.as_str(), tab.start, tab.end), ("chr1", 10, 19));

        let region = parse_region_line("chr1:11-20").unwrap().unwrap();
        assert_eq!(
            (region.seq.as_str(), region.start, region.end),
            ("chr1", 10, 19)
        );
    }

    #[test]
    fn region_list_builds_queryable_index() {
        let index = RegionIndex::from_region_list("chr1:11-20,chr2:5", ',').unwrap();
        assert!(index.has_overlap("chr1", 10, 10));
        assert!(index.has_overlap("chr2", 4, 4));
        assert!(!index.has_overlap("chr2", 5, 5));
        assert_eq!(index.iter().count(), 2);
    }
}
