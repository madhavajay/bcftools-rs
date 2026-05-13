//! bcftools ploidy specification helpers.
//!
//! This ports the public behavior of `ploidy.c`: parse records of the form
//! `CHROM FROM TO SEX PLOIDY`, support `* * * SEX PLOIDY` per-sex defaults,
//! query 0-based positions, and preserve sex IDs in first-seen order.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone, PartialEq, Eq)]
struct PloidyRecord {
    seq: String,
    start: u32,
    end: u32,
    sex: usize,
    ploidy: i32,
}

/// Result of a ploidy query.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryResult {
    /// Whether the queried position overlapped at least one explicit region.
    pub listed: bool,
    /// Ploidy by numeric sex ID.
    pub sex2ploidy: Vec<i32>,
    /// Minimum ploidy encountered for this query.
    pub min: i32,
    /// Maximum ploidy encountered for this query.
    pub max: i32,
}

/// Ploidy lookup table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Ploidy {
    default_ploidy: i32,
    min: i32,
    max: i32,
    sex2id: HashMap<String, usize>,
    id2sex: Vec<String>,
    sex_defaults: Vec<Option<i32>>,
    records: Vec<PloidyRecord>,
}

impl Ploidy {
    /// Reads a ploidy specification from a local file.
    pub fn from_path(path: impl AsRef<Path>, default_ploidy: i32) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut ploidy = Self::empty(default_ploidy);
        for line in BufReader::new(file).lines() {
            ploidy.insert_line(&line?)?;
        }
        ploidy.set_defaults();
        Ok(ploidy)
    }

    /// Parses a complete ploidy specification from a string.
    pub fn from_string(s: &str, default_ploidy: i32) -> io::Result<Self> {
        let mut ploidy = Self::empty(default_ploidy);
        for line in s.lines() {
            ploidy.insert_line(line)?;
        }
        ploidy.set_defaults();
        Ok(ploidy)
    }

    /// Builds an empty table with no recognised sexes.
    pub fn empty(default_ploidy: i32) -> Self {
        Self {
            default_ploidy,
            min: -1,
            max: -1,
            sex2id: HashMap::new(),
            id2sex: Vec::new(),
            sex_defaults: Vec::new(),
            records: Vec::new(),
        }
    }

    /// Adds a sex name and returns its numeric ID.
    pub fn add_sex(&mut self, sex: &str) -> usize {
        self.intern_sex(sex, Some(self.default_ploidy))
    }

    fn intern_sex(&mut self, sex: &str, default: Option<i32>) -> usize {
        if let Some(&id) = self.sex2id.get(sex) {
            return id;
        }
        let id = self.id2sex.len();
        self.id2sex.push(sex.to_string());
        self.sex2id.insert(sex.to_string(), id);
        self.sex_defaults.push(default);
        id
    }

    /// Number of recognised sex names.
    pub fn nsex(&self) -> usize {
        self.id2sex.len()
    }

    /// Maps numeric sex ID to sex name.
    pub fn id2sex(&self, id: usize) -> Option<&str> {
        self.id2sex.get(id).map(String::as_str)
    }

    /// Maps sex name to numeric sex ID.
    pub fn sex2id(&self, sex: &str) -> Option<usize> {
        self.sex2id.get(sex).copied()
    }

    /// Minimum recognised ploidy, including defaults.
    pub fn min(&self) -> i32 {
        self.default_ploidy.min(self.min)
    }

    /// Maximum recognised ploidy, including defaults.
    pub fn max(&self) -> i32 {
        self.default_ploidy.max(self.max)
    }

    /// Query ploidy for all registered sexes at a 0-based position.
    pub fn query(&self, seq: &str, pos: u32) -> QueryResult {
        let overlaps: Vec<_> = self
            .records
            .iter()
            .filter(|record| record.seq == seq && record.start <= pos && pos <= record.end)
            .collect();

        if overlaps.is_empty() {
            let sex2ploidy = self
                .sex_defaults
                .iter()
                .map(|ploidy| ploidy.unwrap_or(self.default_ploidy))
                .collect();
            return QueryResult {
                listed: false,
                sex2ploidy,
                min: self.default_ploidy,
                max: self.default_ploidy,
            };
        }

        let mut sex2ploidy = vec![self.default_ploidy; self.nsex()];
        let mut min = i32::MAX;
        let mut max = -1;
        for record in overlaps {
            if record.ploidy != self.default_ploidy {
                sex2ploidy[record.sex] = record.ploidy;
                min = min.min(record.ploidy);
                max = max.max(record.ploidy);
            }
        }
        if max == -1 {
            min = self.default_ploidy;
            max = self.default_ploidy;
        }

        QueryResult {
            listed: true,
            sex2ploidy,
            min,
            max,
        }
    }

    /// Formats a parseable ploidy specification for debugging.
    pub fn format_spec(&self) -> String {
        let mut out = String::new();
        for record in &self.records {
            out.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                record.seq,
                record.start + 1,
                record.end + 1,
                self.id2sex[record.sex],
                record.ploidy
            ));
        }
        for (id, sex) in self.id2sex.iter().enumerate() {
            out.push_str(&format!(
                "*\t*\t*\t{}\t{}\n",
                sex,
                self.sex_defaults[id].unwrap_or(self.default_ploidy)
            ));
        }
        out
    }

    fn insert_line(&mut self, line: &str) -> io::Result<()> {
        let line = line.trim();
        if line.is_empty() {
            return Ok(());
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != 5 {
            return Err(invalid_line(line));
        }

        let sex = self.intern_sex(fields[3], None);
        let ploidy = fields[4].parse::<i32>().map_err(|_| invalid_line(line))?;
        self.observe_ploidy(ploidy);

        if fields[0] == "*" {
            if fields[1] != "*" || fields[2] != "*" {
                return Err(invalid_line(line));
            }
            self.sex_defaults[sex] = Some(ploidy);
            return Ok(());
        }

        let from = fields[1].parse::<u32>().map_err(|_| invalid_line(line))?;
        let to = fields[2].parse::<u32>().map_err(|_| invalid_line(line))?;
        if from == 0 || to == 0 || to < from {
            return Err(invalid_line(line));
        }

        self.records.push(PloidyRecord {
            seq: fields[0].to_string(),
            start: from - 1,
            end: to - 1,
            sex,
            ploidy,
        });
        Ok(())
    }

    fn set_defaults(&mut self) {
        if let Some(star) = self.sex2id("*").and_then(|id| self.sex_defaults[id]) {
            self.default_ploidy = star;
        }
        for default in &mut self.sex_defaults {
            if default.is_none() {
                *default = Some(self.default_ploidy);
            }
        }
        self.observe_ploidy(self.default_ploidy);
    }

    fn observe_ploidy(&mut self, ploidy: i32) {
        if self.min < 0 || ploidy < self.min {
            self.min = ploidy;
        }
        if self.max < 0 || ploidy > self.max {
            self.max = ploidy;
        }
    }
}

fn invalid_line(line: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("could not parse ploidy line: {line}"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPEC: &str = "\
X 1 60000 M 1
X 2699521 154931043 M 1
Y 1 59373566 M 1
Y 1 59373566 F 0
MT 1 16569 M 1
MT 1 16569 F 1
* * * M 2
* * * F 2
";

    #[test]
    fn parses_sexes_defaults_and_min_max() {
        let ploidy = Ploidy::from_string(SPEC, 2).unwrap();
        assert_eq!(ploidy.nsex(), 2);
        assert_eq!(ploidy.id2sex(0), Some("M"));
        assert_eq!(ploidy.id2sex(1), Some("F"));
        assert_eq!(ploidy.sex2id("M"), Some(0));
        assert_eq!(ploidy.sex2id("missing"), None);
        assert_eq!(ploidy.min(), 0);
        assert_eq!(ploidy.max(), 2);
    }

    #[test]
    fn query_uses_region_ploidy_and_defaults() {
        let ploidy = Ploidy::from_string(SPEC, 2).unwrap();

        let chr_x = ploidy.query("X", 59_999);
        assert!(chr_x.listed);
        assert_eq!(chr_x.sex2ploidy, [1, 2]);
        assert_eq!((chr_x.min, chr_x.max), (1, 1));

        let chr_y = ploidy.query("Y", 0);
        assert!(chr_y.listed);
        assert_eq!(chr_y.sex2ploidy, [1, 0]);
        assert_eq!((chr_y.min, chr_y.max), (0, 1));

        let chr1 = ploidy.query("1", 123);
        assert!(!chr1.listed);
        assert_eq!(chr1.sex2ploidy, [2, 2]);
        assert_eq!((chr1.min, chr1.max), (2, 2));
    }

    #[test]
    fn star_sex_overrides_global_default() {
        let ploidy = Ploidy::from_string("* * * * 1\nX 1 10 F 0\n", 2).unwrap();
        assert_eq!(ploidy.min(), 0);
        assert_eq!(ploidy.max(), 1);
        assert_eq!(ploidy.query("1", 0).sex2ploidy, [1, 1]);
    }

    #[test]
    fn add_sex_registers_default_ploidy() {
        let mut ploidy = Ploidy::from_string("* * * M 2\n", 2).unwrap();
        let id = ploidy.add_sex("F");
        assert_eq!(id, 1);
        assert_eq!(ploidy.query("1", 0).sex2ploidy, [2, 2]);
    }

    #[test]
    fn format_spec_round_trips() {
        let ploidy = Ploidy::from_string(SPEC, 2).unwrap();
        let formatted = ploidy.format_spec();
        let reparsed = Ploidy::from_string(&formatted, 2).unwrap();
        assert_eq!(reparsed, ploidy);
    }

    #[test]
    fn rejects_malformed_lines() {
        assert!(Ploidy::from_string("X 10 1 M 2\n", 2).is_err());
        assert!(Ploidy::from_string("* 1 * M 2\n", 2).is_err());
        assert!(Ploidy::from_string("X 1 10 M nope\n", 2).is_err());
    }
}
