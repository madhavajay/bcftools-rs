//! TSV-to-VCF helper logic ported from bcftools `tsv2vcf.c`.
//!
//! The upstream helper maps a comma-separated column declaration to per-column
//! setters and then walks whitespace-separated input lines. This module keeps
//! that dispatcher shape and provides a small record model for the common
//! `CHROM`, `POS`, `ID`, `REF`, and `ALT` columns used by `convert --tsv2vcf`.

use std::{collections::HashMap, io};

/// Parsed column dispatcher.
#[derive(Clone, Debug)]
pub struct Tsv {
    columns: Vec<Option<String>>,
}

impl Tsv {
    /// Parses a comma-separated column list. A `-` column is skipped.
    pub fn new(columns: &str) -> Self {
        Self {
            columns: columns
                .split(',')
                .map(|column| {
                    let column = column.trim();
                    if column == "-" || column.is_empty() {
                        None
                    } else {
                        Some(column.to_string())
                    }
                })
                .collect(),
        }
    }

    /// Returns the configured column names. Skipped columns are `None`.
    pub fn columns(&self) -> &[Option<String>] {
        &self.columns
    }

    /// Walks one whitespace-separated line and calls registered setters for
    /// matching columns.
    ///
    /// Returns the number of setters that were called, matching upstream's
    /// "status" count. A line with no registered setters is an error.
    pub fn parse_with<T>(
        &self,
        record: &mut T,
        line: &str,
        setters: &SetterMap<T>,
    ) -> io::Result<usize> {
        let mut status = 0usize;
        for (idx, field) in line.split_whitespace().take(self.columns.len()).enumerate() {
            let Some(name) = self.columns[idx].as_deref() else {
                continue;
            };
            let Some(setter) = setters.get(&name.to_ascii_uppercase()) else {
                continue;
            };
            setter(record, field)?;
            status += 1;
        }

        if status == 0 {
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no TSV setters matched",
            ))
        } else {
            Ok(status)
        }
    }

    /// Parses one line into a simple VCF-shaped record.
    pub fn parse_record(&self, line: &str) -> io::Result<TsvRecord> {
        let mut record = TsvRecord::default();
        self.parse_with(&mut record, line, &default_setters())?;
        Ok(record)
    }
}

/// Setter callback registry.
pub type SetterMap<T> = HashMap<String, Box<dyn Fn(&mut T, &str) -> io::Result<()> + Send + Sync>>;

/// Minimal VCF-shaped record populated by TSV setters.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct TsvRecord {
    /// Chromosome / contig name.
    pub chrom: Option<String>,
    /// 0-based position.
    pub pos: Option<i64>,
    /// ID field.
    pub id: Option<String>,
    /// REF allele.
    pub ref_allele: Option<String>,
    /// ALT alleles.
    pub alt_alleles: Vec<String>,
    /// Other named columns not handled by core setters.
    pub extras: HashMap<String, String>,
}

/// Returns default setters for common VCF columns.
pub fn default_setters() -> SetterMap<TsvRecord> {
    let mut setters: SetterMap<TsvRecord> = HashMap::new();
    setters.insert(
        "CHROM".to_string(),
        Box::new(|record, field| {
            record.chrom = Some(field.to_string());
            Ok(())
        }),
    );
    setters.insert(
        "POS".to_string(),
        Box::new(|record, field| {
            let pos = field.parse::<i64>().map_err(|e| {
                io::Error::new(io::ErrorKind::InvalidData, format!("invalid POS: {e}"))
            })?;
            if pos < 1 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "POS must be 1-based",
                ));
            }
            record.pos = Some(pos - 1);
            Ok(())
        }),
    );
    setters.insert(
        "ID".to_string(),
        Box::new(|record, field| {
            record.id = Some(field.to_string());
            Ok(())
        }),
    );
    setters.insert(
        "REF".to_string(),
        Box::new(|record, field| {
            record.ref_allele = Some(field.to_string());
            Ok(())
        }),
    );
    setters.insert(
        "ALT".to_string(),
        Box::new(|record, field| {
            record.alt_alleles = field
                .split(',')
                .filter(|allele| {
                    *allele != "."
                        && record
                            .ref_allele
                            .as_deref()
                            .is_none_or(|ref_allele| *allele != ref_allele)
                })
                .map(ToOwned::to_owned)
                .collect();
            Ok(())
        }),
    );
    setters.insert(
        "AA".to_string(),
        Box::new(|record, field| {
            record.extras.insert("AA".to_string(), field.to_string());
            Ok(())
        }),
    );
    setters
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_columns_and_skips_dash() {
        let tsv = Tsv::new("ID,-,CHROM,POS,REF,ALT");
        assert_eq!(
            tsv.columns(),
            &[
                Some("ID".to_string()),
                None,
                Some("CHROM".to_string()),
                Some("POS".to_string()),
                Some("REF".to_string()),
                Some("ALT".to_string())
            ]
        );
    }

    #[test]
    fn parses_simple_vcf_shaped_record() {
        let tsv = Tsv::new("ID,CHROM,POS,REF,ALT");
        let record = tsv.parse_record("rs1 chr1 11 A C,G").unwrap();
        assert_eq!(record.id.as_deref(), Some("rs1"));
        assert_eq!(record.chrom.as_deref(), Some("chr1"));
        assert_eq!(record.pos, Some(10));
        assert_eq!(record.ref_allele.as_deref(), Some("A"));
        assert_eq!(record.alt_alleles, vec!["C".to_string(), "G".to_string()]);
    }

    #[test]
    fn skipped_columns_still_consume_input_fields() {
        let tsv = Tsv::new("-,CHROM,POS");
        let record = tsv.parse_record("ignored chr2 5 trailing").unwrap();
        assert_eq!(record.chrom.as_deref(), Some("chr2"));
        assert_eq!(record.pos, Some(4));
    }

    #[test]
    fn unknown_columns_are_ignored_unless_registered() {
        #[derive(Default)]
        struct Rec {
            value: Option<String>,
        }

        let tsv = Tsv::new("X");
        let mut rec = Rec::default();
        let mut setters: SetterMap<Rec> = HashMap::new();
        setters.insert(
            "X".to_string(),
            Box::new(|record, field| {
                record.value = Some(field.to_string());
                Ok(())
            }),
        );

        assert_eq!(tsv.parse_with(&mut rec, "abc", &setters).unwrap(), 1);
        assert_eq!(rec.value.as_deref(), Some("abc"));
    }

    #[test]
    fn no_matching_setters_is_error() {
        let tsv = Tsv::new("MISSING");
        let err = tsv.parse_record("value").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn invalid_position_is_error() {
        let tsv = Tsv::new("CHROM,POS");
        assert!(tsv.parse_record("chr1 x").is_err());
        assert!(tsv.parse_record("chr1 0").is_err());
    }
}
