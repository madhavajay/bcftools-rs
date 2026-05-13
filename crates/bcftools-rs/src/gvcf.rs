//! gVCF block helpers.
//!
//! This ports the state-machine behavior of `gvcf.c`: parse depth ranges,
//! decide whether reference records can collapse into a block, maintain block
//! aggregate DP/PL/MIN_DP state, and flush block records when contiguity or
//! range constraints break. Concrete VCF/BCF record mutation is represented by
//! [`GvcfRecord`] so callers can bridge it to their writer API.

use std::io;

/// Header lines added by `gvcf_update_header`.
pub const HEADER_LINES: [&str; 2] = [
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
    "##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum per-sample depth in this gVCF block\">",
];

/// Minimal record model consumed and produced by [`Gvcf`].
#[derive(Debug, Clone, PartialEq)]
pub struct GvcfRecord {
    /// Numeric reference ID.
    pub rid: i32,
    /// 0-based start position.
    pub pos: i32,
    /// REF plus ALT alleles.
    pub alleles: Vec<String>,
    /// Per-sample FORMAT/DP values.
    pub format_dp: Option<Vec<i32>>,
    /// Flat per-sample FORMAT/PL values. The current bcftools algorithm expects
    /// three values per sample for collapsible reference blocks.
    pub format_pl: Option<Vec<i32>>,
    /// Optional INFO/QS values retained from the first record in a block.
    pub info_qs: Option<Vec<f32>>,
    /// Optional encoded genotypes retained from the first record in a block.
    pub genotypes: Option<Vec<i32>>,
    /// Optional 1-based INFO/END. If absent, `pos` is used as the block end.
    pub info_end: Option<i32>,
    /// INFO/MIN_DP value set on non-collapsed reference records.
    pub min_dp: Option<i32>,
}

/// Action returned from [`Gvcf::write`].
#[derive(Debug, Clone, PartialEq)]
pub enum GvcfAction {
    /// No record should be emitted yet.
    None,
    /// Emit an input record, possibly with `MIN_DP` set.
    EmitRecord(GvcfRecord),
    /// Emit a collapsed gVCF block.
    EmitBlock(GvcfRecord),
    /// Emit a block first, then emit the current non-collapsible record.
    EmitBlockAndRecord(GvcfRecord, GvcfRecord),
}

/// gVCF collapse state.
#[derive(Debug, Clone, PartialEq)]
pub struct Gvcf {
    dp_ranges: Vec<i32>,
    prev_range: Option<usize>,
    dp: Vec<i32>,
    pl: Option<Vec<i32>>,
    info_qs: Option<Vec<f32>>,
    genotypes: Option<Vec<i32>>,
    rid: i32,
    start: i32,
    end: i32,
    min_dp: i32,
    alleles: Vec<String>,
}

impl Gvcf {
    /// Parses comma-separated DP thresholds, matching `gvcf_init`.
    pub fn new(dp_ranges: &str) -> io::Result<Self> {
        let mut ranges = Vec::new();
        for part in dp_ranges.split(',') {
            if part.is_empty() {
                return Err(invalid_ranges(dp_ranges));
            }
            ranges.push(part.parse::<i32>().map_err(|_| invalid_ranges(dp_ranges))?);
        }
        if ranges.is_empty() {
            return Err(invalid_ranges(dp_ranges));
        }

        Ok(Self {
            dp_ranges: ranges,
            prev_range: None,
            dp: Vec::new(),
            pl: None,
            info_qs: None,
            genotypes: None,
            rid: -1,
            start: 0,
            end: 0,
            min_dp: 0,
            alleles: Vec::new(),
        })
    }

    /// Processes a record. Pass `None` to flush at end of stream.
    pub fn write(&mut self, record: Option<GvcfRecord>, is_ref: bool) -> io::Result<GvcfAction> {
        let mut can_collapse = is_ref && record.is_some();
        let mut dp_range = 0;
        let mut min_dp = 0;
        let mut needs_flush = !can_collapse;

        if let Some(record) = record.as_ref()
            && can_collapse
        {
            match record.format_dp.as_ref() {
                Some(dp) if !dp.is_empty() => {
                    min_dp = *dp.iter().min().expect("non-empty DP");
                    dp_range = self
                        .dp_ranges
                        .iter()
                        .position(|threshold| min_dp < *threshold)
                        .unwrap_or(self.dp_ranges.len());
                    if dp_range == 0 {
                        needs_flush = true;
                        can_collapse = false;
                    }
                }
                _ => needs_flush = true,
            }
        }

        if self.prev_range.is_some_and(|prev| prev != dp_range) {
            needs_flush = true;
        }
        if let Some(record) = record.as_ref() {
            if self.prev_range.is_some() && (self.rid != record.rid || record.pos > self.end + 1) {
                needs_flush = true;
            }
        } else {
            needs_flush = true;
        }

        let flushed = if self.prev_range.is_some() && needs_flush {
            Some(self.flush(record.as_ref()))
        } else {
            None
        };

        if can_collapse {
            let record = record.expect("record exists when can_collapse");
            self.extend_block(&record, dp_range, min_dp)?;
            return Ok(flushed.map_or(GvcfAction::None, GvcfAction::EmitBlock));
        }

        let current = record.map(|mut record| {
            if is_ref && min_dp != 0 {
                record.min_dp = Some(min_dp);
            }
            record
        });

        Ok(match (flushed, current) {
            (Some(block), Some(record)) => GvcfAction::EmitBlockAndRecord(block, record),
            (Some(block), None) => GvcfAction::EmitBlock(block),
            (None, Some(record)) => GvcfAction::EmitRecord(record),
            (None, None) => GvcfAction::None,
        })
    }

    fn extend_block(
        &mut self,
        record: &GvcfRecord,
        dp_range: usize,
        min_dp: i32,
    ) -> io::Result<()> {
        if self.prev_range.is_none() {
            self.dp = record
                .format_dp
                .clone()
                .expect("DP checked before collapse");
            self.pl = record.format_pl.clone();
            self.info_qs = record.info_qs.clone();
            self.genotypes = record.genotypes.clone();
            self.rid = record.rid;
            self.start = record.pos;
            self.alleles = record.alleles.clone();
            self.min_dp = min_dp;
        } else {
            self.min_dp = self.min_dp.min(min_dp);
            let current_dp = record
                .format_dp
                .as_ref()
                .expect("DP checked before collapse");
            for (stored, current) in self.dp.iter_mut().zip(current_dp) {
                *stored = (*stored).min(*current);
            }
            match (&mut self.pl, &record.format_pl) {
                (Some(stored), Some(current)) => merge_pl(stored, current)?,
                (Some(_), None) => self.pl = None,
                _ => {}
            }
        }

        self.prev_range = Some(dp_range);
        self.end = record.info_end.map(|end| end - 1).unwrap_or(record.pos);
        Ok(())
    }

    fn flush(&mut self, next: Option<&GvcfRecord>) -> GvcfRecord {
        if let Some(next) = next
            && next.rid == self.rid
            && next.pos == self.end
        {
            self.end -= 1;
        }

        let end_1based = self.end + 1;
        let mut record = GvcfRecord {
            rid: self.rid,
            pos: self.start,
            alleles: self.alleles.clone(),
            format_dp: Some(self.dp.clone()),
            format_pl: self.pl.clone(),
            info_qs: self.info_qs.clone(),
            genotypes: self.genotypes.clone(),
            info_end: None,
            min_dp: Some(self.min_dp),
        };
        if self.start + 1 < end_1based {
            record.info_end = Some(end_1based);
        }

        self.prev_range = None;
        self.rid = -1;
        self.pl = None;
        self.info_qs = None;
        self.genotypes = None;
        record
    }
}

fn merge_pl(stored: &mut [i32], current: &[i32]) -> io::Result<()> {
    if stored.len() != current.len() || !stored.len().is_multiple_of(3) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "unexpected number of PL fields",
        ));
    }
    for (stored, current) in stored.chunks_exact_mut(3).zip(current.chunks_exact(3)) {
        if stored[1] > current[1] {
            stored[1] = current[1];
            stored[2] = current[2];
        } else if stored[1] == current[1] && stored[2] > current[2] {
            stored[2] = current[2];
        }
    }
    Ok(())
}

fn invalid_ranges(ranges: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("invalid gVCF DP ranges: {ranges}"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(pos: i32, dp: &[i32], pl: &[i32]) -> GvcfRecord {
        GvcfRecord {
            rid: 0,
            pos,
            alleles: vec!["A".into(), "<*>".into()],
            format_dp: Some(dp.to_vec()),
            format_pl: Some(pl.to_vec()),
            info_qs: Some(vec![1.0, 2.0]),
            genotypes: Some(vec![0, 0]),
            info_end: None,
            min_dp: None,
        }
    }

    #[test]
    fn parses_dp_ranges() {
        assert_eq!(Gvcf::new("1,5,10").unwrap().dp_ranges, [1, 5, 10]);
        assert!(Gvcf::new("1,").is_err());
        assert!(Gvcf::new("abc").is_err());
    }

    #[test]
    fn collapses_contiguous_reference_records_and_flushes_at_end() {
        let mut gvcf = Gvcf::new("1,5,10").unwrap();
        assert_eq!(
            gvcf.write(Some(rec(0, &[6, 7], &[0, 5, 9, 0, 4, 8])), true)
                .unwrap(),
            GvcfAction::None
        );
        assert_eq!(
            gvcf.write(Some(rec(1, &[5, 8], &[0, 3, 7, 0, 4, 6])), true)
                .unwrap(),
            GvcfAction::None
        );

        let GvcfAction::EmitBlock(block) = gvcf.write(None, false).unwrap() else {
            panic!("expected flushed block");
        };
        assert_eq!(block.pos, 0);
        assert_eq!(block.info_end, Some(2));
        assert_eq!(block.min_dp, Some(5));
        assert_eq!(block.format_dp, Some(vec![5, 7]));
        assert_eq!(block.format_pl, Some(vec![0, 3, 7, 0, 4, 6]));
        assert_eq!(block.info_qs, Some(vec![1.0, 2.0]));
        assert_eq!(block.genotypes, Some(vec![0, 0]));
    }

    #[test]
    fn low_depth_reference_record_is_emitted_with_min_dp() {
        let mut gvcf = Gvcf::new("5,10").unwrap();
        let input = rec(0, &[2, 3], &[0, 1, 2, 0, 1, 2]);
        let GvcfAction::EmitRecord(output) = gvcf.write(Some(input), true).unwrap() else {
            panic!("expected original record");
        };
        assert_eq!(output.min_dp, Some(2));
    }

    #[test]
    fn range_change_flushes_then_starts_new_block() {
        let mut gvcf = Gvcf::new("1,5,10").unwrap();
        assert_eq!(
            gvcf.write(Some(rec(0, &[6], &[0, 5, 9])), true).unwrap(),
            GvcfAction::None
        );
        let action = gvcf.write(Some(rec(1, &[12], &[0, 2, 4])), true).unwrap();
        let GvcfAction::EmitBlock(block) = action else {
            panic!("expected flushed block");
        };
        assert_eq!(block.pos, 0);
        assert_eq!(block.info_end, None);

        let GvcfAction::EmitBlock(block) = gvcf.write(None, false).unwrap() else {
            panic!("expected second block");
        };
        assert_eq!(block.pos, 1);
    }

    #[test]
    fn non_reference_record_after_block_emits_both() {
        let mut gvcf = Gvcf::new("1,5,10").unwrap();
        assert_eq!(
            gvcf.write(Some(rec(0, &[6], &[0, 5, 9])), true).unwrap(),
            GvcfAction::None
        );
        let mut variant = rec(2, &[6], &[0, 5, 9]);
        variant.alleles = vec!["A".into(), "C".into()];
        let GvcfAction::EmitBlockAndRecord(block, output) =
            gvcf.write(Some(variant.clone()), false).unwrap()
        else {
            panic!("expected block and variant");
        };
        assert_eq!(block.pos, 0);
        assert_eq!(output, variant);
    }

    #[test]
    fn adjacent_non_reference_at_block_end_trims_end() {
        let mut gvcf = Gvcf::new("1,5,10").unwrap();
        let mut first = rec(0, &[6], &[0, 5, 9]);
        first.info_end = Some(2);
        assert_eq!(gvcf.write(Some(first), true).unwrap(), GvcfAction::None);
        let variant = rec(1, &[6], &[0, 5, 9]);
        let GvcfAction::EmitBlockAndRecord(block, _) = gvcf.write(Some(variant), false).unwrap()
        else {
            panic!("expected block and record");
        };
        assert_eq!(block.pos, 0);
        assert_eq!(block.info_end, None);
    }
}
