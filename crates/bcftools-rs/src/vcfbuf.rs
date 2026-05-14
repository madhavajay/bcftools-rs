//! Windowed VCF record buffer modelled after upstream `vcfbuf.c`.
//!
//! The C implementation stores `bcf1_t` records and is used by commands that
//! need a short look-ahead window, such as `norm` and overlap-removal plugins.
//! This module keeps the buffering logic independent from any concrete VCF
//! record representation: callers provide a span and payload, then decide how
//! to translate flushed payloads back to their writer API.

use std::collections::VecDeque;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VcfSpan {
    pub contig: String,
    /// 0-based inclusive start coordinate.
    pub start: i64,
    /// 0-based exclusive end coordinate.
    pub end: i64,
}

impl VcfSpan {
    pub fn new(contig: impl Into<String>, start: i64, end: i64) -> Self {
        let end = end.max(start + 1);

        Self {
            contig: contig.into(),
            start,
            end,
        }
    }

    pub fn overlaps(&self, other: &Self) -> bool {
        self.contig == other.contig && self.start < other.end && other.start < self.end
    }

    fn is_before(&self, contig: &str, start: i64) -> bool {
        self.contig.as_str() < contig || (self.contig == contig && self.end <= start)
    }

    fn is_outside_window(&self, contig: &str, start: i64, window: i64) -> bool {
        self.contig.as_str() < contig || (self.contig == contig && self.end + window <= start)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BufferedRecord<T> {
    pub span: VcfSpan,
    pub payload: T,
}

#[derive(Debug, Clone)]
pub struct VcfBuffer<T> {
    records: VecDeque<BufferedRecord<T>>,
    window: i64,
}

impl<T> Default for VcfBuffer<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> VcfBuffer<T> {
    pub fn new() -> Self {
        Self {
            records: VecDeque::new(),
            window: 0,
        }
    }

    pub fn with_window(window: i64) -> Self {
        Self {
            records: VecDeque::new(),
            window: window.max(0),
        }
    }

    pub fn window(&self) -> i64 {
        self.window
    }

    pub fn set_window(&mut self, window: i64) {
        self.window = window.max(0);
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn clear(&mut self) {
        self.records.clear();
    }

    pub fn push(&mut self, span: VcfSpan, payload: T) {
        let index = self
            .records
            .iter()
            .position(|record| compare_span(&span, &record.span).is_lt())
            .unwrap_or(self.records.len());

        self.records.insert(index, BufferedRecord { span, payload });
    }

    pub fn iter(&self) -> impl Iterator<Item = &BufferedRecord<T>> {
        self.records.iter()
    }

    pub fn overlaps<'a>(
        &'a self,
        span: &'a VcfSpan,
    ) -> impl Iterator<Item = &'a BufferedRecord<T>> {
        self.records
            .iter()
            .filter(move |record| record.span.overlaps(span))
    }

    /// Removes records that end before `contig:start`.
    pub fn drain_before(&mut self, contig: &str, start: i64) -> Vec<BufferedRecord<T>> {
        let mut out = Vec::new();

        while self
            .records
            .front()
            .is_some_and(|record| record.span.is_before(contig, start))
        {
            out.push(self.records.pop_front().expect("front checked"));
        }

        out
    }

    /// Removes records that cannot overlap records at `contig:start` after
    /// accounting for the configured look-ahead window.
    pub fn drain_outside_window(&mut self, contig: &str, start: i64) -> Vec<BufferedRecord<T>> {
        let mut out = Vec::new();

        while self
            .records
            .front()
            .is_some_and(|record| record.span.is_outside_window(contig, start, self.window))
        {
            out.push(self.records.pop_front().expect("front checked"));
        }

        out
    }

    pub fn drain_all(&mut self) -> Vec<BufferedRecord<T>> {
        self.records.drain(..).collect()
    }
}

fn compare_span(a: &VcfSpan, b: &VcfSpan) -> std::cmp::Ordering {
    a.contig
        .cmp(&b.contig)
        .then_with(|| a.start.cmp(&b.start))
        .then_with(|| a.end.cmp(&b.end))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn span(contig: &str, start: i64, end: i64) -> VcfSpan {
        VcfSpan::new(contig, start, end)
    }

    #[test]
    fn span_end_is_normalized_to_non_empty_interval() {
        assert_eq!(span("chr1", 10, 10).end, 11);
        assert_eq!(span("chr1", 10, 5).end, 11);
    }

    #[test]
    fn buffer_keeps_records_sorted_by_coordinate() {
        let mut buf = VcfBuffer::new();
        buf.push(span("chr2", 1, 2), "d");
        buf.push(span("chr1", 20, 21), "c");
        buf.push(span("chr1", 10, 11), "a");
        buf.push(span("chr1", 10, 12), "b");

        let payloads: Vec<_> = buf.iter().map(|record| record.payload).collect();

        assert_eq!(payloads, ["a", "b", "c", "d"]);
    }

    #[test]
    fn overlaps_reports_half_open_interval_intersections() {
        let mut buf = VcfBuffer::new();
        buf.push(span("chr1", 5, 10), "left");
        buf.push(span("chr1", 10, 15), "right");
        buf.push(span("chr2", 9, 12), "other");

        let query = span("chr1", 9, 11);
        let payloads: Vec<_> = buf.overlaps(&query).map(|record| record.payload).collect();

        assert_eq!(payloads, ["left", "right"]);
    }

    #[test]
    fn drain_before_removes_records_that_cannot_overlap_position() {
        let mut buf = VcfBuffer::new();
        buf.push(span("chr1", 1, 5), "a");
        buf.push(span("chr1", 5, 8), "b");
        buf.push(span("chr1", 8, 9), "c");

        let drained: Vec<_> = buf
            .drain_before("chr1", 8)
            .into_iter()
            .map(|record| record.payload)
            .collect();
        let remaining: Vec<_> = buf.iter().map(|record| record.payload).collect();

        assert_eq!(drained, ["a", "b"]);
        assert_eq!(remaining, ["c"]);
    }

    #[test]
    fn drain_outside_window_keeps_nearby_records() {
        let mut buf = VcfBuffer::with_window(5);
        buf.push(span("chr1", 1, 5), "old");
        buf.push(span("chr1", 8, 10), "near");
        buf.push(span("chr1", 14, 16), "current");

        let drained: Vec<_> = buf
            .drain_outside_window("chr1", 14)
            .into_iter()
            .map(|record| record.payload)
            .collect();
        let remaining: Vec<_> = buf.iter().map(|record| record.payload).collect();

        assert_eq!(drained, ["old"]);
        assert_eq!(remaining, ["near", "current"]);
    }

    #[test]
    fn contig_advance_flushes_previous_contigs() {
        let mut buf = VcfBuffer::new();
        buf.push(span("chr1", 100, 101), "old-contig");
        buf.push(span("chr2", 1, 2), "new-contig");

        let drained: Vec<_> = buf
            .drain_before("chr2", 0)
            .into_iter()
            .map(|record| record.payload)
            .collect();

        assert_eq!(drained, ["old-contig"]);
        assert_eq!(buf.len(), 1);
    }
}
