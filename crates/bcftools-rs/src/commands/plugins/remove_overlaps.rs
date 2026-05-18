//! `bcftools +remove-overlaps` (upstream `bcftools/plugins/remove-overlaps.c`
//! plus the `MARK_OVERLAP`/`MARK_DUP` paths of `bcftools/vcfbuf.c`).
//!
//! Faithful port of the streaming overlap/duplicate mark state machine: a
//! FIFO record buffer with a parallel mark buffer, the `overlap_rid`/
//! `overlap_end` running span, the left-aligned-indel `imin` shared-prefix
//! adjustment, and the `can_flush` drain logic. Records are emitted oldest
//! first exactly as upstream does, so output order and marking match
//! byte-for-byte.
//!
//! Implemented modes: `-m overlap`, `-m dup`, `-m 'min(QUAL)'`,
//! `-M TAG`, `--reverse`, `-O t` (plain `chr<TAB>pos` site list), and
//! both `--missing N` (scalar) and `--missing DP` (the max-QUAL/DP
//! coverage-scaling heuristic) for `min(QUAL)`. `-m 'min(QUAL)'`
//! greedily marks the lowest-QUAL record in each connected overlap
//! component until no overlaps remain (upstream `vcfbuf` `MARK_EXPR`),
//! with the mark-tag `##INFO` appended last per `bcf_hdr_printf`. All
//! upstream `remove-overlaps.*.out` fixtures pass byte-for-byte; only
//! `-i`/`-e` record filtering remains (tracked in `TODO.md`).

use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mark {
    Overlap,
    Dup,
    /// `-m 'min(QUAL)'`: within each overlap cluster, mark the
    /// lowest-QUAL records for removal until no overlaps remain
    /// (upstream `vcfbuf` `MARK_EXPR`).
    MinQual,
}

impl Mark {
    /// Parses the `-m`/`--mark` value.
    pub fn parse(expr: &str) -> Result<Mark, String> {
        let e = expr.replace(' ', "");
        if expr.eq_ignore_ascii_case("overlap") {
            Ok(Mark::Overlap)
        } else if expr.eq_ignore_ascii_case("dup") {
            Ok(Mark::Dup)
        } else if e.eq_ignore_ascii_case("min(QUAL)") {
            Ok(Mark::MinQual)
        } else {
            Err(format!(
                "remove-overlaps -m '{expr}' is not supported in this local slice (only 'overlap', 'dup', 'min(QUAL)')"
            ))
        }
    }
}

/// A buffered record: the raw line plus the parsed fields the state machine
/// needs (everything else is preserved verbatim on output).
struct Rec {
    line: String,
    chrom: String,
    pos0: i64,
    rlen: i64,
    /// Minimum shared REF/ALT prefix across non-symbolic alleles.
    imin: i64,
    /// QUAL (`None` when `.`); used by `-m 'min(QUAL)'`.
    qual: Option<f32>,
    /// `INFO/DP` (single int); used by `--missing DP`.
    dp: Option<i64>,
}

fn common_prefix_len_ci(a: &[u8], b: &[u8]) -> i64 {
    let mut k = 0i64;
    let n = a.len().min(b.len());
    while (k as usize) < n && a[k as usize].eq_ignore_ascii_case(&b[k as usize]) {
        k += 1;
    }
    k
}

fn parse_rec(line: &str) -> Option<Rec> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return None;
    }
    let chrom = f[0].to_owned();
    let pos1: i64 = f[1].parse().ok()?;
    let pos0 = pos1 - 1;
    let reference = f[3].as_bytes();
    let rlen = reference.len() as i64;

    // imin starts at rlen (== REF vs REF common prefix); the REF allele and
    // every non-symbolic ALT pull it down to the minimum shared prefix.
    let mut imin = rlen;
    let alt_field = f[4];
    if alt_field != "." {
        for allele in std::iter::once(f[3]).chain(alt_field.split(',')) {
            if allele.starts_with('<') {
                continue; // ignore symbolic alleles
            }
            let k = common_prefix_len_ci(reference, allele.as_bytes());
            if imin > k {
                imin = k;
            }
        }
    }

    let qual = match f[5] {
        "." | "" => None,
        q => q.parse::<f32>().ok(),
    };
    let dp = f.get(7).and_then(|info| {
        info.split(';')
            .find_map(|kv| kv.strip_prefix("DP="))
            .and_then(|v| v.parse::<i64>().ok())
    });

    Some(Rec {
        line: line.to_owned(),
        chrom,
        pos0,
        rlen,
        imin,
        qual,
        dp,
    })
}

#[derive(PartialEq, Eq)]
enum Status {
    Clean,
    Dirty,
}

/// The `vcfbuf` MARK state machine (overlap + dup paths only).
struct VcfBuf {
    mode: Mark,
    buf: VecDeque<Rec>,
    marks: VecDeque<u8>,
    status: Status,
    overlap_rid: Option<String>,
    overlap_end: i64,
    last_mark: u8,
}

impl VcfBuf {
    fn new(mode: Mark) -> Self {
        VcfBuf {
            mode,
            buf: VecDeque::new(),
            marks: VecDeque::new(),
            status: Status::Clean,
            overlap_rid: None,
            overlap_end: 0,
            last_mark: 0,
        }
    }

    fn push(&mut self, rec: Rec) {
        debug_assert!(self.status != Status::Dirty);
        self.status = Status::Dirty;
        self.buf.push_back(rec);
    }

    /// Port of `mark_overlap_helper_`. Returns the resulting `flush` flag.
    fn mark_overlap_helper(&mut self, flush_all: bool) -> bool {
        if self.status != Status::Dirty {
            return flush_all;
        }
        let mut flush = flush_all;
        self.status = Status::Clean;

        self.marks.push_back(0);
        let last = self.buf.back().unwrap();

        if self.overlap_rid.as_deref() != Some(last.chrom.as_str()) {
            self.overlap_end = 0;
        }
        let mut beg_pos = last.pos0;
        let mut end_pos = last.pos0 + last.rlen - 1;
        let imin = last.imin;

        if beg_pos <= self.overlap_end {
            beg_pos += imin;
            if beg_pos > end_pos {
                end_pos = beg_pos;
            }
        }
        if self.buf.len() == 1 {
            self.overlap_rid = Some(last.chrom.clone());
            self.overlap_end = end_pos;
            return flush;
        }
        if beg_pos <= self.overlap_end {
            if self.overlap_end < end_pos {
                self.overlap_end = end_pos;
            }
            let n = self.marks.len();
            self.marks[n - 1] = 1;
            self.marks[n - 2] = 1;
        } else {
            if self.overlap_end < end_pos {
                self.overlap_end = end_pos;
            }
            flush = true;
        }
        flush
    }

    /// Port of `mark_overlap_can_flush_`.
    fn mark_overlap_can_flush(&mut self, flush_all: bool) -> bool {
        let flush = if self.status == Status::Dirty {
            self.mark_overlap_helper(flush_all)
        } else if self.buf.len() > 1 {
            true
        } else {
            flush_all
        };
        if !flush {
            return false;
        }
        self.last_mark = self.marks.pop_front().unwrap_or(0);
        true
    }

    /// Port of `mark_dup_can_flush_`.
    fn mark_dup_can_flush(&mut self, flush_all: bool) -> bool {
        let mut flush = flush_all;
        if self.status == Status::Dirty {
            self.marks.push_back(0);
            if self.buf.len() == 1 {
                // fall through to flush check (flush == flush_all)
            } else {
                let r1 = &self.buf[self.buf.len() - 1];
                let r2 = &self.buf[self.buf.len() - 2];
                let is_dup = r1.chrom == r2.chrom && r1.pos0 == r2.pos0;
                if is_dup {
                    let n = self.marks.len();
                    self.marks[n - 1] = 1;
                    self.marks[n - 2] = 1;
                } else {
                    flush = true;
                }
            }
        } else if self.buf.len() > 1 {
            flush = true;
        }
        if !flush {
            return false;
        }
        self.last_mark = self.marks.pop_front().unwrap_or(0);
        true
    }

    /// Port of `vcfbuf_flush`. Returns the next record to emit, if any.
    fn flush(&mut self, flush_all: bool) -> Option<Rec> {
        if self.buf.is_empty() {
            return None;
        }
        let can_flush = match self.mode {
            Mark::Overlap => self.mark_overlap_can_flush(flush_all),
            Mark::Dup => self.mark_dup_can_flush(flush_all),
            // MinQual is handled by `min_qual_process`, never the buffer.
            Mark::MinQual => flush_all,
        };
        self.status = Status::Clean;
        if !can_flush {
            return None;
        }
        self.buf.pop_front()
    }
}

/// Reads the input VCF/BCF and returns the processed output. When
/// `text_list` is set (`-O t`), the result is a plain `chr<TAB>pos` list of
/// the surviving sites; otherwise it is the rewritten VCF text.
pub fn run(
    input: &Path,
    mode: Mark,
    mark_tag: Option<&str>,
    reverse: bool,
    text_list: bool,
    missing_qual: f32,
    missing_dp: bool,
) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(process(
        &text,
        mode,
        mark_tag,
        reverse,
        text_list,
        missing_qual,
        missing_dp,
    ))
}

/// Upstream `vcfbuf` `MARK_EXPR` (`min(QUAL)`): within each connected
/// overlap component, repeatedly mark the lowest-QUAL record and drop
/// its overlap edges until the component has no overlaps left. A
/// non-overlapping record (singleton component) is never marked, which
/// is exactly upstream's per-dirty-window behavior.
fn min_qual_process(
    lines: &[&str],
    mark_tag: Option<&str>,
    reverse: bool,
    text_list: bool,
    missing_qual: f32,
    missing_dp: bool,
    out: &mut String,
) {
    let recs: Vec<Rec> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|l| parse_rec(l))
        .collect();
    let n = recs.len();

    // Symmetric overlap adjacency (upstream `records_overlap`: same
    // chrom and the earlier record's end >= the later record's start).
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = (&recs[i], &recs[j]);
            if a.chrom != b.chrom {
                continue;
            }
            let (e, s) = if a.pos0 <= b.pos0 {
                (a.pos0 + a.rlen - 1, b.pos0)
            } else {
                (b.pos0 + b.rlen - 1, a.pos0)
            };
            if e >= s {
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }

    let mut mark = vec![0u8; n];
    let mut seen = vec![false; n];
    for start in 0..n {
        if seen[start] || adj[start].is_empty() {
            seen[start] = true;
            continue;
        }
        // Collect the connected component (BFS).
        let mut comp = Vec::new();
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(v) = stack.pop() {
            comp.push(v);
            for &w in &adj[v] {
                if !seen[w] {
                    seen[w] = true;
                    stack.push(w);
                }
            }
        }
        // Per-component overlap value (upstream `mark_expr` `overlap_t.value`).
        // `--missing DP`: missing QUAL is scaled by coverage using the
        // highest-QUAL record in the component (`mark_expr_missing_*_`);
        // the earliest such record's DP is `max_qual_dp`.
        let (mut max_qual, mut max_qual_dp) = (0.0f32, 0i64);
        if missing_dp {
            let mut by_idx = comp.clone();
            by_idx.sort_unstable();
            for &v in &by_idx {
                if let (Some(q), Some(d)) = (recs[v].qual, recs[v].dp)
                    && max_qual < q
                {
                    max_qual = q;
                    max_qual_dp = d;
                }
            }
        }
        let value = |v: usize| -> f32 {
            if let Some(q) = recs[v].qual {
                return q;
            }
            if missing_dp
                && max_qual_dp != 0
                && let Some(d) = recs[v].dp
            {
                return max_qual * d as f32 / max_qual_dp as f32;
            }
            missing_qual
        };

        // Live edge sets within the component.
        use std::collections::BTreeSet;
        let mut edges: std::collections::HashMap<usize, BTreeSet<usize>> = comp
            .iter()
            .map(|&v| (v, adj[v].iter().copied().collect()))
            .collect();
        let mut nolap: usize = comp.iter().map(|&v| edges[&v].len()).sum::<usize>() / 2;
        // Process members ascending by (value, original index).
        let mut order = comp.clone();
        order.sort_by(|&a, &b| {
            value(a)
                .partial_cmp(&value(b))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.cmp(&b))
        });
        for &oi in &order {
            if nolap == 0 {
                break;
            }
            let neighbours: Vec<usize> = edges[&oi].iter().copied().collect();
            for j in neighbours {
                edges.get_mut(&oi).unwrap().remove(&j);
                edges.get_mut(&j).unwrap().remove(&oi);
                nolap -= 1;
            }
            mark[oi] = 1;
        }
    }

    let mut idx = 0;
    for line in lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if parse_rec(line).is_none() {
            continue;
        }
        let r = &recs[idx];
        let m = mark[idx];
        idx += 1;
        let mut keep = m == 0;
        if reverse {
            keep = !keep;
        }
        if text_list {
            if keep {
                out.push_str(&r.chrom);
                out.push('\t');
                out.push_str(&(r.pos0 + 1).to_string());
                out.push('\n');
            }
            continue;
        }
        if keep {
            out.push_str(&r.line);
            out.push('\n');
        } else if let Some(tag) = mark_tag {
            out.push_str(&set_info_flag(&r.line, tag));
            out.push('\n');
        }
    }
}

fn process(
    text: &str,
    mode: Mark,
    mark_tag: Option<&str>,
    reverse: bool,
    text_list: bool,
    missing_qual: f32,
    missing_dp: bool,
) -> String {
    let lines: Vec<&str> = text.lines().collect();

    let mut out = String::with_capacity(text.len() + 256);

    if mode == Mark::MinQual {
        // Upstream `bcf_hdr_printf` appends the mark-tag INFO line at the
        // *end* of the header (just before `#CHROM`), not grouped after
        // the existing `##INFO` lines.
        if !text_list {
            emit_header_append_last(&lines, mark_tag, &mut out);
        }
        min_qual_process(
            &lines,
            mark_tag,
            reverse,
            text_list,
            missing_qual,
            missing_dp,
            &mut out,
        );
        return out;
    }

    if !text_list {
        emit_header(&lines, mark_tag, &mut out);
    }

    let mut vbuf = VcfBuf::new(mode);
    let emit = |rec: &Rec, last_mark: u8, out: &mut String| {
        let mut keep = last_mark == 0;
        if reverse {
            keep = !keep;
        }
        let mut line = rec.line.clone();
        if !keep {
            match mark_tag {
                None => return, // removed, not emitted
                Some(tag) => line = set_info_flag(&rec.line, tag),
            }
        }
        if text_list {
            out.push_str(&rec.chrom);
            out.push('\t');
            out.push_str(&(rec.pos0 + 1).to_string());
            out.push('\n');
        } else {
            out.push_str(&line);
            out.push('\n');
        }
    };

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let Some(rec) = parse_rec(line) else {
            continue;
        };
        vbuf.push(rec);
        while let Some(r) = vbuf.flush(false) {
            let m = vbuf.last_mark;
            emit(&r, m, &mut out);
        }
    }
    while let Some(r) = vbuf.flush(true) {
        let m = vbuf.last_mark;
        emit(&r, m, &mut out);
    }
    out
}

/// Like [`emit_header`] but appends the mark-tag `##INFO` line as the
/// last header line (immediately before `#CHROM`), matching upstream
/// `bcf_hdr_printf` + `bcf_hdr_sync` append order for `-m 'min(QUAL)'`.
fn emit_header_append_last(lines: &[&str], mark_tag: Option<&str>, out: &mut String) {
    let info_header = mark_tag.map(|tag| {
        format!("##INFO=<ID={tag},Type=Flag,Number=0,Description=\"Marked by +remove-overlaps\">")
    });
    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));
    for (idx, line) in lines.iter().enumerate() {
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("#CHROM") {
            if let Some(h) = &info_header {
                out.push_str(h);
                out.push('\n');
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        out.push_str(line);
        out.push('\n');
        if Some(idx) == fileformat && !has_pass {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
        }
    }
}

/// Emits the header with htslib-style normalization: `##FILTER=<ID=PASS>`
/// injected right after `##fileformat` when absent, and the mark-tag INFO
/// definition placed after the last `##INFO=` line (or just before
/// `#CHROM`), matching upstream `bcf_hdr_printf`.
fn emit_header(lines: &[&str], mark_tag: Option<&str>, out: &mut String) {
    let info_header = mark_tag.map(|tag| {
        format!("##INFO=<ID={tag},Type=Flag,Number=0,Description=\"Marked by +remove-overlaps\">")
    });
    let last_info = lines
        .iter()
        .rposition(|l| l.starts_with("##INFO="))
        .or_else(|| lines.iter().position(|l| l.starts_with("#CHROM")));
    let fileformat = lines.iter().position(|l| l.starts_with("##fileformat="));
    let has_pass = lines.iter().any(|l| l.starts_with("##FILTER=<ID=PASS,"));

    for (idx, line) in lines.iter().enumerate() {
        if !line.starts_with('#') {
            break;
        }
        if line.starts_with("#CHROM") {
            if Some(idx) == last_info
                && let Some(h) = &info_header
            {
                out.push_str(h);
                out.push('\n');
            }
            out.push_str(line);
            out.push('\n');
            continue;
        }
        out.push_str(line);
        out.push('\n');
        if Some(idx) == fileformat && !has_pass {
            out.push_str("##FILTER=<ID=PASS,Description=\"All filters passed\">");
            out.push('\n');
        }
        if line.starts_with("##INFO=")
            && Some(idx) == last_info
            && let Some(h) = &info_header
        {
            out.push_str(h);
            out.push('\n');
        }
    }
}

/// Sets an INFO flag (`bcf_update_info_flag`): appends `;TAG`, or replaces a
/// bare `.` INFO column with `TAG`.
fn set_info_flag(line: &str, tag: &str) -> String {
    let mut f: Vec<&str> = line.split('\t').collect();
    if f.len() < 8 {
        return line.to_owned();
    }
    let info = f[7];
    let new_info = if info == "." || info.is_empty() {
        tag.to_owned()
    } else if info.split(';').any(|kv| kv == tag) {
        info.to_owned()
    } else {
        format!("{info};{tag}")
    };
    f[7] = new_info.as_str();
    f.join("\t")
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
        ".bcftools-rs-remove-overlaps-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.2\n\
##reference=file:///ref.fa\n\
##contig=<ID=1,length=248956422>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t100000\t.\tCC\tG\t.\t.\t.\n\
1\t100001\t.\tC\tG\t.\t.\t.\n\
1\t789241\t.\tC\tG\t.\t.\t.\n\
1\t789242\t.\tC\tG\t.\t.\t.\n\
1\t789242\t.\tC\tA\t.\t.\t.\n\
1\t789243\t.\tC\tCA\t.\t.\t.\n\
1\t789243\t.\tC\tCCA\t.\t.\t.\n\
1\t790000\t.\tC\tG\t.\t.\t.\n\
1\t900000\t.\tC\tG\t.\t.\t.\n";

    fn data_lines(s: &str) -> Vec<&str> {
        s.lines().filter(|l| !l.starts_with('#')).collect()
    }

    // Hand-traced expectations for the truncated VCF above (no 789244/789245
    // run, so the second 789243 record is not pulled into an overlap). Full
    // upstream-fixture parity is covered in tests/plugin_remove_overlaps.rs.

    #[test]
    fn overlap_remove_keeps_non_overlapping() {
        let out = process(VCF, Mark::Overlap, None, false, false, 0.0, false);
        assert_eq!(
            data_lines(&out),
            vec![
                "1\t789241\t.\tC\tG\t.\t.\t.",
                "1\t789243\t.\tC\tCA\t.\t.\t.",
                "1\t789243\t.\tC\tCCA\t.\t.\t.",
                "1\t790000\t.\tC\tG\t.\t.\t.",
                "1\t900000\t.\tC\tG\t.\t.\t.",
            ]
        );
    }

    #[test]
    fn overlap_mark_tags_overlapping_sites() {
        let out = process(
            VCF,
            Mark::Overlap,
            Some("overlap"),
            false,
            false,
            0.0,
            false,
        );
        let d = data_lines(&out);
        assert_eq!(d[0], "1\t100000\t.\tCC\tG\t.\t.\toverlap");
        assert_eq!(d[1], "1\t100001\t.\tC\tG\t.\t.\toverlap");
        assert_eq!(d[2], "1\t789241\t.\tC\tG\t.\t.\t.");
        assert_eq!(d[3], "1\t789242\t.\tC\tG\t.\t.\toverlap");
        assert_eq!(d[4], "1\t789242\t.\tC\tA\t.\t.\toverlap");
        assert_eq!(d[5], "1\t789243\t.\tC\tCA\t.\t.\t.");
        assert_eq!(d[6], "1\t789243\t.\tC\tCCA\t.\t.\t.");
        assert!(out.contains("##FILTER=<ID=PASS,Description=\"All filters passed\">\n"));
        assert!(out.contains(
            "##INFO=<ID=overlap,Type=Flag,Number=0,Description=\"Marked by +remove-overlaps\">\n"
        ));
    }

    #[test]
    fn overlap_reverse_keeps_only_overlapping() {
        let out = process(VCF, Mark::Overlap, None, true, false, 0.0, false);
        let d = data_lines(&out);
        assert_eq!(
            d,
            vec![
                "1\t100000\t.\tCC\tG\t.\t.\t.",
                "1\t100001\t.\tC\tG\t.\t.\t.",
                "1\t789242\t.\tC\tG\t.\t.\t.",
                "1\t789242\t.\tC\tA\t.\t.\t.",
            ]
        );
    }

    #[test]
    fn dup_marks_same_position_only() {
        let out = process(VCF, Mark::Dup, Some("DUP"), false, false, 0.0, false);
        let d = data_lines(&out);
        // Different positions are not duplicates even if spans overlap.
        assert_eq!(d[0], "1\t100000\t.\tCC\tG\t.\t.\t.");
        assert_eq!(d[1], "1\t100001\t.\tC\tG\t.\t.\t.");
        assert_eq!(d[2], "1\t789241\t.\tC\tG\t.\t.\t.");
        assert_eq!(d[3], "1\t789242\t.\tC\tG\t.\t.\tDUP");
        assert_eq!(d[4], "1\t789242\t.\tC\tA\t.\t.\tDUP");
        assert_eq!(d[5], "1\t789243\t.\tC\tCA\t.\t.\tDUP");
        assert_eq!(d[6], "1\t789243\t.\tC\tCCA\t.\t.\tDUP");
    }

    #[test]
    fn text_list_emits_kept_sites() {
        let out = process(VCF, Mark::Overlap, None, false, true, 0.0, false);
        assert_eq!(
            out,
            "1\t789241\n1\t789243\n1\t789243\n1\t790000\n1\t900000\n"
        );
    }
}
