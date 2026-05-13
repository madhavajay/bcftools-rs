//! Sample-list helpers mirroring bcftools `smpl_ilist.c`.
//!
//! The C helper is parameterized by `bcf_hdr_t`; this Rust version accepts a
//! sample-name slice so callers can source names from whichever header wrapper
//! they are using.

use std::{collections::HashMap, fs, io, path::Path};

pub const SMPL_NONE: u32 = 0;
pub const SMPL_STRICT: u32 = 1;
pub const SMPL_SINGLE: u32 = 2;
pub const SMPL_PAIR1: u32 = 4;
pub const SMPL_PAIR2: u32 = 8;
pub const SMPL_VERBOSE: u32 = 16;
pub const SMPL_REORDER: u32 = 32;

/// Parsed sample index list.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleIndexList {
    /// Indexes into the input header sample list.
    pub idx: Vec<usize>,
    /// Optional paired sample names, aligned with `idx`.
    pub pair: Option<Vec<Option<String>>>,
}

impl SampleIndexList {
    pub fn len(&self) -> usize {
        self.idx.len()
    }

    pub fn is_empty(&self) -> bool {
        self.idx.is_empty()
    }
}

/// Parses `--samples` or `--samples-file` style input.
pub fn init(
    header_samples: &[impl AsRef<str>],
    sample_list: Option<&str>,
    is_file: bool,
    mut flags: u32,
) -> io::Result<SampleIndexList> {
    let names: Vec<&str> = header_samples.iter().map(AsRef::as_ref).collect();
    if sample_list.is_none() {
        return Ok(SampleIndexList {
            idx: (0..names.len()).collect(),
            pair: None,
        });
    }

    let sample_list = sample_list.unwrap();
    let negate = sample_list.starts_with('^');
    let raw_list = if negate {
        &sample_list['^'.len_utf8()..]
    } else {
        sample_list
    };
    let entries = read_list(raw_list, is_file)?;
    if negate && flags & SMPL_REORDER != 0 {
        flags &= !SMPL_REORDER;
    }

    let name_to_idx: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(idx, name)| (*name, idx))
        .collect();

    if flags & SMPL_REORDER != 0 {
        let mut idx = Vec::new();
        for entry in &entries {
            let (sample1, sample2) = split_pair(entry);
            let sample_name = if flags & SMPL_PAIR2 != 0 {
                sample2.unwrap_or(sample1)
            } else {
                sample1
            };
            if let Some(sample_idx) = name_to_idx.get(sample_name) {
                idx.push(*sample_idx);
            } else if flags & SMPL_STRICT != 0 {
                return Err(no_such_sample(sample_name));
            }
        }
        return Ok(SampleIndexList { idx, pair: None });
    }

    let mut selected = vec![false; names.len()];
    let mut pair_by_header_idx: Option<Vec<Option<String>>> = None;
    let mut found = 0usize;

    for entry in &entries {
        let (sample1, sample2) = split_pair(entry);
        let sample_name = if flags & SMPL_PAIR2 != 0 {
            sample2.unwrap_or(sample1)
        } else {
            sample1
        };
        let Some(&sample_idx) = name_to_idx.get(sample_name) else {
            if flags & SMPL_STRICT != 0 {
                return Err(no_such_sample(sample_name));
            }
            continue;
        };

        if !selected[sample_idx] {
            found += 1;
        }
        selected[sample_idx] = true;

        if let Some(other) = sample2 {
            let pairs = pair_by_header_idx.get_or_insert_with(|| vec![None; names.len()]);
            if flags & SMPL_PAIR2 != 0 {
                pairs[sample_idx] = Some(sample1.to_string());
            } else if flags & SMPL_PAIR1 != 0 {
                pairs[sample_idx] = Some(other.to_string());
            }
        }
    }

    let mut idx = Vec::new();
    let mut pair = pair_by_header_idx.as_ref().map(|_| Vec::new());
    for sample_idx in 0..names.len() {
        let keep = if negate {
            !selected[sample_idx]
        } else {
            selected[sample_idx]
        };
        if keep {
            idx.push(sample_idx);
            if let (Some(src), Some(dst)) = (&pair_by_header_idx, pair.as_mut()) {
                dst.push(src[sample_idx].clone());
            }
        }
    }

    if negate {
        debug_assert_eq!(idx.len(), names.len().saturating_sub(found));
    }
    Ok(SampleIndexList { idx, pair })
}

/// Maps samples in `header_a` to their indexes in `header_b`.
pub fn map(
    header_a: &[impl AsRef<str>],
    header_b: &[impl AsRef<str>],
    flags: u32,
) -> io::Result<SampleIndexList> {
    if flags & SMPL_STRICT != 0 && header_a.len() != header_b.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "Different number of samples: {} vs {}",
                header_a.len(),
                header_b.len()
            ),
        ));
    }

    let b_names: Vec<&str> = header_b.iter().map(AsRef::as_ref).collect();
    let b_lookup: HashMap<&str, usize> = b_names
        .iter()
        .enumerate()
        .map(|(idx, name)| (*name, idx))
        .collect();

    let mut idx = Vec::with_capacity(header_a.len());
    for name in header_a.iter().map(AsRef::as_ref) {
        if let Some(sample_idx) = b_lookup.get(name) {
            idx.push(*sample_idx);
        } else if flags & SMPL_STRICT != 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("The sample {name} is not present in the second file"),
            ));
        }
    }

    Ok(SampleIndexList { idx, pair: None })
}

fn read_list(raw: &str, is_file: bool) -> io::Result<Vec<String>> {
    if is_file {
        let text = fs::read_to_string(Path::new(raw))?;
        Ok(text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(ToOwned::to_owned)
            .collect())
    } else {
        Ok(raw
            .split(',')
            .map(str::trim)
            .filter(|item| !item.is_empty())
            .map(ToOwned::to_owned)
            .collect())
    }
}

fn split_pair(entry: &str) -> (&str, Option<&str>) {
    let bytes = entry.as_bytes();
    for (idx, byte) in bytes.iter().enumerate() {
        if !byte.is_ascii_whitespace() || escaped_space(bytes, idx) {
            continue;
        }
        let left = &entry[..idx];
        let right = entry[idx..].trim_start();
        return (left, (!right.is_empty()).then_some(right));
    }
    (entry, None)
}

fn escaped_space(bytes: &[u8], space_idx: usize) -> bool {
    let mut n = 0usize;
    let mut idx = space_idx;
    while idx > 0 {
        idx -= 1;
        if bytes[idx] != b'\\' {
            break;
        }
        n += 1;
    }
    n % 2 == 1
}

fn no_such_sample(name: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("No such sample: \"{name}\""),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn none_selects_all_samples_in_header_order() {
        let list = init(&["a", "b", "c"], None, false, SMPL_NONE).unwrap();
        assert_eq!(list.idx, vec![0, 1, 2]);
        assert_eq!(list.pair, None);
    }

    #[test]
    fn include_preserves_header_order_by_default() {
        let list = init(&["a", "b", "c"], Some("c,a"), false, SMPL_NONE).unwrap();
        assert_eq!(list.idx, vec![0, 2]);
    }

    #[test]
    fn reorder_preserves_requested_order() {
        let list = init(&["a", "b", "c"], Some("c,a"), false, SMPL_REORDER).unwrap();
        assert_eq!(list.idx, vec![2, 0]);
    }

    #[test]
    fn exclude_keeps_unlisted_samples() {
        let list = init(&["a", "b", "c"], Some("^b"), false, SMPL_NONE).unwrap();
        assert_eq!(list.idx, vec![0, 2]);
    }

    #[test]
    fn strict_rejects_unknown_samples() {
        let err = init(&["a"], Some("missing"), false, SMPL_STRICT).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(err.to_string().contains("No such sample"));
    }

    #[test]
    fn pairs_can_match_first_or_second_column() {
        let pair1 = init(&["a", "b"], Some("a x,b y"), false, SMPL_PAIR1).unwrap();
        assert_eq!(pair1.idx, vec![0, 1]);
        assert_eq!(
            pair1.pair,
            Some(vec![Some("x".to_string()), Some("y".to_string())])
        );

        let pair2 = init(&["x", "y"], Some("a x,b y"), false, SMPL_PAIR2).unwrap();
        assert_eq!(pair2.idx, vec![0, 1]);
        assert_eq!(
            pair2.pair,
            Some(vec![Some("a".to_string()), Some("b".to_string())])
        );
    }

    #[test]
    fn escaped_space_does_not_split_pair() {
        let list = init(&["a\\ b", "c"], Some("a\\ b,c"), false, SMPL_NONE).unwrap();
        assert_eq!(list.idx, vec![0, 1]);
    }

    #[test]
    fn file_input_reads_one_entry_per_line() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "c").unwrap();
        writeln!(file, "a").unwrap();
        let list = init(
            &["a", "b", "c"],
            Some(file.path().to_str().unwrap()),
            true,
            SMPL_REORDER,
        )
        .unwrap();
        assert_eq!(list.idx, vec![2, 0]);
    }

    #[test]
    fn map_projects_first_header_into_second() {
        let list = map(&["b", "a"], &["a", "b", "c"], SMPL_NONE).unwrap();
        assert_eq!(list.idx, vec![1, 0]);

        let err = map(&["b", "a"], &["a", "b", "c"], SMPL_STRICT).unwrap_err();
        assert!(err.to_string().contains("Different number of samples"));
    }
}
