//! Allele atomization helpers modelled after upstream `abuf.c`.
//!
//! Upstream `abuf` is used by `norm --atomize` to split complex alleles into
//! primitive output records and to build the allele-translation table used when
//! projecting INFO/FORMAT values. This module keeps that core record-independent
//! so the VCF command layer can apply it to parsed records.

use std::cmp::Ordering;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AbufMode {
    Split,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomOverlap {
    Identical,
    NonOverlappingRecord,
    Conflicting,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Atom {
    pub ref_allele: String,
    pub alt_allele: String,
    pub original_alt_index: usize,
    pub beg: usize,
    pub end: usize,
    pub prefix_len: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitAtomRow {
    pub beg: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub needs_star_allele: bool,
    pub original_alt_map: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SplitAlleles {
    Unchanged,
    Split(Vec<SplitAtomRow>),
}

pub fn split_alleles(
    ref_allele: &str,
    alt_alleles: &[&str],
    use_star_allele: bool,
) -> SplitAlleles {
    if alt_alleles.is_empty()
        || !is_acgtn(ref_allele)
        || alt_alleles.iter().any(|allele| !is_acgtn(allele))
    {
        return SplitAlleles::Unchanged;
    }

    let mut atoms = Vec::new();
    for (i, alt) in alt_alleles.iter().enumerate() {
        atomize_allele(ref_allele, alt, i + 1, &mut atoms);
    }
    atoms.sort_by(compare_atoms);
    let rows = build_split_rows(&atoms, alt_alleles.len(), use_star_allele);
    SplitAlleles::Split(rows)
}

pub fn atomize_allele(
    ref_allele: &str,
    alt_allele: &str,
    original_alt_index: usize,
    out: &mut Vec<Atom>,
) {
    let ref_chars: Vec<_> = ref_allele.chars().collect();
    let alt_chars: Vec<_> = alt_allele.chars().collect();
    let mut rlen = ref_chars.len();
    let mut alen = alt_chars.len();
    while rlen > 1 && alen > 1 && ref_chars[rlen - 1].eq_ignore_ascii_case(&alt_chars[alen - 1]) {
        rlen -= 1;
        alen -= 1;
    }

    let max_len = rlen.max(alen);
    let mut current: Option<usize> = None;
    for i in 0..max_len {
        let refb = ref_chars.get(i).copied().unwrap_or('-');
        let altb = alt_chars.get(i).copied().unwrap_or('-');
        if !refb.eq_ignore_ascii_case(&altb) {
            if refb == '-' || altb == '-' {
                let idx = current.expect("indel extension requires a current atom");
                if altb != '-' {
                    out[idx].alt_allele.push(altb);
                }
                if refb != '-' {
                    out[idx].ref_allele.push(refb);
                    out[idx].end += 1;
                }
                continue;
            }

            out.push(Atom {
                ref_allele: refb.to_string(),
                alt_allele: altb.to_string(),
                original_alt_index,
                beg: i,
                end: i,
                prefix_len: 0,
            });
            current = Some(out.len() - 1);

            if rlen != alen && (i + 1 >= rlen || i + 1 >= alen) {
                out.push(Atom {
                    ref_allele: refb.to_string(),
                    alt_allele: refb.to_string(),
                    original_alt_index,
                    beg: i,
                    end: i,
                    prefix_len: 1,
                });
                current = Some(out.len() - 1);
            }
            continue;
        }

        if i + 1 >= rlen || i + 1 >= alen {
            out.push(Atom {
                ref_allele: refb.to_string(),
                alt_allele: altb.to_string(),
                original_alt_index,
                beg: i,
                end: i,
                prefix_len: 0,
            });
            current = Some(out.len() - 1);
        }
    }
}

pub fn atoms_inconsistent(a: &Atom, b: &Atom) -> Ordering {
    a.beg
        .cmp(&b.beg)
        .then_with(|| ascii_cmp(&a.ref_allele, &b.ref_allele))
        .then_with(|| ascii_cmp(&a.alt_allele, &b.alt_allele))
}

pub fn atoms_overlap(a: &Atom, b: &Atom) -> AtomOverlap {
    if a.beg != b.beg {
        return AtomOverlap::Conflicting;
    }
    if a.prefix_len != 0 && a.prefix_len >= b.ref_allele.len() {
        return AtomOverlap::NonOverlappingRecord;
    }
    if b.prefix_len != 0 && b.prefix_len >= a.ref_allele.len() {
        return AtomOverlap::NonOverlappingRecord;
    }
    if !a.ref_allele.eq_ignore_ascii_case(&b.ref_allele) {
        return AtomOverlap::Conflicting;
    }
    if a.prefix_len != 0 && a.prefix_len >= b.alt_allele.len() {
        return AtomOverlap::NonOverlappingRecord;
    }
    if b.prefix_len != 0 && b.prefix_len >= a.alt_allele.len() {
        return AtomOverlap::NonOverlappingRecord;
    }
    if !a.alt_allele.eq_ignore_ascii_case(&b.alt_allele) {
        return AtomOverlap::Conflicting;
    }
    AtomOverlap::Identical
}

fn build_split_rows(
    atoms: &[Atom],
    n_original_alts: usize,
    use_star_allele: bool,
) -> Vec<SplitAtomRow> {
    let mut row_atoms: Vec<&Atom> = Vec::new();
    for atom in atoms {
        if row_atoms
            .last()
            .is_some_and(|last| atoms_inconsistent(last, atom) == Ordering::Equal)
        {
            continue;
        }
        row_atoms.push(atom);
    }

    let mut rows = Vec::with_capacity(row_atoms.len());
    for atom in &row_atoms {
        let mut map = vec![0; n_original_alts];
        map[atom.original_alt_index - 1] = 1;
        rows.push(SplitAtomRow {
            beg: atom.beg,
            ref_allele: atom.ref_allele.clone(),
            alt_allele: atom.alt_allele.clone(),
            needs_star_allele: false,
            original_alt_map: map,
        });
    }

    for atom in atoms {
        for (i, out_atom) in row_atoms.iter().enumerate() {
            if std::ptr::eq(*out_atom, atom) {
                continue;
            }
            if atom.beg > out_atom.end {
                continue;
            }
            if atom.end < out_atom.beg {
                break;
            }
            match atoms_overlap(atom, out_atom) {
                AtomOverlap::Identical => rows[i].original_alt_map[atom.original_alt_index - 1] = 1,
                AtomOverlap::NonOverlappingRecord => {
                    rows[i].original_alt_map[atom.original_alt_index - 1] = 1;
                    rows[i].needs_star_allele = true;
                }
                AtomOverlap::Conflicting => {
                    rows[i].original_alt_map[atom.original_alt_index - 1] = 2;
                    rows[i].needs_star_allele = true;
                }
            }
        }
    }

    if !use_star_allele {
        for row in &mut rows {
            row.needs_star_allele = false;
        }
    }
    rows
}

fn compare_atoms(a: &Atom, b: &Atom) -> Ordering {
    atoms_inconsistent(a, b).then_with(|| a.original_alt_index.cmp(&b.original_alt_index))
}

fn is_acgtn(seq: &str) -> bool {
    seq.chars()
        .all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T' | 'N'))
}

fn ascii_cmp(a: &str, b: &str) -> Ordering {
    a.to_ascii_uppercase().cmp(&b.to_ascii_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atomize_splits_complex_substitutions_without_alignment() {
        let mut atoms = Vec::new();
        atomize_allele("GCGT", "GTGA", 1, &mut atoms);
        assert_eq!(
            atoms,
            [
                Atom {
                    ref_allele: "C".into(),
                    alt_allele: "T".into(),
                    original_alt_index: 1,
                    beg: 1,
                    end: 1,
                    prefix_len: 0
                },
                Atom {
                    ref_allele: "T".into(),
                    alt_allele: "A".into(),
                    original_alt_index: 1,
                    beg: 3,
                    end: 3,
                    prefix_len: 0
                }
            ]
        );
    }

    #[test]
    fn atomize_tracks_indel_prefix_rows_like_upstream_abuf() {
        let mut atoms = Vec::new();
        atomize_allele("C", "CA", 1, &mut atoms);
        assert_eq!(
            atoms,
            [Atom {
                ref_allele: "C".into(),
                alt_allele: "CA".into(),
                original_alt_index: 1,
                beg: 0,
                end: 0,
                prefix_len: 0
            }]
        );

        atoms.clear();
        atomize_allele("C", "GGG", 1, &mut atoms);
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[1].prefix_len, 1);
    }

    #[test]
    fn split_rows_deduplicate_identical_atoms_and_map_source_alts() {
        let SplitAlleles::Split(rows) = split_alleles("CC", &["TC", "TC"], true) else {
            panic!("expected split");
        };
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].beg, 0);
        assert_eq!(rows[0].ref_allele, "C");
        assert_eq!(rows[0].alt_allele, "T");
        assert_eq!(rows[0].original_alt_map, [1, 1]);
        assert!(!rows[0].needs_star_allele);
    }

    #[test]
    fn split_rows_mark_conflicting_overlaps_with_star_allele() {
        let SplitAlleles::Split(rows) = split_alleles("A", &["AT", "C"], true) else {
            panic!("expected split");
        };
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].original_alt_map, [1, 2]);
        assert!(rows[0].needs_star_allele);
        assert_eq!(rows[1].original_alt_map, [2, 1]);
        assert!(rows[1].needs_star_allele);
    }

    #[test]
    fn split_rows_can_disable_star_allele() {
        let SplitAlleles::Split(rows) = split_alleles("A", &["AT", "C"], false) else {
            panic!("expected split");
        };
        assert!(rows.iter().all(|row| !row.needs_star_allele));
    }

    #[test]
    fn non_acgtn_alleles_are_left_unchanged() {
        assert_eq!(
            split_alleles("A", &["<DEL>"], true),
            SplitAlleles::Unchanged
        );
    }
}
