//! FASTA reference helpers shared by bcftools subcommands.
//!
//! This is a bcftools-facing facade over `htslib-rs::faidx_compat`, mirroring
//! the common `fai_load` / `faidx_fetch_seq` workflow without exposing callers
//! to the lower-level index plumbing.

use std::{
    fs::File,
    io::{self, BufReader},
    path::{Path, PathBuf},
};

pub use htslib_rs::faidx_compat::Index as FastaIndex;

/// Indexed FASTA reference.
#[derive(Clone, Debug)]
pub struct FastaReference {
    path: PathBuf,
    index: FastaIndex,
}

impl FastaReference {
    /// Opens a FASTA reference and its `.fai` if present, otherwise builds an
    /// in-memory index from the FASTA contents.
    pub fn open<P>(path: P) -> io::Result<Self>
    where
        P: AsRef<Path>,
    {
        let path = path.as_ref();
        let index = match File::open(fai_path(path)) {
            Ok(file) => htslib_rs::faidx_compat::read_index(BufReader::new(file))?,
            Err(error) if error.kind() == io::ErrorKind::NotFound => {
                let file = File::open(path)?;
                htslib_rs::faidx_compat::build_index(BufReader::new(file))?
            }
            Err(error) => return Err(error),
        };

        Ok(Self {
            path: path.to_path_buf(),
            index,
        })
    }

    /// Returns the loaded FASTA index.
    pub fn index(&self) -> &FastaIndex {
        &self.index
    }

    /// Returns whether the reference contains `name`.
    pub fn has_sequence(&self, name: &str) -> bool {
        htslib_rs::faidx_compat::has_sequence(&self.index, name)
    }

    /// Returns the sequence length for `name`.
    pub fn sequence_len(&self, name: &str) -> Option<u64> {
        htslib_rs::faidx_compat::sequence_len(&self.index, name)
    }

    /// Fetches a sequence for an HTSlib-style region string such as
    /// `chr1:10-20`.
    pub fn fetch_region(&self, region: &str) -> io::Result<Vec<u8>> {
        let file = File::open(&self.path)?;
        let mut reader = BufReader::new(file);
        htslib_rs::faidx_compat::fetch_region_sequence(&mut reader, &self.index, region)
    }
}

/// Returns the conventional FASTA index path for `path`.
pub fn fai_path(path: &Path) -> PathBuf {
    let mut raw = path.as_os_str().to_os_string();
    raw.push(".fai");
    PathBuf::from(raw)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as _;

    fn write_fasta() -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("ref.fa");
        std::fs::write(&path, b">chr1\nACGTACGT\n>chr2 description\nNNNN\n").unwrap();
        (dir, path)
    }

    #[test]
    fn opens_fasta_and_builds_index_when_fai_is_absent() {
        let (_dir, path) = write_fasta();
        let reference = FastaReference::open(&path).unwrap();

        assert!(reference.has_sequence("chr1"));
        assert_eq!(reference.sequence_len("chr1"), Some(8));
        assert_eq!(reference.fetch_region("chr1:2-5").unwrap(), b"CGTA");
    }

    #[test]
    fn prefers_existing_fai_index() {
        let (_dir, path) = write_fasta();
        let index =
            htslib_rs::faidx_compat::build_index(BufReader::new(File::open(&path).unwrap()))
                .unwrap();
        let mut out = File::create(fai_path(&path)).unwrap();
        htslib_rs::faidx_compat::write_index(&mut out, &index).unwrap();
        out.flush().unwrap();

        let reference = FastaReference::open(&path).unwrap();
        assert_eq!(reference.sequence_len("chr2"), Some(4));
        assert_eq!(reference.fetch_region("chr2:1-4").unwrap(), b"NNNN");
    }

    #[test]
    fn missing_region_errors() {
        let (_dir, path) = write_fasta();
        let reference = FastaReference::open(&path).unwrap();
        let err = reference.fetch_region("missing:1-1").unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::NotFound);
    }
}
