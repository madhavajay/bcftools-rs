//! Pure Rust port of `bcftools`.
//!
//! Each subcommand is a module under [`commands`] exposing
//! `pub fn main(args: &[std::ffi::OsString]) -> std::process::ExitCode`.
//! The CLI crate dispatches on `argv[1]` exactly like the upstream
//! `main.c` program. HTSlib-shaped behavior is delegated to `htslib-rs`.

pub mod abuf;
pub mod commands;
pub mod convert;
pub mod diagnostics;
pub mod filter;
pub mod getopt;
pub mod gff;
pub mod gvcf;
pub mod header_version;
pub mod hmm;
pub mod io;
pub mod numerics;
pub mod ploidy;
pub mod reference;
pub mod regidx;
pub mod smpl_ilist;
pub mod synced;
pub mod tsv2vcf;
pub mod vcfbuf;
pub mod version;

pub use diagnostics::{error, error_errno};
