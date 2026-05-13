//! Pure Rust port of `bcftools`.
//!
//! Each subcommand is a module under [`commands`] exposing
//! `pub fn main(args: &[std::ffi::OsString]) -> std::process::ExitCode`.
//! The CLI crate dispatches on `argv[1]` exactly like the upstream
//! `main.c` program. HTSlib-shaped behavior is delegated to `htslib-rs`.

pub mod commands;
pub mod diagnostics;
pub mod getopt;
pub mod header_version;
pub mod io;
pub mod version;

pub use diagnostics::{error, error_errno};
