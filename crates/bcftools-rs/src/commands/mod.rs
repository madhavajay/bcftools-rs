//! Subcommand implementations.
//!
//! Each module exposes:
//!
//! ```ignore
//! pub fn main(args: &[std::ffi::OsString]) -> std::process::ExitCode
//! ```
//!
//! `args` is the full argv slice the subcommand receives, where `args[0]` is
//! the subcommand name (matching upstream's `main_<name>(argc, argv)` calling
//! convention where `argv[0]` is e.g. `"view"`).

pub mod head;
pub mod index;
pub mod query;
pub mod reheader;
pub mod sort;
pub mod tabix;
pub mod view;
