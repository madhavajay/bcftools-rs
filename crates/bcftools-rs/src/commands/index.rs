//! Port of `bcftools index` (upstream `vcfindex.c`).
//!
//! Currently supported:
//!
//! - `bcftools index FILE.bcf` — build a CSI index for a BCF file (default).
//! - `bcftools index FILE.vcf.gz` — build a CSI index for a bgzf-compressed
//!   VCF file (does not rewrite the input; uses
//!   [`htslib_rs::index_compat::build_vcf_csi_from_path_with_min_shift`]).
//! - `-c, --csi` — generate CSI (default).
//! - `-t, --tbi` — generate TBI for VCF.gz (uses
//!   [`htslib_rs::index_compat::build_vcf_tbi_from_path`]).
//! - `-m, --min-shift INT` — CSI min_shift (default 14).
//! - `-o, --output FILE` — write index to a custom path.
//! - `-f, --force` — overwrite existing index.
//! - `-v, --verbosity INT` — verbosity passthrough.
//!
//! Not yet supported (yields an explicit error pointing at the gap):
//!
//! - `-s, --stats`, `-n, --nrecords`, `-a, --all` — index introspection
//!   subcommands; tracked separately.
//! - `--threads INT` — accepted as a no-op for now.

use std::ffi::OsString;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    build_bcf_csi_with_min_shift, build_vcf_csi_from_path_with_min_shift, build_vcf_tbi_from_path,
    write_csi, write_tbi,
};

use crate::diagnostics::fmt_etag;
use crate::getopt::{Getopt, HasArg, LongOpt};
use crate::io::apply_verbosity;

const USAGE: &str = "\n\
About:   Index bgzip compressed VCF/BCF files for random access.\n\
Usage:   bcftools index [options] <in.bcf>|<in.vcf.gz>\n\
\n\
Indexing options:\n\
    -c, --csi                generate CSI-format index for VCF/BCF files [default]\n\
    -f, --force              overwrite index if it already exists\n\
    -m, --min-shift INT      set minimal interval size for CSI indices to 2^INT [14]\n\
    -o, --output FILE        optional output index file name\n\
    -t, --tbi                generate TBI-format index for VCF files\n\
        --threads INT        use multithreading with INT worker threads [0]\n\
    -v, --verbosity INT      verbosity level\n\
\n\
Stats options:\n\
    -a, --all            with --stats, print stats for all contigs even when zero\n\
    -n, --nrecords       print number of records based on existing index file\n\
    -s, --stats          print per contig stats based on existing index file\n\
\n";

const OPT_THREADS: i32 = 9;

/// Subcommand entry point. `argv[0]` is `"index"`.
pub fn main(argv: &[OsString]) -> ExitCode {
    let long_opts = [
        LongOpt::new("all", HasArg::None, b'a' as i32),
        LongOpt::new("csi", HasArg::None, b'c' as i32),
        LongOpt::new("tbi", HasArg::None, b't' as i32),
        LongOpt::new("force", HasArg::None, b'f' as i32),
        LongOpt::new("min-shift", HasArg::Required, b'm' as i32),
        LongOpt::new("stats", HasArg::None, b's' as i32),
        LongOpt::new("nrecords", HasArg::None, b'n' as i32),
        LongOpt::new("threads", HasArg::Required, OPT_THREADS),
        LongOpt::new("verbosity", HasArg::Required, b'v' as i32),
        LongOpt::new("output-file", HasArg::Required, b'o' as i32),
        LongOpt::new("output", HasArg::Required, b'o' as i32),
    ];

    let mut force = false;
    let mut tbi = false;
    let mut min_shift: i32 = 14;
    let mut outfn: Option<String> = None;
    let mut stats_flags: u32 = 0;
    let _n_threads: i32 = 0;

    const STATS_PER_CONTIG: u32 = 1;
    const STATS_ALL_CONTIGS: u32 = 2;
    const STATS_TOTAL: u32 = 4;

    let mut g = Getopt::new("ctfm:snao:v:", &long_opts, argv);
    loop {
        match g.next() {
            Ok(Some(m)) => match m.code {
                v if v == b'v' as i32 => {
                    if apply_verbosity(m.value.as_deref().unwrap_or("0")).is_err() {
                        eprintln!(
                            "Could not parse argument: --verbosity {}",
                            m.value.as_deref().unwrap_or("")
                        );
                        return ExitCode::from(255);
                    }
                }
                v if v == b'c' as i32 => tbi = false,
                v if v == b't' as i32 => {
                    tbi = true;
                    min_shift = 0;
                }
                v if v == b'f' as i32 => force = true,
                v if v == b'm' as i32 => {
                    let raw = m.value.as_deref().unwrap_or("");
                    match raw.parse::<i32>() {
                        Ok(n) => min_shift = n,
                        Err(_) => {
                            eprintln!("Could not parse argument: --min-shift {raw}");
                            return ExitCode::from(255);
                        }
                    }
                }
                v if v == b's' as i32 => stats_flags |= STATS_PER_CONTIG,
                v if v == b'n' as i32 => stats_flags |= STATS_TOTAL,
                v if v == b'a' as i32 => stats_flags |= STATS_ALL_CONTIGS,
                v if v == OPT_THREADS => {
                    // Accepted but currently a no-op.
                }
                v if v == b'o' as i32 => outfn = m.value,
                _ => {
                    eprint!("{USAGE}");
                    return ExitCode::FAILURE;
                }
            },
            Ok(None) => break,
            Err(_) => {
                eprint!("{USAGE}");
                return ExitCode::FAILURE;
            }
        }
    }

    if min_shift == 0 {
        tbi = true;
    }
    if stats_flags > STATS_TOTAL {
        eprintln!(
            "{}",
            fmt_etag(
                "main_vcfindex",
                "expected only one of --stats or --nrecords options"
            )
        );
        return ExitCode::FAILURE;
    }
    if tbi && min_shift > 0 {
        eprintln!(
            "{}",
            fmt_etag(
                "main_vcfindex",
                "min-shift option only expected for CSI indices "
            )
        );
        return ExitCode::FAILURE;
    }
    if !(0..=30).contains(&min_shift) {
        eprintln!(
            "{}",
            fmt_etag(
                "main_vcfindex",
                &format!("expected min_shift in range [0,30] ({min_shift})")
            )
        );
        return ExitCode::FAILURE;
    }

    let positional = g.rest();
    if positional.is_empty() {
        eprint!("{USAGE}");
        return ExitCode::FAILURE;
    }
    let fname = &positional[0];
    let path = Path::new(fname);

    if stats_flags != 0 {
        eprintln!(
            "{}",
            fmt_etag(
                "main_vcfindex",
                "--stats / --nrecords are not yet implemented"
            )
        );
        return ExitCode::FAILURE;
    }

    let idx_path: PathBuf = match outfn {
        Some(o) => PathBuf::from(o),
        None => {
            let mut p = path.as_os_str().to_owned();
            p.push(if tbi { ".tbi" } else { ".csi" });
            PathBuf::from(p)
        }
    };

    if !force && idx_path.exists() {
        let data_meta = fs::metadata(path);
        let idx_meta = fs::metadata(&idx_path);
        if let (Ok(d), Ok(i)) = (&data_meta, &idx_meta)
            && let (Ok(dt), Ok(it)) = (d.modified(), i.modified())
            && dt <= it
        {
            eprintln!(
                "{}",
                fmt_etag(
                    "main_vcfindex",
                    &format!(
                        "the index file exists. Please use '-f' to overwrite {}",
                        idx_path.display()
                    )
                )
            );
            return ExitCode::FAILURE;
        }
    }

    let detected = match format::detect_path(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!(
                "{}",
                fmt_etag(
                    "main_vcfindex",
                    &format!("index: failed to open {}: {e}", path.display())
                )
            );
            return ExitCode::FAILURE;
        }
    };

    if detected.exact == Exact::Bcf {
        if min_shift == 0 {
            eprintln!(
                "{}",
                fmt_etag("main_vcfindex", "BCF requires CSI (min_shift > 0)")
            );
            return ExitCode::FAILURE;
        }
        match build_bcf_csi_with_min_shift(path, min_shift as u8) {
            Ok(idx) => match write_csi(&idx_path, &idx) {
                Ok(()) => ExitCode::SUCCESS,
                Err(e) => {
                    eprintln!(
                        "{}",
                        fmt_etag(
                            "main_vcfindex",
                            &format!("failed to write index {}: {e}", idx_path.display())
                        )
                    );
                    ExitCode::FAILURE
                }
            },
            Err(e) => {
                eprintln!(
                    "{}",
                    fmt_etag(
                        "main_vcfindex",
                        &format!("failed to create index for \"{}\": {e}", path.display())
                    )
                );
                ExitCode::FAILURE
            }
        }
    } else if detected.exact == Exact::Vcf
        && (detected.compression == Compression::Bgzf || detected.compression == Compression::Gzip)
    {
        if tbi {
            match build_vcf_tbi_from_path(path) {
                Ok(idx) => match write_tbi(&idx_path, &idx) {
                    Ok(()) => ExitCode::SUCCESS,
                    Err(e) => {
                        eprintln!(
                            "{}",
                            fmt_etag(
                                "main_vcfindex",
                                &format!("failed to write index {}: {e}", idx_path.display())
                            )
                        );
                        ExitCode::FAILURE
                    }
                },
                Err(e) => {
                    eprintln!(
                        "{}",
                        fmt_etag(
                            "main_vcfindex",
                            &format!("failed to create index for \"{}\": {e}", path.display())
                        )
                    );
                    ExitCode::FAILURE
                }
            }
        } else {
            match build_vcf_csi_from_path_with_min_shift(path, min_shift as u8) {
                Ok(idx) => match write_csi(&idx_path, &idx) {
                    Ok(()) => ExitCode::SUCCESS,
                    Err(e) => {
                        eprintln!(
                            "{}",
                            fmt_etag(
                                "main_vcfindex",
                                &format!("failed to write index {}: {e}", idx_path.display())
                            )
                        );
                        ExitCode::FAILURE
                    }
                },
                Err(e) => {
                    eprintln!(
                        "{}",
                        fmt_etag(
                            "main_vcfindex",
                            &format!("failed to create index for \"{}\": {e}", path.display())
                        )
                    );
                    ExitCode::FAILURE
                }
            }
        }
    } else {
        eprintln!(
            "{}",
            fmt_etag(
                "main_vcfindex",
                &format!(
                    "\"{}\" is in a format that cannot be usefully indexed",
                    path.display()
                )
            )
        );
        ExitCode::FAILURE
    }
}
