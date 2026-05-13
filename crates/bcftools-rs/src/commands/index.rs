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
//! - `-s, --stats`, `-n, --nrecords`, `-a, --all` — index introspection
//!   based on existing CSI/TBI metadata.
//!
//! Not yet supported (yields an explicit error pointing at the gap):
//!
//! - `--threads INT` — accepted as a no-op for now.

use std::ffi::OsString;
use std::fs;
use std::io::{self, BufReader, Read as _};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::index_compat::{
    IndexFormat, build_bcf_csi_with_min_shift, build_vcf_csi_from_path_with_min_shift,
    build_vcf_tbi_from_path, locate_associated_index, read_csi, read_tbi, write_csi, write_tbi,
};
use htslib_rs::{csi::BinningIndex, vcf};

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
                    let raw = m.value.as_deref().unwrap_or("");
                    if raw.parse::<i32>().is_err() {
                        eprintln!("Could not parse argument: --threads {raw}");
                        return ExitCode::from(255);
                    }
                    // Accepted but currently a no-op; indexing is single-threaded
                    // until htslib-rs exposes BGZF worker pools for index builds.
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
        let Some(idx_path) = outfn.as_deref().map(PathBuf::from) else {
            eprint!("{USAGE}");
            return ExitCode::FAILURE;
        };
        return match index_stdin(&idx_path, tbi, min_shift) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("main_vcfindex", &format!("{e}")));
                ExitCode::FAILURE
            }
        };
    }
    let fname = &positional[0];
    let path = Path::new(fname);

    if stats_flags != 0 {
        return match print_index_stats(path, stats_flags) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("{}", fmt_etag("vcf_index_stats", &format!("{e}")));
                ExitCode::FAILURE
            }
        };
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

fn index_stdin(idx_path: &Path, tbi: bool, min_shift: i32) -> io::Result<()> {
    let mut data = Vec::new();
    io::stdin().lock().read_to_end(&mut data)?;
    if data.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "no input data on stdin",
        ));
    }

    let tmp_path = stdin_tmp_path(idx_path);
    fs::write(&tmp_path, &data)?;
    let result = build_index_to_path(&tmp_path, idx_path, tbi, min_shift);
    let _ = fs::remove_file(&tmp_path);
    result
}

fn stdin_tmp_path(idx_path: &Path) -> PathBuf {
    let parent = idx_path.parent().unwrap_or_else(|| Path::new("."));
    let stem = idx_path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("stdin");
    parent.join(format!(".{stem}.stdin.{}.tmp", std::process::id()))
}

fn build_index_to_path(path: &Path, idx_path: &Path, tbi: bool, min_shift: i32) -> io::Result<()> {
    let detected = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;

    if detected.exact == Exact::Bcf {
        if min_shift == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "BCF requires CSI (min_shift > 0)",
            ));
        }
        let idx = build_bcf_csi_with_min_shift(path, min_shift as u8)?;
        write_csi(idx_path, &idx)
    } else if detected.exact == Exact::Vcf
        && (detected.compression == Compression::Bgzf || detected.compression == Compression::Gzip)
    {
        if tbi {
            let idx = build_vcf_tbi_from_path(path)?;
            write_tbi(idx_path, &idx)
        } else {
            let idx = build_vcf_csi_from_path_with_min_shift(path, min_shift as u8)?;
            write_csi(idx_path, &idx)
        }
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "\"{}\" is in a format that cannot be usefully indexed",
                path.display()
            ),
        ))
    }
}

#[derive(Debug)]
struct ContigStat {
    name: Option<String>,
    length: Option<usize>,
    records: u64,
}

fn print_index_stats(path: &Path, stats_flags: u32) -> io::Result<()> {
    const STATS_PER_CONTIG: u32 = 1;
    const STATS_ALL_CONTIGS: u32 = 2;
    const STATS_TOTAL: u32 = 4;

    let header = read_variant_header_for_stats(path).ok();
    let index = read_index_for_stats(path)?;
    let contigs = collect_contig_stats(index.as_ref(), header.as_ref());
    let total: u64 = contigs.iter().map(|stat| stat.records).sum();

    if stats_flags & STATS_PER_CONTIG != 0 {
        for stat in &contigs {
            if stat.records == 0 && stats_flags & STATS_ALL_CONTIGS == 0 {
                continue;
            }
            let name = stat.name.as_deref().unwrap_or("n/a");
            let length = stat
                .length
                .map(|n| n.to_string())
                .unwrap_or_else(|| ".".to_string());
            println!("{name}\t{length}\t{}", stat.records);
        }
    }

    if stats_flags & STATS_TOTAL != 0 {
        println!("{total}");
    }

    Ok(())
}

fn read_index_for_stats(path: &Path) -> io::Result<Box<dyn BinningIndex>> {
    if path_has_suffix(path, ".tbi") {
        return read_tbi(path).map(|idx| Box::new(idx) as Box<dyn BinningIndex>);
    }
    if path_has_suffix(path, ".csi") {
        return read_csi(path).map(|idx| Box::new(idx) as Box<dyn BinningIndex>);
    }

    let detected = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if detected.exact == Exact::Bcf {
        let located = locate_associated_index(path, IndexFormat::Csi)
            .ok_or_else(|| missing_index_error(path))?;
        read_csi(located.path).map(|idx| Box::new(idx) as Box<dyn BinningIndex>)
    } else {
        let located = locate_associated_index(path, IndexFormat::Tbi)
            .ok_or_else(|| missing_index_error(path))?;
        match located.format {
            IndexFormat::Tbi => read_tbi(located.path).map(|idx| Box::new(idx) as _),
            IndexFormat::Csi => read_csi(located.path).map(|idx| Box::new(idx) as _),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "unsupported index format for VCF/BCF stats",
            )),
        }
    }
}

fn read_variant_header_for_stats(path: &Path) -> io::Result<vcf::Header> {
    let data_path = data_path_for_index_path(path);
    let detected = format::detect_path(&data_path).map_err(|e| io::Error::other(e.to_string()))?;
    if detected.exact == Exact::Bcf {
        htslib_rs::variant_io_compat::read_bcf_header_from_path(&data_path)
    } else if detected.exact == Exact::Vcf
        && (detected.compression == Compression::Bgzf || detected.compression == Compression::Gzip)
    {
        let file = fs::File::open(&data_path)?;
        let reader = flate2::read::MultiGzDecoder::new(file);
        htslib_rs::variant_io_compat::read_vcf_header(BufReader::new(reader))
    } else {
        htslib_rs::variant_io_compat::read_vcf_header_from_path(&data_path)
    }
}

fn collect_contig_stats(index: &dyn BinningIndex, header: Option<&vcf::Header>) -> Vec<ContigStat> {
    let reference_sequences: Vec<_> = index.reference_sequences().collect();
    let index_names: Vec<String> = index
        .header()
        .map(|h| {
            h.reference_sequence_names()
                .iter()
                .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
                .collect()
        })
        .unwrap_or_default();
    let can_use_header_order = index_names.is_empty()
        && header
            .map(|h| h.contigs().len() == reference_sequences.len())
            .unwrap_or(false);

    reference_sequences
        .into_iter()
        .enumerate()
        .map(|(tid, reference_sequence)| {
            let indexed_name = index_names.get(tid).cloned();
            let (name, length) = if let (Some(h), Some(indexed_name)) = (header, indexed_name) {
                let length = h
                    .contigs()
                    .get(&indexed_name)
                    .and_then(|contig| contig.length());
                (Some(indexed_name), length)
            } else if can_use_header_order {
                header
                    .and_then(|h| h.contigs().get_index(tid))
                    .map(|(name, contig)| (Some(name.clone()), contig.length()))
                    .unwrap_or((None, None))
            } else {
                (None, None)
            };
            ContigStat {
                name,
                length,
                records: reference_sequence
                    .metadata()
                    .map(|m| m.mapped_record_count())
                    .unwrap_or(0),
            }
        })
        .collect()
}

fn data_path_for_index_path(path: &Path) -> PathBuf {
    let raw = path.as_os_str().to_string_lossy();
    let data = raw
        .strip_suffix(".tbi")
        .or_else(|| raw.strip_suffix(".csi"))
        .unwrap_or(&raw);
    PathBuf::from(data)
}

fn path_has_suffix(path: &Path, suffix: &str) -> bool {
    path.as_os_str()
        .to_string_lossy()
        .to_ascii_lowercase()
        .ends_with(suffix)
}

fn missing_index_error(path: &Path) -> io::Error {
    io::Error::new(
        io::ErrorKind::NotFound,
        format!(
            "No index file could be found for '{}'. Use 'bcftools index' to create one",
            path.display()
        ),
    )
}
