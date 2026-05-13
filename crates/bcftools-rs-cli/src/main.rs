//! `bcftools` binary entry point. Port of upstream `main.c`.
//!
//! Dispatches to per-subcommand `bcftools_rs::commands::<name>::main(argv)`.
//! `+name` is a shorthand for `plugin name` exactly like upstream
//! (`main.c:289-296`). The version/help dispatch shapes match upstream
//! verbatim where feasible.

use std::env;
use std::ffi::OsString;
use std::io::Write;
use std::process::ExitCode;

use bcftools_rs::commands;
use bcftools_rs::version::{version_block, version_only_string};

fn unsupported(argv: &[OsString]) -> ExitCode {
    let name = argv
        .first()
        .map(|s| s.to_string_lossy())
        .unwrap_or_else(|| "<unknown>".into());
    eprintln!("[E::main] command '{name}' is not yet implemented");
    ExitCode::FAILURE
}

struct Cmd {
    func: Option<fn(&[OsString]) -> ExitCode>,
    alias: &'static str,
    help: &'static str,
}

const CMDS: &[Cmd] = &[
    Cmd {
        func: None,
        alias: "Indexing",
        help: "",
    },
    Cmd {
        func: Some(commands::index::main),
        alias: "index",
        help: "index VCF/BCF files",
    },
    Cmd {
        func: Some(commands::tabix::main),
        alias: "tabix",
        help: "-tabix for BGZF'd BED, GFF, SAM, VCF and more",
    },
    Cmd {
        func: None,
        alias: "VCF/BCF manipulation",
        help: "",
    },
    Cmd {
        func: Some(unsupported),
        alias: "annotate",
        help: "annotate and edit VCF/BCF files",
    },
    Cmd {
        func: Some(unsupported),
        alias: "concat",
        help: "concatenate VCF/BCF files from the same set of samples",
    },
    Cmd {
        func: Some(unsupported),
        alias: "convert",
        help: "convert VCF/BCF files to different formats and back",
    },
    Cmd {
        func: Some(commands::head::main),
        alias: "head",
        help: "view VCF/BCF file headers",
    },
    Cmd {
        func: Some(unsupported),
        alias: "isec",
        help: "intersections of VCF/BCF files",
    },
    Cmd {
        func: Some(unsupported),
        alias: "merge",
        help: "merge VCF/BCF files files from non-overlapping sample sets",
    },
    Cmd {
        func: Some(unsupported),
        alias: "norm",
        help: "left-align and normalize indels",
    },
    Cmd {
        func: Some(commands::query::main),
        alias: "query",
        help: "transform VCF/BCF into user-defined formats",
    },
    Cmd {
        func: Some(commands::reheader::main),
        alias: "reheader",
        help: "modify VCF/BCF header, change sample names",
    },
    Cmd {
        func: Some(commands::sort::main),
        alias: "sort",
        help: "sort VCF/BCF file",
    },
    Cmd {
        func: Some(commands::view::main),
        alias: "view",
        help: "VCF/BCF conversion, view, subset and filter VCF/BCF files",
    },
    Cmd {
        func: None,
        alias: "VCF/BCF analysis",
        help: "",
    },
    Cmd {
        func: Some(unsupported),
        alias: "call",
        help: "SNP/indel calling",
    },
    Cmd {
        func: Some(unsupported),
        alias: "consensus",
        help: "create consensus sequence by applying VCF variants",
    },
    Cmd {
        func: Some(unsupported),
        alias: "cnv",
        help: "HMM CNV calling",
    },
    Cmd {
        func: Some(unsupported),
        alias: "csq",
        help: "call variation consequences",
    },
    Cmd {
        func: Some(unsupported),
        alias: "filter",
        help: "filter VCF/BCF files using fixed thresholds",
    },
    Cmd {
        func: Some(unsupported),
        alias: "gtcheck",
        help: "check sample concordance, detect sample swaps and contamination",
    },
    Cmd {
        func: Some(unsupported),
        alias: "mpileup",
        help: "multi-way pileup producing genotype likelihoods",
    },
    Cmd {
        func: Some(unsupported),
        alias: "roh",
        help: "identify runs of autozygosity (HMM)",
    },
    Cmd {
        func: Some(unsupported),
        alias: "stats",
        help: "produce VCF/BCF stats",
    },
    Cmd {
        func: Some(unsupported),
        alias: "som",
        help: "-filter using Self-Organized Maps (experimental)",
    },
    Cmd {
        func: None,
        alias: "Plugins",
        help: "",
    },
    Cmd {
        func: Some(unsupported),
        alias: "plugin",
        help: "user-defined plugins",
    },
];

fn usage<W: Write>(out: &mut W) -> std::io::Result<()> {
    writeln!(out)?;
    writeln!(
        out,
        "Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)"
    )?;
    writeln!(
        out,
        "Version: {} (using htslib {})",
        bcftools_rs::version::BCFTOOLS_VERSION,
        bcftools_rs::version::HTSLIB_RS_VERSION
    )?;
    writeln!(out)?;
    writeln!(
        out,
        "Usage:   bcftools [--version|--version-only] [--help] <command> <argument>"
    )?;
    writeln!(out)?;
    writeln!(out, "Commands:")?;
    let mut last_section: Option<&str> = None;
    for c in CMDS {
        if c.func.is_none() {
            last_section = Some(c.alias);
            writeln!(out)?;
            writeln!(out, " -- {}", c.alias)?;
            continue;
        }
        if last_section.is_some() && !c.help.starts_with('-') {
            writeln!(out, "    {:<12} {}", c.alias, c.help)?;
        }
    }
    writeln!(out)?;
    writeln!(
        out,
        " Most commands accept VCF, bgzipped VCF, and BCF with the file type detected"
    )?;
    writeln!(
        out,
        " automatically even when streaming from a pipe. Indexed VCF and BCF will work"
    )?;
    writeln!(
        out,
        " in all situations. Un-indexed VCF and BCF and streams will work in most but"
    )?;
    writeln!(out, " not all situations.")?;
    writeln!(out)?;
    Ok(())
}

fn main() -> ExitCode {
    let argv: Vec<OsString> = env::args_os().collect();
    if argv.len() < 2 {
        let mut err = std::io::stderr().lock();
        let _ = usage(&mut err);
        return ExitCode::FAILURE;
    }

    let arg1 = argv[1].to_string_lossy().into_owned();
    match arg1.as_str() {
        "version" | "--version" | "-v" => {
            print!("{}", version_block());
            return ExitCode::SUCCESS;
        }
        "--version-only" => {
            println!("{}", version_only_string());
            return ExitCode::SUCCESS;
        }
        "help" | "--help" | "-h" => {
            if argv.len() == 2 {
                let mut out = std::io::stdout().lock();
                let _ = usage(&mut out);
                return ExitCode::SUCCESS;
            }
            // `bcftools help COMMAND [...]` → `bcftools COMMAND` (subcommand
            // prints its usage when called without args). We pass an
            // argv slice consisting of just the subcommand name.
            let cmd = argv[2].clone();
            let sub_argv = vec![cmd];
            return dispatch(&sub_argv);
        }
        _ if arg1.starts_with('+') => {
            // `bcftools +name ...` → `bcftools plugin name ...`
            let mut sub_argv = Vec::with_capacity(argv.len() + 1);
            sub_argv.push(OsString::from("plugin"));
            sub_argv.push(OsString::from(&arg1[1..]));
            sub_argv.extend(argv.iter().skip(2).cloned());
            return dispatch(&sub_argv);
        }
        _ => {}
    }

    dispatch(&argv[1..])
}

fn dispatch(sub_argv: &[OsString]) -> ExitCode {
    let name = sub_argv[0].to_string_lossy().into_owned();
    for c in CMDS {
        if let Some(func) = c.func
            && c.alias == name
        {
            return func(sub_argv);
        }
    }
    eprintln!("[E::main] unrecognized command '{name}'");
    ExitCode::FAILURE
}
