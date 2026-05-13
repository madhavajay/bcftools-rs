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
        func: None,
        alias: "VCF/BCF manipulation",
        help: "",
    },
    Cmd {
        func: Some(commands::head::main),
        alias: "head",
        help: "view VCF/BCF file headers",
    },
    Cmd {
        func: Some(commands::sort::main),
        alias: "sort",
        help: "sort VCF/BCF files",
    },
    Cmd {
        func: Some(commands::view::main),
        alias: "view",
        help: "VCF/BCF conversion, view, subset and filter VCF/BCF files",
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
        "Version: {} (using htslib-rs {})",
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
    writeln!(out, " automatically even when streaming from a pipe.")?;
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
            // Plugin support is deferred; emit a clear unsupported message.
            eprintln!(
                "[E::main] plugin dispatch '+{}' is not yet implemented",
                &arg1[1..]
            );
            return ExitCode::FAILURE;
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
