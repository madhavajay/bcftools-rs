//! Minimal `bgzip`-compatible helper for the upstream Perl test harness.
//!
//! This is intentionally not advertised as a bcftools subcommand. The harness
//! expects sibling `bgzip` and `tabix` executables, and invokes `bgzip -c` to
//! BGZF-compress raw fixture bytes before indexing them.

use std::ffi::OsString;
use std::fs::File;
use std::io::{self, Read as _, Write as _};
use std::path::Path;
use std::process::ExitCode;

const USAGE: &str = "\n\
Usage: bgzip -c [FILE]\n\
\n\
";

pub fn main(argv: &[OsString]) -> ExitCode {
    match run(argv) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}

fn run(argv: &[OsString]) -> io::Result<()> {
    let mut write_stdout = false;
    let mut input: Option<&Path> = None;

    for arg in argv.iter().skip(1) {
        match arg.to_string_lossy().as_ref() {
            "-c" | "--stdout" => write_stdout = true,
            "-h" | "--help" => {
                print!("{USAGE}");
                return Ok(());
            }
            raw if raw.starts_with('-') => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unsupported bgzip option: {raw}"),
                ));
            }
            _ if input.is_none() => input = Some(Path::new(arg)),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "bgzip expects at most one input path",
                ));
            }
        }
    }

    if !write_stdout {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "this bgzip helper currently supports only -c",
        ));
    }

    let mut data = Vec::new();
    match input {
        Some(path) => File::open(path)?.read_to_end(&mut data)?,
        None => io::stdin().lock().read_to_end(&mut data)?,
    };

    let stdout = io::stdout().lock();
    let mut writer = htslib_rs::bgzf::io::Writer::new(stdout);
    writer.write_all(&data)?;
    let _stdout = writer.finish()?;
    Ok(())
}
