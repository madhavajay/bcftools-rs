//! Static in-process plugin registry.
//!
//! The upstream command discovers shared objects through `BCFTOOLS_PLUGINS`.
//! This Rust port builds plugins into the binary, so this module provides the
//! registry/listing surface first while individual plugin algorithms are ported
//! behind the same names.

use std::ffi::OsString;
use std::io::{self, Write};
use std::path::Path;
use std::process::ExitCode;

use crate::diagnostics::fmt_etag;

const USAGE: &str = "\n\
About:   Run user defined plugin\n\
Usage:   bcftools plugin <name> [OPTIONS] <file> [-- PLUGIN_OPTIONS]\n\
         bcftools +name [OPTIONS] <file>  [-- PLUGIN_OPTIONS]\n\
\n\
VCF input options:\n\
   -e, --exclude EXPR             Exclude sites for which the expression is true\n\
   -i, --include EXPR             Select sites for which the expression is true\n\
   -r, --regions REGION           Restrict to comma-separated list of regions\n\
   -R, --regions-file FILE        Restrict to regions listed in a file\n\
       --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n\
   -t, --targets REGION           Similar to -r but streams rather than index-jumps\n\
   -T, --targets-file FILE        Similar to -R but streams rather than index-jumps\n\
       --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n\
VCF output options:\n\
       --no-version               Do not append version and command line to the header\n\
   -o, --output FILE              Write output to a file [standard output]\n\
   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n\
       --threads INT              Use multithreading with <int> worker threads [0]\n\
Plugin options:\n\
   -h, --help                     List plugin's options\n\
   -l, --list-plugins             List available plugins. See BCFTOOLS_PLUGINS environment variable and man page for details\n\
   -v, --verbosity INT            Verbosity level\n\
   -V, --version                  Print version string and exit\n\
   -W, --write-index[=FMT]        Automatically index the output files [off]\n\
\n";

#[derive(Clone, Copy, PartialEq, Eq)]
enum OutKind {
    VcfText,
    VcfGz,
    Bcf,
}

fn parse_out_kind(raw: &str) -> OutKind {
    match raw.as_bytes().first().copied() {
        Some(b'z') => OutKind::VcfGz,
        Some(b'u' | b'b') => OutKind::Bcf,
        _ => OutKind::VcfText,
    }
}

#[derive(Clone, Copy)]
struct Plugin {
    name: &'static str,
    about: &'static str,
}

const PLUGINS: &[Plugin] = &[
    Plugin {
        name: "GTisec",
        about: "Count genotype intersections between sample groups.",
    },
    Plugin {
        name: "GTsubset",
        about: "Output positions where selected samples have exclusive genotypes.",
    },
    Plugin {
        name: "ad-bias",
        about: "Detect allele-depth strand and position bias.",
    },
    Plugin {
        name: "add-variantkey",
        about: "Add VariantKey INFO annotations.",
    },
    Plugin {
        name: "af-dist",
        about: "Calculate allele-frequency distribution diagnostics.",
    },
    Plugin {
        name: "allele-length",
        about: "Calculate allele length statistics.",
    },
    Plugin {
        name: "check-ploidy",
        about: "Check sex/ploidy from genotype data.",
    },
    Plugin {
        name: "check-sparsity",
        about: "Check sparse VCF/BCF genotype representation.",
    },
    Plugin {
        name: "color-chrs",
        about: "Color chromosome names in VCF output.",
    },
    Plugin {
        name: "contrast",
        about: "Compare allele counts between groups of samples.",
    },
    Plugin {
        name: "counts",
        about: "Count samples, SNPs, indels, and total sites.",
    },
    Plugin {
        name: "dosage",
        about: "Print genotype dosage values.",
    },
    Plugin {
        name: "fill-AN-AC",
        about: "Fill INFO/AN and INFO/AC fields. Deprecated in favor of fill-tags.",
    },
    Plugin {
        name: "fill-from-fasta",
        about: "Fill REF/ALT or INFO values from a FASTA reference.",
    },
    Plugin {
        name: "fill-tags",
        about: "Fill INFO tags from FORMAT/genotype fields.",
    },
    Plugin {
        name: "fixploidy",
        about: "Set or fix genotype ploidy.",
    },
    Plugin {
        name: "fixref",
        about: "Check and fix reference allele orientation.",
    },
    Plugin {
        name: "frameshifts",
        about: "Annotate frameshift indels.",
    },
    Plugin {
        name: "guess-ploidy",
        about: "Guess sample ploidy from genotype data.",
    },
    Plugin {
        name: "gvcfz",
        about: "Compress gVCF blocks.",
    },
    Plugin {
        name: "impute-info",
        about: "Add imputation INFO metrics.",
    },
    Plugin {
        name: "indel-stats",
        about: "Calculate indel statistics.",
    },
    Plugin {
        name: "isecGT",
        about: "Collect genotype intersection counts.",
    },
    Plugin {
        name: "mendelian2",
        about: "Find Mendelian inconsistency sites.",
    },
    Plugin {
        name: "missing2ref",
        about: "Set missing genotypes to reference alleles.",
    },
    Plugin {
        name: "parental-origin",
        about: "Infer parental origin of alleles.",
    },
    Plugin {
        name: "prune",
        about: "Prune variants by linkage or distance.",
    },
    Plugin {
        name: "remove-overlaps",
        about: "Remove overlapping variants.",
    },
    Plugin {
        name: "scatter",
        about: "Scatter variants into genomic chunks.",
    },
    Plugin {
        name: "setGT",
        about: "Set genotypes matching a query.",
    },
    Plugin {
        name: "smpl-stats",
        about: "Calculate per-sample statistics.",
    },
    Plugin {
        name: "split",
        about: "Split VCF by sample groups or annotations.",
    },
    Plugin {
        name: "split-vep",
        about: "Extract fields from VEP CSQ annotations.",
    },
    Plugin {
        name: "tag2tag",
        about: "Convert between related FORMAT/INFO tags.",
    },
    Plugin {
        name: "trio-dnm2",
        about: "Find de novo mutations in trios.",
    },
    Plugin {
        name: "trio-stats",
        about: "Calculate trio statistics.",
    },
    Plugin {
        name: "trio-switch-rate",
        about: "Estimate trio switch error rates.",
    },
    Plugin {
        name: "variant-distance",
        about: "Calculate distances between variants.",
    },
    Plugin {
        name: "variantkey-hex",
        about: "Generate VariantKey lookup tables in hexadecimal form.",
    },
    Plugin {
        name: "vcf2table",
        about: "Convert VCF records to tabular output.",
    },
    Plugin {
        name: "vrfs",
        about: "Calculate variant read frequency statistics.",
    },
];

pub fn count_plugins() -> usize {
    PLUGINS.len()
}

pub fn main(argv: &[OsString]) -> ExitCode {
    match run(argv) {
        Ok(code) => code,
        Err(e) => {
            eprintln!("{}", fmt_etag("main_plugin", &e.to_string()));
            ExitCode::FAILURE
        }
    }
}

fn run(argv: &[OsString]) -> io::Result<ExitCode> {
    let mut list = false;
    let mut verbose = 0usize;
    let mut help = false;
    let mut version = false;
    let mut plugin_name: Option<String> = None;
    let mut input: Option<String> = None;
    let mut output: Option<String> = None;
    let mut output_kind = OutKind::VcfText;
    // Plugin-specific options consumed for the plugins ported so far.
    let mut direction: Option<String> = None;
    let mut tag_name: Option<String> = None;
    let mut use_missing = false;
    let mut past_separator = false;
    let mut replace = false;
    let mut threshold: f64 = 0.1;
    let mut conversion: Option<&'static str> = None;

    let mut iter = argv.iter().skip(1).peekable();
    while let Some(arg) = iter.next() {
        let raw = arg.to_string_lossy();
        match raw.as_ref() {
            "-l" | "--list-plugins" => list = true,
            "-lv" => {
                list = true;
                verbose += 1;
            }
            "-lvv" => {
                list = true;
                verbose += 2;
            }
            "-lvvv" => {
                list = true;
                verbose += 3;
            }
            "-h" | "--help" | "-?" => help = true,
            "-V" | "--version" => version = true,
            "-v" | "--verbose" => verbose += 1,
            "-vv" => verbose += 2,
            "-vvv" => verbose += 3,
            "-W" | "--write-index" => {}
            "-o" | "--output" => {
                output = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-O" | "--output-type" => {
                if let Some(v) = iter.next() {
                    output_kind = parse_out_kind(&v.to_string_lossy());
                }
            }
            _ if raw.starts_with("--output=") => {
                output = Some(raw["--output=".len()..].to_owned());
            }
            _ if raw.starts_with("--output-type=") => {
                output_kind = parse_out_kind(&raw["--output-type=".len()..]);
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => {
                output = Some(raw[2..].to_owned());
            }
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                output_kind = parse_out_kind(&raw[2..]);
            }
            "-i" | "--include" | "-e" | "--exclude" | "--regions" | "-R" | "--regions-file"
            | "--targets" | "-T" | "--targets-file" | "--regions-overlap" | "--targets-overlap"
            | "--threads" => {
                let _ = iter.next();
            }
            // `-r` is `--regions` (value) before `--`, `--replace` (flag) after.
            "-r" => {
                if past_separator {
                    replace = true;
                } else {
                    let _ = iter.next();
                }
            }
            "--replace" => replace = true,
            // `-t` is `--targets` (value) before `--`, `--threshold` after.
            "-t" | "--threshold" => {
                if past_separator || raw == "--threshold" {
                    if let Some(v) = iter.next() {
                        threshold = v.to_string_lossy().parse().unwrap_or(threshold);
                    }
                } else {
                    let _ = iter.next();
                }
            }
            "--gl-to-pl" => conversion = Some("gl-to-pl"),
            "--gp-to-gt" => conversion = Some("gp-to-gt"),
            "--gl-to-gp" => conversion = Some("gl-to-gp"),
            "-d" | "--direction" => {
                direction = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-n" | "--tag-name" => {
                tag_name = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-m" | "--use-missing" => use_missing = true,
            _ if raw.starts_with("--direction=") => {
                direction = Some(raw["--direction=".len()..].to_owned());
            }
            _ if raw.starts_with("--tag-name=") => {
                tag_name = Some(raw["--tag-name=".len()..].to_owned());
            }
            _ if raw.starts_with("-d") && raw.len() > 2 => {
                direction = Some(raw[2..].to_owned());
            }
            _ if raw.starts_with("-n") && raw.len() > 2 => {
                tag_name = Some(raw[2..].to_owned());
            }
            "--" => past_separator = true,
            "--no-version" => {}
            _ if raw.starts_with("--verbosity=") => {
                verbose = raw["--verbosity=".len()..].parse().unwrap_or(verbose);
            }
            _ if raw.starts_with("--verbose=") => {
                verbose = raw["--verbose=".len()..].parse().unwrap_or(verbose);
            }
            _ if raw.starts_with("-v") && raw.len() > 2 => {
                verbose = raw[2..].parse().unwrap_or(raw.len() - 1);
            }
            _ if raw.starts_with("--") => {}
            _ if raw.starts_with('-') => {}
            _ if plugin_name.is_none() => plugin_name = Some(raw.into_owned()),
            _ if input.is_none() => input = Some(raw.into_owned()),
            _ => {}
        }
    }

    if list {
        let mut out = io::stdout().lock();
        list_plugins(&mut out, verbose > 0)?;
        return Ok(ExitCode::SUCCESS);
    }

    let Some(name) = plugin_name else {
        let mut err = io::stderr().lock();
        err.write_all(USAGE.as_bytes())?;
        return Ok(ExitCode::FAILURE);
    };

    let Some(plugin) = find_plugin(&name) else {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("plugin '{name}' is not registered"),
        ));
    };

    if version {
        println!(
            "bcftools  {} using htslib {}",
            crate::version::BCFTOOLS_VERSION,
            crate::version::HTSLIB_RS_VERSION
        );
        println!(
            "plugin at {} using htslib {}\n",
            crate::version::BCFTOOLS_VERSION,
            crate::version::HTSLIB_RS_VERSION
        );
        return Ok(ExitCode::SUCCESS);
    }

    if help {
        let mut err = io::stderr().lock();
        writeln!(err, "\nAbout:   {}", plugin.about)?;
        writeln!(
            err,
            "Usage:   bcftools +{} [General Options] -- [Plugin Options]",
            plugin.name
        )?;
        writeln!(
            err,
            "\nThe '{}' plugin is registered but its record-processing implementation is not yet ported.",
            plugin.name
        )?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "counts" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let report = crate::commands::plugins::counts::run(Path::new(&input))?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "allele-length" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let report = crate::commands::plugins::allele_length::run(Path::new(&input))?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "check-ploidy" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let report = crate::commands::plugins::check_ploidy::run(Path::new(&input), use_missing)?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "tag2tag" {
        use crate::commands::plugins::tag2tag::{self, Conversion};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let conv = match conversion {
            Some("gl-to-pl") => Conversion::GlToPl,
            Some("gp-to-gt") => Conversion::GpToGt,
            Some("gl-to-gp") => {
                return Err(io::Error::new(
                    io::ErrorKind::Unsupported,
                    "tag2tag --gl-to-gp is not supported in this local slice",
                ));
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "tag2tag requires one of --gl-to-pl or --gp-to-gt in this local slice",
                ));
            }
        };
        let vcf = tag2tag::run(Path::new(&input), conv, replace, threshold)?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "missing2ref" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let vcf = crate::commands::plugins::missing2ref::run(Path::new(&input))?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "fill-AN-AC" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let vcf = crate::commands::plugins::fill_an_ac::run(Path::new(&input))?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "variant-distance" {
        use crate::commands::plugins::variant_distance::{self, Direction};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let dir = match direction.as_deref() {
            None => Direction::Nearest,
            Some(d) => Direction::parse(d).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown -d direction '{d}' (expected nearest|fwd|rev|both)"),
                )
            })?,
        };
        let tag = tag_name.as_deref().unwrap_or("DIST");
        let vcf = variant_distance::run(Path::new(&input), dir, tag)?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    Err(io::Error::new(
        io::ErrorKind::Unsupported,
        format!(
            "plugin '{}' is registered but not yet implemented",
            plugin.name
        ),
    ))
}

fn write_plugin_output(bytes: &[u8], output: Option<&str>, kind: OutKind) -> io::Result<()> {
    use std::fs::File;
    match output {
        Some(path) if path != "-" => write_kind(bytes, kind, File::create(path)?),
        _ => write_kind(bytes, kind, io::stdout().lock()),
    }
}

fn write_kind<W: Write>(bytes: &[u8], kind: OutKind, out: W) -> io::Result<()> {
    match kind {
        OutKind::VcfText => {
            let mut w = io::BufWriter::new(out);
            w.write_all(bytes)
        }
        OutKind::VcfGz => {
            let mut bgzf = htslib_rs::bgzf::io::Writer::new(out);
            bgzf.write_all(bytes)?;
            bgzf.finish().map(|_| ())
        }
        OutKind::Bcf => {
            use htslib_rs::vcf::variant::io::Write as _;
            let mut reader = htslib_rs::vcf::io::Reader::new(std::io::BufReader::new(bytes));
            let header = reader.read_header()?;
            let mut writer = htslib_rs::bcf::io::Writer::new(out);
            writer.write_variant_header(&header)?;
            for result in reader.records() {
                writer.write_variant_record(&header, &result?)?;
            }
            writer.try_finish()
        }
    }
}

fn list_plugins<W: Write>(out: &mut W, verbose: bool) -> io::Result<()> {
    for plugin in PLUGINS {
        if verbose {
            writeln!(out, "\n-- {} --\n{}", plugin.name, plugin.about)?;
        } else {
            writeln!(out, "{}", plugin.name)?;
        }
    }
    Ok(())
}

fn find_plugin(name: &str) -> Option<Plugin> {
    PLUGINS.iter().copied().find(|plugin| plugin.name == name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn registry_tracks_upstream_plugin_count() {
        assert_eq!(count_plugins(), 41);
        assert!(find_plugin("fill-tags").is_some());
        assert!(find_plugin("trio-dnm2").is_some());
    }
}
