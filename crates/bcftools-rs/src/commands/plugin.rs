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
    let mut extra: Option<String> = None;
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
    // af-dist options.
    let mut af_tag: Option<String> = None;
    let mut dev_bins: Option<String> = None;
    let mut prob_bins: Option<String> = None;
    // remove-overlaps options.
    let mut mark_expr: Option<String> = None;
    let mut mark_tag: Option<String> = None;
    let mut missing_expr: Option<String> = None;
    let mut reverse = false;
    let mut out_type_raw: Option<String> = None;
    // ad-bias options.
    let mut samples_file: Option<String> = None;
    let mut clean_vcf = false;
    let mut min_dp: Option<i32> = None;
    let mut min_alt_dp: Option<i32> = None;
    let mut ad_threshold: Option<f64> = None;
    // prune options.
    let mut window: Option<String> = None;
    let mut nsites: Option<i32> = None;
    let mut nsites_mode: Option<String> = None;
    let mut prune_af_tag: Option<String> = None;
    let mut prune_annot: Option<String> = None;
    let mut prune_max: Option<String> = None;
    let mut prune_set_filter: Option<String> = None;
    // contrast options.
    let mut contrast_annots: Option<String> = None;
    let mut contrast_control: Option<String> = None;
    let mut contrast_case: Option<String> = None;
    let mut force_samples = false;
    // fixref options.
    let mut fixref_fasta: Option<String> = None;
    let mut fixref_mode: Option<String> = None;
    let mut fixref_discard = false;
    // PED-driven plugins (trio-switch-rate, trio-stats, ...).
    let mut ped_file: Option<String> = None;
    let mut trio_stats_alt: Option<i32> = None;
    let mut trio_stats_debug: Option<String> = None;
    // dosage options.
    let mut tags_list: Option<String> = None;
    // guess-ploidy options.
    let mut gp_region: Option<String> = None;

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
                    let v = v.to_string_lossy();
                    out_type_raw = Some(v.clone().into_owned());
                    output_kind = parse_out_kind(&v);
                }
            }
            _ if raw.starts_with("--output=") => {
                output = Some(raw["--output=".len()..].to_owned());
            }
            _ if raw.starts_with("--output-type=") => {
                out_type_raw = Some(raw["--output-type=".len()..].to_owned());
                output_kind = parse_out_kind(&raw["--output-type=".len()..]);
            }
            _ if raw.starts_with("-o") && raw.len() > 2 => {
                output = Some(raw[2..].to_owned());
            }
            _ if raw.starts_with("-O") && raw.len() > 2 => {
                out_type_raw = Some(raw[2..].to_owned());
                output_kind = parse_out_kind(&raw[2..]);
            }
            // guess-ploidy region restriction (`-r X` / `-rX` / `-R file`).
            "-r" | "--regions" | "-R" | "--regions-file"
                if plugin_name.as_deref() == Some("guess-ploidy") =>
            {
                gp_region = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("-r")
                && raw.len() > 2
                && plugin_name.as_deref() == Some("guess-ploidy") =>
            {
                gp_region = Some(raw[2..].to_owned());
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
            // `-t` is `--af-tag` for af-dist, `--threshold` after `--` for
            // tag2tag, otherwise `--targets` (value, ignored).
            "--tags" => {
                tags_list = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--tags=") => {
                tags_list = Some(raw["--tags=".len()..].to_owned());
            }
            "-t" | "--threshold" | "--af-tag" => {
                if plugin_name.as_deref() == Some("dosage") {
                    tags_list = iter.next().map(|s| s.to_string_lossy().into_owned());
                } else if raw == "--af-tag" || plugin_name.as_deref() == Some("af-dist") {
                    af_tag = iter.next().map(|s| s.to_string_lossy().into_owned());
                } else if past_separator || raw == "--threshold" {
                    if let Some(v) = iter.next() {
                        if plugin_name.as_deref() == Some("ad-bias") {
                            ad_threshold = v.to_string_lossy().parse().ok();
                        } else {
                            threshold = v.to_string_lossy().parse().unwrap_or(threshold);
                        }
                    }
                } else {
                    let _ = iter.next();
                }
            }
            _ if raw.starts_with("--af-tag=") => {
                af_tag = Some(raw["--af-tag=".len()..].to_owned());
            }
            "--dev-bins" => {
                dev_bins = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--dev-bins=") => {
                dev_bins = Some(raw["--dev-bins=".len()..].to_owned());
            }
            "--prob-bins" => {
                prob_bins = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--prob-bins=") => {
                prob_bins = Some(raw["--prob-bins=".len()..].to_owned());
            }
            "-d" if plugin_name.as_deref() == Some("af-dist") => {
                dev_bins = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-d" if plugin_name.as_deref() == Some("ad-bias") => {
                min_dp = iter.next().and_then(|s| s.to_string_lossy().parse().ok());
            }
            "-d" if plugin_name.as_deref() == Some("fixref") => {
                fixref_discard = true;
            }
            "-d" | "--debug" if plugin_name.as_deref() == Some("trio-stats") => {
                trio_stats_debug = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-a" | "--alt-trios" if plugin_name.as_deref() == Some("trio-stats") => {
                trio_stats_alt = iter.next().and_then(|s| s.to_string_lossy().parse().ok());
            }
            "-p" if plugin_name.as_deref() == Some("af-dist") => {
                prob_bins = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-p" | "--ped" => {
                ped_file = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "--gl-to-pl" => conversion = Some("gl-to-pl"),
            "--gp-to-gt" => conversion = Some("gp-to-gt"),
            "--gl-to-gp" => conversion = Some("gl-to-gp"),
            "-d" | "--direction" => {
                direction = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            // `-n` is `--nsites-per-win` for prune, `--tag-name` otherwise.
            "-n" | "--tag-name" | "--nsites-per-win" => {
                if raw == "--nsites-per-win" || plugin_name.as_deref() == Some("prune") {
                    nsites = iter.next().and_then(|s| s.to_string_lossy().parse().ok());
                } else {
                    tag_name = iter.next().map(|s| s.to_string_lossy().into_owned());
                }
            }
            "-N" | "--nsites-per-win-mode" => {
                nsites_mode = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-w" | "--window" => {
                window = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--window=") => {
                window = Some(raw["--window=".len()..].to_owned());
            }
            "--AF-tag" => {
                prune_af_tag = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--AF-tag=") => {
                prune_af_tag = Some(raw["--AF-tag=".len()..].to_owned());
            }
            // ad-bias options.
            "-s" | "--samples" => {
                samples_file = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--samples=") => {
                samples_file = Some(raw["--samples=".len()..].to_owned());
            }
            "-c" | "--clean-vcf" => clean_vcf = true,
            "--min-dp" => {
                min_dp = iter.next().and_then(|s| s.to_string_lossy().parse().ok());
            }
            "--min-alt-dp" => {
                min_alt_dp = iter.next().and_then(|s| s.to_string_lossy().parse().ok());
            }
            "-a" | "--annots" if plugin_name.as_deref() == Some("contrast") => {
                contrast_annots = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-0" | "--control-samples" | "--bg-samples" => {
                contrast_control = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-1" | "--case-samples" | "--novel-samples" => {
                contrast_case = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "--force-samples" => force_samples = true,
            "-a" if plugin_name.as_deref() == Some("ad-bias") => {
                min_alt_dp = iter.next().and_then(|s| s.to_string_lossy().parse().ok());
            }
            "-a" | "--annotate" if plugin_name.as_deref() == Some("prune") => {
                prune_annot = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-f" | "--set-filter" if plugin_name.as_deref() == Some("prune") => {
                prune_set_filter = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "-f" | "--fasta-ref" if plugin_name.as_deref() == Some("fixref") => {
                fixref_fasta = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "--mode" if plugin_name.as_deref() == Some("fixref") => {
                fixref_mode = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            "--discard" if plugin_name.as_deref() == Some("fixref") => {
                fixref_discard = true;
            }
            "--max" if plugin_name.as_deref() == Some("prune") => {
                prune_max = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            // `-m` is `--max` for prune, `--mark EXPR` for remove-overlaps,
            // the boolean `--use-missing` flag otherwise.
            "-m" => {
                if plugin_name.as_deref() == Some("prune") {
                    prune_max = iter.next().map(|s| s.to_string_lossy().into_owned());
                } else if plugin_name.as_deref() == Some("remove-overlaps") {
                    mark_expr = iter.next().map(|s| s.to_string_lossy().into_owned());
                } else if plugin_name.as_deref() == Some("fixref") {
                    fixref_mode = iter.next().map(|s| s.to_string_lossy().into_owned());
                } else {
                    use_missing = true;
                }
            }
            "--mark" => {
                mark_expr = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--mark=") => {
                mark_expr = Some(raw["--mark=".len()..].to_owned());
            }
            "--use-missing" => use_missing = true,
            "-M" | "--mark-tag" => {
                mark_tag = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--mark-tag=") => {
                mark_tag = Some(raw["--mark-tag=".len()..].to_owned());
            }
            "--reverse" => reverse = true,
            "--missing" => {
                missing_expr = iter.next().map(|s| s.to_string_lossy().into_owned());
            }
            _ if raw.starts_with("--missing=") => {
                missing_expr = Some(raw["--missing=".len()..].to_owned());
            }
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
            // Trailing positional, e.g. `+variantkey-hex in.vcf <dir>`.
            _ if extra.is_none() => extra = Some(raw.into_owned()),
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

    if plugin.name == "add-variantkey" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let vcf = crate::commands::plugins::add_variantkey::run(Path::new(&input))?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "variantkey-hex" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        // Upstream's optional output directory positional; defaults to "./".
        let dir = extra.unwrap_or_else(|| "./".to_owned());
        let summary = crate::commands::plugins::variantkey_hex::run(Path::new(&input), &dir)?;
        io::stdout().lock().write_all(summary.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "remove-overlaps" {
        use crate::commands::plugins::remove_overlaps::{self, Mark};
        let input = input.unwrap_or_else(|| "-".to_owned());
        if let Some(m) = missing_expr.as_deref() {
            // --missing only applies to the `min(QUAL)` expression mode,
            // which needs the not-yet-ported filter engine.
            let _ = m;
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "remove-overlaps --missing requires the min(QUAL) expression mode (filter engine not yet ported)",
            ));
        }
        let mode = Mark::parse(mark_expr.as_deref().unwrap_or("overlap"))
            .map_err(|e| io::Error::new(io::ErrorKind::Unsupported, e))?;
        let text_list = out_type_raw.as_deref().is_some_and(|o| o.starts_with('t'));
        let vcf = remove_overlaps::run(
            Path::new(&input),
            mode,
            mark_tag.as_deref(),
            reverse,
            text_list,
        )?;
        if text_list {
            write_plugin_output(vcf.as_bytes(), output.as_deref(), OutKind::VcfText)?;
        } else {
            write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        }
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "trio-stats" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let Some(ped) = ped_file.as_deref() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "trio-stats requires -p/--ped",
            ));
        };
        let dbg = trio_stats_debug.unwrap_or_default();
        let dbg_mendel = dbg
            .split(',')
            .any(|t| t.eq_ignore_ascii_case("mendel-errors"));
        let dbg_tr = dbg
            .split(',')
            .any(|t| t.eq_ignore_ascii_case("transmitted"));
        let report = crate::commands::plugins::trio_stats::run(
            Path::new(&input),
            Path::new(ped),
            trio_stats_alt.unwrap_or(0),
            dbg_mendel,
            dbg_tr,
        )?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "trio-switch-rate" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let Some(ped) = ped_file.as_deref() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "trio-switch-rate requires -p/--ped",
            ));
        };
        let report =
            crate::commands::plugins::trio_switch_rate::run(Path::new(&input), Path::new(ped))?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "fixref" {
        use crate::commands::plugins::fixref::{self, Mode};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let Some(fa) = fixref_fasta.as_deref() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "fixref requires -f/--fasta-ref",
            ));
        };
        let mode = Mode::parse(fixref_mode.as_deref().unwrap_or("stats"))
            .map_err(|e| io::Error::new(io::ErrorKind::Unsupported, e))?;
        let vcf = fixref::run(
            Path::new(&input),
            Path::new(fa),
            mode,
            fixref_discard,
            "FIXREF",
        )?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "contrast" {
        use crate::commands::plugins::contrast::{self, Annots};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let annots = Annots::parse(contrast_annots.as_deref().unwrap_or("PASSOC,FASSOC"))
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
        let (Some(ctrl), Some(case)) = (contrast_control.as_deref(), contrast_case.as_deref())
        else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "contrast requires -0/--control-samples and -1/--case-samples",
            ));
        };
        let vcf = contrast::run(Path::new(&input), annots, ctrl, case, force_samples)?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "prune" {
        use crate::commands::plugins::prune::{self, LdAnnot, Mode};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let win = match window.as_deref() {
            Some(w) => prune::parse_window(w)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?,
            None => 0,
        };

        // LD path: `-a`/`--annotate` or `-m`/`--max`.
        if prune_annot.is_some() || prune_max.is_some() {
            let mut annot = LdAnnot::default();
            if let Some(a) = prune_annot.as_deref() {
                for t in a.split(',') {
                    match t.to_ascii_uppercase().as_str() {
                        "R2" => annot.annot[0] = true,
                        "LD" => annot.annot[1] = true,
                        "RD" | "HD" => annot.annot[2] = true,
                        "COUNT" => {
                            return Err(io::Error::new(
                                io::ErrorKind::Unsupported,
                                "prune -a count is not supported in this slice",
                            ));
                        }
                        other => {
                            return Err(io::Error::new(
                                io::ErrorKind::InvalidInput,
                                format!("The tag \"{other}\" is not supported"),
                            ));
                        }
                    }
                }
            }
            let mut max: [Option<f64>; 3] = [None; 3];
            if let Some(m) = prune_max.as_deref() {
                let (idx, num) = if let Some(v) = m.strip_prefix("R2=") {
                    (0, v)
                } else if let Some(v) = m.strip_prefix("LD=") {
                    (1, v)
                } else if let Some(v) = m.strip_prefix("RD=").or_else(|| m.strip_prefix("HD=")) {
                    (2, v)
                } else if m.starts_with("count=") {
                    return Err(io::Error::new(
                        io::ErrorKind::Unsupported,
                        "prune -m count= is not supported in this slice",
                    ));
                } else {
                    (0, m)
                };
                max[idx] = Some(num.parse().map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("Could not parse: --max {m}"),
                    )
                })?);
            }
            let vcf = prune::run_ld(
                Path::new(&input),
                win,
                annot,
                max,
                prune_set_filter.as_deref(),
            )?;
            write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
            return Ok(ExitCode::SUCCESS);
        }

        let Some(n) = nsites else {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "prune in this slice requires -n/--nsites-per-win or -a/-m",
            ));
        };
        let mode = Mode::parse(nsites_mode.as_deref().unwrap_or("maxAF"))
            .map_err(|e| io::Error::new(io::ErrorKind::Unsupported, e))?;
        let vcf = prune::run(Path::new(&input), win, n, mode, prune_af_tag.as_deref())?;
        write_plugin_output(vcf.as_bytes(), output.as_deref(), output_kind)?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "ad-bias" {
        if clean_vcf {
            return Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "ad-bias -c/--clean-vcf (VCF allele-removal output) is not supported in this local slice",
            ));
        }
        let input = input.unwrap_or_else(|| "-".to_owned());
        let Some(samples) = samples_file.as_deref() else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "ad-bias requires -s/--samples <file>",
            ));
        };
        let report = crate::commands::plugins::ad_bias::run(
            Path::new(&input),
            Path::new(samples),
            ad_threshold.unwrap_or(1e-3),
            min_dp.unwrap_or(0),
            min_alt_dp.unwrap_or(1),
        )?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "guess-ploidy" {
        use crate::commands::plugins::guess_ploidy::{self, Tag};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let report = guess_ploidy::run(
            Path::new(&input),
            Tag::Pl, // default; auto-switches PL->GL->GT on header
            gp_region.as_deref(),
            1e-3,
            0.5,
            false,
            verbose as u32,
        )?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "dosage" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let tags: Vec<String> = tags_list
            .as_deref()
            .unwrap_or("PL,GL,GT")
            .split(',')
            .map(|s| s.to_owned())
            .collect();
        let report = crate::commands::plugins::dosage::run(Path::new(&input), &tags)?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "indel-stats" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let report = crate::commands::plugins::indel_stats::run(Path::new(&input))?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "smpl-stats" {
        let input = input.unwrap_or_else(|| "-".to_owned());
        let report = crate::commands::plugins::smpl_stats::run(Path::new(&input))?;
        io::stdout().lock().write_all(report.as_bytes())?;
        return Ok(ExitCode::SUCCESS);
    }

    if plugin.name == "af-dist" {
        use crate::commands::plugins::af_dist::{self, DEFAULT_BINS};
        let input = input.unwrap_or_else(|| "-".to_owned());
        let af = af_tag.as_deref().unwrap_or("AF");
        let dev = dev_bins.as_deref().unwrap_or(DEFAULT_BINS);
        let prob = prob_bins.as_deref().unwrap_or(DEFAULT_BINS);
        let report = af_dist::run(Path::new(&input), af, dev, prob)?;
        io::stdout().lock().write_all(report.as_bytes())?;
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
