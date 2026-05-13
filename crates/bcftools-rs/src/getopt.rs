//! Minimal `getopt_long`-style argument parser for `bcftools` subcommands.
//!
//! Upstream uses POSIX `getopt_long` with shapes like `"h:n:s:v:"` plus a long
//! options table that maps to the same chars or to extra integer codes. We
//! reproduce that here, including:
//!
//! - attached short values (`-Oz` ≡ `-O z` ≡ `-O=z`)
//! - bundled boolean shorts (`-vh` ≡ `-v -h`)
//! - `--name=value` and `--name value` for long options
//! - `--` ends option processing; remaining argv is positional
//! - returns the matched option's char (or its long-only int code) plus its
//!   string value if the option takes one
//!
//! It does NOT replicate libc's global state (`optarg`, `optind`, `optopt`);
//! the parser is owned by the caller.

use std::ffi::OsString;

/// Whether a long option takes an argument.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HasArg {
    /// No argument (`no_argument` in libc).
    None,
    /// Required argument (`required_argument` in libc).
    Required,
    /// Optional argument — accepted only as `--name=value`, never `--name value`.
    Optional,
}

/// One long option entry, parallel to libc's `struct option`.
#[derive(Debug, Clone)]
pub struct LongOpt {
    /// Long option name (without the leading `--`).
    pub name: &'static str,
    /// Whether the option takes an argument.
    pub has_arg: HasArg,
    /// Integer code returned when this option matches. By convention, ASCII
    /// codes in `0..128` map to a corresponding short option; values outside
    /// that range are long-only sentinels.
    pub code: i32,
}

impl LongOpt {
    /// Convenience constructor.
    pub const fn new(name: &'static str, has_arg: HasArg, code: i32) -> Self {
        Self {
            name,
            has_arg,
            code,
        }
    }
}

/// One entry of the short-option spec parsed from a getopt-style string like
/// `"h:n:s:v:"`. The trailing colon means "required argument".
#[derive(Debug, Clone, Copy)]
struct ShortOpt {
    ch: char,
    has_arg: HasArg,
}

/// One option-match emitted by the parser.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Match {
    /// `LongOpt::code` for long matches, or the ASCII code of the short letter.
    pub code: i32,
    /// String form of the option's argument, if it took one.
    pub value: Option<String>,
}

/// Parser state.
#[derive(Debug)]
pub struct Getopt<'a> {
    short: Vec<ShortOpt>,
    long: &'a [LongOpt],
    /// Argv with `args[0]` being the subcommand name (e.g. `"head"`).
    argv: &'a [OsString],
    /// Index of the next argv slot to look at.
    cursor: usize,
    /// When parsing bundled shorts like `-vh`, the byte offset within
    /// `argv[cursor]` after the leading `-`.
    bundled_offset: Option<usize>,
}

impl<'a> Getopt<'a> {
    /// Build a new parser from a getopt-style short spec, a long-option table,
    /// and the argv slice received by the subcommand. `argv[0]` is treated as
    /// the program name (skipped).
    pub fn new(short_spec: &str, long: &'a [LongOpt], argv: &'a [OsString]) -> Self {
        let short = parse_short_spec(short_spec);
        Self {
            short,
            long,
            argv,
            cursor: 1,
            bundled_offset: None,
        }
    }

    /// Index of the first non-option argument once iteration completes.
    pub fn optind(&self) -> usize {
        self.cursor
    }

    /// Yield the next option match. Returns `Ok(None)` when option processing
    /// is over (the next argv slot is positional or there are none left).
    /// Returns `Err(msg)` for unknown options or missing required arguments;
    /// the caller is expected to print usage and exit non-zero.
    ///
    /// Not implemented as `Iterator::next` because the state-machine semantics
    /// (bundled short clusters, error propagation) don't compose naturally
    /// with the standard iterator protocol.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Result<Option<Match>, String> {
        // Continue consuming a bundled short cluster like `-abc`.
        if let Some(off) = self.bundled_offset {
            let arg = self.argv[self.cursor].to_string_lossy().into_owned();
            let bytes = arg.as_bytes();
            if off >= bytes.len() {
                self.bundled_offset = None;
                self.cursor += 1;
            } else {
                let ch = bytes[off] as char;
                let spec = self
                    .lookup_short(ch)
                    .ok_or_else(|| format!("unrecognized option `-{ch}`"))?;
                if spec.has_arg == HasArg::None {
                    self.bundled_offset = Some(off + 1);
                    return Ok(Some(Match {
                        code: ch as i32,
                        value: None,
                    }));
                }
                // Required-arg short: the rest of the cluster (or next argv) is
                // the value, exactly like libc getopt.
                let rest = &arg[off + 1..];
                let value = if !rest.is_empty() {
                    rest.trim_start_matches('=').to_string()
                } else {
                    self.cursor += 1;
                    if self.cursor >= self.argv.len() {
                        return Err(format!("option `-{ch}` requires an argument"));
                    }
                    self.argv[self.cursor].to_string_lossy().into_owned()
                };
                self.bundled_offset = None;
                self.cursor += 1;
                return Ok(Some(Match {
                    code: ch as i32,
                    value: Some(value),
                }));
            }
        }

        if self.cursor >= self.argv.len() {
            return Ok(None);
        }
        let arg = self.argv[self.cursor].to_string_lossy().into_owned();

        // `--` ends option processing.
        if arg == "--" {
            self.cursor += 1;
            return Ok(None);
        }
        // Bare `-` is a positional (means stdin).
        if arg == "-" {
            return Ok(None);
        }
        // Long option.
        if let Some(rest) = arg.strip_prefix("--") {
            let (name, eq_value) = match rest.find('=') {
                Some(i) => (&rest[..i], Some(rest[i + 1..].to_string())),
                None => (rest, None),
            };
            let lo = self
                .lookup_long(name)
                .ok_or_else(|| format!("unrecognized option `--{name}`"))?;
            self.cursor += 1;
            let value = match (lo.has_arg, eq_value) {
                (HasArg::None, Some(_)) => {
                    return Err(format!("option `--{name}` does not take an argument"));
                }
                (HasArg::None, None) => None,
                (HasArg::Required, Some(v)) => Some(v),
                (HasArg::Required, None) => {
                    if self.cursor >= self.argv.len() {
                        return Err(format!("option `--{name}` requires an argument"));
                    }
                    let v = self.argv[self.cursor].to_string_lossy().into_owned();
                    self.cursor += 1;
                    Some(v)
                }
                (HasArg::Optional, v) => v,
            };
            return Ok(Some(Match {
                code: lo.code,
                value,
            }));
        }
        // Short option(s).
        if let Some(rest) = arg.strip_prefix('-') {
            if rest.is_empty() {
                return Ok(None);
            }
            let bytes = rest.as_bytes();
            let ch = bytes[0] as char;
            let spec = self
                .lookup_short(ch)
                .ok_or_else(|| format!("unrecognized option `-{ch}`"))?;
            if spec.has_arg == HasArg::None {
                if bytes.len() == 1 {
                    self.cursor += 1;
                } else {
                    // Bundle: `-vh` → continue consuming after the first char.
                    self.bundled_offset = Some(2);
                }
                return Ok(Some(Match {
                    code: ch as i32,
                    value: None,
                }));
            }
            // Required-arg short.
            let value = if bytes.len() > 1 {
                let rest_str = &rest[1..];
                rest_str.trim_start_matches('=').to_string()
            } else {
                self.cursor += 1;
                if self.cursor >= self.argv.len() {
                    return Err(format!("option `-{ch}` requires an argument"));
                }
                self.argv[self.cursor].to_string_lossy().into_owned()
            };
            self.cursor += 1;
            return Ok(Some(Match {
                code: ch as i32,
                value: Some(value),
            }));
        }
        // Positional.
        Ok(None)
    }

    /// All argv entries from `optind` onwards. Useful for collecting positional
    /// arguments after parsing options.
    pub fn rest(&self) -> &'a [OsString] {
        &self.argv[self.cursor..]
    }

    fn lookup_short(&self, ch: char) -> Option<ShortOpt> {
        self.short.iter().copied().find(|s| s.ch == ch)
    }

    fn lookup_long(&self, name: &str) -> Option<&'a LongOpt> {
        self.long.iter().find(|l| l.name == name)
    }
}

fn parse_short_spec(spec: &str) -> Vec<ShortOpt> {
    let mut out = Vec::new();
    let mut chars = spec.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch == ':' || ch == '+' || ch == '-' {
            continue;
        }
        let has_arg = if chars.peek() == Some(&':') {
            chars.next();
            if chars.peek() == Some(&':') {
                chars.next();
                HasArg::Optional
            } else {
                HasArg::Required
            }
        } else {
            HasArg::None
        };
        out.push(ShortOpt { ch, has_arg });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(s: &[&str]) -> Vec<OsString> {
        s.iter().map(|x| OsString::from(*x)).collect()
    }

    #[test]
    fn simple_required_short() {
        let argv = argv(&["head", "-n", "5", "in.vcf"]);
        let mut g = Getopt::new("n:", &[], &argv);
        let m = g.next().unwrap().unwrap();
        assert_eq!(m.code, b'n' as i32);
        assert_eq!(m.value.as_deref(), Some("5"));
        assert_eq!(g.next().unwrap(), None);
        assert_eq!(g.rest(), [OsString::from("in.vcf")]);
    }

    #[test]
    fn attached_short_value() {
        let argv = argv(&["view", "-Oz", "in.vcf"]);
        let mut g = Getopt::new("O:", &[], &argv);
        let m = g.next().unwrap().unwrap();
        assert_eq!(m.code, b'O' as i32);
        assert_eq!(m.value.as_deref(), Some("z"));
    }

    #[test]
    fn long_with_eq() {
        let argv = argv(&["index", "--threads=4", "in.bcf"]);
        let long = [LongOpt::new("threads", HasArg::Required, 9)];
        let mut g = Getopt::new("", &long, &argv);
        let m = g.next().unwrap().unwrap();
        assert_eq!(m.code, 9);
        assert_eq!(m.value.as_deref(), Some("4"));
    }

    #[test]
    fn long_with_space() {
        let argv = argv(&["index", "--min-shift", "14", "in.bcf"]);
        let long = [LongOpt::new("min-shift", HasArg::Required, b'm' as i32)];
        let mut g = Getopt::new("m:", &long, &argv);
        let m = g.next().unwrap().unwrap();
        assert_eq!(m.code, b'm' as i32);
        assert_eq!(m.value.as_deref(), Some("14"));
    }

    #[test]
    fn bundled_booleans() {
        let argv = argv(&["index", "-cf", "in.bcf"]);
        let mut g = Getopt::new("cf", &[], &argv);
        let m1 = g.next().unwrap().unwrap();
        let m2 = g.next().unwrap().unwrap();
        assert_eq!(m1.code, b'c' as i32);
        assert_eq!(m2.code, b'f' as i32);
        assert_eq!(g.next().unwrap(), None);
    }

    #[test]
    fn missing_required_arg() {
        let argv = argv(&["head", "-n"]);
        let mut g = Getopt::new("n:", &[], &argv);
        assert!(g.next().is_err());
    }

    #[test]
    fn unknown_option() {
        let argv = argv(&["head", "-Z"]);
        let mut g = Getopt::new("n:", &[], &argv);
        assert!(g.next().is_err());
    }

    #[test]
    fn double_dash_terminates() {
        let argv = argv(&["view", "-Oz", "--", "-not-a-flag.vcf"]);
        let mut g = Getopt::new("O:", &[], &argv);
        let _ = g.next().unwrap().unwrap();
        assert_eq!(g.next().unwrap(), None);
        assert_eq!(g.rest(), [OsString::from("-not-a-flag.vcf")]);
    }

    #[test]
    fn long_only_codes() {
        let argv = argv(&["index", "--threads", "8"]);
        let long = [LongOpt::new("threads", HasArg::Required, 9)];
        let mut g = Getopt::new("", &long, &argv);
        let m = g.next().unwrap().unwrap();
        assert_eq!(m.code, 9);
        assert_eq!(m.value.as_deref(), Some("8"));
    }
}
