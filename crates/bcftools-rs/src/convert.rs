//! Format-string parsing foundation for bcftools `convert.c`.
//!
//! Upstream `convert.c` implements the large `query -f` mini-language. This
//! module starts with a reusable syntax layer: split a format string into
//! literals, percent tokens, and sample loops while preserving upstream escape
//! behavior for common control characters.

use std::io;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FormatItem {
    Literal(String),
    Token(String),
    SampleLoop(Vec<FormatItem>),
}

pub fn parse_format(format: &str) -> io::Result<Vec<FormatItem>> {
    let mut chars = format.char_indices().peekable();
    parse_until(format, &mut chars, None)
}

pub trait FormatRecord {
    fn sample_count(&self) -> usize;
    fn value(&self, token: &str, sample_index: Option<usize>) -> Option<String>;
}

pub fn render_format(format: &str, record: &impl FormatRecord) -> io::Result<String> {
    let items = parse_format(format)?;
    render_items(&items, record, None)
}

pub fn render_items(
    items: &[FormatItem],
    record: &impl FormatRecord,
    sample_index: Option<usize>,
) -> io::Result<String> {
    let mut out = String::new();

    for item in items {
        match item {
            FormatItem::Literal(value) => out.push_str(value),
            FormatItem::Token(token) if token == "%" => out.push('%'),
            FormatItem::Token(token) => {
                out.push_str(render_token(token, record, sample_index)?.as_str());
            }
            FormatItem::SampleLoop(children) => {
                for i in 0..record.sample_count() {
                    out.push_str(render_items(children, record, Some(i))?.as_str());
                }
            }
        }
    }

    Ok(out)
}

fn render_token(
    token: &str,
    record: &impl FormatRecord,
    sample_index: Option<usize>,
) -> io::Result<String> {
    if let Some((function, argument)) = split_function(token) {
        return render_function(function, argument, record, sample_index);
    }

    if let Some((base, index)) = split_indexed_token(token) {
        let value = record
            .value(base, sample_index)
            .unwrap_or_else(|| ".".to_string());
        return Ok(value
            .split(',')
            .nth(index)
            .filter(|value| !value.is_empty())
            .unwrap_or(".")
            .to_string());
    }

    Ok(record
        .value(token, sample_index)
        .unwrap_or_else(|| ".".to_string()))
}

fn render_function(
    function: &str,
    argument: &str,
    record: &impl FormatRecord,
    sample_index: Option<usize>,
) -> io::Result<String> {
    let function_key = function.to_ascii_uppercase();
    let numbers = numeric_function_values(function, argument, record, sample_index)?;

    match function_key.as_str() {
        "SUM" | "SSUM" | "SMPL_SUM" => Ok(format_number(numbers.iter().sum())),
        "AVG" | "SAVG" | "SMPL_AVG" => {
            if numbers.is_empty() {
                Ok(".".into())
            } else {
                Ok(format_number(
                    numbers.iter().sum::<f64>() / numbers.len() as f64,
                ))
            }
        }
        "MIN" | "SMIN" | "SMPL_MIN" => Ok(numbers
            .into_iter()
            .reduce(f64::min)
            .map(format_number)
            .unwrap_or_else(|| ".".into())),
        "MAX" | "SMAX" | "SMPL_MAX" => Ok(numbers
            .into_iter()
            .reduce(f64::max)
            .map(format_number)
            .unwrap_or_else(|| ".".into())),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported format function %{function}"),
        )),
    }
}

fn numeric_function_values(
    function: &str,
    argument: &str,
    record: &impl FormatRecord,
    sample_index: Option<usize>,
) -> io::Result<Vec<f64>> {
    let argument = argument.trim();
    let rendered_values = if sample_index.is_none() && is_format_argument(argument) {
        (0..record.sample_count())
            .map(|i| render_token(argument, record, Some(i)))
            .collect::<io::Result<Vec<_>>>()?
    } else {
        vec![render_token(argument, record, sample_index)?]
    };

    rendered_values
        .iter()
        .flat_map(|values| values.split(','))
        .filter(|value| !value.is_empty() && *value != ".")
        .map(|value| {
            value.parse::<f64>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("non-numeric value {value:?} in %{function}"),
                )
            })
        })
        .collect()
}

fn is_format_argument(argument: &str) -> bool {
    argument.starts_with("FORMAT/") || argument.starts_with("FMT/")
}

fn split_function(token: &str) -> Option<(&str, &str)> {
    let open = token.find('(')?;
    token
        .ends_with(')')
        .then_some((&token[..open], &token[open + 1..token.len() - 1]))
        .filter(|(function, _)| !function.is_empty())
}

fn split_indexed_token(token: &str) -> Option<(&str, usize)> {
    let open = token.rfind('{')?;
    let index = token[open + 1..token.len() - 1].parse().ok()?;
    let base = &token[..open];
    (!base.is_empty()).then_some((base, index))
}

fn format_number(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{}", value as i64)
    } else {
        format!("{value}")
    }
}

fn parse_until(
    format: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    terminator: Option<char>,
) -> io::Result<Vec<FormatItem>> {
    let mut items = Vec::new();
    let mut literal = String::new();

    while let Some((_, ch)) = chars.next() {
        if Some(ch) == terminator {
            flush_literal(&mut items, &mut literal);
            return Ok(items);
        }

        match ch {
            '\\' => literal.push(parse_escape(chars)?),
            '%' => {
                flush_literal(&mut items, &mut literal);
                items.push(FormatItem::Token(parse_token(format, chars)?));
            }
            '[' => {
                flush_literal(&mut items, &mut literal);
                items.push(FormatItem::SampleLoop(parse_until(
                    format,
                    chars,
                    Some(']'),
                )?));
            }
            ']' if terminator.is_none() => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "unmatched closing bracket in format string",
                ));
            }
            _ => literal.push(ch),
        }
    }

    if terminator.is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unterminated sample loop in format string",
        ));
    }

    flush_literal(&mut items, &mut literal);
    Ok(items)
}

fn parse_escape(chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>) -> io::Result<char> {
    let Some((_, ch)) = chars.next() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "trailing backslash in format string",
        ));
    };

    Ok(match ch {
        'n' => '\n',
        't' => '\t',
        'r' => '\r',
        '\\' => '\\',
        other => other,
    })
}

fn parse_token(
    format: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
) -> io::Result<String> {
    let Some((start, ch)) = chars.next() else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "trailing percent in format string",
        ));
    };

    if ch == '%' {
        return Ok("%".into());
    }

    if ch == '{' {
        let token_start = start + ch.len_utf8();
        for (i, ch) in chars.by_ref() {
            if ch == '}' {
                return Ok(format[token_start..i].to_string());
            }
        }

        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unterminated braced token in format string",
        ));
    }

    let token_start = start;
    let mut token_end = start + ch.len_utf8();

    while let Some(&(i, next)) = chars.peek() {
        if is_token_char(next) {
            token_end = i + next.len_utf8();
            chars.next();
        } else {
            break;
        }
    }

    parse_token_suffix(format, chars, token_start, token_end)
}

fn is_token_char(ch: char) -> bool {
    ch.is_ascii_alphanumeric() || matches!(ch, '_' | '/' | '-' | '.')
}

fn parse_token_suffix(
    format: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    token_start: usize,
    mut token_end: usize,
) -> io::Result<String> {
    let Some(&(suffix_start, suffix_open)) = chars.peek() else {
        return Ok(format[token_start..token_end].to_string());
    };

    let suffix_close = match suffix_open {
        '{' => '}',
        '(' => ')',
        _ => return Ok(format[token_start..token_end].to_string()),
    };

    chars.next();
    let mut depth = 1usize;
    let mut in_string = false;

    for (idx, ch) in chars.by_ref() {
        match ch {
            '"' => in_string = !in_string,
            c if c == suffix_open && !in_string => depth += 1,
            c if c == suffix_close && !in_string => {
                depth -= 1;
                if depth == 0 {
                    token_end = idx + ch.len_utf8();
                    return Ok(format[token_start..token_end].to_string());
                }
            }
            _ => {}
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!(
            "unterminated {} suffix in format token",
            &format[suffix_start..suffix_start + suffix_open.len_utf8()]
        ),
    ))
}

fn flush_literal(items: &mut Vec<FormatItem>, literal: &mut String) {
    if !literal.is_empty() {
        items.push(FormatItem::Literal(std::mem::take(literal)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    #[derive(Default)]
    struct MockRecord {
        values: BTreeMap<String, String>,
        sample_values: Vec<BTreeMap<String, String>>,
    }

    impl MockRecord {
        fn with_value(mut self, key: &str, value: &str) -> Self {
            self.values.insert(key.into(), value.into());
            self
        }

        fn with_sample(mut self, values: &[(&str, &str)]) -> Self {
            self.sample_values.push(
                values
                    .iter()
                    .map(|(key, value)| ((*key).into(), (*value).into()))
                    .collect(),
            );
            self
        }
    }

    impl FormatRecord for MockRecord {
        fn sample_count(&self) -> usize {
            self.sample_values.len()
        }

        fn value(&self, token: &str, sample_index: Option<usize>) -> Option<String> {
            match sample_index {
                Some(i) if token == "SAMPLE" => Some(format!("S{}", i + 1)),
                Some(i) => self.sample_values.get(i).and_then(|values| {
                    values
                        .get(token)
                        .or_else(|| {
                            token
                                .strip_prefix("FORMAT/")
                                .and_then(|key| values.get(key))
                        })
                        .or_else(|| token.strip_prefix("FMT/").and_then(|key| values.get(key)))
                        .cloned()
                }),
                None => self.values.get(token).cloned(),
            }
        }
    }

    #[test]
    fn parses_literals_tokens_and_escapes() {
        let items = parse_format("%CHROM\\t%POS\\n").unwrap();

        assert_eq!(
            items,
            [
                FormatItem::Token("CHROM".into()),
                FormatItem::Literal("\t".into()),
                FormatItem::Token("POS".into()),
                FormatItem::Literal("\n".into()),
            ]
        );
    }

    #[test]
    fn parses_sample_loops() {
        let items = parse_format("%CHROM[\\t%SAMPLE=%GT]\\n").unwrap();

        assert_eq!(
            items,
            [
                FormatItem::Token("CHROM".into()),
                FormatItem::SampleLoop(vec![
                    FormatItem::Literal("\t".into()),
                    FormatItem::Token("SAMPLE".into()),
                    FormatItem::Literal("=".into()),
                    FormatItem::Token("GT".into()),
                ]),
                FormatItem::Literal("\n".into()),
            ]
        );
    }

    #[test]
    fn parses_braced_and_forced_namespace_tokens() {
        let items = parse_format("%{INFO/CSQ}\\t%/DP").unwrap();

        assert_eq!(
            items,
            [
                FormatItem::Token("INFO/CSQ".into()),
                FormatItem::Literal("\t".into()),
                FormatItem::Token("/DP".into()),
            ]
        );
    }

    #[test]
    fn parses_vector_indexes_and_function_tokens() {
        let items = parse_format("%AC{1}\\t%SUM(INFO/AD)\\n").unwrap();

        assert_eq!(
            items,
            [
                FormatItem::Token("AC{1}".into()),
                FormatItem::Literal("\t".into()),
                FormatItem::Token("SUM(INFO/AD)".into()),
                FormatItem::Literal("\n".into()),
            ]
        );
    }

    #[test]
    fn parses_literal_percent() {
        let items = parse_format("rate=%%\\n").unwrap();

        assert_eq!(
            items,
            [
                FormatItem::Literal("rate=".into()),
                FormatItem::Token("%".into()),
                FormatItem::Literal("\n".into()),
            ]
        );
    }

    #[test]
    fn rejects_unbalanced_constructs() {
        assert!(parse_format("%").is_err());
        assert!(parse_format("[%GT").is_err());
        assert!(parse_format("%{INFO").is_err());
        assert!(parse_format("%GT]").is_err());
    }

    #[test]
    fn renders_record_tokens_literals_and_missing_values() {
        let record = MockRecord::default()
            .with_value("CHROM", "chr1")
            .with_value("POS", "7");

        assert_eq!(
            render_format("%CHROM:%POS:%ID:%%\\n", &record).unwrap(),
            "chr1:7:.:%\n"
        );
    }

    #[test]
    fn renders_sample_loops() {
        let record = MockRecord::default()
            .with_value("CHROM", "chr1")
            .with_sample(&[("GT", "0/1"), ("DP", "5")])
            .with_sample(&[("GT", "1/1"), ("DP", "9")]);

        assert_eq!(
            render_format("%CHROM[\\t%SAMPLE=%GT:%DP]\\n", &record).unwrap(),
            "chr1\tS1=0/1:5\tS2=1/1:9\n"
        );
    }

    #[test]
    fn renders_vector_indexes_and_numeric_functions() {
        let record = MockRecord::default()
            .with_value("AC", "2,5,8")
            .with_value("INFO/AD", "3,4,5")
            .with_sample(&[("AD", "1,2"), ("DP", "7")])
            .with_sample(&[("AD", "4,6"), ("DP", "11")]);

        assert_eq!(
            render_format("%AC{1}\\t%SUM(INFO/AD)\\t%AVG(INFO/AD)\\n", &record).unwrap(),
            "5\t12\t4\n"
        );
        assert_eq!(
            render_format("[%SAMPLE:%AD{1}:%sSUM(AD)\\n]", &record).unwrap(),
            "S1:2:3\nS2:6:10\n"
        );
    }

    #[test]
    fn renders_case_insensitive_and_sample_vector_numeric_functions() {
        let record = MockRecord::default()
            .with_sample(&[("AD", "1,2"), ("DP", "7")])
            .with_sample(&[("AD", "4,6"), ("DP", "11")]);

        assert_eq!(
            render_format(
                "%sum(FORMAT/AD)\\t%SMPL_MAX(FMT/DP)\\t%smpl_avg(FORMAT/DP)\\n",
                &record,
            )
            .unwrap(),
            "13\t11\t9\n"
        );
    }

    #[test]
    fn renders_missing_for_out_of_range_vector_indexes() {
        let record = MockRecord::default().with_value("AC", "2");

        assert_eq!(render_format("%AC{3}\\n", &record).unwrap(), ".\n");
    }

    #[test]
    fn renders_parsed_item_trees_directly() {
        let record = MockRecord::default().with_value("INFO/CSQ", "missense");
        let items = parse_format("%{INFO/CSQ}\\n").unwrap();

        assert_eq!(render_items(&items, &record, None).unwrap(), "missense\n");
    }
}
