//! Lexer foundation for the bcftools filter-expression engine.
//!
//! bcftools' `filter.c` is separate from HTSlib's lighter expression parser.
//! This module starts the native port with a token stream for record/sample
//! expressions used by `-i` and `-e`.

use std::collections::BTreeMap;
use std::io;

use regex::Regex;

#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    Identifier(String),
    Number(String),
    String(String),
    Operator(Operator),
    LeftParen,
    RightParen,
    LeftBracket,
    RightBracket,
    Comma,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Expr {
    Identifier(String),
    Number(String),
    String(String),
    Unary {
        op: UnaryOp,
        expr: Box<Expr>,
    },
    Binary {
        op: BinaryOp,
        lhs: Box<Expr>,
        rhs: Box<Expr>,
    },
    Call {
        function: String,
        args: Vec<Expr>,
    },
    Index {
        expr: Box<Expr>,
        index: Box<Expr>,
    },
    Wildcard,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Missing,
    Bool(bool),
    Number(f64),
    String(String),
    List(Vec<Value>),
}

impl Value {
    pub fn truthy(&self) -> bool {
        match self {
            Self::Missing => false,
            Self::Bool(value) => *value,
            Self::Number(value) => *value != 0.0 && !value.is_nan(),
            Self::String(value) => !value.is_empty() && value != ".",
            Self::List(values) => values.iter().any(Self::truthy),
        }
    }

    fn as_number(&self) -> Option<f64> {
        match self {
            Self::Number(value) => Some(*value),
            Self::Bool(value) => Some(if *value { 1.0 } else { 0.0 }),
            Self::String(value) => value.parse().ok(),
            _ => None,
        }
    }

    fn scalar_eq(&self, other: &Self) -> bool {
        match (self.as_number(), other.as_number()) {
            (Some(lhs), Some(rhs)) => lhs == rhs,
            _ => self.as_string() == other.as_string(),
        }
    }

    fn as_string(&self) -> String {
        match self {
            Self::Missing => ".".into(),
            Self::Bool(value) => value.to_string(),
            Self::Number(value) => {
                if value.fract() == 0.0 {
                    format!("{value:.0}")
                } else {
                    value.to_string()
                }
            }
            Self::String(value) => value.clone(),
            Self::List(values) => values
                .iter()
                .map(Self::as_string)
                .collect::<Vec<_>>()
                .join(","),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct EvalContext {
    values: BTreeMap<String, Value>,
    sample_values: Vec<BTreeMap<String, Value>>,
    active_sample: Option<usize>,
}

impl EvalContext {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with(mut self, name: impl Into<String>, value: Value) -> Self {
        self.values.insert(name.into(), value);
        self
    }

    pub fn with_sample(
        mut self,
        values: impl IntoIterator<Item = (impl Into<String>, Value)>,
    ) -> Self {
        self.sample_values.push(
            values
                .into_iter()
                .map(|(name, value)| (name.into(), value))
                .collect(),
        );
        self
    }

    pub fn get(&self, name: &str) -> Value {
        if let Some(sample) = self.active_sample
            && let Some(value) = self.sample_value(sample, name)
        {
            return value;
        }

        self.values.get(name).cloned().unwrap_or(Value::Missing)
    }

    pub fn sample_count(&self) -> usize {
        self.sample_values.len()
    }

    fn sample_value(&self, sample: usize, name: &str) -> Option<Value> {
        let values = self.sample_values.get(sample)?;
        values
            .get(name)
            .or_else(|| name.strip_prefix("FMT/").and_then(|key| values.get(key)))
            .or_else(|| name.strip_prefix("FORMAT/").and_then(|key| values.get(key)))
            .cloned()
    }

    fn for_sample(&self, sample: usize) -> Self {
        let mut context = self.clone();
        context.active_sample = Some(sample);
        context
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvalTrace {
    pub status: EvalStatus,
    pub lookups: Vec<LookupRequest>,
    pub short_circuits: usize,
}

impl Default for EvalTrace {
    fn default() -> Self {
        Self {
            status: EvalStatus::NotStarted,
            lookups: Vec::new(),
            short_circuits: 0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvalStatus {
    NotStarted,
    Completed,
    Error,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LookupRequest {
    pub name: String,
    pub sample_index: Option<usize>,
    pub source: LookupSource,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LookupSource {
    SampleContext,
    RecordContext,
    External,
    Missing,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnaryOp {
    Not,
    Neg,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BinaryOp {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Regex,
    NotRegex,
    And,
    Or,
    Add,
    Sub,
    Mul,
    Div,
    Mod,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operator {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Regex,
    NotRegex,
    And,
    Or,
    Not,
    Plus,
    Minus,
    Star,
    Slash,
    Percent,
}

pub fn lex(expression: &str) -> io::Result<Vec<Token>> {
    let mut chars = expression.char_indices().peekable();
    let mut tokens = Vec::new();

    while let Some((i, ch)) = chars.next() {
        match ch {
            ch if ch.is_whitespace() => {}
            '(' => tokens.push(Token::LeftParen),
            ')' => tokens.push(Token::RightParen),
            '[' => tokens.push(Token::LeftBracket),
            ']' => tokens.push(Token::RightBracket),
            ',' => tokens.push(Token::Comma),
            '"' | '\'' => tokens.push(Token::String(read_quoted(expression, &mut chars, ch)?)),
            ch if is_identifier_start(ch) => {
                tokens.push(Token::Identifier(read_identifier(
                    expression, &mut chars, i,
                )));
            }
            ch if ch.is_ascii_digit()
                || (ch == '.' && chars.peek().is_some_and(|&(_, next)| next.is_ascii_digit())) =>
            {
                tokens.push(Token::Number(read_number(expression, &mut chars, i)));
            }
            '=' => match chars.peek().copied() {
                Some((_, '=')) => {
                    chars.next();
                    tokens.push(Token::Operator(Operator::Eq));
                }
                Some((_, '~')) => {
                    chars.next();
                    tokens.push(Token::Operator(Operator::Regex));
                }
                _ => tokens.push(Token::Operator(Operator::Eq)),
            },
            '!' => match chars.peek().copied() {
                Some((_, '=')) => {
                    chars.next();
                    tokens.push(Token::Operator(Operator::Ne));
                }
                Some((_, '~')) => {
                    chars.next();
                    tokens.push(Token::Operator(Operator::NotRegex));
                }
                _ => tokens.push(Token::Operator(Operator::Not)),
            },
            '<' => match chars.peek().copied() {
                Some((_, '=')) => {
                    chars.next();
                    tokens.push(Token::Operator(Operator::Le));
                }
                _ => tokens.push(Token::Operator(Operator::Lt)),
            },
            '>' => match chars.peek().copied() {
                Some((_, '=')) => {
                    chars.next();
                    tokens.push(Token::Operator(Operator::Ge));
                }
                _ => tokens.push(Token::Operator(Operator::Gt)),
            },
            '&' if chars.peek().is_some_and(|&(_, next)| next == '&') => {
                chars.next();
                tokens.push(Token::Operator(Operator::And));
            }
            '|' if chars.peek().is_some_and(|&(_, next)| next == '|') => {
                chars.next();
                tokens.push(Token::Operator(Operator::Or));
            }
            '+' => tokens.push(Token::Operator(Operator::Plus)),
            '-' => tokens.push(Token::Operator(Operator::Minus)),
            '~' => tokens.push(Token::Operator(Operator::Regex)),
            '*' => tokens.push(Token::Operator(Operator::Star)),
            '/' => tokens.push(Token::Operator(Operator::Slash)),
            '%' => tokens.push(Token::Operator(Operator::Percent)),
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unexpected filter-expression character {ch:?}"),
                ));
            }
        }
    }

    Ok(tokens)
}

pub fn parse_expression(expression: &str) -> io::Result<Expr> {
    let tokens = lex(expression)?;
    let mut parser = Parser { tokens, pos: 0 };
    let expr = parser.parse_expr(0)?;

    if parser.peek().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "unexpected trailing token in filter expression",
        ));
    }

    Ok(expr)
}

pub fn eval_expression(expression: &str, context: &EvalContext) -> io::Result<Value> {
    let expr = parse_expression(expression)?;
    eval(&expr, context)
}

pub fn eval_expression_with(
    expression: &str,
    context: &EvalContext,
    mut resolver: impl FnMut(&str, Option<usize>) -> Option<Value>,
) -> io::Result<Value> {
    let expr = parse_expression(expression)?;
    eval_with(&expr, context, &mut resolver)
}

pub fn eval_expression_traced(
    expression: &str,
    context: &EvalContext,
    mut resolver: impl FnMut(&str, Option<usize>) -> Option<Value>,
) -> (io::Result<Value>, EvalTrace) {
    let mut trace = EvalTrace {
        status: EvalStatus::NotStarted,
        ..EvalTrace::default()
    };
    let result = parse_expression(expression)
        .and_then(|expr| eval_with_trace(&expr, context, &mut resolver, Some(&mut trace)));
    trace.status = if result.is_ok() {
        EvalStatus::Completed
    } else {
        EvalStatus::Error
    };
    (result, trace)
}

pub fn eval(expr: &Expr, context: &EvalContext) -> io::Result<Value> {
    eval_with(expr, context, &mut |_, _| None)
}

pub fn eval_with(
    expr: &Expr,
    context: &EvalContext,
    resolver: &mut impl FnMut(&str, Option<usize>) -> Option<Value>,
) -> io::Result<Value> {
    eval_with_trace(expr, context, resolver, None)
}

fn eval_with_trace(
    expr: &Expr,
    context: &EvalContext,
    resolver: &mut impl FnMut(&str, Option<usize>) -> Option<Value>,
    mut trace: Option<&mut EvalTrace>,
) -> io::Result<Value> {
    match expr {
        Expr::Identifier(name) => Ok(resolve_identifier(
            name,
            context,
            resolver,
            trace.as_deref_mut(),
        )),
        Expr::Number(value) => value
            .parse::<f64>()
            .map(Value::Number)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e)),
        Expr::String(value) => Ok(Value::String(value.clone())),
        Expr::Wildcard => Ok(Value::String("*".into())),
        Expr::Unary { op, expr } => {
            let value = eval_with_trace(expr, context, resolver, trace.as_deref_mut())?;
            match op {
                UnaryOp::Not => Ok(Value::Bool(!value.truthy())),
                UnaryOp::Neg => Ok(Value::Number(-number(&value)?)),
            }
        }
        Expr::Binary { op, lhs, rhs } => {
            eval_binary(*op, lhs, rhs, context, resolver, trace.as_deref_mut())
        }
        Expr::Call { function, args } => {
            eval_call(function, args, context, resolver, trace.as_deref_mut())
        }
        Expr::Index { expr, index } => {
            let value = eval_with_trace(expr, context, resolver, trace.as_deref_mut())?;
            let index = eval_with_trace(index, context, resolver, trace)?;
            match (value, index) {
                (Value::List(values), Value::String(index)) if index == "*" => {
                    Ok(Value::List(values))
                }
                (Value::List(values), index) => {
                    let index = number(&index)? as usize;
                    Ok(values.get(index).cloned().unwrap_or(Value::Missing))
                }
                (value, Value::String(index)) if index == "*" => Ok(value),
                _ => Ok(Value::Missing),
            }
        }
    }
}

fn resolve_identifier(
    name: &str,
    context: &EvalContext,
    resolver: &mut impl FnMut(&str, Option<usize>) -> Option<Value>,
    trace: Option<&mut EvalTrace>,
) -> Value {
    let mut source = LookupSource::Missing;
    let mut value = None;

    if let Some(sample) = context.active_sample
        && let Some(value) = context.sample_value(sample, name)
    {
        source = LookupSource::SampleContext;
        if let Some(trace) = trace {
            trace.lookups.push(LookupRequest {
                name: name.to_string(),
                sample_index: context.active_sample,
                source,
            });
        }
        return value;
    }

    if let Some(value) = context.values.get(name).cloned() {
        source = LookupSource::RecordContext;
        if let Some(trace) = trace {
            trace.lookups.push(LookupRequest {
                name: name.to_string(),
                sample_index: context.active_sample,
                source,
            });
        }
        return value;
    }

    if let Some(sample) = context.active_sample
        && let Some(value) = resolver(name, Some(sample))
    {
        source = LookupSource::External;
        if let Some(trace) = trace {
            trace.lookups.push(LookupRequest {
                name: name.to_string(),
                sample_index: Some(sample),
                source,
            });
        }
        return value;
    }

    if let Some(external) = resolver(name, None) {
        source = LookupSource::External;
        value = Some(external);
    }

    if let Some(trace) = trace {
        trace.lookups.push(LookupRequest {
            name: name.to_string(),
            sample_index: context.active_sample,
            source,
        });
    }

    value.unwrap_or(Value::Missing)
}

fn eval_binary(
    op: BinaryOp,
    lhs: &Expr,
    rhs: &Expr,
    context: &EvalContext,
    resolver: &mut impl FnMut(&str, Option<usize>) -> Option<Value>,
    mut trace: Option<&mut EvalTrace>,
) -> io::Result<Value> {
    if op == BinaryOp::And {
        let lhs = eval_with_trace(lhs, context, resolver, trace.as_deref_mut())?;
        if !lhs.truthy() {
            if let Some(trace) = trace {
                trace.short_circuits += 1;
            }
            return Ok(Value::Bool(false));
        }
        return Ok(Value::Bool(
            eval_with_trace(rhs, context, resolver, trace.as_deref_mut())?.truthy(),
        ));
    }

    if op == BinaryOp::Or {
        let lhs = eval_with_trace(lhs, context, resolver, trace.as_deref_mut())?;
        if lhs.truthy() {
            if let Some(trace) = trace {
                trace.short_circuits += 1;
            }
            return Ok(Value::Bool(true));
        }
        return Ok(Value::Bool(
            eval_with_trace(rhs, context, resolver, trace.as_deref_mut())?.truthy(),
        ));
    }

    let lhs = eval_with_trace(lhs, context, resolver, trace.as_deref_mut())?;
    let rhs = eval_with_trace(rhs, context, resolver, trace)?;

    match op {
        BinaryOp::Eq => Ok(Value::Bool(compare_any(&lhs, &rhs, Value::scalar_eq))),
        BinaryOp::Ne => Ok(Value::Bool(!compare_any(&lhs, &rhs, Value::scalar_eq))),
        BinaryOp::Lt => compare_numbers(&lhs, &rhs, |lhs, rhs| lhs < rhs),
        BinaryOp::Le => compare_numbers(&lhs, &rhs, |lhs, rhs| lhs <= rhs),
        BinaryOp::Gt => compare_numbers(&lhs, &rhs, |lhs, rhs| lhs > rhs),
        BinaryOp::Ge => compare_numbers(&lhs, &rhs, |lhs, rhs| lhs >= rhs),
        BinaryOp::Regex => regex_match(&lhs, &rhs),
        BinaryOp::NotRegex => {
            let Value::Bool(matches) = regex_match(&lhs, &rhs)? else {
                unreachable!("regex_match returns bool")
            };
            Ok(Value::Bool(!matches))
        }
        BinaryOp::Add => Ok(Value::Number(number(&lhs)? + number(&rhs)?)),
        BinaryOp::Sub => Ok(Value::Number(number(&lhs)? - number(&rhs)?)),
        BinaryOp::Mul => Ok(Value::Number(number(&lhs)? * number(&rhs)?)),
        BinaryOp::Div => Ok(Value::Number(number(&lhs)? / number(&rhs)?)),
        BinaryOp::Mod => Ok(Value::Number(number(&lhs)? % number(&rhs)?)),
        BinaryOp::And | BinaryOp::Or => unreachable!("handled before eager evaluation"),
    }
}

fn eval_call(
    function: &str,
    args: &[Expr],
    context: &EvalContext,
    resolver: &mut impl FnMut(&str, Option<usize>) -> Option<Value>,
    mut trace: Option<&mut EvalTrace>,
) -> io::Result<Value> {
    match function.to_ascii_uppercase().as_str() {
        "COUNT" => {
            require_arity(function, args, 1)?;
            let value = eval_with_trace(&args[0], context, resolver, trace.as_deref_mut())?;
            Ok(Value::Number(count_truthy(&value) as f64))
        }
        "N_PASS" => {
            require_arity(function, args, 1)?;
            if context.sample_count() == 0 {
                let value = eval_with_trace(&args[0], context, resolver, trace.as_deref_mut())?;
                return Ok(Value::Number(count_truthy(&value) as f64));
            }

            let mut count = 0usize;
            for sample in 0..context.sample_count() {
                if eval_with_trace(
                    &args[0],
                    &context.for_sample(sample),
                    resolver,
                    trace.as_deref_mut(),
                )?
                .truthy()
                {
                    count += 1;
                }
            }
            Ok(Value::Number(count as f64))
        }
        "MIN" => {
            require_arity(function, args, 1)?;
            let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
            values
                .into_iter()
                .reduce(f64::min)
                .map(Value::Number)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "MIN requires numeric values")
                })
        }
        "SUM" | "SSUM" | "SMPL_SUM" => {
            require_arity(function, args, 1)?;
            let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
            Ok(Value::Number(values.iter().sum()))
        }
        "AVG" | "MEAN" | "SAVG" | "SMEAN" | "SMPL_AVG" | "SMPL_MEAN" => {
            require_arity(function, args, 1)?;
            let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
            if values.is_empty() {
                Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("{function} requires numeric values"),
                ))
            } else {
                Ok(Value::Number(
                    values.iter().sum::<f64>() / values.len() as f64,
                ))
            }
        }
        "STDEV" | "SSTDEV" | "SMPL_STDEV" => {
            require_arity(function, args, 1)?;
            let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
            if values.is_empty() {
                Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("{function} requires numeric values"),
                ))
            } else {
                Ok(Value::Number(sample_standard_deviation(&values)))
            }
        }
        "ABS" => {
            require_arity(function, args, 1)?;
            Ok(Value::Number(
                number(&eval_with_trace(&args[0], context, resolver, trace)?)?.abs(),
            ))
        }
        "PHRED" => {
            require_arity(function, args, 1)?;
            let value = number(&eval_with_trace(&args[0], context, resolver, trace)?)?;
            Ok(Value::Number(phred_score(value)))
        }
        "BINOM" => {
            if args.len() == 1 {
                let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
                if values.len() < 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "binom requires at least two numeric values",
                    ));
                }
                Ok(Value::Number(binom_two_sided(values[0], values[1], 0.5)))
            } else if args.len() == 2 {
                let lhs = number(&eval_with_trace(
                    &args[0],
                    context,
                    resolver,
                    trace.as_deref_mut(),
                )?)?;
                let rhs = number(&eval_with_trace(&args[1], context, resolver, trace)?)?;
                Ok(Value::Number(binom_two_sided(lhs, rhs, 0.5)))
            } else {
                Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "binom expects one or two argument(s)",
                ))
            }
        }
        "FISHER" => {
            if args.len() == 1 {
                let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
                if values.len() != 4 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "fisher with one argument requires four numeric values",
                    ));
                }
                Ok(Value::Number(fisher_two_sided(
                    values[0], values[1], values[2], values[3],
                )))
            } else if args.len() == 2 {
                let lhs = numeric_values(&eval_with_trace(
                    &args[0],
                    context,
                    resolver,
                    trace.as_deref_mut(),
                )?);
                let rhs = numeric_values(&eval_with_trace(&args[1], context, resolver, trace)?);
                if lhs.len() < 2 || rhs.len() < 2 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "fisher with two arguments requires two numeric values in each argument",
                    ));
                }
                Ok(Value::Number(fisher_two_sided(
                    lhs[0], rhs[0], lhs[1], rhs[1],
                )))
            } else {
                Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "fisher expects one or two argument(s)",
                ))
            }
        }
        "MAX" => {
            require_arity(function, args, 1)?;
            let values = numeric_values(&eval_with_trace(&args[0], context, resolver, trace)?);
            values
                .into_iter()
                .reduce(f64::max)
                .map(Value::Number)
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "MAX requires numeric values")
                })
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported filter function {function}"),
        )),
    }
}

fn sample_standard_deviation(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }

    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / values.len() as f64;
    variance.sqrt()
}

fn binom_two_sided(a: f64, b: f64, probability: f64) -> f64 {
    if a < 0.0 || b < 0.0 {
        return 0.0;
    }
    let a = a.round() as i32;
    let b = b.round() as i32;
    if a == 0 && b == 0 {
        return -1.0;
    }
    if a == b {
        return 1.0;
    }

    let n = (a + b) as usize;
    let limit = a.min(b) as usize;
    let mut term = (1.0 - probability).powi(n as i32);
    let mut cdf = term;

    for k in 0..limit {
        term *= (n - k) as f64 / (k + 1) as f64 * probability / (1.0 - probability);
        cdf += term;
    }

    (2.0 * cdf).min(1.0)
}

fn fisher_two_sided(n11: f64, n12: f64, n21: f64, n22: f64) -> f64 {
    if [n11, n12, n21, n22].iter().any(|value| *value < 0.0) {
        return 0.0;
    }

    let n11 = n11.round() as i32;
    let n12 = n12.round() as i32;
    let n21 = n21.round() as i32;
    let n22 = n22.round() as i32;
    let row1 = n11 + n12;
    let row2 = n21 + n22;
    let col1 = n11 + n21;
    let total = row1 + row2;
    if total == 0 {
        return 1.0;
    }

    let min_a = 0.max(col1 - row2);
    let max_a = row1.min(col1);
    let observed = hypergeom_table_probability(n11, row1, col1, total);
    let epsilon = observed * 1e-12 + 1e-15;

    (min_a..=max_a)
        .map(|a| hypergeom_table_probability(a, row1, col1, total))
        .filter(|probability| *probability <= observed + epsilon)
        .sum::<f64>()
        .min(1.0)
}

fn hypergeom_table_probability(a: i32, row1: i32, col1: i32, total: i32) -> f64 {
    let b = row1 - a;
    let c = col1 - a;
    let d = total - row1 - c;
    (ln_factorial(row1)
        + ln_factorial(total - row1)
        + ln_factorial(col1)
        + ln_factorial(total - col1)
        - ln_factorial(total)
        - ln_factorial(a)
        - ln_factorial(b)
        - ln_factorial(c)
        - ln_factorial(d))
    .exp()
}

fn ln_factorial(value: i32) -> f64 {
    (1..=value).map(|i| (i as f64).ln()).sum()
}

fn phred_score(probability: f64) -> f64 {
    if probability == 0.0 {
        return 99.0;
    }

    (-4.3429 * probability.ln()).min(99.0)
}

fn require_arity(function: &str, args: &[Expr], expected: usize) -> io::Result<()> {
    if args.len() == expected {
        Ok(())
    } else {
        Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("{function} expects {expected} argument(s)"),
        ))
    }
}

fn count_truthy(value: &Value) -> usize {
    match value {
        Value::List(values) => values.iter().filter(|value| value.truthy()).count(),
        value => usize::from(value.truthy()),
    }
}

fn numeric_values(value: &Value) -> Vec<f64> {
    match value {
        Value::List(values) => values.iter().flat_map(numeric_values).collect(),
        value => value.as_number().into_iter().collect(),
    }
}

fn compare_numbers(
    lhs: &Value,
    rhs: &Value,
    compare: impl Fn(f64, f64) -> bool + Copy,
) -> io::Result<Value> {
    Ok(Value::Bool(compare_any(lhs, rhs, |lhs, rhs| {
        match (lhs.as_number(), rhs.as_number()) {
            (Some(lhs), Some(rhs)) => compare(lhs, rhs),
            _ => false,
        }
    })))
}

fn compare_any(lhs: &Value, rhs: &Value, compare: impl Fn(&Value, &Value) -> bool + Copy) -> bool {
    match (lhs, rhs) {
        (Value::List(lhs), Value::List(rhs)) => lhs
            .iter()
            .any(|lhs| rhs.iter().any(|rhs| compare(lhs, rhs))),
        (Value::List(values), rhs) => values.iter().any(|lhs| compare(lhs, rhs)),
        (lhs, Value::List(values)) => values.iter().any(|rhs| compare(lhs, rhs)),
        _ => compare(lhs, rhs),
    }
}

fn regex_match(lhs: &Value, rhs: &Value) -> io::Result<Value> {
    let pattern = Regex::new(rhs.as_string().as_str())
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;

    Ok(Value::Bool(match lhs {
        Value::List(values) => values
            .iter()
            .any(|value| pattern.is_match(value.as_string().as_str())),
        value => pattern.is_match(value.as_string().as_str()),
    }))
}

fn number(value: &Value) -> io::Result<f64> {
    value.as_number().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("expected numeric filter value, got {value:?}"),
        )
    })
}

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn parse_expr(&mut self, min_bp: u8) -> io::Result<Expr> {
        let mut lhs = self.parse_prefix()?;

        loop {
            if self.consume(&Token::LeftBracket) {
                let index = self.parse_expr(0)?;
                self.expect(Token::RightBracket)?;
                lhs = Expr::Index {
                    expr: Box::new(lhs),
                    index: Box::new(index),
                };
                continue;
            }

            let Some(op) = self.peek_binary_op() else {
                break;
            };
            let (left_bp, right_bp) = infix_binding_power(op);
            if left_bp < min_bp {
                break;
            }

            self.pos += 1;
            let rhs = self.parse_expr(right_bp)?;
            lhs = Expr::Binary {
                op,
                lhs: Box::new(lhs),
                rhs: Box::new(rhs),
            };
        }

        Ok(lhs)
    }

    fn parse_prefix(&mut self) -> io::Result<Expr> {
        match self.next() {
            Some(Token::Identifier(name)) => {
                if self.consume(&Token::LeftParen) {
                    let mut args = Vec::new();
                    if !self.consume(&Token::RightParen) {
                        loop {
                            args.push(self.parse_expr(0)?);
                            if self.consume(&Token::RightParen) {
                                break;
                            }
                            self.expect(Token::Comma)?;
                        }
                    }
                    Ok(Expr::Call {
                        function: name,
                        args,
                    })
                } else {
                    Ok(Expr::Identifier(name))
                }
            }
            Some(Token::Number(value)) => Ok(Expr::Number(value)),
            Some(Token::String(value)) => Ok(Expr::String(value)),
            Some(Token::Operator(Operator::Star)) => Ok(Expr::Wildcard),
            Some(Token::Operator(Operator::Not)) => Ok(Expr::Unary {
                op: UnaryOp::Not,
                expr: Box::new(self.parse_expr(13)?),
            }),
            Some(Token::Operator(Operator::Minus)) => Ok(Expr::Unary {
                op: UnaryOp::Neg,
                expr: Box::new(self.parse_expr(13)?),
            }),
            Some(Token::LeftParen) => {
                let expr = self.parse_expr(0)?;
                self.expect(Token::RightParen)?;
                Ok(expr)
            }
            Some(token) => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unexpected filter-expression token {token:?}"),
            )),
            None => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "expected filter-expression operand",
            )),
        }
    }

    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<Token> {
        let token = self.tokens.get(self.pos).cloned();
        if token.is_some() {
            self.pos += 1;
        }
        token
    }

    fn consume(&mut self, expected: &Token) -> bool {
        if self.peek() == Some(expected) {
            self.pos += 1;
            true
        } else {
            false
        }
    }

    fn expect(&mut self, expected: Token) -> io::Result<()> {
        if self.consume(&expected) {
            Ok(())
        } else {
            Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("expected filter-expression token {expected:?}"),
            ))
        }
    }

    fn peek_binary_op(&self) -> Option<BinaryOp> {
        match self.peek()? {
            Token::Operator(Operator::Eq) => Some(BinaryOp::Eq),
            Token::Operator(Operator::Ne) => Some(BinaryOp::Ne),
            Token::Operator(Operator::Lt) => Some(BinaryOp::Lt),
            Token::Operator(Operator::Le) => Some(BinaryOp::Le),
            Token::Operator(Operator::Gt) => Some(BinaryOp::Gt),
            Token::Operator(Operator::Ge) => Some(BinaryOp::Ge),
            Token::Operator(Operator::Regex) => Some(BinaryOp::Regex),
            Token::Operator(Operator::NotRegex) => Some(BinaryOp::NotRegex),
            Token::Operator(Operator::And) => Some(BinaryOp::And),
            Token::Operator(Operator::Or) => Some(BinaryOp::Or),
            Token::Operator(Operator::Plus) => Some(BinaryOp::Add),
            Token::Operator(Operator::Minus) => Some(BinaryOp::Sub),
            Token::Operator(Operator::Star) => Some(BinaryOp::Mul),
            Token::Operator(Operator::Slash) => Some(BinaryOp::Div),
            Token::Operator(Operator::Percent) => Some(BinaryOp::Mod),
            _ => None,
        }
    }
}

fn infix_binding_power(op: BinaryOp) -> (u8, u8) {
    match op {
        BinaryOp::Or => (1, 2),
        BinaryOp::And => (3, 4),
        BinaryOp::Eq
        | BinaryOp::Ne
        | BinaryOp::Lt
        | BinaryOp::Le
        | BinaryOp::Gt
        | BinaryOp::Ge
        | BinaryOp::Regex
        | BinaryOp::NotRegex => (5, 6),
        BinaryOp::Add | BinaryOp::Sub => (7, 8),
        BinaryOp::Mul | BinaryOp::Div | BinaryOp::Mod => (9, 10),
    }
}

fn read_identifier(
    expression: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    start: usize,
) -> String {
    let mut end = start;

    while let Some(&(i, ch)) = chars.peek() {
        if is_identifier_continue(ch) {
            end = i + ch.len_utf8();
            chars.next();
        } else {
            break;
        }
    }

    if end == start {
        end = start
            + expression[start..]
                .chars()
                .next()
                .expect("start char")
                .len_utf8();
    }

    expression[start..end].to_string()
}

fn read_number(
    expression: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    start: usize,
) -> String {
    let mut end = start;
    let mut allow_sign = false;

    while let Some(&(i, ch)) = chars.peek() {
        if ch.is_ascii_digit() || ch == '.' {
            end = i + ch.len_utf8();
            chars.next();
            allow_sign = false;
        } else if matches!(ch, 'e' | 'E') {
            end = i + ch.len_utf8();
            chars.next();
            allow_sign = true;
        } else if allow_sign && matches!(ch, '+' | '-') {
            end = i + ch.len_utf8();
            chars.next();
            allow_sign = false;
        } else {
            break;
        }
    }

    if end == start {
        end = start
            + expression[start..]
                .chars()
                .next()
                .expect("start char")
                .len_utf8();
    }

    expression[start..end].to_string()
}

fn read_quoted(
    expression: &str,
    chars: &mut std::iter::Peekable<std::str::CharIndices<'_>>,
    quote: char,
) -> io::Result<String> {
    let mut out = String::new();

    while let Some((_, ch)) = chars.next() {
        if ch == quote {
            return Ok(out);
        }

        if ch == '\\' {
            let Some((_, escaped)) = chars.next() else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "trailing escape in quoted filter string",
                ));
            };
            out.push(match escaped {
                'n' => '\n',
                't' => '\t',
                'r' => '\r',
                other => other,
            });
        } else {
            out.push(ch);
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidInput,
        format!("unterminated quoted filter string in {expression:?}"),
    ))
}

fn is_identifier_start(ch: char) -> bool {
    ch.is_ascii_alphabetic() || ch == '_' || ch == '.'
}

fn is_identifier_continue(ch: char) -> bool {
    ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | '/' | ':' | '-')
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lexes_core_comparisons_and_boolean_operators() {
        let tokens = lex(r#"QUAL>=30 && FILTER="PASS""#).unwrap();

        assert_eq!(
            tokens,
            [
                Token::Identifier("QUAL".into()),
                Token::Operator(Operator::Ge),
                Token::Number("30".into()),
                Token::Operator(Operator::And),
                Token::Identifier("FILTER".into()),
                Token::Operator(Operator::Eq),
                Token::String("PASS".into()),
            ]
        );
    }

    #[test]
    fn lexes_format_paths_functions_and_regex() {
        let tokens = lex(r#"N_PASS(FMT/DP>10) > 0 || INFO/CSQ~"missense""#).unwrap();

        assert_eq!(
            tokens,
            [
                Token::Identifier("N_PASS".into()),
                Token::LeftParen,
                Token::Identifier("FMT/DP".into()),
                Token::Operator(Operator::Gt),
                Token::Number("10".into()),
                Token::RightParen,
                Token::Operator(Operator::Gt),
                Token::Number("0".into()),
                Token::Operator(Operator::Or),
                Token::Identifier("INFO/CSQ".into()),
                Token::Operator(Operator::Regex),
                Token::String("missense".into()),
            ]
        );
    }

    #[test]
    fn lexes_indexing_not_regex_and_missing_value() {
        let tokens = lex(r#"TAG[0]!="." && ALT[*]!~"^<""#).unwrap();

        assert_eq!(
            tokens,
            [
                Token::Identifier("TAG".into()),
                Token::LeftBracket,
                Token::Number("0".into()),
                Token::RightBracket,
                Token::Operator(Operator::Ne),
                Token::String(".".into()),
                Token::Operator(Operator::And),
                Token::Identifier("ALT".into()),
                Token::LeftBracket,
                Token::Operator(Operator::Star),
                Token::RightBracket,
                Token::Operator(Operator::NotRegex),
                Token::String("^<".into()),
            ]
        );
    }

    #[test]
    fn lexes_single_quoted_strings_and_escapes() {
        let tokens = lex("'a\\tb'").unwrap();

        assert_eq!(tokens, [Token::String("a\tb".into())]);
    }

    #[test]
    fn parses_precedence_for_boolean_and_arithmetic_expressions() {
        let expr = parse_expression("QUAL + 1 >= 30 && DP < 10").unwrap();

        assert_eq!(
            expr,
            Expr::Binary {
                op: BinaryOp::And,
                lhs: Box::new(Expr::Binary {
                    op: BinaryOp::Ge,
                    lhs: Box::new(Expr::Binary {
                        op: BinaryOp::Add,
                        lhs: Box::new(Expr::Identifier("QUAL".into())),
                        rhs: Box::new(Expr::Number("1".into())),
                    }),
                    rhs: Box::new(Expr::Number("30".into())),
                }),
                rhs: Box::new(Expr::Binary {
                    op: BinaryOp::Lt,
                    lhs: Box::new(Expr::Identifier("DP".into())),
                    rhs: Box::new(Expr::Number("10".into())),
                }),
            }
        );
    }

    #[test]
    fn parses_functions_indexes_regex_and_wildcards() {
        let expr = parse_expression(r#"N_PASS(FMT/DP[0] > 10) || ALT[*] ~ "^<""#).unwrap();

        assert_eq!(
            expr,
            Expr::Binary {
                op: BinaryOp::Or,
                lhs: Box::new(Expr::Call {
                    function: "N_PASS".into(),
                    args: vec![Expr::Binary {
                        op: BinaryOp::Gt,
                        lhs: Box::new(Expr::Index {
                            expr: Box::new(Expr::Identifier("FMT/DP".into())),
                            index: Box::new(Expr::Number("0".into())),
                        }),
                        rhs: Box::new(Expr::Number("10".into())),
                    }],
                }),
                rhs: Box::new(Expr::Binary {
                    op: BinaryOp::Regex,
                    lhs: Box::new(Expr::Index {
                        expr: Box::new(Expr::Identifier("ALT".into())),
                        index: Box::new(Expr::Wildcard),
                    }),
                    rhs: Box::new(Expr::String("^<".into())),
                }),
            }
        );
    }

    #[test]
    fn parses_unary_not_and_negation() {
        let expr = parse_expression("!(DP < -1.5e-2)").unwrap();

        assert_eq!(
            expr,
            Expr::Unary {
                op: UnaryOp::Not,
                expr: Box::new(Expr::Binary {
                    op: BinaryOp::Lt,
                    lhs: Box::new(Expr::Identifier("DP".into())),
                    rhs: Box::new(Expr::Unary {
                        op: UnaryOp::Neg,
                        expr: Box::new(Expr::Number("1.5e-2".into())),
                    }),
                }),
            }
        );
    }

    #[test]
    fn parse_rejects_trailing_or_unbalanced_tokens() {
        assert!(parse_expression("DP >").is_err());
        assert!(parse_expression("(DP > 1").is_err());
        assert!(parse_expression("DP > 1 2").is_err());
    }

    #[test]
    fn evaluates_scalar_boolean_and_arithmetic_expressions() {
        let context = EvalContext::new()
            .with("QUAL", Value::Number(31.0))
            .with("DP", Value::Number(8.0))
            .with("FILTER", Value::String("PASS".into()));

        assert_eq!(
            eval_expression(r#"QUAL + 1 >= 30 && DP < 10 && FILTER="PASS""#, &context).unwrap(),
            Value::Bool(true)
        );
    }

    #[test]
    fn evaluates_lists_indexes_and_count_functions() {
        let context = EvalContext::new().with(
            "FMT/DP",
            Value::List(vec![
                Value::Number(4.0),
                Value::Number(12.0),
                Value::Missing,
            ]),
        );

        assert_eq!(
            eval_expression("FMT/DP[1] > 10", &context).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            eval_expression("COUNT(FMT/DP[*] > 10)", &context).unwrap(),
            Value::Number(1.0)
        );
    }

    #[test]
    fn evaluates_n_pass_against_sample_contexts() {
        let context = EvalContext::new()
            .with("QUAL", Value::Number(60.0))
            .with_sample([
                ("GT", Value::String("0/1".into())),
                ("DP", Value::Number(12.0)),
            ])
            .with_sample([
                ("GT", Value::String("0/0".into())),
                ("DP", Value::Number(22.0)),
            ])
            .with_sample([
                ("GT", Value::String("./.".into())),
                ("DP", Value::Number(3.0)),
            ]);

        assert_eq!(
            eval_expression(r#"N_PASS(FMT/DP > 10)"#, &context).unwrap(),
            Value::Number(2.0)
        );
        assert_eq!(
            eval_expression(r#"N_PASS(GT != "0/0" && DP >= 10)"#, &context).unwrap(),
            Value::Number(1.0)
        );
        assert_eq!(
            eval_expression(r#"QUAL = 60 && N_PASS(FORMAT/DP > 10) = 2"#, &context).unwrap(),
            Value::Bool(true)
        );
    }

    #[test]
    fn evaluates_external_value_injection() {
        let context = EvalContext::new()
            .with("QUAL", Value::Number(60.0))
            .with_sample([("GT", Value::String("0/1".into()))])
            .with_sample([("GT", Value::String("0/0".into()))]);

        let value = eval_expression_with(
            r#"QUAL = 60 && EXT = "ok" && N_PASS(EXT_DP > 10 && GT != "0/0") = 1"#,
            &context,
            |name, sample| match (name, sample) {
                ("EXT", None) => Some(Value::String("ok".into())),
                ("EXT_DP", Some(0)) => Some(Value::Number(12.0)),
                ("EXT_DP", Some(1)) => Some(Value::Number(20.0)),
                _ => None,
            },
        )
        .unwrap();

        assert_eq!(value, Value::Bool(true));
    }

    #[test]
    fn local_context_values_take_precedence_over_external_values() {
        let context = EvalContext::new().with("DP", Value::Number(8.0));

        let value = eval_expression_with("DP = 8", &context, |name, _| {
            (name == "DP").then_some(Value::Number(99.0))
        })
        .unwrap();

        assert_eq!(value, Value::Bool(true));
    }

    #[test]
    fn evaluates_min_max_and_short_circuit_logic() {
        let context = EvalContext::new().with(
            "AD",
            Value::List(vec![Value::Number(3.0), Value::Number(9.0)]),
        );

        assert_eq!(
            eval_expression("MIN(AD) = 3 && MAX(AD) = 9", &context).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            eval_expression("0 && UNKNOWN > 1", &context).unwrap(),
            Value::Bool(false)
        );
    }

    #[test]
    fn evaluates_common_numeric_functions() {
        let context = EvalContext::new()
            .with(
                "AD",
                Value::List(vec![Value::Number(3.0), Value::Number(9.0)]),
            )
            .with("DELTA", Value::Number(-4.0))
            .with("PVALUE", Value::Number(0.01));

        assert_eq!(
            eval_expression("SMPL_SUM(AD)", &context).unwrap(),
            Value::Number(12.0)
        );
        assert_eq!(
            eval_expression("sAVG(AD)", &context).unwrap(),
            Value::Number(6.0)
        );
        assert_eq!(
            eval_expression("MEAN(AD)", &context).unwrap(),
            Value::Number(6.0)
        );
        assert_eq!(
            eval_expression("sSTDEV(AD)", &context).unwrap(),
            Value::Number(3.0)
        );
        assert_eq!(
            eval_expression("ABS(DELTA)", &context).unwrap(),
            Value::Number(4.0)
        );

        let Value::Number(phred) = eval_expression("PHRED(PVALUE)", &context).unwrap() else {
            panic!("PHRED should return a number");
        };
        assert!((phred - 20.0).abs() < 0.001);
    }

    #[test]
    fn evaluates_simple_binomial_tail_function() {
        let context = EvalContext::new().with(
            "AD",
            Value::List(vec![Value::Number(10.0), Value::Number(2.0)]),
        );

        let Value::Number(pvalue) = eval_expression("binom(AD)", &context).unwrap() else {
            panic!("binom should return a number");
        };
        assert!((pvalue - 0.03857421875).abs() < f64::EPSILON);

        let Value::Number(phred) = eval_expression("phred(binom(10,2))", &context).unwrap() else {
            panic!("phred(binom()) should return a number");
        };
        assert!((phred - 14.137028).abs() < 0.001);
    }

    #[test]
    fn evaluates_simple_fisher_exact_function() {
        let context = EvalContext::new()
            .with(
                "DP4",
                Value::List(vec![
                    Value::Number(1.0),
                    Value::Number(9.0),
                    Value::Number(11.0),
                    Value::Number(3.0),
                ]),
            )
            .with(
                "ADF",
                Value::List(vec![Value::Number(1.0), Value::Number(11.0)]),
            )
            .with(
                "ADR",
                Value::List(vec![Value::Number(9.0), Value::Number(3.0)]),
            );

        let Value::Number(pvalue) = eval_expression("fisher(DP4)", &context).unwrap() else {
            panic!("fisher should return a number");
        };
        assert!((pvalue - 0.0027594561852200836).abs() < 1e-12);

        let Value::Number(two_arg) = eval_expression("fisher(ADF,ADR)", &context).unwrap() else {
            panic!("fisher with two arguments should return a number");
        };
        assert!((two_arg - pvalue).abs() < 1e-12);
    }

    #[test]
    fn traced_eval_records_lookup_sources_and_short_circuits() {
        let context = EvalContext::new()
            .with("QUAL", Value::Number(60.0))
            .with_sample([("GT", Value::String("0/1".into()))])
            .with_sample([("GT", Value::String("0/0".into()))]);

        let (value, trace) = eval_expression_traced(
            r#"QUAL > 50 || MISSING = 1 && N_PASS(GT = "0/1" && EXT_DP > 10) > 0"#,
            &context,
            |name, sample| match (name, sample) {
                ("EXT_DP", Some(0)) => Some(Value::Number(12.0)),
                ("EXT_DP", Some(1)) => Some(Value::Number(3.0)),
                _ => None,
            },
        );

        assert_eq!(value.unwrap(), Value::Bool(true));
        assert_eq!(trace.status, EvalStatus::Completed);
        assert_eq!(trace.short_circuits, 1);
        assert_eq!(
            trace.lookups,
            [LookupRequest {
                name: "QUAL".into(),
                sample_index: None,
                source: LookupSource::RecordContext
            }]
        );
    }

    #[test]
    fn traced_eval_records_sample_external_and_missing_lookups() {
        let context = EvalContext::new()
            .with_sample([("GT", Value::String("0/1".into()))])
            .with_sample([("GT", Value::String("0/0".into()))]);

        let (value, trace) = eval_expression_traced(
            r#"N_PASS(GT = "0/1" && EXT_DP > 10 && UNKNOWN = 1)"#,
            &context,
            |name, sample| match (name, sample) {
                ("EXT_DP", Some(0)) => Some(Value::Number(12.0)),
                ("EXT_DP", Some(1)) => Some(Value::Number(3.0)),
                _ => None,
            },
        );

        assert_eq!(value.unwrap(), Value::Number(0.0));
        assert_eq!(trace.status, EvalStatus::Completed);
        assert!(trace.lookups.contains(&LookupRequest {
            name: "GT".into(),
            sample_index: Some(0),
            source: LookupSource::SampleContext
        }));
        assert!(trace.lookups.contains(&LookupRequest {
            name: "EXT_DP".into(),
            sample_index: Some(0),
            source: LookupSource::External
        }));
        assert!(trace.lookups.contains(&LookupRequest {
            name: "UNKNOWN".into(),
            sample_index: Some(0),
            source: LookupSource::Missing
        }));
        assert!(trace.short_circuits > 0);
    }

    #[test]
    fn traced_eval_marks_errors() {
        let (value, trace) =
            eval_expression_traced("MIN(UNKNOWN)", &EvalContext::new(), |_, _| None);

        assert!(value.is_err());
        assert_eq!(trace.status, EvalStatus::Error);
        assert_eq!(
            trace.lookups,
            [LookupRequest {
                name: "UNKNOWN".into(),
                sample_index: None,
                source: LookupSource::Missing
            }]
        );
    }

    #[test]
    fn evaluates_regex_with_real_patterns_and_lists() {
        let context = EvalContext::new()
            .with("INFO/CSQ", Value::String("missense_variant".into()))
            .with(
                "ALT",
                Value::List(vec![
                    Value::String("<DEL>".into()),
                    Value::String("A".into()),
                ]),
            );

        assert_eq!(
            eval_expression(r#"INFO/CSQ ~ "^missense_""#, &context).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            eval_expression(r#"INFO/CSQ !~ "synonymous""#, &context).unwrap(),
            Value::Bool(true)
        );
        assert_eq!(
            eval_expression(r#"ALT ~ "^<""#, &context).unwrap(),
            Value::Bool(true)
        );
        assert!(eval_expression(r#"ALT ~ "[""#, &context).is_err());
    }

    #[test]
    fn rejects_unexpected_or_unterminated_input() {
        assert!(lex("@").is_err());
        assert!(lex("\"unterminated").is_err());
        assert!(lex("'trailing\\").is_err());
    }
}
