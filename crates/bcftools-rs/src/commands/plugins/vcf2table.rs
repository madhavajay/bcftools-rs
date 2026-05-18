//! `bcftools +vcf2table` (upstream `bcftools/plugins/vcf2table.c`).
//!
//! Renders each VCF record as a set of ASCII box tables: a `<<<`/`>>>`
//! delimited block containing the `# Variant`, `# INFO`,
//! `# GENOTYPE TYPES`, and `# GENOTYPES` tables. This is the filter-free
//! first slice: it targets the no-argument upstream fixture
//! (`vcf2table.1.out`) and the non-tty ASCII rendering path (upstream
//! forces `ascii=1` when stdout is not a tty, which is always the case
//! when output is captured). The shared HTSlib `kputd` formatter renders
//! the numeric IDX / percentage cells, exactly as upstream's
//! `CellSetD` -> `kputd` does.
//!
//! Deferred to later slices (tracked under "Remaining" in TODO.md):
//! VEP/CSQ, BCSQ, ANN/SNPEFF, LOF, SpliceAI and HYPERLINKS tables, the
//! Unicode/color (tty) rendering path, genome-build hyperlink
//! generation, and the `-x`/`--hide` option (its fixture,
//! `vcf2table.2.out`, needs the VEP/BCSQ tables).

use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};

use super::prune::kputd;
use crate::vcf_compat::normalize_vcf_text;

#[derive(Clone, Copy, PartialEq)]
enum Build {
    Undefined,
    Hg19,
    Hg38,
    RotavirusRf,
}

impl Build {
    /// Upstream `PRINT_HEADER` build prefix (note the surrounding spaces).
    fn prefix(self) -> &'static str {
        match self {
            Build::Hg19 => " GRCh37 : ",
            Build::Hg38 => " GRCh38 : ",
            Build::RotavirusRf => " Rotavirus : ",
            Build::Undefined => "",
        }
    }
}

/// A simple table mirroring upstream's `TablePtr`: a header row plus body
/// rows, rendered in ASCII mode with `+`/`-`/`|` borders and one space of
/// padding on each side of every cell.
struct Table {
    header: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl Table {
    fn new(header: &[&str]) -> Self {
        Self {
            header: header.iter().map(|s| (*s).to_owned()).collect(),
            rows: Vec::new(),
        }
    }

    fn ncols(&self) -> usize {
        self.header.len()
    }

    /// Push a row, padding/truncating to the column count (upstream rows
    /// are allocated with exactly `TableNCols` empty cells).
    fn push_row(&mut self, mut cells: Vec<String>) {
        cells.resize(self.ncols(), String::new());
        self.rows.push(cells);
    }

    /// Faithful port of upstream `TablePrint` in ASCII mode.
    fn render(&self, out: &mut String) {
        let ncols = self.ncols();
        let mut widths = vec![0usize; ncols];
        for (x, w) in widths.iter_mut().enumerate() {
            *w = self.header[x].len();
        }
        for row in &self.rows {
            for (x, cell) in row.iter().enumerate() {
                if cell.len() > widths[x] {
                    widths[x] = cell.len();
                }
            }
        }
        let empty = self.rows.is_empty();

        // line 1: top border
        for w in &widths {
            out.push('+');
            for _ in 0..(2 + w) {
                out.push('-');
            }
        }
        out.push_str("+\n");

        // line 2: header cells
        for (x, title) in self.header.iter().enumerate() {
            out.push_str("| ");
            out.push_str(title);
            for _ in 0..(widths[x] - title.len()) {
                out.push(' ');
            }
            out.push(' ');
        }
        out.push_str("|\n");

        // line 3: header/body separator (`+` in ASCII regardless of
        // upstream's ├┼┤ vs └┴┘ choice; the dash run is identical).
        for w in &widths {
            out.push('+');
            for _ in 0..(2 + w) {
                out.push('-');
            }
        }
        out.push_str("+\n");

        // body
        for row in &self.rows {
            for (x, cell) in row.iter().enumerate() {
                out.push_str("| ");
                out.push_str(cell);
                for _ in 0..(widths[x] - cell.len()) {
                    out.push(' ');
                }
                out.push(' ');
            }
            out.push_str("|\n");
        }

        // last line (only when there are body rows)
        if !empty {
            for w in &widths {
                out.push('+');
                for _ in 0..(2 + w) {
                    out.push('-');
                }
            }
            out.push_str("+\n");
        }
        out.push('\n');
    }
}

/// Reads the input VCF/BCF and returns the rendered table text.
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(render(&text))
}

fn render(text: &str) -> String {
    let mut samples: Vec<&str> = Vec::new();
    let build = detect_build(text);
    let mut out = String::with_capacity(text.len() * 4);
    let mut n_variants: u64 = 0;

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("#CHROM") {
            // Sample names are the columns after FORMAT (10th onward).
            let cols: Vec<&str> = rest.split('\t').collect();
            // `rest` begins right after "#CHROM"; full header is
            // CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <samples...>
            let full: Vec<&str> = line.split('\t').collect();
            if full.len() > 9 {
                samples = full[9..].to_vec();
            }
            let _ = cols;
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        n_variants += 1;
        render_record(line, &samples, build, n_variants, &mut out);
    }
    out
}

fn header_line(marker: &str, build: Build, t: &[&str], n: u64, out: &mut String) {
    out.push_str(marker);
    out.push_str(build.prefix());
    // Upstream: fprintf(" %s:%s:%s (n. %ld)\n", CHROM, POS, REF, n).
    let _ = writeln!(out, " {}:{}:{} (n. {n})", t[0], t[1], t[3]);
}

fn render_record(line: &str, samples: &[&str], build: Build, n: u64, out: &mut String) {
    let t: Vec<&str> = line.split('\t').collect();
    if t.len() < 8 {
        return;
    }

    header_line("<<<", build, &t, n, out);
    out.push('\n');

    // # Variant
    let mut vc = Table::new(&["KEY", "VALUE"]);
    vc.push_row(vec!["CHROM".into(), t[0].into()]);
    vc.push_row(vec!["POS".into(), t[1].into()]);
    // end/length only when the variant spans more than one base.
    let pos0: i64 = t[1].parse::<i64>().map(|p| p - 1).unwrap_or(0);
    let end1 = info_end(t[7]).unwrap_or_else(|| pos0 + t[3].len() as i64);
    if pos0 + 1 != end1 {
        vc.push_row(vec!["end".into(), end1.to_string()]);
        vc.push_row(vec!["length".into(), (end1 - pos0).to_string()]);
    }
    vc.push_row(vec!["ID".into(), t[2].into()]);
    vc.push_row(vec!["REF".into(), t[3].into()]);
    vc.push_row(vec!["ALT".into(), t[4].into()]);
    vc.push_row(vec!["QUAL".into(), t[5].into()]);
    vc.push_row(vec!["FILTER".into(), t[6].into()]);
    out.push_str("# Variant\n");
    vc.render(out);

    // # INFO
    if t.len() > 7 && t[7] != "." {
        let mut info = Table::new(&["KEY", "IDX", "VALUE"]);
        for entry in t[7].split(';') {
            let Some(eq) = entry.find('=') else {
                continue; // flag INFO (no '='): skipped upstream
            };
            if eq == 0 {
                continue;
            }
            let key = &entry[..eq];
            let vals: Vec<&str> = entry[eq + 1..].split(',').collect();
            for (j, v) in vals.iter().enumerate() {
                let idx = if vals.len() > 1 {
                    kputd((j + 1) as f64)
                } else {
                    String::new()
                };
                info.push_row(vec![key.to_owned(), idx, (*v).to_owned()]);
            }
        }
        if !info.rows.is_empty() {
            out.push_str("# INFO\n");
            info.render(out);
        }
    }

    // genotypes
    if t.len() > 9 {
        let formats: Vec<&str> = t[8].split(':').collect();
        let gt_col = formats.iter().position(|f| *f == "GT");
        let ft_col = formats.iter().position(|f| *f == "FT");
        let _ = ft_col; // FT coloring is tty-only (deferred)

        let mut gcols: Vec<&str> = vec!["SAMPLE", "GTYPE"];
        gcols.extend_from_slice(&formats);
        let mut gtable = Table::new(&gcols);

        let (mut c_ref, mut c_het, mut c_var, mut c_mis, mut c_other) = (0, 0, 0, 0, 0);

        for (si, raw) in t[9..].iter().enumerate() {
            // vcf_format pads short FORMAT samples to the full key count
            // with ".", which is what process() observes.
            let mut vals: Vec<String> = raw.split(':').map(|s| s.to_owned()).collect();
            if vals.len() < formats.len() {
                vals.resize(formats.len(), ".".to_owned());
            }

            let mut gtype = String::new();
            if let Some(gc) = gt_col
                && gc < vals.len()
            {
                let gt = vals[gc].replace('|', "/");
                let alleles: Vec<&str> = gt.split('/').collect();
                let (mut a0, mut a1, mut amiss, mut aother) = (0, 0, 0, 0);
                for a in &alleles {
                    match *a {
                        "0" => a0 += 1,
                        "1" => a1 += 1,
                        "." => amiss += 1,
                        _ => aother += 1,
                    }
                }
                match alleles.len() {
                    2 => {
                        if a0 == 0 && a1 == 0 && aother == 0 {
                            gtype = "NO_CALL".into();
                            c_mis += 1;
                        } else if a0 == 2 {
                            gtype = "HOM_REF".into();
                            c_ref += 1;
                        } else if amiss == 0 && alleles[0] == alleles[1] {
                            gtype = "HOM_VAR".into();
                            c_var += 1;
                        } else if amiss == 0 && alleles[0] != alleles[1] {
                            gtype = "HET".into();
                            c_het += 1;
                        } else {
                            c_other += 1;
                        }
                    }
                    1 => {
                        if a0 == 1 {
                            gtype = "REF".into();
                            c_ref += 1;
                        } else if a1 == 1 {
                            gtype = "ALT".into();
                            c_var += 1;
                        } else if amiss == 1 {
                            gtype = "NO_CALL".into();
                            c_mis += 1;
                        } else {
                            c_other += 1;
                        }
                    }
                    n => {
                        if a0 == n {
                            gtype = "HOM_REF".into();
                            c_ref += 1;
                        } else if a1 == n {
                            gtype = "HOM_VAR".into();
                            c_ref += 1; // upstream increments count_hom_ref here
                        } else if amiss == n {
                            gtype = "NO_CALL".into();
                            c_mis += 1;
                        } else {
                            c_other += 1;
                        }
                    }
                }
            }

            let mut row = vec![samples.get(si).copied().unwrap_or("").to_owned(), gtype];
            row.extend(vals);
            gtable.push_row(row);
        }

        // # GENOTYPE TYPES
        let total = c_ref + c_het + c_var + c_mis + c_other;
        let mut gt = Table::new(&["Type", "Count", "%"]);
        let add = |label: &str, count: i64, gt: &mut Table| {
            if count > 0 && total > 0 {
                gt.push_row(vec![
                    label.to_owned(),
                    count.to_string(),
                    kputd(100.0 * (count as f32 / total as f32) as f64),
                ]);
            }
        };
        add("REF only ", c_ref, &mut gt);
        add("HET", c_het, &mut gt);
        add("ALT only", c_var, &mut gt);
        add("MISSING", c_mis, &mut gt);
        add("OTHER", c_other, &mut gt);
        if !gt.rows.is_empty() {
            out.push_str("# GENOTYPE TYPES\n");
            gt.render(out);
        }

        if !gtable.rows.is_empty() {
            out.push_str("# GENOTYPES\n");
            gtable.render(out);
        }
    }

    header_line(">>>", build, &t, n, out);
    out.push('\n');
}

/// Parses `INFO/END=` (1-based); upstream uses it to set `rlen`.
fn info_end(info: &str) -> Option<i64> {
    if info == "." {
        return None;
    }
    info.split(';')
        .find_map(|kv| kv.strip_prefix("END="))
        .and_then(|v| v.parse::<i64>().ok())
}

/// Mirrors upstream `findContigs`: a build matches when both reference
/// contigs are present with their canonical lengths (bare or `chr`-prefixed).
fn detect_build(text: &str) -> Build {
    let mut contigs: Vec<(String, u64)> = Vec::new();
    for line in text.lines() {
        if !line.starts_with("##contig=") {
            if line.starts_with("#CHROM") || !line.starts_with('#') {
                break;
            }
            continue;
        }
        let id = extract_attr(line, "ID=");
        let len = extract_attr(line, "length=").and_then(|s| s.parse::<u64>().ok());
        if let (Some(id), Some(len)) = (id, len) {
            contigs.push((id, len));
        }
    }
    let has = |name: &str, len: u64| {
        let chr = format!("chr{name}");
        contigs
            .iter()
            .any(|(c, l)| *l == len && (c == name || c == &chr))
    };
    if has("1", 249250621) && has("2", 243199373) {
        Build::Hg19
    } else if has("1", 248956422) && has("2", 242193529) {
        Build::Hg38
    } else if has("RF01", 3302) && has("RF02", 2687) {
        Build::RotavirusRf
    } else {
        Build::Undefined
    }
}

fn extract_attr(line: &str, key: &str) -> Option<String> {
    let start = line.find(key)? + key.len();
    let rest = &line[start..];
    let end = rest.find([',', '>']).unwrap_or(rest.len());
    Some(rest[..end].to_owned())
}

fn read_vcf_text(path: &Path) -> io::Result<String> {
    if path == Path::new("-") {
        let tmp = stdin_tmp_path();
        let mut data = Vec::new();
        io::stdin().lock().read_to_end(&mut data)?;
        fs::write(&tmp, data)?;
        let result = read_vcf_text(&tmp);
        let _ = fs::remove_file(&tmp);
        return result;
    }

    let fmt = format::detect_path(path).map_err(|e| io::Error::other(e.to_string()))?;
    if fmt.exact == Exact::Bcf {
        return htslib_rs::variant_io_compat::view_bcf_as_vcf_text_from_path_with_limit(path, None);
    }

    let mut text = String::new();
    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let file = File::open(path)?;
        let mut dec = MultiGzDecoder::new(file);
        dec.read_to_string(&mut text)?;
    } else {
        text = fs::read_to_string(path)?;
    }
    normalize_vcf_text(&mut text);
    Ok(text)
}

fn stdin_tmp_path() -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_nanos())
        .unwrap_or(0);
    std::env::temp_dir().join(format!(
        ".bcftools-rs-vcf2table-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_simple_variant_block() {
        let vcf = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=100>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tD\n\
1\t3000000\t.\tC\tA\t59.2\tPASS\tAN=4;AC=2\tGT:GQ\t0/1:245\t0/1:245\n";
        let o = render(vcf);
        assert!(o.starts_with("<<< 1:3000000:C (n. 1)\n\n# Variant\n"));
        assert!(o.contains("| CHROM  | 1       |\n"));
        assert!(o.contains("| QUAL   | 59.2    |\n"));
        assert!(o.contains("# INFO\n+-----+-----+-------+\n| KEY | IDX | VALUE |\n"));
        assert!(o.contains("| AN  |     | 4     |\n"));
        assert!(o.contains("# GENOTYPE TYPES\n"));
        assert!(o.contains("| HET  | 2     | 100 |\n"));
        assert!(o.contains("| C      | HET   | 0/1 | 245 |\n"));
        assert!(o.trim_end().ends_with(">>> 1:3000000:C (n. 1)"));
    }

    #[test]
    fn multi_value_info_uses_one_based_idx() {
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tC\tA\t.\t.\tTR=1,2;TA=1;TG=1,2,3\n";
        let o = render(vcf);
        assert!(o.contains("| TR  | 1   | 1     |\n"));
        assert!(o.contains("| TR  | 2   | 2     |\n"));
        assert!(o.contains("| TA  |     | 1     |\n"));
        assert!(o.contains("| TG  | 3   | 3     |\n"));
    }

    #[test]
    fn short_format_sample_padded_with_dots() {
        // FORMAT has 5 keys; sample C provides 2 -> padded XR/XG/XA = "."
        let vcf = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tC\tD\n\
1\t30\tid3\tC\tA\t59.2\tPASS\tAN=4\tGT:GQ:XR:XG:XA\t0/1:245\t0/1:245:1,2:1,2,3:2\n";
        let o = render(vcf);
        assert!(o.contains("| C      | HET   | 0/1 | 245 | .   | .     | .  |\n"));
        assert!(o.contains("| D      | HET   | 0/1 | 245 | 1,2 | 1,2,3 | 2  |\n"));
    }

    #[test]
    fn build_detection_defaults_to_undefined() {
        // chr2 length differs from hg19 -> undefined -> no prefix.
        let vcf = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=249250621>\n\
##contig=<ID=2,length=249250621>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tC\tA\t.\t.\t.\n";
        let o = render(vcf);
        assert!(o.starts_with("<<< 1:10:C (n. 1)\n"));
    }

    #[test]
    fn hg19_build_prefix() {
        let vcf = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=249250621>\n\
##contig=<ID=2,length=243199373>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t10\t.\tC\tA\t.\t.\t.\n";
        let o = render(vcf);
        assert!(o.starts_with("<<< GRCh37 :  1:10:C (n. 1)\n"));
    }
}
