//! bcftools-shaped synced reader helpers.
//!
//! This module is a facade over the currently exposed `htslib-rs` synced
//! pairing primitives. It is intentionally small: it gives command ports a
//! stable place to parse `--collapse` modes, collect per-position variant
//! groups across inputs, and ask for paired allele rows.

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{self, BufReader};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use htslib_rs::core::Position;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant_io_compat::{SyncedPairLogic, SyncedVariantGroup};
use htslib_rs::vcf;

pub type SyncedPairRows = Vec<Vec<Option<String>>>;
pub type SyncedPositionRows = Vec<(SyncedSiteKey, SyncedPairRows)>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CollapseMode {
    None,
    Snps,
    Indels,
    Both,
    Any,
    Some,
    Id,
    All,
}

impl CollapseMode {
    pub fn pair_logic(self) -> SyncedPairLogic {
        match self {
            Self::None | Self::Id => SyncedPairLogic::Exact,
            Self::Snps => SyncedPairLogic::Snps,
            Self::Indels => SyncedPairLogic::Indels,
            Self::Both => SyncedPairLogic::Both,
            Self::Any | Self::All => SyncedPairLogic::All,
            Self::Some => SyncedPairLogic::Some,
        }
    }
}

impl FromStr for CollapseMode {
    type Err = ParseCollapseModeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "none" => Ok(Self::None),
            "snps" => Ok(Self::Snps),
            "indels" => Ok(Self::Indels),
            "both" => Ok(Self::Both),
            "any" => Ok(Self::Any),
            "some" => Ok(Self::Some),
            "id" => Ok(Self::Id),
            "all" => Ok(Self::All),
            _ => Err(ParseCollapseModeError),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ParseCollapseModeError;

impl std::fmt::Display for ParseCollapseModeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("unrecognized collapse mode")
    }
}

impl std::error::Error for ParseCollapseModeError {}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SyncedSiteKey {
    pub reference_sequence_order: usize,
    pub reference_sequence_name: String,
    pub position: Position,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SyncedSite {
    pub key: SyncedSiteKey,
    pub groups: Vec<SyncedVariantGroup>,
}

impl SyncedSite {
    pub fn paired_rows(&self, collapse: CollapseMode) -> Vec<Vec<Option<String>>> {
        htslib_rs::variant_io_compat::pair_synced_variant_groups(
            &self.groups,
            collapse.pair_logic(),
        )
    }
}

#[derive(Debug, Clone)]
pub struct SyncedReader {
    inputs: Vec<PathBuf>,
    collapse: CollapseMode,
    regions: Vec<String>,
    targets: Vec<String>,
}

impl Default for SyncedReader {
    fn default() -> Self {
        Self {
            inputs: Vec::new(),
            collapse: CollapseMode::None,
            regions: Vec::new(),
            targets: Vec::new(),
        }
    }
}

impl SyncedReader {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_path<P>(&mut self, path: P)
    where
        P: Into<PathBuf>,
    {
        self.inputs.push(path.into());
    }

    pub fn set_collapse(&mut self, collapse: CollapseMode) {
        self.collapse = collapse;
    }

    pub fn set_regions<I, S>(&mut self, regions: I)
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.regions = regions.into_iter().map(Into::into).collect();
    }

    pub fn set_targets<I, S>(&mut self, targets: I)
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.targets = targets.into_iter().map(Into::into).collect();
    }

    pub fn regions(&self) -> &[String] {
        &self.regions
    }

    pub fn targets(&self) -> &[String] {
        &self.targets
    }

    pub fn read_sites(&self) -> io::Result<Vec<SyncedSite>> {
        let sites = collect_synced_sites_from_paths(&self.inputs)?;
        if self.regions.is_empty() && self.targets.is_empty() {
            return Ok(sites);
        }

        Ok(sites
            .into_iter()
            .filter(|site| site_matches_text_filters(site, &self.regions, &self.targets))
            .collect())
    }

    pub fn paired_rows(&self) -> io::Result<SyncedPositionRows> {
        self.read_sites().map(|sites| {
            sites
                .into_iter()
                .map(|site| {
                    let rows = site.paired_rows(self.collapse);
                    (site.key, rows)
                })
                .collect()
        })
    }
}

pub fn collect_synced_sites_from_paths<P>(paths: &[P]) -> io::Result<Vec<SyncedSite>>
where
    P: AsRef<Path>,
{
    let mut sites: BTreeMap<SyncedSiteKey, HashMap<String, Vec<usize>>> = BTreeMap::new();

    for (input_index, path) in paths.iter().enumerate() {
        let fmt =
            format::detect_path(path.as_ref()).map_err(|e| io::Error::other(e.to_string()))?;
        let (contigs, records) = read_variant_records(path.as_ref(), fmt)?;
        for record in records {
            let Some(position) = record.variant_start() else {
                continue;
            };
            let reference_sequence_name = record.reference_sequence_name().to_string();
            let reference_sequence_order = contigs
                .get(&reference_sequence_name)
                .copied()
                .unwrap_or(usize::MAX);
            let variant = variant_summary(&record);
            let key = SyncedSiteKey {
                reference_sequence_order,
                reference_sequence_name,
                position,
            };
            sites
                .entry(key)
                .or_default()
                .entry(variant)
                .or_default()
                .push(input_index);
        }
    }

    Ok(sites
        .into_iter()
        .map(|(key, variants)| {
            let mut groups = variants
                .into_iter()
                .map(|(variant, input_indexes)| SyncedVariantGroup {
                    variants: vec![variant],
                    input_indexes,
                })
                .collect::<Vec<_>>();
            groups.sort_by(|a, b| a.variants.cmp(&b.variants));
            SyncedSite { key, groups }
        })
        .collect())
}

fn read_variant_records(
    path: &Path,
    fmt: format::Format,
) -> io::Result<(HashMap<String, usize>, Vec<vcf::variant::RecordBuf>)> {
    use htslib_rs::bcf;

    if fmt.exact == Exact::Bcf {
        let mut reader = File::open(path).map(bcf::io::Reader::new)?;
        let header = reader.read_header()?;
        let contigs = contig_order(&header);
        let records = reader
            .record_bufs(&header)
            .collect::<io::Result<Vec<_>>>()?;
        return Ok((contigs, records));
    }

    if fmt.compression == Compression::Bgzf || fmt.compression == Compression::Gzip {
        let f = File::open(path)?;
        let dec = flate2::read::MultiGzDecoder::new(f);
        let mut reader = vcf::io::Reader::new(BufReader::new(dec));
        let header = reader.read_header()?;
        let contigs = contig_order(&header);
        let records = reader
            .records()
            .map(|result| {
                let record = result?;
                vcf::variant::RecordBuf::try_from_variant_record(&header, &record)
            })
            .collect::<io::Result<Vec<_>>>()?;
        return Ok((contigs, records));
    }

    let mut reader = File::open(path)
        .map(BufReader::new)
        .map(vcf::io::Reader::new)?;
    let header = reader.read_header()?;
    let contigs = contig_order(&header);
    let records = reader
        .records()
        .map(|result| {
            let record = result?;
            vcf::variant::RecordBuf::try_from_variant_record(&header, &record)
        })
        .collect::<io::Result<Vec<_>>>()?;
    Ok((contigs, records))
}

fn contig_order(header: &vcf::Header) -> HashMap<String, usize> {
    header
        .contigs()
        .iter()
        .enumerate()
        .map(|(index, (name, _))| (name.to_string(), index))
        .collect()
}

fn variant_summary(record: &vcf::variant::RecordBuf) -> String {
    let reference = record.reference_bases();
    let alternates = record.alternate_bases().as_ref();
    if alternates.is_empty() {
        return format!("{reference}>.");
    }
    alternates
        .iter()
        .map(|alternate| format!("{reference}>{alternate}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn site_matches_text_filters(site: &SyncedSite, regions: &[String], targets: &[String]) -> bool {
    (regions.is_empty()
        || regions
            .iter()
            .any(|region| site_matches_region(site, region)))
        && (targets.is_empty()
            || targets
                .iter()
                .any(|target| site_matches_region(site, target)))
}

fn site_matches_region(site: &SyncedSite, raw: &str) -> bool {
    let (name, range) = raw.split_once(':').unwrap_or((raw, ""));
    if name != site.key.reference_sequence_name {
        return false;
    }
    if range.is_empty() {
        return true;
    }

    let (start, end) = range.split_once('-').unwrap_or((range, ""));
    let position = usize::from(site.key.position);
    let start = start.replace(',', "").parse::<usize>().unwrap_or(1);
    let end = if end.is_empty() {
        start
    } else {
        end.replace(',', "").parse::<usize>().unwrap_or(usize::MAX)
    };
    start <= position && position <= end
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "##fileformat=VCFv4.2\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    #[test]
    fn collapse_modes_parse_to_pairing_logic() {
        assert_eq!("none".parse::<CollapseMode>().unwrap(), CollapseMode::None);
        assert_eq!("snps".parse::<CollapseMode>().unwrap(), CollapseMode::Snps);
        assert!("bad".parse::<CollapseMode>().is_err());
        assert_eq!(CollapseMode::Any.pair_logic(), SyncedPairLogic::All);
    }

    #[test]
    fn synced_reader_pairs_records_across_inputs() {
        let dir = tempfile::TempDir::new().unwrap();
        let a = dir.path().join("a.vcf");
        let b = dir.path().join("b.vcf");
        std::fs::write(
            &a,
            format!("{HEADER}1\t10\t.\tA\tC\t.\tPASS\t.\n1\t20\t.\tG\tT\t.\tPASS\t.\n"),
        )
        .unwrap();
        std::fs::write(
            &b,
            format!("{HEADER}1\t10\t.\tA\tG\t.\tPASS\t.\n1\t30\t.\tC\tCA\t.\tPASS\t.\n"),
        )
        .unwrap();

        let mut reader = SyncedReader::new();
        reader.add_path(a);
        reader.add_path(b);
        reader.set_collapse(CollapseMode::Snps);
        let rows = reader.paired_rows().unwrap();
        assert_eq!(rows.len(), 3);

        let pos10 = &rows[0];
        assert_eq!(usize::from(pos10.0.position), 10);
        assert_eq!(pos10.1, vec![vec![Some("C".into()), Some("G".into())]]);
    }

    #[test]
    fn synced_reader_filters_regions_and_targets() {
        let dir = tempfile::TempDir::new().unwrap();
        let a = dir.path().join("a.vcf");
        std::fs::write(
            &a,
            format!("{HEADER}1\t10\t.\tA\tC\t.\tPASS\t.\n1\t20\t.\tG\tT\t.\tPASS\t.\n"),
        )
        .unwrap();

        let mut reader = SyncedReader::new();
        reader.add_path(a);
        reader.set_regions(["1:15-25"]);
        let sites = reader.read_sites().unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(usize::from(sites[0].key.position), 20);
    }
}
