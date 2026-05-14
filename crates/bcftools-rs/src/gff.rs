//! GFF/GTF parsing helpers for the future `csq` port.
//!
//! Upstream `gff.c` carries a large transcript model on top of line parsing.
//! This module starts with the format-stable foundation: parse records into
//! normalized half-open coordinates and decode GFF3/GTF attributes.

use std::collections::BTreeMap;
use std::io;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
    NotStranded,
}

impl Strand {
    fn parse(raw: &str) -> io::Result<Self> {
        match raw {
            "+" => Ok(Self::Forward),
            "-" => Ok(Self::Reverse),
            "." => Ok(Self::NotStranded),
            "?" => Ok(Self::Unknown),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid GFF strand {raw:?}"),
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Phase {
    Zero,
    One,
    Two,
    None,
}

impl Phase {
    fn parse(raw: &str) -> io::Result<Self> {
        match raw {
            "0" => Ok(Self::Zero),
            "1" => Ok(Self::One),
            "2" => Ok(Self::Two),
            "." => Ok(Self::None),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid GFF phase {raw:?}"),
            )),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GffRecord {
    pub seqid: String,
    pub source: Option<String>,
    pub feature_type: String,
    /// 0-based inclusive start coordinate.
    pub start: i64,
    /// 0-based exclusive end coordinate.
    pub end: i64,
    pub score: Option<String>,
    pub strand: Strand,
    pub phase: Phase,
    pub attributes: BTreeMap<String, Vec<String>>,
}

impl GffRecord {
    pub fn id(&self) -> Option<&str> {
        self.attributes
            .get("ID")
            .and_then(|values| values.first())
            .map(String::as_str)
    }

    pub fn parent_ids(&self) -> impl Iterator<Item = &str> {
        self.attributes
            .get("Parent")
            .into_iter()
            .flat_map(|values| values.iter().map(String::as_str))
    }

    pub fn attr_first(&self, key: &str) -> Option<&str> {
        self.attributes
            .get(key)
            .and_then(|values| values.first())
            .map(String::as_str)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct GffModel {
    pub genes: BTreeMap<String, Gene>,
    pub transcripts: BTreeMap<String, Transcript>,
    pub unplaced: Vec<GffRecord>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub id: String,
    pub record: GffRecord,
    pub transcripts: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Transcript {
    pub id: String,
    pub gene_id: Option<String>,
    pub record: Option<GffRecord>,
    pub exons: Vec<GffRecord>,
    pub cds: Vec<GffRecord>,
    pub other_features: Vec<GffRecord>,
}

impl Transcript {
    fn new(id: String, gene_id: Option<String>) -> Self {
        Self {
            id,
            gene_id,
            record: None,
            exons: Vec::new(),
            cds: Vec::new(),
            other_features: Vec::new(),
        }
    }
}

pub fn build_model(records: impl IntoIterator<Item = GffRecord>) -> GffModel {
    let mut model = GffModel::default();

    for record in records {
        if is_gene_feature(&record.feature_type) {
            if let Some(id) = gene_id(&record) {
                model.genes.entry(id.clone()).or_insert(Gene {
                    id,
                    record,
                    transcripts: Vec::new(),
                });
            } else {
                model.unplaced.push(record);
            }
            continue;
        }

        if is_transcript_feature(&record.feature_type) {
            if let Some(id) = transcript_id(&record) {
                let gene_id = transcript_gene_id(&record);
                let transcript = model
                    .transcripts
                    .entry(id.clone())
                    .or_insert_with(|| Transcript::new(id.clone(), gene_id.clone()));
                if transcript.gene_id.is_none() {
                    transcript.gene_id = gene_id.clone();
                }
                transcript.record = Some(record);
                if let Some(gene_id) = gene_id {
                    link_gene_transcript(&mut model, &gene_id, &id);
                }
            } else {
                model.unplaced.push(record);
            }
            continue;
        }

        let parents = feature_parent_transcripts(&record);
        if parents.is_empty() {
            model.unplaced.push(record);
            continue;
        }

        for parent in parents {
            let transcript = model
                .transcripts
                .entry(parent.clone())
                .or_insert_with(|| Transcript::new(parent.clone(), None));
            match record.feature_type.as_str() {
                "exon" => transcript.exons.push(record.clone()),
                "CDS" => transcript.cds.push(record.clone()),
                _ => transcript.other_features.push(record.clone()),
            }
        }
    }

    let gene_links = model
        .transcripts
        .values()
        .filter_map(|transcript| {
            transcript
                .gene_id
                .as_ref()
                .map(|gene_id| (gene_id.clone(), transcript.id.clone()))
        })
        .collect::<Vec<_>>();
    for (gene_id, transcript_id) in gene_links {
        link_gene_transcript(&mut model, &gene_id, &transcript_id);
    }

    for transcript in model.transcripts.values_mut() {
        sort_transcript_features(transcript);
    }

    model
}

fn is_gene_feature(feature_type: &str) -> bool {
    feature_type == "gene"
}

fn is_transcript_feature(feature_type: &str) -> bool {
    matches!(
        feature_type,
        "mRNA" | "transcript" | "primary_transcript" | "lnc_RNA" | "ncRNA" | "rRNA" | "tRNA"
    ) || feature_type.ends_with("RNA")
}

fn gene_id(record: &GffRecord) -> Option<String> {
    record
        .id()
        .or_else(|| record.attr_first("gene_id"))
        .map(str::to_string)
}

fn transcript_id(record: &GffRecord) -> Option<String> {
    record
        .id()
        .or_else(|| record.attr_first("transcript_id"))
        .map(str::to_string)
}

fn transcript_gene_id(record: &GffRecord) -> Option<String> {
    record
        .parent_ids()
        .next()
        .or_else(|| record.attr_first("gene_id"))
        .map(str::to_string)
}

fn feature_parent_transcripts(record: &GffRecord) -> Vec<String> {
    let parents = record
        .parent_ids()
        .map(str::to_string)
        .collect::<Vec<String>>();
    if !parents.is_empty() {
        return parents;
    }

    record
        .attr_first("transcript_id")
        .map(|id| vec![id.to_string()])
        .unwrap_or_default()
}

fn link_gene_transcript(model: &mut GffModel, gene_id: &str, transcript_id: &str) {
    if let Some(gene) = model.genes.get_mut(gene_id)
        && !gene.transcripts.iter().any(|id| id == transcript_id)
    {
        gene.transcripts.push(transcript_id.to_string());
    }
}

fn sort_transcript_features(transcript: &mut Transcript) {
    let reverse = transcript
        .record
        .as_ref()
        .is_some_and(|record| record.strand == Strand::Reverse);
    sort_features(&mut transcript.exons, reverse);
    sort_features(&mut transcript.cds, reverse);
    sort_features(&mut transcript.other_features, reverse);
}

fn sort_features(features: &mut [GffRecord], reverse: bool) {
    features.sort_by(|a, b| {
        let ordering = a
            .seqid
            .cmp(&b.seqid)
            .then_with(|| a.start.cmp(&b.start))
            .then_with(|| a.end.cmp(&b.end))
            .then_with(|| a.feature_type.cmp(&b.feature_type));
        if reverse {
            ordering.reverse()
        } else {
            ordering
        }
    });
}

pub fn parse_line(line: &str) -> io::Result<Option<GffRecord>> {
    let line = line.trim_end_matches(['\r', '\n']);

    if line.is_empty() || line.starts_with('#') {
        return Ok(None);
    }

    let mut fields = line.split('\t');
    let seqid = required_field(&mut fields, "seqid")?;
    let source = required_field(&mut fields, "source")?;
    let feature_type = required_field(&mut fields, "type")?;
    let start = parse_one_based_start(required_field(&mut fields, "start")?)?;
    let end = parse_one_based_end(required_field(&mut fields, "end")?, start)?;
    let score = none_if_dot(required_field(&mut fields, "score")?);
    let strand = Strand::parse(required_field(&mut fields, "strand")?)?;
    let phase = Phase::parse(required_field(&mut fields, "phase")?)?;
    let attributes = parse_attributes(required_field(&mut fields, "attributes")?);

    if fields.next().is_some() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "too many GFF fields",
        ));
    }

    Ok(Some(GffRecord {
        seqid: percent_decode(seqid),
        source: none_if_dot(source).map(percent_decode),
        feature_type: percent_decode(feature_type),
        start,
        end,
        score: score.map(ToOwned::to_owned),
        strand,
        phase,
        attributes,
    }))
}

pub fn parse_attributes(raw: &str) -> BTreeMap<String, Vec<String>> {
    let mut attributes = BTreeMap::new();

    if raw == "." {
        return attributes;
    }

    for field in raw
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
    {
        let (key, value) = field
            .split_once('=')
            .or_else(|| field.split_once(char::is_whitespace))
            .map(|(key, value)| (key.trim(), value.trim()))
            .unwrap_or((field, ""));

        if key.is_empty() {
            continue;
        }

        let values = value
            .trim_matches('"')
            .split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(percent_decode)
            .collect::<Vec<_>>();

        attributes.entry(percent_decode(key)).or_insert(values);
    }

    attributes
}

fn required_field<'a>(
    fields: &mut impl Iterator<Item = &'a str>,
    name: &str,
) -> io::Result<&'a str> {
    fields.next().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("missing GFF {name} field"),
        )
    })
}

fn parse_one_based_start(raw: &str) -> io::Result<i64> {
    let start = raw
        .parse::<i64>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    if start < 1 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GFF start must be 1-based and positive",
        ));
    }

    Ok(start - 1)
}

fn parse_one_based_end(raw: &str, start: i64) -> io::Result<i64> {
    let end = raw
        .parse::<i64>()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    if end <= start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GFF end must be greater than or equal to start",
        ));
    }

    Ok(end)
}

fn none_if_dot(raw: &str) -> Option<&str> {
    (raw != ".").then_some(raw)
}

fn percent_decode(raw: &str) -> String {
    let bytes = raw.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;

    while i < bytes.len() {
        if bytes[i] == b'%'
            && i + 2 < bytes.len()
            && let Some(value) = decode_hex_pair(bytes[i + 1], bytes[i + 2])
        {
            out.push(value);
            i += 3;
            continue;
        }

        out.push(bytes[i]);
        i += 1;
    }

    String::from_utf8_lossy(&out).into_owned()
}

fn decode_hex_pair(high: u8, low: u8) -> Option<u8> {
    Some(hex_value(high)? << 4 | hex_value(low)?)
}

fn hex_value(byte: u8) -> Option<u8> {
    match byte {
        b'0'..=b'9' => Some(byte - b'0'),
        b'a'..=b'f' => Some(byte - b'a' + 10),
        b'A'..=b'F' => Some(byte - b'A' + 10),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_line_skips_comments_and_blank_lines() {
        assert_eq!(parse_line("#gff-version 3").unwrap(), None);
        assert_eq!(parse_line("").unwrap(), None);
    }

    #[test]
    fn parse_gff3_record_normalizes_coordinates_and_attributes() {
        let record = parse_line(
            "chr1\tsrc\tmRNA\t10\t20\t.\t+\t.\tID=tx1;Parent=gene1,gene2;Name=hello%20world",
        )
        .unwrap()
        .unwrap();

        assert_eq!(record.seqid, "chr1");
        assert_eq!(record.source.as_deref(), Some("src"));
        assert_eq!(record.feature_type, "mRNA");
        assert_eq!(record.start, 9);
        assert_eq!(record.end, 20);
        assert_eq!(record.score, None);
        assert_eq!(record.strand, Strand::Forward);
        assert_eq!(record.phase, Phase::None);
        assert_eq!(record.id(), Some("tx1"));
        assert_eq!(record.parent_ids().collect::<Vec<_>>(), ["gene1", "gene2"]);
        assert_eq!(record.attributes["Name"], ["hello world"]);
    }

    #[test]
    fn parse_gtf_style_attributes() {
        let record =
            parse_line("chr1\t.\texon\t1\t3\t42\t-\t0\tgene_id \"g1\"; transcript_id \"t1\";")
                .unwrap()
                .unwrap();

        assert_eq!(record.source, None);
        assert_eq!(record.start, 0);
        assert_eq!(record.end, 3);
        assert_eq!(record.score.as_deref(), Some("42"));
        assert_eq!(record.strand, Strand::Reverse);
        assert_eq!(record.phase, Phase::Zero);
        assert_eq!(record.attributes["gene_id"], ["g1"]);
        assert_eq!(record.attributes["transcript_id"], ["t1"]);
    }

    #[test]
    fn parse_rejects_malformed_coordinates() {
        assert!(parse_line("chr1\t.\tgene\t0\t1\t.\t.\t.\tID=g").is_err());
        assert!(parse_line("chr1\t.\tgene\t5\t4\t.\t.\t.\tID=g").is_err());
    }

    #[test]
    fn parse_rejects_bad_strand_and_phase() {
        assert!(parse_line("chr1\t.\tgene\t1\t2\t.\tx\t.\tID=g").is_err());
        assert!(parse_line("chr1\t.\tCDS\t1\t2\t.\t+\t3\tID=cds").is_err());
    }

    #[test]
    fn build_model_groups_gff3_genes_transcripts_and_features() {
        let records = parse_records(
            "\
chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1
chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1
chr1\t.\texon\t40\t60\t.\t+\t.\tID=ex2;Parent=tx1
chr1\t.\texon\t10\t20\t.\t+\t.\tID=ex1;Parent=tx1
chr1\t.\tCDS\t50\t60\t.\t+\t2\tID=cds2;Parent=tx1
chr1\t.\tCDS\t10\t18\t.\t+\t0\tID=cds1;Parent=tx1
",
        );

        let model = build_model(records);
        let gene = &model.genes["gene1"];
        assert_eq!(gene.transcripts, ["tx1"]);

        let transcript = &model.transcripts["tx1"];
        assert_eq!(transcript.gene_id.as_deref(), Some("gene1"));
        assert_eq!(
            transcript.exons.iter().map(|r| r.start).collect::<Vec<_>>(),
            [9, 39]
        );
        assert_eq!(
            transcript.cds.iter().map(|r| r.start).collect::<Vec<_>>(),
            [9, 49]
        );
        assert!(model.unplaced.is_empty());
    }

    #[test]
    fn build_model_supports_gtf_ids_and_reverse_ordering() {
        let records = parse_records(
            "\
chr1\t.\texon\t1\t10\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";
chr1\t.\tCDS\t3\t8\t.\t-\t0\tgene_id \"g1\"; transcript_id \"t1\";
chr1\t.\ttranscript\t1\t50\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";
chr1\t.\texon\t30\t50\t.\t-\t.\tgene_id \"g1\"; transcript_id \"t1\";
",
        );

        let model = build_model(records);
        let transcript = &model.transcripts["t1"];
        assert_eq!(transcript.gene_id.as_deref(), Some("g1"));
        assert_eq!(
            transcript.exons.iter().map(|r| r.start).collect::<Vec<_>>(),
            [29, 0]
        );
        assert_eq!(
            transcript.cds.iter().map(|r| r.start).collect::<Vec<_>>(),
            [2]
        );
    }

    #[test]
    fn build_model_keeps_unparented_features_unplaced() {
        let records = parse_records("chr1\t.\texon\t1\t10\t.\t+\t.\tID=orphan\n");
        let model = build_model(records);

        assert_eq!(model.unplaced.len(), 1);
        assert!(model.transcripts.is_empty());
    }

    fn parse_records(raw: &str) -> Vec<GffRecord> {
        raw.lines()
            .map(|line| parse_line(line).unwrap().unwrap())
            .collect()
    }
}
