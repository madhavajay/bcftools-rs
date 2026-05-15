//! `bcftools +counts` (upstream `bcftools/plugins/counts.c`).
//!
//! A minimal plugin that counts samples, SNPs, indels, MNPs, "others", and
//! total sites. Record classification routes through
//! `htslib_rs::variant::classify_variant`, the HTSlib `bcf_set_variant_type`
//! port, and is OR-combined across every ALT allele exactly like upstream's
//! `bcf_get_variant_types`.

use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::MultiGzDecoder;
use htslib_rs::format::{self, Compression, Exact};
use htslib_rs::variant::{VariantType, classify_variant};

use crate::vcf_compat::normalize_vcf_text;

#[derive(Debug, Default, PartialEq, Eq)]
pub struct Counts {
    pub samples: u64,
    pub snps: u64,
    pub indels: u64,
    pub mnps: u64,
    pub others: u64,
    pub sites: u64,
}

impl Counts {
    pub fn report(&self) -> String {
        format!(
            "Number of samples: {}\n\
Number of SNPs:    {}\n\
Number of INDELs:  {}\n\
Number of MNPs:    {}\n\
Number of others:  {}\n\
Number of sites:   {}\n",
            self.samples, self.snps, self.indels, self.mnps, self.others, self.sites
        )
    }
}

/// Reads the input VCF/BCF and returns the upstream-shaped count report.
pub fn run(input: &Path) -> io::Result<String> {
    let text = read_vcf_text(input)?;
    Ok(count_text(&text).report())
}

fn count_text(text: &str) -> Counts {
    let mut counts = Counts::default();
    for line in text.lines() {
        if line.starts_with("#CHROM") {
            counts.samples = sample_count(line);
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let (Some(_chrom), Some(_pos), Some(_id), Some(reference), Some(alt)) = (
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
            fields.next(),
        ) else {
            continue;
        };
        let mut record_type = VariantType::REF;
        for allele in alt.split(',') {
            if allele == "." {
                continue;
            }
            record_type |= classify_variant(reference, allele).variant_type;
        }
        if record_type.contains(VariantType::SNP) {
            counts.snps += 1;
        }
        if record_type.contains(VariantType::INDEL) {
            counts.indels += 1;
        }
        if record_type.contains(VariantType::MNP) {
            counts.mnps += 1;
        }
        if record_type.contains(VariantType::OTHER) {
            counts.others += 1;
        }
        counts.sites += 1;
    }
    counts
}

fn sample_count(chrom_header: &str) -> u64 {
    let cols: Vec<&str> = chrom_header.split('\t').collect();
    if cols.len() > 9 {
        (cols.len() - 9) as u64
    } else {
        0
    }
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
        ".bcftools-rs-counts-{}-{nanos}.tmp",
        std::process::id()
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const VCF: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1/1\n\
1\t20\t.\tA\tAT\t.\t.\t.\tGT\t0/1\t0/0\n\
1\t30\t.\tAC\tGT\t.\t.\t.\tGT\t0/1\t0/0\n\
1\t40\t.\tA\t<DEL>\t.\t.\t.\tGT\t0/1\t0/0\n\
1\t50\t.\tA\tG,AT\t.\t.\t.\tGT\t1/2\t0/0\n";

    #[test]
    fn counts_match_expected_classification() {
        let c = count_text(VCF);
        assert_eq!(c.samples, 2);
        // SNP rows: pos 10 (A>C) and pos 50 (multiallelic A>G,AT — SNP+INDEL).
        assert_eq!(c.snps, 2);
        // INDEL rows: pos 20 (A>AT), pos 50 (A>G,AT).
        assert_eq!(c.indels, 2);
        // MNP row: pos 30 (AC>GT).
        assert_eq!(c.mnps, 1);
        // OTHER row: pos 40 (<DEL> symbolic).
        assert_eq!(c.others, 1);
        assert_eq!(c.sites, 5);
    }

    #[test]
    fn report_matches_upstream_layout() {
        let c = Counts {
            samples: 3,
            snps: 7,
            indels: 2,
            mnps: 0,
            others: 1,
            sites: 10,
        };
        assert_eq!(
            c.report(),
            "Number of samples: 3\n\
Number of SNPs:    7\n\
Number of INDELs:  2\n\
Number of MNPs:    0\n\
Number of others:  1\n\
Number of sites:   10\n"
        );
    }

    #[test]
    fn no_samples_sites_only() {
        let sites_only = "##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t1\t.\tA\tT\t.\t.\t.\n";
        let c = count_text(sites_only);
        assert_eq!(c.samples, 0);
        assert_eq!(c.snps, 1);
        assert_eq!(c.sites, 1);
    }
}
