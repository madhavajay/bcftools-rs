//! In-process plugin record-processing implementations.
//!
//! Each module mirrors one upstream `bcftools/plugins/<name>.c` and is invoked
//! through `bcftools plugin <name>` / `bcftools +<name>`. The static
//! registry/listing surface lives in [`super::plugin`]; this module holds the
//! actual per-plugin algorithms as they are ported.

pub mod ad_bias;
pub mod add_variantkey;
pub mod af_dist;
pub mod allele_length;
pub mod check_ploidy;
pub mod check_sparsity;
pub mod contrast;
pub mod counts;
pub mod dosage;
pub mod fill_an_ac;
pub mod fill_from_fasta;
pub mod fixploidy;
pub mod fixref;
pub mod frameshifts;
pub mod gtisec;
pub mod gtsubset;
pub mod guess_ploidy;
pub mod gvcfz;
pub mod impute_info;
pub mod indel_stats;
pub mod isecgt;
pub mod mendelian2;
pub mod missing2ref;
pub mod parental_origin;
pub mod prune;
pub mod remove_overlaps;
pub mod scatter;
pub mod setgt;
pub mod smpl_stats;
pub mod split;
pub mod split_vep;
pub mod tag2tag;
pub mod trio_stats;
pub mod trio_switch_rate;
pub mod variant_distance;
pub mod variantkey;
pub mod variantkey_hex;
pub mod vcf2table;
