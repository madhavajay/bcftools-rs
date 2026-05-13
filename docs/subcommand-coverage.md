# Subcommand and Plugin Coverage Map

This is the Phase 2 surface map for the bcftools-rs port. It maps each
bcftools command/plugin group to upstream source files, the main HTSlib APIs it
uses, current `htslib-rs` coverage, and the bcftools-rs internal module that
owns the bcftools-specific dependency.

Coverage status values:

- `exposed`: `htslib-rs` already has a Rust API or compatibility helper.
- `partial`: `htslib-rs` has a related helper but not the full bcftools call
  surface.
- `needed`: missing or not broad enough for the bcftools call sites.
- `out-of-scope`: intentionally not part of the Rust-only target.

Internal owner values refer to Phase 1 modules in `TODO.md`.

## Shared HTSlib API Matrix

| HTSlib / bcftools API family | `htslib-rs` coverage | bcftools-rs owner | Downstream users |
| --- | --- | --- | --- |
| VCF/BCF read/write (`bcf_hdr_*`, `bcf_read`, `bcf_write`, VCF text I/O) | partial | command-specific I/O plus `io.rs` | `view`, `query`, `filter`, `annotate`, `merge`, `norm`, `call`, `mpileup`, plugins |
| VCF/BCF record getters (`bcf_get_info_*`, `bcf_get_format_*`, genotypes) | partial | filter engine, convert formatter, command ports | `view`, `query`, `stats`, `call`, `annotate`, plugins |
| VCF/BCF mutation (`bcf_update_*`, alleles, FILTER/INFO/FORMAT/ID/QUAL/POS) | partial | command ports plus htslib-rs gap `bcf_update_*` | `annotate`, `norm`, `filter`, `call`, many plugins |
| Header translation (`bcf_translate`, merged/per-input header dictionaries) | partial | synced reader wrapper plus htslib-rs gap `bcf_translate` | `merge`, `concat`, `isec`, plugins |
| Synced reader (`bcf_sr_*`, regions, targets, collapse modes, paired readers) | partial | `synced.rs` plus htslib-rs gap `synced_bcf_reader` | `view`, `query`, `isec`, `merge`, `stats`, `annotate` |
| Filter expressions (`filter_init`, `filter_test`, sample-vector semantics) | needed | `filter/` | `view`, `filter`, `query`, `annotate`, `norm`, `stats`, `call`, `mpileup`, plugins |
| Region indexes (`regidx_*`) | exposed | `regidx.rs` | `view`, `filter`, `annotate`, `isec`, plugins |
| Sample indexes (`smpl_ilist_*`) | exposed | `smpl_ilist.rs` | `view`, `call`, `stats`, plugins |
| FASTA/FAI (`faidx_*`) | exposed | `reference.rs` | `csq`, `consensus`, `mpileup`, `norm`, `+fill-from-fasta` |
| Tabix/TBI/CSI and generic indexes (`tbx_*`, `bcf_index_*`, `hts_idx_*`) | exposed for current ports, partial for all writer-thread cases | `io.rs`, command ports | `index`, `tabix`, `view`, `sort`, writers with `-W` |
| BGZF threading (`hts_set_threads`, writer worker pools) | partial | htslib-rs gap `hts_set_threads` | `view`, `merge`, `norm`, `concat`, `sort` |
| Pileup (`bam_mplp_*`, BAQ/probaln, multi-input pileup) | partial | htslib-rs gaps `Pileup iterator surface`, `BAQ/probaln_glocal` | `mpileup`, `call` |
| Convert formatter (`convert_init`, `convert_line`, query format tokens) | needed | `convert/` | `query`, `convert`, plugins |
| Allele/window buffers (`abuf`, `vcfbuf`) | needed | `abuf.rs`, `vcfbuf.rs` | `norm`, `merge`, `+remove-overlaps`, `+prune`, `+scatter` |
| gVCF helpers | needed | `gvcf.rs` | `call`, `convert`, `+gvcfz` |
| GFF parser | needed | `gff.rs` | `csq` |
| Ploidy spec | needed | `ploidy.rs` | `call`, `+fixploidy`, `+guess-ploidy` |
| HMM/numerics | needed | `hmm.rs`, `numerics.rs` | `roh`, `cnv`, `gtcheck`, `polysomy`, stats plugins |

## Command Surface

| Command | Upstream source files | Main APIs called | `htslib-rs` status | Internal owner / blocker |
| --- | --- | --- | --- | --- |
| `index` | `vcfindex.c`, `version.c` | `bcf_index_build3`, `tbx_index_build`, `hts_idx_load`, stats over CSI/TBI | exposed for current VCF/BCF/TBI/CSI paths | complete command; `io.rs`, `index_compat` |
| `tabix` | `tabix.c` | `tbx_index_build`, `tbx_index_load`, `tbx_itr_querys`, BGZF getline | exposed for BED/GFF/SAM/VCF presets | complete command; `tabix_compat` |
| `head` | `vcfhead.c` | VCF/BCF header read, raw VCF header text | exposed | complete command |
| `view` | `vcfview.c` | VCF/BCF I/O, `bcf_sr_*`, `filter_*`, sample subsets, regions/targets, `bcf_translate`, `bcf_update_*` | partial | blocked on `filter/`, `synced.rs`, sample/region semantics |
| `query` | `vcfquery.c`, `convert.c` | `convert_*`, `filter_*`, `bcf_sr_*`, sample lists | partial | blocked on `convert/`, `filter/`, `synced.rs` |
| `stats` | `vcfstats.c`, `smpl_ilist.c`, `filter.c` | typed INFO/FORMAT/GT getters, sample subsets, filters, synced pairs | partial | blocked on `filter/`, `synced.rs`, numerics |
| `isec` | `vcfisec.c`, `vcmp.c` | `bcf_sr_*`, collapse modes, header translation, output writers | partial | blocked on `synced.rs`, `bcf_translate` gap |
| `sort` | `vcfsort.c`, `extsort.c`, `version.c` | VCF/BCF read/write, external sort temp files, optional `init_index2` | partial | in-memory subset complete; full item blocked on `extsort` disk-backed path |
| `norm` | `vcfnorm.c`, `abuf.c`, `vcfbuf.c`, `vcmp.c`, `gvcf.c` | allele mutation, header/record updates, FASTA, buffers, filters | partial | blocked on `abuf.rs`, `vcfbuf.rs`, `bcf_update_*`, `filter/` |
| `concat` | `vcfconcat.c`, `vcmp.c` | header compatibility, naive BCF concat, synced regions, translation | partial | blocked on `bcf_translate`, `synced.rs`, writer indexing |
| `merge` | `vcfmerge.c`, `vcmp.c`, `abuf.c` | multi-reader sync, header merge/translate, INFO merge rules, record mutation | partial | blocked on `synced.rs`, `bcf_translate`, `bcf_update_*`, `abuf.rs` |
| `reheader` | `reheader.c` | raw header replacement, sample rename, FAI contig fill, BCF in-place temp prefix | partial | blocked on BCF header rewrite and `reference.rs` integration |
| `convert` | `vcfconvert.c`, `convert.c`, `gvcf.c`, `tsv2vcf.c` | format conversions, FORMAT/INFO token formatter, gVCF, TSV helpers | partial | blocked on `convert/`, `gvcf.rs`; TSV helper exists |
| `filter` | `vcffilter.c`, `filter.c` | expression engine, FILTER/GT mutation, masks, gap filters | partial | blocked on `filter/`, `bcf_update_*` |
| `annotate` | `vcfannotate.c`, `regidx.c`, `filter.c`, `tsv2vcf.c` | annotation readers, region indexes, column mapping, `bcf_update_*`, filters | partial | blocked on `bcf_update_*`, `synced.rs`, `filter/` |
| `mpileup` | `mpileup.c`, `bam2bcf*.c`, `bam_sample.c`, `read_consensus.c`, `cigar_state.h`, `mw.h` | BAM pileup, BAQ/probaln, FASTA, BCF writes | partial | blocked on pileup/BAQ gaps, `reference.rs`, numerics |
| `call` | `vcfcall.c`, `mcall.c`, `ccall.c`, `prob1.c`, `em.c`, `ploidy.c`, `gvcf.c` | genotype likelihoods, ploidy, constrained alleles, gVCF, FORMAT updates | partial | blocked on `ploidy.rs`, `gvcf.rs`, `numerics.rs`, `bcf_update_*` |
| `consensus` | `consensus.c`, `filter.c`, `read_consensus.c` | FASTA fetch, VCF record application, sample/haplotype selection, filters | partial | blocked on `reference.rs` wiring, filter/sample semantics |
| `csq` | `csq.c`, `gff.c` | GFF parse, FASTA, phased genotypes, record annotation | partial | blocked on `gff.rs`, `reference.rs`, `bcf_update_*` |
| `roh` | `vcfroh.c`, `HMM.c`, `peakfit.c` | HMM, FORMAT/GT/AF getters, filters | partial | blocked on `hmm.rs`, `numerics.rs`, `filter/` |
| `cnv` | `vcfcnv.c`, `HMM.c`, `peakfit.c` | HMM, BAF/LRR parsing, numerics | partial | blocked on `hmm.rs`, `numerics.rs` |
| `gtcheck` | `vcfgtcheck.c`, `prob1.c`, `filter.c` | genotype likelihoods, sample concordance, filters | partial | blocked on `numerics.rs`, `filter/` |
| `polysomy` | `polysomy.c`, `peakfit.c`, GSL upstream | chromosomal-copy stats, peak fitting | needed | blocked on `numerics.rs`; replace GSL with Rust/statrs |
| `som` | `vcfsom.c` | experimental SOM filter, VCF I/O | partial | deferred/out-of-scope unless tests require it |

## Plugin Surface

| Plugin group | Upstream source files | Main APIs called | `htslib-rs` status | Internal owner / blocker |
| --- | --- | --- | --- | --- |
| Tag fixers | `plugins/fill-AN-AC.c`, `fill-tags.c`, `missing2ref.c`, `tag2tag.c`, `setGT.c`, `add-variantkey.c`, `variantkey-hex.c`, `allele-length.c`, `impute-info.c`, `counts.c`, `dosage.c`, `frameshifts.c`, `remove-overlaps.c`, `fill-from-fasta.c` | typed INFO/FORMAT/GT getters, `bcf_update_*`, FASTA, buffers, filters | partial | blocked on `bcf_update_*`, `filter/`, `abuf.rs`, `vcfbuf.rs`, `reference.rs` |
| Reference fixers | `plugins/fixref.c`, `fixploidy.c` | FASTA, allele mutation, ploidy/sample handling | partial | blocked on `reference.rs`, `ploidy.rs`, `bcf_update_*` |
| Subset/split | `plugins/split.c`, `scatter.c`, `GTsubset.c`, `GTisec.c`, `isecGT.c` | sample subsets, synced readers, header translation, writer indexing | partial | blocked on `synced.rs`, `bcf_translate`, `smpl_ilist.rs` |
| Stats/reports | `plugins/smpl-stats.c`, `indel-stats.c`, `trio-stats.c`, `variant-distance.c`, `ad-bias.c`, `af-dist.c`, `check-ploidy.c`, `check-sparsity.c`, `vcf2table.c`, `vrfs.c` | typed getters, sample lists, numerics, convert formatting | partial | blocked on `numerics.rs`, `convert/`, `smpl_ilist.rs` |
| VEP-aware | `plugins/split-vep.c` | CSQ/VEP INFO parser, convert-style output, filters | partial | blocked on `convert/`, `filter/`, INFO parsing/mutation |
| Trio/pedigree | `plugins/mendelian2.c`, `trio-dnm2.c`, `trio-switch-rate.c`, `parental-origin.c` | pedigree logic, genotype likelihoods, HMM/numerics | partial | blocked on `numerics.rs`, `hmm.rs`, filter/sample helpers |
| Sample inference | `plugins/guess-ploidy.c`, `contrast.c` | ploidy, sample subsets, numeric tests | partial | blocked on `ploidy.rs`, `numerics.rs`, `smpl_ilist.rs` |
| Misc | `plugins/color-chrs.c`, `gvcfz.c`, `prune.c` | terminal color output, gVCF, VCF buffer/windowing | partial; `color-chrs` terminal behavior out-of-scope for core library | blocked on `gvcf.rs`, `vcfbuf.rs`; terminal color may stay CLI-only |

## htslib-rs Gap Rollup

The downstream gaps from this map are tracked in `htslib-rs/TODO.md` under
the bcftools-rs rollup section:

- `synced_bcf_reader` full API parity.
- `bcf_translate` beyond the current synthetic translation fixture.
- Complete `bcf_update_*` mutation primitives for all command call sites.
- Pileup iterator surface needed by bcftools `mpileup`.
- BAQ and `probaln_glocal` wiring for `bam2bcf*.c`.
- BGZF writer thread-pool wiring via `hts_set_threads`.
- Region/target arithmetic differences for `-r`/`-R` vs `-t`/`-T`.

