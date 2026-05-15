# Upstream Perl Test Status

Status values:

- `enabled`: runs in CI through `scripts/run-bcftools-test-pl.sh`.
- `rust-covered`: Rust-native tests cover a meaningful local slice, but the
  upstream Perl function is not enabled in CI yet.
- `not-yet-ported`: upstream test exists, but the required Rust subcommand or
  shared infrastructure is not complete enough to enable it.
- `partial-rust-covered`: the Rust command has meaningful unit/integration
  coverage, but the full upstream Perl slice is not enabled yet.
- `blocked-local-tool`: the Rust subcommand has partial coverage, but this
  test function also requires external tools not yet provided by CI.

| `test.pl` function | Status | Notes |
| --- | --- | --- |
| `test_vcf_head` | enabled | Plain VCF `bcftools head`; CI runs it in the enabled parity slice regex. |
| `test_vcf_head2` | enabled | Compressed VCF and BCF-pipe `head` coverage runs with Rust-backed staged `bgzip`/`tabix` harness helpers. |
| `test_index` | enabled | VCF/BCF indexing, explicit output path, and streaming index creation pass in CI. |
| `test_vcf_idxstats` | enabled | `bcftools index -s/-n` over TBI/CSI, VCF.gz, BCF, and direct CSI paths passes in CI. |
| `test_vcf_sort` | rust-covered | Rust integration coverage covers coordinate sorting, output formats, write-index, temp-prefix/compression-level shape, and Kestrel header compatibility; full Perl slice still depends on broader `query` parity. |
| `test_vcf_view` | rust-covered | Rust integration coverage includes VCF/VCF.gz/BCF I/O, sample/region/target filtering, common site filters, genotype/phasing filters, expressions over core/INFO fields, output formats, threading, and Kestrel headers; full Perl parity remains incomplete for advanced FORMAT/sample expression and structured writer semantics. |
| `test_csq` | not-yet-ported | `csq` not ported. |
| `test_csq_real` | not-yet-ported | `csq` not ported. |
| `test_gtcheck` | not-yet-ported | `gtcheck` not ported. |
| `test_mpileup` | not-yet-ported | `mpileup` not ported. |
| `test_naive_concat` | rust-covered | Rust integration coverage includes `--naive` and `--naive-force`; Perl slice not enabled yet. |
| `test_plugin_scatter` | not-yet-ported | Plugin not ported. |
| `test_plugin_split` | not-yet-ported | Plugin not ported. |
| `test_plugin_vrfs` | not-yet-ported | Plugin not ported. |
| `test_rename_chrs` | not-yet-ported | Depends on `annotate`/`query`. |
| `test_roh` | not-yet-ported | `roh` not ported. |
| `test_tabix` | enabled | VCF BGZF indexing/querying runs with Rust-backed staged `bgzip`/`tabix` harness helpers. |
| `test_trio_dnm3` | not-yet-ported | Plugin not ported. |
| `test_usage` | not-yet-ported | Harness requires `IO::Pty`; full command table is not ported. |
| `test_vcf_64bit` | not-yet-ported | Depends on broader VCF command coverage. |
| `test_vcf_annotate` | not-yet-ported | `annotate` not ported. |
| `test_vcf_call` | not-yet-ported | `call` not ported. |
| `test_vcf_call_cAls` | not-yet-ported | `call` not ported. |
| `test_vcf_check` | rust-covered | Rust `stats` has integration coverage for the local stats slice; the upstream check Perl flow is not enabled yet. |
| `test_vcf_check_merge` | not-yet-ported | Depends on `stats`/`merge`. |
| `test_vcf_concat` | rust-covered | Rust integration coverage includes same-sample concat, regions, duplicate removal, output formats, indexing, headers, threads, and Kestrel reads; ligation and full synced-reader edge cases remain. |
| `test_vcf_consensus` | not-yet-ported | `consensus` not ported. |
| `test_vcf_consensus_chain` | not-yet-ported | `consensus` not ported. |
| `test_vcf_convert` | rust-covered | Rust integration coverage includes TSV/23andMe-style conversion, gVCF expansion, Oxford GEN/HAP/HAP-LEGEND forward and reverse paths, sample selection, filters, output formats, indexing, and many upstream fixtures; full edge-case parity remains. |
| `test_vcf_convert_gvcf` | rust-covered | Rust integration coverage includes VCF/VCF.gz/BCF gVCF expansion and filter-gated expansion behavior; Perl slice not enabled yet. |
| `test_vcf_convert_hls2vcf` | rust-covered | Rust integration coverage includes HAP/LEGEND/SAMPLE back-conversion, text/VCF.gz/BCF output, indexing, and fixture parity for the current slice; haploid-missing BCF serialization remains dependency-blocked. |
| `test_vcf_convert_hs2vcf` | rust-covered | Rust integration coverage includes HAP/SAMPLE back-conversion, sample selection, sex/haploid2diploid handling, output formats, indexing, and fixture parity for the current slice; haploid-missing BCF serialization remains dependency-blocked. |
| `test_vcf_convert_tsv2vcf` | rust-covered | Rust integration coverage includes explicit-column TSV, AA/reference-derived alleles, GT sample fields, skipped malformed rows, diagnostics counters, output formats, and indexing; full 23andMe edge-case parity remains. |
| `test_vcf_filter` | rust-covered | Rust integration coverage includes VCF/VCF.gz/BCF I/O, expression filtering over the current shared filter slice, soft/hard FILTER rewriting, masks, gap filters, set-GTs, regions/targets, output formats, indexing, headers, threads, and Kestrel reads; full FORMAT/sample-vector parity remains. |
| `test_vcf_isec` | rust-covered | Rust integration coverage includes pairwise intersections, collapse modes, record/target filters, directory output, output formats, indexing, and Kestrel reads; full synced-reader parity remains. |
| `test_vcf_isec2` | rust-covered | Rust integration coverage covers the current multi-file/directory-output slice; full upstream parity remains. |
| `test_vcf_merge` | not-yet-ported | `merge` not ported. |
| `test_vcf_merge_big` | not-yet-ported | `merge` not ported. |
| `test_vcf_norm` | not-yet-ported | `norm` not ported. |
| `test_vcf_plugin` | partial-rust-covered | Static plugin registry/listing is in; `+counts`, `+missing2ref`, `+fill-AN-AC`, `+allele-length`, `+variant-distance`, `+check-ploidy`, `+tag2tag` (gl-to-pl/gp-to-gt), `+add-variantkey`, `+variantkey-hex`, `+remove-overlaps`, `+af-dist`, `+smpl-stats`, `+indel-stats`, `+ad-bias`, `+prune`, `+dosage`, `+guess-ploidy`, `+contrast`, `+fixref`, `+trio-switch-rate`, `+trio-stats`, `+mendelian2`, `+parental-origin`, `+fixploidy`, and `+GTsubset` have Rust implementations (25 of 41). `+missing2ref`/`+fill-AN-AC`/`+allele-length`, all four `+variant-distance` modes, all three `+check-ploidy` fixtures, the two integer `+tag2tag` conversions, both VariantKey plugins (full `query.add-variantkey.vcf` / `variantkey-hex.out` fixtures), `+remove-overlaps` (all six `remove-overlaps.1.*` overlap/dup/`-O t`/`--reverse` fixtures), `+af-dist` (`af-dist.out`, HWE prob + AF-deviation histograms), `+smpl-stats` (`smpl-stats.1.out`, per-sample/per-site genotype stats), `+indel-stats` (`indel-stats.1.out`, SN/DVAF/DLEN/DFRAC/NFRAC), `+ad-bias` (`ad-bias.out`, Fisher exact test on FMT/AD, both inputs), `+prune` (`prune.1.{1,2,3,4,6}.out` + `prune.2.1.out`, windowed maxAF/1st pruning **and** `calc_ld` r2/LD'/RD `-a`/`-m`/`-f` modes), `+dosage` (`dosage.{1,2,3}.out`, PL/GL/GT dosages), `+guess-ploidy` (`guess-ploidy.{PL,GL}.out`, sex inference from PL/GL likelihoods), `+contrast` (`contrast.out`/`.1.out`/`.1.1.out`/`.1.2.out`, control/case association + novel allele/genotype), `+fixref` (`fixref.{4,5,6,7}.out`, FASTA REF/ALT ref-alt/swap/flip/flip-all), `+trio-switch-rate` (`trio.out`, PED-trio phase-switch rate), `+trio-stats` (`trio-stats.out`/`trio-stats.2.out`, Mendelian/DNM/transmitted trio stats), `+mendelian2` (`mendelian.{1,3,4,6,7,8}.out`, built-in GRCh37 ruleset → per-record ploidy/inheritance, count/annotate/delete/list modes), `+parental-origin` (`parental-origin.{1,2,3,4,5}.out`, DEL/DUP genotype-likelihood parental-origin models + local `kf_betai` binomial filters), `+fixploidy` (`fixploidy.out`, `ploidy.c`-subset region/sex ploidy query + bcf GT expand/in-place rewrite), and `+GTsubset` (`view.GTsubset.{NA1,NA1NA2,NA1NA2NA3}.out`, exclusive-shared-genotype site filter on bcf-encoded allele/phase ints) match their upstream `*.out` fixtures byte-for-byte. The umbrella Perl slice stays disabled until enough plugins land. |
| `test_vcf_query` | rust-covered | Rust integration coverage includes list-samples, sample selection, headers, regions/targets, record and sample filters, and many formatter tokens/functions; full shared formatter parity remains. |
| `test_vcf_regions` | rust-covered | Rust `query`/`view`/`filter`/`isec` have local region and target coverage, but full upstream region/target semantics are not complete. |
| `test_vcf_reheader` | rust-covered | Rust integration coverage includes header replacement, sample renaming, FAI contig updates, stdin, BCF output, in-place BCF, and threads; `test_rename_chrs` still depends on `annotate`/full `query`. |
| `test_vcf_stats` | rust-covered | Rust integration coverage includes the current stats report slice, region filtering, sample selection, two-input comparison, computed TYPE expressions, and upstream fixture checks; full stats parity remains. |
