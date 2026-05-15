# Upstream Perl Test Status

Status values:

- `enabled`: runs in CI through `scripts/run-bcftools-test-pl.sh`.
- `not-yet-ported`: upstream test exists, but the required Rust subcommand or
  shared infrastructure is not complete enough to enable it.
- `partial-rust-covered`: the Rust command has meaningful unit/integration
  coverage, but the full upstream Perl slice is not enabled yet.
- `blocked-local-tool`: the Rust subcommand has partial coverage, but this
  test function also requires external tools not yet provided by CI.

| `test.pl` function | Status | Notes |
| --- | --- | --- |
| `test_vcf_head` | enabled | Plain VCF `bcftools head`; CI runs `-f '^test_vcf_head$'`. |
| `test_vcf_head2` | blocked-local-tool | Exercises compressed VCF/BCF paths and currently needs external `bgzip`/`tabix` availability in the Perl harness. |
| `test_index` | partial-rust-covered | Rust integration coverage exists for BCF/VCF.gz CSI/TBI builds, stdin indexing, stats, large-coordinate CSI queries, and option validation; Perl slice not enabled yet. |
| `test_vcf_idxstats` | partial-rust-covered | Rust integration coverage exists for `index -s`/`-n`; Perl slice not enabled yet. |
| `test_vcf_sort` | partial-rust-covered | Rust integration coverage exists for coordinate sorting, compressed writes, indexing, temp-run spill, Kestrel headers, and threading; full Perl slice still depends on broader command parity. |
| `test_vcf_view` | partial-rust-covered | Rust integration coverage exists for VCF/VCF.gz/BCF I/O, sample/region/target filtering, many simple filters, expressions, Kestrel headers, and threaded writes; full Perl slice still needs complete expression and structured-path parity. |
| `test_csq` | not-yet-ported | `csq` not ported. |
| `test_csq_real` | not-yet-ported | `csq` not ported. |
| `test_gtcheck` | not-yet-ported | `gtcheck` not ported. |
| `test_mpileup` | not-yet-ported | `mpileup` not ported. |
| `test_naive_concat` | partial-rust-covered | Rust `concat` has naive concat coverage, but full Perl slice is not enabled yet. |
| `test_plugin_scatter` | not-yet-ported | Plugin not ported. |
| `test_plugin_split` | not-yet-ported | Plugin not ported. |
| `test_plugin_vrfs` | not-yet-ported | Plugin not ported. |
| `test_rename_chrs` | not-yet-ported | Depends on `annotate`/`query`. |
| `test_roh` | not-yet-ported | `roh` not ported. |
| `test_tabix` | blocked-local-tool | Rust `tabix` has integration coverage for preset TBI/CSI builds and queries; Perl slice requires external `bgzip` on `PATH`. |
| `test_trio_dnm3` | not-yet-ported | Plugin not ported. |
| `test_usage` | not-yet-ported | Harness requires `IO::Pty`; full command table is not ported. |
| `test_vcf_64bit` | not-yet-ported | Depends on broader VCF command coverage. |
| `test_vcf_annotate` | not-yet-ported | `annotate` not ported. |
| `test_vcf_call` | not-yet-ported | `call` not ported. |
| `test_vcf_call_cAls` | not-yet-ported | `call` not ported. |
| `test_vcf_check` | partial-rust-covered | Depends on `stats`; Rust `stats` has substantial local coverage but full Perl slice is not enabled yet. |
| `test_vcf_check_merge` | not-yet-ported | Depends on full `stats`/`merge`; `merge` is not ported. |
| `test_vcf_concat` | partial-rust-covered | Rust `concat` covers same-sample vertical concat, naive concat, duplicate handling, region restriction, indexing, Kestrel headers, and threading; full Perl slice is not enabled yet. |
| `test_vcf_consensus` | not-yet-ported | `consensus` not ported. |
| `test_vcf_consensus_chain` | not-yet-ported | `consensus` not ported. |
| `test_vcf_convert` | partial-rust-covered | Rust `convert` has broad GEN/SAMPLE, HAP/SAMPLE, HAP/LEGEND/SAMPLE, TSV, gVCF, BCF stdin, output/indexing, and fixture parity coverage; full Perl slice is not enabled yet. |
| `test_vcf_convert_gvcf` | partial-rust-covered | Rust gVCF conversion has FASTA-backed reference-block expansion and fixture-pipe coverage; advanced expression parity remains. |
| `test_vcf_convert_hls2vcf` | partial-rust-covered | Rust HAP/LEGEND/SAMPLE back-conversion has text and BCF output coverage; haploid missing `GT=.` BCF serialization remains dependency-blocked. |
| `test_vcf_convert_hs2vcf` | partial-rust-covered | Rust HAP/SAMPLE back-conversion has text and BCF output coverage; haploid missing `GT=.` BCF serialization remains dependency-blocked. |
| `test_vcf_convert_tsv2vcf` | partial-rust-covered | Rust TSV/23andMe conversion has fixture parity and `-Ou | view` coverage; full diagnostics and edge-case parity remain. |
| `test_vcf_filter` | partial-rust-covered | Rust `filter` has text-mode expression, mask, soft-filter, set-GTs, region/target, index, Kestrel, and threading coverage; full FORMAT/sample-vector and structured-path parity remain. |
| `test_vcf_isec` | partial-rust-covered | Rust `isec` has text-backed set intersection, complement, prefix output, collapse modes, targets, regions, and VCF/BCF record-output coverage; full synced-reader parity remains. |
| `test_vcf_isec2` | partial-rust-covered | Same Rust `isec` coverage as `test_vcf_isec`; full upstream multi-file edge parity remains. |
| `test_vcf_merge` | not-yet-ported | `merge` not ported. |
| `test_vcf_merge_big` | not-yet-ported | `merge` not ported. |
| `test_vcf_norm` | not-yet-ported | `norm` not ported. |
| `test_vcf_plugin` | not-yet-ported | Plugin registry and plugin implementations not ported. |
| `test_vcf_query` | partial-rust-covered | Rust `query` has list-samples, sample selection, region/target, expression, formatter, numeric function, and PBINOM fixture coverage; full formatter and sample-vector parity remain. |
| `test_vcf_regions` | partial-rust-covered | Depends on `query` and full region/target semantics; Rust `view`/`query` have substantial POS-based region/target coverage but not full indexed overlap parity. |
| `test_vcf_reheader` | partial-rust-covered | Rust `reheader` has header replacement, sample rename, FAI contig update, stdin, BCF, in-place, and threading coverage; rename-chrs still depends on `annotate`/full `query`. |
| `test_vcf_stats` | partial-rust-covered | Rust `stats` has substantial single-input and pairwise text-backed section coverage; full indexed synced-reader and edge-case parity remain. |
