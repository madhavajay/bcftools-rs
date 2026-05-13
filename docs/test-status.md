# Upstream Perl Test Status

Status values:

- `enabled`: runs in CI through `scripts/run-bcftools-test-pl.sh`.
- `not-yet-ported`: upstream test exists, but the required Rust subcommand or
  shared infrastructure is not complete enough to enable it.
- `blocked-local-tool`: the Rust subcommand has partial coverage, but this
  test function also requires external tools not yet provided by CI.

| `test.pl` function | Status | Notes |
| --- | --- | --- |
| `test_vcf_head` | enabled | Plain VCF `bcftools head`; CI runs `-f '^test_vcf_head$'`. |
| `test_vcf_head2` | blocked-local-tool | Exercises compressed VCF/BCF paths and currently needs external `bgzip`/`tabix` availability in the Perl harness. |
| `test_index` | not-yet-ported | Rust integration coverage exists; Perl slice not enabled yet. |
| `test_vcf_idxstats` | not-yet-ported | Rust integration coverage exists; Perl slice not enabled yet. |
| `test_vcf_sort` | not-yet-ported | Rust integration coverage exists for the current subset; full Perl slice depends on `query`. |
| `test_vcf_view` | not-yet-ported | Current Rust `view` lacks expression filtering and sample subsetting. |
| `test_csq` | not-yet-ported | `csq` not ported. |
| `test_csq_real` | not-yet-ported | `csq` not ported. |
| `test_gtcheck` | not-yet-ported | `gtcheck` not ported. |
| `test_mpileup` | not-yet-ported | `mpileup` not ported. |
| `test_naive_concat` | not-yet-ported | `concat` not ported. |
| `test_plugin_scatter` | not-yet-ported | Plugin not ported. |
| `test_plugin_split` | not-yet-ported | Plugin not ported. |
| `test_plugin_vrfs` | not-yet-ported | Plugin not ported. |
| `test_rename_chrs` | not-yet-ported | Depends on `annotate`/`query`. |
| `test_roh` | not-yet-ported | `roh` not ported. |
| `test_tabix` | blocked-local-tool | Rust `tabix` has integration coverage; Perl slice requires external `bgzip` on `PATH`. |
| `test_trio_dnm3` | not-yet-ported | Plugin not ported. |
| `test_usage` | not-yet-ported | Harness requires `IO::Pty`; full command table is not ported. |
| `test_vcf_64bit` | not-yet-ported | Depends on broader VCF command coverage. |
| `test_vcf_annotate` | not-yet-ported | `annotate` not ported. |
| `test_vcf_call` | not-yet-ported | `call` not ported. |
| `test_vcf_call_cAls` | not-yet-ported | `call` not ported. |
| `test_vcf_check` | not-yet-ported | Depends on `stats`. |
| `test_vcf_check_merge` | not-yet-ported | Depends on `stats`/`merge`. |
| `test_vcf_concat` | not-yet-ported | `concat` not ported. |
| `test_vcf_consensus` | not-yet-ported | `consensus` not ported. |
| `test_vcf_consensus_chain` | not-yet-ported | `consensus` not ported. |
| `test_vcf_convert` | not-yet-ported | `convert` not ported. |
| `test_vcf_convert_gvcf` | not-yet-ported | `convert`/gVCF helpers not ported. |
| `test_vcf_convert_hls2vcf` | not-yet-ported | `convert` not ported. |
| `test_vcf_convert_hs2vcf` | not-yet-ported | `convert` not ported. |
| `test_vcf_convert_tsv2vcf` | not-yet-ported | `convert`/TSV helpers not ported. |
| `test_vcf_filter` | not-yet-ported | Filter expression engine not ported. |
| `test_vcf_isec` | not-yet-ported | `isec` not ported. |
| `test_vcf_isec2` | not-yet-ported | `isec` not ported. |
| `test_vcf_merge` | not-yet-ported | `merge` not ported. |
| `test_vcf_merge_big` | not-yet-ported | `merge` not ported. |
| `test_vcf_norm` | not-yet-ported | `norm` not ported. |
| `test_vcf_plugin` | not-yet-ported | Plugin registry and plugin implementations not ported. |
| `test_vcf_query` | not-yet-ported | `query`/convert formatter not ported. |
| `test_vcf_regions` | not-yet-ported | Depends on `query` and full region/target semantics. |
| `test_vcf_reheader` | not-yet-ported | `reheader` not ported. |
| `test_vcf_stats` | not-yet-ported | `stats` not ported. |
