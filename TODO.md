# bcftools-rs Fixes Needed For BioScript VNtyper

This file tracks `bcftools-rs` changes currently needed by the BioScript VNtyper
port. These should be fixed in `bcftools-rs` rather than worked around in
BioScript unless we explicitly decide otherwise.

## 1. Accept Kestrel Java-style `VCF4.2` headers ✅

Resolved in `bcftools-rs` (no submodule changes). A `Read` adapter and a
companion in-place text normalizer in `crates/bcftools-rs/src/vcf_compat.rs`
rewrite a non-canonical `##fileformat=VCF<x>.<y>` first line to
`##fileformat=VCFv<x>.<y>` and emit the upstream-style warning
(`[W::bcf_get_version] Couldn't get VCF version, considering as <ver>`).
Wired into:
- `bcftools-rs sort` text + BGZF/gzip read paths
  (`crates/bcftools-rs/src/commands/sort.rs`).
- `bcftools-rs head` `-s/-n` record-emitting paths
  (`crates/bcftools-rs/src/commands/head.rs::write_n_records`). Header-text
  output already preserves raw bytes via line-by-line read.
- `bcftools-rs view` text passthrough, filtered passthrough, structured
  VCF/BCF writer paths, and `read_header`
  (`crates/bcftools-rs/src/commands/view.rs`).
- `bcftools-rs query` `-l/-s/-S` sample-list parser path
  (`crates/bcftools-rs/src/commands/query.rs::header_sample_names_from_path`).
  Query's text-mode formatter never validates the fileformat line, so
  Kestrel input naturally passes through.

`bcftools-rs reheader` reads VCF text line-by-line without invoking the
strict parser, so Kestrel headers already round-trip and no wrapper was
needed there.

Covered by:
- `crates/bcftools-rs/src/vcf_compat.rs` unit tests (`tests` module).
- `crates/bcftools-rs/tests/sort.rs::sort_accepts_kestrel_non_canonical_fileformat_header`.
- `crates/bcftools-rs/tests/sort.rs::sort_accepts_kestrel_header_with_compressed_write_index`.
- `crates/bcftools-rs/tests/sort.rs::sort_does_not_warn_for_canonical_fileformat_header`.
- `crates/bcftools-rs/tests/head.rs::head_with_s_accepts_kestrel_non_canonical_fileformat_header`.
- `crates/bcftools-rs/tests/view.rs::view_accepts_kestrel_non_canonical_fileformat_header`.

### Original notes (kept for context)

### Problem

Java Kestrel emits VCF files with this first line:

```text
##fileformat=VCF4.2
```

That is not the canonical VCF spelling, which is:

```text
##fileformat=VCFv4.2
```

However, current upstream `bcftools 1.23.1` still accepts the Kestrel form with
a warning and treats it as VCF 4.2.

### Upstream bcftools behavior

Verified with:

```bash
ports/vntyper/test-data/tools/local/bin/bcftools --version
```

Output:

```text
bcftools 1.23.1
Using htslib 1.23.1
```

Reproduction:

```bash
tmpdir=$(mktemp -d)
cat > "$tmpdir/kestrel-style.vcf" <<'VCF'
##fileformat=VCF4.2
##contig=<ID=chr1,length=10>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	C	.	PASS	.
VCF

ports/vntyper/test-data/tools/local/bin/bcftools sort \
  "$tmpdir/kestrel-style.vcf" \
  -o "$tmpdir/kestrel-style.sorted.vcf.gz" \
  -W \
  -O z
```

Observed upstream behavior:

```text
[W::bcf_get_version] Couldn't get VCF version, considering as 4.2
```

The command exits `0`.

### Current bcftools-rs behavior

Reproduction:

```bash
tmpdir=$(mktemp -d)
cat > "$tmpdir/kestrel-style.vcf" <<'VCF'
##fileformat=VCF4.2
##contig=<ID=chr1,length=10>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	C	.	PASS	.
VCF

cd vendor/rust/bcftools-rs
CC=cc AR=ar cargo run -q -p bcftools-rs-cli -- sort \
  "$tmpdir/kestrel-style.vcf" \
  -o "$tmpdir/kestrel-style.sorted.vcf" \
  -O v
```

Observed current behavior:

```text
[E::main_vcfsort] invalid record
```

The command exits `1`.

### Expected bcftools-rs behavior

`bcftools-rs sort` should match upstream `bcftools` for this case:

- Accept `##fileformat=VCF4.2`.
- Warn that the version could not be parsed, if warning parity is feasible.
- Treat the file as VCF 4.2.
- Continue sorting/compressing/indexing normally.
- Exit `0` for otherwise valid input.

### Likely fix location

This probably belongs in the VCF reader/version handling layer used by
`bcftools-rs sort`, not in BioScript.

The Rust implementation appears stricter than HTSlib here, likely because it is
parsing through a strict VCF parser path instead of reproducing HTSlib's
`bcf_get_version` fallback behavior.

Potential approaches:

- Add a compatibility path in `bcftools-rs sort` before strict VCF parsing that
  treats `##fileformat=VCF4.2` as `VCFv4.2`.
- Preferably, add the fallback in the shared `htslib-rs`/VCF compatibility layer
  if other commands also parse VCF headers through the same path.

### Test to add

Add a reduced test in `bcftools-rs` that writes the VCF shown above and asserts:

- `bcftools-rs sort input.vcf -o output.vcf -O v` exits successfully.
- The output VCF exists and contains the record.
- If warning capture is practical, the warning matches upstream intent:
  unable to parse version, considering as 4.2.

Also keep an existing canonical-header test:

```text
##fileformat=VCFv4.2
```

so the compatibility path does not regress normal VCF parsing.

### BioScript VNtyper impact

BioScript native VNtyper runs:

```text
kestrel.run_native -> bcftools.sort_native -> bcftools.index_native
```

Kestrel Java-compatible output currently uses `##fileformat=VCF4.2`, so
`bcftools-rs` rejecting this blocks sorting/indexing raw Kestrel VCF output.

BioScript could normalize the header before calling `bcftools-rs`, but that
would hide a real upstream parity gap. The correct parity behavior is for
`bcftools-rs` to accept the file the same way upstream `bcftools` does.

# TODO: Port bcftools to Pure Rust

Goal: build a pure Rust replacement for the `bcftools` C program with full subcommand and plugin parity, then port and pass the upstream `test/test.pl` suite plus add Rust-native unit/integration tests. Implementation routes through `htslib-rs` (sibling submodule) for HTSlib-shaped APIs and may use `noodles` only where there is no HTSlib analogue.

Current goal: keep momentum inside `bcftools-rs` only. If a TODO item requires changes to underlying libraries (`htslib-rs`, `noodles`, or their submodules), move that dependency work to the end of this file under the rolling dependency-blocker list, continue with other `bcftools-rs` items that can be completed locally, and then stop once the remaining work is blocked. Do not change the underlying libraries during this goal.

PR workflow (locked in 2026-05-15): land one PR at a time. Open a single
focused branch, run the Rust gate (`cargo fmt --all --check`,
`cargo clippy --workspace --all-targets -- -D warnings`,
`cargo test --workspace`) locally, push, wait for `Rust tests` and
`bcftools Perl parity tests` to go green on GitHub, then squash-merge with
`gh pr merge <N> --squash --delete-branch`, pull `main`, and only then start
the next slice. Do **not** stack multiple open PRs against `main` again — the
stack landed 2026-05-15 generated cascading `TODO.md`/`docs/test-status.md`/
`commands/mod.rs` conflicts that all had to be hand-resolved.

Latest landed progress:

- 2026-05-17: PR #107 (`progress/merge-format-vector-text-fixture`, merge
  commit `22b9403`) added FORMAT-key union plus Number=A/R vector remapping in
  the text-backed `merge` slice, with byte-for-byte `merge.6.out` coverage.
- 2026-05-17: PR #105 (`progress/merge-star-ref-extension-text-fixture`,
  merge commit `9f81a10`) added sampled star-allele REF-extension merging in
  the text-backed `merge` slice, with byte-for-byte `merge.5.out` coverage.
- 2026-05-17: PR #103 (`progress/merge-multiallelic-subset-text-fixture`,
  merge commit `208e263`) added sampled same-ref ALT-subset merging for
  `merge --merge none|both`, with byte-for-byte
  `merge.multiallelics.1.1.out` coverage.
- 2026-05-17: PR #101 (`progress/merge-fileformat-max-text-fixture`, merge
  commit `d7338d4`) preserved the highest input VCF fileformat in the
  text-backed merge header, with byte-for-byte `merge.symbolic.1.1.out`
  coverage.
- 2026-05-17: PR #99 (`progress/merge-info-af-join-text-fixture`, merge
  commit `e1f4e91`) added narrow `-i AF:join` INFO-rule joining for
  sites-only ALT-union records, with byte-for-byte `merge.join.1.out`
  coverage.
- 2026-05-17: PR #97 (`progress/merge-nonref-symbolic-text-fixture`, merge
  commit `6d1c20e`) added narrow `<NON_REF>` symbolic allele union for
  `merge --merge none|both`, with byte-for-byte `merge.12.1.out` coverage.
- 2026-05-17: PR #95 (`progress/merge-ad-vector-text-fixture`, merge commit
  `a8c0c20`) added default-mode sampled same-position allele union with
  FORMAT/AD vector expansion and byte-for-byte `merge.11.1.out` coverage.
- 2026-05-17: PR #93 (`progress/merge-mode-both-text-fixtures`, merge commit
  `b1e7cb0`) added narrow sampled same-position allele unions for
  `merge -m both` and `merge -m snp-ins-del`, including GT allele remapping
  and byte-for-byte `merge.10.{2,3}.out` coverage.
- 2026-05-17: PR #91 (`progress/merge-mode-none-text-fixture`, merge commit
  `78ce1f5`) added a narrow `merge -m none` text-backed slice that keeps
  same-position conflicting records as separate rows, with byte-for-byte
  `merge.10.1.out` coverage.
- 2026-05-17: PR #89 (`progress/merge-force-samples-prefix-chain`, merge
  commit `ae51059`) matched upstream repeated duplicate-sample prefixing under
  `merge --force-samples`, with byte-for-byte `merge.9.out` fixture coverage.
- 2026-05-17: PR #87 (`progress/merge-info-rules-alt-union`, merge commit
  `b151b0c`) added narrow `AC:sum`/`AN:sum` INFO-rule handling for the
  sampled plus sites-only same-locus ALT-union `merge` slice, with
  byte-for-byte `merge.9.2.out` text coverage.
- 2026-05-17: PR #85 (`progress/merge-sampled-sites-only-alt`, merge commit
  `a10673b`) added sampled plus sites-only same-locus ALT union in `merge`,
  preserving sampled columns and upstream default AC/AN behavior with
  byte-for-byte `merge.9.1.out` text coverage.
- 2026-05-17: PR #83 (`progress/merge-sites-only-alt-union`, merge commit
  `1290900`) added sites-only same-locus ALT union in `merge`, including
  AC/AN summing, `-i AN:sum,AC:sum` command-shape compatibility for the
  current slice, and byte-for-byte `merge.8.out` text coverage.
- 2026-05-17: PR #81 (`progress/merge-force-single`, merge commit
  `d4491ca`) added `merge --force-single` command-shape compatibility and
  byte-for-byte `merge.LPL.0.out` text coverage while preserving default
  single-input rejection.
- 2026-05-17: PR #79 (`progress/merge-missing-sites-slice`, merge commit
  `8fb3c4d`) extended the local text-backed `merge` slice to simple sorted
  exact-site unions with missing FORMAT sample synthesis, `--no-index`
  command-shape compatibility, upstream-style PASS header injection, and
  byte-for-byte `merge.noidx.abc.out` coverage.
- 2026-05-17: PR #77 (`progress/filter-format-aggregates`, merge commit
  `5acf26b`) landed FORMAT/AO aggregate `-S .` text parity in `filter` —
  site-level `MIN`/`MAX`/`SUM`/`AVG`/`MEDIAN`/`STDEV` versus sample-level
  `SMPL_*`/`s*` behavior, with byte-for-byte `filter.{30,31,32,33,34,36}.out`
  coverage.
- 2026-05-17: PR #75 (`progress/filter-missing-subscripts`, merge commit
  `41c7b73`) tightened text-backed `filter` parity for missing subscript /
  QUAL comparisons and SnpGap type-list matching — byte-for-byte
  `filter.{24,25,26}.out` plus `filter.29.out` text coverage, including
  symbolic deletion spans from `INFO/END`.
- 2026-05-17: PR #73 (`progress/filter-format-indexes`, merge commit
  `f4803d0`) landed FORMAT sample/value subscripts in shared filter
  evaluation — `TAG[:value]`, `TAG[sample:]`, `TAG[sample:value]`, and
  existing `TAG[*]` behavior over text-backed FORMAT vectors, with
  byte-for-byte `filter.{20,21,22,23}.out` coverage.
- 2026-05-17: PR #72 (`progress/filter-gt-literals`, merge commit
  `a08cf9f`) landed bcftools GT class literals in shared filter evaluation —
  `A`, `R`, `RR`, `RA`, `AR`, `AA`, `aA`, `Aa`, `HOM`, `HET`, and `HAP`,
  with byte-for-byte `filter.{12,13,14,15,16,17,18,19}.out` coverage.
- 2026-05-17: PR #71 (`progress/filter-engine-slice`, merge commit
  `0707e1c`) landed a shared filter-engine / `filter` slice — upstream
  single `&`/`|` boolean operators, top-level `FMT/`/`FORMAT/` sample-vector
  predicates, joined filter short options, and byte-for-byte `filter.2.out`
  coverage.
- 2026-05-17: PR #70 (`progress/check-sparsity`, merge commit `f7a5395`)
  landed `+check-sparsity` — local text-backed sparse-sample reporting by
  contig or requested region, with `-n`, `-r`, and `-R` support plus
  synthetic Rust coverage because upstream has no dedicated fixture row.
- 2026-05-17: PR #69 (`progress/frameshifts`, merge commit `e7be6be`)
  landed `+frameshifts` — local OOF annotation for simple indel alleles
  using simple BED-like exon files, with synthetic Rust integration coverage
  because upstream has no dedicated fixture row.
- 2026-05-17: PR #68 (`progress/isecgt`, merge commit `549ad56`)
  landed `+isecGT` — local text-backed two-input genotype comparison,
  sample-name remapping, and missing-GT rewrite, with synthetic Rust
  integration coverage because upstream has no dedicated fixture row.
- 2026-05-17: PR #67 (`progress/split`, merge commit `1611123`)
  landed `+split` — per-sample / per-group VCF splitting by default,
  `-S` sample files, `-G` group files, and `-k` INFO/FORMAT projection,
  with VCF text and BGZF VCF output; byte-for-byte against
  `split.1.{1,2,3,7}.out` and `split.2.1.out` (sort-dir / `query -l` /
  `view -H` harness).
- 2026-05-17: PR #66 (`progress/scatter`, merge commit `42ae63d`)
  landed `+scatter` — split a VCF into multiple VCFs by `-n` chunks or
  `-s`/`-S` regions (`-x` extra, `-p` prefix), byte-for-byte against
  `scatter.1.{1,2,3}.out` (sort-dir / `cat` / `grep -v ^##` harness).
- 2026-05-16: PR #65 (`progress/fill-from-fasta`, merge commit
  `c4cc5cd`) landed `+fill-from-fasta` (`-c REF` modes) — REF fill
  from a FASTA reference, byte-for-byte against `ref.out` / `aa.2.out`
  (the `-i` `aa.out` row deferred to the filter engine).
- 2026-05-16: PR #64 (`progress/gtisec`, merge commit `68cc665`)
  landed `+GTisec` — genotype-intersection subset counts in
  banker's-sequence order, byte-for-byte against
  `view.GTisec.{,H,Hm,Hmv,Hv,m,mv,v}.out`.
- 2026-05-16: PR #63 (`progress/gtsubset`, merge commit `7a0f7e0`)
  landed `+GTsubset` — exclusive-shared-genotype site filter on
  bcf-encoded allele/phase ints, byte-for-byte against
  `view.GTsubset.{NA1,NA1NA2,NA1NA2NA3}.out`.
- 2026-05-16: PR #62 (`progress/fixploidy`, merge commit `b4bc138`)
  landed `+fixploidy` — `ploidy.c`-subset region/sex ploidy query +
  bcf GT expand/in-place rewrite, byte-for-byte against
  `fixploidy.out`.
- 2026-05-16: PR #61 (`progress/parental-origin`, merge commit
  `bbae5f4`) landed `+parental-origin` — DEL/DUP genotype-likelihood
  parental-origin models + a local `kf_betai` port for the binomial
  parental-het filters, byte-for-byte against
  `parental-origin.{1,2,3,4,5}.out`.
- 2026-05-16: PR #60 (`progress/mendelian2`, merge commit `3dde694`)
  landed `+mendelian2` — built-in default GRCh37 ruleset → per-record
  region ploidy/inheritance, haploid + diploid consistency branches,
  count/annotate/delete/list modes, byte-for-byte against
  `mendelian.{1,3,4,6,7,8}.out`.
- 2026-05-16: PR #59 (`progress/trio-stats`, merge commit `ca2221f`)
  landed `+trio-stats` — PED trios, `bcf_calc_ac`, per-trio
  Mendelian-error / DNM / transmitted-doubleton classification, the
  `-a` deferral and interleaved `-d` debug dump, byte-for-byte against
  `trio-stats.out` / `trio-stats.2.out`.
- 2026-05-16: PR #58 (`progress/trio-switch-rate`, merge commit
  `5446b0e`) landed `+trio-switch-rate` + a reusable PED parser,
  byte-for-byte against `trio.out`.
- 2026-05-16: PR #57 (`progress/fixref`, merge commit `8806ceb`) landed
  `+fixref` — FASTA REF/ALT ref-alt/swap/flip/flip-all, byte-for-byte
  against `fixref.{4,5,6,7}.out`.
- 2026-05-16: PR #56 (`progress/contrast`, merge commit `8b2b65a`) landed
  `+contrast` — control/case PASSOC/FASSOC/NASSOC/NOVELAL/NOVELGT,
  byte-for-byte against `contrast.out`/`.1.out`/`.1.1.out`/`.1.2.out`.
- 2026-05-16: PR #55 (`progress/prune-ld`, merge commit `f4996b0`) landed
  the `+prune` LD `-a`/`-m`/`-f` modes (`calc_ld` + HTSlib `kputd`),
  byte-for-byte against `prune.1.{1,2,3}.out` / `prune.2.1.out`.
- 2026-05-15: PR #54 (`progress/guess-ploidy`, merge commit `901d5a1`)
  landed `+guess-ploidy` — PL/GL/GT haploid/diploid log-likelihood sex
  inference in `f64`, byte-for-byte against `guess-ploidy.{PL,GL}.out`.
- 2026-05-15: PR #53 (`progress/dosage`, merge commit `6e35df2`) landed
  `+dosage` — PL/GL/GT diploid likelihood/genotype dosages in `f32`,
  byte-for-byte against `dosage.{1,2,3}.out`.
- 2026-05-15: PR #52 (`progress/prune`, merge commit `c4ecd2e`) landed
  the `+prune` window subset — the `vcfbuf` windowed `_prune_sites`
  `1st`/`maxAF` modes, byte-for-byte against `prune.1.4.out`/
  `prune.1.6.out`.
- 2026-05-15: PR #51 (`progress/ad-bias`, merge commit `9be9c42`) landed
  `+ad-bias` (report mode) — Fisher's exact test on FORMAT/AD via the
  HTSlib `kfunc.c` port, byte-for-byte against `ad-bias.out` for two
  inputs.
- 2026-05-15: PR #50 (`progress/indel-stats`, merge commit `9139169`)
  landed `+indel-stats` (no-PED default) — SN/DVAF/DLEN/DFRAC/NFRAC with
  the FORMAT/AD VAF + minor-allele-fraction analysis, byte-for-byte
  against `indel-stats.1.out`.
- 2026-05-15: PR #49 (`progress/smpl-stats`, merge commit `3864e03`)
  landed `+smpl-stats` (default "all" filter) — per-sample/per-site
  genotype stats with `bcf_calc_ac` + the `bcf_acgt2int` ts/tv walk,
  byte-for-byte against `smpl-stats.1.out`.
- 2026-05-15: PR #48 (`progress/af-dist`, merge commit `c14c442`) landed
  `+af-dist` with the `bin.c` histogram port (`f32` binning), byte-for-byte
  against `af-dist.out`.
- 2026-05-15: PR #47 (`progress/remove-overlaps`, merge commit `25ecebf`)
  landed `+remove-overlaps` — a faithful port of the `vcfbuf`
  `MARK_OVERLAP`/`MARK_DUP` streaming state machine plus the
  `remove-overlaps.c` driver (`-m overlap|dup`, `-M TAG`, `--reverse`,
  `-O t`), byte-for-byte against all six `remove-overlaps.1.*` fixtures.
- 2026-05-15: PR #46 (`progress/variantkey-plugins`, merge commit
  `07cebd2`) landed `+add-variantkey` and `+variantkey-hex` with a shared
  faithful MIT VariantKey algorithm port, byte-for-byte against the full
  `query.add-variantkey.vcf` / `variantkey-hex.out` fixtures.
- 2026-05-15: PR #45 (`progress/todo-batch`, merge commit `3572038`)
  squash-landed the 7-plugin batch onto `main`: in-process
  `counts`, `missing2ref`, `fill-AN-AC`, `allele-length`,
  `variant-distance`, `check-ploidy`, and `tag2tag` (gl-to-pl/gp-to-gt),
  every one with an upstream `*.out` fixture byte-for-byte verified, plus
  the plugin output writer and dispatcher wiring.
- 2026-05-15: PR #41 (`progress/merge-same-site-slice`, merge commit
  `7543a42`) added the first local-text `bcftools merge` slice: same-position
  record merging across VCF/VCF.gz/BCF inputs with identical fixed site fields,
  duplicate-sample-name rejection unless `--force-samples` prefixes later
  inputs, `-l`/`--file-list`, `-o`, `-O u|b|v|z`, `-m TYPE` accepted for
  command-shape compat, `--no-version`, BGZF and BCF write paths, Kestrel-
  tolerant text reads, and CLI dispatcher wiring; 7 integration tests in
  `crates/bcftools-rs/tests/merge.rs`.
- 2026-05-15: PRs #10-#40 (30 PRs) landed in a single batch squash-merge onto
  `main`. Coverage spans `concat` overlap guard, `filter` mask-file/median/
  sample-fraction functions, `isec`/`stats` class-aware collapse modes,
  `query` allow-undef-tags/TGT/`FMT/AD[:GT]`, `reheader` `--samples-list`,
  `head` BGZF VCF records, `convert` attached FASTA reference/unsupported-tag
  diagnostics, `sort` attached write-index + Perl parity options (`-m 0`,
  `-Ob`, `view` pipe), `tabix` long aliases, `index` extra-input rejection +
  Perl parity, `view` numeric `-O0`-`-O9` + options-after-input + 64-bit
  text-output parity, `concat --naive` Perl parity, `stats` computed TYPE
  filters, the static plugin registry listing, first `consensus` FASTA-ALT
  slice, `annotate --rename-chrs` slice, `norm -d`/`--rm-dup` slice, plus
  refreshed README/test-status docs.
- 2026-05-14: PR #9 (`progress/convert-fixture-parity-2`, merge commit
  `05f3c18`) landed another convert parity slice after PR #8: more upstream
  GEN/SAMPLE, HAP/SAMPLE, and HAP/LEGEND/SAMPLE fixture-output parity, the
  upstream `-h` alias for HAP/LEGEND/SAMPLE output, single-precision PL/GL
  probability normalization, haploid missing HAP output parity, harness-style
  BCF stdin input for forward GEN/SAMPLE, HAP/SAMPLE, and HAP/LEGEND/SAMPLE
  output modes, upstream-style `--tsv2vcf -Ou | view` fixture pipes,
  upstream-style reverse GEN/SAMPLE `-Ou | view` fixture coverage, the
  whole-project progress estimate, and the BCF serialization blocker note for
  HAP-family reverse `-Ou` pipes.
- 2026-05-14: PR #8 (`progress/todo-local-bcftools-parity`, merge commit
  `8742124`) landed the first broad local-only parity batch for `concat`,
  `convert`, `filter`, `isec`, and `stats`, plus the dispatcher exports,
  command integration tests, TSV-to-VCF ALT normalization, and the snapshot
  coverage notes below.
- Validation before each merge: `cargo fmt --all --check`,
  `cargo clippy --workspace --all-targets -- -D warnings`,
  `cargo test --workspace`, plus GitHub CI Rust tests and bcftools Perl parity
  tests. Run these as **separate, fully-completed commands** (clippy, then
  fmt, then test) — concurrent `cargo` invocations share the `target/` lock
  and report stale green results that fail CI. Per-suite test counts are kept
  current in each command/plugin snapshot bullet rather than enumerated here
  (that enumeration drifted repeatedly); the workspace is green as of the
  latest merged commit on `main` (`22b9403`) (~220 lib unit tests plus per-command
  and per-plugin integration suites).
- No code slice is in flight after PR #107. The next branch should pick one
  local-only item from the queue below, keep the one-branch rule, run the
  full local gate, and wait for both required GitHub checks before merge.
- Next local-only queue:
  continue extending the `merge` slice toward full synced-reader alignment,
  allele unification, and `-m none|snps|indels|both|all|id`; deepen the
  `consensus`, `annotate`, and `norm` first slices; continue tightening `concat`, `filter`, `stats`,
  `isec`, `query`, `view`, `reheader`, and `convert` edge cases that do not
  require changes in `htslib-rs`, `noodles`, or their submodules. Remaining
  `convert` HAP/SAMPLE and gVCF `-Ou` pipe gaps are blocked on BCF writer
  support for haploid missing `GT=.` and the out-of-range/missing typed-value
  blockers listed at the end of this file.

Subcommand coverage at a glance (CLI dispatcher state on `main`):

| Subcommand | Status | Module / notes |
| --- | --- | --- |
| `annotate` | first slice | `commands/annotate.rs` — `--rename-chrs` only |
| `call` | not started | dispatched to `unsupported` |
| `cnv` | not started | dispatched to `unsupported` |
| `concat` | broad slice | `commands/concat.rs` — `-a`/`-l` ligate remain |
| `consensus` | first slice | `commands/consensus.rs` — simple ALT application |
| `convert` | broad slice | `commands/convert.rs` — gVCF/HAP `-Ou` pipes blocked |
| `csq` | not started | dispatched to `unsupported`; `gff.rs` partial |
| `filter` | broad slice | `commands/filter.rs` |
| `gtcheck` | not started | dispatched to `unsupported` |
| `head` | complete enough | `commands/head.rs` + Perl `test_vcf_head`/`head2` enabled |
| `index` | complete enough | `commands/index.rs` + Perl `test_index`/`idxstats` enabled |
| `isec` | broad text slice | `commands/isec.rs` — full synced-reader pending |
| `merge` | first slice | `commands/merge.rs` — same-site only |
| `mpileup` | not started | dispatched to `unsupported` |
| `norm` | first slice | `commands/norm.rs` — `-d`/`--rm-dup` only |
| `plugin` | registry + 32 impls | `commands/plugin.rs` registry of 41 names; `commands/plugins/` implements `counts`, `missing2ref`, `fill-AN-AC`, `allele-length`, `variant-distance`, `check-ploidy`, `tag2tag`, `add-variantkey`, `variantkey-hex`, `remove-overlaps`, `af-dist`, `smpl-stats`, `indel-stats`, `ad-bias`, `prune`, `dosage`, `guess-ploidy`, `contrast`, `fixref`, `trio-switch-rate`, `trio-stats`, `mendelian2`, `parental-origin`, `fixploidy`, `GTsubset`, `GTisec`, `fill-from-fasta`, `scatter`, `split`, `isecGT`, `frameshifts`, `check-sparsity` |
| `query` | broad slice | `commands/query.rs` |
| `reheader` | broad slice | `commands/reheader.rs` |
| `roh` | not started | dispatched to `unsupported`; HMM kernel ready |
| `som` | out of scope | dispatched to `unsupported` (deferred) |
| `sort` | VNtyper-ready | `commands/sort.rs` + Perl `test_vcf_sort` options |
| `stats` | broad slice | `commands/stats.rs` |
| `tabix` | complete enough | `commands/tabix.rs` + Perl `test_tabix` enabled |
| `view` | broad slice | `commands/view.rs` — 64-bit BCF pipe parity pending |
| `bgzip` (helper) | Perl harness | `commands/bgzip.rs` — staged bgzip/tabix for `test.pl` |

32 of 41 plugin record-processing implementations done (see Wave F);
13 remain.

Current whole-project estimate:

- 2026-05-17 (post FORMAT subscript filter slice, PR #73 landed; no open PR):
  approximately 50-53% complete toward the full stated goal. Movement since
  the prior estimate includes `+split`, `+isecGT`, `+frameshifts`,
  `+check-sparsity`, and three shared filter-engine / `filter` slices covering
  FORMAT vector predicates, GT class literals, and FORMAT sample/value
  subscripts. 32 of 41 plugins done.
- 2026-05-16 (post `+trio-switch-rate`, PR #58 landed): approximately
  38-41% complete toward the full stated goal. Movement since the prior
  estimate is `+trio-switch-rate` (PED-trio phase-switch rate) verified
  byte-for-byte against `trio.out`, plus a reusable PED parser. 20 of 41
  plugins done.
- 2026-05-16 (post `+fixref`, PR #57 landed): approximately
  37-40% complete toward the full stated goal. Movement since the prior
  estimate is `+fixref` (FASTA REF/ALT strand fixing:
  ref-alt/swap/flip/flip-all) verified byte-for-byte against
  `fixref.{4,5,6,7}.out`. 19 of 41 plugins done.
- 2026-05-16 (post `+contrast`, PR #56 landed): approximately
  36-39% complete toward the full stated goal. Movement since the prior
  estimate is `+contrast` (control/case association: PASSOC Fisher exact,
  FASSOC, NASSOC, NOVELAL/NOVELGT) verified byte-for-byte against
  `contrast.out`/`.1.out`/`.1.1.out`/`.1.2.out`. 18 of 41 plugins done.
- 2026-05-16 (post `+prune -a/-m` LD modes, PR #55 landed): approximately
  35-38% complete toward the full stated goal of a pure Rust
  bcftools replacement with full subcommand, plugin, upstream `test.pl`,
  Rust integration-test, and parity-polishing coverage. Movement since the
  prior estimate is the full `vcfbuf` `calc_ld` (r2 / Lewontin's D' /
  Ragsdale's hd) + the HTSlib `kputd` float formatter ported, completing
  `+prune` (`-a`/`-m`/`-f`) byte-for-byte against `prune.1.1/1.2/1.3.out`
  and `prune.2.1.out`, on top of `guess-ploidy`, `dosage`, `ad-bias`,
  `indel-stats`, `smpl-stats`, `af-dist`, the `vcfbuf` overlap/dup state
  machine, the VariantKey pair, the PR #45 7-plugin batch and the
  PRs #10-#41 command slices. The raw checklist is well past two-thirds
  checked, but the estimate still weights the unfinished large
  subcommands (`mpileup`, `call`, `csq`, full `merge`/`annotate`/`norm`),
  the 24 remaining plugins (most coupled to the
  filter-engine/FASTA/PED infra still in progress), full upstream
  byte-for-byte parity, exit-code parity, and performance triage more
  heavily than scaffolding. The narrower BioScript VNtyper-useful local
  parity slice is roughly 75%+.
- 2026-05-15 (post `+guess-ploidy`): approximately 34-37% (kept for trend).
- 2026-05-15 (post `+dosage`): approximately 33-36% (kept for trend).
- 2026-05-15 (post `+prune`): approximately 32-35% (kept for trend).
- 2026-05-15 (post `+ad-bias`): approximately 31-34% (kept for trend).
- 2026-05-15 (post `+indel-stats`): approximately 30-33% (kept for trend).
- 2026-05-15 (post `+smpl-stats`): approximately 29-32% (kept for trend).
- 2026-05-15 (post `+af-dist`): approximately 28-31% (kept for trend).
- 2026-05-15 (post `+remove-overlaps`): approximately 27-30% (kept for trend).
- 2026-05-15 (post VariantKey pair): approximately 26-29% (kept for trend).
- 2026-05-15 (post 7-plugin batch): approximately 25-28% (kept for trend).
- 2026-05-15 (prior): approximately 22-25% (kept for trend).
- 2026-05-14: approximately 20% complete (prior estimate, kept for trend).

## Current Inputs

- `bcftools/`: upstream C bcftools source and test suite. 59 C files (~60k LOC) plus 41 plugin `.c` files under `plugins/`. ~28 built-in subcommands dispatched from `main.c:73-201`. A 2406-line Perl test harness (`bcftools/test/test.pl`) with ~1098 `run_test()` invocations and expected-output fixtures under `bcftools/test/<subcommand-or-plugin>/`.
- `htslib-rs/`: sibling pure-Rust HTSlib compatibility workspace. Re-exports `noodles` and ships HTSlib-shaped adapters under `crates/htslib-rs/src/*_compat.rs`. VCF/BCF coverage already includes header IDs, typed FORMAT/INFO adapters, allele removal, variant classification, sweep, synced-reader pairing, region index, and HTSlib expression evaluation.

## Pinned Scope Decisions

The following are decided up front and shape every phase below:

- **Subcommands**: target full parity with all upstream subcommands except those explicitly deferred (see *Out of Scope*).
- **Plugins**: all 41 upstream plugins are in scope and ported as in-process subcommands. There is no `dlopen`. Plugins are invoked via `bcftools +<name>` (alias of `bcftools plugin <name>`) exactly like upstream.
- **Layout**: workspace mirroring `htslib-rs`:
  - `crates/bcftools-rs` — library, one module per built-in subcommand under `src/commands/`, one module per plugin under `src/commands/plugins/`, shared infra under `src/`
  - `crates/bcftools-rs-cli` — the `bcftools` binary (dispatch + main, including the `+name` plugin-name shortcut)
- **HTSlib/noodles API gaps**: when bcftools-rs needs an HTSlib-shaped API that `htslib-rs` does not yet expose, or a `noodles` change is required, do not change those underlying libraries during the current bcftools-rs pass. Move the item to the end-of-file dependency-blocker list, keep working on other bcftools-rs tasks that do not need underlying-library changes, and stop once the remaining work is blocked on those dependencies. Do not bypass `htslib-rs` for HTSlib-shaped APIs. (Direct `noodles` use from bcftools-rs is acceptable only for code that has no HTSlib analogue.)
- **Tests — two gates**:
  1. **Parity gate**: upstream `bcftools/test/test.pl` run against the Rust binary. Expected outputs are the checked-in files under `bcftools/test/`. Used as a regression gate in CI.
  2. **Rust unit/integration tests**: per-subcommand `tests/` under `crates/bcftools-rs` using `cargo test`. Used for fine-grained development feedback and Rust-native edge cases.
- **Parity level**:
  - **Strict (byte-for-byte)**: VCF/BCF binary outputs, TBI/CSI index bytes, sort order, FASTA output from `consensus`, TSV/text outputs from `stats`, `query`, `gtcheck`, `roh`, `cnv`, `csq`, `isec`, exit codes.
  - **Semantic**: `##bcftools_<command>_Version`/`##bcftools_<command>_Command` header lines (match upstream's `ID:bcftools VN:<version> CL:<...>` shape — see *Header-versioning strategy* below), stderr error messages (same key information, wording may differ), usage/help text. `--no-version` must suppress the header line exactly like upstream (heavily used in `test.pl`).
- **C oracle**: local dev only. Devs MAY build upstream `bcftools` (in `bcftools/`) and use it to regenerate expected fixtures. CI does NOT build or run C bcftools — it only diffs against the checked-in fixtures.
- **Binary name**: `bcftools`. `test.pl` reads the binary path from `$opts->{bin}` (set via the harness `-b` flag); we pass that to point at our Rust build.
- **Expression engine**: bcftools ships its own filter-expression compiler in `filter.c` (171k, ~4500 LOC, 41+ operator tokens, sample-vector evaluation, lazy AC/AN/genotype caching, external value injection via `filter_test_ext`). This is **distinct** from `htslib-rs::expr` (which mirrors HTSlib's lighter `hts_expr.c`). The bcftools filter engine is ported as its own module in `bcftools-rs` — not by extending `htslib-rs::expr`.

## Out of Scope (deferred)

- `bcftools som` — flagged "experimental, do not advertise" in `main.c:194`. Port last or skip; the `test.pl` harness does not test it heavily.
- `bcftools polysomy` — only built under `USE_GPL` (links GNU Scientific Library). Defer or replace its GSL calls with `statrs`/native code; track separately.
- Remote I/O backends: `https://`, `s3://`, `ftp://`, `gs://`. Local-file paths only. (`htslib-rs` also defers these.)
- The C plugin ABI (`dlopen`, `init`/`process`/`destroy`, `BCFTOOLS_PLUGINS` env-var lookup of `.so` files). Plugins are in-process Rust modules; the env var is honored only for listing/help purposes.
- Windows-specific build paths (`_CRT_glob` in `main.c:260`, MSYS wildcard handling).
- Plugin makefiles and `plugins/*.mk` config — plugins are built into the binary.
- C ABI exposure. bcftools-rs is a Rust binary, not a library callable from C.

## Porting Principles

- Stay pure Rust. No `bindgen`, no `cc` crate, no linking to HTSlib C or to bcftools C.
- Default to `htslib-rs` for HTSlib-shaped helpers (header manipulation, INFO/FORMAT typed access, synced reader, region parsing, format detection, BGZF, index I/O). When `htslib-rs` lacks the API, or when a `noodles` change is needed, record that dependency at the end of this file, continue with other bcftools-rs-local work, and stop without making underlying-library changes once only dependency-blocked work remains.
- Preserve observable behavior under the parity rules above. Treat each `test.pl` test case as an acceptance test; do not mark a subcommand complete until both its `test.pl` cases and its Rust integration tests pass.
- Each subcommand is one module under `crates/bcftools-rs/src/commands/<name>.rs`, exposing `pub fn main(args: &[OsString]) -> ExitCode`. The CLI crate dispatches on `argv[1]` exactly like `main.c:298-306` and translates `+name` → `plugin name` exactly like `main.c:289-296`.
- Use `clap` for arg parsing but configure it to accept upstream's flag forms (short flags, long flags, value layout). Both `-Oz`/`-O z`/`--output-type z` must work — upstream uses `getopt_long` with attached short-arg values throughout.
- Errors: prefer `Result<T, E>` internally with a bcftools-rs error type; surface via `error` / `error_errno` equivalents that match upstream's "[E::funcname] message" stderr format (see `bcftools.h:54-60`).

## Header-versioning strategy

Upstream bcftools writes a per-command header line into output VCF/BCF, e.g. `##bcftools_viewVersion=1.21+htslib-1.21` and `##bcftools_viewCommand=view -Oz file.vcf.gz; Date=...` (see `bcf_hdr_append_version` in `bcftools.c`). To stay close to byte parity:

- Emit `##bcftools_<cmd>Version=<bcftools-upstream-version>+htslib-<htslib-rs-version>` and `##bcftools_<cmd>Command=<reconstructed argv>; Date=<RFC-2822 date>` where the VN matches the upstream bcftools version we are tracking (pin this in `version.rs`).
- `--no-version` must suppress both lines. `test.pl` passes `--no-version` to most invocations; that path **must** be exact.
- The `Date=` field uses HTSlib's `hts_time_with_tz` format; reproduce it bit-for-bit for tests that don't pass `--no-version`. Where the date makes a test non-deterministic, document it and either inject a fixed timestamp via env var or expect the test to set `--no-version`.

## Phase 0: Workspace Skeleton

- [x] Create root `Cargo.toml` workspace mirroring `htslib-rs/Cargo.toml`:
  - members: `crates/bcftools-rs`, `crates/bcftools-rs-cli`
  - workspace deps include `htslib-rs = { path = "../htslib-rs/crates/htslib-rs" }`
  - shared deps: `clap`, `anyhow`, `bstr`, `bytes`, `flate2`, `libdeflater`, `memchr`, `regex`, `noodles` (only if escape hatch needed). For HMM/stats: `statrs` (replaces GSL usage in `polysomy`, `cnv`, `roh`).
  - rust-version + edition matched to htslib-rs (`1.89.0`, `2024`)
- [x] Create the two crate skeletons with empty `lib.rs` / `main.rs` and a placeholder dispatcher.
- [x] Wire up `cargo fmt`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace` as the Rust gate.
- [x] Add a top-level GitHub Actions workflow with two jobs:
  - Rust gate (fmt + clippy + test).
  - Parity gate: build release binary, then run `cd bcftools/test && perl test.pl -b $WORKSPACE/target/release/bcftools` (or the harness's equivalent binary-override flag — confirm in Phase 4).
- [x] Document the project goal, scope, and CI gates in `README.md`.

## Phase 1: Shared Infrastructure

These are used by nearly every subcommand and must exist before subcommands can land.

- [x] **CLI dispatcher** (`bcftools-rs-cli/src/main.rs`): port `main.c:263-309`. Subcommand table grouped by section ("Indexing", "VCF/BCF manipulation", "VCF/BCF analysis", "Plugins"). `bcftools version` / `--version` / `--version-only` / `help [cmd]` exactly mirror `main.c:267-288`. `+name` plugin shortcut as in `main.c:289-296`. Unknown-subcommand error matches `main.c:307` wording.
- [x] **Version + feature string** (`bcftools-rs/src/version.rs`): export `BCFTOOLS_VERSION` constant (tracking upstream tag — see Submodule Pinning). Helper that prints both bcftools and htslib-rs versions.
- [x] **Common error helpers** (`bcftools-rs/src/diagnostics.rs`): port `error` / `error_errno` from `bcftools.h`/`bcftools.c` — `[E::funcname] message[: strerror(errno)]` then exit non-zero. Used uniformly.
- [x] **I/O & format helpers** (`bcftools-rs/src/io.rs`): port the `hts_bcf_wmode` / `hts_bcf_wmode2` / `set_wmode` family, `init_index` / `init_index2`, `init_tmp_prefix`, `write_index_parse`, `parse_overlap_option`, `apply_verbosity` (`bcftools.h:62-78`). Output-type dispatch (`-O v|z|u|b`) shared between every writer.
- [x] **Header-version writer** (`bcftools-rs/src/header_version.rs`): port `bcf_hdr_append_version`. Reconstructs argv with HTSlib-compatible quoting and produces `##bcftools_<cmd>Version` / `##bcftools_<cmd>Command` lines. `--no-version` short-circuits this.
- [ ] **Filter expression engine** (`bcftools-rs/src/filter/`): port `filter.c` (~4500 LOC) and `filter.h`. This is the bcftools `-i`/`-e` expression compiler/evaluator with sample-vector semantics, lazy AC/AN/genotype caching, `filter_test_ext` external-value injection, and `filter_max_unpack`/`filter_status` instrumentation. Used by `view`, `filter`, `query`, `isec`, `annotate`, `norm`, `stats`, `call`, `mpileup`, and many plugins. **This is the single largest porting task in the project.**
  - [x] Snapshot coverage: lexer for identifiers, INFO/FORMAT paths, numeric literals, quoted strings/escapes, comparison/regex/boolean/arithmetic operators including upstream single `&`/`|` aliases, function punctuation, vector index brackets, and FORMAT sample/value subscript colons; Pratt parser/AST for unary, binary, function-call, index, wildcard, and simple FORMAT subscript expressions; scalar evaluator for booleans, arithmetic, comparisons, regex matching, list indexing, FORMAT sample/value indexing (`TAG[:value]`, `TAG[sample:]`, `TAG[sample:value]`), simple list comparisons, top-level `FMT/`/`FORMAT/` sample-vector site predicates, `COUNT`/`MIN`/`MAX`/`SUM`/`AVG`/`MEAN`/`MEDIAN`/`STDEV`/`ABS`/`PHRED`/simple `binom`/simple `fisher`, plus `s*`/`SMPL_*` aliases for simple numeric aggregations, sample-context `N_PASS`/`F_PASS` over `FMT/`/`FORMAT/`/bare sample fields, no-argument `F_MISSING` missing-genotype fraction, limited GT special-literal equality for sample-count expressions (`GT="mis"` and related forms) plus bcftools GT class literals (`A`, `R`, `RR`, `RA`, `AR`, `AA`, `aA`, `Aa`, `HOM`, `HET`, `HAP`), external value injection callbacks for record/sample lookups, and evaluation tracing for lookup source/status plus short-circuit counts.
  - [ ] Remaining: full bcftools type system, exact regex/case-sensitivity parity, complete sample-vector semantics, lazy AC/AN/genotype caching, full `filter_max_unpack` parity, and integration into `view`/`query`/dependent commands.
- [x] **Synced reader wrapper** (`bcftools-rs/src/synced.rs`): bcftools-shaped facade over `htslib-rs::variant_io_compat::SyncedVariantGroup`/`pair_synced_variant_groups`. Exposes the `bcf_sr_t`-style API surface bcftools subcommands expect (add inputs, set regions/targets, iterate paired groups, `--collapse` modes). Where htslib-rs lacks a needed mode, extend it.
- [x] **Sample-list helpers** (`bcftools-rs/src/smpl_ilist.rs`): port `smpl_ilist.c` (sample subset, `^` exclusion, file-input form). Used by `view -s`, `call -s`, `stats -s`, many plugins.
- [x] **Region/target index** (`bcftools-rs/src/regidx.rs`): thin wrapper over `htslib-rs::regidx::RegionIndex` with bcftools-specific BED/region parsing helpers. Used by `view -R/-T`, `filter -R/-T`, `annotate`, `isec`, etc.
- [ ] **VCF buffer** (`bcftools-rs/src/vcfbuf.rs`): port `vcfbuf.c` (windowed buffer of `bcf1_t` records with overlap/window controls). Used by `+prune`, `+remove-overlaps`, `norm`, `+scatter`.
  - [x] Snapshot coverage: record-shape-independent window buffer with sorted insertion, half-open span overlap queries, contig-aware flushing, and configurable look-ahead window flushing.
  - [x] `MARK_OVERLAP`/`MARK_DUP` streaming state machine ported faithfully (FIFO + parallel mark buffer, `overlap_rid`/`overlap_end` span, `imin` left-aligned-indel prefix adjustment, `can_flush` drain) and wired into `+remove-overlaps` with full six-fixture byte parity.
  - [x] Windowed `_prune_sites` (`1st`/`maxAF`) and the `buf->win` site/bp window flush, wired into `+prune -n`/`-N`.
  - [x] `_calc_r2_ld` (genotype-dosage r2 / Lewontin's D' / Ragsdale's hd, `f64`) + the `vcfbuf_ld` window driver + the HTSlib `kputd` float formatter, wired into `+prune -a`/`-m`/`-f` with full `prune.1.1/1.2/1.3.out` + `prune.2.1.out` byte parity.
  - [ ] Remaining: the `cluster` mode (`-a count`/`-m count=`), wiring to concrete record-mutation paths in `norm`/`+scatter`, and the `MARK_EXPR` (`min(QUAL)`) / `prune -i,-e` paths (blocked on the filter engine).
- [ ] **`abuf` allele buffer** (`bcftools-rs/src/abuf.rs`): port `abuf.c` (allele-aware comparison buffer). Used by `norm`, `merge`, `+remove-overlaps`.
  - [x] Snapshot coverage: record-independent allele atomization for complex substitutions and indels, split-row construction, duplicate atom collapsing, source-ALT translation maps, and star-allele overlap marking.
  - [ ] Remaining: full allele comparison buffer behavior, concrete VCF/BCF record rewrite integration for `norm`/`merge`/plugins, INFO/FORMAT projection across atomized rows, and upstream `abuf.c` edge-case parity.
- [ ] **`convert` formatter** (`bcftools-rs/src/convert/`): port `convert.c` (76k). The `-f` format-string mini-language used by `query`, `convert`, and several plugins. Decisively non-trivial: token grammar, FORMAT/INFO tag expansion, sample iteration, GT special forms.
  - [x] Snapshot coverage: reusable format-string syntax parser for literals, common escapes, percent tokens, braced tokens, literal percent, sample loops, vector-index suffixes like `%TAG{1}`, and function-call tokens like `%SUM(TAG)`; record-agnostic renderer trait for scalar tokens, vector indexing, sample-loop expansion, case-insensitive numeric functions (`SUM`/`AVG`/`MEAN`/`MIN`/`MAX`/`ABS` plus `s*`/`SMPL_*` forms), simple `FORMAT/` sample-vector aggregation outside sample loops, and limited `%PBINOM(TAG)` support in sample loops using diploid text `GT` plus integer allele-indexed values.
  - [ ] Remaining: evaluation against real VCF/BCF records, full INFO/FORMAT expansion, complete GT special forms, exact sample iteration semantics, typed VCF/BCF `%PBINOM` parity, and wiring `query`/`convert` to the shared parser.
- [x] **gVCF helpers** (`bcftools-rs/src/gvcf.rs`): port `gvcf.c`. Used by `call`, `convert --gvcf2vcf`, `+gvcfz`.
- [x] **Reference helpers** (`bcftools-rs/src/reference.rs`): FASTA + FAI handling shared by `csq`, `consensus`, `mpileup`, `norm -c`, `+fill-from-fasta`. Routes through `htslib-rs::faidx_compat`.
- [ ] **GFF parser** (`bcftools-rs/src/gff.rs`): port `gff.c` (45k). Used only by `csq` but large and self-contained.
  - [x] Snapshot coverage: GFF3/GTF line parser with 1-based to 0-based half-open coordinate normalization, strand/phase validation, GFF3 comma attributes, GTF quoted attributes, and percent decoding; gene/transcript model grouping for GFF3 `ID`/`Parent` and GTF `gene_id`/`transcript_id`, including exon/CDS ordering by transcript strand and unplaced-feature tracking.
  - [ ] Remaining: full upstream transcript validation, CDS phase/frame reconciliation, SO consequence-specific structures, FASTA/reference integration, and full `csq` parity with upstream `gff.c`.
- [x] **Ploidy specification** (`bcftools-rs/src/ploidy.rs`): port `ploidy.c`. Used by `call`, `+fixploidy`, `+guess-ploidy`.
- [x] **HMM kernel** (`bcftools-rs/src/hmm.rs`): port `HMM.c`. Used by `roh`, `cnv`, `+parental-origin`.
- [x] **Math/numerics** (`bcftools-rs/src/numerics.rs`): port `kmin.c`, `peakfit.c`, `hclust.c`, `dist.c`, `em.c`, `prob1.c`. Used by `call`, `cnv`, `polysomy`, `gtcheck`, `+contrast`, `+af-dist`.
- [x] **TSV → VCF helper** (`bcftools-rs/src/tsv2vcf.rs`): port `tsv2vcf.c`. Used by `convert --tsv2vcf`, `+impute-info`, `+vrfs`.
- [x] **Logging passthrough**: bridge to `htslib-rs::log_compat` so `--verbosity` flows correctly to BGZF and synced-reader.

## Phase 2: Subcommand & Plugin Surface Mapping

Before implementing in waves, build `docs/subcommand-coverage.md` enumerating:

- [x] For each subcommand and plugin, list the upstream C source files it spans and every HTSlib/bcftools-internal API it calls (`bcf_sr_init`, `bcf_hdr_*`, `bcf_get_*`, `bcf_update_*`, `bcf_translate`, `filter_init`, `regidx_*`, etc.).
- [x] For each HTSlib API: column for `htslib-rs` coverage status (already exposed / needs to be added / out of scope).
- [x] For each bcftools-internal API: column for which Phase 1 module owns it.
- [x] Resulting gap list rolls up into `htslib-rs/TODO.md` extensions to be done before the dependent bcftools-rs subcommand can land.

## Phase 3: Subcommand Implementation Waves

Each subcommand below maps to: (a) one Rust module under `crates/bcftools-rs/src/commands/`, (b) the corresponding `test_*` cases in `bcftools/test/test.pl` passing against the Rust binary, (c) at least one Rust integration test under `crates/bcftools-rs/tests/<name>.rs`.

The waves are ordered to land foundational machinery first (read/write/index, the filter engine, the synced reader wrapper) and unblock the rest.

### Wave A — Read/Write/Index Foundation

- [ ] `view` (`vcfview.c`, 41k) — VCF↔BCF conversion, filtering (`-i`/`-e`/`-f`/`-G`/`-m`/`-M`/`-q`/`-Q`/`-v`/`-V`), sample/region restriction, `--no-version`. Anchor subcommand for parity testing. Covered by `test_vcf_view`. Depends on Phase 1 filter engine + synced reader wrapper.
  - [x] Snapshot coverage: VCF/VCF.gz/BCF read paths, VCF text/BGZF/BCF write paths including numeric `-O0`-`-O9` BGZF shorthand, stdin spooling, raw `--no-version` VCF passthrough, raw-header BCF VCF-text output, header-only/no-header modes, upstream-style option parsing when options appear after the input path, simple text VCF passthrough with textual version-header injection, HTSlib-style text normalization for Integer/Float INFO/FORMAT values outside BCF output, simple positional region filtering including `-r`/`-R` and braced contig names, text VCF region overlap modes via `--regions-overlap 0|1|2`, simple target filtering including `-t`/`-T` and `^` exclusion, text VCF target overlap modes via `--targets-overlap 0|1|2`, text VCF sample subsetting via `-s`/`-S` including BCF input to VCF output, BCF-output sample subsetting via the VCF projection path, text VCF `-G` genotype-column dropping, simple text VCF FILTER-list filtering via `-f`, limited text VCF expression filtering via `-i`/`-e` for core fields plus scalar and indexed INFO fields, simple text VCF type filtering via `-v`/`-V`, simple text VCF allele-count filtering via `-m`/`-M`, simple text VCF allele count/frequency filtering via `-c`/`-C`/`-q`/`-Q`, simple text VCF known/novel filtering via `-k`/`-n`, simple text VCF uncalled-site filtering via `-u`/`-U`, simple text VCF genotype-class filtering via `-g`, simple text VCF phased-site filtering via `-p`/`-P`, and threaded BGZF VCF/BCF writes. 52 integration tests in `crates/bcftools-rs/tests/view.rs`.
  - [ ] Remaining: full filter expression handling including FORMAT/sample expressions and advanced vector/sample slicing, complete FILTER/frequency/count/allele/type/genotype/phasing/known-novel/uncalled filter semantics across structured VCF/BCF writer paths, overlap-aware indexed region semantics, structured BCF/VCF writer overlap filtering, BCF output parity for 64-bit/out-of-range integer and missing INFO/FORMAT values, and full upstream `test_vcf_view` parity.
- [x] `head` (`vcfhead.c`) — header-only output, `-n N` line cap, `-s N` records-after-header cap. Covered by `test_vcf_head`, `test_vcf_head2`. Snapshot coverage: VCF/VCF.gz/BCF input paths, stdin handling for VCF/BCF, Kestrel-tolerant non-canonical VCF headers, BGZF VCF record-tail coverage, and dispatcher version/help/plugin shortcut behavior. 16 integration tests in `crates/bcftools-rs/tests/head.rs`.
- [x] `index` (`vcfindex.c`) — TBI/CSI build, `-s/--stats`, `-n/--nrecords`, `-c/--csi`, `--threads`. Covered by `test_index`, `test_vcf_idxstats`.
  - [x] Snapshot coverage: BCF CSI indexing, BGZF VCF CSI/TBI indexing,
        stdin indexing with explicit `-o`, custom output paths, overwrite
        protection, per-contig and total-record stats from data or index paths,
        large-coordinate CSI fixture queries through `view`, option validation,
        and rejection of extra input paths. 12 integration tests in
        `crates/bcftools-rs/tests/index.rs`.
- [x] `tabix` (`tabix.c`) — generic BGZF index/query for BED/GFF/SAM/VCF. Marked "do not advertise" upstream (`main.c:85`) but kept for tests. Covered by `test_tabix`.
  - [x] Snapshot coverage: BGZF VCF TBI build and query, CSI build and query,
        BED/GFF/SAM preset builds, `-a` streaming, existing-index refusal,
        attached `-pTYPE`/`-mINT` forms, and long aliases including
        `--preset`, `--force`, `--csi`, `--tbi`, and `--min-shift`.
        6 integration tests in `crates/bcftools-rs/tests/tabix.rs`.
- [ ] `query` (`vcfquery.c`, 20k) — `-f` format-string output, `-l/--list-samples`, region/target restriction, `--include`/`--exclude`. Depends on Phase 1 `convert` formatter + filter engine. Covered by `test_vcf_query`.
  - [x] Snapshot coverage: `-l/--list-samples` for VCF/VCF.gz/BCF, `-s`/`-S` sample selection including `^` exclusion, `-H`/`-HH` column headers for simple formats, `-u`/`--allow-undef-tags` compatibility for unknown format tokens rendered as `.`, POS-based `-r`/`-R`/`-t`/`-T` region and target restriction including braced contig names, limited record-level `-i`/`-e` filtering for core/INFO expressions including missing INFO values, comma-separated string element matches, `@file` string membership, semicolon-separated `ID` exact/regex/file membership, `strlen(...)`, indexed `AC`, and computed `AF`/`MAC`/`MAF`, simple ALT/INFO vector predicates (`ALT[*]~"..."`, `ALT="*"`, `TAG[*]="."`, `TAG[*]!="."`), `FILTER` ID and semicolon-set comparisons, simple FORMAT/sample predicates (`GT="."`, `GT="0|1"`, `GT="hom"`, `FMT/AD[:N]`, `FMT/AD[GT]`, `FMT/AD[:GT]`, `FMT/AD[0:GT]`, `sSUM(FMT/AD[GT])`, `binom(FMT/AD)`, `phred(binom(FMT/AD))`, `binom(FMT/AD[:N],FMT/AD[:N])`, simple `AD[:N]/sum(AD[*])`, and `FMT/`/`FORMAT/` tags), single-pipe sample masking vs double-pipe record OR for simple FORMAT predicates, simple sample-count filters (`N_PASS(...)`, `COUNT(...)`, `smpl_count(...)`), modulo comparisons, simple computed fields (`N_ALT`, `N_SAMPLES`, exact/regex/negated `TYPE`, `%ILEN`), core-field predicates (`CHROM`, `REF`, large `POS`), negative integer range predicates, and native bcftools-rs fallback expression evaluation for text records including simple `phred(binom(FMT/AD))`, `binom(INFO/TAG[N],INFO/TAG[N])`, and `phred(fisher(INFO/DP4))`; small text-backed `-f` formatter for core fields, implicit record newlines, `%LINE`, `%FORMAT`, INFO lookups, brace vector indexes (`%TAG{N}`), scientific-notation numeric output normalization, `%SAMPLE`, forced record namespace `%/TAG`, `%N_PASS(...)` sample counts, simple FORMAT/sample loops, `%TGT` allele-name genotype formatting, upstream-backed `query.func.1` numeric formatter fixtures, `%smpl_count(FMT/TAG)` sample-loop formatting, limited numeric functions (`SUM`/`AVG`/`MEAN`/`MIN`/`MAX`/`ABS`) over INFO and FORMAT values, and upstream-backed sample-loop `%PBINOM(TAG)`. 77 integration tests in `crates/bcftools-rs/tests/query.rs`.
  - [ ] Remaining: full `convert.c` formatter grammar, complete functions and GT special forms, indexed/overlap-aware region and target semantics, and full bcftools filter expression/sample-vector semantics.
- [ ] `stats` (`vcfstats.c`, 87k) — single-input and pairwise stats, depth/INFO/FORMAT histograms, sample-level stats, `-s` selection, `--af-bins`, `-i`/`-e`. The largest "report" subcommand. Covered by `test_vcf_stats`, `test_vcf_check`, `test_vcf_check_merge`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/stats.rs`): single-input `# SN`, `# TSTV`, fixed-row `# ST`, `# AF`, `# QUAL`, `# IDD`, `# HWE`, `# DP`, core `# PSC` sections with genotype-derived selected-sample singleton counts, selected-sample `# PSI` indel het/hom counters, genotype-derived `# SiS` singleton stats, selected-sample `# VAF` distributions from FORMAT/AD, basic `-E`/`--exons` indel frame-shift counts in `# FS` and `# PSI`, and basic `-F`/`--fasta-ref` indel-context sections (`# ICS`/`# ICL`); basic text-backed two-input pairwise A-only/B-only/shared reporting for the same sections; `-c`/`--collapse none|snps|indels|both|any|all|some` for pairwise text grouping including class-aware `both` separation of same-position SNP and indel records; `-f`/`--apply-filters`, `-i`/`-e` expression filtering via the shared filter engine for core fields, INFO tags, and computed `TYPE`; `-1`/`--1st-allele-only`, `-I`/`--split-by-ID`, `--af-bins`, `--af-tag`, `-u`/`--user-tstv TAG[:min:max:n]` including indexed INFO tags, `-d`/`--depth min,max,step` distribution from INFO/DP and FORMAT/DP or FORMAT/AD sample depths, `-s`/`-S` sample selection including `^` exclusion and `-s -` all-sample form, `-r`/`-R`/`-t`/`-T` POS-based region/target restriction. 25 integration tests in `crates/bcftools-rs/tests/stats.rs`.
  - [ ] Remaining: full indexed synced-reader pairwise parity, exact `--collapse` edge-case parity across multiallelic records, exact indel context, exon boundary, PSC/PSI AC/AN/frame-shift and edge-case parity.
- [ ] `isec` (`vcfisec.c`, 31k) — multi-input intersections, `-n`, `-w`, `-c`, `-C`, prefix output, `-p` directory output. Depends on synced reader. Covered by `test_vcf_isec`, `test_vcf_isec2`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/isec.rs`): text-backed VCF/VCF.gz/BCF set intersections and complements, bitmap stdout summaries to stdout or `-o FILE`, `-n [+-=]INT` and exact bitmask `-n~101` forms, `-C`, `-w LIST` VCF/BCF record output, `-c none|exact|any|all|some|both|snps|indels|id` including class-aware `snps`/`indels`/`both` collapse behavior for same-position mixed variant classes, single-input target-file VCF filtering, `-i`/`-e` record-level expression filtering through the shared filter engine, POS-based `-r`/`-R`/`-t`/`-T`, `-p DIR` directory output with `README.txt`, `sites.txt`, numbered VCF/VCF.gz/BCF files, automatic TBI/CSI indexing for `-p -O z|b` numbered outputs, and two-input default Venn layout (`0000`/`0001` private, `0002`/`0003` shared), Kestrel-tolerant text reads, `--no-version`, VCF.gz/BCF record output with `-O z|b`. 13 integration tests in `crates/bcftools-rs/tests/isec.rs`.
  - [ ] Remaining: full synced-reader multi-file iteration and overlap-aware indexed region/target semantics, exact upstream collapse-mode parity across multiallelic edge cases, structured header translation, and full upstream `test_vcf_isec*` parity.

### Wave B — File-Level Manipulation

- [ ] `norm` (`vcfnorm.c`, 116k) — left-align indels, split/join multiallelics (`-m -/+any/+snps/+indels/+both`), `-c` reference-check modes, `--rm-dup`, `--atomize`, `-N`. Depends on Phase 1 `abuf`, `vcfbuf`, reference. Covered by `test_vcf_norm`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/norm.rs`): first local command slice for duplicate-record removal with `-d`/`--rm-dup none|exact|snps|indels|both|any|all`, VCF/VCF.gz/BCF input, VCF/BGZF VCF/BCF output via `-O v|z|u|b`, `-o` file output, upstream-style PASS filter header insertion for normalized VCF text, and `--no-version` command-shape compatibility. 4 integration tests in `crates/bcftools-rs/tests/norm.rs`.
  - [ ] Remaining: left alignment, reference-check modes (`-c`), split/join multiallelics (`-m`), atomization, old-record tags, keep-sum INFO/FORMAT projection, overlap handling, right alignment with GFF, symbolic allele edge cases, and full upstream `test_vcf_norm` parity.
- [x] `sort` (`vcfsort.c`) — coordinate sort with disk-backed external-sort fallback (`extsort.c`). Covered by `test_vcf_sort`.
  - [x] **VNtyper compatibility target**: support the exact command shape used by upstream VNtyper's Kestrel post-processing:
        `bcftools sort <output_indel.vcf> -o <output_indel.vcf.gz> -W -O z`.
        This means coordinate-sorting VCF records, writing BGZF-compressed VCF
        output for `-O z`, and honoring `-W` by creating the matching VCF index.
        Full external-sort parity can come later, but this small-file path
        unblocks the BioScript VNtyper port.
  - [x] Snapshot coverage: coordinate/ref/ALT sorting, compressed VCF output,
        CSI/TBI indexing including attached `-Wtbi`, BGZF writer threading,
        Kestrel-tolerant header normalization, disk-backed temp-run spill, and
        upstream Perl `test_vcf_sort` option forms including `-m 0`,
        tiny `-m 1000` external-sort spills, attached `-Ob` BCF output, and
        BCF stdout piping through `view`. 11 integration tests in
        `crates/bcftools-rs/tests/sort.rs`.
- [ ] `concat` (`vcfconcat.c`, 52k) — vertical concat (`-a`, `-d`, `-l`, `--naive`, `--ligate`, `--regions`). Covered by `test_vcf_concat`, `test_naive_concat`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/concat.rs`): same-sample vertical concat for VCF/VCF.gz/BCF inputs, header preservation from first file, sample-column verification across inputs, default adjacent-input overlap rejection plus `-a`/`--allow-overlaps`, `-o`/`--output`, `-O u|b|v|z[0-9]`, `-f`/`--file-list`, `-G`/`--drop-genotypes`, `-D`/`--remove-duplicates`, `-d`/`--rm-dups snps|indels|both|all|exact`, `-n`/`--naive` VCF/VCF.gz body concatenation and `--naive-force`, upstream Perl `test_naive_concat` coverage for generated VCF.gz metadata-header differences and BCF inputs piped through `view`, `-r`/`-R` POS-based region restriction including BED coordinate conversion, `--regions-overlap 0|1|2` with record-span matching for 1/2, `-W`/`--write-index[=csi|tbi]` for VCF.gz/BCF outputs, `--threads` for VCF.gz/BCF file outputs, full `##bcftools_concat{Version,Command}` header line emission with `--no-version` suppression, Kestrel-tolerant text reads. 29 integration tests in `crates/bcftools-rs/tests/concat.rs`.
  - [ ] Remaining: full `-a`/`--allow-overlaps` edge-case parity with synced-reader overlap semantics, `-l`/`--ligate` and ligate-force/warn variants, true BCF block-level `--naive` concat for byte-level BCF output parity, `-c`/`--compact-PS`, `-q`/`--min-PQ`.
- [ ] `merge` (`vcfmerge.c`, 155k — largest single file in bcftools) — multi-sample merge across files, `-m none/snps/indels/both/all/id`, `--info-rules`, `-l`, `--regions`. Covered by `test_vcf_merge`, `test_vcf_merge_big`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/merge.rs`): same-site sample-column concatenation across VCF/VCF.gz/BCF inputs with identical fixed site fields, simple sorted text-mode union of exact sites that are absent from some inputs with upstream-style missing FORMAT sample values, sites-only same-locus ALT union with AC/AN summing and byte-for-byte `merge.8.out` text coverage, sampled plus sites-only same-locus ALT union with byte-for-byte `merge.9.1.out` text coverage, narrow `-i AN:sum,AC:sum` INFO-rule summing for that sampled/sites-only ALT-union slice with byte-for-byte `merge.9.2.out` text coverage, narrow `-i AF:join` INFO-rule joining for sites-only ALT-union records with byte-for-byte `merge.join.1.out` text coverage, repeated `--force-samples` duplicate-prefixing with byte-for-byte `merge.9.out` coverage, narrow `-m none` same-position conflict retention with byte-for-byte `merge.10.1.out` coverage, narrow sampled same-position allele unions for `-m both` and `-m snp-ins-del` with byte-for-byte `merge.10.{2,3}.out` coverage, sampled same-ref ALT-subset merge for `--merge none|both` with byte-for-byte `merge.multiallelics.1.1.out` coverage, sampled star-allele REF-extension merge with byte-for-byte `merge.5.out` coverage, sampled same-position INFO Number=A/R/G remapping and summing with byte-for-byte `merge.3.out` coverage, FORMAT-key union plus Number=A/R/G remapping with byte-for-byte `merge.6.out` coverage, default-mode sampled same-position allele union with FORMAT/AD vector expansion and byte-for-byte `merge.11.1.out` coverage, narrow `<NON_REF>` symbolic allele union for `--merge none|both` with byte-for-byte `merge.12.1.out` coverage, highest-input VCF fileformat preservation with byte-for-byte `merge.symbolic.1.1.out` coverage, byte-for-byte `merge.noidx.abc.out` text fixture coverage including `--no-index`, single-input `--force-single` command-shape compatibility with byte-for-byte `merge.LPL.0.out` coverage, upstream-style PASS header injection, header preservation from the first input, duplicate-sample-name rejection unless `--force-samples` prefixes later inputs, `-l`/`--file-list` input enumeration, `-o`/`--output`, `-O u|b|v|z`, other `-m TYPE` values accepted for command-shape compatibility, `--no-version` accepted, BGZF and BCF write paths, Kestrel-tolerant text reads. 23 integration tests in `crates/bcftools-rs/tests/merge.rs`.
  - [ ] Remaining: full synced-reader multi-input alignment beyond the simple sorted text union, overlapping sample-set merging with general allele unification, complete `-m none|snps|indels|both|all|id` semantics, `--info-rules`, `-l`/`--regions` index-driven streaming, `--missing-to-ref`, gVCF merging, and full upstream `test_vcf_merge*` parity.
- [ ] `reheader` (`reheader.c`, 27k) — header replacement, sample rename, FAI-driven contig fill, `--in-place` for BCF. Covered by `test_vcf_reheader`, `test_rename_chrs`.
  - [x] Snapshot coverage: VCF/BGZF VCF header replacement, sample rename via file/list including `--samples-list`, FAI contig updates with upstream-style attribute ordering, stdin handling, BCF output, BCF `--in-place`, and threaded BGZF/BCF output. 14 integration tests in `crates/bcftools-rs/tests/reheader.rs`.
  - [ ] Remaining: BCF header serialization order/quoting parity and `test_rename_chrs` dependencies on `annotate`/full `query`.
- [ ] `convert` (`vcfconvert.c`, 76k) — VCF ↔ {HAP/LEGEND/SAMPLE, GEN, HAPS-SAMPLE, TSV, gVCF, 23andMe}. Plus `--tsv2vcf`. Covered by `test_vcf_convert*`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/convert.rs`): explicit-column `--tsv2vcf` for TSV inputs with `CHROM,POS,ID,REF,ALT`, ignored `-` input columns, REF-matching ALT normalization to sites-only `.`, checked text-output parity for upstream `convert.tsv.vcf` and `convert.23andme.vcf` fixtures, upstream-style `--tsv2vcf -Ou | view` fixture pipes, optional `-s`/`--samples` and `-S`/`--samples-file` GT output from trailing VCF-style GT fields or allele-letter pairs, FASTA-backed `AA` columns for common 23andMe-style SNP genotypes with reference-derived REF/ALT alleles and contig headers including attached `-fFILE` reference parsing, skipped AA insertion/deletion rows, malformed data-row recovery with skip warnings and row/site counters, upstream-style row/site and genotype-class stderr counters, FASTA-backed `--gvcf2vcf` expansion of VCF/VCF.gz/BCF reference blocks with `INFO/END`, upstream mode-flag `--gvcf2vcf` argument shape plus legacy `--gvcf2vcf FILE`, BCF stdin input for gVCF conversion, upstream-style `--gvcf2vcf` filter-as-expansion-gate behavior where failing records are emitted unchanged, checked text-output parity for upstream `convert.gvcf.out`, basic VCF/BCF to `--gensample`, `--hapsample`, and `--haplegendsample` output with `.gen.gz`/`.hap.gz`/`.legend.gz`/`.samples`, upstream harness-style BCF stdin input for forward GEN/SAMPLE, HAP/SAMPLE, and HAP/LEGEND/SAMPLE output modes, stdout output for GEN/SAMPLE/HAP/LEGEND sinks named `-`, checked text-output parity for upstream `convert.gs.gt.gen`, `convert.gs.gt.ids.gen`, `convert.gs.gt.ids.gen6`, `convert.gs.gt.samples`, `convert.gs.pl.gen`, `convert.gs.pl.samples`, `check.gs.vcfids.gen`, `check.gs.vcfids.samples`, `check.gs.chrom.gen`, `check.gs.chrom.samples`, `check.gs.vcfids_chrom.gen`, `check.gs.vcfids_chrom.samples`, `convert.hs.hap`, `convert.hs.ids.hap`, `convert.hs.sample`, `convert.hls.haps`, `convert.hls.legend`, `convert.hls.ids.legend`, `convert.hls.samples`, and `convert.hap-missing.haps` fixtures, upstream-style single-precision PL/GL likelihood normalization and haploid missing HAP output (`? -`), upstream `-h` alias for HAP/LEGEND/SAMPLE output, `--gensample2vcf` back-conversion with GT/GP reconstruction, first-max GT tie handling, upstream-style GP number formatting, reversed marker/ID GEN columns, `--3N6`, VCF IDs, checked text-output parity for upstream `convert.gs.vcf` and `convert.gs.noids.vcf` fixtures, upstream-style reverse GEN/SAMPLE `-Ou | view` fixture pipe, VCF/VCF.gz/BCF output and indexing, basic `--hapsample2vcf` back-conversion with GT reconstruction, VCF IDs, checked text-output parity for upstream `convert.gt.noHead.vcf` and `convert.gt.noHead.ids.vcf` fixtures, VCF/VCF.gz/BCF output and indexing, basic `--haplegendsample2vcf` back-conversion with GT reconstruction, checked text-output parity for the upstream HAP/LEGEND/SAMPLE `convert.gt.noHead.vcf` fixture, VCF/VCF.gz/BCF output and indexing, GT/GP/PL/GL-backed GEN probability triples with clean unsupported-tag diagnostics, deprecated `--chrom` diagnostics, `--vcf-ids`, `--keep-duplicates`, `--haploid2diploid`, `--sex`, `-i`/`-e`, and `-s`/`-S` sample selection including `^` exclusion, VCF/VCF.gz/BCF output, `-o`/`-O u|b|v|z[0-9]`, `-W` indexing for VCF.gz/BCF outputs, `--threads`, and `--no-version`. 50 integration tests in `crates/bcftools-rs/tests/convert.rs`.
  - [ ] Remaining: VCF/BCF to TSV, 23andMe and full GEN/SAMPLE/HAPS/SAMPLE/HAP-LEGEND-SAMPLE edge-case parity; advanced gVCF2VCF filter-expression parity, full `--tsv2vcf` 23andMe edge-case parity and exact diagnostics; full upstream `test_vcf_convert*` parity.

### Wave C — Filtering & Annotation

- [ ] `filter` (`vcffilter.c`) — apply expression-based soft/hard filtering, set FILTER tags, `--mask`, `--SnpGap`, `--IndelGap`, `--set-GTs`. Heavily depends on Phase 1 filter engine. Covered by `test_vcf_filter`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/filter.rs`): VCF/VCF.gz/BCF read paths and VCF/VCF.gz/BCF write paths, `-i`/`-e` text-mode expression filtering via the shared filter engine for core fields and INFO tags plus simple FORMAT/sample contexts, upstream-backed sample-fraction functions `F_PASS(...)` and `F_MISSING` for text VCF filters, upstream single `&`/`|` FORMAT-vector site predicates with byte-for-byte `filter.2.out` coverage, bcftools GT class literals with byte-for-byte `filter.{12,13,14,15,16,17,18,19}.out` coverage, FORMAT sample/value subscripts with byte-for-byte `filter.{20,21,22,23,24,25}.out` coverage, missing QUAL comparison with byte-for-byte `filter.26.out` coverage, FORMAT/AO aggregate `-S .` behavior (`MIN`/`MAX`/`SUM`/`AVG`/`MEDIAN`/`STDEV` site-level functions vs `SMPL_*`/`s*` sample-level functions) with byte-for-byte `filter.{30,31,32,33,34,36}.out` text coverage, `-s`/`--soft-filter` re-tagging plus auto `##FILTER` header injection, joined short options used by upstream fixtures (`-sTAG`, `-S.`, `-m...`, `-g...`, `-G...`, `-iEXPR`, `-eEXPR`), `-m +` additive and `-m x` reset-pass modes, `--mask`/`-M` soft-filter masks including mask files and `^` negation, `--mask-overlap 0|1|2` POS/span matching, `-S`/`--set-GTs .|0` site-level failed-record genotype rewriting plus simple per-sample rewrites for FORMAT-scoped expressions, with existing INFO/AC and INFO/AN recalculation, `-g`/`--SnpGap` and `-G`/`--IndelGap` local text-mode gap filters including `--SnpGap` type-list matching for symbolic deletion spans with byte-for-byte `filter.29.out` text coverage and `--IndelGap` QUAL/AC/first-record tie-breaking, `-r`/`-R`/`-t`/`-T` POS-based region/target restriction, `-W`/`--write-index[=csi|tbi]` for VCF.gz/BCF outputs, `--threads` for VCF.gz/BCF file outputs, full `##bcftools_filter{Version,Command}` header line emission with `--no-version` suppression, Kestrel-tolerant text reads, shared `record_lookup` helper reused by `stats`. 37 integration tests in `crates/bcftools-rs/tests/filter.rs`.
  - [ ] Remaining: exact buffered gap-filter edge-case parity, full filter-expression FORMAT/sample-vector semantics, structured BCF write path that round-trips through the soft-filter rewrite without re-parsing.
- [ ] `annotate` (`vcfannotate.c`, 180k — single largest file in bcftools) — INFO/FORMAT/FILTER/ID column transfer from VCF/BCF/TAB sources, rename chrs, `-x` removal, header injection, `-c CHROM,POS,REF,ALT,…` column mapping, `--columns-file`, `--single-overlaps`, `--regions-overlap`. Covered by `test_vcf_annotate`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/annotate.rs`): `--rename-chrs` with two-column chromosome maps (contig-header `ID=` rewriting + CHROM-column rewriting), `-x`/`--remove` tag removal for `ID`, `QUAL`, `FILTER`, `FILTER/<ID>` (substituting `PASS` when the FILTER set empties), `INFO`, and `INFO/<ID>` (dropping the matching `##INFO`/`##FILTER` header lines), combined `--rename-chrs` + `-x`, VCF/VCF.gz/BCF input, VCF/BGZF VCF/BCF output via `-O v|z|u|b`, `-o` file output, and `--no-version` command-shape compatibility. 9 integration tests in `crates/bcftools-rs/tests/annotate.rs`.
  - [ ] Remaining: INFO/FORMAT/FILTER/ID transfer from VCF/BCF/TAB sources, `FORMAT/<ID>` removal, the `^`-keep-only `-x` form, header injection, `-c`/`--columns-file` mapping, sample-aware FORMAT annotation, overlap modes, merge logic, mark/missing modifiers, and full upstream `test_vcf_annotate` / `test_rename_chrs` parity.

### Wave D — Calling & Consequence

- [ ] `mpileup` (`mpileup.c`, 84k + `bam2bcf*.c`, `bam_sample.c`, `read_consensus.c`, `cigar_state.h`, `mw.h`) — multi-way pileup producing genotype likelihoods as BCF. Distinct from `samtools mpileup` (which produces text/VCF). Depends on `htslib-rs::alignment_compat` pileup iterators and `htslib-rs::probaln` for BAQ. Covered by `test_mpileup`.
- [ ] `call` (`vcfcall.c` + `mcall.c`, 65k + `ccall.c`, `prob1.c`, `em.c`) — multi-allelic (`-m`) and consensus (`-c`) callers, `--ploidy`, `--variants-only`, `--annotate FORMAT/PV4`, `-S` constrained allele set, gVCF mode. Covered by `test_vcf_call`, `test_vcf_call_cAls`.
- [ ] `consensus` (`consensus.c`, 55k) — apply VCF variants to FASTA reference, chain-file mode, `--missing`, `--mark-del/-ins/-snv`, `-H A/R/I/L` for haplotype selection, sample filters. Covered by `test_vcf_consensus`, `test_vcf_consensus_chain`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/consensus.rs`): first local command slice for FASTA consensus generation with plain VCF, VCF.gz/BGZF, and BCF inputs; required `-f`/`--fasta-ref`; simple REF-to-ALT SNP/indel application; FASTA record header/order preservation including region-style headers such as `1:2-501`; empty-VCF passthrough against the upstream `consensus.5.out` fixture; `-s` accepted for command-shape compatibility; `-H` numeric ALT-index selection; reference mismatch diagnostics. 4 integration tests in `crates/bcftools-rs/tests/consensus.rs`.
  - [ ] Remaining: sample-aware genotype/haplotype selection, IUPAC/marking modes, masks, absent/missing handling, chain output, overlap policy, symbolic allele edge cases, indexed region semantics, and full upstream `test_vcf_consensus*` parity.
- [ ] `csq` (`csq.c`, 166k) — variant consequence annotation given a GFF, supports phased calls, compound variants, splice consequences. Depends on Phase 1 `gff.rs` and reference helpers. Covered by `test_csq`, `test_csq_real`.

### Wave E — HMM / Stats / Trio

- [ ] `roh` (`vcfroh.c`, 52k) — HMM for runs of homozygosity, `--AF-dflt`, `--GTs-only`, `--estimate-AF`, viterbi/fwd-bwd. Depends on Phase 1 HMM. Covered by `test_roh`.
- [ ] `cnv` (`vcfcnv.c`, 60k) — HMM CNV calling from BAF/LRR. Depends on Phase 1 HMM and `peakfit`. Has no upstream `test.pl` coverage but the Rust gate must include synthetic integration tests.
- [ ] `gtcheck` (`vcfgtcheck.c`, 67k) — sample-concordance / contamination checks, `--pairwise`, `--dry-run`, `-e` for per-sample error rate. Covered by `test_gtcheck`.
- [ ] `polysomy` (`polysomy.c`, 34k) — chromosomal-copy detection. GPL-only upstream (uses GSL). Replace GSL with `statrs`/native; track as a separate milestone.
- [x] `som` (`vcfsom.c`, 25k) — experimental SOM-based filter. Defer (out of scope unless tests demand).

### Wave F — Plugins

All 41 plugins are in scope as in-process Rust implementations rather than
`dlopen`-loaded shared objects. They are invoked through `bcftools plugin
<name>` and `bcftools +<name>`. The `+name` dispatch lives in the CLI crate;
the `plugin` command's listing/help (`-l`, `-lv`, `-h`) walks a static plugin
registry rather than scanning `BCFTOOLS_PLUGINS` for `.so` files.

Implemented so far (PRs #45–#70): 32 plugins
under `crates/bcftools-rs/src/commands/plugins/` —
`counts`, `missing2ref`, `fill-AN-AC`, `allele-length`, `variant-distance`,
`check-ploidy`, `tag2tag` (gl-to-pl/gp-to-gt), `add-variantkey`,
`variantkey-hex`, `remove-overlaps`, `af-dist`, `smpl-stats`,
`indel-stats`, `ad-bias`, `prune`, `dosage`, `guess-ploidy`, `contrast`,
`fixref`, `trio-switch-rate`, `trio-stats`, `mendelian2`,
`parental-origin`, `fixploidy`, `GTsubset`, `GTisec`,
`fill-from-fasta`, `scatter`, `split`, `isecGT`, `frameshifts`,
`check-sparsity`. Every implemented plugin with an upstream `*.out`
fixture is byte-for-byte verified;
`variant-distance`/`check-ploidy` pass their entire `test_vcf_plugin`
slices, the two VariantKey plugins match the full
`query.add-variantkey.vcf` / `variantkey-hex.out` fixtures (66 records, 3
hash/non-reversible), `remove-overlaps` matches all six
`remove-overlaps.1.*` fixtures (overlap/dup/`-O t`/`--reverse`), `af-dist`
matches `af-dist.out` (HWE prob + AF-deviation histograms, `f32` binning),
`smpl-stats` matches `smpl-stats.1.out` (per-sample/per-site genotype
stats), `indel-stats` matches `indel-stats.1.out` (SN/DVAF/DLEN/DFRAC/
NFRAC), `ad-bias` matches `ad-bias.out` for both inputs (Fisher exact test
on FORMAT/AD), `prune` matches `prune.1.{1,2,3,4,6}.out` and
`prune.2.1.out` (windowed `_prune_sites` maxAF/1st **and** the
`calc_ld` r2/LD'/RD `-a`/`-m`/`-f` LD modes), `dosage` matches
`dosage.{1,2,3}.out` (PL/GL/GT likelihood/genotype dosages, `f32`), and
`guess-ploidy` matches `guess-ploidy.{PL,GL}.out` (PL/GL/GT
haploid/diploid log-likelihood sex inference, `f64`), and `contrast`
matches `contrast.out`/`.1.out`/`.1.1.out`/`.1.2.out` (control/case
PASSOC/FASSOC/NASSOC/NOVELAL/NOVELGT), `fixref` matches
`fixref.{4,5,6,7}.out` (FASTA REF/ALT ref-alt/swap/flip/flip-all), and
`trio-switch-rate` matches `trio.out` (PED-trio phase-switch rate +
per-population averages), and `trio-stats` matches `trio-stats.out`/
`trio-stats.2.out` (Mendelian/DNM/transmitted classification + debug
dump). The 9 remaining unimplemented plugins and many still-open plugin
subfeatures are heavier and coupled to shared infra still in progress: the
bcftools filter engine (`+setGT`, `+split-vep` expressions,
`remove-overlaps -m 'min(QUAL)'`, `smpl-stats`/`indel-stats`/`prune
-i/-e`), `hts_drand48` parity (`prune -N rand`), FASTA/reference
(`+fixref`, `+fill-from-fasta`), PED/trio handling (`+trio-stats`,
`+mendelian2`, `+trio-dnm3`, `indel-stats -p`), or VCF-rewrite/convert
(`ad-bias -c/-f`). The bcftools **filter engine** (`filter.c`, the
single largest task — unblocks the most plugins and subcommands) is the
preferred next pick.

Current local slice:

- [x] Static registry/listing surface in `crates/bcftools-rs/src/commands/plugin.rs`:
  `bcftools plugin -l`, `bcftools plugin -lv`, `bcftools +<name> --help`,
  and the `+name` shortcut cover all 41 upstream plugin names from
  `bcftools/plugins/*.c` without `dlopen` or `BCFTOOLS_PLUGINS` scanning.
  Per-plugin record-processing implementations live under
  `crates/bcftools-rs/src/commands/plugins/<name>.rs` and are dispatched from
  `plugin.rs` once ported.
- [x] `+counts` (`crates/bcftools-rs/src/commands/plugins/counts.rs`): counts
  samples/SNPs/indels/MNPs/others/sites with per-ALT classification routed
  through `htslib_rs::variant::classify_variant` (HTSlib `bcf_set_variant_type`
  port) OR-combined across ALTs like upstream `bcf_get_variant_types`; VCF/
  VCF.gz/BCF and stdin input via `bcftools +counts` / `bcftools plugin counts`;
  upstream-shaped six-line report. 5 integration tests in
  `crates/bcftools-rs/tests/plugin_counts.rs` + 3 unit tests. No upstream
  `test.pl` case exists for `+counts`.
- [x] `+fill-AN-AC` (`crates/bcftools-rs/src/commands/plugins/fill_an_ac.rs`):
  fills `INFO/AN` (total called alleles) and per-ALT `INFO/AC` from
  `FORMAT/GT`; existing AN/AC stripped first; `AC` omitted when there are no
  ALT alleles; `##INFO` lines for `AC` then `AN` inserted after the last
  existing `##INFO` line (HTSlib `bcf_hdr_append` grouping). VCF/VCF.gz/BCF +
  stdin input, `-o`/`-O u|b|v|z` output. Byte-for-byte parity with the
  upstream `plugin1.vcf` -> `fill-AN-AC.out` fixture. 4 integration tests +
  6 unit tests. Remaining: `+fill-tags` superset semantics.
- [x] `+tag2tag` (`crates/bcftools-rs/src/commands/plugins/tag2tag.rs`):
  exact integer conversions `--gl-to-pl` (`PL = lround(-10*GL)`, missing
  preserved) and `--gp-to-gt` (hard-call from normalized `GP`,
  `-t`/`--threshold`, call iff max posterior >= 1 - threshold, alleles via
  the HTSlib `bcf_gt2alleles` layout); `-r`/`--replace` drops the source
  FORMAT tag and its `##FORMAT` header and appends the dst header as the
  last `##` line. Byte-for-byte parity with upstream `view.GL.vcf`->`view.PL.vcf`
  and `view.GP.vcf`->`view.GT.vcf` (`test.pl` lines 676, 678). 4 integration
  tests + 4 unit tests. Remaining: float `--gl-to-gp` (`%g`) and the
  localized `--LXX-to-XX` family (`test.pl` 677, 679-681).
- [x] `+check-ploidy` (`crates/bcftools-rs/src/commands/plugins/check_ploidy.rs`):
  per-sample contiguous constant-ploidy regions
  (`Sample Chrom Start End Ploidy`); default ignores genotypes with any
  missing allele, `-m`/`--use-missing` counts missing slots; faithful
  upstream flush timing (chrom-change flush uses the previous chrom name,
  ploidy-change flush the current). Byte-for-byte parity with the upstream
  `checkploidy{,.2}.vcf` -> `checkploidy.{1,2,3}.out` fixtures covering
  `test.pl` lines 646-648. 4 integration tests + 4 unit tests.
- [x] `+variant-distance` (`crates/bcftools-rs/src/commands/plugins/variant_distance.rs`):
  annotates `INFO/<tag>` (default `DIST`) with distance to the nearest
  variant; `-d nearest|fwd|rev|both` (`both` is a Number=2 `<rev>,<fwd>` tag
  with `0` for a missing side), `-n`/`--tag-name`; same-POS records are
  duplicates sharing one distance; injects the implicit
  `##FILTER=<ID=PASS>` after `##fileformat` (HTSlib write behavior) and the
  `##INFO` tag line after the last `##INFO` (or before `#CHROM`). Byte-for-byte
  parity with all four upstream fixtures `variant-distance.{1,2,3,4}.out`
  covering `test.pl` lines 873-877. 7 integration tests + 5 unit tests.
- [x] `+allele-length` (`crates/bcftools-rs/src/commands/plugins/allele_length.rs`):
  REF / first-ALT / REF+ALT length histograms (MAXLEN=512, clamped) plus a
  non-base (`[^ACGTacgt]`) tally; first ALT only, matching upstream's
  `rec->d.allele[1]`. Text report. Byte-for-byte parity with the upstream
  `query.nucleotide.vcf` -> `query.allele-length.tsv` fixture. 2 integration
  tests + 3 unit tests.
- [x] `+missing2ref` (`crates/bcftools-rs/src/commands/plugins/missing2ref.rs`):
  default missing-genotype-to-ref behavior — every `.` allele token inside the
  `GT` FORMAT subfield becomes `0` while phase/unphase separators and all other
  FORMAT subfields are byte-preserved; GT located by FORMAT key index (not
  position). VCF/VCF.gz/BCF and stdin input; `-o`/`-O u|b|v|z` output via a
  shared `write_plugin_output` helper in `plugin.rs`. Byte-for-byte parity
  with the upstream `plugin1.vcf` -> `missing2ref.out` fixture
  (`test_vcf_plugin` / `+missing2ref --no-version`). 5 integration tests in
  `crates/bcftools-rs/tests/plugin_missing2ref.rs` + 5 unit tests. Remaining:
  `-e`/expression-gated and major-allele set modes.
- [x] `+add-variantkey` (`crates/bcftools-rs/src/commands/plugins/add_variantkey.rs`,
  shared algorithm in `plugins/variantkey.rs`): appends `VKX` (16-hex 64-bit
  VariantKey over CHROM, 0-based POS, REF, first ALT) and `RSX` (8-hex of the
  numeric `rs` ID) INFO fields, injecting the two `##INFO` lines immediately
  before `#CHROM` to match upstream `bcf_hdr_append` ordering after the
  harness `grep -v ^##bcftools_`. Faithful port of the MIT VariantKey
  reference (reversible ACGT encoding + MurmurHash3-like non-reversible
  path, exact `uint8_t`/`uint32_t` wrapping). VCF/VCF.gz/BCF and stdin input;
  `-o`/`-O u|b|v|z` via `write_plugin_output`. Byte-for-byte parity with
  `query.add-variantkey.vcf` (66 records, 3 hash/non-reversible).
  1 integration test in `crates/bcftools-rs/tests/plugin_add_variantkey.rs`
  + 5 unit tests.
- [x] `+variantkey-hex` (`crates/bcftools-rs/src/commands/plugins/variantkey_hex.rs`):
  suppresses VCF output and generates the three unsorted lookup files
  (`vkrs.unsorted.hex`, `rsvk.unsorted.hex`, `nrvk.unsorted.tsv`) under the
  optional output-directory positional (raw `strcat`-style concatenation,
  default `./`); `destroy()` summary (`VariantKeys:`/`Non-reversible
  VariantKeys:`) to stdout. Byte-for-byte parity with `variantkey-hex.out`
  plus regenerated lookup files. 1 integration test in
  `crates/bcftools-rs/tests/plugin_variantkey_hex.rs` + 1 unit test.
- [x] `+remove-overlaps` (`crates/bcftools-rs/src/commands/plugins/remove_overlaps.rs`):
  faithful port of the `vcfbuf` `MARK_OVERLAP`/`MARK_DUP` streaming state
  machine (FIFO record buffer + parallel mark buffer, `overlap_rid`/
  `overlap_end` running span, left-aligned-indel `imin` shared-prefix
  adjustment, `can_flush` drain) plus the `remove-overlaps.c` driver:
  `-m overlap`, `-m dup`, `-M TAG` (INFO flag injection with htslib-style
  `##FILTER=<ID=PASS>`/`##INFO` header normalization), `--reverse`, and
  `-O t` plain `chr<TAB>pos` site list. VCF/VCF.gz/BCF and stdin input;
  `-o`/`-O u|b|v|z` via `write_plugin_output`. Byte-for-byte parity with
  all six `remove-overlaps.1.{1..6}.out` fixtures. 6 integration tests in
  `crates/bcftools-rs/tests/plugin_remove_overlaps.rs` + 5 unit tests.
  Remaining: `-m 'min(QUAL)'` expression mode, `--missing`, and `-i`/`-e`
  filtering (all blocked on the bcftools filter engine port).
- [x] `+af-dist` (`crates/bcftools-rs/src/commands/plugins/af_dist.rs`):
  port of `af-dist.c` + the `bin.c` histogram (`bin_init`/`bin_get_idx`/
  `bin_get_value` for the `0..1` boundary case). Computes the HWE
  genotype-probability distribution (`2*AF*(1-AF)` for RA, `AF**2` for AA;
  `dosage==1`/`==2` only) and the AF-deviation distribution
  (`|AF - nALT/nALL|`), with all binning arithmetic in `f32` to match
  upstream's `float` edge sensitivity. Skips records with no INFO/AF and
  samples that are not fully called (vector_end/missing). VCF/VCF.gz/BCF
  and stdin input; `-t`/`--af-tag`, `-d`/`--dev-bins`, `-p`/`--prob-bins`.
  Byte-for-byte parity with `af-dist.out` after the harness
  `grep -v bcftools`. 1 integration test in
  `crates/bcftools-rs/tests/plugin_af_dist.rs` + 4 unit tests. Remaining:
  the `-l`/`--list` debug genotype dump.
- [x] `+smpl-stats` (`crates/bcftools-rs/src/commands/plugins/smpl_stats.rs`,
  default "all" filter): port of `smpl-stats.c` `process_record`/`destroy`.
  Per-sample stats (npass, non-ref, homRR, homAA, het, hemi, SNV, indel,
  singleton, missing, ts, tv, ts/tv) and the per-site rollup; `parse_genotype`
  hemizygous/vector_end semantics; `bcf_calc_ac` allele counts (INFO/AC+AN
  when present, else tallied from FORMAT/GT across all samples) for
  singleton detection; the upstream per-base `bcf_acgt2int` ts/tv walk;
  `classify_variant` for SNV/MNP-vs-indel typing; `ntv==0` → `inf` ts/tv.
  Emits the verbatim comment block + `CMD` line (harness strips `^CMD`).
  Byte-for-byte parity with `smpl-stats.1.out`. 1 integration test in
  `crates/bcftools-rs/tests/plugin_smpl_stats.rs` + 4 unit tests.
  Remaining: `-i`/`-e` filter-threshold scanning (curly-brace expansion +
  per-sample filter), blocked on the bcftools filter engine port.
- [x] `+indel-stats` (`crates/bcftools-rs/src/commands/plugins/indel_stats.rs`,
  no-PED default): port of `indel-stats.c` `process_record`/
  `update_indel_stats`/`destroy`. Record-level `VCF_INDEL` prefilter
  (`bcf_get_variant_types`), SN summary (nsites/npass/npass_gt/nins/ndel/
  nframeshift/ninframe), the FORMAT/AD variant-allele-frequency histogram
  (DVAF, `vaf2bin`), the indel-length histogram (DLEN, `len2bin`,
  het-of-two-indels both-allele recording), and the mean minor-allele
  fraction at HET indel genotypes vs length (DFRAC/NFRAC); the
  more-frequent-indel-allele selection from FORMAT/AD; CSQ
  `inframe`/`frameshift` substring detection; `var.n = len(ALT)-len(REF)`.
  Verbatim comment block + `CMD` line (harness strips `^CMD`). Byte-for-byte
  parity with `indel-stats.1.out`. 1 integration test in
  `crates/bcftools-rs/tests/plugin_indel_stats.rs` + 3 unit tests.
  Remaining: `-p` trio/de-novo mode, `-i`/`-e` filter scanning, and
  `--max-len`/`--nvaf` overrides (blocked on PED/filter infra).
- [x] `+ad-bias` (`crates/bcftools-rs/src/commands/plugins/ad_bias.rs`,
  report mode): port of `ad-bias.c`'s report path. Parses the
  sample/control pair file against the `#CHROM` order (skipping pairs not
  in the VCF), runs the upstream stateful two-most-frequent-allele search
  over FORMAT/AD (sample loop then control loop, with the carry-over
  `ibig`/`ismall`/`nbig`/`nsmall` state), applies `-d`/`-a` depth gates,
  and computes Fisher's exact test two-tail via
  `htslib_rs::math::kt_fisher_exact` (the HTSlib `kfunc.c` port). Emits
  `FT` lines below `-t` (default 1e-3) and the `SN` summary with C-style
  `%e` formatting (signed two-digit exponent). Byte-for-byte parity with
  `ad-bias.out` for both `ad-bias.vcf` and `ad-bias.2.vcf`. 2 integration
  tests in `crates/bcftools-rs/tests/plugin_ad_bias.rs` + 3 unit tests.
  Remaining: `-c`/`--clean-vcf` (VCF allele removal), `-v` variant-type
  filtering, and `-f` convert format (need convert/VCF-rewrite infra).
- [x] `+prune` (`crates/bcftools-rs/src/commands/plugins/prune.rs`): port
  of the `vcfbuf` windowed flush (`buf->win`, bp/site span) + `_prune_sites`
  (`1ST`/`MAX_AF`) for `-n`/`-N`, **and** the LD path: `_calc_r2_ld`
  (genotype-dosage r2 / Lewontin's D' / Ragsdale's hd, `f64`), the
  `vcfbuf_ld` window driver (per-metric max + position, `-m` early-exit),
  the HTSlib `kstring.c:kputd` float formatter, and the `prune.c`
  `-a`/`-m`/`-f` driver (`POS_*`/`R2`/`LD`/`RD` INFO + header injection,
  hard-drop vs soft-`FILTER`). The streaming push/flush dynamics, the
  `nbuf`/`nprune`/`eoff` removal arithmetic, the `cmpvrec` ordering, and
  the `kputd` formatter were traced against the fixtures before coding.
  htslib-style `##FILTER=<ID=PASS>` header injection. Byte-for-byte parity
  with `prune.1.1.out` (`-a r2,LD,HD`), `prune.1.2.out` (`-m 0.5 -f
  MaxR2`), `prune.1.3.out` (`-m 0.5`), `prune.1.4.out` (maxAF
  `--AF-tag`), `prune.1.6.out` (1st), and `prune.2.1.out` (20-sample).
  6 integration tests in `crates/bcftools-rs/tests/plugin_prune.rs` + 3
  unit tests. Remaining: `-a count`/`-m count=` cluster mode, `-N rand`
  (`hts_drand48` parity), and `-i`/`-e` filtering (filter engine).
- [x] `+dosage` (`crates/bcftools-rs/src/commands/plugins/dosage.rs`):
  port of `dosage.c`. `-t PL,GL,GT` ordered handlers (first applicable
  wins, header-gated for PL/GL); PL/GL dosages from diploid GL-ordered
  likelihoods (`10^(-0.1*PL)` / `10^GL`, normalized, accumulated per
  allele via the upstream `j/k/l` triangular loop) all in `f32` to match
  upstream `float` precision; GT alt-allele-count dosage; missing/short
  vector → `-1`; the `#[1]CHROM…[5]<sample>` header and per-record
  `CHROM/POS/REF/ALT` table with `%f` (PL/GL) and `%.1f` (GT) formatting.
  Byte-for-byte parity with `dosage.1.out` (`-t PL`), `dosage.2.out`
  (`-t GL`), and `dosage.3.out` (`-t GT`). 3 integration tests in
  `crates/bcftools-rs/tests/plugin_dosage.rs` + 4 unit tests.
- [x] `+guess-ploidy` (`crates/bcftools-rs/src/commands/plugins/guess_ploidy.rs`):
  port of `guess-ploidy.c`. `-r`/`-R` region restriction, SNP-only sites,
  the PL/GL/GT per-site genotype-probability derivation (PL via
  `pl2p[i]=10^(-i/10)` with the `<0||>=256→pl2p[255]` clamp, GL via
  `10^GL`, GT via the `-e` error model), per-site observed AF from the
  summed probabilities, then per-sample `log P(haploid)` /
  `log P(diploid)` accumulation (all `f64` to match upstream `double`);
  the diploid/all-haploid record split, the `vector_end`/missing/
  non-informative skips, the PL→GL→GT header auto-switch, and the
  verbose `SEX` report (`%f`, score computed at full precision). 
  Byte-for-byte parity with `guess-ploidy.PL.out` and
  `guess-ploidy.GL.out` (identical, exercising the PL→GL auto-switch).
  2 integration tests in `crates/bcftools-rs/tests/plugin_guess_ploidy.rs`
  + 2 unit tests. Remaining: `-g` genome shortcut begin-end sub-region,
  `--AF-tag`, and `-i`/`-e` filtering (filter engine).
- [x] `+contrast` (`crates/bcftools-rs/src/commands/plugins/contrast.rs`):
  port of `contrast.c`. `-0`/`-1` control/case sample groups (comma
  list or one-per-line file, sample-name precedence, `--force-samples`
  drops absent names); per-record `control_als`/`gt` allele bitmasks and
  `nals[4]` (ctrl-ref/ctrl-alt/case-ref/case-alt); `PASSOC`
  (`kt_fisher_exact` two-tail), `FASSOC` (`f32` non-REF proportions, `.`
  when undefined), `NASSOC` (4 ints), `NOVELAL` (case samples with an
  allele absent from controls), `NOVELGT` (novel genotype bitmask vs the
  control genotype set; `else if` after NOVELAL exactly as upstream); the
  requested `##INFO` defs + htslib `##FILTER=<ID=PASS>` header injection;
  every record written (skipped ones verbatim); floats via the shared
  `kputd`. Byte-for-byte parity with `contrast.out` (PASSOC,FASSOC,
  NOVELAL,NOVELGT; list **and** file `-0`/`-1`), `contrast.1.out`
  (NASSOC, `--force-samples` with an absent case sample), `contrast.1.1.out`
  (NOVELAL,NOVELGT) and `contrast.1.2.out` (NOVELGT). 5 integration tests
  in `crates/bcftools-rs/tests/plugin_contrast.rs` + 2 unit tests.
  Remaining: `-f` rare-allele enrichment (`max_AC`) and `-i`/`-e`
  filtering (filter engine).
- [x] `+fixref` (`crates/bcftools-rs/src/commands/plugins/fixref.rs`):
  port of `fixref.c` FASTA-reference strand fixing — `ref-alt` & `swap`
  (REF/ALT column changes only), `flip` & `flip-all` (also flip + swap
  genotypes 0<->1). `nt2int`/`int2nt`/`revint` complement, the
  single-base FASTA lookup via `crate::reference::FastaReference`
  (`faidx_compat`, builds the `.fai` on the fly), the exact `ir`/`ia`/`ib`
  decision tree per mode, the `FIXREF` dirty-bit INFO annotation
  (err/skip/none/flip/swap/GT order), non-SNP/non-biallelic/non-ACGT →
  `skip` (record written verbatim unless `-d`), whole-sequence suppress
  when the contig is absent from the FASTA, and `##INFO`/PASS header
  injection. Byte-for-byte parity with `fixref.4.out` (ref-alt),
  `fixref.5.out` (flip), `fixref.6.out` (flip-all), `fixref.7.out`
  (swap). 4 integration tests in `crates/bcftools-rs/tests/plugin_fixref.rs`
  + 3 unit tests. Remaining: `-m top` (Illumina TOP-strand sequence
  walking), `-i`/`-m id` (dbSNP second-VCF), and `-m stats` reporting
  (stats go to stderr, not exercised by `test.pl`).
- [x] `+trio-switch-rate` (`crates/bcftools-rs/src/commands/plugins/trio_switch_rate.rs`):
  port of `trio-switch-rate.c`. PED parser (`familyID sampleID
  paternalID maternalID sex phenotype [population]`, whitespace-split)
  resolving father/mother/child to header sample indices (rows with a
  parent absent from the VCF skipped); the phased-het-child phase-switch
  detection per trio — child must be a phased het, not both parents het,
  Mendelian-error when parent dosages tie, `test_phase` from the
  homozygous parent, `prev`/`nswitch` tracking with per-chromosome
  reset; the `TRIO` rows (`%.2f` switch rate) and the per-population
  averaged `POP` rows (`%.0f`/`%.2f`). Byte-for-byte parity with
  `trio.out`. 1 integration test in
  `crates/bcftools-rs/tests/plugin_trio_switch_rate.rs` + 2 unit tests.
  The PED parser here is the reusable basis for `+trio-stats`,
  `+mendelian2`, and `indel-stats -p`.
- [x] `+trio-stats` (`crates/bcftools-rs/src/commands/plugins/trio_stats.rs`,
  default "all" filter): port of `trio-stats.c` — the largest plugin so
  far. PED trios (dedup by `child father mother`), `bcf_calc_ac`
  (INFO/AC+AN else GT tally), per-trio `parse_genotype` (haploid → hom),
  `ac_trio`/star/non-ref handling, `bcf_acgt2int` per-site ts/tv, the
  Mendelian-error decision (`a0F`/`a1M`/`a0M`/`a1F`) with `ndnm_hom` and
  the `ndnm_recurrent` culprit selection via global `ac[culprit]`, the
  novel-singleton / untransmitted-singleton / transmitted-doubleton
  classification, and the `-a` max-alt-trios per-site cross-trio
  deferral (only counted when `nalt ≤ -a`). Interleaved
  `MERR`/`TRANSMITTED` debug lines (`-d mendel-errors,transmitted`)
  emitted during processing, then the verbatim 15-line comment header +
  `DEF`/`FLT0` summary (`%.2f` ts/tv, `inf` when `ntv==0`). Byte-for-byte
  parity with `trio-stats.out` (`-a 1`) and `trio-stats.2.out` (no `-a`).
  2 integration tests in `crates/bcftools-rs/tests/plugin_trio_stats.rs`
  + 2 unit tests. Remaining: `-i`/`-e` filter-threshold scanning and
  `-P` pfm trio (filter engine / single-trio mode).
- [x] `+mendelian2` (`crates/bcftools-rs/src/commands/plugins/mendelian2.rs`):
  port of `mendelian2.c`. Single `-p [1X:|2X:]P,F,M` trio; the built-in
  default `GRCh37` ruleset (`init_rules(args, NULL)` → alias `"GRCh37"`)
  is reproduced as a region table so each record resolves its
  `(sex_id, inherits, ploidy)` — chrX/Y/MT haploid (ploidy 1, M/F/.
  inheritance), all other regions MF/ploidy 2. `parse_gt` allele
  bitmasks; the haploid-kid branch (compare against the single
  inheriting parent, `ngood_alt` unless both ref), the diploid
  consistency branch (phase-consistent GOOD, parent-missing guards,
  else MERR), and the `c`/`a`/`d`/`e`/`g`/`m`/`E`/`M`/`S` modes +
  summary/per-trio count table. Byte-for-byte parity with
  `mendelian.{1,3,4,6,7,8}.out`. 6 integration tests in
  `crates/bcftools-rs/tests/plugin_mendelian2.rs` + 4 unit tests.
  Remaining: explicit `--rules`/`--rules-file` (other assemblies /
  custom ploidy) and `-i`/`-e` filtering (filter engine).
- [x] `+parental-origin` (`crates/bcftools-rs/src/commands/plugins/parental_origin.rs`):
  port of `parental-origin.c`. `-p P,F,M` trio, `-r REGION`,
  `-t del|dup`. Per-SNP FORMAT/PL→normalized GL, FORMAT/GT ALT dosage,
  FORMAT/AD; the DEL and DUP genotype-likelihood origin models (with
  the upstream observed-vs-deleted-allele accumulator swap for DEL),
  `-g` greedy ambiguous-site inclusion, and `-b` skewed-parental-het
  exclusion. Includes a local port of HTSlib `kfunc.c` `kf_betai`
  (modified Lentz continued fraction, reusing `htslib_rs::math::
  kf_lgamma`) backing `calc_binom_one_sided`/`calc_binom_two_sided`.
  Emits the `type/predicted_origin/quality/nmarkers` summary.
  Byte-for-byte parity with `parental-origin.{1,2,3,4,5}.out`. 5
  integration tests in `crates/bcftools-rs/tests/
  plugin_parental_origin.rs` + 3 unit tests. Remaining: `-i`/`-e`
  filtering and the `-d` informative-site debug listing (filter
  engine).
- [x] `+fixploidy` (`crates/bcftools-rs/src/commands/plugins/fixploidy.rs`):
  port of `fixploidy.c` plus the needed `ploidy.c` subset. Parses the
  `CHROM FROM TO SEX PLOIDY` regions file (`-p`) and the `NAME SEX`
  sex file (`-s`); supports `-f` forced ploidy, `-d` default ploidy,
  and the built-in default ruleset when neither `-p` nor `-f` is
  given. Mirrors `ploidy_query` (overlap vs default, non-default
  `_max`) and the bcf GT expand / in-place rewrite, preserving
  per-sample separators and phasing and producing exactly the
  per-sample ploidy (single `.` for ploidy 0). Byte-for-byte parity
  with `fixploidy.out`. 1 integration test in
  `crates/bcftools-rs/tests/plugin_fixploidy.rs` + 3 unit tests.
  Remaining: the `* * *` default-ploidy line form and `-t` tags
  other than GT (upstream also only supports GT).
- [x] `+GTsubset` (`crates/bcftools-rs/src/commands/plugins/gtsubset.rs`):
  port of `GTsubset.c`. Outputs only sites where every `-s` sample
  exclusively shares a genotype — all selected samples equal, no other
  sample equal — with missing GTs always passing. The genotype
  comparison is on the bcf-encoded allele/phase ints (so phasing is
  significant), `a1/a2` seeded from the first selected sample with
  both alleles set, vector-end padding for haploid samples, and the
  lone-`.`-sample → `.:.:…` expansion bcftools applies on write.
  Byte-for-byte parity with `view.GTsubset.{NA1,NA1NA2,NA1NA2NA3}.out`.
  3 integration tests in `crates/bcftools-rs/tests/plugin_gtsubset.rs`
  + 1 unit test.
- [x] `+GTisec` (`crates/bcftools-rs/src/commands/plugins/gtisec.rs`):
  port of `GTisec.c`. Counts genotype intersections across all
  possible sample subsets, emitting subset counts in banker's-sequence
  order — the upstream `compute_bankers` / `choose` recursion is
  ported verbatim (bitwise-symmetry shortcut included). Per record,
  samples are grouped by their `bcf_alleles2gt` key (haploid/short
  genotypes use a vector-end second allele so they never collide with
  diploid keys); each distinct genotype increments its sample-set
  bucket. Supports `-m` (per-sample missing counts), `-v` (annotate
  with sample lists) and `-H` (human-readable, implies `-v`, repeats
  each subset per contained sample) plus the verbatim comment header.
  Byte-for-byte parity with `view.GTisec.{,H,Hm,Hmv,Hv,m,mv,v}.out`.
  8 integration tests in `crates/bcftools-rs/tests/plugin_gtisec.rs`
  + 2 unit tests.
- [x] `+fill-from-fasta` (`crates/bcftools-rs/src/commands/plugins/fill_from_fasta.rs`):
  port of `fill-from-fasta.c`. Fills the REF allele (`-c REF`) or an
  INFO string tag from a FASTA reference. Parses the FASTA locally
  (name = first token after `>`, no `.fai` needed), fetches
  `[pos, pos+ref_len)`, applies the upstream uppercasing
  (`c>96 ? c-32`) and `-N` non-ACGTN→N, and supports the `-h`
  header-line append. Byte-for-byte parity with `ref.out` (`-c REF`)
  and `aa.2.out` (`-c REF -N`). 2 integration tests in
  `crates/bcftools-rs/tests/plugin_fill_from_fasta.rs` + 1 unit test.
  Remaining: `-i`/`-e` filtering (the `aa.out` row needs
  `-i 'TYPE="snp"'`) — filter engine.
- [x] `+scatter` (`crates/bcftools-rs/src/commands/plugins/scatter.rs`):
  port of `scatter.c`. Splits a VCF into multiple VCFs by fixed-size
  chunks (`-n N`, sequential integer-named files) or a region list
  (`-s REGIONS` / `-S FILE`, each region's records routed to its named
  file, `-x EXTRA` collecting records that match no region, `-p`
  prefix), every output file getting the full input header. Region
  parsing mirrors `regidx_parse_reg_name` (chrom-only ⇒ whole
  contig). Byte-for-byte parity with `scatter.1.{1,2,3}.out` via the
  sort-dir / `cat` / `grep -v ^##` harness. 3 integration tests in
  `crates/bcftools-rs/tests/plugin_scatter.rs` + 2 unit tests.
  Remaining: `-i`/`-e` filtering (filter engine).
- [x] `+split` (`crates/bcftools-rs/src/commands/plugins/split.rs`):
  filter-free slice of `split.c`. Splits VCF text into per-sample output
  VCFs by default, `-S` samples file, or `-G` groups file, including
  upstream filename sanitization / suffix collision behavior and sample
  renaming. Also supports `-k`/`--keep-tags` INFO/FORMAT projection in
  the text path and `-Oz` BGZF VCF output. Byte-for-byte parity with
  `split.1.{1,2,3,7}.out` and `split.2.1.out` through the harness
  shape that sorts output files, runs `query -l`, then `view -H`.
  7 integration tests in
  `crates/bcftools-rs/tests/plugin_split.rs` + 3 unit tests.
  Remaining: `-i`/`-e` filtering (rows 882–884), region/target
  restriction, BCF output, `-W` indexing, and output threading.
- [x] `+isecGT` (`crates/bcftools-rs/src/commands/plugins/isecgt.rs`):
  local text-backed slice of `isecGT.c`. Compares two VCF/BCF inputs
  by matching `CHROM/POS/REF/ALT`, maps samples from file A to file B
  by name with strict sample-count checks, and sets non-identical
  genotypes in file A to missing while preserving other FORMAT fields.
  Supports VCF text output and `-Oz` BGZF VCF output through the shared
  plugin writer. Upstream has no dedicated `test.pl` fixture for
  `+isecGT`; covered by 2 synthetic integration tests in
  `crates/bcftools-rs/tests/plugin_isecgt.rs` + 2 unit tests.
  Remaining: true `bcf_sr_t` synced-reader pairing/collapse behavior,
  `-r/-R/-t/-T` restriction, BCF output parity edge cases, `-W`
  indexing, and output threading.
- [x] `+frameshifts` (`crates/bcftools-rs/src/commands/plugins/frameshifts.rs`):
  local slice of `frameshifts.c`. Requires `-- -e EXONS`, parses simple
  BED-like exon files, and annotates simple insertion/deletion alleles
  with `INFO/OOF` values (out-of-frame `1`, in-frame `0`, not-applicable
  `-1`) while preserving existing INFO fields. Supports VCF text output
  and `-Oz` BGZF VCF output through the shared plugin writer. Upstream
  has no dedicated `test.pl` fixture row for `+frameshifts`; covered by
  2 synthetic integration tests in
  `crates/bcftools-rs/tests/plugin_frameshifts.rs` + 2 unit tests.
  Remaining: exact `bcf_sr_regions_t` parser parity for compressed/indexed
  exon lists, complex/symbolic allele edge cases, BCF output parity edge
  cases, `-W` indexing, and output threading.
- [x] `+check-sparsity` (`crates/bcftools-rs/src/commands/plugins/check_sparsity.rs`):
  local text-backed slice of `check-sparsity.c`. Reports samples with fewer
  than `-n` non-missing genotype calls per contig by default, or per
  requested `-r` / `-R` region. VCF/VCF.gz/BCF input is normalized through
  the existing text-view path. Upstream has no dedicated `test.pl` fixture
  row for `+check-sparsity`; covered by 2 synthetic integration tests in
  `crates/bcftools-rs/tests/plugin_check_sparsity.rs` + 3 unit tests.
  Remaining: true indexed tabix/BCF iterator behavior, exact upstream
  region-list label/ordering parity, GT storage-type edge cases, and
  stdin-with-region diagnostics.

Grouped roughly by complexity / shared dependencies:

- [ ] **Tag fixers** — `+fill-AN-AC`, `+fill-tags` (45k — heaviest of this group), `+missing2ref`, `+tag2tag`, `+setGT`, `+add-variantkey`, `+variantkey-hex`, `+allele-length`, `+impute-info`, `+counts`, `+dosage`, `+frameshifts`, `+remove-overlaps`, `+fill-from-fasta`.
- [ ] **Reference fixers** — `+fixref`, `+fixploidy`.
- [ ] **Subset/split** — `+split` (30k), `+scatter`, `+GTsubset`, `+GTisec`, `+isecGT`.
- [ ] **Stats / reports** — `+smpl-stats`, `+indel-stats`, `+trio-stats`, `+variant-distance`, `+ad-bias`, `+af-dist`, `+check-ploidy`, `+check-sparsity`, `+vcf2table` (46k), `+vrfs` (38k).
- [ ] **VEP-aware** — `+split-vep` (74k — the heaviest plugin by far).
- [ ] **Trio / pedigree** — `+mendelian2` (37k), `+trio-dnm3` (105k — the single largest plugin; has its own `test/trio-dnm3/test.sh` fixture), `+trio-switch-rate`, `+parental-origin`.
- [ ] **Sample inference** — `+guess-ploidy`, `+contrast`.
- [ ] **Misc** — `+color-chrs` (curses-style colored output), `+gvcfz`, `+prune`.

For each plugin: at least one `test.pl` case (most are covered by `test_vcf_plugin`, with named cases like `test_plugin_vrfs`, `test_plugin_split`, `test_plugin_scatter`, `test_trio_dnm3`) plus a Rust integration test.

## Phase 4: Test Harness Integration

- [x] **Parity gate setup**: confirm `bcftools/test/test.pl` can be driven against our Rust binary. Identify the harness's binary-override mechanism (read the `$$opts{bin}` setup at the top of `test.pl`) and document the exact invocation in `README.md`.
- [x] **`##bcftools_*` handling**: where upstream expected outputs include `##bcftools_<cmd>{Version,Command,Date}` lines we cannot reproduce, prefer adding `--no-version` to the test invocation (already used pervasively in `test.pl`). For tests that intentionally exercise the version line, pin the date via env var.
- [x] **Run progressively**: as each subcommand or plugin lands in Phase 3, enable its `test_*` in CI. Disabled tests should be tracked in `docs/test-status.md` as `not-yet-ported` (NOT just commented out).
- [x] **Rust integration tests per subcommand**: every implemented subcommand and plugin has its own `crates/bcftools-rs/tests/<name>.rs` (e.g. `tests/view.rs`, `tests/plugin_variant_distance.rs`) covering happy paths, error paths, regions/`-i`/`-e`/`-R`/`-T`/`-s` variants, threaded writes where applicable, and upstream `*.out` fixture parity where a fixture exists. The authoritative per-suite counts live in each command/plugin snapshot bullet (they drift if duplicated here). This item stays done; extend alongside each new slice rather than re-opening it.
- [ ] **`trio-dnm3` harness**: `bcftools/test/trio-dnm3/test.sh` is shelled out from `test.pl`. Confirm it works against the Rust binary unchanged, or port it to a Rust integration test.
- [x] **`csq` and `mpileup` fixtures**: `bcftools/test/csq/` and `bcftools/test/mpileup/` have nested fixture directories. The Rust gate must locate them via the `--path` form `test.pl` already uses.

## Phase 5: Parity Polishing

- [ ] **Diff every `test_*` output byte-for-byte** against the C bcftools outputs on a known fixture corpus (locally, dev-only). For each diff: classify (real bug / acceptable cosmetic / `##bcftools_*` only) and either fix or document.
- [ ] **Threads**: verify `--threads N` propagates to BGZF worker pools in writers (`view`, `merge`, `norm`, etc.) and matches upstream's parallelism behavior.
- [ ] **Exit codes**: confirm exit code matches upstream for invalid inputs, missing files, truncated BGZF, malformed records, header mismatches in `merge`/`isec`, etc.
- [ ] **Performance triage**: measure each subcommand on a representative dataset vs C bcftools. Goal is "within 2x" initially; performance fixes come after parity. Focus areas: `merge`, `annotate`, `csq`, `mpileup`, the filter engine.
- [x] **Bench harness**: criterion or custom timing harness under `benches/` for `view`, `merge`, `norm`, `annotate`, the filter engine.

## Dependency Blockers: htslib-rs/noodles Extensions Needed (rolling list)

This end-of-file list is filled as the subcommand surface mapping uncovers gaps that require changes outside `bcftools-rs`. During the current goal, move any `htslib-rs`, `noodles`, or submodule dependency work here, continue with independent bcftools-rs tasks, and stop when the remaining work depends on one of these blockers. Do not edit, patch, commit, or push the underlying libraries from this bcftools-rs pass.

- [ ] **`synced_bcf_reader` parity** — `htslib-rs::variant_io_compat` exposes pairing logic and no-index summaries, but bcftools depends on the full `bcf_sr_t` API surface: streaming iteration across multiple inputs, region/target restriction, `--collapse` modes (`none/snps/indels/both/any/some/id`), per-reader allele translation. Audit and extend.
- [x] **`bcf_translate`** — header-translation table between merged header and per-input header. Used in `merge`, `concat`, plugins. Confirm htslib-rs covers it beyond the existing translation fixture.
- [x] **`bcf_update_*` mutation API** — INFO/FORMAT/FILTER/ID/QUAL/POS/alleles mutation primitives. Partially in `VcfRecordAdapter`; extend to cover all upstream call sites.
- [x] **Pileup iterator surface for `mpileup`** — bcftools `mpileup` exercises far more of the HTSlib pileup API than samtools, including multi-input synchronized pileup. Audit `htslib-rs::alignment_compat` and extend.
- [x] **BAQ and `probaln_glocal`** — exposed by `htslib-rs::probaln`; verify wiring for `mpileup` (called from `bam2bcf*.c`).
- [ ] **`hts_set_threads` for BGZF writers** — wire to BCF/VCF.gz writers used by `view`, `merge`, `norm`, `concat`.
- [ ] **Custom tabix text presets** — `bcftools tabix -s/-b/-e/-0/-S/-c`
  needs an htslib-rs `tabix_compat` equivalent of a custom text-index format
  with configurable sequence/start/end columns, coordinate base, skip count,
  and comment character. Current `TextFormat` only exposes fixed
  BED/GFF/SAM/VCF presets, so bcftools-rs accepts and discards these arguments
  only for preset paths.
- [x] **CSI index 64-bit coordinate support** — `large_chrom_csi_limit` test in `test.pl:39` asserts the 2^31-1 boundary. Confirm htslib-rs CSI handles it.
- [x] **`hts_expr` vs bcftools filter** — `htslib-rs::expr` is the HTSlib expression language. bcftools has its own. Decide whether to also expose helpers from htslib-rs that bcftools's filter engine can reuse (numeric helpers, token utilities) or keep them fully separate. Document the decision.
- [x] **Region-with-target arithmetic** — `htslib-rs::region` covers HTSlib's grammar; confirm `-r`/`-R`/`-t`/`-T` semantics including the difference between regions (index-driven) and targets (streaming-filter) match upstream.
- [ ] **BCF serialization of haploid missing `GT=.`** — reverse
  `convert --hapsample2vcf -Ou` and `convert -H -Ou` hit
  `[E::main_vcfconvert] invalid input parameter` on the upstream Oxford
  fixtures when a haploid missing genotype (`GT=.`) is serialized through the
  current text-VCF-to-BCF writer path. The upstream gVCF pipe
  `convert --gvcf2vcf ... -Ou | view` hits the same class of failure on the
  first expanded record (`GT=.`). Text VCF parity is correct, GEN/SAMPLE
  `-G -Ou | view` passes, and TSV2VCF `-Ou | view` passes; the remaining
  HAP/SAMPLE, HAP/LEGEND/SAMPLE, and gVCF upstream pipe fixtures need
  `htslib-rs`/writer support for this genotype shape.
- [ ] **BCF serialization of 64-bit/out-of-range typed VCF values** —
  upstream `test_vcf_64bit` now passes all simple text-output cases in
  `bcftools-rs`, but the `view input.vcf -Ou | view -H` cases still fail while
  encoding or re-reading BCF. Failures include out-of-range Integer INFO/FORMAT
  values that HTSlib maps to missing, `INFO` missing-value arrays currently
  reaching an unimplemented `noodles-bcf` encoder path, and invalid `END`
  position handling. This needs `htslib-rs`/writer support for HTSlib's integer
  sentinel and missing-vector BCF semantics.

## Submodule Pinning

- [x] Pin `bcftools/` to a specific upstream release tag once Phase 0 lands (record tag + commit in `README.md` and `version.rs`). The pinned VN is what we emit in `##bcftools_*Version`.
- [x] Pin `htslib-rs/` to a known-green commit when Phase 0 lands. Bump deliberately when new `htslib-rs` extensions are required.

## Repository Map (target end state)

- `crates/bcftools-rs/` — library with one module per subcommand and one per plugin plus shared infra.
- `crates/bcftools-rs-cli/` — the `bcftools` binary.
- `bcftools/` — upstream C source + tests, used as fixture and reference only.
- `htslib-rs/` — sibling Rust HTSlib compatibility workspace consumed via path dep.
- `docs/subcommand-coverage.md` — per-subcommand/plugin HTSlib API surface and `htslib-rs` coverage status.
- `docs/test-status.md` — per-test pass/skip/not-yet-ported status.
- `TODO.md` — this file.
- `README.md` — project overview, scope decisions, build/test instructions.

## Development Workflow

```sh
# clone with submodules
git clone --recurse-submodules <repo>

# rust gate
cargo fmt --all
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace

# parity gate (against checked-in expected outputs)
cargo build --release
cd bcftools/test && perl test.pl   # binary path override per Phase 4

# optional: refresh expected outputs from C bcftools (local dev only)
cd bcftools && autoreconf -i && ./configure && make
cd test && perl test.pl --redo-outputs   # confirm flag name in Phase 4
```
