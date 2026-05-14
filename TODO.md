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

Latest landed progress:

- 2026-05-14: PR #8 (`progress/todo-local-bcftools-parity`, merge commit
  `8742124`) landed the first broad local-only parity batch for `concat`,
  `convert`, `filter`, `isec`, and `stats`, plus the dispatcher exports,
  command integration tests, TSV-to-VCF ALT normalization, and the snapshot
  coverage notes below.
- Validation before merge: `cargo fmt --all --check`,
  `cargo clippy --workspace --all-targets -- -D warnings`, and
  `cargo test --workspace`.
- Next local-only queue: continue upstream fixture parity for `convert`, then
  tighten `concat`, `filter`, `stats`, and `isec` edge cases that do not
  require changes in `htslib-rs`, `noodles`, or their submodules.

Current in-flight local progress:

- 2026-05-14: `progress/convert-fixture-parity-2` adds another convert parity
  slice after PR #8: more upstream GEN/SAMPLE, HAP/SAMPLE, and
  HAP/LEGEND/SAMPLE fixture-output parity, the upstream `-h` alias for
  HAP/LEGEND/SAMPLE output, single-precision PL/GL probability normalization,
  haploid missing HAP output parity, and harness-style BCF stdin input for
  forward GEN/SAMPLE, HAP/SAMPLE, and HAP/LEGEND/SAMPLE output modes,
  upstream-style `--tsv2vcf -Ou | view` fixture pipes, and upstream-style
  reverse GEN/SAMPLE `-Ou | view` fixture coverage. Local validation has
  passed with `cargo fmt --all --check`,
  `cargo clippy --workspace --all-targets -- -D warnings`, and
  `cargo test --workspace`.

Current whole-project estimate:

- 2026-05-14: approximately 20% complete toward the full stated goal of a pure
  Rust bcftools replacement with full subcommand, plugin, upstream `test.pl`,
  Rust integration-test, and parity-polishing coverage. The raw checklist is
  roughly half checked, but setup/foundation items are smaller than the
  remaining implementation surface. The estimate weights the unfinished large
  subcommands (`annotate`, `merge`, `norm`, `mpileup`, `call`, `csq`), most
  plugins, full upstream byte-for-byte parity, exit-code parity, and performance
  triage more heavily than scaffolding. The narrower BioScript VNtyper-useful
  local parity slice is much further along, roughly 60-70%+.

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
  - [x] Snapshot coverage: lexer for identifiers, INFO/FORMAT paths, numeric literals, quoted strings/escapes, comparison/regex/boolean/arithmetic operators, function punctuation, and vector index brackets; Pratt parser/AST for unary, binary, function-call, index, and wildcard expressions; scalar evaluator for booleans, arithmetic, comparisons, regex matching, list indexing, simple list comparisons, `COUNT`/`MIN`/`MAX`/`SUM`/`AVG`/`MEAN`/`STDEV`/`ABS`/`PHRED`/simple `binom`/simple `fisher`, plus `s*`/`SMPL_*` aliases for simple numeric aggregations, sample-context `N_PASS` over `FMT/`/`FORMAT/`/bare sample fields, external value injection callbacks for record/sample lookups, and evaluation tracing for lookup source/status plus short-circuit counts.
  - [ ] Remaining: full bcftools type system, exact regex/case-sensitivity parity, complete sample-vector semantics, lazy AC/AN/genotype caching, full `filter_max_unpack` parity, and integration into `view`/`query`/dependent commands.
- [x] **Synced reader wrapper** (`bcftools-rs/src/synced.rs`): bcftools-shaped facade over `htslib-rs::variant_io_compat::SyncedVariantGroup`/`pair_synced_variant_groups`. Exposes the `bcf_sr_t`-style API surface bcftools subcommands expect (add inputs, set regions/targets, iterate paired groups, `--collapse` modes). Where htslib-rs lacks a needed mode, extend it.
- [x] **Sample-list helpers** (`bcftools-rs/src/smpl_ilist.rs`): port `smpl_ilist.c` (sample subset, `^` exclusion, file-input form). Used by `view -s`, `call -s`, `stats -s`, many plugins.
- [x] **Region/target index** (`bcftools-rs/src/regidx.rs`): thin wrapper over `htslib-rs::regidx::RegionIndex` with bcftools-specific BED/region parsing helpers. Used by `view -R/-T`, `filter -R/-T`, `annotate`, `isec`, etc.
- [ ] **VCF buffer** (`bcftools-rs/src/vcfbuf.rs`): port `vcfbuf.c` (windowed buffer of `bcf1_t` records with overlap/window controls). Used by `+prune`, `+remove-overlaps`, `norm`, `+scatter`.
  - [x] Snapshot coverage: record-shape-independent window buffer with sorted insertion, half-open span overlap queries, contig-aware flushing, and configurable look-ahead window flushing.
  - [ ] Remaining: wire to concrete VCF/BCF record mutation paths in `norm`/plugins and audit behavior against upstream `vcfbuf.c` edge cases.
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
  - [x] Snapshot coverage: VCF/VCF.gz/BCF read paths, VCF text/BGZF/BCF write paths, stdin spooling, raw `--no-version` VCF passthrough, raw-header BCF VCF-text output, header-only/no-header modes, simple positional region filtering including `-r`/`-R` and braced contig names, text VCF region overlap modes via `--regions-overlap 0|1|2`, simple target filtering including `-t`/`-T` and `^` exclusion, text VCF target overlap modes via `--targets-overlap 0|1|2`, text VCF sample subsetting via `-s`/`-S` including BCF input to VCF output, BCF-output sample subsetting via the VCF projection path, text VCF `-G` genotype-column dropping, simple text VCF FILTER-list filtering via `-f`, limited text VCF expression filtering via `-i`/`-e` for core fields plus scalar and indexed INFO fields, simple text VCF type filtering via `-v`/`-V`, simple text VCF allele-count filtering via `-m`/`-M`, simple text VCF allele count/frequency filtering via `-c`/`-C`/`-q`/`-Q`, simple text VCF known/novel filtering via `-k`/`-n`, simple text VCF uncalled-site filtering via `-u`/`-U`, simple text VCF genotype-class filtering via `-g`, simple text VCF phased-site filtering via `-p`/`-P`, and threaded BGZF VCF/BCF writes.
  - [ ] Remaining: full filter expression handling including FORMAT/sample expressions and advanced vector/sample slicing, complete FILTER/frequency/count/allele/type/genotype/phasing/known-novel/uncalled filter semantics across structured VCF/BCF writer paths, overlap-aware indexed region semantics, structured BCF/VCF writer overlap filtering, and full upstream `test_vcf_view` parity.
- [x] `head` (`vcfhead.c`) — header-only output, `-n N` line cap, `-s N` records-after-header cap. Covered by `test_vcf_head`, `test_vcf_head2`.
- [x] `index` (`vcfindex.c`) — TBI/CSI build, `-s/--stats`, `-n/--nrecords`, `-c/--csi`, `--threads`. Covered by `test_index`, `test_vcf_idxstats`.
- [x] `tabix` (`tabix.c`) — generic BGZF index/query for BED/GFF/SAM/VCF. Marked "do not advertise" upstream (`main.c:85`) but kept for tests. Covered by `test_tabix`.
- [ ] `query` (`vcfquery.c`, 20k) — `-f` format-string output, `-l/--list-samples`, region/target restriction, `--include`/`--exclude`. Depends on Phase 1 `convert` formatter + filter engine. Covered by `test_vcf_query`.
  - [x] Snapshot coverage: `-l/--list-samples` for VCF/VCF.gz/BCF, `-s`/`-S` sample selection including `^` exclusion, `-H`/`-HH` column headers for simple formats, `-u`/`--allow-undef-tags` compatibility for unknown format tokens rendered as `.`, POS-based `-r`/`-R`/`-t`/`-T` region and target restriction including braced contig names, limited record-level `-i`/`-e` filtering for core/INFO expressions including missing INFO values, comma-separated string element matches, `@file` string membership, semicolon-separated `ID` exact/regex/file membership, `strlen(...)`, indexed `AC`, and computed `AF`/`MAC`/`MAF`, simple ALT/INFO vector predicates (`ALT[*]~"..."`, `ALT="*"`, `TAG[*]="."`, `TAG[*]!="."`), `FILTER` ID and semicolon-set comparisons, simple FORMAT/sample predicates (`GT="."`, `GT="0|1"`, `GT="hom"`, `FMT/AD[:N]`, `FMT/AD[GT]`, `FMT/AD[0:GT]`, `sSUM(FMT/AD[GT])`, `binom(FMT/AD)`, `phred(binom(FMT/AD))`, `binom(FMT/AD[:N],FMT/AD[:N])`, simple `AD[:N]/sum(AD[*])`, and `FMT/`/`FORMAT/` tags), single-pipe sample masking vs double-pipe record OR for simple FORMAT predicates, simple sample-count filters (`N_PASS(...)`, `COUNT(...)`, `smpl_count(...)`), modulo comparisons, simple computed fields (`N_ALT`, `N_SAMPLES`, exact/regex/negated `TYPE`, `%ILEN`), core-field predicates (`CHROM`, `REF`, large `POS`), negative integer range predicates, and native bcftools-rs fallback expression evaluation for text records including simple `phred(binom(FMT/AD))`, `binom(INFO/TAG[N],INFO/TAG[N])`, and `phred(fisher(INFO/DP4))`; small text-backed `-f` formatter for core fields, implicit record newlines, `%LINE`, `%FORMAT`, INFO lookups, brace vector indexes (`%TAG{N}`), scientific-notation numeric output normalization, `%SAMPLE`, forced record namespace `%/TAG`, `%N_PASS(...)` sample counts, simple FORMAT/sample loops, upstream-backed `query.func.1` numeric formatter fixtures, `%smpl_count(FMT/TAG)` sample-loop formatting, limited numeric functions (`SUM`/`AVG`/`MEAN`/`MIN`/`MAX`/`ABS`) over INFO and FORMAT values, and upstream-backed sample-loop `%PBINOM(TAG)`. 74 integration tests in `crates/bcftools-rs/tests/query.rs`.
  - [ ] Remaining: full `convert.c` formatter grammar, complete functions and GT special forms, indexed/overlap-aware region and target semantics, and full bcftools filter expression/sample-vector semantics.
- [ ] `stats` (`vcfstats.c`, 87k) — single-input and pairwise stats, depth/INFO/FORMAT histograms, sample-level stats, `-s` selection, `--af-bins`, `-i`/`-e`. The largest "report" subcommand. Covered by `test_vcf_stats`, `test_vcf_check`, `test_vcf_check_merge`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/stats.rs`): single-input `# SN`, `# TSTV`, fixed-row `# ST`, `# AF`, `# QUAL`, `# IDD`, `# HWE`, `# DP`, core `# PSC` sections with genotype-derived selected-sample singleton counts, selected-sample `# PSI` indel het/hom counters, genotype-derived `# SiS` singleton stats, selected-sample `# VAF` distributions from FORMAT/AD, basic `-E`/`--exons` indel frame-shift counts in `# FS` and `# PSI`, and basic `-F`/`--fasta-ref` indel-context sections (`# ICS`/`# ICL`); basic text-backed two-input pairwise A-only/B-only/shared reporting for the same sections; `-c`/`--collapse none|snps|indels|both|any|all|some` for pairwise text grouping; `-f`/`--apply-filters`, `-i`/`-e` expression filtering via the shared filter engine, `-1`/`--1st-allele-only`, `-I`/`--split-by-ID`, `--af-bins`, `--af-tag`, `-u`/`--user-tstv TAG[:min:max:n]` including indexed INFO tags, `-d`/`--depth min,max,step` distribution from INFO/DP and FORMAT/DP or FORMAT/AD sample depths, `-s`/`-S` sample selection including `^` exclusion and `-s -` all-sample form, `-r`/`-R`/`-t`/`-T` POS-based region/target restriction. 22 integration tests in `crates/bcftools-rs/tests/stats.rs`.
  - [ ] Remaining: full indexed synced-reader pairwise parity, exact `--collapse` edge-case parity across multiallelic records, exact indel context, exon boundary, PSC/PSI AC/AN/frame-shift and edge-case parity.
- [ ] `isec` (`vcfisec.c`, 31k) — multi-input intersections, `-n`, `-w`, `-c`, `-C`, prefix output, `-p` directory output. Depends on synced reader. Covered by `test_vcf_isec`, `test_vcf_isec2`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/isec.rs`): text-backed VCF/VCF.gz/BCF set intersections and complements, bitmap stdout summaries to stdout or `-o FILE`, `-n [+-=]INT` and exact bitmask `-n~101` forms, `-C`, `-w LIST` VCF/BCF record output, `-c none|exact|any|all|some|both|snps|indels|id` for the common upstream fixture cases, single-input target-file VCF filtering, `-i`/`-e` record-level expression filtering through the shared filter engine, POS-based `-r`/`-R`/`-t`/`-T`, `-p DIR` directory output with `README.txt`, `sites.txt`, numbered VCF/VCF.gz/BCF files, automatic TBI/CSI indexing for `-p -O z|b` numbered outputs, and two-input default Venn layout (`0000`/`0001` private, `0002`/`0003` shared), Kestrel-tolerant text reads, `--no-version`, VCF.gz/BCF record output with `-O z|b`. 12 integration tests in `crates/bcftools-rs/tests/isec.rs`.
  - [ ] Remaining: full synced-reader multi-file iteration and overlap-aware indexed region/target semantics, exact upstream collapse-mode parity across multiallelic edge cases, structured header translation, and full upstream `test_vcf_isec*` parity.

### Wave B — File-Level Manipulation

- [ ] `norm` (`vcfnorm.c`, 116k) — left-align indels, split/join multiallelics (`-m -/+any/+snps/+indels/+both`), `-c` reference-check modes, `--rm-dup`, `--atomize`, `-N`. Depends on Phase 1 `abuf`, `vcfbuf`, reference. Covered by `test_vcf_norm`.
- [x] `sort` (`vcfsort.c`) — coordinate sort with disk-backed external-sort fallback (`extsort.c`). Covered by `test_vcf_sort`.
  - [x] **VNtyper compatibility target**: support the exact command shape used by upstream VNtyper's Kestrel post-processing:
        `bcftools sort <output_indel.vcf> -o <output_indel.vcf.gz> -W -O z`.
        This means coordinate-sorting VCF records, writing BGZF-compressed VCF
        output for `-O z`, and honoring `-W` by creating the matching VCF index.
        Full external-sort parity can come later, but this small-file path
        unblocks the BioScript VNtyper port.
- [ ] `concat` (`vcfconcat.c`, 52k) — vertical concat (`-a`, `-d`, `-l`, `--naive`, `--ligate`, `--regions`). Covered by `test_vcf_concat`, `test_naive_concat`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/concat.rs`): same-sample vertical concat for VCF/VCF.gz/BCF inputs, header preservation from first file, sample-column verification across inputs, `-o`/`--output`, `-O u|b|v|z[0-9]`, `-f`/`--file-list`, `-G`/`--drop-genotypes`, `-D`/`--remove-duplicates`, `-d`/`--rm-dups snps|indels|both|all|exact`, `-n`/`--naive` VCF body concatenation and `--naive-force`, `-r`/`-R` POS-based region restriction including BED coordinate conversion, `--regions-overlap 0|1|2` with record-span matching for 1/2, `-W`/`--write-index[=csi|tbi]` for VCF.gz/BCF outputs, `--threads` for VCF.gz/BCF file outputs, full `##bcftools_concat{Version,Command}` header line emission with `--no-version` suppression, Kestrel-tolerant text reads. 26 integration tests in `crates/bcftools-rs/tests/concat.rs`.
  - [ ] Remaining: `-a`/`--allow-overlaps` (depends on synced reader parity), `-l`/`--ligate` and ligate-force/warn variants, true BCF block-level `--naive` concat, `-c`/`--compact-PS`, `-q`/`--min-PQ`.
- [ ] `merge` (`vcfmerge.c`, 155k — largest single file in bcftools) — multi-sample merge across files, `-m none/snps/indels/both/all/id`, `--info-rules`, `-l`, `--regions`. Covered by `test_vcf_merge`, `test_vcf_merge_big`.
- [ ] `reheader` (`reheader.c`, 27k) — header replacement, sample rename, FAI-driven contig fill, `--in-place` for BCF. Covered by `test_vcf_reheader`, `test_rename_chrs`.
  - [x] Snapshot coverage: VCF/BGZF VCF header replacement, sample rename via file/list, FAI contig updates with upstream-style attribute ordering, stdin handling, BCF output, BCF `--in-place`, and threaded BGZF/BCF output.
  - [ ] Remaining: BCF header serialization order/quoting parity and `test_rename_chrs` dependencies on `annotate`/full `query`.
- [ ] `convert` (`vcfconvert.c`, 76k) — VCF ↔ {HAP/LEGEND/SAMPLE, GEN, HAPS-SAMPLE, TSV, gVCF, 23andMe}. Plus `--tsv2vcf`. Covered by `test_vcf_convert*`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/convert.rs`): explicit-column `--tsv2vcf` for TSV inputs with `CHROM,POS,ID,REF,ALT`, ignored `-` input columns, REF-matching ALT normalization to sites-only `.`, checked text-output parity for upstream `convert.tsv.vcf` and `convert.23andme.vcf` fixtures, upstream-style `--tsv2vcf -Ou | view` fixture pipes, optional `-s`/`--samples` and `-S`/`--samples-file` GT output from trailing VCF-style GT fields or allele-letter pairs, FASTA-backed `AA` columns for common 23andMe-style SNP genotypes with reference-derived REF/ALT alleles and contig headers, skipped AA insertion/deletion rows, malformed data-row recovery with skip warnings and row/site counters, upstream-style row/site and genotype-class stderr counters, FASTA-backed `--gvcf2vcf` expansion of VCF/VCF.gz/BCF reference blocks with `INFO/END`, upstream mode-flag `--gvcf2vcf` argument shape plus legacy `--gvcf2vcf FILE`, BCF stdin input for gVCF conversion, upstream-style `--gvcf2vcf` filter-as-expansion-gate behavior where failing records are emitted unchanged, checked text-output parity for upstream `convert.gvcf.out`, basic VCF/BCF to `--gensample`, `--hapsample`, and `--haplegendsample` output with `.gen.gz`/`.hap.gz`/`.legend.gz`/`.samples`, upstream harness-style BCF stdin input for forward GEN/SAMPLE, HAP/SAMPLE, and HAP/LEGEND/SAMPLE output modes, stdout output for GEN/SAMPLE/HAP/LEGEND sinks named `-`, checked text-output parity for upstream `convert.gs.gt.gen`, `convert.gs.gt.ids.gen`, `convert.gs.gt.ids.gen6`, `convert.gs.gt.samples`, `convert.gs.pl.gen`, `convert.gs.pl.samples`, `check.gs.vcfids.gen`, `check.gs.vcfids.samples`, `check.gs.chrom.gen`, `check.gs.chrom.samples`, `check.gs.vcfids_chrom.gen`, `check.gs.vcfids_chrom.samples`, `convert.hs.hap`, `convert.hs.ids.hap`, `convert.hs.sample`, `convert.hls.haps`, `convert.hls.legend`, `convert.hls.ids.legend`, `convert.hls.samples`, and `convert.hap-missing.haps` fixtures, upstream-style single-precision PL/GL likelihood normalization and haploid missing HAP output (`? -`), upstream `-h` alias for HAP/LEGEND/SAMPLE output, `--gensample2vcf` back-conversion with GT/GP reconstruction, first-max GT tie handling, upstream-style GP number formatting, reversed marker/ID GEN columns, `--3N6`, VCF IDs, checked text-output parity for upstream `convert.gs.vcf` and `convert.gs.noids.vcf` fixtures, upstream-style reverse GEN/SAMPLE `-Ou | view` fixture pipe, VCF/VCF.gz/BCF output and indexing, basic `--hapsample2vcf` back-conversion with GT reconstruction, VCF IDs, checked text-output parity for upstream `convert.gt.noHead.vcf` and `convert.gt.noHead.ids.vcf` fixtures, VCF/VCF.gz/BCF output and indexing, basic `--haplegendsample2vcf` back-conversion with GT reconstruction, checked text-output parity for the upstream HAP/LEGEND/SAMPLE `convert.gt.noHead.vcf` fixture, VCF/VCF.gz/BCF output and indexing, GT/GP/PL/GL-backed GEN probability triples, deprecated `--chrom` diagnostics, `--vcf-ids`, `--keep-duplicates`, `--haploid2diploid`, `--sex`, `-i`/`-e`, and `-s`/`-S` sample selection including `^` exclusion, VCF/VCF.gz/BCF output, `-o`/`-O u|b|v|z[0-9]`, `-W` indexing for VCF.gz/BCF outputs, `--threads`, and `--no-version`. 49 integration tests in `crates/bcftools-rs/tests/convert.rs`.
  - [ ] Remaining: VCF/BCF to TSV, 23andMe and full GEN/SAMPLE/HAPS/SAMPLE/HAP-LEGEND-SAMPLE edge-case parity; advanced gVCF2VCF filter-expression parity, full `--tsv2vcf` 23andMe edge-case parity and exact diagnostics; full upstream `test_vcf_convert*` parity.

### Wave C — Filtering & Annotation

- [ ] `filter` (`vcffilter.c`) — apply expression-based soft/hard filtering, set FILTER tags, `--mask`, `--SnpGap`, `--IndelGap`, `--set-GTs`. Heavily depends on Phase 1 filter engine. Covered by `test_vcf_filter`.
  - [x] Snapshot coverage (`crates/bcftools-rs/src/commands/filter.rs`): VCF/VCF.gz/BCF read paths and VCF/VCF.gz/BCF write paths, `-i`/`-e` text-mode expression filtering via the shared filter engine for core fields and INFO tags, `-s`/`--soft-filter` re-tagging plus auto `##FILTER` header injection, `-m +` additive and `-m x` reset-pass modes, `--mask`/`-M` soft-filter masks with `^` negation and `--mask-overlap 0|1|2` POS/span matching, `-S`/`--set-GTs .|0` site-level failed-record genotype rewriting plus simple per-sample rewrites for FORMAT-scoped expressions, with existing INFO/AC and INFO/AN recalculation, `-g`/`--SnpGap` and `-G`/`--IndelGap` local text-mode gap filters including `--IndelGap` QUAL/AC/first-record tie-breaking, `-r`/`-R`/`-t`/`-T` POS-based region/target restriction, `-W`/`--write-index[=csi|tbi]` for VCF.gz/BCF outputs, `--threads` for VCF.gz/BCF file outputs, full `##bcftools_filter{Version,Command}` header line emission with `--no-version` suppression, Kestrel-tolerant text reads, shared `record_lookup` helper reused by `stats`. 28 integration tests in `crates/bcftools-rs/tests/filter.rs`.
  - [ ] Remaining: exact buffered gap-filter edge-case parity, full filter-expression FORMAT/sample-vector semantics, structured BCF write path that round-trips through the soft-filter rewrite without re-parsing.
- [ ] `annotate` (`vcfannotate.c`, 180k — single largest file in bcftools) — INFO/FORMAT/FILTER/ID column transfer from VCF/BCF/TAB sources, rename chrs, `-x` removal, header injection, `-c CHROM,POS,REF,ALT,…` column mapping, `--columns-file`, `--single-overlaps`, `--regions-overlap`. Covered by `test_vcf_annotate`.

### Wave D — Calling & Consequence

- [ ] `mpileup` (`mpileup.c`, 84k + `bam2bcf*.c`, `bam_sample.c`, `read_consensus.c`, `cigar_state.h`, `mw.h`) — multi-way pileup producing genotype likelihoods as BCF. Distinct from `samtools mpileup` (which produces text/VCF). Depends on `htslib-rs::alignment_compat` pileup iterators and `htslib-rs::probaln` for BAQ. Covered by `test_mpileup`.
- [ ] `call` (`vcfcall.c` + `mcall.c`, 65k + `ccall.c`, `prob1.c`, `em.c`) — multi-allelic (`-m`) and consensus (`-c`) callers, `--ploidy`, `--variants-only`, `--annotate FORMAT/PV4`, `-S` constrained allele set, gVCF mode. Covered by `test_vcf_call`, `test_vcf_call_cAls`.
- [ ] `consensus` (`consensus.c`, 55k) — apply VCF variants to FASTA reference, chain-file mode, `--missing`, `--mark-del/-ins/-snv`, `-H A/R/I/L` for haplotype selection, sample filters. Covered by `test_vcf_consensus`, `test_vcf_consensus_chain`.
- [ ] `csq` (`csq.c`, 166k) — variant consequence annotation given a GFF, supports phased calls, compound variants, splice consequences. Depends on Phase 1 `gff.rs` and reference helpers. Covered by `test_csq`, `test_csq_real`.

### Wave E — HMM / Stats / Trio

- [ ] `roh` (`vcfroh.c`, 52k) — HMM for runs of homozygosity, `--AF-dflt`, `--GTs-only`, `--estimate-AF`, viterbi/fwd-bwd. Depends on Phase 1 HMM. Covered by `test_roh`.
- [ ] `cnv` (`vcfcnv.c`, 60k) — HMM CNV calling from BAF/LRR. Depends on Phase 1 HMM and `peakfit`. Has no upstream `test.pl` coverage but the Rust gate must include synthetic integration tests.
- [ ] `gtcheck` (`vcfgtcheck.c`, 67k) — sample-concordance / contamination checks, `--pairwise`, `--dry-run`, `-e` for per-sample error rate. Covered by `test_gtcheck`.
- [ ] `polysomy` (`polysomy.c`, 34k) — chromosomal-copy detection. GPL-only upstream (uses GSL). Replace GSL with `statrs`/native; track as a separate milestone.
- [x] `som` (`vcfsom.c`, 25k) — experimental SOM-based filter. Defer (out of scope unless tests demand).

### Wave F — Plugins

All 41 plugins are ported as in-process modules under `crates/bcftools-rs/src/commands/plugins/<name>.rs`. They are invoked through `bcftools plugin <name>` and `bcftools +<name>`. The `+name` dispatch lives in the CLI crate; the `plugin` command's listing/help (`-l`, `-lv`, `-h`) walks a static plugin registry rather than scanning `BCFTOOLS_PLUGINS` for `.so` files.

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
- [ ] **Rust integration tests per subcommand**: under `crates/bcftools-rs/tests/<name>.rs`, write at least: happy path, error path, region/`-i`/`-e`/`-R`/`-T`/`-s` variants, threaded variant where applicable. These run on every PR independently of the Perl gate.
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
- [x] **CSI index 64-bit coordinate support** — `large_chrom_csi_limit` test in `test.pl:39` asserts the 2^31-1 boundary. Confirm htslib-rs CSI handles it.
- [x] **`hts_expr` vs bcftools filter** — `htslib-rs::expr` is the HTSlib expression language. bcftools has its own. Decide whether to also expose helpers from htslib-rs that bcftools's filter engine can reuse (numeric helpers, token utilities) or keep them fully separate. Document the decision.
- [x] **Region-with-target arithmetic** — `htslib-rs::region` covers HTSlib's grammar; confirm `-r`/`-R`/`-t`/`-T` semantics including the difference between regions (index-driven) and targets (streaming-filter) match upstream.
- [ ] **BCF serialization of haploid missing `GT=.`** — reverse
  `convert --hapsample2vcf -Ou` and `convert -H -Ou` hit
  `[E::main_vcfconvert] invalid input parameter` on the upstream Oxford
  fixtures when a haploid missing genotype (`GT=.`) is serialized through the
  current text-VCF-to-BCF writer path. Text VCF parity is correct, and
  GEN/SAMPLE `-G -Ou | view` passes; the remaining HAP/SAMPLE and
  HAP/LEGEND/SAMPLE upstream pipe fixtures need `htslib-rs`/writer support for
  this genotype shape.

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
