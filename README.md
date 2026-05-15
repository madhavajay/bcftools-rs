# bcftools-rs

`bcftools-rs` is a pure-Rust port of [`bcftools`](https://github.com/samtools/bcftools).
The goal is full subcommand and plugin parity with the upstream C program,
verified against the upstream `test/test.pl` suite plus Rust-native unit and
integration tests.

It does not link to HTSlib C, does not use `bindgen` or `cc`, and does not
ship the upstream C plugin ABI. All HTSlib-shaped behavior routes through the
sibling [`htslib-rs`](htslib-rs/) workspace.

## Repository layout

| Path | Purpose |
| --- | --- |
| `crates/bcftools-rs/` | Library: one module per subcommand (`commands/`) and one per plugin (`commands/plugins/`), plus shared infrastructure |
| `crates/bcftools-rs-cli/` | The `bcftools` binary — dispatches `argv[1]` to the matching subcommand exactly like upstream `main.c` |
| `bcftools/` | Upstream C source + tests (submodule). Used only as fixture/reference; never built or linked at runtime |
| `htslib-rs/` | Sibling pure-Rust HTSlib compatibility workspace (submodule); consumed via path dep |
| `TODO.md` | The phased porting plan |

## Status

Phase 0 (workspace skeleton), shared infrastructure, and the first local-only
parity batches are in. The project is still far from full upstream parity, but
these commands have meaningful Rust-native coverage:

| Subcommand | Coverage |
| --- | --- |
| `bcftools head` | Full upstream surface: `-h N`, `-n N`, `-s N`, `-v N`. VCF + VCF.gz + BCF input, including stdin. |
| `bcftools index` | BCF CSI and VCF.gz CSI/TBI builds, stdin indexing with `-o`, overwrite protection, `--stats`, `--nrecords`, and large-coordinate CSI fixture coverage. |
| `bcftools tabix` | Preset BGZF indexes and queries for VCF/BED/GFF/SAM, plus all-record streaming. Custom `-s/-b/-e/-0/-S/-c` text layouts are dependency-blocked on `htslib-rs` API surface. |
| `bcftools view` | VCF/VCF.gz/BCF reads and VCF/BGZF/BCF writes, stdin spooling, header modes, `--no-version`, sample subsetting, simple region/target restriction, many simple site filters, limited expression filtering, Kestrel-compatible VCF headers, and threaded writes. |
| `bcftools query` | Sample listing and selection, POS-based regions/targets, a text-backed subset of record/sample expressions, and a growing subset of the `convert.c` formatter including sample loops, numeric functions, `%N_PASS(...)`, and `%PBINOM(...)`. |
| `bcftools sort` | Coordinate sorting with disk-backed temp-run spill, VCF/BGZF output, automatic indexing, Kestrel-compatible VCF headers, and threaded BGZF writes. |
| `bcftools concat` | Same-sample vertical concat for VCF/VCF.gz/BCF inputs, file lists, genotype dropping, duplicate removal, naive concat, region restriction, indexing, Kestrel-compatible headers, and threaded writes. |
| `bcftools convert` | Focused TSV/23andMe, gVCF, GEN/SAMPLE, HAP/SAMPLE, and HAP/LEGEND/SAMPLE conversion paths with fixture-backed text parity, BCF stdin/output paths, indexing, filtering hooks, and sample selection. |
| `bcftools filter` | Text-backed expression filtering, soft-filter tagging, masks, gap filters, simple genotype rewrites, region/target restriction, indexing, Kestrel-compatible headers, and threaded writes. |
| `bcftools isec` | Text-backed set intersections/complements, collapse modes for common fixtures, target/region filtering, prefix output, directory output, record-output selection, indexing, and Kestrel-compatible reads. |
| `bcftools reheader` | VCF/BGZF VCF header replacement, sample rename, FAI contig updates, stdin handling, BCF output, BCF `--in-place`, and threaded output. |
| `bcftools stats` | Substantial single-input and pairwise text-backed stats sections, sample selection, AF/depth/user-TSTV options, expression filtering, regions/targets, and selected indel-context/exon summaries. |

The remaining large subcommands (`annotate`, `merge`, `norm`, `mpileup`,
`call`, `consensus`, `csq`, `roh`, `cnv`, `gtcheck`) and most plugins are not
ported yet. Detailed parity status is tracked in [`TODO.md`](TODO.md), and
upstream Perl harness enablement is tracked in
[`docs/test-status.md`](docs/test-status.md).

## Build and test

Clone with submodules:

```sh
git clone --recurse-submodules git@github.com:madhavajay/bcftools-rs.git
```

Run the Rust gate:

```sh
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
```

This produces a `bcftools` binary at `target/debug/bcftools` (or `target/release/bcftools`).

Run selected upstream Perl parity tests against the Rust binary:

```sh
cargo build --release -p bcftools-rs-cli
scripts/run-bcftools-test-pl.sh -f '^(test_vcf_head|test_vcf_head2|test_tabix|test_index|test_vcf_idxstats)$'
```

The upstream harness does not expose a direct `--bin` option for `bcftools`
itself. It derives `$$opts{bin}` from the parent of `bcftools/test`, so the
wrapper stages the expected layout with a real `test/` directory containing
symlinked upstream fixtures and `bcftools` symlinked to the Rust binary. Pass
`-f` with one or more `test.pl` function names or regular expressions to run
only the subcommands currently ported; anchor names when needed because the
harness matches function names as regexes.

Rust integration tests locate upstream fixtures relative to
`CARGO_MANIFEST_DIR` at `../../bcftools/test/<name>`, including nested fixture
sets such as `csq/` and `mpileup/`. This mirrors `test.pl`'s `$$opts{path}`
substitution for arguments containing `{PATH}`.

## CI

Two jobs run in CI ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)),
mirroring the split used by `htslib-rs`:

- **Rust gate** — `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
  and `cargo test --workspace`. Includes integration tests that run the built
  `bcftools` binary as a subprocess against fixtures from `bcftools/test/`.
- **bcftools Perl parity tests** — builds the Rust release binary and runs the
  enabled upstream `bcftools/test/test.pl` slices through
  `scripts/run-bcftools-test-pl.sh`. Enabled and disabled Perl functions are
  tracked in [`docs/test-status.md`](docs/test-status.md).

## Submodule pinning

- **`bcftools/`** is pinned at upstream tag `1.23.1`, commit
  `ff1604d4622dc715a921f8e21e0e5d88438d10d1`. This is the version emitted in
  `##bcftools_<cmd>Version=...+htslib-...` header lines.
- **`htslib-rs/`** is pinned at commit
  `6bd6fb051ee7898c2afa4e619bf99dbad5f60dd7` on `main`. The two
  consumer-driven extensions added for bcftools-rs are merged upstream:
  - `index_compat::build_vcf_csi_from_path` /
    `build_vcf_csi_from_path_with_min_shift` /
    `build_vcf_tbi_from_path` — index existing BGZF-compressed VCFs without
    rewriting the data (analogous to HTSlib's `tbx_index_build3`).
  - `header_compat::append_other_record` / `append_other_records` /
    `append_line` — string-keyed `##key=value` injection (analogous to
    HTSlib's `bcf_hdr_append`).

## Design rules

These are decided and shape every subcommand and plugin:

- **Stay pure Rust.** No `bindgen`, no `cc` crate, no linking to HTSlib C or
  to bcftools C.
- **Route HTSlib-shaped APIs through `htslib-rs`.** When `htslib-rs` lacks
  what's needed, extend it in a feature branch first and pin the submodule
  to that branch — do not bypass with direct `noodles` calls just because
  it's faster (escape hatch is allowed only for behavior with no HTSlib
  analogue).
- **Two test gates.** Rust unit/integration tests under
  `crates/bcftools-rs/tests/` run on every PR. The Perl parity gate
  (`bcftools/test/test.pl` driven against the Rust binary) lands per
  subcommand as Wave A → F lands.
- **Strict byte parity** for VCF/BCF binary output, index bytes, sort order,
  text output from `stats`/`query`/`gtcheck`/`roh`/`cnv`/`csq`/`isec`, and
  exit codes. **Semantic parity** for `##bcftools_<cmd>Version` /
  `##bcftools_<cmd>Command` lines, stderr text, and usage/help text.
- **`--no-version`** suppresses the `##bcftools_<cmd>{Version,Command}` lines
  exactly like upstream — `bcftools/test/test.pl` uses this pervasively, so
  the suppression path is a hard requirement.
- Tests that intentionally keep `##bcftools_<cmd>Command` lines can pin the
  generated `Date=` field with `BCFTOOLS_RS_FIXED_TIME=<unix-seconds>`.

## License

MIT. The vendored upstream `bcftools/` source remains under its upstream
license (MIT/Expat unless `USE_GPL` is enabled).
