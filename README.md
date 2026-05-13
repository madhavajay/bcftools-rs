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

Phase 0 (workspace skeleton) and the Wave A foundation are in. Implemented
subcommands:

| Subcommand | Coverage |
| --- | --- |
| `bcftools head` | Full upstream surface: `-h N`, `-n N`, `-s N`, `-v N`. VCF + VCF.gz + BCF input, including stdin. |
| `bcftools index` | BCF→CSI and VCF.gz→CSI/TBI. `-c`/`-t`/`-m`/`-o`/`-f`/`-v`, plus `--stats`/`--nrecords` from existing CSI/TBI metadata. |
| `bcftools view` | I/O backbone: `-O v\|z\|u\|b`, `-o`, `-h` (header-only), `-H` (no-header), `--no-version`, positional `CHROM` / `CHROM:START-END` region filtering. Filtering expressions and sample subsetting not yet wired. |

Subcommand parity beyond Wave A is tracked in [`TODO.md`](TODO.md).

## Build and test

Clone with submodules:

```sh
git clone --recurse-submodules git@github.com:madhavajay/bcftools-rs.git
```

Run the Rust gate:

```sh
cargo fmt --all -- --check
cargo clippy -p bcftools-rs -p bcftools-rs-cli --all-targets -- -D warnings
cargo test -p bcftools-rs -p bcftools-rs-cli
```

This produces a `bcftools` binary at `target/debug/bcftools` (or `target/release/bcftools`).

Run selected upstream Perl parity tests against the Rust binary:

```sh
cargo build --release -p bcftools-rs-cli
scripts/run-bcftools-test-pl.sh -f '^test_vcf_head$'
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
  `56ddf62df73efe96a3a906081ca50fbc3a350b70` on `main`. The two
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
