# Design Decisions

## bcftools Filter Expressions vs HTSlib Expressions

Decision: keep the bcftools filter expression engine as a separate
`bcftools-rs` module instead of adapting `htslib-rs::expr` as the primary
implementation.

Rationale:

- Upstream bcftools uses `filter.c`, not HTSlib's `hts_expr.c`, for `-i` and
  `-e` expressions.
- The bcftools language has observable behavior that is wider than
  `htslib-rs::expr`: sample-vector evaluation, genotype and AC/AN lazy
  caching, external values via `filter_test_ext`, and bcftools-specific status
  and unpack tracking.
- Subcommands such as `view`, `filter`, `query`, `isec`, `annotate`, `norm`,
  `stats`, `call`, and `mpileup` rely on those bcftools semantics for parity
  with `bcftools/test/test.pl`.

`htslib-rs::expr` can still be used as a reference for shared low-level
tokenization or numeric utilities if a future implementation finds a compatible
piece, but it is not the acceptance target and should not define bcftools-rs
filter behavior.

## Experimental `som` Subcommand

Decision: do not port `bcftools som` for the initial parity target unless an
upstream test or downstream consumer starts requiring it.

Rationale:

- Upstream marks `som` with a help string beginning with `-`, which keeps it
  out of the advertised command list.
- `som` is experimental, has a custom map-training/classification workflow,
  and is not a dependency of the current `bcftools/test/test.pl` parity slices.
- Keeping the dispatch-table alias as an unsupported hidden command preserves
  the upstream command-table shape without making this experimental path a
  release blocker for the core VCF/BCF commands.
