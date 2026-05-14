#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
binary="${BCFTOOLS_RS_BIN:-$repo_root/target/release/bcftools}"

if [[ ! -x "$binary" ]]; then
    echo "missing executable: $binary" >&2
    echo "build it with: cargo build --release -p bcftools-rs-cli" >&2
    exit 2
fi

stage="${TMPDIR:-/tmp}/bcftools-rs-test-pl.$$"
cleanup() {
    rm -rf "$stage"
}
trap cleanup EXIT

mkdir -p "$stage"
mkdir -p "$stage/test"
find "$repo_root/bcftools/test" -mindepth 1 -maxdepth 1 -exec ln -s {} "$stage/test/" \;
rm -f "$stage/test/test.pl"
cp "$repo_root/bcftools/test/test.pl" "$stage/test/test.pl"
ln -s "$binary" "$stage/bcftools"
ln -s "$binary" "$stage/bgzip"
cat > "$stage/tabix" <<'SH'
#!/usr/bin/env bash
set -euo pipefail

exec "$(dirname -- "$0")/bcftools" tabix "$@"
SH
chmod +x "$stage/tabix"

cd "$stage/test"
exec perl ./test.pl -e "bgzip=$stage/bgzip" -e "tabix=$stage/tabix" "$@"
