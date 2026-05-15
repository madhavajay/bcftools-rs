//! VariantKey core algorithm (upstream `bcftools/variantkey.h`).
//!
//! Faithful port of the MIT-licensed VariantKey reference implementation by
//! Nicola Asuni / GENOMICS plc, used by the `+add-variantkey` and
//! `+variantkey-hex` plugins. It produces a 64-bit key from CHROM, 0-based
//! POS, REF and ALT, fully reversible for variants with up to 11 ACGT bases
//! and hash-encoded otherwise (low bit set on the hash path). All integer
//! arithmetic mirrors the C `uint8_t`/`uint32_t` wrapping semantics exactly.

const VKSHIFT_CHROM: u32 = 59;
const VKSHIFT_POS: u32 = 31;
const MAXUINT32: u32 = 0xFFFF_FFFF;

/// Numeric chromosome encoding (`encode_numeric_chrom`). Returns 0 if any
/// non-digit byte is found; the running value wraps like the C `uint8_t`.
fn encode_numeric_chrom(chrom: &[u8]) -> u8 {
    let mut v = chrom[0].wrapping_sub(b'0');
    for &c in &chrom[1..] {
        if !c.is_ascii_digit() {
            return 0;
        }
        v = v.wrapping_mul(10).wrapping_add(c - b'0');
    }
    v
}

/// True if `chrom` has a case-insensitive `chr` prefix and is longer than 3
/// bytes (`has_chrom_chr_prefix`).
fn has_chrom_chr_prefix(chrom: &[u8]) -> bool {
    chrom.len() > 3
        && (chrom[0] == b'c' || chrom[0] == b'C')
        && (chrom[1] == b'h' || chrom[1] == b'H')
        && (chrom[2] == b'r' || chrom[2] == b'R')
}

fn onecharmap(c: u8) -> u8 {
    match c {
        b'X' | b'x' => 23,
        b'Y' | b'y' => 24,
        b'M' | b'm' => 25,
        _ => 0,
    }
}

/// Chromosome numerical encoding (`encode_chrom`). Returns 0 on invalid input.
pub fn encode_chrom(mut chrom: &[u8]) -> u8 {
    if has_chrom_chr_prefix(chrom) {
        chrom = &chrom[3..];
    }
    if chrom.is_empty() {
        return 0;
    }
    if chrom[0].is_ascii_digit() {
        return encode_numeric_chrom(chrom);
    }
    if chrom.len() == 1 || (chrom.len() == 2 && (chrom[1] == b'T' || chrom[1] == b't')) {
        return onecharmap(chrom[0]);
    }
    0
}

/// Encode a single base: A=0, C=1, G=2, T=3 (case-insensitive), else 4.
fn encode_base(c: u8) -> u32 {
    match c {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

/// Encode an allele into `h` starting at `*bitpos`. Returns false on a
/// non-ACGT base, matching the C `-1` error return.
fn encode_allele(h: &mut u32, bitpos: &mut u8, s: &[u8]) -> bool {
    for &c in s {
        let v = encode_base(c);
        if v > 3 {
            return false;
        }
        *bitpos -= 2;
        *h |= v << *bitpos;
    }
    true
}

/// Reversible REF+ALT encoding (`encode_refalt_rev`); `MAXUINT32` on error.
fn encode_refalt_rev(reference: &[u8], alt: &[u8]) -> u32 {
    let mut h: u32 = 0;
    h |= (reference.len() as u32) << 27;
    h |= (alt.len() as u32) << 23;
    let mut bitpos: u8 = 23;
    if !encode_allele(&mut h, &mut bitpos, reference) || !encode_allele(&mut h, &mut bitpos, alt) {
        return MAXUINT32;
    }
    h
}

/// MurmurHash3-like 32-bit mix (`muxhash`).
fn muxhash(mut k: u32, mut h: u32) -> u32 {
    k = k.wrapping_mul(0xcc9e_2d51);
    k = k.rotate_left(15); // == (k >> 17) | (k << 15)
    k = k.wrapping_mul(0x1b87_3593);
    h ^= k;
    h = h.rotate_left(13); // == (h >> 19) | (h << 13)
    h.wrapping_mul(5).wrapping_add(0xe654_6b64)
}

fn encode_packchar(c: u8) -> u32 {
    if c < b'A' {
        return 27;
    }
    if c >= b'a' {
        return (c - b'a' + 1) as u32;
    }
    (c - b'A' + 1) as u32
}

/// Pack the trailing 1..=5 characters (`pack_chars_tail`).
fn pack_chars_tail(s: &[u8]) -> u32 {
    let size = s.len();
    let mut h: u32 = 0;
    // pos walks backwards from the last byte, matching the C `switch`
    // fall-through which consumes the tail right-to-left.
    let mut idx = size;
    if size >= 5 {
        idx -= 1;
        h ^= encode_packchar(s[idx]) << (1 + 5);
    }
    if size >= 4 {
        idx -= 1;
        h ^= encode_packchar(s[idx]) << (1 + 5 * 2);
    }
    if size >= 3 {
        idx -= 1;
        h ^= encode_packchar(s[idx]) << (1 + 5 * 3);
    }
    if size >= 2 {
        idx -= 1;
        h ^= encode_packchar(s[idx]) << (1 + 5 * 4);
    }
    if size >= 1 {
        idx -= 1;
        h ^= encode_packchar(s[idx]) << (1 + 5 * 5);
    }
    h
}

/// Pack a full block of 6 characters (`pack_chars`).
fn pack_chars(s: &[u8]) -> u32 {
    (encode_packchar(s[5]) << 1)
        ^ (encode_packchar(s[4]) << (1 + 5))
        ^ (encode_packchar(s[3]) << (1 + 5 * 2))
        ^ (encode_packchar(s[2]) << (1 + 5 * 3))
        ^ (encode_packchar(s[1]) << (1 + 5 * 4))
        ^ (encode_packchar(s[0]) << (1 + 5 * 5))
}

/// 32-bit hash of a nucleotide string (`hash32`).
fn hash32(s: &[u8]) -> u32 {
    let mut h: u32 = 0;
    let mut rest = s;
    while rest.len() >= 6 {
        h = muxhash(pack_chars(&rest[..6]), h);
        rest = &rest[6..];
    }
    if !rest.is_empty() {
        h = muxhash(pack_chars_tail(rest), h);
    }
    h
}

/// Non-reversible hash REF+ALT encoding (`encode_refalt_hash`).
fn encode_refalt_hash(reference: &[u8], alt: &[u8]) -> u32 {
    let mut h = muxhash(hash32(alt), muxhash(0x3, hash32(reference)));
    h ^= h >> 16;
    h = h.wrapping_mul(0x85eb_ca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2_ae35);
    h ^= h >> 16;
    (h >> 1) | 0x1
}

/// REF+ALT numerical encoding (`encode_refalt`).
fn encode_refalt(reference: &[u8], alt: &[u8]) -> u32 {
    if reference.len() + alt.len() <= 11 {
        let h = encode_refalt_rev(reference, alt);
        if h != MAXUINT32 {
            return h;
        }
    }
    encode_refalt_hash(reference, alt)
}

fn encode_variantkey(chrom: u8, pos: u32, refalt: u32) -> u64 {
    ((chrom as u64) << VKSHIFT_CHROM) | ((pos as u64) << VKSHIFT_POS) | (refalt as u64)
}

/// 64-bit VariantKey from CHROM, 0-based POS, REF and ALT.
pub fn variantkey(chrom: &[u8], pos: u32, reference: &[u8], alt: &[u8]) -> u64 {
    encode_variantkey(encode_chrom(chrom), pos, encode_refalt(reference, alt))
}

/// Lowercase, zero-padded 16-hex-digit VariantKey string (`variantkey_hex`).
pub fn variantkey_hex(vk: u64) -> String {
    format!("{vk:016x}")
}

/// Parse an `rs<digits>` ID into the upstream 32-bit RSX value: drop the
/// first two bytes (`rs`), then `strtoul(..., 10)` truncated to `uint32_t`.
pub fn rsid_u32(id: &str) -> u32 {
    let bytes = id.as_bytes();
    if bytes.len() < 2 {
        return 0;
    }
    let mut rest = &bytes[2..];
    while let Some((&c, tail)) = rest.split_first() {
        if c.is_ascii_whitespace() {
            rest = tail;
        } else {
            break;
        }
    }
    let mut acc: u64 = 0;
    let mut saw_digit = false;
    for &c in rest {
        if !c.is_ascii_digit() {
            break;
        }
        saw_digit = true;
        acc = acc.wrapping_mul(10).wrapping_add((c - b'0') as u64);
    }
    if !saw_digit {
        return 0;
    }
    acc as u32
}

/// Lowercase, zero-padded 8-hex-digit RSX string.
pub fn rsx_hex(rs: u32) -> String {
    format!("{rs:08x}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reversible_key_matches_upstream_fixture() {
        // query.add-variantkey.vcf record 1: chr "1", POS 10019 (0-based
        // 10018), REF "TA", ALT "T" -> VKX=0800139110e60000 RSX=2e3deb1d.
        let vk = variantkey(b"1", 10018, b"TA", b"T");
        assert_eq!(variantkey_hex(vk), "0800139110e60000");
        assert_eq!(rsx_hex(rsid_u32("rs775809821")), "2e3deb1d");
        assert_eq!(vk & 1, 0, "short ACGT variant must be reversible");
    }

    #[test]
    fn second_fixture_record() {
        // POS 10055 -> 0-based 10054, REF "T", ALT "TA".
        let vk = variantkey(b"1", 10054, b"T", b"TA");
        assert_eq!(variantkey_hex(vk), "080013a309780000");
        assert_eq!(rsx_hex(rsid_u32("rs768019142")), "2dc70ac6");
    }

    #[test]
    fn non_reversible_hash_records_match_fixture() {
        // REF+ALT length > 11 -> hash mode, low bit set (non-reversible).
        let vk = variantkey(b"1", 10227, b"TAACCCCTAACCCTAACCCTAAACCCTA", b"T");
        assert_eq!(variantkey_hex(vk), "080013f9a00e1d03");
        assert_eq!(rsx_hex(rsid_u32("rs200462216")), "0bf2cf88");
        assert_eq!(vk & 1, 1);

        let vk = variantkey(b"1", 10327, b"AACCCCTAACCCTAACCCTAACCCT", b"A");
        assert_eq!(variantkey_hex(vk), "0800142b90367897");
        assert_eq!(rsx_hex(rsid_u32("rs201106462")), "0bfca41e");

        let vk = variantkey(b"1", 10615, b"CCGCCGTTGCAAAGGCGCGCCG", b"C");
        assert_eq!(variantkey_hex(vk), "080014bb8ad3d64f");
        assert_eq!(rsx_hex(rsid_u32("rs376342519")), "166e87f7");
    }

    #[test]
    fn chrom_encoding_special_cases() {
        assert_eq!(encode_chrom(b"1"), 1);
        assert_eq!(encode_chrom(b"22"), 22);
        assert_eq!(encode_chrom(b"X"), 23);
        assert_eq!(encode_chrom(b"Y"), 24);
        assert_eq!(encode_chrom(b"MT"), 25);
        assert_eq!(encode_chrom(b"chr7"), 7);
        assert_eq!(encode_chrom(b"chrX"), 23);
        assert_eq!(encode_chrom(b""), 0);
        assert_eq!(encode_chrom(b"GL000"), 0);
    }

    #[test]
    fn rsid_edge_cases() {
        assert_eq!(rsid_u32("."), 0);
        assert_eq!(rsid_u32("rs"), 0);
        assert_eq!(rsid_u32("rs0"), 0);
        assert_eq!(rsid_u32("rs123abc"), 123);
    }
}
