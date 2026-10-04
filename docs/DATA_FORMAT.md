# Atlas storage and commitments

All binary integers are little-endian. SHA256 commitments use the exact bytes
specified below. The original read-only arithmetic codec is extracted into
`scripts/atlas_codec.py`; its source identities are in `CODE_EXTRACTION.json`.

## Scan packs

The 62 `facts-*.pack` files concatenate frames. A scan header uses Python's
`struct.Struct('<8sQQQQII32s')`, in this order:

| Field | Meaning |
| --- | --- |
| 8 bytes | Magic `VDWSM1` followed by two zero bytes |
| uint64 | Compressed payload length |
| uint64 | Uncompressed payload length |
| uint64, uint64 | Interval lo, hi (hi exclusive) |
| uint32 | Number of prime records |
| uint32 | Witness flag, equal to 1 in this atlas |
| 32 bytes | SHA256 of uncompressed payload |

The payload is zlib-compressed. Each record contains a positive prime delta,
one admissibility mask byte, and M,L,witness for every set bit in ascending
r = 2,…,9 order. Bit r−2 is set exactly when r divides p−1. The first prime is
relative to lo−1; later primes are relative to the preceding prime. Unsigned
integers use canonical base-128 varints, at most 64 bits. K is derived by the
scientific contract rather than redundantly stored.

The prime commitment hashes `b'vdw-primes-u64le'` followed by each p as uint64.
The profile commitment hashes `b'vdw-p-profile-v1-u64le'`, then for each prime
`struct.pack('<QB', p, mask)` followed by each admissible `(M,L,K)` encoded with
`struct.pack('<QQQ', ...)`. Witnesses are protected by frame and pack hashes;
they are not part of that profile commitment.

## A4 packs

The eight `A4-*.pack` files concatenate frames with header
`struct.Struct('<8sQQIIII32s32s32s')`:

```text
magic = b'VDWA4V1\0', lo, hi, number_of_primes,
compressed_bytes, raw_bytes, flags = 0,
prime_hash, profile_hash, values_hash
```

The zlib payload contains one uint32 A4 value per prime in the corresponding
scan block. `raw_bytes = 4 * number_of_primes`; `values_hash` is its SHA256.

## Publication indices

Each group has `ledger.sqlite`. Use read-only SQLite connections. The `blocks`
tables retain the original block id, lo, hi, n, pack basename, byte offset,
frame length and scientific commitments. Both have zlib-compressed JSON block
summaries. The `audits` tables retain the mathematical evidence. The scan index
also retains `best` and `milestones`; both retain scientific `meta` entries.

These two databases were rebuilt into fresh files, including a SQLite integrity
check. The science columns were compared row by row with the originals. The
manifest lists their original and publication hashes and omitted summary/meta
keys. Audit wall-clock durations are also omitted. No old SQLite pages are
copied into a publication database. Pack lengths and hashes are unchanged.

## Chains

Both chains begin with 32 zero bytes. For every block in ascending id order:

```text
scan = SHA256(previous_scan || prime_hash || profile_hash || uint64(lo) || uint64(hi))
A4   = SHA256(previous_A4 || prime_hash || profile_hash || values_hash || uint64(lo) || uint64(hi))
```

Their final hexadecimal values are:

```text
scan: f1b8309688bf345077c3373dff44e875498ea1cc748146711a48e03e20209433
A4:   cc45f668ee4c5a78dfc71969b94efda123da661b2f6690bb283a4ea7032c5b56
```

The original sealed archive manifest identity is retained as a provenance
anchor. The publication subset is not the original full archive and does not
claim to reproduce that original manifest byte for byte.
