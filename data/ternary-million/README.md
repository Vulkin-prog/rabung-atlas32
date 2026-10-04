# Fixed million ternary profiles: counts and long-run positions

This is the existing descriptive supplement sealed on 2 October 2026, curated
for redistribution. It contains **1,000,000 exact ternary profiles** and their
long-run counts in **[2^33, 2^34)**. It is a probability sample of primes
p = 1 mod 6, separate from the complete atlas below 2^32. It does not extend
the domain of that complete census.

All scientific packs and prime lists are committed here; no release download,
Git LFS, Python dependency or GPU is needed to read them. The scientific data
occupy about 35.5 MiB, including their parent profiles. [manifest.json](manifest.json)
records every file's SHA256 and byte length, original identities and publication
operations. The original scientific archives have not been modified.

| Stored quantity | Exact total |
| --- | ---: |
| Sampled primes / ternary profiles | 1,000,000 |
| Count and position blocks | 248 |
| Canonical maximal runs of length at least 18 | 33,257,058 |
| Retained canonical run positions of length at least 21 | 1,235,473 |
| Criterion (a) at k = 21: M < 21 | 299,150 |
| Criteria (a)+(b) at k = 21: K <= 21 | 293,872 |
| Primes with witnessed doubling edges | 390 |
| Directed doubling edges | 398 |

**Positions of runs of lengths 18, 19 and 20 are not retained for the full
million.** Their counts, colour totals and position-bin totals are retained.
The separate 512 diagnostic complete catalogues were used in the historical
qualification; those primes are included in this million. This supplement is
descriptive after exploration, not an independent validation or a reopened
primary test. No coefficient is refitted here.

## Read, verify and export

Run these commands from the repository root:

```bash
python3 -B scripts/verify_ternary.py
python3 -B scripts/query_ternary.py 8611114591
python3 -B scripts/query_ternary.py --stratum 0 --index 1241
```

Verification covers every pack hash, all one million records, complete coverage
of the four fixed lists, all maxima and prefixes against their parent profiles,
counts against positions, the integer count and factorial-moment tables, and
all doubling graphs against their stored witnesses. It repeats full modular
character and maximality checks for **1,189 runs** at four fixed indices in
each block. It takes about 13 seconds on the preparation machine. This audit
does not enumerate all p-1 character values for every sampled prime; full
maximum completeness rests on the historical qualified scientific computation.

The portable reader checks at most 4,096 profiles in a block. Query indices are
zero-based within a stratum. A prime outside the fixed lists is reported as
absent, rather than computed. Query output includes M, L, K, counts at 18..27,
colours, echoes, position bins, and the retained start/length/colour triples.

For ordinary CSV files, export to an ignored working directory:

```bash
mkdir -p build/ternary_exports
python3 -B scripts/query_ternary.py --export counts > build/ternary_exports/counts.csv
python3 -B scripts/query_ternary.py --export positions > build/ternary_exports/positions21.csv
python3 -B scripts/query_ternary.py --stratum 0 --export positions --minimum 24 > build/ternary_exports/stratum0_positions24.csv
```

The full count export has one million data rows. The full position export has
1,235,473 data rows. Exports stream through bounded blocks and can be much
larger than the compressed source. `--minimum` filters retained positions and
must be at least 21; it cannot reconstruct shorter positions. CSV colours 1
and 2 follow the per-prime worker root convention described below.

[summary.json](summary.json) and [threshold_counts.csv](threshold_counts.csv)
give regenerated exact totals, histograms and factorial-moment sums. The
historical descriptive comparisons and graph witness evidence are in
[analysis/](analysis/). Old source hashes inside redacted analysis JSON are
historical anchors. The hashes of the published files are separately recorded
in the manifest; they must not be substituted for those original anchors.

## Four sampling strata

There are 250,000 primes per stratum, sorted within each list:

| Stratum | Band, with upper endpoint excluded | Is 2 a cube modulo p? |
| --- | --- | --- |
| 0 | [8,589,934,592, 12,884,901,888) | Yes |
| 1 | [8,589,934,592, 12,884,901,888) | No |
| 2 | [12,884,901,888, 17,179,869,184) | Yes |
| 3 | [12,884,901,888, 17,179,869,184) | No |

The original strata were sampled after excluding the conservative prior
profile inventory. [SAMPLING.json](SAMPLING.json) retains the sampling rule,
band seeds, list hashes and parent protocol identity. The accepted lists are
fully included. Historical sampler traces and the earlier exclusion inventory
are not bundled, so this publication subset does not claim to regenerate the
selection from the seeds alone. The million is not an exhaustive prime scan.

## What a run and a position mean

The ternary character is x^((p-1)/3) modulo p, for x = 1,...,p-1.
A run is a maximal interval of consecutive integers on which this value is
constant. Because p = 1 mod 6, reflection x -> p-x preserves the character.
Each reflected pair is stored once, with the smaller start x in the first
half-period. A central run crossing (p-1)/2 is stored once with its full even
length; its start satisfies 2x + length = p+1. The reflected start of an
interval is p-x-length+1. All positions are one-based integer residues.

An **ordinary** run excludes the initial run starting at 1 and the merged
central run. The counts N18,...,N27 in the fixed record count ordinary canonical
maximal runs of length **at least** the threshold, not overlapping windows or
individual starts within the same run. Canonical counts add the initial and
central runs when their lengths reach the threshold. Counts for the full
unfolded interval can be recovered by doubling canonical ordinary and initial
counts and adding the central run once.

M is the exact full maximum, L the exact initial length. For these primes the
Rabung threshold is K = max(M+1, 2L+2). The stored central field is the full
central length when it is at least 18, and **zero if it is below 18**; zero does
not assert that the central run is absent. The position catalogue includes
initial, ordinary and central canonical runs of length at least 21.

Colour 0 denotes the value 1. Colours 1 and 2 denote z and z^2 for the primitive
cubic root chosen by the historical worker for that prime. That root is not
stored in these compact count packs. For any retained nonneutral run it can
be recovered exactly from its start: if t=x^((p-1)/3), then z=t for colour 1
and z=t^2 for colour 2. An orientation of colours 1 and 2 must not be compared
across primes without recovering their roots. Numeric character values are
unambiguous; neutrality and equality within a profile are unambiguous.

The 16 bins count ordinary runs at 18 using
min(15, floor(16*x/((p-1)/2))). Echo counts at 18 are the historical
first-18-window containment echoes under multiplication by 2 and 3, with the
source fold restricted to one monotone branch. They differ from the graph
edges on full ordinary runs at 21: such an edge has an actual integer point
whose doubled folded image lies in the target full run. Self-edges are kept.
The graph records enumerate directed simple cycles of lengths 1, 2 and 3;
absence of longer cycles is not claimed.

## Binary format and commitments

`primes/stratum-s.u64` contains 250,000 little-endian unsigned 64-bit integers.
The pack `blocks/sS-NNNNNN.pack` is one zlib stream. Its decompressed payload
is a little-endian uint32 fixed-array byte length, the fixed array, then the
concatenated position stream. Each prime has a **68-byte `<34H`** record:

| Zero-based fields | Meaning |
| --- | --- |
| 0, 1 | M, L |
| 2 | Full central length if >=18; otherwise 0 |
| 3..12 | Ordinary counts N18..N27 |
| 13..15 | Ordinary colour counts at 18 |
| 16, 17 | First-window echo counts at 18 under 2 and 3 |
| 18..33 | Sixteen ordinary position-bin counts at 18 |

For each prime, the position stream begins with a canonical unsigned varint
run count. Each run then stores a positive varint start difference (previous
start initially zero), a varint `length-21`, and a one-byte colour 0..2.
Varints use seven payload bits per byte and the high continuation bit. The
reader rejects noncanonical encodings, overflow, overlapping intervals,
trailing bytes and extra zlib streams.

Metadata retains the original compressed-pack, fixed-count, position-stream,
prime and scientific hashes. The scientific commitment is SHA256 of:

```text
ASCII("vdw-count-catalogue-extension-v1") || original_protocol_sha256_bytes
|| prime_list_slice_u64le || fixed_counts || position_stream
```

The prime commitment is SHA256 of
`ASCII("vdw-compact-primes-u64le-v1") || prime_list_slice_u64le`.
[PROTOCOLE.json](PROTOCOLE.json) is the original unchanged protocol, preserving
the scientific commitment identity. Its dates and execution bounds document
the historical computation; reading it authorizes no new run.

The 248 byte-identical `parents/*.pack` contain the original `<HHIBII>`
records (M, L, normalized cubic signature, flags, rational crossing masks at
18 and 21). Their header is `<8sIIIII32s32s16s>`, with magic `VDWCSP1\0`,
stratum, list start, record count, compressed and raw byte lengths, prime
hash, exact-field hash and the first 16 bytes of the frozen parent identity.
The exact-field hash is SHA256 of
`ASCII("vdw-sample-fields-v1") || prime_list_slice_u64le || parent_raw`.
Flag bits 0 and 1 give the central eligibility gates at 18 and 21; bit 2
is the cube class of 2. The extra signature and masks are retained for parent
integrity; their original encoding is described by the collector sources.

## Collector sources and scope

The [frozen sparse collector](../../src/ternary_catalogue/) includes its exact
58 source files, compile definitions, coverage argument and historical
qualification evidence. The count/capture reduction and doubling graph
functions were extracted unchanged into `scripts/ternary_capture.py` and
`scripts/ternary_graph.py`. The portable reader and verification driver are
new publication code, with negative corruption tests and CI coverage.

CUDA compilation is optional and explicit:

```bash
python3 -B scripts/build_ternary_worker.py
```

This compiles an SM120 worker and starts no calculation. Compilation and a
fresh million-profile regeneration have **not** been rerun for this repository
addition. Historical worker binaries, host logs, launchers and received review
reports are excluded. The original qualification covers the fixed-threshold
ternary catalogue below 2^34; it is not a capacity guarantee for arbitrary
future inputs or a universal statistical law. Rights and visibility follow
the repository's [RIGHTS.md](../../RIGHTS.md).
