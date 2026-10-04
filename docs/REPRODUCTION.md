# Reproduction routes

## Reading and checking the published data

Use `download_data.py` to retrieve immutable, pinned release assets. It verifies
every size and SHA256 before atomically installing a file. The manifest is the
commitment to all scientific pack bytes and to the derived indices.

`verify_atlas.py --data data/atlas` checks all 72 hashes, complete interval
coverage, block cardinalities and byte extents, paired scan/A4 identities,
chains, histogram aggregates, total A4 and periodic audit counts. It decodes
the two endpoints as a bounded extra check. It does not independently repeat
primality tests and the GPU algorithm for 203 million records.

`reproduce_analysis.py --data data/atlas` performs complete frame decoding with
the original C++ analysis and compares all descriptive aggregates. Expected
compiler dependencies on Ubuntu are `g++`, `zlib1g-dev` and `libssl-dev`.
The script never downloads dependencies or runs a GPU worker automatically.
Other original readers for tail, overlap, joint and phase investigations are
provided next to their scientific evidence.

## Independent bounded profile and A4 regeneration

Build with `python3 -B scripts/build_workers.py`. The CPU demonstration
`recompute_small.py --count 256` compares independently computed dense profiles,
prefixes, reflected witnesses and exact A4 values with the first committed block.
Counts can range from 1 to 4,096. There is a 30-second subprocess limit.

The persistent dense profile worker (`build/small_cpu`) accepts:

```text
0 2 3 5 7 11
1 13 17 19
QUIT
```

Batch ids start at zero. Each batch contains 1..128 primes; dense CPU input is
limited to p ≤ 2,000,000. Output is line-delimited JSON (`ready`, `measurement`,
`batch_done`, `bye`). Each measurement has p and the admissible profiles.

The exact A4 worker (`build/moment_worker`) accepts a different protocol:

```text
batch_id audit_index number_of_primes ordered_primes...
QUIT
```

The audit index is −1 for no certificate, or a zero-based prime index. Batches
contain at most 4,096 strictly increasing primes below 2^32. After the initial
JSON ready line, output is binary: `struct.Struct('<8sQIIQQIIQQ')`, then uint32
A4 values, then certificate JSON. The header fields are `A4BATCH1`, batch id,
count, certificate length, compute nanoseconds, sum A4, maximum search attempts,
small dense prime count, first prime, last prime. The certificate checker
`src/census/theorie/certify_elliptic.py` supplies the independent arithmetic.

## Historical GPU sources

The frozen source tree is copied with its original scientific source hashes.
Historical binaries are deliberately not redistributed. `GPU_CONFIG.json`
records the original SM120 macros, compiler configuration and the later
capacity policy. CUDA headers and CUB are external dependencies, not bundled
copies. Rebuilt executables may have different binary hashes across toolchains.

`build_workers.py --gpu` compiles the baseline bridge and persistent worker for
SM120 using the recorded macros. This optional build has not been retested on
every compiler or GPU; the publication validation covers CPU workers and the
data, not a new full CUDA census.

The GPU worker accepts the same line protocol as the dense profile worker, with
prime domain 1000 ≤ p < 2^32 and at most 128 primes per batch. Original arguments
are in `GPU_CONFIG.json`. For the final campaign capacity policy, replace
`--capacity` with 4194304 and `--epoch` with 2048. A recorded capacity refusal
required replaying the exact same entire batch with epoch halved, down to 32,
within a single 30-second total deadline. Other failures were not retried and
partial batches were not accepted. Earlier blocks used the original policy;
no policy changes altered the mathematical outputs.

The source directory includes inherited project snapshots under nested
`vendor/research/` and `baseline/` paths. These names identify this project's
historical engine layers, not an installed third-party software distribution.
Machine-specific launch, sensor and thermal-override scripts are not part of
this publication. A new large census therefore requires a separately reviewed
execution supervisor and explicit resource planning; this repository does not
silently schedule or restart one.

## Exact theoretical certificates

The three Python verification entry points in the main README need only the
Python standard library. They verify finite character matrices, exact square
and fourth-power identities in number fields, elementary Galois certificates,
norm spans, density formulae, code counts and bounded residual investigations.
Assertions must be enabled: use ordinary `python3 -B`, never `python3 -O`.

For independent regeneration of local maximum distributions:

```bash
g++ -O3 -std=c++17 reproducibility/proofs/local_maximum.cpp -o build/local_maximum
build/local_maximum reproducibility/evidence/new/local_relations_input.txt build/local_maximum.csv
cmp build/local_maximum.csv reproducibility/evidence/new/LOCAL_MAXIMUM_CUBIC.csv
```

PARI/GP scripts and their mathematical inputs are included. PARI is an optional
external dependency for regenerating number-field certificates, not a
prerequisite for exact certificate verification. Historical command transcripts
and installed PARI binaries are not included.

## Validation limits

Exact arithmetic audit is not a proof-assistant formalization. Byte commitments
support archive integrity, not ownership or novelty claims. Small independent
recomputations and periodic audits support the complete census, but do not
constitute a new complete independent computation of every prime. Post-census
models do not reopen closed primary tests or establish universal laws.
