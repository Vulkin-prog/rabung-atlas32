# Historical ternary long-run collector

The `fast_source/` tree consists of the 58 exact scientific source files listed
by the historical `FAST_BUILD.json`. Their hashes, compile definitions and
architecture are preserved in [SOURCE_IDENTITY.json](SOURCE_IDENTITY.json).
They produced the fixed ternary-million supplement in [2^33,2^34), separate
from the complete atlas below 2^32. No binary or CUDA installation is bundled.

The fixed minimum is 18 for the million-profile counts. The collector retains
all complete canonical maximal runs above that minimum, including d=1, disables
early stopping, completes every surviving sparse-cover anchor and deduplicates
exact intervals with duplicate agreement. The historical complete small-prime
and 512 wide-prime comparisons are in
[FAST_QUALIFICATION.json](../../data/ternary-million/analysis/FAST_QUALIFICATION.json).
Scientific comparisons are unchanged; execution timings are omitted in that
publication copy. The original qualification hash remains a historical anchor
in the unchanged count protocol.

From the repository root, `python3 -B scripts/build_ternary_worker.py` explicitly
compiles the source for CUDA SM120. It never runs the worker. This optional
compilation has not been requalified for the publication addition. The
standard-library stored-data verifier needs neither compilation nor a GPU.

The worker takes `--minimum 18`, then line-oriented input containing a sequence
number and at most 128 admissible prime integers. It emits full long-run
catalogues and a batch completion line. Capture overflow or incomplete
completion refuses with an error; output must not be accepted as a partial
catalogue. See the source for the complete protocol. No automatic launcher or
hardware-parameter adjustment is provided.

Unchanged scientific functions `vu` and `process` from the historical
`count_million.py`, and `graph` from `critical_graphs.py`, are extracted in the
repository's `scripts/` directory. [CODE_EXTRACTION.json](CODE_EXTRACTION.json)
records their source identities and extraction operations. The original host
supervisors are omitted. Other source provenance and rights follow the main
repository's curation and third-party documents.
