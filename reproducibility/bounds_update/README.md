# Separate addendum to the published Rabung certificates paper

The revised draft is `paper/RABUNG_BOUNDS_EN.tex` (v1.1, 9 October 2026).
The published v1.0.0 remains at <https://doi.org/10.5281/zenodo.22980965>.
Its original repository snapshot is
`Vulkin-prog/vdw-rabung-bounds@704b4125b274bab3df45454adc9c421f8068b8a1`.
`PUBLISHED_IMPORTS.json` binds nine selected original files to their Git blob
identities, lengths and SHA256 values. The thirteen-claim gate, original source
program and original closure are copied unchanged, not rewritten or reopened.
The complete original publication artifact remains at its original location.

`results/` retains nine new complete CPU executions with exact source/binary
identities, commands and literal output. Six concern atlas target primes;
three concern stronger ternary constructions outside `2^32`. A full scan visits
every nonzero residue in ordinary order, with no mirror or sparse-cover
assumption. This is author-controlled separate-source verification, not an
external replication. No new GPU census or benchmark was run.

The byte-table checker accepts primes below `2^32` and orders 2 through 9.
The two-bit checker accepts orders 2 or 3, primes below `2^35`, target length
at most 256, and an explicit product-capacity guard. It uses exact 128-bit
products for powering and checked 64-bit products for its multiplicative walk.
These are CPU domains, not an extension of the GPU worker domain. At the
largest retained prime, two-bit storage needs approximately 4.3 GB.

From the repository root:

```bash
mkdir -p build
python3 -B scripts/verify_bounds_update.py
g++ -O3 -std=c++17 -Wall -Wextra -Wconversion -o build/independent_profile reproducibility/bounds_update/independent_profile.cpp
python3 -B reproducibility/bounds_update/check_small.py build/independent_profile
g++ -O3 -std=c++17 -Wall -Wextra -Wconversion -o build/two_bit_profile reproducibility/bounds_update/two_bit_profile.cpp
python3 -B reproducibility/bounds_update/check_small.py build/two_bit_profile --max-order 3
```

The default verifier checks the retained records, primality and root orders,
local run/prefix powers, atlas binding, all strict/weak integer conversions,
the unchanged old closure and the new closure. It does not re-establish a
global maximum from local powers: that upper bound uses the separately retained
complete CPU execution. Small controls compare the checker with literal power
colours and explicit arithmetic progressions, including refusal cases.

`--generate` regenerates `CLAIMS_UPDATE.json`, `CLOSURE_UPDATE.json` and the two
LaTeX tables. The updated admitted scope is lengths 12 through 28 and up to six
colours; the reported three-colour table starts at 17. Published inputs retain
their source attribution. The stronger Landman–Robertson published statement
keeps its unresolved proof status and is displayed separately. No new global
literature priority claim is made.

To repeat a complete large check, run the appropriate rebuilt checker with
`p r k`. This allocates gigabytes and scans the entire period; it is an optional
pointwise reproduction, separate from the bounded CI and default verification.
The succinct triple, Rabung proof and CPU source suffice without any timing
claim, law of maxima or full search campaign.
