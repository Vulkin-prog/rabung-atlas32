# Four articles, including the published certificates paper

This revision follows the accepted organization: one published article, revised
as a new draft, and three companion articles. It does not create a fifth paper.

| Article | Version in this PR | Main question and proof ownership |
| --- | --- | --- |
| Succinct Rabung certificates and recurrence-closed lower bounds for van der Waerden numbers | Revised draft v1.1; published v1.0.0 remains at DOI 10.5281/zenodo.22980965 | Which succinct triples satisfy Rabung's criterion, and which bounds follow after closing an explicit source corpus? Owns the full boundary equivalence, pointwise checks, fixed-prime monotonicity, exact-modulus products and bound comparisons. |
| Certified rational-anchor filtering for exact GPU scans of power-residue runs | v0.3 | How can a scanner preserve exact maxima, prefixes and canonical witnesses while avoiding a dense expensive character evaluation? Owns the integer cover, anchor selection, conservative rejection, completion, arithmetic and controlled timings. |
| Arithmetic structure in a complete atlas of power-residue runs below 2^32 | v1.5 | What does the complete finite population reveal, and which selected structures have exact algebraic explanations? Keeps census, Kummer densities, norm dependencies, local codes and ternary diagnostics together. |
| Cubic-residue quadruplets: elliptic formulas and arithmetic fluctuations | v1.3 | How can a consecutive cubic four-window count be computed exactly, and what are its first four moments over primes? Owns the assembled Greene–Lennon count, boundary corrections, geometric proof, six curves, CM twists and analytic moment argument. |

The atlas retains a short Rabung specialization for readability and cites the
published full boundary proof. The cubic proofs use neither the census nor a
GPU speed measurement. Short-pattern counts and local density statements do
not imply a logarithmic-length maximum law. Census times are not inserted into
the methods paper's two controlled performance comparisons.

## Conventions and references

All four current manuscripts put colours first: `W(r,k)`. A strict stored
integer `B` means `W(r,k) > B`, equivalently `W(r,k) >= B+1`. Historical CSV
fields, seals, registries and closed tests keep their original names and bytes.
The atlas field `W_length_first_lower_bound` remains a historical schema name;
its numerical entries are weak integer bounds, one above construction length.

Article bibliographies cite published work, versioned public preprints,
scientific software/data artifacts and the companion manuscripts intended for
publication. Supplied reviews, private working notes and conversation histories
are not scientific references. Existing internal review records remain work
records and are not redistributed as received PDFs.

## Changes in this revision

The GPU revision counts in-table scalar completion steps explicitly as
`C_tab`, distinguishes launched grids from resident blocks, separates 1024-thread
generators from 128-thread completion blocks, uses "worker slot" for concurrent
prime work and "warp lane" inside a warp, and follows one small window through
the full cover/selection/completion chain. The measured algorithms and timings
are preserved. The selected reproduction package curates host-path and disk-only
metadata explicitly; original file identities remain recorded.

The bounds revision keeps the published thirteen claims and their gate intact.
Its new register has six primitive certificates giving fourteen direct rows,
supported by nine complete CPU profile rechecks. Three extra checks bind the
finite-atlas ternary maxima; the stronger outside-atlas ternary witnesses are
used for the current bound table. The finite recurrence scope extends the
published four-colour corpus to six colours, with no fresh global priority
audit. The unresolved Landman–Robertson statement remains a separate comparator.

The original published edition is not overwritten. This PR provides a revised
draft in the atlas repository; it does not update the original deposit, change
repository visibility, submit an article or merge the PR. A public immutable
artifact and specialist review of the cubic analytic interfaces remain needed
for submission. The four-moment argument is not formally machine-verified.
