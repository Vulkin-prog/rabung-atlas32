# Cubic quadruplets: curated certificates

The [companion manuscript](../../paper/CUBIC_QUADRUPLETS_EN.tex) gives the
four-position elliptic reduction, exact consecutive count and four vertical
moments. These directories retain the scientific inputs, algorithms and
finite results from the two selected elliptic research archives. Archive
identities and source hashes are in
[EDITORIAL_PROVENANCE.json](../../docs/EDITORIAL_PROVENANCE.json).

From the repository root:

    python3 -B reproducibility/cubic_quadruplets/elliptic/src/verify.py
    python3 -B reproducibility/cubic_quadruplets/moments/src/verify.py
    python3 -B scripts/cubic_quadruplets.py 1000033

These commands use standard-library Python. The count oracle validates six
curve traces with point-order and Hasse certificates, checks the CM relations,
and evaluates the exact integer count. It ships certificates for the three
large examples in the manuscript. For another admissible prime below 10^10,
explicit --generate requests a bounded certificate search; failure to find
a witness is reported rather than guessed. It does not enumerate a full
character word. Optional historical generators are retained, but are not run
by the default verification command.

Finite identity checks and case enumeration do not prove automorphy or a
prime number theorem. The independent analytic review status is described in
[ELLIPTIC_REVIEW.md](../../docs/ELLIPTIC_REVIEW.md). These certificates do not
establish a full limiting distribution or a law for critical maxima.
