# Rabung atlas below 2^32

Reproducible data, source code, exact algebraic certificates and the manuscript
**Arithmetic structure in a complete atlas of power-residue runs below 2^32**,
by Brice Pouly. The manuscript was prepared with AI assistance; it is a research
manuscript, not a peer-reviewed publication. See [the paper](paper/PAPER_EN.pdf)
and [its standalone LaTeX source](paper/PAPER_EN.tex).

The complete census covers every prime in **[1, 4,294,967,296)** and every
admissible order **r = 2, …, 9** (r divides p−1).

| Quantity | Exact value |
| --- | ---: |
| Primes, including p = 2 | 203,280,221 |
| Admissible power-residue profiles | 677,591,249 |
| Committed blocks, each at most 4,096 primes | 49,901 |
| Sum of exact quadratic A4 counts over odd primes | 53,206,216,723,199,456 |
| Scientific pack files | 70 |

Each profile stores the longest monochromatic run **M**, initial residue run
**L**, and canonical reflected witness. The Rabung threshold **K** incorporates
both the interior and boundary conditions. Quadratic **A4** counts are exact
counts of monochromatic four-windows on 1,…,p−1; p = 2 has a zero placeholder.
Definitions and the distinction between criterion (a) and (a)+(b) are in
[METHODS.md](docs/METHODS.md).

## Start without downloading the full atlas

Python 3.10 or newer is sufficient for the portable checks. No Python package,
network connection, PARI installation, or GPU is required.

```bash
python3 -B scripts/verify_atlas.py
python3 -B scripts/query.py 1009
python3 -B -m unittest discover -s tests
python3 -B reproducibility/proofs/verify_synthesis.py
python3 -B reproducibility/proofs/verify_report4.py
python3 -B reproducibility/anteriorities/theorie_densites_residuels_20261004/preuves/verify_results.py
```

The two committed original blocks contain 6,409 primes near the two ends of the
domain. **They are non-contiguous examples, not a smaller complete census.**
The theory verification checks 445 fourth-power identities and 808,064 finite
character values, norm dependencies, local code laws and bounded negative
certificates using exact arithmetic. Typical time is about one minute on an
ordinary desktop; hardware affects timing.

## Complete atlas data

The [pinned release `atlas32-v1.0.0`](https://github.com/Vulkin-prog/rabung-atlas32/releases/tag/atlas32-v1.0.0)
contains 72 assets, totaling **4,852,996,533 bytes**: the 70 original scientific
packs and two publication-derived SQLite indices. The asset hashes, lengths,
destination paths and original seal identity are in
[data/atlas/manifest.json](data/atlas/manifest.json).

While the repository is private, downloading requires authorized GitHub access.
Authenticate with [GitHub CLI](https://cli.github.com/) using `gh auth login`,
then run:

```bash
python3 -B scripts/download_data.py
python3 -B scripts/verify_atlas.py --data data/atlas
python3 -B scripts/query.py 2940141229 --data data/atlas
```

The downloader checks each SHA256, resumes by reusing verified files and keeps
at least 5 GiB of free disk space. The complete verification checks every asset
hash, interval coverage, both scientific chains, indexed counts and aggregates,
and endpoint frames. This is a stored-data check; it does not recompute every
prime. `--indices-only` on the verifier omits the whole-pack hash checks.
Downloading only the indices is possible with `download_data.py --indices-only`,
but querying or verifying endpoint frames also requires their packs.

Raw packs remain byte-identical to the sealed source. The SQLite indices are
new, explicitly documented publication copies: mathematical data and indexed
commitments are unchanged, while execution metadata is omitted. They must not
be described as byte-identical copies of the original indices. See
[DATA_FORMAT.md](docs/DATA_FORMAT.md) and [CURATION.json](docs/CURATION.json).
Release assets keep the multi-gigabyte dataset outside Git history, in accordance
with [GitHub's large-file guidance](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github).

## Recompute results

With a C++17 compiler, bounded CPU regeneration is:

```bash
python3 -B scripts/build_workers.py
python3 -B scripts/recompute_small.py --count 256
```

The unchanged dense profile worker and exact A4 worker regenerate the first
256 primes and compare every result with the original committed data. This is
a bounded check, not a GPU benchmark or a new campaign.

To regenerate **all descriptive aggregates** from the complete atlas, install
the zlib and OpenSSL development headers and run:

```bash
python3 -B scripts/reproduce_analysis.py --data data/atlas
```

This streams every frame, verifies frame hashes, checks structural inequalities
and the elliptic Hasse bound, and compares the numerical aggregates with the
committed analysis. Machine-dependent timings are excluded from comparison;
integer fields are compared exactly and floating summaries with relative
tolerance 10^−10. Output is written under ignored `build/reanalysis/`.

The qualified CUDA SM120 engine sources, macros and amended capacity policy are
included for researchers wishing to rebuild the original workers. Compilation
is explicit: `python3 -B scripts/build_workers.py --gpu` requires a compatible
NVIDIA toolchain. No GPU calculation is started by any default command or CI.
Historical binaries, machine-specific monitors and launch scripts are omitted.
See [REPRODUCTION.md](docs/REPRODUCTION.md) for worker protocols and limitations.

## Fixed million ternary profiles and long-run positions

The [ternary-million supplement](data/ternary-million/) is committed directly
in Git: **1,000,000 sampled profiles** in **[2^33, 2^34)**, split into four
fixed strata, with exact M, L and K, ordinary and canonical counts at
thresholds 18..27, and **1,235,473 canonical maximal-run positions at length
at least 21**. It contains 33,257,058 canonical maximal runs at length at least
18. Positions at lengths 18..20 were not retained for the full million.

This existing descriptive supplement is separate from the complete atlas
below 2^32. The 248 count blocks, four prime lists, parent profiles, doubling
graph witnesses and collector sources are included with original hash anchors
and publication provenance. No extra download or GPU calculation is required.

```bash
python3 -B scripts/verify_ternary.py
python3 -B scripts/query_ternary.py 8611114591
```

The verifier checks the entire stored million and all long positions, repeats
1,189 exact full-run CPU checks and reconstructs all 398 doubling edges.
The [supplement documentation](data/ternary-million/README.md) explains
canonical reflection, colour labels, binary formats and CSV exports. This
sample includes the earlier 512 diagnostic primes; it is not an independent
validation and does not alter any closed primary test.

The atlas revision also includes the fixed-coordinate catalogue and the exact
iid comparator using the same full-run mask. These diagnostics and their
continuation budgets are in [reproducibility/ternary_diagnostics/](reproducibility/ternary_diagnostics/).
They remain descriptive analyses of this fixed sample.

## Cubic quadruplets and elliptic fluctuations

The separate research manuscript
[Cubic-residue quadruplets: elliptic formulas and arithmetic fluctuations](paper/CUBIC_QUADRUPLETS_EN.pdf)
has a [standalone LaTeX source](paper/CUBIC_QUADRUPLETS_EN.tex). It assembles
the Greene–Lennon four-position identity into an exact consecutive cubic count,
including the boundary corrections, with six-curve trace certificates and
the first four moments over primes. The direct geometric proof cites the
preceding family and isogeny of Freeman–Satoh and Guillevic–Vergnaud. The two
arithmetic regimes depend on whether three is a cube. This cubic statistic
is distinct from the historical quadratic A4 column of the atlas.

From the repository root:

    python3 -B scripts/cubic_quadruplets.py 1000033
    python3 -B reproducibility/proofs/residual_irreducibility.py
    python3 -B scripts/verify_editorial.py
    python3 -B scripts/verify_editorial.py --catalogue

The count command validates stored elliptic certificates. The irreducibility
command binds all 100 residual models to their character certificates and
replays 200 finite witnesses over Q and Q(i). The editorial command additionally
checks the hypergeometric derivation against 2,360 literal four-position sums,
then replays the curated elliptic, moment and ternary diagnostics.
Its `--catalogue` option requires a C++17 compiler and
reclassifies all retained ordinary long-run positions. No command rescans
the full million character words.

The manuscript gives an argument without GRH. Independent review of the
analytic passage is still required before publication; see
[ELLIPTIC_REVIEW.md](docs/ELLIPTIC_REVIEW.md). Four moments do not establish
a complete limiting distribution or a critical maximum law. The editorial
scope and separate provenance manifests are in
[EDITORIAL_INTEGRATION.md](docs/EDITORIAL_INTEGRATION.md).
The atlas v1.4 and companion v1.2 incorporate the two reports received on
8 October 2026. The [response to the reviews](docs/REVIEW_RESPONSE_2026-10-08_FR.md)
records the corrections and remaining editorial tasks. The
[reproduction map](docs/PUBLICATION_REPRODUCIBILITY.md) identifies the frozen
data and scanner separately from the manuscript revisions.
The article bibliographies exclude private reviews and internal working notes;
scientific statements are supported by the supplied proofs and certificates.

## Manuscript and theory

The paper combines the atlas analysis with exact character-class certification,
Kummer and Chebotarev densities, quadratic propagation, a length-53 cubic code,
norm-generated dependencies and a dihedral quartic. Eight order-four density
laws are certified. The remaining order-four witness is still open; its bounded
negative evidence is included without claiming an impossibility theorem.

- [Record catalogue, K = 3,…,64](data/derived/ATLAS_CATALOGUE_K3_K64.csv)
- [Order-four density table](data/derived/ORDER_FOUR_DENSITIES.csv)
- [Descriptive analysis and exploration](analysis/)
- [New exact certificates](reproducibility/evidence/new/)
- [Earlier mathematical certificates](reproducibility/anteriorities/)

The theoretical models were selected after observing the atlas. Old primary
tests remain closed; these examples are not prospective statistical tests,
universal growth laws, or claims of established bibliographic priority.

## Visibility, rights and provenance

This repository was created **private**. Making it public is a separate owner
decision. Host logs, sensor histories, local absolute paths, credentials,
installed third-party software and externally received review reports are not
included. Provenance records distinguish exact copies from redacted metadata
and editorial packaging changes; original research archives are preserved.

No new open-source or open-content licence has been granted by depositing these
files. Rights are reserved pending an explicit owner choice. Consult
[RIGHTS.md](RIGHTS.md), [THIRD_PARTY.md](docs/THIRD_PARTY.md),
[PUBLICATION_AUDIT.json](docs/PUBLICATION_AUDIT.json) and [CITATION.cff](CITATION.cff)
before a future public release. Automated checks support this curation; they do
not establish legal ownership or mathematical novelty.

To compile both manuscripts independently:

```bash
mkdir -p build/paper
for source in paper/PAPER_EN.tex paper/CUBIC_QUADRUPLETS_EN.tex; do
  for pass in 1 2; do
    pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/paper "$source"
  done
done
```
