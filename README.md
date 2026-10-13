# Rabung atlas below 2^32

Research code, data and exact certificates accompanying four papers by
**Brice Pouly** on power-residue colourings, Rabung lower bounds and arithmetic
fluctuations. The repository contains a complete census below $2^{32}$, a
separate sampled ternary supplement, CPU and CUDA sources, and reproducibility
packages for the finite computations reported in the papers.

## Papers

| Article | PDF | LaTeX source | Reproducibility |
| --- | --- | --- | --- |
| Succinct Rabung certificates and recurrence-closed lower bounds for van der Waerden numbers | [v2.0](paper/RABUNG_BOUNDS_EN.pdf) | [Source](paper/RABUNG_BOUNDS_EN.tex) | [Certificates and recurrence closure](reproducibility/bounds_update/README.md) |
| Certified rational-anchor filtering for exact GPU scans of power-residue runs | [v1.0](paper/GPU_METHOD_EN.pdf) | [Source](paper/GPU_METHOD_EN.tex) | [GPU method and measurements](reproducibility/gpu_method/README.md) |
| Arithmetic structure in a complete atlas of power-residue runs below 2^32 | [v1.0](paper/PAPER_EN.pdf) | [Source](paper/PAPER_EN.tex) | [Atlas reproduction guide](docs/REPRODUCTION.md) |
| Cubic-residue quadruplets: elliptic formulas and arithmetic fluctuations | [v1.0](paper/CUBIC_QUADRUPLETS_EN.pdf) | [Source](paper/CUBIC_QUADRUPLETS_EN.tex) | [Elliptic counts and moments](reproducibility/cubic_quadruplets/README.md) |

All manuscripts use $W(r,k)$ with the number of colours first. They acknowledge
AI assistance. The certificates paper is published as
[v2.0, DOI 10.5281/zenodo.23337824](https://doi.org/10.5281/zenodo.23337824).
Its [first edition, v1.0.0](https://doi.org/10.5281/zenodo.22980965), remains
available for the original evidence and historical comparisons. The GPU
method paper is published as
[v1.0, DOI 10.5281/zenodo.23338707](https://doi.org/10.5281/zenodo.23338707).
The atlas paper is published as
[v1.0, DOI 10.5281/zenodo.23339668](https://doi.org/10.5281/zenodo.23339668).

## GPU hardware

**The GPU program was designed specifically for NVIDIA SM120 (`sm_120`,
compute capability 12.0) and has only been used on a single NVIDIA GeForce
RTX 5090.** All GPU census runs, qualification checks and performance
measurements reported in this repository were obtained on that card.
Other GPU architectures and models have not been validated.

The CUDA build targets `compute_120` / `sm_120`. Its configuration, compiler
options and source hashes are recorded in
[GPU_CONFIG.json](src/census/GPU_CONFIG.json). An optional build requires a
C++17 host compiler and a CUDA toolkit with SM120 support:

```bash
python3 -B scripts/build_workers.py --gpu
```

This command compiles the workers. GPU execution is explicit; the default
verification commands and CI do not start GPU calculations. The reported
performance figures apply to the RTX 5090 and the configurations documented in
the GPU paper and its reproducibility package.

## Complete atlas

The census covers every prime $p < 2^{32}$ and every admissible order
$r = 2,\ldots,9$ with $r \mid p-1$.

| Quantity | Exact value |
| --- | ---: |
| Primes, including $p = 2$ | 203,280,221 |
| Admissible power-residue profiles | 677,591,249 |
| Blocks, each containing at most 4,096 primes | 49,901 |
| Sum of exact quadratic $A_4$ counts over odd primes | 53,206,216,723,199,456 |
| Scientific pack files | 70 |

Each profile stores the longest monochromatic run $M$, the initial residue run
$L$, and a canonical reflected witness. The Rabung threshold $K$ combines the
interior and boundary conditions. The quadratic $A_4$ statistic counts
monochromatic four-windows on $1,\ldots,p-1$; $p = 2$ has a zero placeholder.
Definitions are given in [METHODS.md](docs/METHODS.md).

### Quick start

Run these commands from the repository root with Python 3.10 or newer.
The sample checks use only the Python standard library and require no
download or GPU:

```bash
python3 -B scripts/verify_atlas.py
python3 -B scripts/query.py 1009
python3 -B -m unittest discover -s tests
```

The two included example blocks contain 6,409 primes near the two ends of the
domain. They are non-contiguous samples; completeness applies to the full
atlas described above.

### Download the full dataset

The [data release `atlas32-v1.0.0`](https://github.com/Vulkin-prog/rabung-atlas32/releases/tag/atlas32-v1.0.0)
contains **72 assets totaling 4,852,996,533 bytes**: the 70 scientific packs
and two SQLite indices. The
[manifest](data/atlas/manifest.json) records each asset's destination path,
size and SHA256.

Download the assets from the release page, or use the checksum-verifying
downloader with [GitHub CLI](https://cli.github.com/):

```bash
gh auth login
python3 -B scripts/download_data.py
python3 -B scripts/verify_atlas.py --data data/atlas
python3 -B scripts/query.py 2940141229 --data data/atlas
```

The downloader reuses files whose hashes already match and maintains at
least 5 GiB of free disk space. Full verification checks asset hashes,
interval coverage, scientific hash chains, indexed counts, aggregates and
endpoint frames. It verifies stored data rather than recomputing every prime.

The raw packs preserve the original bytes. The SQLite indices are curated
publication copies that preserve mathematical data and indexed commitments
while omitting execution metadata. See [DATA_FORMAT.md](docs/DATA_FORMAT.md)
and [CURATION.json](docs/CURATION.json) for the formats and provenance.

## Reproduce the computations

### CPU regeneration and atlas analysis

The CPU build uses `g++` with C++17 support:

```bash
python3 -B scripts/build_workers.py
python3 -B scripts/recompute_small.py --count 256
```

The profile and exact quadratic $A_4$ workers regenerate the first 256 primes
and compare the results with the reference data.

With the complete dataset, zlib and OpenSSL development headers, regenerate
the descriptive aggregates by streaming every frame:

```bash
python3 -B scripts/reproduce_analysis.py --data data/atlas
```

Integer fields are compared exactly; floating summaries use a relative
tolerance of $10^{-10}$. Machine-dependent timings are excluded from that
comparison. Worker protocols and requirements are documented in
[REPRODUCTION.md](docs/REPRODUCTION.md).

### Article checks

Each article has a dedicated verification route:

```bash
python3 -B scripts/verify_bounds_update.py
python3 -B scripts/verify_gpu_method.py
python3 -B scripts/verify_editorial.py --scope atlas --catalogue
python3 -B scripts/verify_editorial.py --scope cubic
python3 -B scripts/cubic_quadruplets.py 1000033
```

The GPU-method verifier and the atlas catalogue replay require a C++17
compiler. These checks use retained measurements, exact CPU references and
finite certificates; they do not launch a GPU campaign.

The atlas includes certified character classes, Kummer and Chebotarev
densities, norm dependencies and local code laws. Useful entry points are the
[record catalogue for $K = 3,\ldots,64$](data/derived/ATLAS_CATALOGUE_K3_K64.csv),
the [order-four density table](data/derived/ORDER_FOUR_DENSITIES.csv), and the
[exact proof checks](reproducibility/proofs/).

### Ternary supplement

The [ternary-million dataset](data/ternary-million/README.md) contains
**1,000,000 sampled profiles** in $[2^{33},2^{34})$, split into four fixed
strata. Its 248 blocks include exact $M$, $L$ and $K$, run counts for thresholds
18 through 27, and **1,235,473 canonical maximal-run positions of length at
least 21**. There are 33,257,058 canonical maximal runs of length at least 18;
positions of lengths 18 through 20 were not retained for the full sample.

The supplement is included in Git and can be checked without downloading the
complete atlas:

```bash
python3 -B scripts/verify_ternary.py
python3 -B scripts/query_ternary.py 8611114591
```

It is a stratified sample separate from the complete census below $2^{32}$.
The [ternary diagnostics](reproducibility/ternary_diagnostics/) describe the
fixed-coordinate catalogue and the comparator using the same full-run mask.

### Scientific scope

The papers distinguish the complete finite census, the sampled ternary
supplement, exact algebraic certificates and descriptive statistical models.
Models selected after observing the data are exploratory. The remaining
order-four witness is open; its bounded negative evidence is not an
impossibility theorem.

The cubic four-window statistic is distinct from the atlas's quadratic
$A_4$ column. Its finite verification checks identities, elliptic trace
certificates and moment bookkeeping. Four moments do not determine a complete
limiting distribution or a law for critical maxima, and finite programs do
not prove the analytic theorems cited in the papers.

## Build article archives and PDFs

From a Git checkout, build one ZIP per article and verify its contents after
extraction. The replay requires Python 3.10 or newer and a C++17 compiler:

```bash
python3 -B scripts/build_publication_archives.py --verify --source-commit "$(git rev-parse HEAD)"
```

The archives are written to `build/publication/`. Each contains the article's
PDF and standalone LaTeX source, an abstract with inline LaTeX, the relevant
code and evidence, rights notices, and a SHA256 manifest. The certificates
archive also contains a concise v2 version note. To include the complete
unchanged first-edition deposit, add
`--legacy-archive /path/to/vdw_rabung_publication_package.zip`.
The full atlas release remains a separate data download.

To rebuild all four PDFs with a LaTeX installation, run the following in a
Bash-compatible shell:

```bash
mkdir -p build/paper
for source in paper/RABUNG_BOUNDS_EN.tex paper/GPU_METHOD_EN.tex paper/PAPER_EN.tex paper/CUBIC_QUADRUPLETS_EN.tex; do
  for pass in 1 2; do
    pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error -output-directory=build/paper "$source"
  done
done
```

## Citation and rights

When using these results, cite the relevant paper and identify the
`atlas32-v1.0.0` data release where applicable. Repository citation metadata
is provided in [CITATION.cff](CITATION.cff).

Copyright 2026 Brice Pouly for project-authored material, subject to
pre-existing rights in individual files. All rights reserved. No general
open-source or open-content licence is granted for the code, manuscripts or
data. See [RIGHTS.md](RIGHTS.md) and [THIRD_PARTY.md](docs/THIRD_PARTY.md) for
rights and dependency notices.
