# Version identifiers and reproduction map

The existing draft PR is [#1](https://github.com/Vulkin-prog/rabung-atlas32/pull/1).
The reviewed editorial baseline is commit
`59a9f9c501ab6781996222470049940ddef2e691`; the sealed-census manuscript
baseline is `452907920dc623d4d1e6ae7675428fe02ed840c4`.
The current four-paper revision is dated 9 October 2026: certificates v1.2,
GPU methods v0.4, atlas v1.6 and cubic companion v1.4. The published certificates
edition v1.0.0 remains at DOI 10.5281/zenodo.22980965 and is not overwritten.
The source-policy correction follows commit
`94d586095c1b12edf185c31f20f6f7c72c786c4a` (atlas v1.3 / companion v1.1).
Private reviews and internal working notes are excluded from the manuscript
bibliographies and prose references. The scientific reproduction artifacts
are identified directly in the articles.

The publication sequence is revised certificates, GPU method, atlas, cubic
quadruplets. `PUBLICATION_SEQUENCE.json` distinguishes current public editions
from planned insertions. `scripts/verify_publication_references.py` checks the
actual series references and exact edition locators. Its offline check does
not establish public availability or publish anything. See
[the per-stage plan](SUCCESSIVE_PUBLICATION_PLAN_FR.md).
The current compilation, bibliography and scoped reproduction checks are in
`SUCCESSIVE_PUBLICATION_VALIDATION.json`. `FOUR_ARTICLE_VALIDATION.json` retains
the scientific checks performed on the preceding four-article revision.

For a checked-out revision, record `git rev-parse HEAD` and use that full SHA
in a GitHub `blob/<SHA>/<path>` permalink. The commit pins the complete tree;
the editorial SHA256 manifest additionally identifies the new and edited
files. A branch URL or a PR number alone is not an immutable version identifier.
The commit containing this document is the version anchor for its revised
manuscripts and additional certificates; the baseline SHA above identifies
the versions supplied to the reviews.

| Result package | Versioned location / identity | Minimal reproduction |
| --- | --- | --- |
| Revised atlas | `paper/PAPER_EN.tex` and `.pdf`; `docs/EDITORIAL_SHA256.json` | Two pdflatex passes with `-no-shell-escape`. |
| Revised cubic companion | `paper/CUBIC_QUADRUPLETS_EN.tex` and `.pdf`; same editorial manifest | Two pdflatex passes; `python3 -B scripts/cubic_quadruplets.py 1000033`. |
| Revised certificates paper | `paper/RABUNG_BOUNDS_EN.tex` and `.pdf`; published edition separately pinned | Two pdflatex passes; `python3 -B scripts/verify_bounds_update.py`. |
| GPU methods paper | `paper/GPU_METHOD_EN.tex` and `.pdf`; source edition f00cdcb95f12d7578bdb37a87db97ff233c86a22 | Two pdflatex passes; `python3 -B scripts/verify_gpu_method.py`. |
| Nine pointwise CPU rechecks | `reproducibility/bounds_update/results/`; source and output bound in each record | Default verification checks retained executions and literal local powers; optional large full-period reruns use the separate CPU checker. |
| Published input preservation | `reproducibility/bounds_update/PUBLISHED_IMPORTS.json`; nine original Git blob identities | Default bounds verifier also reproduces the original closure exactly. |
| GPU selected packet | `reproducibility/gpu_method/PROVENANCE.json` and `SOURCE_MAP.json`; 34 selected imports, not all 272 original imports | Benchmark regeneration and finite host source bridges; no CUDA run or new qualification. |
| Frozen scientific scanner | `src/census/`; census baseline commit above; scientific data seal | `python3 -B scripts/build_workers.py`, then `python3 -B scripts/recompute_small.py --count 256`. |
| Complete atlas data | Release `atlas32-v1.0.0`; asset lengths and hashes in `data/atlas/manifest.json` | `python3 -B scripts/verify_atlas.py` checks the committed samples; full-release checking requires the separately downloaded assets. |
| Original mathematical certification | `reproducibility/anteriorities/` and `reproducibility/evidence/new/` | `python3 -B reproducibility/proofs/verify_synthesis.py`; `python3 -B reproducibility/proofs/verify_report4.py`. |
| Residual irreducibility addition | `RESIDUAL_IRREDUCIBILITY.json` in the evidence directory; 100 model bindings and 200 finite witnesses | `python3 -B reproducibility/proofs/residual_irreducibility.py`. `--generate` regenerates the small witness file. |
| Greene–Lennon exposition check | `scripts/verify_hypergeometric.py`; archived geometric verifier kept intact | `python3 -B scripts/verify_hypergeometric.py`. |
| Elliptic count and moments | `reproducibility/cubic_quadruplets/`; copied-source identities in `docs/EDITORIAL_PROVENANCE.json` | `python3 -B scripts/verify_editorial.py --scope cubic`. |
| Fixed million ternary population | `data/ternary-million/`; original commitments and publication provenance | `python3 -B scripts/verify_ternary.py`. |
| Same-filter diagnostics and retained coordinates | `reproducibility/ternary_diagnostics/`; complete ordinary catalogue retained | `python3 -B scripts/verify_editorial.py --scope atlas --catalogue` with a C++17 compiler. |

The complete original atlas manifest SHA256 is
`762dc0aadc9beaa8a6bf2cb1a82604ad0c3dfcbb9f67c1a852ca7a2bb4468c20`.
It is separate from the editorial manifest. No revised mathematical claim
is presented as a new full census. None of these default checks starts a GPU
campaign. Finite checks do not formally verify the automorphic input theorems.

Repository access currently serves authorized reviewers. A public immutable
archive, persistent identifier, confirmed redistribution rights and specialist
human review remain publication tasks. No public availability or first public
disclosure date is inferred from a private commit, an internal note date, or a
hash. This revision does not change visibility, publish an archive, or submit
any paper to a journal. Supplied external review PDFs are not redistributed.
