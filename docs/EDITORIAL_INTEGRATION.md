# Editorial integration, 8 October 2026

The atlas manuscript is updated in place. Its already integrated Kummer
densities, norm dependencies, order-four models and dihedral quartic are
retained, rather than announced again as new results. Its current unresolved
catalogue has 22 binary profiles and one order-four profile.

| Destination | Added contribution |
| --- | --- |
| Atlas manuscript | Fixed million ternary population; ordinary/canonical accounting; coordinate catalogue; continuation diagnostics; exact iid comparator with the same mask; precise scope of the critical conjecture. |
| Cubic quadruplets companion | Exact consecutive count assembled from the Greene–Lennon identity; direct proof with the preceding Freeman–Satoh and Guillevic–Vergnaud geometry; six-curve certificates; two CM fluctuation regimes; first four moments over primes; signed support-four appendix. |
| Reproducibility packages | Curated elliptic and moment certificates, fixed-coordinate diagnostics, same-filter arithmetic and portable verification commands. |
| Review record | Internal review plus two supplied reports; attribution and local corrections; explicit analytic interfaces; specialist human review still needed before submission. |

The atlas is now v1.3 and the companion v1.1. See the
[response to the two reports](REVIEW_RESPONSE_2026-10-08_FR.md). The historical
residual verifier is extended to check polynomial binding and 200 new
irreducibility witnesses; its earlier character matrices remain unchanged.
A separate small-prime check verifies the hypergeometric exposition against
literal four-position sums and direct elliptic point counts.

See [ELLIPTIC_REVIEW.md](ELLIPTIC_REVIEW.md). Four moments do not prove a
complete limiting law. The ternary diagnostic counts are descriptive, and
neither manuscript establishes the first critical intensity or the critical
void probability. No new GRH assumption is introduced for the four moments.

## Provenance and scope

The integration uses the atlas at commit
452907920dc623d4d1e6ae7675428fe02ed840c4 and four selected research archives.
[EDITORIAL_PROVENANCE.json](EDITORIAL_PROVENANCE.json) records the source
archive hashes and the individual files used. Private consolidated archives,
conversation histories, external review reports and host logs are excluded.

The release atlas32-v1.0.0, its scientific seals, the fixed million data and
earlier primary tests are unchanged. PUBLICATION_SHA256.json remains a
historical publication snapshot. EDITORIAL_SHA256.json separately anchors the
current new or edited files, excluding that manifest itself. No original
manifest is silently rewritten to suggest a new scientific campaign.

The coordinate exporter has one packaging adaptation: its classifier source
is found beside the exporter, rather than under the original temporary research
directory. The local source checksum is updated accordingly; original source
hashes remain in the provenance record. Its complete HITS.csv and SUMMARY.csv
reproduce byte for byte against the retained ordinary-position catalogue.

## Reproduction

From the repository root, with standard-library Python 3.10 or newer:

    python3 -B scripts/cubic_quadruplets.py 1000033
    python3 -B reproducibility/proofs/residual_irreducibility.py
    python3 -B scripts/verify_hypergeometric.py
    python3 -B -m unittest discover -s tests
    python3 -B scripts/verify_editorial.py

To also reclassify every retained ordinary long-run position, use a C++17
compiler and run:

    python3 -B scripts/verify_editorial.py --catalogue

This reads the retained million-profile supplement and its 1,234,808 ordinary
positions. It does not rescan the full character words. The verifier rejects
Python optimization because inherited coordinate checks use assertions.
No default command or CI job starts a GPU calculation or downloads the full
multi-gigabyte atlas.

Both LaTeX sources are standalone. Compile each twice with pdflatex and
no-shell-escape. CI retains the earlier checks, adds the curated reproductions
and catalogue comparison, and compiles both manuscripts, rejecting overfull
boxes and unresolved references. Numerical checks do not machine-verify the
analytic proof.

Repository visibility, rights and licence status are unchanged. The existing
rights documentation applies to the new manuscript and curated packages too.
The [reproduction map](PUBLICATION_REPRODUCIBILITY.md) distinguishes current
version identifiers from the public immutable archive still to be arranged.
