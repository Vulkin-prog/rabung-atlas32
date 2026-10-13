# Dependencies, references and withheld material

No third-party runtime, downloaded source distribution or executable is included
in this repository. The frozen engine uses external CUDA and CUB headers from
the user's toolchain. CPU re-analysis links external zlib and OpenSSL libraries.
PARI/GP is an optional external tool for certificate regeneration. LaTeX and
Poppler are external document tools. Python scripts use its standard library.
Each dependency retains the licence and terms of its own distribution.

The optional PARI documentation reference is pinned to the official
[pari-2.19.0.tar.gz release](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz),
announced on 28 September 2026, including its shipped manual. The official
distribution SHA256 is
`f317b9722eb5d9094a60303774f066f3a83e3ec1f170be8546c44d7583f30b6d`.
The [official download page](https://pari.math.u-bordeaux.fr/download.html)
provides that checksum. The distribution and manual are not redistributed
here; this pinned reference replaces the moving `html-stable` citation.

Historical paths called `vendor/research/`, `vendor/` and `baseline/` are project
snapshots inherited by the census engine. File identities were checked against
the frozen campaign source record. No external licence notice was found in the
selected source files; this observation is not proof of ownership. Their
redistribution is subject to the owner's rights confirmation described in
`RIGHTS.md`. No new licence is applied to inherited material by this deposit.

Externally supplied French reports are **not redistributed**. Their rights were
not established. Editorial reports and internal working notes are not cited as
scientific sources in any manuscript. The mathematical arguments are given
in the articles and their reproducibility packages. The review texts are not
included in source, Git history or release assets.
References to papers and official documentation are bibliographic citations and
links, not copies of those documents.

Machine-specific logs, sensor readings, process identifiers, account credentials,
absolute home paths, software installations and executable binaries are omitted.
Compressed research records are audited after decompression as well as on disk.
Publication-derived indices contain mathematical records and commitments only;
they are rebuilt rather than copied with deleted SQLite pages.

The certificates addendum imports nine author-owned files from the already
published v1.0.0 artifact. Their original Git identities and SHA256 values are
listed in `reproducibility/bounds_update/PUBLISHED_IMPORTS.json`; the original
artifact's component licences remain applicable to those copies. This does
not grant a new blanket licence to the assembled repository.

The GPU selected package imports author-owned scientific source extracts,
configuration and measured data from the pinned methods packet. Source code
remains byte-identical. Its thirteen explicit JSON publication copies remove
host-path and disk-only metadata; timings, profiles and mathematical values
are unchanged. The full original packet and its build environment are not
redistributed here. Source paths containing `vendor/research` have the same
rights-confirmation boundary as the existing census snapshots.
