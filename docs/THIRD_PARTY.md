# Dependencies, references and withheld material

No third-party runtime, downloaded source distribution or executable is included
in this repository. The frozen engine uses external CUDA and CUB headers from
the user's toolchain. CPU re-analysis links external zlib and OpenSSL libraries.
PARI/GP is an optional external tool for certificate regeneration. LaTeX and
Poppler are external document tools. Python scripts use its standard library.
Each dependency retains the licence and terms of its own distribution.

Historical paths called `vendor/research/`, `vendor/` and `baseline/` are project
snapshots inherited by the census engine. File identities were checked against
the frozen campaign source record. No external licence notice was found in the
selected source files; this observation is not proof of ownership. Their
redistribution is subject to the owner's rights confirmation described in
`RIGHTS.md`. No new licence is applied to inherited material by this deposit.

Externally supplied French reports are **not redistributed**. Their rights were
not established. The manuscript may cite a privately supplied review as a
source; the review text is not included in source, Git history or release assets.
References to papers and official documentation are bibliographic citations and
links, not copies of those documents.

Machine-specific logs, sensor readings, process identifiers, account credentials,
absolute home paths, software installations and executable binaries are omitted.
Compressed research records are audited after decompression as well as on disk.
Publication-derived indices contain mathematical records and commitments only;
they are rebuilt rather than copied with deleted SQLite pages.
