# Selected GPU methods reproduction package

The source article is `paper/GPU_METHOD_EN.tex` (v0.3, 9 October 2026).
The complete original experimental packet remains pinned at
`Vulkin-prog/vdw-gpu-starter@f2d515c383776d3f915b87423ea101029d041bc5`;
the v0.2 editorial source is pinned at `f00cdcb95f12d7578bdb37a87db97ff233c86a22`.
These private-repository pins identify versions; they are not public availability
claims. The full packet has 272 imports. This selected package has 34 imports,
not a second full packet or a CUDA build tree.

`PROVENANCE.json` records original paths and hashes as well as current bytes.
Measured source extracts and the configuration remain byte-identical. Thirteen
JSON publication copies replace historical host paths and omit disk-only
metadata. Their timings, profiles and mathematical values are unchanged.
No binaries, host logs, supplied review texts or software installations are
included. `SOURCE_MAP.json` contains the selected source/function and
qualification bridges; its original review-history fields are omitted.

From the repository root:

```bash
python3 -B scripts/verify_gpu_method.py
```

This regenerates both benchmark tables from per-prime medians and six repeated
batch sums, checks all 6,144 final prime identities and 20,512 admissible
profiles, and verifies the selected import hashes. A C++17 compiler checks the
extracted anchor policies for every `q=3..256`. Finite scalar checks compare
the recorded oracle arithmetic with exact Python integers. A literal character
word checks the small `p=101` worked example.

The wrappers adapt the original v0.2 offline checkers to the selected package
and current manuscript paths. They do not alter measured algorithm sources.
Neither default command compiles CUDA, uses a GPU, times an engine, rechecks
full billion-position character words or reproduces the complete frozen
qualification. `check_source_bridges.py --report PATH` writes an optional
current CPU-check report; ordinary checks do not mutate retained inputs.

For complete CUDA build commands, use `paper-gpu/rebuild_workers.py` in the
full pinned packet. The 32-bit and wider GPU domains remain separate and
guarded. The new full CPU checks of larger-prime bounds have their own source
and domain in `../bounds_update/`; they do not extend the methods article's
GPU qualification.
