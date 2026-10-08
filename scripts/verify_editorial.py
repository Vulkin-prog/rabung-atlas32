#!/usr/bin/env python3
"""Replay the curated cubic certificates and ternary same-filter diagnostics.

The external automorphy and prime-number theorems are mathematical inputs.
--catalogue additionally reclassifies the complete retained ordinary catalogue
on CPU; it never enumerates the full character words or uses a GPU.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def run(path, *args):
    subprocess.run([sys.executable, "-B", str(ROOT / path), *map(str, args)],
                   cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", action="store_true")
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("Run without -O: inherited coordinate checks use assertions.")
    run("reproducibility/proofs/residual_irreducibility.py")
    run("scripts/verify_hypergeometric.py")
    run("reproducibility/cubic_quadruplets/elliptic/src/verify.py")
    run("reproducibility/cubic_quadruplets/moments/src/verify.py")
    run("reproducibility/ternary_diagnostics/coordinates/src/analyse.py")
    run("reproducibility/ternary_diagnostics/coordinates/tests/verify_results.py")
    run("reproducibility/ternary_diagnostics/same_filter/src/analyse.py")
    run("reproducibility/ternary_diagnostics/same_filter/tests/verify.py")
    if args.catalogue:
        pack = ROOT / "reproducibility/ternary_diagnostics/coordinates"
        with tempfile.TemporaryDirectory(prefix="rabung-coordinates-") as folder:
            output = Path(folder) / "export"
            run("reproducibility/ternary_diagnostics/coordinates/src/export_classified.py",
                "--repo", ROOT, "--output-dir", output, "--seconds", 450)
            for name in ["HITS.csv", "SUMMARY.csv"]:
                if (output / name).read_bytes() != (pack / "input/coordinate_export" / name).read_bytes():
                    raise ArithmeticError(f"Full catalogue differs: {name}")
            report = json.loads((output / "EXPORT.json").read_text())
            if report["input_tsv_sha256"] != "79862cfc07141191eba77c2c4a75d3534274576a61c8b6f969d12db40fdf4541":
                raise ArithmeticError("Retained ordinary catalogue differs.")
    manifest = json.loads((ROOT / "docs/EDITORIAL_SHA256.json").read_text())
    for path, expected in manifest.items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
            raise ArithmeticError(f"Editorial file changed: {path}")
    print(json.dumps({"status": "PASS", "editorial_files": len(manifest),
                      "full_retained_catalogue_reclassified": args.catalogue,
                      "scope": "finite certificates and archived diagnostics; analytic proofs are not machine-verified"}))


if __name__ == "__main__":
    main()
