#!/usr/bin/env python3
"""Small CPU-only checks of v0.2's source bridges, not a CUDA qualification.

Compiles only the extracted scalar anchor policy with a host C++ compiler.
No nvcc, worker, profiler, campaign launcher, or GPU API is called.
"""
from pathlib import Path
import hashlib
import json
import random
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
import sys
import argparse

if sys.flags.optimize:
    raise SystemExit("Run without -O: inherited finite checks use assertions.")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PAPER = ROOT / "paper"
DATA = HERE


def check(report_path=None):
    source_map = json.loads((HERE / "SOURCE_MAP.json").read_text())
    for name, meta in source_map["qualification_records"].items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest() == meta["sha256"]
    oracle_record = source_map["oracle_build_record"]
    assert hashlib.sha256((HERE/oracle_record["file"]).read_bytes()).hexdigest() == oracle_record["sha256"]
    for name, meta in source_map["unchanged_sources"].items():
        buf = (DATA / name).read_bytes()
        assert hashlib.sha256(buf).hexdigest() == meta["sha256"], name
        lines = buf.decode().splitlines()
        for symbol, indices in meta["symbols"].items():
            assert indices, (name, symbol)
            assert all(symbol in lines[i-1] for i in indices)
    compiler = shutil.which("g++") or shutil.which("c++")
    if compiler is None:
        raise RuntimeError("A host C++ compiler is needed; CUDA is not needed")
    policies = []
    for bits in (32, 33):
        text = (DATA / f"source/engine{bits}/src/gpu_pipeline.cuh").read_text()
        start = text.index("u32 anchor_t(")
        end = text.index("\n__global__", start)
        function = text[start:end]
        prefix = "#include <algorithm>\n#include <cstdint>\n#include <iostream>\n"
        prefix += "using u32=std::uint32_t;using std::min;using std::max;\n"
        prefix += "int __clz(u32 x){return __builtin_clz(x);}\n"
        driver = "\nint main(){for(u32 q=3;q<=256;++q)std::cout<<q<<' '<<anchor_t(q,13)<<'\\n';}\n"
        with tempfile.TemporaryDirectory(prefix="vdw-paper-policy-") as directory:
            src, exe = Path(directory) / "policy.cpp", Path(directory) / "policy"
            src.write_text(prefix + function + driver)
            subprocess.run([compiler, "-std=c++17", "-O0", str(src), "-o", str(exe)],
                           check=True, timeout=15, capture_output=True)
            output = subprocess.run([str(exe)], check=True, timeout=2,
                                    capture_output=True, text=True).stdout
        pairs = [tuple(map(int, row.split())) for row in output.splitlines()]
        assert len(pairs) == 254
        assert dict(pairs)[12] == 6 and dict(pairs)[13] == 7
        fallbacks = []
        for q, t in pairs:
            assert 3 <= t <= q
            half = (q+1)//2
            space = q-half+1
            power = 1 << (space.bit_length()-1)
            proposed = q-power+1
            if proposed < 3:
                assert t == 3
            elif proposed > half+2:
                assert t == max(3, half)
                fallbacks.append(q)
            else:
                assert q-t+1 == power
        # Both fixed t (32-bit) and monotone t (wider) stay within a rising q.
        epochs = 0
        for old_q, old_proposal in pairs:
            for new_q, new_proposal in pairs:
                if new_q < old_q:
                    continue
                assert max(old_q, new_proposal) <= new_q
                assert 3 <= old_proposal <= new_q
                assert 3 <= max(old_proposal, new_proposal) <= new_q
                epochs += 1
        policies.append({"bits": bits, "tested_q": [3, 256], "values": len(pairs),
                         "general_period_fallbacks": fallbacks,
                         "monotone_q_pairs": epochs})

    # Check the reciprocal algorithm's scalar IEEE-754 interpretation against
    # exact integer residues. This is a finite arithmetic check, not a run scan.
    rng = random.Random(20261003)
    moduli = [1009, 257000011, 4294967291, 4294992241, 8589931561]
    cases = 0
    errors = set()
    wrapped = 0
    for p in moduli:
        edges = [0, 1, 2, p//2, p-2, p-1]
        inputs = [(a,b) for a in edges for b in edges]
        inputs += [(rng.randrange(p), rng.randrange(p)) for _ in range(512)]
        for a, b in inputs:
            quotient = int((float(a)*float(b))*(1.0/float(p)))
            errors.add(quotient-(a*b//p))
            low = (a*b-quotient*p) & ((1<<64)-1)
            remainder = low if low < (1<<63) else low-(1<<64)
            assert -p <= remainder < 2*p
            if remainder < 0:
                remainder += p
            if remainder >= p:
                remainder -= p
            assert remainder == a*b % p
            wrapped += a*b >= 1<<64
            cases += 1
    assert errors <= {-1,0,1}
    assert wrapped > 0

    # Component counters are sums over concurrent work, not wall percentages.
    diagnostics = []
    for bits, filename, variant, lo in [
        (32, "holdout-32-mr_selected.json", "selected", 2971000001),
        (33, "holdout-33-selected_selected_binary48.json", "selected_binary48", 8041000021),
    ]:
        obj = json.loads((DATA / "throughput" / filename).read_text())
        window = min(obj["windows"], key=lambda w:abs(w["lo"]-lo))
        rows = [r for r in window["samples"] if r["variant"] == variant]
        assert len(rows) == 48
        fields = ["gpu_prepare_ms", "queue_ms", "prepare_ms", "setup_ms", "verify_ms"]
        means = {field:sum(r[field] for r in rows)/6 for field in fields}
        tex = (PAPER / "GPU_METHOD_EN.tex").read_text()
        assert all(f"{v:.2f}" in tex for v in means.values())
        diagnostics.append({"bits":bits, "lo":window["lo"],
                            "mean_per_window_summed_work_ms":means,
                            "overlap_prevents_wall_time_partition":True})

    report = {"status":"PASS", "version":"0.3",
              "operation":"CPU_ONLY_REVISION_BRIDGE_CHECK",
              "utc":datetime.now(timezone.utc).isoformat(),
              "source_map_sha256":hashlib.sha256((HERE/"SOURCE_MAP.json").read_bytes()).hexdigest(),
              "host_compiler":compiler, "anchor_policies":policies,
              "oracle_scalar_interpretation":{"cases":cases, "moduli":moduli,
                 "full_products_exceeding_64_bits":wrapped,
                 "observed_quotient_errors":sorted(errors)},
              "diagnostics":diagnostics, "new_CUDA_executions":0,
              "limitations":["Finite host checks do not verify CUDA code generation, a full run maximum, or all inputs.",
                 "The oracle check models round-to-nearest binary64 and the recorded two's-complement interpretation.",
                 "These finite checks add no CUDA qualification or maximum-run oracle."]}
    if report_path:
        report_path.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({"status":"PASS", "anchor_values_per_engine":254,
                      "oracle_scalar_cases":cases, "new_CUDA_executions":0}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    check(parser.parse_args().report)

