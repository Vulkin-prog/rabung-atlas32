#!/usr/bin/env python3
"""Exact cubic four-window count from six certified elliptic traces.

The count is cubic and is separate from the historical quadratic A4 column.
No Jacobi evaluator, statistical estimate, GPU, or long word is used.
"""
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "reproducibility/cubic_quadruplets"
sys.path.insert(0, str(PACK / "elliptic/src"))
from elliptic import Curve, prime_exact, verify_trace

E_CURVES = ((-99, 414), (-21, 37), (-2835, 58158))
C_CURVES = ((0, -27), (0, -243), (0, -2187))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def admissible(p):
    require(type(p) is int and 3 < p < 10**10 and p % 3 == 1,
            "A prime p > 3, p = 1 mod 3, p < 10^10 is required.")
    require(prime_exact(p), "Composite modulus.")


def certified_traces(p, certificates, curves):
    require(len(certificates) == len(curves), "Three certificates are required.")
    traces = []
    for cert, (a, b) in zip(certificates, curves):
        require(cert["p"] == p, "Certificate modulus differs from the query.")
        require(cert["a"] == a % p and cert["b"] == b % p,
                "Certificate belongs to a different curve.")
        traces.append(verify_trace(cert, check_p=False))
    return traces


def evaluate(p, noncm_certificates, cm_certificates):
    admissible(p)
    t = certified_traces(p, noncm_certificates, E_CURVES)
    c = certified_traces(p, cm_certificates, C_CURVES)
    cube2 = pow(2, (p - 1) // 3, p) == 1
    cube3 = pow(3, (p - 1) // 3, p) == 1
    require(c[0] == c[1] == c[2] if cube3 else sum(c) == 0,
            "The certified CM traces violate the cubic twist identity.")
    correction = 26 + 18 * int(cube2) + 18 * int(cube2 and cube3)
    beta3 = 2 if cube3 else -1
    fluctuation_numerator = -2*c[0] - 2*c[1] - 2*t[0] - 2*t[1] - beta3*t[2]
    numerator = p - correction + fluctuation_numerator
    require(numerator % 27 == 0 and 0 <= numerator // 27 <= p - 4,
            "Nonintegral or out-of-range four-window count.")
    return {"p": p, "A4_cubic": numerator // 27,
            "CM_traces": c, "nonCM_traces": t,
            "cube2": cube2, "cube3": cube3,
            "finite_correction": correction,
            "corrected_fluctuation_numerator": fluctuation_numerator,
            "method": "six point-order and Hasse certificates"}


def stored_certificates(p):
    a = json.loads((PACK / "elliptic/results/A4_CUBIQUE.json").read_text())
    b = json.loads((PACK / "moments/results/CERTIFICATS_CM.json").read_text())
    old = next((r for r in a["cases"] if r["p"] == p), None)
    new = next((r for r in b["cases"] if r["p"] == p), None)
    require(old is not None and new is not None,
            "No stored six-curve certificate set for this prime. Use --generate explicitly.")
    ec = [dict(c) for c in old["elliptic_certificates"]]
    # The inherited third certificate uses the rational E_(10/3) model.
    # Move its point and coefficients through (X,Y) = (9x,27y).
    last = ec[2]
    require(last.get("parameter") == {"numerator": 10, "denominator": 3},
            "Unexpected inherited third-curve parameter.")
    last["a"] = last["a"] * 81 % p
    last["b"] = last["b"] * 729 % p
    last["point"] = [last["point"][0] * 9 % p, last["point"][1] * 27 % p]
    return ec, new["CM_certificates"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p", type=int)
    parser.add_argument("--generate", action="store_true",
                        help="Search bounded point-order certificates on six fixed curves.")
    args = parser.parse_args()
    try:
        admissible(args.p)
        if args.generate:
            ec = [Curve(args.p, a, b).certify_order() for a, b in E_CURVES]
            cc = [Curve(args.p, a, b).certify_order() for a, b in C_CURVES]
        else:
            ec, cc = stored_certificates(args.p)
        result = evaluate(args.p, ec, cc)
        if args.generate:
            result["nonCM_certificates"] = ec
            result["CM_certificates"] = cc
        print(json.dumps(result, indent=2, sort_keys=True))
    except (ValueError, ArithmeticError, RuntimeError, KeyError) as exc:
        parser.exit(2, f"cubic_quadruplets: {exc}\n")


if __name__ == "__main__":
    main()
