#!/usr/bin/env python3
"""Cross-check the Greene--Lennon derivation against literal character sums.

This checks the new exposition at small primes, including A=0, the cubic
phase, infinity and opposite parameters. It is not a proof of the cited
theorems. The archived geometric verifier is retained unchanged.
"""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / "reproducibility/cubic_quadruplets/elliptic"
sys.path.insert(0, str(PACK / "src"))
sys.path.insert(0, str(PACK / "input"))
from core import char_labels
from elliptic import primes_upto
from jacobi_sector_inherited import ROOTS, add, conj, mul


def main():
    cases = zero_parameters = nontrivial_phases = 0
    for p in primes_upto(97):
        if p % 3 != 1:
            continue
        labels = char_labels(p)

        def chi(x):
            x %= p
            return ROOTS[labels[x]] if x else (0, 0)

        def hypergeometric(s):
            value = (0, 0)
            for u in range(p):
                value = add(value, mul(mul(conj(chi(u)), chi(1-u)),
                                       conj(chi(1-s*u))))
            return value  # p times the normalized function, in Z[omega].

        for c in range(2, min(p, 7)):
            for d in range(c+1, p):
                a, b = 0, 1
                denominator = (a-b)*(c-d) % p
                A = (4*(a*b+c*d)-2*(a+b)*(c+d))*pow(denominator, -1, p) % p
                rho = (c-d)*pow(a-b, -1, p) % p
                s = (2-A)*pow(4, -1, p) % p
                if s in (0, 1):
                    raise ArithmeticError("A distinct quadruple gave a singular parameter.")
                original = (1, 0)  # The projective source point at infinity.
                for z in range(p):
                    original = add(original, chi((z+a)*(z+b)*(z+c)**2*(z+d)**2))
                    if (z+b) % p:
                        u = (a-b)*pow(a-d, -1, p)*(z+d)*pow(z+b, -1, p) % p
                        if (z+c)*(z+d) % p:
                            left = (z+a)*(z+b)*pow((z+c)*(z+d), -1, p) % p
                            right = pow(rho, -1, p)*(1-u)*pow(u*(1-s*u), -1, p) % p
                            if left != right:
                                raise ArithmeticError("The fractional linear identity failed.")
                value = hypergeometric(s)
                # A direct Legendre count on the short Weierstrass model.
                trace = 0
                for x in range(p):
                    rhs = (x**3-(6*A+15)*x+A*A+14*A+22) % p
                    trace -= 0 if not rhs else (1 if pow(rhs, (p-1)//2, p) == 1 else -1)
                if (original != mul(conj(chi(rho)), value)
                        or value != (-trace, 0)
                        or value != hypergeometric((1-s) % p)):
                    raise ArithmeticError("Phase, trace or Greene transformation failed.")
                # Complete the square after the translated 3-torsion model.
                for x in range(3):
                    y = x+1
                    X, Y = x+3, y-3*x+A-2
                    short = (Y*Y-X**3+(6*A+15)*X-(A*A+14*A+22)) % p
                    torsion = (y*y-6*x*y+2*(A-2)*y-x**3) % p
                    if short != torsion:
                        raise ArithmeticError("The 3-torsion coordinate change failed.")
                cases += 1
                zero_parameters += A == 0
                nontrivial_phases += chi(rho) != (1, 0)
    if not zero_parameters or not nontrivial_phases:
        raise ArithmeticError("Degenerate-parameter and nontrivial-phase cases were not covered.")
    print(json.dumps({"status": "PASS", "four_position_cases": cases,
                      "zero_parameter_cases": zero_parameters,
                      "nontrivial_phase_cases": nontrivial_phases,
                      "maximum_prime": 97}, sort_keys=True))


if __name__ == "__main__":
    main()
