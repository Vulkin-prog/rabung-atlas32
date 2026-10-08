"""Exact algebra for the first four moments; no prime equidistribution is assumed in code.
The justification for replacing prime averages by these functionals is in RAPPORT_FR.md.
Python standard library only.
"""
from fractions import Fraction
from itertools import product
from math import comb, factorial

E_CURVES = ((-99, 414), (-21, 37), (-2835, 58158))
C_CURVES = ((0, -27), (0, -243), (0, -2187))


def compositions(n: int, parts: int):
    if parts == 1:
        yield (n,)
        return
    for v in range(n + 1):
        for tail in compositions(n - v, parts - 1):
            yield (v,) + tail


def one_moment(n: int, cm: bool) -> int:
    if n < 0:
        raise ValueError('nonnegative order required')
    if n % 2:
        return 0
    return comb(n, n // 2) if cm else comb(n, n // 2) // (n // 2 + 1)


def multinomial_moment(coefficients, degree: int) -> int:
    if len(coefficients) != 4 or not 0 <= degree <= 4:
        raise ValueError('four variables; certified moment range 0..4')
    out = 0
    for v in compositions(degree, 4):
        term = factorial(degree)
        for r in v:
            term //= factorial(r)
        for j, r in enumerate(v):
            term *= coefficients[j] ** r * one_moment(r, j == 0)
        out += term
    return out


def laurent_moment(coefficients, degree: int) -> int:
    """Independent computation: constant term against Weyl weights, exact integers."""
    if len(coefficients) != 4 or not 0 <= degree <= 4:
        raise ValueError('four variables; certified moment range 0..4')
    pol = {(0, 0, 0, 0): 1}
    for _ in range(degree):
        new = {}
        for exp, value in pol.items():
            for j, c in enumerate(coefficients):
                for direction in (-1, 1):
                    e = list(exp)
                    e[j] += direction
                    e = tuple(e)
                    new[e] = new.get(e, 0) + c * value
        pol = new
    ans = 0
    for bits in product((0, 1), repeat=3):
        exp = (0,) + tuple(-2 * b for b in bits)
        ans += (-1) ** sum(bits) * pol.get(exp, 0)
    return ans


def route(h: int, m: tuple[int, int, int], u: int) -> dict:
    """Bookkeeping of the analytic proof, NOT an automorphy or PNT verifier."""
    if len(m) != 3 or any(v < 0 for v in m) or abs(h) + sum(m) > 4:
        raise ValueError('outside the proved total-degree-four domain')
    if u not in (0, 1, 2):
        raise ValueError('finite character exponent must be reduced modulo three')
    indices = [j for j, v in enumerate(m) if v]
    if not indices:
        return {'type': 'Hecke', 'constant': h == 0 and u == 0}
    if len(indices) == 1:
        return {'type': 'Rankin-Selberg', 'dimensions': [m[indices[0]] + 1, 1],
                'pole_check': 'different dimensions'}
    if len(indices) == 2:
        a, b = indices
        return {'type': 'Rankin-Selberg', 'dimensions': [m[a] + 1, m[b] + 1],
                'pole_check': 'witness p=31 if equal dimensions'}
    if sum(m) == 3:
        return {'type': 'tensor then Rankin-Selberg', 'dimensions': [4, 2],
                'pole_check': 'different dimensions'}
    if sorted(m) == [1, 1, 2]:
        return {'type': 'tensor then Rankin-Selberg', 'dimensions': [3, 4],
                'pole_check': 'different dimensions'}
    raise ArithmeticError('unhandled admissible monomial')


def moment_table() -> dict:
    weights = {'cube': (-4, -2, -2, -2), 'noncube': (2, -2, -2, 1)}
    table = {}
    for name, coeff in weights.items():
        vals = [multinomial_moment(coeff, n) for n in range(5)]
        if vals != [laurent_moment(coeff, n) for n in range(5)]:
            raise ArithmeticError('moment implementations disagree')
        table[name] = {'coefficients_CM_then_E': list(coeff), 'moments_0_to_4': vals}
    table['all'] = {'moments_0_to_4': [
        int((Fraction(table['cube']['moments_0_to_4'][r], 3)
          + Fraction(2 * table['noncube']['moments_0_to_4'][r], 3)))
        for r in range(5)]}
    return table
