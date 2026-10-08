"""Exact elementary elliptic arithmetic. No statistical assertions.

A trace certificate proves a point order larger than the Hasse interval width.
Its verification does not repeat the baby-step/giant-step search.
"""
from __future__ import annotations
from math import isqrt
from functools import lru_cache
from typing import Optional

Point = Optional[tuple[int, int]]

def primes_upto(n: int) -> list[int]:
    if n < 2:
        return []
    a = bytearray(b'\x01') * (n + 1)
    a[0:2] = b'\x00\x00'
    for q in range(2, isqrt(n) + 1):
        if a[q]:
            start = q*q
            a[start:n+1:q] = b'\x00' * ((n-start)//q+1)
    return [i for i, b in enumerate(a) if b]

PRIMES = primes_upto(100_000)

@lru_cache(maxsize=512)
def prime_exact(n: int) -> bool:
    """Trial division, with an explicit domain below 10^10."""
    if not isinstance(n, int) or n < 2 or n >= 10**10:
        return False
    for q in PRIMES:
        if q*q > n:
            return True
        if n % q == 0:
            return n == q
    return True  # All primes up to sqrt(n) were tested in the declared domain.

def factor_exact(n: int) -> list[list[int]]:
    if not isinstance(n, int) or not 1 <= n < 10**10:
        raise ValueError('factor domain')
    out = []
    for q in PRIMES:
        if q*q > n:
            break
        if n % q == 0:
            e = 0
            while n % q == 0:
                n //= q
                e += 1
            out.append([q, e])
    if n > 1:
        out.append([n, 1])
    return out

class Curve:
    def __init__(self, p: int, a: int, b: int):
        if p <= 3 or p >= 10**10 or not prime_exact(p):
            raise ValueError('p is not a prime in the implemented domain')
        self.p, self.a, self.b = p, a % p, b % p
        if (4*self.a**3 + 27*self.b**2) % p == 0:
            raise ValueError('singular curve')
        z = 2
        while pow(z, (p-1)//2, p) != p-1:
            z += 1
        q, s = p-1, 0
        while q % 2 == 0:
            q //= 2
            s += 1
        self.sqrt_data = (q, s, pow(z, q, p))

    def on_curve(self, P: Point) -> bool:
        return P is None or ((P[1]*P[1] - (P[0]**3 + self.a*P[0] + self.b)) % self.p == 0)

    def neg(self, P: Point) -> Point:
        return None if P is None else (P[0], (-P[1]) % self.p)

    def add(self, P: Point, Q: Point) -> Point:
        if P is None:
            return Q
        if Q is None:
            return P
        p = self.p
        x, y = P
        u, v = Q
        if x == u:
            if (y+v) % p == 0:
                return None
            slope = (3*x*x+self.a)*pow(2*y, -1, p) % p
        else:
            slope = (v-y)*pow(u-x, -1, p) % p
        w = (slope*slope-x-u) % p
        return (w, (slope*(x-w)-y) % p)

    def mul(self, n: int, P: Point) -> Point:
        if n < 0:
            return self.mul(-n, self.neg(P))
        R = None
        while n:
            if n & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            n >>= 1
        return R

    def sqrt(self, a: int) -> Optional[int]:
        p = self.p
        a %= p
        if a == 0:
            return 0
        if pow(a, (p-1)//2, p) != 1:
            return None
        q, s, c = self.sqrt_data
        x, t, m = pow(a, (q+1)//2, p), pow(a, q, p), s
        while t != 1:
            z, i = t, 0
            while z != 1 and i < m:
                z = z*z % p
                i += 1
            if i == m:
                raise ArithmeticError('Tonelli failure')
            b = pow(c, 1 << (m-i-1), p)
            x = x*b % p
            c = b*b % p
            t = t*c % p
            m = i
        if x*x % p != a:
            raise ArithmeticError('square root check')
        return min(x, p-x)

    def literal_order(self) -> int:
        p = self.p
        out = p+1
        for x in range(p):
            v = (x*x*x+self.a*x+self.b) % p
            if v:
                out += 1 if pow(v, (p-1)//2, p) == 1 else -1
        return out

    def hasse_candidates(self, P: Point) -> list[int]:
        """All integers n in the Hasse interval such that [n]P=O."""
        p = self.p
        r = isqrt(4*p)
        lo, hi = p+1-r, p+1+r
        width = hi-lo
        step = isqrt(width)+1
        baby: dict[Point, list[int]] = {}
        Q = None
        for j in range(step):
            baby.setdefault(Q, []).append(j)
            Q = self.add(Q, P)
        jump = self.neg(Q)
        target = self.neg(self.mul(lo, P))
        candidates = set()
        for i in range(width//step+1):
            for j in baby.get(target, []):
                n = lo+i*step+j
                if n <= hi:
                    candidates.add(n)
            target = self.add(target, jump)
        return sorted(candidates)

    def certify_order(self, max_x: int = 256) -> dict:
        """Bounded deterministic search. Refuses rather than guesses on failure."""
        p = self.p
        width = 2*isqrt(4*p)
        for x in range(min(max_x, p)):
            y = self.sqrt(x*x*x+self.a*x+self.b)
            if y is None:
                continue
            P = (x, y)
            cand = self.hasse_candidates(P)
            if len(cand) != 1:
                continue
            N = cand[0]
            factors = factor_exact(N)
            order = N
            for q, _ in factors:
                while order % q == 0 and self.mul(order//q, P) is None:
                    order //= q
            if order <= width:
                continue
            cert = {'p': p, 'a': self.a, 'b': self.b, 'point': [x,y],
                    'group_order': N, 'trace': p+1-N,
                    'point_order': order, 'factorization': factor_exact(order)}
            verify_trace(cert, check_p=False)
            return cert
        raise RuntimeError(f'No large-order certificate for p={p},a={self.a},b={self.b}; no trace accepted')

def verify_trace(c: dict, check_p: bool = True) -> int:
    p = c['p']
    if check_p and not prime_exact(p):
        raise ValueError('not a certified prime')
    E = Curve(p, c['a'], c['b'])
    P = tuple(c['point'])
    if not (0 <= P[0] < p and 0 <= P[1] < p and E.on_curve(P)):
        raise ValueError('bad point')
    r = c['point_order']
    fac = c['factorization']
    product, last = 1, 1
    for q, e in fac:
        if q <= last or e < 1 or not prime_exact(q):
            raise ValueError('bad factorization')
        product *= q**e
        last = q
    if product != r or r <= 2*isqrt(4*p):
        raise ValueError('order not large enough')
    if E.mul(r,P) is not None:
        raise ValueError('not an order multiple')
    for q, _ in fac:
        if E.mul(r//q, P) is None:
            raise ValueError('not exact order')
    N = c['group_order']
    if abs(N-p-1) > isqrt(4*p) or N % r != 0:
        raise ValueError('not unique in Hasse interval')
    if c['trace'] != p+1-N:
        raise ValueError('wrong trace')
    return p+1-N
