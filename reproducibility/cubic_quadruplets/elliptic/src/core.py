"""Exact signed fourth-support reduction; no claim on the remaining supports."""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb
from elliptic import prime_exact


def partitions4(t):
    a,b,c,d=t
    return [((a,b),(c,d)),((a,c),(b,d)),((a,d),(b,c))]

def parameters(pair1, pair2):
    a,b=pair1; c,d=pair2
    if len({a,b,c,d}) != 4:
        raise ValueError('four distinct offsets required')
    den=(a-b)*(c-d)
    num=4*(a*b+c*d)-2*(a+b)*(c+d)
    return Fraction(num,den), Fraction(c-d,a-b)

def mod_fraction(x, p):
    if x.denominator % p == 0:
        raise ValueError('denominator vanishes')
    return x.numerator*pow(x.denominator,-1,p) % p

def elliptic_coefficients(A, p):
    A=mod_fraction(A,p) if isinstance(A,Fraction) else A % p
    if p <= 3 or A in (2,p-2):
        raise ValueError('bad reduction')
    return (-6*A-15)%p, (A*A+14*A+22)%p

def char_labels(p, length=None):
    """Canonical orientation fixed by the first nontrivial Euler value."""
    e=(p-1)//3
    if p % 3 != 1 or not prime_exact(p):
        raise ValueError('not cubic admissible')
    g=2
    while pow(g,e,p)==1:
        g+=1
    rho=pow(g,e,p)
    lookup={1:0,rho:1,rho*rho%p:2}
    n=p-1 if length is None else length
    return [None]+[lookup[pow(a,e,p)] for a in range(1,n+1)]

def aggregate(k, p):
    if not 4 <= k < p-1 or p % 3 != 1:
        raise ValueError('invalid k or p')
    labels=char_labels(p,k)
    groups=defaultdict(int)
    holes=0
    pairs=0
    raw_abs=0
    csum=0
    for offsets in combinations(range(k+1),4):
        c=-1 if offsets[0]==0 else 2
        outside=[r for r in range(k+1) if r not in offsets]
        for pair1,pair2 in partitions4(offsets):
            A,rho=parameters(pair1,pair2)
            phase=(labels[abs(pair2[1]-pair2[0])]-labels[abs(pair1[1]-pair1[0])])%3
            beta=2 if phase==0 else -1
            groups[abs(A)]+=c*beta
            pairs+=1; csum+=c; raw_abs+=abs(c*beta)
            a,b=pair1; d,e=pair2
            for r in outside:
                exp=(labels[abs(a-r)]+labels[abs(b-r)]+2*labels[abs(d-r)]+2*labels[abs(e-r)])%3
                holes+=c*(2 if exp==0 else -1)
    return dict(groups), {'p':p,'k':k,'conjugate_pairs':pairs,
        'raw_monomials':2*pairs,'rational_parameters':len(groups),
        'nonzero_parameters':sum(v!=0 for v in groups.values()),
        'raw_absolute_weight':raw_abs,'grouped_absolute_weight':sum(abs(v) for v in groups.values()),
        'coefficient_sum':sum(groups.values()),'pair_coefficient_sum':csum,
        'constant_term_numerator':2*csum,'holes_numerator':holes}

def fourth_sector(k, p, traces):
    groups, meta=aggregate(k,p)
    packet=sum(w*traces[A] for A,w in groups.items() if w)
    value=Fraction(-meta['constant_term_numerator']-packet-meta['holes_numerator'],3**k)
    meta.update({'elliptic_packet':packet,'fourth_sector':frac_json(value)})
    return value,meta

def frac_json(q):
    q=Fraction(q)
    return {'numerator':q.numerator,'denominator':q.denominator}
