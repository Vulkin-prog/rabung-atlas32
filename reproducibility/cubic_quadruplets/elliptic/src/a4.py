"""Exact cubic four-letter count from three elliptic traces and one Jacobi sum.
The inherited Jacobi normalization is kept separate with its provenance.
"""
from pathlib import Path
from fractions import Fraction
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'input'))
from jacobi_sector_inherited import jacobi_fast, mul, add, conj, twice_real, ROOTS
from core import char_labels
from elliptic import prime_exact

PARAMETERS=(Fraction(14),Fraction(1),Fraction(10,3))

def a4_from_traces(p:int,traces:dict) -> tuple[int,dict]:
    if p%3!=1 or not prime_exact(p) or p<=3:
        raise ValueError('prime p=1 mod3 required, p<10^10')
    J,z,steps=jacobi_fast(p,check_prime=False)
    labels=char_labels(p,6)
    # Same orientation: both routines use the first nontrivial Euler value.
    lookup={1:0,z:1,z*z%p:2}
    if any(labels[a] != lookup[pow(a,(p-1)//3,p)] for a in range(1,7)):
        raise ArithmeticError('Jacobi/character orientation mismatch')
    chi2,chi3,chi6=(ROOTS[labels[a]] for a in (2,3,6))
    jacobi_term=2*twice_real(mul(J,add(conj(chi2),conj(chi6))))
    I2=int(labels[2]==0)
    I23=int(labels[2]==labels[3]==0)
    beta3=2 if labels[3]==0 else -1
    elliptic_term=-2*traces[PARAMETERS[0]]-2*traces[PARAMETERS[1]]-beta3*traces[PARAMETERS[2]]
    num=p-26-18*I2-18*I23+jacobi_term+elliptic_term
    if num%27 or not 0 <= num//27 <= p-4:
        raise ArithmeticError('nonintegral or out-of-range cubic A4')
    return num//27, {'p':p,'jacobi':list(J),'zeta_mod_p':z,'euclidean_steps':steps,
        'chi2_label':labels[2],'chi3_label':labels[3], 'I2':I2,'I23':I23,
        'beta3':beta3,'jacobi_term':jacobi_term,'elliptic_term':elliptic_term,
        'numerator':num,'A4_cubic':num//27}
