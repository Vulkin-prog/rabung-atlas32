"""Exact support-three sector for ternary monochromatic-run correlations.

Integers (a,b) encode a+b*omega, omega**2+omega+1=0.
No GPU/network and no full word is built by the fast evaluator.
Input p must be a prime >3, congruent to 1 modulo 3. A separately
certified prime may skip trial primality; all certificates are checked.
"""
from __future__ import annotations
from collections import Counter
from fractions import Fraction
from math import isqrt
from typing import Iterable

E = tuple[int,int]
ROOTS: tuple[E,...]=((1,0),(0,1),(-1,-1))

def need(ok: bool, message: str) -> None:
    if not ok: raise ValueError(message)

def add(x:E,y:E)->E:return (x[0]+y[0],x[1]+y[1])
def neg(x:E)->E:return (-x[0],-x[1])
def sub(x:E,y:E)->E:return add(x,neg(y))
def scale(n:int,x:E)->E:return n*x[0],n*x[1]
def mul(x:E,y:E)->E:
    a,b=x;c,d=y
    return a*c-b*d,a*d+b*c-b*d

def conj(x:E)->E:return x[0]-x[1],-x[1]
def norm(x:E)->int:return x[0]*x[0]-x[0]*x[1]+x[1]*x[1]
def twice_real(x:E)->int:return 2*x[0]-x[1]

def nearest(n:int,d:int)->int:
    need(d>0,'Positive denominator required')
    return (2*n+d)//(2*d)

def divrem(x:E,y:E)->tuple[E,E]:
    n=norm(y);need(n>0,'Division by zero')
    a,b=mul(x,conj(y));q=(nearest(a,n),nearest(b,n))
    r=sub(x,mul(q,y));need(4*norm(r)<=3*n,'Euclidean norm bound failed')
    return q,r

def gcd_e(x:E,y:E)->tuple[E,int]:
    steps=0
    while norm(y):
        _,z=divrem(x,y);x,y=y,z;steps+=1
    return x,steps

def prime_trial(p:int)->bool:
    if p<2:return False
    if p%2==0:return p==2
    return all(p%d for d in range(3,isqrt(p)+1,2))

def primes_to(n:int)->list[int]:
    if n<2:return []
    sieve=bytearray(b'\x01')*(n+1);sieve[:2]=b'\x00\x00'
    for q in range(2,isqrt(n)+1):
        if sieve[q]:sieve[q*q:n+1:q]=b'\x00'*(((n-q*q)//q)+1)
    return [q for q in range(2,n+1) if sieve[q]]

def primitive_cube_root(p:int,limit:int=256)->int:
    """Bounded deterministic search, not a worst-case polylog guarantee.

A supplied nontrivial cubic root avoids this search entirely. Failure is
explicit; it never substitutes a wrong root or an incomplete answer.
"""
    need(p>3 and p%3==1,'p must be prime and 1 modulo 3')
    for a in range(2,min(p,limit+2)):
        z=pow(a,(p-1)//3,p)
        if z!=1:
            need((z*z+z+1)%p==0,'Not a cube root')
            return z
    raise ValueError('No cubic nonresidue found within the declared search budget')

def jacobi_fast(p:int,zeta:int|None=None,check_prime:bool=True)->tuple[E,int,int]:
    need(p>3 and p%3==1,'Inadmissible modulus')
    if check_prime:need(prime_trial(p),'Composite modulus')
    z=primitive_cube_root(p) if zeta is None else zeta%p
    need(z!=1 and (z*z+z+1)%p==0,'Invalid orientation root')
    g,steps=gcd_e((p,0),(-z,1))
    need(norm(g)==p,'GCD is not a prime element of norm p')
    units=ROOTS+tuple(neg(t) for t in ROOTS)
    candidates=[mul(g,u) for u in units if mul(g,u)[0]%3==2 and mul(g,u)[1]%3==0]
    need(len(candidates)==1,'Primary normalization not unique')
    J=candidates[0]
    need((J[0]+J[1]*z)%p==0,'Wrong Jacobi conjugate')
    return J,z,steps

def char_labels(p:int,z:int,values:Iterable[int])->dict[int,int]:
    mapping={1:0,z:1,z*z%p:2}
    need(len(mapping)==3,'Invalid cubic roots')
    out={}
    for n in values:
        need(n%p!=0,'Character evaluated at zero')
        out[n]=mapping[pow(n%p,(p-1)//3,p)]
    return out

def small_labels(p:int,z:int,bound:int)->tuple[list[int],list[int],tuple[int,...]]:
    qs=primes_to(bound);sig=char_labels(p,z,qs)
    c=[-1,0]+[0]*max(0,bound-1)
    for n in range(2,bound+1):
        q=next(q for q in qs if n%q==0)
        c[n]=(c[n//q]+sig[q])%3
    return c,qs,tuple(sig[q] for q in qs)

def convolution(c:list[int],dmax:int)->list[E]:
    out=[(0,0)]*(dmax+1)
    for d in range(2,dmax+1):
        counts=[0,0,0]
        for a in range(1,d):counts[(c[d]+c[a]+c[d-a])%3]+=1
        out[d]=(counts[0]-counts[2],counts[1]-counts[2])
    return out

def weights_from_labels(c:list[int],k:int)->E:
    C=convolution(c,k+1)
    total=(0,0)
    for d in range(2,k+1):total=add(total,scale(2,C[d]))
    return sub(total,C[k+1])

def support_three(p:int,k:int,zeta:int|None=None,check_prime:bool=True)->dict:
    need(k>=2 and k+2<p,'Window must not identify branch points')
    J,z,steps=jacobi_fast(p,zeta,check_prime)
    c,qs,sig=small_labels(p,z,k+1)
    W=weights_from_labels(c,k)
    value=Fraction(twice_real(mul(J,conj(W))),2*3**k)
    return {'p':p,'k':k,'zeta_mod_p':z,'J':list(J),'norm_J':norm(J),
            'W':list(W),'norm_W':norm(W),'primes':qs,'signature':list(sig),
            'I_support_three':str(value),'gcd_steps':steps}

def factor_vector(n:int,qs:list[int])->tuple[int,...]:
    v=[]
    for q in qs:
        e=0
        while n%q==0:n//=q;e+=1
        v.append(e%3)
    need(n==1,'Incomplete small prime basis')
    return tuple(v)

def weight_polynomial(k:int)->tuple[list[int],dict[tuple[int,...],int]]:
    """W_k = sum weight[e] chi(prod q**e); cubes already cancelled."""
    qs=primes_to(k+1);fac={n:factor_vector(n,qs) for n in range(1,k+2)}
    out=Counter()
    for d in range(2,k+2):
        weight=2 if d<=k else -1
        for a in range(1,d):
            v=tuple((x+y+z)%3 for x,y,z in zip(fac[d],fac[a],fac[d-a]))
            out[v]+=weight
    return qs,{e:w for e,w in out.items() if w}

def evaluate_polynomial(poly:dict[tuple[int,...],int],sig:tuple[int,...])->E:
    counts=[0,0,0]
    for e,w in poly.items():counts[sum(x*y for x,y in zip(e,sig))%3]+=w
    return counts[0]-counts[2],counts[1]-counts[2]

def labels_from_signature(bound:int,qs:list[int],sig:tuple[int,...])->list[int]:
    c=[-1,0]+[0]*max(0,bound-1)
    for n in range(2,bound+1):
        i=next(i for i,q in enumerate(qs) if n%q==0)
        c[n]=(c[n//qs[i]]+sig[i])%3
    return c

def direct_jacobi(p:int,z:int)->E:
    c=char_labels(p,z,range(1,p));counts=[0,0,0]
    for x in range(2,p):counts[(c[x]+c[(1-x)%p])%3]+=1
    return counts[0]-counts[2],counts[1]-counts[2]
