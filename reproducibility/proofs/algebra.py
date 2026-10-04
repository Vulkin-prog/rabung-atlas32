"""Small exact linear algebra helpers (standard library only)."""
from fractions import Fraction

def echelon(vectors):
    piv={}
    for v in vectors:
        while v:
            k=v.bit_length()-1
            if k in piv:v^=piv[k]
            else:piv[k]=v;break
    return piv

def reduce(v,piv):
    for k in sorted(piv,reverse=True):
        if v>>k&1:v^=piv[k]
    return v

def nullspace(rows,n):
    piv=echelon(rows);answers=[]
    for k in range(n):
        if k in piv:continue
        v=1<<k
        for j in sorted(piv):
            if (piv[j]&v).bit_count()%2:v^=1<<j
        assert all((v&r).bit_count()%2==0 for r in rows)
        answers.append(v)
    return answers

def coordinates(v,basis):
    piv={}
    for i,b in enumerate(basis):
        t=1<<i
        while b:
            k=b.bit_length()-1
            if k in piv:w,u=piv[k];b^=w;t^=u
            else:piv[k]=(b,t);break
    answer=0
    while v:
        k=v.bit_length()-1
        if k not in piv:raise ValueError('outside span')
        w,t=piv[k];v^=w;answer^=t
    return answer

def indexes(v,n):return [i for i in range(n) if v>>i&1]
def bits(I):return sum(1<<i for i in I)
def add(F,a,b):
    aa,ad=a;bb,bd=b
    return F.normal([x*bd+y*ad for x,y in zip(aa,bb)],ad*bd)
def scale(F,a,c):
    c=Fraction(c);return F.normal([x*c.numerator for x in a[0]],a[1]*c.denominator)

