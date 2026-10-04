"""Exact quotient-field arithmetic and elementary irreducibility certificates.

The implementation uses only Python integers and fractions. Elements have
one common denominator; this avoids repeated fractional gcds in products.
"""
from fractions import Fraction
import math


def prime(q):
    return q>=2 and all(q%d for d in range(2,math.isqrt(q)+1))


def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def modpoly(a,P,q):
    a=trim([x%q for x in a]);d=len(P)-1
    while a!=[0] and len(a)>d:
        k=len(a)-d-1;b=a[-1]*pow(P[-1],-1,q)%q
        for j in range(d+1):a[k+j]=(a[k+j]-b*P[j])%q
        a=trim(a)
    return a


def pmul(a,b,P,q):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return modpoly(z,P,q)


def ppow(a,e,P,q):
    v=[1]
    while e:
        if e&1:v=pmul(v,a,P,q)
        e//=2
        if e:a=pmul(a,a,P,q)
    return trim(v)


def pgcd(a,b,q):
    a=trim([x%q for x in a]);b=trim([x%q for x in b])
    while b!=[0]:a,b=b,modpoly(a,b,q)
    return trim([x*pow(a[-1],-1,q)%q for x in a])


def irreducible_mod(P,q):
    d=len(P)-1
    if P[-1]%q==0:return False
    divisors=[k for k in range(2,d+1) if d%k==0 and prime(k)]
    required={d//k for k in divisors}
    x=modpoly([0,1],P,q);f=x
    for j in range(1,d+1):
        f=ppow(f,q,P,q)
        if j in required:
            diff=[0]*max(len(f),len(x))
            for i,a in enumerate(f):diff[i]+=a
            for i,a in enumerate(x):diff[i]-=a
            if len(pgcd(P,diff,q))!=1:return False
    return trim(f)==trim(x)


def irreducibility_prime(P):
    for q in range(2,5001):
        if prime(q) and irreducible_mod(P,q):return q
    raise ValueError('No irreducible specialization below 5001')


class Field:
    def __init__(self,P):
        self.P=[int(x) for x in P];assert self.P[-1]==1
        self.d=len(P)-1
        self.zero=((0,)*self.d,1)
        self.one=self.integer(1)

    def normal(self,a,b):
        a=list(a)+[0]*max(0,self.d-len(a));assert len(a)==self.d and b
        if b<0:a=[-x for x in a];b=-b
        g=math.gcd(b,*a)
        return tuple(x//g for x in a),b//g

    def element(self,a):
        v=[Fraction(x) for x in a];den=math.lcm(*(x.denominator for x in v))
        return self.normal([x.numerator*(den//x.denominator) for x in v],den)

    def integer(self,n):
        return self.normal([n],1)

    def add_integer(self,a,n):
        v,b=a;v=list(v);v[0]+=n*b
        return self.normal(v,b)

    def multiply(self,a,b):
        aa,ad=a;bb,bd=b;z=[0]*(2*self.d-1)
        for i,x in enumerate(aa):
            if x:
                for j,y in enumerate(bb):
                    if y:z[i+j]+=x*y
        for k in range(2*self.d-2,self.d-1,-1):
            v=z[k]
            if v:
                for j in range(self.d):z[k-self.d+j]-=v*self.P[j]
        return self.normal(z[:self.d],ad*bd)

    def power(self,a,n):
        assert n>=0;v=self.one
        while n:
            if n&1:v=self.multiply(v,a)
            n//=2
            if n:a=self.multiply(a,a)
        return v

    def evaluate(self,coeff_ascending,a):
        v=self.zero
        for c in reversed(coeff_ascending):v=self.add_integer(self.multiply(v,a),int(c))
        return v


def verify_decompositions(field,units,matrix,targets):
    assert len(matrix)==len(units) and all(len(row)==len(targets) for row in matrix)
    assert all(u!=field.zero for u in units)
    cache={}
    for j,target in enumerate(targets):
        left,right=field.one,target
        for i,row in enumerate(matrix):
            n=row[j]
            if not n:continue
            key=(i,abs(n))
            if key not in cache:cache[key]=field.power(units[i],abs(n))
            if n>0:left=field.multiply(left,cache[key])
            else:right=field.multiply(right,cache[key])
        assert left==right,('exact multiplicative identity',j)
    return len(targets)
