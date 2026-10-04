"""Frozen scientific functions; no GPU, launch or host monitoring side effects."""
import struct
from collections import Counter
from bisect import bisect_left, bisect_right

ROW = struct.Struct('<34H')


def vu(out,n):
    assert 0<=n<2**64
    while n>=128:out.append((n&127)|128);n>>=7
    out.append(n)


def process(value,prior):
    p=value['p'];half=(p-1)//2;e=(p-1)//3;z=value['character_root'];roots=[1,z,z*z%p]
    M,L,sig,flags,m18,m21=prior
    assert value['M_bound']==M and value['L']==L and M>=18
    phase2=roots.index(pow(2,e,p));assert bool(flags&4)==(phase2==0)
    runs=value['runs'];ordinary=[];middle=0;critical=[];last=0
    for x,length,c in runs:
        assert x>last and length>=18 and 0<=c<=2;last=x
        colour=roots[c];assert pow(x,e,p)==colour and pow(x+length-1,e,p)==colour
        if x>1:assert pow(x-1,e,p)!=colour
        if x+length<p:assert pow(x+length,e,p)!=colour
        if length>=21:
            for t in range(x+1,x+length-1):assert pow(t,e,p)==colour
            critical.append((x,length,c))
        if x+length-1>half:
            assert 2*x+length==p+1;middle=length
        elif x!=1:ordinary.append((x,length,c))
    Ns=[sum(l>=k for x,l,c in ordinary)for k in range(18,28)]
    eligible=bool(flags&2)
    assert (M<21)==(eligible and Ns[3]==0)
    colourcounts=Counter(c for x,l,c in ordinary);bins=[0]*16
    for x,l,c in ordinary:bins[min(15,x*16//half)]+=1
    starts=[x for x,l,c in ordinary];echo=[]
    for m in [2,3]:
        phase=roots.index(pow(m,e,p));count=0;k=18
        for j,(b,l,c)in enumerate(ordinary):
            lower=m*b;upper=m*(b+k-1);branch=2*lower//p
            if branch!=2*upper//p:continue
            def fold(v):v%=p;return min(v,p-v)
            a,h=sorted((fold(lower),fold(upper)))
            for i in range(bisect_left(starts,a),bisect_right(starts,h-k+1)):
                if i==j:continue
                assert (ordinary[i][2]-c)%3==phase and (ordinary[i][0]-a)%m!=1
                count+=1
        echo.append(count)
    fields=[M,L,middle,*Ns,*[colourcounts[c]for c in range(3)],*echo,*bins]
    assert len(fields)==34 and all(0<=x<=65535 for x in fields)
    positions=bytearray();vu(positions,len(critical));previous=0
    for x,l,c in critical:vu(positions,x-previous);vu(positions,l-21);positions.append(c);previous=x
    return ROW.pack(*fields),positions,len(runs),len(critical)
