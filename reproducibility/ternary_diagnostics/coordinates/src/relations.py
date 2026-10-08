"""Elementary sufficient character certificates in Q[X]/(X^2+uX+v).
The code does not compute complete Kummer groups. All implications use
only identities involving two factors on the left and at most two on the right.
"""
from functools import lru_cache
from math import gcd

@lru_cache(None)
def factor_mod3(a: int) -> dict[int,int]:
    a=abs(a)
    if a==0: raise ValueError('zero coefficient')
    out={};q=2
    while q*q<=a:
        e=0
        while a%q==0:a//=q;e+=1
        if e%3:out[q]=e%3
        q+=1
    if a>1:out[a]=1
    return out

def row_of(rel):
    out={0:(len(rel['left'])-len(rel['right']))%3}
    for val,sign in [(rel['num'],-1),(rel['den'],1)]:
        for p,e in factor_mod3(val).items():out[p]=(out.get(p,0)+sign*e)%3
    return {p:e for p,e in out.items() if e}

def combine(a,b,factor):
    out=a.copy()
    for k,e in b.items():
        v=(out.get(k,0)+factor*e)%3
        if v:out[k]=v
        elif k in out:del out[k]
    return out

def generate(u,v,L):
    rels=[];directions={}
    for i in range(L):
        for j in range(i,L):
            a=i+j-u;b=i*j-v
            if a and b%a==0 and 0<=b//a<L:
                rels.append({'left':[i,j],'right':[b//a],'num':a,'den':1})
            elif a==0 and b:
                rels.append({'left':[i,j],'right':[],'num':b,'den':1})
            g=gcd(a,b)
            if not g:raise ValueError('reducible input')
            if a<0 or (a==0 and b<0):g=-g
            key=(a//g,b//g)
            if key in directions:
                t,z= directions[key]
                rels.append({'left':[i,j],'right':t,'num':g,'den':z})
            else:directions[key]=([i,j],g)
    rels.sort(key=lambda rel:(len(row_of(rel)),abs(rel['num'])+abs(rel['den']),len(rel['left'])+len(rel['right']),rel['left'],rel['right']))
    return rels

@lru_cache(None)
def certificate(u:int,v:int,L:int):
    rels=generate(u,v,L);basis={}
    for i,rel in enumerate(rels):
        row=row_of(rel);comb={i:1}
        while row:
            k=min(row)
            if k not in basis:
                inv=pow(row[k],-1,3)
                basis[k]=({p:inv*e%3 for p,e in row.items()},{p:inv*e%3 for p,e in comb.items()});break
            scale=-row[k];rr,cc=basis[k];row=combine(row,rr,scale);comb=combine(comb,cc,scale)
    def solve(k):
        row={k:1};result={}
        while row:
            p=min(row)
            if p not in basis:return None
            scale=row[p];rr,cc=basis[p];row=combine(row,rr,-scale);result=combine(result,cc,scale)
        return result
    phase=solve(0);two=solve(2)
    used=sorted(set((phase or {}).keys())|set((two or {}).keys()));remap={old:new for new,old in enumerate(used)}
    conv=lambda proof:None if proof is None else [[remap[i],e] for i,e in sorted(proof.items())]
    return {'u':u,'v':v,'length':L,'D':u*u-4*v,'generated_relations':len(rels),'relation_rank_mod3':len(basis),'relations':[rels[i] for i in used],'phase_proof':conv(phase),'two_proof':conv(two)}

def product(indices,u,v):
    a,b=1,0
    for j in indices:a,b=a*j-b*v,a+b*j-b*u
    return (a,b)

def verify_certificate(cert):
    u,v,L=cert['u'],cert['v'],cert['length'];rels=cert['relations']
    assert cert['D']==u*u-4*v
    for rel in rels:
        assert rel['num'] and rel['den']
        assert all(0<=j<L for j in rel['left']+rel['right'])
        left=product(rel['left'],u,v);right=product(rel['right'],u,v)
        assert tuple(rel['den']*a for a in left)==tuple(rel['num']*b for b in right)
    for name,key in [('phase_proof',0),('two_proof',2)]:
        if cert[name] is None:continue
        row={}
        for i,e in cert[name]:
            assert 0<=i<len(rels) and e in (1,2)
            row=combine(row,row_of(rels[i]),e)
        assert row=={key:1}
    return True
