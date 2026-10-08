#!/usr/bin/env python3
"""Exact degree-preserving irreducibility certificates over Q(i).

Inputs have descending integer coefficients. An irreducible reduction at a
prime q=1 mod 4 proves irreducibility over Q(i). Generation is bounded;
verification uses explicit rejections and works also under python -O.
"""
from __future__ import annotations
import argparse, copy, itertools, json
from math import isqrt
from pathlib import Path

def require(c, message):
    if not c: raise ValueError(message)

def prime(q):
    return isinstance(q,int) and q>=2 and all(q%d for d in range(2,isqrt(q)+1))

def trim(a):
    a=list(a)
    while a and a[-1]==0: a.pop()
    return a

def remainder(a,b,q):
    a=trim([x%q for x in a]); b=trim([x%q for x in b])
    require(bool(b),'zero polynomial divisor'); inv=pow(b[-1],-1,q)
    while len(a)>=len(b):
        c=a[-1]*inv%q; shift=len(a)-len(b)
        for j,x in enumerate(b): a[j+shift]=(a[j+shift]-c*x)%q
        a=trim(a)
    return a

def multiply(a,b,f,q):
    if not a or not b: return []
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]=(c[i+j]+x*y)%q
    return remainder(c,f,q)

def power(a,n,f,q):
    r=[1]; a=remainder(a,f,q)
    while n:
        if n&1: r=multiply(r,a,f,q)
        a=multiply(a,a,f,q); n//=2
    return r

def gcd(a,b,q):
    while b: a,b=b,remainder(a,b,q)
    if not a: return []
    inv=pow(a[-1],-1,q)
    return [x*inv%q for x in a]

def factors(n):
    return [l for l in range(2,n+1) if n%l==0 and prime(l)]

def frobenius_minus_x(n,f,q):
    a=power([0,1],q**n,f,q)
    a += [0]*max(0,2-len(a)); a[1]=(a[1]-1)%q
    return trim(a)

def irreducible(poly,q):
    if not prime(q) or len(poly)<2 or poly[0]%q==0: return False
    f=list(reversed([int(x)%q for x in poly])); d=len(f)-1
    inv=pow(f[-1],-1,q); f=[x*inv%q for x in f]
    if remainder(frobenius_minus_x(d,f,q),f,q): return False
    return all(gcd(f,frobenius_minus_x(d//l,f,q),q)==[1] for l in factors(d))

def check_binding(model,data,witness,degree):
    poly=model['polynomial']
    require(poly==data['polynomial'],'character/list polynomial mismatch')
    require(poly==witness['polynomial'],'irreducibility/list polynomial mismatch')
    require(witness['certificate']==model['certificate'],'certificate identity mismatch')
    require(witness['degree']==degree==len(poly)-1,'degree mismatch')
    q=witness['q']
    require(prime(q) and q%4==1,'witness is not a split Gaussian prime')
    require(irreducible(poly,q),'reduction is not degree-preserving and irreducible')

def generate(root):
    records={}; summary={}
    for d in range(2,7):
        models=json.loads((root/f'RESIDUAL_QUALIFICATION_D{d}.json').read_text())['models']; qs=[]
        for model in models:
            name=model['certificate']; data=json.loads((root/name).read_text())
            q=next((q for q in range(5,10001,4) if irreducible(model['polynomial'],q)),None)
            require(q is not None,'bounded witness search failed')
            record={'certificate':name,'polynomial':model['polynomial'],'degree':d,'q':q}
            check_binding(model,data,record,d); require(name not in records,'duplicate certificate name')
            records[name]=record; qs.append(q)
        summary[str(d)]={'count':len(qs),'max_q':max(qs)}
    require(len(records)==100,'expected exactly 100 candidates')
    return {'schema':1,'base_field':'Q(i)','criterion':'degree-preserving irreducible reduction at q=1 mod 4','source_order':'integer coefficients descending','models':records,'summary':summary}

def self_test():
    total=0
    for q in [2,3,5]:
        for d in range(2,5):
            for lower in itertools.product(range(q),repeat=d):
                poly=[1]+list(lower); f=list(reversed(poly)); trial=True
                for r in range(1,d//2+1):
                    for coeff in itertools.product(range(q),repeat=r):
                        if not remainder(f,list(coeff)+[1],q): trial=False; break
                    if not trial: break
                require(irreducible(poly,q)==trial,'criterion/trial factorization disagreement'); total+=1
    model={'certificate':'example.json','polynomial':[1,0,2]}; data={'polynomial':[1,0,2]}
    w={'certificate':'example.json','polynomial':[1,0,2],'degree':2,'q':5}
    check_binding(model,data,w,2); negatives=[]
    wrong=copy.deepcopy(data); wrong['polynomial']=[1,1,2]; negatives.append((model,wrong,w,2))
    for key,value in [('polynomial',[1,1,2]),('certificate','other'),('degree',3),('q',3),('q',9)]:
        wrong=copy.deepcopy(w); wrong[key]=value; negatives.append((model,data,wrong,2))
    for poly in [[1,0,-1],[5,0,2]]:
        m={'certificate':'example.json','polynomial':poly}; dt={'polynomial':poly}; w1={**w,'polynomial':poly}
        negatives.append((m,dt,w1,2))
    for args in negatives:
        try: check_binding(*args)
        except ValueError: continue
        raise ValueError('altered certificate accepted')
    return {'status':'PASS','independent_trial_factorizations':total,'negative_binding_checks':len(negatives)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--generate',type=Path); p.add_argument('--output',type=Path); args=p.parse_args()
    report=self_test()
    if args.generate:
        data=generate(args.generate); require(args.output is not None,'--output is required')
        args.output.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n'); report['summary']=data['summary']
    print(json.dumps(report,sort_keys=True))
if __name__=='__main__': main()
