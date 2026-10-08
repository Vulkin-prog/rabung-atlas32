"""Portable exact verifier. All finite fields and all error decisions use integers.
Large traces are checked by point orders + Hasse, NOT by the trace search.
Small sums and the million-letter support check use different constructions.
"""
from __future__ import annotations
from collections import Counter
from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from math import comb, isqrt
from pathlib import Path
import hashlib, json, sys
from core import (parameters, partitions4, mod_fraction, elliptic_coefficients,
                  char_labels, aggregate, fourth_sector, frac_json)
from elliptic import Curve, primes_upto, prime_exact, factor_exact, verify_trace
from a4 import a4_from_traces, PARAMETERS
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'input'))
from jacobi_sector_inherited import ROOTS, add, mul, scale, conj, jacobi_fast, direct_jacobi, norm
ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()

def need(ok, text):
    if not ok: raise ArithmeticError(text)

def reject(fn):
    try: fn()
    except (ValueError,ArithmeticError,RuntimeError,KeyError,TypeError,IndexError):
        COUNTS['negative_tests']+=1
        return
    raise ArithmeticError('altered object was accepted')

def lam(x,p):
    x%=p
    return 0 if not x else (1 if pow(x,(p-1)//2,p)==1 else -1)

def zchi(c,x,p):
    return (0,0) if x%p==0 else ROOTS[c[x%p]]

# Polynomials over Z in two formal variables x,A.
def pa(x,y):
    z=dict(x)
    for e,v in y.items():z[e]=z.get(e,0)+v
    return {e:v for e,v in z.items() if v}
def ps(x,k):return {e:k*v for e,v in x.items() if k*v}
def pm(x,y):
    z={}
    for (i,j),v in x.items():
        for (a,b),w in y.items():z[i+a,j+b]=z.get((i+a,j+b),0)+v*w
    return {e:v for e,v in z.items() if v}
def pp(x,n):
    z={(0,0):1}
    for _ in range(n):z=pm(z,x)
    return z

def symbolic():
    one={(0,0):1}; x={(1,0):1}; A={(0,1):1}
    h=pa(x,ps(one,-3)); am2=pa(A,ps(one,-2))
    aa=pa(ps(A,-6),ps(one,-15)); bb=pa(pa(pp(A,2),ps(A,14)),ps(one,22))
    t=ps(am2,-12); v=ps(pp(am2,2),4)
    rnum=pa(pa(pm(x,pp(h,2)),pm(t,h)),v)
    dnum=pa(pa(pp(h,3),ps(pm(t,h),-1)),ps(v,-2))
    ca=pa(ps(A,54),ps(one,-135)); cb=pa(pa(ps(pp(A,2),-27),ps(A,378)),ps(one,-594))
    rhs=pa(pa(pp(rnum,3),pm(pm(ca,rnum),pp(h,4))),pm(cb,pp(h,6)))
    lhs=pm(pp(dnum,2),pa(pa(pp(x,3),pm(aa,x)),bb))
    need(lhs==rhs,'formal Velu identity');COUNTS['formal_polynomial_identities']+=1
    disc=ps(pa(ps(pp(aa,3),4),ps(pp(bb,2),27)),-16)
    expected=ps(pm(pp(am2,3),pa(A,ps(one,2))),-432)
    need(disc==expected,'formal discriminant');COUNTS['formal_polynomial_identities']+=1
    # Birational transformation of the quartic quotient; x now denotes w=1/(u+2).
    btrans=pa(pa(pa(one,ps(x,-6)),ps(pp(x,2),9)),pm(am2,pp(x,3)))
    XX=pa(pm(am2,x),ps(one,3))
    rhs=pa(pa(pp(XX,3),pm(aa,XX)),bb)
    need(pm(pp(am2,2),btrans)==rhs,'quartic to Weierstrass');COUNTS['formal_polynomial_identities']+=1
    return {'variables':['x','A'],'velu_identity_nonzero_terms_after_subtraction':0,
            'discriminant':'-432*(A-2)^3*(A+2)'}

def small_curves():
    examples=[]
    for p in primes_upto(199):
        if p%3!=1:continue
        c=char_labels(p)
        for A in range(p):
            if A in (2,p-2):continue
            E=Curve(p,*elliptic_coefficients(A,p))
            t=p+1-E.literal_order()
            opp=p+1-Curve(p,*elliptic_coefficients(-A,p)).literal_order()
            need(t==opp,'opposite parameter traces'); COUNTS['opposite_trace_equalities']+=1
            mixed=(0,0)
            for w in range(p):mixed=add(mixed,scale(lam(w*w+A*w+1,p),zchi(c,w,p)))
            need(mixed==(-t,0),'mixed sum equals -trace');COUNTS['mixed_sum_trace_equalities']+=1
            sex=sum(lam(x**6+A*x**3+1,p) for x in range(p))
            sp=sum(lam((u+2)*(u**3-3*u+A),p) for u in range(p))
            sm=sum(lam((u-2)*(u**3-3*u+A),p) for u in range(p))
            need(-1-sex==(-1-sp)+(-1-sm)==2*t,'two quotient count identity')
            COUNTS['genus_two_quotient_counts']+=1
            Q=(3%p,(A-2)%p)
            need(E.on_curve(Q) and E.mul(3,Q) is None and Q is not None,'order three point')
            COUNTS['order_three_points']+=1
            target=Curve(p,(54*A-135)%p,(-27*A*A+378*A-594)%p)
            tt=-12*(A-2)%p; vv=4*(A-2)**2%p
            for x in range(min(p,12)):
                if x==3:continue
                y=E.sqrt(x**3+E.a*x+E.b)
                if y is None:continue
                inv=pow(x-3,-1,p)
                xx=(x+tt*inv+vv*inv*inv)%p
                yy=(1-tt*inv*inv-2*vv*inv**3)*y%p
                need(target.on_curve((xx,yy)),'Velu image over finite field')
                COUNTS['finite_velu_images']+=1
        # all three partitions of small offset quadruples, actual complex char sums.
        for q4 in combinations(range(min(7,p)),4):
            for ab,cd in partitions4(q4):
                A,rho=parameters(ab,cd);a,b=ab;d,e=cd
                t=p+1-Curve(p,*elliptic_coefficients(A,p)).literal_order()
                got=(1,0)
                for z in range(p):
                    value=((z+a)*(z+b)*(z+d)**2*(z+e)**2)%p
                    got=add(got,zchi(c,value,p))
                phase=zchi(c,mod_fraction(rho,p),p)
                target=scale(-t,conj(phase))
                need(got==target,'four point complex identity');COUNTS['four_point_complex_identities']+=1
                if phase!=(1,0) and t and len(examples)<3:
                    examples.append({'p':p,'offsets':list(q4),'partition':[list(ab),list(cd)],
                       'parameter':frac_json(A),'rho':frac_json(rho),
                       'cubic_phase':list(phase),'trace':t,'S_plus_one':list(got),
                       'untwisted_wrong_value':[-t,0]})
    return examples

def literal_fourth(p,k,c):
    total=0
    for q4 in combinations(range(k+1),4):
        coeff=-1 if q4[0]==0 else 2
        for ab,cd in partitions4(q4):
            a,b=ab;d,e=cd
            for z in range(1,p-k):
                phase=(c[z+a]+c[z+b]+2*c[z+d]+2*c[z+e])%3
                total+=coeff*(2 if phase==0 else -1)
    return Fraction(total,3**k)

def small_sectors():
    for p in primes_upto(109):
        if p%3!=1:continue
        c=char_labels(p)
        for k in range(4,min(8,p-2)+1):
            weights,_=aggregate(k,p)
            traces={A:p+1-Curve(p,*elliptic_coefficients(A,p)).literal_order() for A,w in weights.items() if w}
            value,meta=fourth_sector(k,p,traces)
            need(value==literal_fourth(p,k,c),'aggregate vs literal monomials')
            need(value==histogram_fourth(c,p,k)[0],'aggregate vs histogram')
            need(meta['pair_coefficient_sum']==3*(2*comb(k,4)-comb(k,3)),'constant coefficient sum')
            COUNTS['literal_fourth_sector_comparisons']+=1

def K(n,k,s):
    def cb(a,b):return comb(a,b) if 0<=b<=a else 0
    return sum(cb(n,r)*2**r*cb(k-n,s-r)*(-1)**(s-r) for r in range(s+1))

def histogram_fourth(c,p,k):
    K3=[K(n,k,3) for n in range(k+1)];K4=[K(n,k,4) for n in range(k+1)]
    totals=[0,0,0]
    for a in c[2:k+2]:totals[a]+=1
    numerator=0; direct_starts=0
    for z in range(1,p-k):
        value=sum(2*K4[totals[a]]-(2 if c[z]==a else -1)*K3[totals[a]] for a in range(3))
        need(value%3==0,'histogram character projection')
        numerator+=value//3
        direct_starts+=int(max(totals)==k and totals[c[z]]==0)
        if z+1<p-k:
            totals[c[z+1]]-=1;totals[c[z+k+1]]+=1
    return Fraction(numerator,3**k),direct_starts

def small_a4():
    for p in primes_upto(3000):
        if p%3!=1:continue
        traces={A:p+1-Curve(p,*elliptic_coefficients(A,p)).literal_order() for A in PARAMETERS}
        val,meta=a4_from_traces(p,traces)
        c=char_labels(p)
        direct=sum(c[x]==c[x+1]==c[x+2]==c[x+3] for x in range(1,p-3))
        need(val==direct,'cubic A4 literal comparison')
        J,z,_=jacobi_fast(p,check_prime=False)
        need(J==direct_jacobi(p,z),'Jacobi orientation direct')
        need(norm(J)==p and J[0]%3==2 and J[1]%3==0,'Jacobi primary normalization')
        COUNTS['literal_cubic_A4_comparisons']+=1
        COUNTS['literal_Jacobi_comparisons']+=1

def verify_big():
    stored=json.loads((ROOT/'results/CERTIFICATS_ELLIPTIQUES.json').read_text())
    summary=[]
    for row in stored['cases']:
        p,k=row['p'],row['k'];need(prime_exact(p),'large prime')
        weights,_=aggregate(k,p)
        traces={}; seen=set()
        for cert in row['trace_certificates']:
            A=Fraction(cert['parameter']['numerator'],cert['parameter']['denominator'])
            need(A not in seen and weights[A]!=0 and cert['weight']==weights[A],'weight or duplicate')
            need((cert['a'],cert['b'])==elliptic_coefficients(A,p) and cert['p']==p,'wrong elliptic curve')
            seen.add(A);traces[A]=verify_trace(cert,check_p=False)
            COUNTS['large_trace_certificates']+=1
        need(seen=={A for A,w in weights.items() if w},'missing trace')
        _,meta=fourth_sector(k,p,traces)
        for key,v in meta.items():need(row[key]==v,'stored fourth result '+key)
        summary.append(meta)
    a4=json.loads((ROOT/'results/A4_CUBIQUE.json').read_text())
    for row in a4['cases']:
        p=row['p'];need(prime_exact(p),'a4 prime')
        traces={}
        for cert in row['elliptic_certificates']:
            A=Fraction(cert['parameter']['numerator'],cert['parameter']['denominator'])
            need((cert['a'],cert['b'])==elliptic_coefficients(A,p) and cert['p']==p,'a4 wrong curve')
            traces[A]=verify_trace(cert,check_p=False);COUNTS['A4_trace_certificate_rechecks']+=1
        _,meta=a4_from_traces(p,traces)
        for key,v in meta.items():need(row[key]==v,'stored A4 '+key)
    return summary

def million():
    p=1000033;k=13
    # Distinct construction: primitive root walk of F_p^*, not Euler exponentiation.
    fac=[q for q,e in factor_exact(p-1)]
    g=2
    while any(pow(g,(p-1)//q,p)==1 for q in fac):g+=1
    J,z,_=jacobi_fast(p,check_prime=False)
    labg={1:0,z:1,z*z%p:2}[pow(g,(p-1)//3,p)]
    c=[None]*p;x=1
    for exponent in range(p-1):
        need(c[x] is None,'primitive root cycle repeated')
        c[x]=(labg*exponent)%3;x=x*g%p
    need(x==1 and all(a is not None for a in c[1:]),'incomplete primitive root cycle')
    # Every entry checked independently by integer modular exponentiation.
    lookup={1:0,z:1,z*z%p:2}
    for x in range(1,p):need(c[x]==lookup[pow(x,(p-1)//3,p)],'Euler/walk discrepancy')
    COUNTS['million_Euler_vs_walk_values']+=p-1
    val,Z=histogram_fourth(c,p,k)
    expected=json.loads((ROOT/'results/CERTIFICATS_ELLIPTIQUES.json').read_text())['cases'][0]['fourth_sector']
    need(frac_json(val)==expected,'million support-four comparison')
    aa=sum(c[x]==c[x+1]==c[x+2]==c[x+3] for x in range(1,p-3))
    ar=json.loads((ROOT/'results/A4_CUBIQUE.json').read_text())['cases'][0]['A4_cubic']
    need(aa==ar,'million A4 comparison')
    COUNTS['million_direct_window_starts']+=p-k-1
    digest=hashlib.sha256(bytes(c[1:])).hexdigest()
    return {'p':p,'k':k,'primitive_root':g,'word_labels_sha256':digest,'fourth_sector':frac_json(val),
            'A4_cubic_direct':aa,'whole_word_start_count':Z,
            'purpose':'finite identity validation, not an estimate of a vertical law'}

def negative_tests():
    row=json.loads((ROOT/'results/CERTIFICATS_ELLIPTIQUES.json').read_text())['cases'][0]
    cert=row['trace_certificates'][0]
    c=deepcopy(cert);c['trace']+=1;reject(lambda:verify_trace(c))
    c=deepcopy(cert);c['group_order']+=1;reject(lambda:verify_trace(c))
    c=deepcopy(cert);c['point'][1]=(c['point'][1]+1)%c['p'];reject(lambda:verify_trace(c))
    c=deepcopy(cert);c['factorization'][0][1]+=1;reject(lambda:verify_trace(c))
    c=deepcopy(cert);c['point_order']*=2;reject(lambda:verify_trace(c))
    c=deepcopy(cert);c['p']=1000005;reject(lambda:verify_trace(c))
    reject(lambda:parameters((0,1),(1,2)))
    reject(lambda:mod_fraction(Fraction(1,7),7))
    reject(lambda:elliptic_coefficients(2,13))
    reject(lambda:elliptic_coefficients(-2,13))
    reject(lambda:aggregate(4,11))
    reject(lambda:Curve(13,0,0))
    reject(lambda:Curve(1000033,*elliptic_coefficients(14,1000033)).certify_order(max_x=0))


def main():
    formal=symbolic()
    twists=small_curves()
    small_sectors();small_a4();big=verify_big();local=million();negative_tests()
    obj={'status':'PASS','scope':'exact finite identities and supplied trace certificates; no universal maximum law',
         'counts':dict(sorted(COUNTS.items())),'formal':formal,'twist_counterexamples':twists,
         'large_signed_sectors':big,'million_independent_check':local,
         'prime_words_near_2pow32_scanned':[],
         'not_checked':['global prime-averaged first intensity','probability of no critical runs',
                        'all supports >=5','completeness of the bounded point search outside supplied certificates',
                        'worldwide novelty']}
    (ROOT/'results/VERIFICATION.json').write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','counts':obj['counts'],'million':local},sort_keys=True,indent=2))

if __name__=='__main__':main()
