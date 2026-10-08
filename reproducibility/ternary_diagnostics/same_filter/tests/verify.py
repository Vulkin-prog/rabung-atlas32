import sys,json,itertools,random,math,hashlib
from pathlib import Path
from fractions import Fraction as F
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from filter_model import *


def explicit_mark(start,end,points,cuts):
    return any(start<=t<=end for t in points) or any(start<=t and t+1<=end for t in cuts)


def runs(word):
    ans=[];a=0
    for i in range(1,len(word)+1):
        if i==len(word) or word[i]!=word[a]:
            if a>0 and i<len(word):ans.append((a+1,i))
            a=i
    return ans


def exact_marked_expectation(h,k,r,comps):
    return F_exact(h-2,k,r)-sum((F_exact(L,k,r) for L in comps),F(0))


def main():
    count={};rng=random.Random(6138415)
    # Closed finite geometric sum versus literal summation.
    n=0
    for r in range(2,10):
        for L in range(0,25):
            for k in range(1,12):
                direct=sum((F((L-ell+1)*(r-1)**2,r**(ell+1)) for ell in range(k,L+1)),F(0))
                require(F_exact(L,k,r)==direct,'geometric sum mismatch');n+=1
    count['finite_geometric_sum_equalities']=n
    words=masks=checks=0
    # Direct colour words; two-sided maximality is tested by literal scanning.
    for h in range(3,10):
        for r in [2,3]:
            configs=[([],[]),([(2,h-1)],[]),([],list(range(2,h-1)))]
            for t in range(4):
                pts=[j for j in range(2,h) if rng.randrange(3)==0]
                cuts=[j for j in range(2,h-1) if rng.randrange(3)==0]
                configs.append(([(j,j) for j in pts],cuts))
            sums=[[0]*h for _ in configs]
            for word in itertools.product(range(r),repeat=h):
                words+=1;rr=runs(word)
                for i,(ints,cuts) in enumerate(configs):
                    pts=[j for a,b in ints for j in range(a,b+1)]
                    for a,b in rr:
                        if explicit_mark(a,b,pts,cuts):
                            for k in range(1,b-a+2):sums[i][k]+=1
            for i,(ints,cuts) in enumerate(configs):
                masks+=1;comps=complement_lengths(h,ints,cuts)
                pts=[j for a,b in ints for j in range(a,b+1)]
                for k in range(1,h):
                    z=exact_marked_expectation(h,k,r,comps)
                    require(z==F(sums[i][k],r**h),'full word baseline mismatch');checks+=1
                    for ell in [k]:
                        W=sum(explicit_mark(a,a+ell-1,pts,cuts) for a in range(2,h-ell+1))
                        require(W==coverage(h,ell,comps),'coverage identity mismatch')
    count.update(words_enumerated=words,masks_tested_against_words=masks,expectation_equalities=checks)
    # Exhaustive site masks and edge masks, with rational-cut effects distinct.
    nm=slopes=innov=0
    for h in range(3,8):
        for bits in range(1<<((h-2)+(h-3))):
            pts=[j for j in range(2,h) if bits>>(j-2)&1]
            cuts=[j for j in range(2,h-1) if bits>>(h-2+j-2)&1]
            ints=merged_intervals([(j,j) for j in pts],2,h-1)
            J=len(ints)+len(cuts);comps=complement_lengths(h,ints,cuts);nm+=1
            Ws=[coverage(h,ell,comps) for ell in range(1,h+2)]
            for x,y in zip(Ws,Ws[1:]):
                require(-1<=y-x<=J,'coverage slope bound');slopes+=1
            for a in range(1,h-1):
                for b in range(a+1,h):
                    z=exact_marked_expectation(h,b,3,comps)-F(1,3**(b-a))*exact_marked_expectation(h,a,3,comps)
                    low,high=innovation_bounds(a,b,3,J)
                    require(low<=z<=high,'signed continuation envelope');innov+=1
            for k in range(1,h):
                z=exact_marked_expectation(h,k,3,comps)
                require(0<=z<=mean_upper(k,3,len(pts),J),'mask first moment bound')
    count.update(exhaustive_masks=nm,coverage_slope_checks=slopes,signed_continuation_bounds=innov)
    # Direct modular enumeration versus root-interval generation on small cases.
    profiles=positions=baselines=0
    for p in [101,127,211,503,997,2003]:
        require(is_prime_trial(p),'small prime verification')
        for A,D,Q in [(2,12,7),(4,16,9),(6,64,16)]:
            if p<=parameters(A,D,Q)['min_p_strict']:continue
            mask=coordinate_mask(p,A,D,Q);h=(p-1)//2
            expanded=set(j for a,b in mask['intervals'] for j in range(a,b+1))
            ds=set(discriminants(D));direct=set()
            for y in range(1,h+1):
                positions+=1
                for a in range(-A,A+1):
                    z=(2*y+a)**2%p
                    opts=[z,z-p]
                    if any(d in ds and (a*a-d)%4==0 for d in opts):
                        direct.add(y);break
            require(expanded==direct,'modular mask disagreement')
            for ell in range(1, min(h-1,12)):
                w=sum(explicit_mark(x,x+ell-1,direct,mask['cut_edges']) for x in range(2,h-ell+1))
                require(w==coverage(h,ell,mask['components']),'modular coverage mismatch')
            for k in [1,3,7]:
                lo,hi=filtered_mean(h,k,3,mask['components'],cutoff=12)
                exact=exact_marked_expectation(h,k,3,mask['components'])
                require(lo<=exact<=hi,'finite tail enclosure');baselines+=1
            profiles+=1
    count.update(modular_mask_cases=profiles,modular_positions=positions,finite_tail_enclosures=baselines)
    # Independently count D and A pairs, and totients for cut budgets.
    for A,D,Q in [(64,256,128),(128,4096,128)]:
        p=parameters(A,D,Q)
        alt=sum(1 for a in range(-A,A+1) for d in range(-D,D+1)
                if d and (a*a-d)%4==0 and (d<0 or math.isqrt(d)**2!=d))
        require(alt==p['L'],'polynomial count budget')
        phi=sum(sum(math.gcd(a,q)==1 for a in range(1,q)) for q in range(3,Q+1))//2
        require(phi==p['R'],'totient budget')
    # Invalid input paths and parameters are rejected.
    bad=[lambda:parameters(-1,4,2),lambda:discriminants(0),lambda:mean_upper(0,3,1,1),
         lambda:innovation_bounds(4,3,3,1),lambda:coordinate_mask(17,128,4096),
         lambda:F_envelope(100,30,3,20),lambda:merged_intervals([(4,3)],1,8),
         lambda:complement_lengths(1,[],[])]
    for f in bad:
        try:f()
        except ValueError:pass
        else:raise AssertionError('invalid input accepted')
    count['negative_tests']=len(bad)
    data=json.loads((ROOT/'results/ANALYSE.json').read_text())
    require(data['parameters']['L']==518112 and data['parameters']['J']==6542,'reported parameters')
    for row in data['envelopes']:
        u=F(row['iid_expected_sum_upper']['fraction'])
        low=F(row['observed_minus_iid_lower']['fraction'])
        require(u+low==row['observed_marked_runs'],'reported finite baseline')
    res={'status':'PASS','checks':count,'scope':'Finite combinatorial models and mask geometry; no re-reading of the million, no new primality or full-run audit of its entries.'}
    (ROOT/'results/VERIFICATION.json').write_text(json.dumps(res,indent=2,sort_keys=True)+'\n')
    print(json.dumps(res,indent=2))

if __name__=='__main__':main()
