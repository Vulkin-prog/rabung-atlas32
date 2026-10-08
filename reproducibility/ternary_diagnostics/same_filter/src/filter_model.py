"""Independent-colour baseline for the fixed coordinate filter.

All count/expectation calculations use integer arithmetic and Fraction.
No arithmetic run scanner, no GPU, no network. Intervals are inclusive.
"""
from fractions import Fraction as F
from math import gcd, isqrt


def require(test, message):
    if not test:
        raise ValueError(message)


def discriminants(bound):
    require(isinstance(bound, int) and bound >= 1, 'positive discriminant bound required')
    return [d for d in range(-bound, bound+1)
            if d and d % 4 in (0, 1) and (d < 0 or isqrt(d)**2 != d)]


def parameters(amax=128, dmax=4096, qmax=128):
    require(isinstance(amax,int) and amax>=0 and isinstance(qmax,int) and qmax>=1,
            'invalid coordinate bounds')
    ds=discriminants(dmax)
    lengths={0:sum(a%2==0 for a in range(-amax,amax+1)),
             1:sum(a%2==1 for a in range(-amax,amax+1))}
    # Each reduced rational in (0,1/2) corresponds to one noninteger cut.
    rational=sum(gcd(a,q)==1 for q in range(3,qmax+1) for a in range(1,(q+1)//2))
    return {'amax':amax,'dmax':dmax,'qmax':qmax,'D':len(ds),
            'even_D':sum(d%4==0 for d in ds),'odd_D':sum(d%4==1 for d in ds),
            'L':sum(lengths[d%2] for d in ds),'R':rational,
            'J':len(ds)+rational,'min_p_strict':(amax+2)**2+dmax}


def mean_upper(k, r, L, J):
    """Uniform upper bound for expected marked, maximal, interior runs >=k."""
    require(k>=1 and r>=2 and L>=0 and J>=0,'invalid mean parameters')
    return F((r-1)*(L+(k-1)*J)+J,r**k)


def innovation_bounds(a,b,r,J):
    """Per-profile E(T_b)-r^(a-b)E(T_a); interval endpoints inclusive."""
    require(1<=a<b and r>=2 and J>=0,'invalid continuation parameters')
    if J==0:
        return F(0),F(0)
    scale=F((r-1)*(b-a),r**b)
    return -scale,J*scale


def merged_intervals(intervals,lo,hi):
    out=[]
    for a,b in sorted(intervals):
        require(a<=b,'reversed mask interval')
        a,b=max(a,lo),min(b,hi)
        if a>b: continue
        if out and a<=out[-1][1]+1:
            out[-1]=(out[-1][0],max(b,out[-1][1]))
        else: out.append((a,b))
    return out


def complement_lengths(h, intervals=(), cut_edges=()):
    """Components of [2,h-1] after marked vertices and crossed edges removed.

    A noninteger cut at c removes the edge (floor(c), floor(c)+1).
    A run avoids the filter exactly when its interior interval lies in one
    of these components. First and last word positions are never eligible
    for strictly interior runs.
    """
    require(isinstance(h,int) and h>=2,'word length must be >=2')
    if h==2:return []
    blocks=merged_intervals(intervals,2,h-1)
    free=[];left=2
    for a,b in blocks:
        if left<a: free.append((left,a-1))
        left=b+1
    if left<=h-1:free.append((left,h-1))
    cuts=sorted(set(t for t in cut_edges if 2<=t<h-1))
    ans=[];ix=0
    for a,b in free:
        while ix<len(cuts) and cuts[ix]<a:ix+=1
        left=a
        while ix<len(cuts) and cuts[ix]<b:
            ans.append(cuts[ix]-left+1)
            left=cuts[ix]+1;ix+=1
        ans.append(b-left+1)
    require(all(l>0 for l in ans),'nonpositive complement length')
    return ans


def coverage(h, ell, components):
    require(ell>=1,'positive run length required')
    return max(0,h-ell-1)-sum(max(0,l-ell+1) for l in components)


def F_exact(L,k,r):
    require(L>=0 and k>=1 and r>=2,'invalid segment parameters')
    if L<k:return F(0)
    return F((r-1)*(L-k+1)-1,r**k)+F(1,r**(L+1))


def F_envelope(L,k,r,cutoff=96):
    require(L>=0 and k>=1 and r>=2 and cutoff>=k,'invalid envelope parameters')
    if L<k:return F(0),F(0)
    main=F((r-1)*(L-k+1)-1,r**k)
    if L<cutoff:
        tail=F(1,r**(L+1));return main+tail,main+tail
    return main,main+F(1,r**(cutoff+1))


def filtered_mean(h,k,r,components,cutoff=96):
    """Rigorous enclosure of E(marked runs); huge powers are never expanded."""
    total_lo,total_hi=F_envelope(h-2,k,r,cutoff)
    for L in components:
        a,b=F_envelope(L,k,r,cutoff)
        total_lo-=b;total_hi-=a
    return max(F(0),total_lo),max(F(0),total_hi)


def is_prime_trial(n):
    if n<2:return False
    if n%2==0:return n==2
    return all(n%d for d in range(3,isqrt(n)+1,2))


def sqrt_mod(n,p):
    """Tonelli--Shanks, assuming p odd prime; None for a nonsquare."""
    n%=p
    if n==0:return 0
    if pow(n,(p-1)//2,p)!=1:return None
    if p%4==3:return pow(n,(p+1)//4,p)
    q=p-1;s=0
    while q%2==0:q//=2;s+=1
    z=2
    while pow(z,(p-1)//2,p)!=p-1:z+=1
    c=pow(z,q,p);x=pow(n,(q+1)//2,p);t=pow(n,q,p);m=s
    while t!=1:
        i=1;tt=t*t%p
        while tt!=1:
            tt=tt*tt%p;i+=1
            require(i<m,'invalid prime/root input')
        b=pow(c,1<<(m-i-1),p)
        x=x*b%p;t=t*b*b%p;c=b*b%p;m=i
    require(x*x%p==n,'square-root verification')
    return x


def coordinate_mask(p,amax=128,dmax=4096,qmax=128):
    """All candidate coordinates, independent of the character word.

    Callers must certify primality separately. Large-p guard guarantees that
    each split D contributes one full interval in each reflected half.
    """
    par=parameters(amax,dmax,qmax)
    require(p%2==1 and p>max(par['min_p_strict'],qmax,5),'prime too small for block geometry')
    h=(p-1)//2;inv2=(p+1)//2;intervals=[];split=0
    for D in discriminants(dmax):
        z=sqrt_mod(D,p)
        if z is None:continue
        split+=1;got=0
        a=amax if amax%2==D%2 else amax-1
        amin=-amax if (-amax)%2==D%2 else -amax+1
        count=(a-amin)//2+1
        if count<=0:continue
        for root in (z,p-z):
            left=(root-a)*inv2%p;right=left+count-1
            require(right<p,'unexpected wrapping anchor interval')
            require(not(left<=h<=right or left<=h+1<=right or left==0),
                    'quadratic interval touches a forbidden boundary')
            if right<h:
                require(left>=2,'quadratic interval reaches first ordinary boundary')
                intervals.append((left,right));got+=1
        require(got==1,'a split discriminant did not give exactly one half-block')
    cuts=[]
    for q in range(3,qmax+1):
        for a in range(1,(q+1)//2):
            if gcd(a,q)==1:
                require(a*p%q!=0,'rational cut is an integer')
                cuts.append(a*p//q)
    require(len(cuts)==par['R'],'rational cut count')
    comps=complement_lengths(h,intervals,cuts)
    return {'p':p,'h':h,'split_discriminants':split,'intervals':intervals,
            'cut_edges':cuts,'components':comps,'parameters':par}
