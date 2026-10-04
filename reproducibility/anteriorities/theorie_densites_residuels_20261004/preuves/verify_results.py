"""Portable exact verifier: standard Python only, no PARI, GPU, or network."""
from pathlib import Path
from fractions import Fraction
import json,sys,time,math,itertools
from exact_fields import Field,prime,irreducible_mod
from finite_modules import echelon,kernel4,intersection

def read(root,name):return json.loads((root/name).read_text())
def add(F,a,b):return F.normal([x*b[1]+y*a[1] for x,y in zip(a[0],b[0])],a[1]*b[1])
def scale(F,a,n):return F.normal([n*x for x in a[0]],a[1])
def evalint(P,x,q=None):
    z=0
    for a in P:z=z*x+a;z=z if q is None else z%q
    return z
def bits(v):return sum((x%2)<<j for j,x in enumerate(v))
def conv(a,b):
    z=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return z
def discriminant(P):
    n=len(P)-1;D=[P[i]*(n-i) for i in range(n)];rows=[]
    for k in range(n-1):rows.append([0]*k+P+[0]*(n-2-k))
    for k in range(n):rows.append([0]*k+D+[0]*(n-1-k))
    A=[[Fraction(a) for a in row] for row in rows];det=Fraction(1)
    for k in range(len(A)):
        i=next(i for i in range(k,len(A)) if A[i][k])
        if i!=k:A[k],A[i]=A[i],A[k];det=-det
        pivot=A[k][k];det*=pivot
        for i in range(k+1,len(A)):
            c=A[i][k]/pivot
            for j in range(k+1,len(A)):A[i][j]-=c*A[k][j]
    value=det*(-1)**(n*(n-1)//2)/P[0];assert value.denominator==1
    return int(value)

def verify_galois(cert):
    P=cert['polynomial'];n=len(P)-1;types=set()
    assert discriminant(P)==cert['discriminant']
    for q,t,fac in cert['factorizations']:
        assert prime(q) and cert['discriminant']%q and P[0]%q
        z=[1]
        for f in fac:
            assert irreducible_mod(f,q)
            z=[x%q for x in conv(z,f)]
        expected=[a*pow(P[0],-1,q)%q for a in reversed(P)]
        assert z==expected and sorted(len(f)-1 for f in fac)==t
        types.add(tuple(t))
    assert (n,) in types
    if n==3:assert (1,2) in types
    if n==4:assert (1,3) in types
    return ('C2' if n==2 else 'S3' if n==3 else 'S4')

def product(F,alpha,positions):
    z=F.one
    for j in positions:z=F.multiply(z,F.add_integer(alpha,j))
    return z

def verify_base(cert):
    P=cert['polynomial'];m=cert['M'];F=Field(list(reversed(P)));alpha=F.element([0,1])
    masks=[]
    for rel in cert['square_relations']:
        target=product(F,alpha,rel['positions'])
        assert F.power(F.element(rel['square_root_power_basis']),2)==target
        masks.append(sum(1<<j for j in rel['positions']))
    rows=[];count=0
    for q,x,row in cert['lower_character_samples']:
        assert prime(q) and evalint(P,x,q)==0 and len(row)==m
        for j,a in enumerate(row):
            assert (x+j)%q and pow((x+j)%q,(q-1)//2,q)==(1 if a==0 else q-1)
            count+=1
        rows.append(row)
    rank=len(echelon([bits(row) for row in rows]));ratiorank=len(echelon([bits([a^row[0] for a in row[1:]]) for row in rows]))
    assert m-len(echelon(masks))==rank==cert['neutral_rank']
    assert ratiorank==cert['equal_phase_rank']==rank and any(v.bit_count()%2 for v in masks)
    for rel in cert['search_relations']:
        assert product(F,alpha,rel['left'])==scale(F,product(F,alpha,rel['right']),rel['coefficient'])
    return count,len(masks)

def verify_number_algebra(cert,field='normal'):
    P=cert['P'];n=cert['n'];m=cert['m'];F=Field([int(x) for x in cert['field_polynomial_ascending']]);R=[F.element(x) for x in cert['roots_power_basis']]
    assert len(R)==n and all(F.evaluate(list(reversed(P)),r)==F.zero for r in R)
    if field=='normal':
        coeff=[F.one]
        for r in R:
            z=[F.zero]*(len(coeff)+1)
            for j,c in enumerate(coeff):z[j]=add(F,z[j],scale(F,F.multiply(c,r),-1));z[j+1]=add(F,z[j+1],c)
            coeff=z
        assert coeff==[F.integer(a) for a in reversed(P)]
        assert F.d==(6 if n==3 else 24)
        # The S_n certificate and this complete, distinct factorization imply
        # that this degree n! quotient algebra is the splitting field.
    else:
        assert F.d==12 and len(R)==2
        delta=add(F,R[0],scale(F,R[1],-1))
        assert F.multiply(delta,F.element(cert['root_difference_inverse_power_basis']))==F.one
        # S4 makes the cofactor cubic irreducible over a marked root field.
    rows=[];count=0;D=discriminant(P)
    for q,roots,row in cert['samples']:
        assert prime(q) and D%q and len(roots)==n and len(set(roots))==n and len(row)==n*m
        assert all(evalint(P,x,q)==0 for x in roots)
        for i,x in enumerate(roots):
            for j in range(m):
                a=row[i*m+j];assert (x+j)%q and a in (0,1)
                assert pow((x+j)%q,(q-1)//2,q)==(1 if a==0 else q-1);count+=1
        rows.append(row)
    masks=[]
    for rel in cert['relations']:
        z=F.one
        for i in rel['factor_indexes']:z=F.multiply(z,F.add_integer(R[i//m],i%m))
        assert F.power(F.element(rel['square_root_power_basis']),2)==z
        masks.append(sum(1<<j for j in rel['factor_indexes']))
    rank=len(echelon([bits(row) for row in rows]));assert rank==cert['rank'] and rank+len(echelon(masks))==n*m
    for result in cert['subset_orders']:
        I=result['root_subset'];k=len(echelon([bits([row[i*m+j] for i in I for j in range(m)]) for row in rows]))
        assert k==result['log2'] and result['order']==2**k
    return F,R,rows,count,len(masks)

def verify_mod4(cert,galois):
    P=cert['polynomial'];m=cert['M'];degree=cert['degree'];F=Field([int(x) for x in cert['FIELD']]);alpha=F.element(cert['ALPHA']);iota=F.element(cert['iota_power_basis'])
    assert F.evaluate(list(reversed(P)),alpha)==F.zero and F.power(iota,2)==F.integer(-1)
    assert F.d==2*degree
    if degree==2:assert cert['DISCRIMINANT']>0 and galois=='C2'
    # A real quadratic field or an odd degree field cannot contain i.
    # An S4 marked-root field has stabilizer S3 and no quadratic subfield.
    # Thus alpha and i generate a field of the entire quotient degree 2*d.
    p=cert['p'];w=cert['witness'];assert prime(p) and p%4==1 and evalint(P,w,p)==0
    assert all(pow((w+j)%p,(p-1)//4,p)==cert['phase'] for j in range(m))
    rows=[];count=m
    for q,x,ii,row in cert['samples']:
        assert prime(q) and q%4==1 and cert['DISCRIMINANT']%q and P[0]%q
        assert evalint(P,x,q)==0 and ii*ii%q==q-1 and len(row)==m
        for j,a in enumerate(row):assert a in range(4) and (x+j)%q and pow((x+j)%q,(q-1)//4,q)==pow(ii,a,q);count+=1
        rows.append(row)
    identities=0
    for kind in ['neutral','ratios']:
        M=rows if kind=='neutral' else [[(a-row[0])%4 for a in row[1:]] for row in rows]
        group=kernel4(M);stored=cert['groups'][kind]
        assert all(group[k]==stored[k] for k in ['C4','C2','log2_order','order'])
        assert [r['exponents'] for r in stored['relations']]==group['kernel_generators']
        for rel in stored['relations']:
            v=rel['exponents'];target=F.one
            for j,a in enumerate(v):target=F.multiply(target,F.power(F.add_integer(alpha,j if kind=='neutral' else j+1),a))
            left=F.power(F.element(rel['fourth_root_power_basis']),4)
            if kind=='ratios':left=F.multiply(left,F.power(alpha,sum(v)))
            assert left==target;identities+=1
    assert cert['phase_order_allowed']==cert['groups']['neutral']['order']//cert['groups']['ratios']['order']
    assert pow(cert['phase'],cert['phase_order_allowed'],p)==1
    for rel in cert['search_relations']:
        assert scale(F,product(F,alpha,rel['left']),rel['a'])==scale(F,product(F,alpha,rel['right']),rel['b'])
    return count,identities

def main(root,report_path=None):
    began=time.perf_counter();inputs=read(root,'INPUTS.json');galois={}
    for cert in read(root,'GALOIS_ELEMENTARY_CERTIFICATES.json'):galois[cert['tag']]=verify_galois(cert)
    character_values=0;squares=0;fourth_powers=0;objects={}
    for kind in ['cubic','quartic']:
        count,k=verify_base(inputs[kind]);character_values+=count;squares+=k
        c=read(root,kind+'_normal_certificate.json');F,R,rows,count,k=verify_number_algebra(c);character_values+=count;squares+=k;objects[kind]=(c,F,R,rows)
    c=read(root,'quartic_pair_certificate.json');_,_,_,count,k=verify_number_algebra(c,'pair');character_values+=count;squares+=k
    simple=read(root,'CUBIC_SIMPLIFICATION.json');cb=inputs['cubic'];F=Field(list(reversed(cb['polynomial'])));a=F.element([0,1])
    target=product(F,a,simple['factored_square_positions'])
    R=scale(F,product(F,a,simple['factored_square_root']['numerator_positions']),4845)
    assert F.multiply(target,F.power(F.add_integer(a,43),2))==F.power(R,2)
    assert simple['coefficient_product']==4845**2==255*285*323
    common=read(root,'CUBIC_COMMON_RATIONAL_CLASSES.json');c,F,R,rows=objects['cubic'];cols=[]
    for rel in common['rational_common_classes']:
        z=scale(F,product(F,R[0],rel['positions']),rel['constant'])
        assert F.power(F.element(rel['square_root_power_basis']),2)==z;squares+=1
        col=sum(int(pow(rel['constant'],(s[0]-1)//2,s[0])!=1)<<i for i,s in enumerate(c['samples']));cols.append(col)
    assert len(echelon(cols))==8
    density=read(root,'EXACT_DENSITIES.json')
    for kind in ['cubic','quartic']:
        c,F,R,rows=objects[kind];n=c['n'];m=c['m'];sp=[[sum(row[i*m+j]<<t for t,row in enumerate(rows)) for j in range(m)] for i in range(n)]
        for result in density[kind]['intersections']:
            I=result['root_subset'];V=sp[I[0]]
            for i in I[1:]:V=intersection(V,sp[i])
            assert len(echelon(V))==result['intersection_dimension']
        delta=(Fraction(1,2**51)-Fraction(1,2**95)+Fraction(1,3*2**122) if kind=='cubic' else Fraction(1,2**44)-Fraction(1,2**86)+Fraction(1,3*2**127)-Fraction(1,3*2**162))
        assert density[kind]['density_fraction']==[delta.numerator,delta.denominator]
        ranks={len(v['root_subset']):v['log2'] for v in c['subset_orders']};co=density[kind]['all_roots_common_dimension']
        assert density[kind]['higher_order_relation_dimension_after_common_quotient']==n*(ranks[1]-co)-(ranks[n]-co)
    quartic=inputs['quartic'];c,F,R,rows=objects['quartic'];classes=[2,3,14]
    for coeff,const in [(32,2),(192,3),(224,14)]:
        relation=next(r for r in quartic['search_relations'] if r['coefficient']==coeff);assert Fraction(coeff,const).denominator==1 and math.isqrt(coeff//const)**2==coeff//const
    assert len(echelon([bits([int(pow(a,(s[0]-1)//2,s[0])!=1) for a in classes]) for s in c['samples']]))==3
    mod4=[]
    for path in sorted(root.glob('r4_*_certificate.json')):
        cert=read(root,path.name);count,k=verify_mod4(cert,galois[f'r4_{cert["p"]}']);character_values+=count;fourth_powers+=k
        mod4.append({'p':cert['p'],'polynomial':cert['polynomial'],'m':cert['M'],'neutral_log2_order':cert['groups']['neutral']['log2_order'],'equal_phase_log2_order':cert['groups']['ratios']['log2_order'],'saved_bits':2*(cert['M']-1)-cert['groups']['ratios']['log2_order']})
    result={'status':'ALL_EXACT_CHECKS_PASS','standard_python_only':True,'exact_square_identities':squares,'exact_fourth_power_identities':fourth_powers,
        'auxiliary_and_origin_character_values_checked':character_values,'elementary_galois_certificates':len(galois),'r4_models':mod4,
        'density_formulas_verified':True,'class_group_completeness_assumed':False,'seconds':time.perf_counter()-began}
    if report_path is not None:Path(report_path).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('directory',nargs='?',default=str(Path(__file__).resolve().parent.parent));p.add_argument('--report')
    args=p.parse_args();main(Path(args.directory).resolve(),args.report)
