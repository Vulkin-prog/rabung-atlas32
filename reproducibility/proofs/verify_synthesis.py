"""Portable verifier of this supplement: exact standard-library arithmetic.

Default layout: this script in proofs/; data in its parent; anteriorities/
contains the four earlier sealed theoretical deliveries. No output is written
unless --report is supplied. No PARI, GPU, network or prime census is needed.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import argparse,json,math,time,itertools,csv,copy
from exact_fields import Field,prime
from finite_modules import kernel4,intersection
from algebra import echelon,reduce,coordinates,nullspace,bits,add,scale

@lru_cache(None)
def isprime(q):return prime(q)
def read(path):return json.loads(path.read_text())
def evaluate(P,x,q=None):
    z=0
    for a in P:z=z*x+a;z=z if q is None else z%q
    return z
def transpose(C,s):return [sum(((v>>i)&1)<<j for j,v in enumerate(C)) for i in range(s)]
def field(P):return Field([int(x) for x in P])

def check_multi(d,old,unit=None):
    F=field(d['field_polynomial_ascending']);R=[F.element(v) for v in d['roots_power_basis']]
    P=d['P'];m=d['m'];n=len(R);degree=len(P)-1
    assert all(F.evaluate(list(reversed(P)),r)==F.zero for r in R)
    if unit:
        iota=F.element(unit['iota_power_basis'])
        inverse=F.element(unit['root_difference_inverse_power_basis'])
        assert F.multiply(add(F,R[0],scale(F,R[1],-1)),inverse)==F.one
    else:iota=F.element(old['iota_power_basis'])
    assert F.power(iota,2)==F.integer(-1)
    if degree>=3:
        D=old['DISCRIMINANT'];assert not(D<0 and math.isqrt(-D)**2==-D)
    if n==degree:
        coeff=[F.one]
        for r in R:
            z=[F.zero]*(len(coeff)+1)
            for j,c in enumerate(coeff):
                z[j]=add(F,z[j],scale(F,F.multiply(c,r),-1));z[j+1]=add(F,z[j+1],c)
            coeff=z
        assert coeff==[scale(F,F.integer(a),Fraction(1,P[0])) for a in reversed(P)]
    assert F.d==({2:4,3:12,4:48}[degree] if n==degree else 24)
    # Earlier elementary Galois proofs establish S_d. All roots and i in a
    # degree 2*d! algebra force the splitting field compositum. S3 and S4
    # have only their discriminant quadratic subfield, distinct from Q(i).
    # For pairs the inverse above distinguishes both roots: S4 transitivity
    # forces root-field degree 12 and adjoining i doubles it.
    rows=[];values=0
    for q,roots,ii,row in d['samples']:
        assert isprime(q) and q%4==1 and old['DISCRIMINANT']%q and P[0]%q
        assert ii*ii%q==q-1 and len(row)==n*m and len(set(roots))==n
        assert all(evaluate(P,w,q)==0 for w in roots)
        for i,w in enumerate(roots):
            for j in range(m):
                e=row[i*m+j];assert e in range(4) and (w+j)%q
                assert pow((w+j)%q,(q-1)//4,q)==pow(ii,e,q);values+=1
        rows.append(row)
    identities=0
    for mode in ['neutral','ratios']:
        L=m if mode=='neutral' else m-1
        M=rows if mode=='neutral' else [[(row[i*m+j]-row[i*m])%4 for i in range(n) for j in range(1,m)] for row in rows]
        lower=kernel4(M);stored=d['groups'][mode]
        assert all(lower[k]==stored[k] for k in ['C4','C2','log2_order','order'])
        assert lower['kernel_generators']==[x['exponents'] for x in stored['relations']]
        for rel in stored['relations']:
            z=F.one;den=F.one
            for i,e in enumerate(rel['exponents']):
                if not e:continue
                r=R[i//L];j=i%L+(mode=='ratios')
                z=F.multiply(z,F.power(F.add_integer(r,j),e))
                if mode=='ratios':den=F.multiply(den,F.power(r,e))
            left=F.power(F.element(rel['fourth_root_power_basis']),4)
            assert F.multiply(left,den)==z;identities+=1
        for subset in stored.get('subset_groups',[]):
            I=subset['root_subset'];g=kernel4([[row[i*L+j] for i in I for j in range(L)] for row in M])
            assert all(g[k]==subset[k] for k in ['C4','C2','order','log2_order'])
    return values,identities

def check_norms(root,previous):
    results={};identities=0
    for kind in ['cubic','quartic']:
        c=read(previous/'theorie_densites_residuels_20261004'/f'{kind}_normal_certificate.json')
        d=read(root/f'{kind}_NORM_DEPENDENCIES.json');short=read(root/f'{kind}_SHORT_NORMS.json')
        F=field(c['field_polynomial_ascending']);alpha=F.element(c['roots_power_basis'][0])
        for rel in read(root/'COMMON_CLASSES_EXACT.json')[kind]:
            z=F.integer(rel['constant'])
            for j in rel['positions']:z=F.multiply(z,F.add_integer(alpha,j))
            assert F.power(F.element(rel['square_root_power_basis']),2)==z;identities+=1
        n,m=c['n'],c['m'];rows=[s[2] for s in c['samples']];s=len(rows)
        C=[sum(int(pow(a,(q-1)//2,q)!=1)<<t for t,(q,_,_) in enumerate(c['samples'])) for a in d['common_rational_classes']]
        cp=echelon(C);assert len(cp)==d['common_dimension']
        columns=[sum(row[j]<<t for t,row in enumerate(rows)) for j in range(n*m)]
        quotient=[[reduce(columns[i*m+j],cp) for j in range(m)] for i in range(n)]
        qdim=len(echelon(quotient[0]));assert qdim==d['single_quotient_dimension']
        for k in range(1,n):
            for I in itertools.combinations(range(n),k):
                assert len(echelon(sum([quotient[i] for i in I],[])))==k*qdim
        relation_basis=[bits(v['factor_indexes']) for v in c['relations']];rp=echelon(relation_basis)
        single=nullspace(transpose(columns[:m],s),m)
        internal=[v<<(i*m) for i in range(n) for v in single]
        common_positions=d['common_class_positions_first_root']
        for i in range(n):
            for a,I in zip(C,common_positions):
                v=0
                for j in I:v^=columns[i*m+j]
                assert v==a
        for i in range(1,n):
            for I in common_positions:internal.append(bits(I)^bits([i*m+j for j in I]))
        assert all(reduce(v,rp)==0 for v in internal)
        generated=list(internal)
        for rel in d['norm_constraints']:
            v=bits(rel['square_factor_indexes']);assert reduce(v,rp)==0
            z=0
            for i in rel['archived_square_relation_combination']:z^=relation_basis[i]
            assert z==v;generated.append(v)
        assert len(echelon(generated))==len(rp)==d['total_square_relation_dimension']
        higher=n*qdim-len(echelon(sum(quotient,[])))
        assert higher==d['higher_relation_dimension']==d['norm_generated_higher_dimension']
        norms=[(-1)**n*evaluate(c['P'],-j) for j in range(m)]
        assert norms==short['norms'] and all(norms)
        shortvectors=[];actual=[]
        for rel in short['short_norm_relation_basis']:
            v=bits(rel['positions']);shortvectors.append(v)
            assert math.prod(norms[j] for j in rel['positions'])==rel['norm_product']
            assert math.prod(rel['constant_factors'])==rel['constant_product']
            assert rel['norm_product']*rel['constant_product']==rel['integer_square_root']**2
            w=0
            for i in range(n):w^=v<<(i*m)
            for a in rel['constant_factors']:
                I=common_positions[d['common_rational_classes'].index(a)];w^=bits(I)
            assert reduce(w,rp)==0;actual.append(w);identities+=1
        assert len(echelon(internal+actual))-len(echelon(internal))==higher
        results[kind]={'common_dimension':len(cp),'single_quotient_dimension':qdim,'higher_norm_relations':higher,'short_weights':[v.bit_count() for v in shortvectors]}
    return results,identities

def check_local(root,previous):
    cert=read(previous/'theorie_prolongements_20261004'/'CUBIQUE_53_RELATIONS.json');law=read(root/'LOCAL_CODE_LAW.json')
    m=cert['m'];R=[bits(x['positions']) for x in cert['relations']];F=Field(list(reversed(cert['polynomial'])));a=F.element([0,1,0])
    for rel in cert['relations']:
        z=F.one
        for j in rel['positions']:z=F.multiply(z,F.add_integer(a,j))
        assert F.power(F.element(rel['square_root_power_basis']),2)==z
    rows=[]
    for q,w,row in cert['lower_character_samples']:
        assert isprime(q) and evaluate(cert['polynomial'],w,q)==0
        assert all((w+j)%q and pow((w+j)%q,(q-1)//2,q)==(1 if e==0 else q-1) for j,e in enumerate(row))
        rows.append(bits([j for j,e in enumerate(row) if e]))
    assert len(echelon(rows))==41 and len(echelon(R))==12
    dual=[0]
    for r in R:dual +=[x^r for x in dual]
    hist={str(w):sum(v.bit_count()==w for v in dual) for w in sorted({v.bit_count() for v in dual})}
    assert hist==law['dual_weight_distribution'] and law['minimum_dual_weight']==4
    weights=[]
    for k in range(m+1):
        s=sum(num*sum((-1)**j*math.comb(int(w),j)*math.comb(m-int(w),k-j) for j in range(max(0,k-(m-int(w))),min(int(w),k)+1)) for w,num in hist.items())
        assert s%4096==0;weights.append(s//4096)
    assert weights==law['code_weight_counts'] and sum(weights)==2**41
    tails=list(csv.DictReader((root/'LOCAL_MAXIMUM_CUBIC.csv').open()))
    assert len(tails)==54
    for row in tails:
        assert row['below_code_fourier']==row['below_code_syndrome']
    assert int(tails[-1]['below_code_fourier'])==2**41
    assert int(tails[-2]['below_code_fourier'])==2**41-2==int(tails[-3]['below_code_fourier'])
    return {'dimension':41,'relations':12,'dual_minimum_distance':4,'checked_maximum_thresholds':54,'monochromatic_enhancement':4096}

def check_residual(root,previous):
    total=0;values=0
    for degree in range(2,7):
        data=read(root/f'RESIDUAL_QUALIFICATION_D{degree}.json')
        for model in data['models']:
            d=read(root/model['certificate']);m=d['M'];rows=[]
            assert d['p']==2940141229 and m==28 and isprime(d['p']) and evaluate(d['polynomial'],d['witness'],d['p'])==0
            assert all(pow((d['witness']+j)%d['p'],(d['p']-1)//4,d['p'])==d['phase'] for j in range(m))
            for q,w,ii,row in d['samples']:
                assert isprime(q) and q%4==1 and d['DISCRIMINANT']%q and d['polynomial'][0]%q
                assert evaluate(d['polynomial'],w,q)==0 and ii*ii%q==q-1 and len(row)==m and all(e in range(4) for e in row)
                assert all((w+j)%q and pow((w+j)%q,(q-1)//4,q)==pow(ii,e,q) for j,e in enumerate(row));values+=m;rows.append(row)
            neutral=kernel4(rows);ratios=kernel4([[(e-row[0])%4 for e in row[1:]] for row in rows])
            assert (neutral['C4'],neutral['C2'],ratios['C4'],ratios['C2'])==(28,0,27,0)
            # An image of maximal possible order proves no fourth-power
            # relation: no unverified high-degree field assumption is needed.
            total+=1
    assert total==100
    return total,values

def check_densities(root,previous):
    cases=read(root/'ORDER4_DENSITIES_EXACT.json')['cases'];assert len(cases)==8
    for case in cases:
        d=case['degree'];p=case['p_origin'];old=read(previous/'theorie_densites_residuels_20261004'/f'r4_{p}_d{d}_certificate.json')
        normal=read(root/(f'order4_pair_{p}_certificate.json' if d==2 else f'order4_{p}_normal_certificate.json'))
        pair=normal if d<=3 else read(root/f'order4_{p}_pair_certificate.json')
        for mode,data in case['modes'].items():
            orders=[old['groups'][mode]['order']]
            for k in range(2,d+1):
                cert=pair if k==2 else normal
                order=cert['groups'][mode]['order'] if d==2 else next(v['order'] for v in cert['groups'][mode]['subset_groups'] if len(v['root_subset'])==k)
                orders.append(order)
            assert orders==[v['order'] for v in data['subset_groups_by_size']]
            moments=[Fraction(1)]+[Fraction(1,2*math.factorial(k)*order) for k,order in enumerate(orders,1)]
            probabilities=[sum((-1)**(j-k)*math.comb(j,k)*moments[j] for j in range(k,d+1)) for k in range(d+1)]
            stored=[Fraction(x['numerator'],x['denominator']) for x in data['root_count_law_among_all_primes']]
            assert probabilities==stored and all(x>=0 for x in probabilities) and sum(probabilities)==1
            assert 1-probabilities[0]==Fraction(data['density_among_all_primes']['numerator'],data['density_among_all_primes']['denominator'])
    for item in read(root/'ORDER4_PAIR_PHASES.json'):
        c=read(root/f"order4_pair_{item['p']}_certificate.json");m=c['m'];gens=[v['exponents'] for v in c['groups']['neutral']['relations']]
        phases=[list(v) for v in itertools.product(range(4),repeat=2) if all(sum(e*v[i//m] for i,e in enumerate(g))%4==0 for g in gens)]
        assert phases==item['possible_constant_phase_exponent_pairs']
        assert len(phases)==c['groups']['neutral']['order']//c['groups']['ratios']['order']
    shifted=read(root/'DIHEDRAL_SHIFTED_63.json');P=shifted['P']
    assert max(abs(evaluate(P,-j)) for j in range(-11,52))==813402
    assert set(range(52))|set(range(-11,41))==set(range(-11,52))
    return {'order_four_density_laws':16,'quadratic_phase_pair_laws':3,'dihedral_reflected_length':63}

def check_small_arithmetic(root):
    d=read(root/'RESIDUAL_LOCAL_CONTEXT.json');p=d['p'];w=d['witness'];m=d['m'];ii=d['iota']
    assert isprime(p) and p%8==5 and ii*ii%p==p-1 and m==28
    assert pow(p-1,(p-1)//4,p)==d['chi4_minus_one']==p-1
    start,end=d['positions'];assert len(d['local_quartic_exponents'])==end-start+1
    for j,e in zip(range(start,end+1),d['local_quartic_exponents']):
        assert e in range(4) and pow((w+j)%p,(p-1)//4,p)==pow(ii,e,p)
    assert all(pow(w+j,(p-1)//4,p)==p-1 for j in range(m))
    reflected=p-w-m+1
    assert reflected==d['reflected_start'] and d['reflected_phase_exponent']==0
    assert all(pow(reflected+j,(p-1)//4,p)==1 for j in range(m))
    run=d['quadratic_neutral_local_run'];assert run['start']==w-1 and run['length']==29
    assert run['left_extra']==1 and run['right_extra']==0 and run['exact_edges_seen']
    a=run['start'];b=a+run['length']
    assert all(pow(x,(p-1)//2,p)==1 for x in range(a,b))
    assert pow(a-1,(p-1)//2,p)==p-1 and pow(b,(p-1)//2,p)==p-1
    checked=[]
    for q in range(5,258,2):
        if not isprime(q):continue
        def chi(x):
            y=pow(x%q,(q-1)//2,q)
            return 0 if y==0 else (1 if y==1 else -1)
        direct=sum(len({chi(x+j) for j in range(4)})==1 for x in range(1,q-3))
        B=4*chi(2)+2*chi(6)+2*chi(3)+2*chi(-1)+2*chi(-2)
        Epoints=q+1+sum(chi(x*(x-1)*(x-9)) for x in range(q))
        trace=q+1-Epoints
        assert 8*direct==q-11-B-trace
        assert sum(chi(x*(x+1)*(x+2)*(x+3)) for x in range(q))==-trace-1
        checked.append(q)
    return {'residual_local_context':'PASS','direct_A4_and_elliptic_checks':len(checked),'small_prime_maximum':max(checked)}

def main(root,previous,report=None):
    began=time.perf_counter();units={x['field_file']:x for x in read(root/'FIELD_UNIT_CERTIFICATES.json')}
    values=0;powers=0;files=[]
    for path in sorted(root.glob('order4_*_certificate.json')):
        d=read(path);old=read(previous/'theorie_densites_residuels_20261004'/d['base_certificate'])
        stem=path.name.replace('_certificate.json','_field.json');u=units.get(stem)
        v,p=check_multi(d,old,u);values+=v;powers+=p;files.append(path.name)
        print(json.dumps({'verified':path.name,'fourth_powers':p}),flush=True)
    assert len(files)==10
    norms,norm_identities=check_norms(root,previous)
    local=check_local(root,previous)
    total,extra=check_residual(root,previous);values+=extra
    density=check_densities(root,previous)
    arithmetic=check_small_arithmetic(root)
    result={'status':'PASS','multi_root_certificates':len(files),'fourth_power_identities':powers,'checked_character_values':values,
        'norm_integer_square_identities':norm_identities,'norms':norms,'local_code':local,'residual_models_with_full_character_group':total,
        'density_laws':density,'small_arithmetic_checks':arithmetic,'seconds':time.perf_counter()-began,'requirements':'Standard Python library only; previous sealed proof bundle.'}
    if report:report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True);return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent/'evidence/new')
    p.add_argument('--previous',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
    main(a.root,a.previous or Path(__file__).resolve().parent.parent/'anteriorities',a.report)
