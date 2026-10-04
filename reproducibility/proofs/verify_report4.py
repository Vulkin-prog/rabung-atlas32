"""Exact independent checks of additions suggested by RAPPORT_FR(4).md.

Standard Python only; uses previously certified fields and finite characters.
No PARI, GPU, network or new global prime census is used.
"""
from pathlib import Path
from fractions import Fraction
import argparse, json, math, time, sys

def run(root):
    sys.path.insert(0,str(root/'proofs'))
    from exact_fields import Field, prime
    from algebra import add, scale, echelon, reduce, nullspace, bits
    from finite_modules import kernel4
    def read(p):return json.loads(p.read_text())
    began=time.perf_counter();quadratics=[]
    models=[(1583868241,19),(2400748661,889),(3318080401,97)]
    for p,D in models:
        cert=read(root/'evidence/new'/f'order4_pair_{p}_certificate.json')
        F=Field(cert['field_polynomial_ascending']);a,b=map(F.element,cert['roots_power_basis']);m=cert['m']
        P=cert['P'];base=read(root/'anteriorities/theorie_densites_residuels_20261004'/cert['base_certificate'])
        assert F.d==4 and math.isqrt(D)**2!=D
        assert F.evaluate(list(reversed(P)),a)==F.zero==F.evaluate(list(reversed(P)),b)
        disc=P[1]**2-4*P[0]*P[2];ds=math.isqrt(disc//D)
        assert disc==D*ds**2
        sqrtD=scale(F,F.add_integer(scale(F,a,2*P[0]),P[1]),Fraction(1,ds))
        assert F.power(sqrtD,2)==F.integer(D)
        assert F.power(F.element(base['iota_power_basis']),2)==F.integer(-1)
        A=lambda j:F.add_integer(a,j)
        B=lambda j:F.add_integer(b,j)
        mul=F.multiply
        identities=[]
        if D==19:
            assert mul(A(9),A(11))==scale(F,mul(A(12),A(26)),-1)
            assert mul(A(12),B(3))==mul(A(9),A(11))
            identities=['A9*A11=-A12*A26','A12*B3=A9*A11']
            assert b==F.add_integer(scale(F,a,-1),-29)
            P=cert['P'];values=[sum(x*(-j)**(len(P)-1-k) for k,x in enumerate(P)) for j in range(30)]
            assert all(values) and max(map(abs,values))==411
            phases=read(root/'evidence/new/ORDER4_PAIR_PHASES.json')
            phase=next(x for x in phases if x['p']==p)
            assert phase['possible_constant_phase_exponent_pairs']==[[0,0],[1,1],[2,2],[3,3]]
            assert cert['groups']['neutral']['order']==4**21
            assert cert['groups']['ratios']['order']==4**20
            # B_j/B_0=A_(29-j)/A_29. In particular A_29/A_0 is
            # (A_26/A_0)/(B_3/B_0); the two ratio groups coincide.
            assert mul(A(29),B(3))==mul(A(26),B(0))
            extended={'length':30,'max_abs_P_negative_positions':411,
                'neutral_group_order':4**21,'constant_group_order':4**20,
                'density_among_admissible_primes_neutral':[1,2**43],
                'density_among_admissible_primes_constant':[1,2**41],
                'density_among_all_primes_neutral':[1,2**44],
                'density_among_all_primes_constant':[1,2**42],
                'extension_given_constant_28':[1,31]}
        elif D==889:
            assert mul(A(0),B(16))==scale(F,A(8),16)
            identities=['A0*B16=16*A8'];extended=None
        else:
            assert mul(A(10),B(11))==A(9)
            identities=['A10*B11=A9'];extended=None
        case=next(x for x in read(root/'evidence/new/ORDER4_DENSITIES_EXACT.json')['cases'] if x['p_origin']==p)
        dens={mode:Fraction(v['density_among_all_primes']['numerator'],v['density_among_all_primes']['denominator']) for mode,v in case['modes'].items()}
        delta=dens['ratios']-dens['neutral']
        expected={19:Fraction(93,2**44),889:Fraction(3,2**35),97:Fraction(1,2**39)}[D]
        assert delta==expected
        if extended:
            assert Fraction(1,2**42)/dens['ratios']==Fraction(1,31)
        criterion=None
        if D==889:
            T=[1,2,3,4,5,6,7,8,9,10,13,19,21,24,26,27]
            base=read(root/'anteriorities/theorie_densites_residuels_20261004'/cert['base_certificate'])
            rows=[[(row[j]-row[0])%4 for j in T] for q,w,ii,row in base['samples']]
            g=kernel4(rows)
            assert (g['C4'],g['C2'],g['order'])==(16,0,4**16)
            assert g['order']==base['groups']['ratios']['order']
            criterion={'reference':0,'comparison_positions':T,'comparisons':16,'distinct_positions':17,
                'ratio_group_order':4**16,'minimal_among_individual_ratio_generators':True}
        P=cert['P'];window=30 if D==19 else m
        norm_values=[sum(x*(-j)**(len(P)-1-k) for k,x in enumerate(P)) for j in range(window)]
        assert all(norm_values) and max(map(abs,norm_values))<=1126
        # The true field is Q(sqrt(D),i); its ramified primes divide 2*D.
        # A prime >1126 also avoids the leading coefficient and every window
        # norm. Power roots of unit products are automatically local units.
        assert 2*P[0]<=1126 and D<=1126
        assert prime(p) and p%4==1
        quadratics.append({'p_origin':p,'D':D,'identities':identities,'extended_window':extended,
            'nonneutral_existence_density_all_primes':[delta.numerator,delta.denominator],
            'comparison_criterion':criterion,'sufficient_good_prime_bound':1126,
            'maximum_absolute_window_norm':max(map(abs,norm_values))})

    normal=[]
    for kind in ['cubic','quartic']:
        cert=read(root/'anteriorities/theorie_densites_residuels_20261004'/f'{kind}_normal_certificate.json')
        norms=read(root/'evidence/new'/f'{kind}_NORM_DEPENDENCIES.json')
        n,m=cert['n'],cert['m'];samples=cert['samples'];s=len(samples)
        C=[sum((pow(a,(q-1)//2,q)!=1)<<i for i,(q,_,_) in enumerate(samples)) for a in norms['common_rational_classes']]
        cp=echelon(C);blocks=[];selected=[]
        for i in range(n):
            values=[reduce(sum(row[i*m+j]<<k for k,(_,_,row) in enumerate(samples)),cp) for j in range(m)]
            basis=[];positions=[]
            for j,v in enumerate(values):
                if reduce(v,echelon(basis)):
                    basis.append(v);positions.append(j)
            blocks.append(basis);selected.append(positions)
        d=norms['single_quotient_dimension'];assert all(len(b)==d for b in blocks)
        columns=sum(blocks,[])
        rows=[sum(((v>>k)&1)<<j for j,v in enumerate(columns)) for k in range(s)]
        R=nullspace(rows,n*d);e=len(R)
        assert e==norms['higher_relation_dimension']==norms['norm_generated_higher_dimension']
        block_changes=[]
        for i in range(n):
            coupled=[(r>>(i*d))&((1<<d)-1) for r in R]
            assert len(echelon(coupled))==e
            complete=coupled[:]
            for k in range(d):
                if reduce(1<<k,echelon(complete)):complete.append(1<<k)
            assert len(complete)==d
            block_changes.append(complete)
        for r in R:
            z=0
            for k,v in enumerate(columns):
                if r>>k&1:z^=v
            assert z==0
        moments=[Fraction(1)]+[Fraction(math.comb(n,k),2**(k*d)) for k in range(1,n)]+[Fraction(1,2**(n*d-e))]
        probs=[sum((-1)**(j-k)*math.comb(j,k)*moments[j] for j in range(k,n+1)) for k in range(n+1)]
        assert all(v>=0 for v in probs) and sum(probs)==1
        normal.append({'family':kind,'roots':n,'block_dimension':d,'coupled_parities':e,
            'remaining_bits_after_other_roots':d-e,'common_quotient_basis_positions':selected,
            'block_coordinate_changes':block_changes,
            'joint_neutral_probability':[1,2**(n*d-e)],
            'independent_marginal_comparator':[1,2**(n*d)],
            'enhancement':2**e,
            'conditional_neutral_root_count_law':[[v.numerator,v.denominator] for v in probs],
            'scope':'Split normal-field primes with common classes neutral and a uniformly chosen place for root labels. Completeness uses the previously verified norm-generated relation certificate.'})
    return {'status':'PASS','quadratic_additions':quadratics,'parity_normal_forms':normal,
        'seconds':time.perf_counter()-began,'requirements':'Standard Python plus previously verified exact proof bundle.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    p.add_argument('--report',type=Path);a=p.parse_args();result=run(a.root)
    if a.report:a.report.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'quadratics':len(result['quadratic_additions']),
        'normal_forms':[{k:v[k] for k in ['family','roots','block_dimension','coupled_parities','remaining_bits_after_other_roots']} for v in result['parity_normal_forms']],
        'seconds':result['seconds']}))
