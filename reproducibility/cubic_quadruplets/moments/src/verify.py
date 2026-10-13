"""Finite checks only: the analytic theorems in the report are external inputs.
Uses exact integers/Fraction/Z[omega], no numerical quadrature or statistical fit.
"""
from pathlib import Path
from fractions import Fraction
from collections import Counter
from itertools import product
from math import isqrt
import json, sys
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'input'))
from elliptic_inherited import Curve, verify_trace, prime_exact
from jacobi_sector_inherited import jacobi_fast, ROOTS, mul, conj, twice_real
from moments import E_CURVES, C_CURVES, moment_table, route, compositions, one_moment

COUNTS = Counter()

def check(test, label):
    if not test:
        raise ArithmeticError(label)
    COUNTS[label] += 1

def frac(x):
    q = Fraction(x)
    return {'numerator': q.numerator, 'denominator': q.denominator}

def primes(n):
    a = bytearray(b'\x01') * (n + 1)
    a[:2] = b'\x00\x00'
    for d in range(2, isqrt(n) + 1):
        if a[d]:
            a[d*d::d] = b'\x00' * (((n-d*d)//d) + 1)
    return [q for q in range(2, n + 1) if a[q]]

def trace_count(p, a, b):
    # Counting square fibers, independent of the inherited pow-based point counter.
    sq = [0] * p
    for y in range(p):
        sq[y*y % p] += 1
    n = 1 + sum(sq[(x*x*x+a*x+b) % p] for x in range(p))
    return p + 1 - n

def cube_data(p):
    J, zeta, _ = jacobi_fast(p, check_prime=False)
    lut = {1: ROOTS[0], zeta: ROOTS[1], zeta*zeta % p: ROOTS[2]}
    ch = [(0, 0)] + [lut[pow(x, (p-1)//3, p)] for x in range(1, p)]
    return J, ch

def cm_from_jacobi(J, chis):
    return [-twice_real(mul(J, conj(chis[D]))) for D in (2,6,18)]

def sym_unscaled(p, t, m):
    a, b = 1, t
    if m == 0:
        return a
    for _ in range(2, m+1):
        a, b = b, t*b - p*a
    return b

def squarefree_part(n):
    sign = -1 if n < 0 else 1
    n = abs(n)
    r = 1
    for d in primes(isqrt(n)+1):
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e % 2:
            r *= d
    return sign*r*n

def small_checks():
    rows = []
    # p=7,13 require reducing the small D=18 modulo p explicitly.
    for p in primes(3000):
        if p <= 3 or p % 3 != 1:
            continue
        COUNTS['small_natural_primes'] += 1
        J, ch = cube_data(p)
        jc = [-twice_real(mul(J, conj(ch[D % p]))) for D in (2,6,18)]
        cm = [trace_count(p, a, b) for a,b in C_CURVES]
        ts = [trace_count(p, a, b) for a,b in E_CURVES]
        for (a,b), t in zip(C_CURVES+E_CURVES, cm+ts):
            check(t == p+1-Curve(p,a,b).literal_order(), 'two_independent_point_counts')
        for j, D in enumerate((2,6,18)):
            check(cm[j] == jc[j], 'CM_Jacobi_trace_identities')
            cubes = [0]*p
            for u in range(p):
                cubes[u*u*u%p] += 1
            nfermat = 3+sum(cubes[r]*cubes[(D-r)%p] for r in range(p))
            check(p+1-nfermat == cm[j], 'Fermat_fiber_vs_elliptic_trace')
            if p <= 100:
                for u in range(p):
                    for v in range(p):
                        if (u**3+v**3-D)%p or (u+v)%p == 0:
                            continue
                        s = pow((u+v)%p,-1,p)
                        x,y = 12*D*s%p,36*D*(u-v)*s%p
                        check((y*y-x**3+432*D*D)%p == 0, 'Fermat_birational_map_points')
        cube3 = ch[3] == (1,0)
        check(cm[0] == cm[1] == cm[2] if cube3 else sum(cm)==0, 'exact_cubic_twist_relations')
        I2 = int(ch[2] == (1,0))
        I23 = I2*int(cube3)
        c = 26+18*I2+18*I23
        beta = 2 if cube3 else -1
        num = p-c-2*cm[0]-2*cm[1]-2*ts[0]-2*ts[1]-beta*ts[2]
        direct = sum(ch[x]==ch[x+1]==ch[x+2]==ch[x+3] for x in range(1,p-3))
        check(num == 27*direct, 'new_A4_vs_natural_word')
        lhs = num-p+c
        rhs = (-4*cm[0]-2*sum(ts)) if cube3 else (2*cm[2]-2*ts[0]-2*ts[1]+ts[2])
        check(lhs == rhs, 'conditional_trace_decompositions')
        rows.append({'p':p,'cube3':cube3,'CM_traces':cm,'E_traces':ts,'A4':direct,
                     'corrected_numerator':lhs})
    return rows

def arithmetic_witnesses():
    proofs = []
    witnesses = ((5,7),(5,7),(7,19))
    for i, ((a,b), ps) in enumerate(zip(E_CURVES, witnesses)):
        vals = []
        for p in ps:
            t = trace_count(p,a,b)
            check(t % p != 0, 'ordinary_reduction_witnesses')
            vals.append({'p':p, 'trace':t, 'Frobenius_field_squarefree':squarefree_part(t*t-4*p)})
        check(vals[0]['Frobenius_field_squarefree'] != vals[1]['Frobenius_field_squarefree'],
              'incompatible_CM_field_witnesses')
        j = Fraction(1728*4*a**3,4*a**3+27*b*b)
        disc = -16*(4*a**3+27*b*b)
        q = 5 if i < 2 else 11
        t = trace_count(q,a,b)
        powers = [sym_unscaled(q,t,m) for m in range(1,5)]
        check(all(v != 0 for v in powers) and q%3==2, 'no_quadratic_basechange_self_twist_witnesses')
        proofs.append({'curve_index':i,'a':a,'b':b,'j_invariant':frac(j),
                       'Weierstrass_discriminant':disc,'ordinary_witnesses':vals,
                       'inert_prime':q,'symmetric_power_integer_traces':powers})
    traces31 = [trace_count(31,*v) for v in E_CURVES]
    check(traces31 == [-4,5,-7], 'trace_at_31')
    check(len(set(abs(t) for t in traces31)) == 3, 'no_GL2_twist_equivalences_witness')
    sym2 = [Fraction(t*t,31)-1 for t in traces31]
    check(len(set(abs(t) for t in sym2)) == 3, 'no_Sym2_cubic_twist_equivalences_witness')
    return {'curves':proofs,'traces31':traces31,'normalized_Sym2_at_31':[frac(v) for v in sym2]}

def large_checks():
    old = json.loads((ROOT/'input/A4_CUBIQUE_PRECEDENT.json').read_text())['cases']
    new = json.loads((ROOT/'results/CERTIFICATS_CM.json').read_text())['cases']
    rows=[]
    for r, s in zip(old, new):
        p=r['p']
        check(p == s['p'], 'large_case_alignment')
        J,zeta,_ = jacobi_fast(p,check_prime=False)
        lut={1:ROOTS[0],zeta:ROOTS[1],zeta*zeta%p:ROOTS[2]}
        ch={D:lut[pow(D,(p-1)//3,p)] for D in (2,3,6,18)}
        ts=[]
        for c in r['elliptic_certificates']:
            verify_trace(c)
            COUNTS['inherited_nonCM_point_order_certificates'] += 1
            ts.append(c['trace'])
        cm=[]
        for c in s['CM_certificates']:
            verify_trace(c)
            j=c['CM_index']
            a,b=C_CURVES[j]
            check(c['a']==a%p and c['b']==b%p,'large_CM_curve_parameters')
            cm.append(c['trace'])
            COUNTS['new_CM_point_order_certificates'] += 1
        check(cm == cm_from_jacobi(J,ch),'large_CM_vs_Jacobi')
        cube3=ch[3]==(1,0)
        I2=int(ch[2]==(1,0)); c=26+18*I2+18*I2*int(cube3)
        z=-2*cm[0]-2*cm[1]-2*ts[0]-2*ts[1]-(2 if cube3 else -1)*ts[2]
        check(27*r['A4_cubic']-p+c==z, 'large_new_A4_identity')
        rows.append({'p':p,'cube3':cube3,'CM_traces':cm,'E_traces':ts,
                     'A4':r['A4_cubic'],'correction':c,'corrected_numerator':z,
                     'corrected_square':frac(Fraction(z*z,p)),
                     'corrected_fourth':frac(Fraction(z**4,p*p))})
    return rows

def formal_checks():
    table=moment_table()
    check(table['cube']['moments_0_to_4']==[1,0,44,0,4224], 'cube_moments')
    check(table['noncube']['moments_0_to_4']==[1,0,17,0,738], 'noncube_moments')
    check(table['all']['moments_0_to_4']==[1,0,26,0,1900], 'mixture_moments')
    # Classify every relevant character monomial, not just the ones appearing in a fitted model.
    routes = Counter()
    for h in range(-4,5):
        for m in product(range(5), repeat=3):
            if abs(h)+sum(m)>4:
                continue
            for u in range(3):
                r=route(h,m,u)
                routes[r['type']] += 1
                COUNTS['analytic_case_bookkeeping'] += 1
    # Independent one-variable constant-term formulas, including every monomial of total degree <=4.
    mono=[]
    for d in range(5):
        for exps in compositions(d,4):
            val=1
            for j,r in enumerate(exps):
                laur={r-2*i:__import__('math').comb(r,i) for i in range(r+1)}
                own=laur.get(0,0)-(laur.get(-2,0) if j else 0)
                check(own==one_moment(r,j==0),'one_variable_Weyl_moments')
                val*=own
            mono.append({'powers':list(exps),'limit':val})
            COUNTS['joint_monomials_total_degree_at_most_four'] += 1
    return {'moment_table':table, 'routes':dict(routes), 'joint_monomials':mono,
            'unconditional_kurtosis':frac(Fraction(1900,26**2)),
            'variance_ratio_cube_over_noncube':frac(Fraction(44,17)),
            'covariance_CM_cube':[[2,2],[2,2]],
            'covariance_CM_noncube':[[2,-1],[-1,2]],
            'not_proved':'Full joint distribution; growing window length; law of maxima.'}

def negative_tests():
    for args in ((0,(2,2,2),0),(5,(0,0,0),0),(0,(-1,0,0),0),(0,(0,0,0),3)):
        try:
            route(*args)
        except ValueError:
            COUNTS['negative_domain_tests'] += 1
        else:
            raise ArithmeticError('expected domain rejection')
    cc=json.loads((ROOT/'results/CERTIFICATS_CM.json').read_text())['cases'][0]['CM_certificates'][0]
    for key in ('trace','group_order','point_order'):
        c=dict(cc); c[key]+=1
        try:
            verify_trace(c)
        except (ValueError,ArithmeticError):
            COUNTS['altered_certificates_rejected']+=1
        else:
            raise ArithmeticError('altered trace certificate accepted')

def main():
    small=small_checks()
    witnesses=arithmetic_witnesses()
    large=large_checks()
    algebra=formal_checks()
    negative_tests()
    objects={
        'CERTIFICATS_FINIS.json':{'status':'PASS','purpose':'exact identities, not statistical validation',
            'small_cases':small,'arithmetic_witnesses':witnesses,'large_cases':large},
        'MOMENTS_EXACTS.json':algebra,
        'VERIFICATION.json':{'status':'PASS','counts':dict(sorted(COUNTS.items())),
            'external_inputs_not_computationally_proved':['CM Hecke theory','symmetric-power automorphy',
             'cyclic base change','GL2 tensor product functoriality','Rankin-Selberg prime number theorem'],
            'no_GPU':True,'no_atlas_rescan':True,'no_full_maximum_law':True}}
    for name,obj in objects.items():
        (ROOT/'results'/name).write_text(json.dumps(obj,indent=2,sort_keys=True)+'\n',encoding='utf-8')
    print(json.dumps(objects['VERIFICATION.json'],indent=2,sort_keys=True))

if __name__ == '__main__':
    main()
