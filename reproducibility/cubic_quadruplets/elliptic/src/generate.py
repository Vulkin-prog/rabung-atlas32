"""Generate exact elliptic trace certificates. Bounded, local CPU computation.
No GPU, network, atlas data or global maxima are accessed.
"""
from pathlib import Path
import json
from fractions import Fraction
from core import aggregate,elliptic_coefficients,fourth_sector,frac_json
from elliptic import Curve,prime_exact
from a4 import a4_from_traces, PARAMETERS

ROOT=Path(__file__).resolve().parents[1]
CASES=[(1000033,13),(4294967197,21),(4294967311,21)]

def main():
    cases=[]; a4cases=[]
    for p,k in CASES:
        if not prime_exact(p): raise ValueError('p not prime')
        groups,meta=aggregate(k,p)
        certs=[];traces={}
        for A in sorted(groups):
            if not groups[A]:continue
            a,b=elliptic_coefficients(A,p)
            c=Curve(p,a,b).certify_order()
            c.update({'parameter':frac_json(A),'weight':groups[A]})
            certs.append(c);traces[A]=c['trace']
        a4certs=[]
        certmap={Fraction(c['parameter']['numerator'],c['parameter']['denominator']):c for c in certs}
        for A in PARAMETERS:
            c=certmap.get(A)
            if c is None:
                c=Curve(p,*elliptic_coefficients(A,p)).certify_order()
                c['parameter']=frac_json(A)
            a4certs.append({key:v for key,v in c.items() if key!='weight'})
        _,a4meta=a4_from_traces(p,{A:c['trace'] for A,c in zip(PARAMETERS,a4certs)})
        a4meta['elliptic_certificates']=a4certs
        a4cases.append(a4meta)
        value,meta=fourth_sector(k,p,traces)
        meta['trace_certificates']=certs
        cases.append(meta)
        print(p,k,'curves',len(certs),'packet',meta['elliptic_packet'],
              'holes',meta['holes_numerator'],'sector',value,float(value),flush=True)
    obj={'status':'exact certificates; not a maximum law',
         'cases_chosen_before_trace_generation':CASES,'cases':cases}
    path=ROOT/'results'/'CERTIFICATS_ELLIPTIQUES.json'
    path.write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
    (ROOT/'results'/'A4_CUBIQUE.json').write_text(json.dumps({'status':'exact cubic four-window counts, not maxima','cases':a4cases},sort_keys=True,indent=2)+'\n')

if __name__=='__main__':main()
