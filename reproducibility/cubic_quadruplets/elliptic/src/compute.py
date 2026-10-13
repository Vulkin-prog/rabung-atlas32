"""Compute and certify cubic A4; optionally the exact fourth-support start sector.
This is NOT a maximum scanner. All output contains auditable certificates.
Implemented prime domain: 3<p<10^10, p=1 modulo3. The point search may refuse.
"""
import argparse,json,sys
from core import aggregate,elliptic_coefficients,fourth_sector,frac_json
from a4 import PARAMETERS,a4_from_traces
from elliptic import Curve,prime_exact

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--p',type=int,required=True)
    ap.add_argument('--k',type=int,help='also compute fourth-support start sector (4<=k<=32)')
    args=ap.parse_args();p=args.p;k=args.k
    if not prime_exact(p) or p<=3 or p%3!=1:raise ValueError('cubic prime p<10^10 required')
    if k is not None and not 4<=k<=min(32,p-2):raise ValueError('declared k budget is 4..min(32,p-2)')
    weights={} if k is None else aggregate(k,p)[0]
    parameters=set(PARAMETERS)|{A for A,w in weights.items() if w}
    certs=[];traces={}
    for A in sorted(parameters):
        cert=Curve(p,*elliptic_coefficients(A,p)).certify_order()
        cert['parameter']=frac_json(A)
        certs.append(cert);traces[A]=cert['trace']
    _,a4=a4_from_traces(p,traces)
    out={'p':p,'A4_cubic':a4,'trace_certificates':certs,'not_computed':'global maximum or critical law'}
    if k is not None:out['fourth_support_start_sector']=fourth_sector(k,p,traces)[1]
    print(json.dumps(out,sort_keys=True,indent=2))

if __name__=='__main__':
    try:main()
    except (ValueError,RuntimeError,ArithmeticError) as exc:
        print('REFUSED: '+str(exc),file=sys.stderr)
        sys.exit(2)
