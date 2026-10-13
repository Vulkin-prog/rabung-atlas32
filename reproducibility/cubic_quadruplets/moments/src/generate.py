"""Generate nine bounded point-order certificates for the fixed CM curves.
No maximum, long-run scan, or statistical sample is generated here.
"""
from pathlib import Path
import sys, json
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'input'))
from elliptic_inherited import Curve, verify_trace
from moments import C_CURVES


def main():
    rows = []
    for p in (1000033, 4294967197, 4294967311):
        certs = []
        for j, (a, b) in enumerate(C_CURVES):
            c = Curve(p, a, b).certify_order(max_x=256)
            verify_trace(c)
            c['CM_index'] = j
            certs.append(c)
        rows.append({'p': p, 'CM_certificates': certs})
    obj = {'status': 'PASS', 'purpose': 'finite identity certificates, not a sample of the limiting law',
           'cases': rows}
    path = ROOT / 'results/CERTIFICATS_CM.json'
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'status': 'PASS', 'CM_point_order_certificates': 9}))

if __name__ == '__main__':
    main()
