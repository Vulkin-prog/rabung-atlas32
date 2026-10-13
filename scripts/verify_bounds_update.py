#!/usr/bin/env python3
"""Check the new pointwise records and regenerate their separate finite closure.

Default checks inspect retained full CPU executions and literal local powers;
they do not rescan billion-position words. --generate updates derived files.
The published v1.0 inputs and their claim gate are kept byte-identical.
"""
import argparse
import copy
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]
PACK = ROOT / 'reproducibility/bounds_update'
LEGACY = PACK / 'published_v1_0'
EXPECTED = [(2072886931, 3, 17), (4291179427, 3, 18), (4294961809, 3, 19),
            (4294967197, 3, 20), (4195501393, 2, 27), (4294720151, 2, 28),
            (5899373413, 3, 18), (14413721599, 3, 19), (17177303971, 3, 20)]

def require(condition, reason):
    if not condition:
        raise ArithmeticError(reason)

def read(path):
    return json.loads(path.read_text())

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prime(p):
    return p >= 2 and all(p % d for d in range(2, math.isqrt(p) + 1))

def factors(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out

def records():
    meta = read(PACK / 'PUBLISHED_IMPORTS.json')
    for item in meta['files']:
        data = (LEGACY / item['path']).read_bytes()
        require(len(data) == item['bytes'] and hashlib.sha256(data).hexdigest() == item['sha256'], item['path'])
        require(hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest() == item['git_blob'], 'published Git identity')
    profiles = {}
    for p, r, k in EXPECTED:
        filename = f'p{p}_r{r}_k{k}.json'
        item = read(PACK / 'results' / filename)
        row = item['profile']
        require(json.loads(item['stdout']) == row, 'raw output binding')
        require(item['exit_code'] == 0 and not item['stderr'] and row['accepted'] and row['Bstar'], filename)
        require((row['p'], row['r'], row['k']) == (p, r, k), 'triple binding')
        checker = 'independent_profile.cpp' if p < 2**32 else 'two_bit_profile.cpp'
        require(item['source_sha256'] == sha(PACK / checker), 'checker source binding')
        require(item['command'][1:] == list(map(str, (p, r, k))), 'command binding')
        require(prime(p) and (p - 1) % r == 0 and p > k, 'prime preconditions')
        g = row['g']
        require(pow(g, p-1, p) == 1 and all(pow(g, (p-1)//q, p) != 1 for q in factors(p-1)), 'primitive root')
        require(row['full_positions'] == p-1 and row['powering_checks'] == row['reflection_checks'] == 2048, 'full execution scope')
        require(row['M'] < k and row['K'] <= k and row['strict_bound'] == (k-1)*p+1, 'criterion or strict integer')
        boundary = 2*row['L']+2 if ((p-1)//2) % r == 0 else row['L']+2
        require(row['K'] == max(row['M']+1, boundary), 'threshold')
        colour = lambda x: pow(x % p, (p-1)//r, p)
        require(all(colour(x) == 1 for x in range(1, row['L']+1)) and colour(row['L']+1) != 1, 'prefix')
        w, M = row['w'], row['M']
        require(1 <= w <= p-M and all(colour(x) == colour(w) for x in range(w, w+M)), 'run witness')
        require((w == 1 or colour(w-1) != colour(w)) and (w+M == p or colour(w+M) != colour(w)), 'run neighbours')
        profiles[(p, r, k)] = {'file': 'results/' + filename, **row}
    # The six atlas rechecks must bind the finite-domain catalogue, not a global optimum.
    catalogue = {(int(x['r']), int(x['k'])): x for x in csv.DictReader((ROOT / 'data/derived/ATLAS_CATALOGUE_K3_K64.csv').open())}
    for p, r, k in EXPECTED[:6]:
        row = catalogue[(r, k)]
        require(int(row['AB_max_p']) == p and int(row['construction_length']) == (k-1)*p+1, 'atlas binding')
    return profiles

def module():
    path = LEGACY / 'tools/close_baselines.py'
    spec = importlib.util.spec_from_file_location('published_closure', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def derive(profiles):
    mod = module()
    old, _, _, _ = mod.closure()
    require(old == read(LEGACY / 'audit/generated/bounds_closure.json'), 'published closure changed')
    current = []
    for p, r, k0 in [(2072886931, 3, 17), (5899373413, 3, 18), (14413721599, 3, 19),
                     (17177303971, 3, 20), (4195501393, 2, 27), (4294720151, 2, 28)]:
        for k in range(k0, 29) if k0 == 20 else [k0]:
            current.append({'id': f'update_w{r}_k{k}_p{p}', 'colors': r, 'length': k,
                            'prime': p, 'lower_bound': (k-1)*p+1,
                            'origin': 'bounds_update_2026_10_09', 'check_record': profiles[(p,r,k0)]['file'],
                            'qualified_at_length': k0,
                            'derivation': 'pointwise' if k == k0 else 'fixed-prime Rabung monotonicity'})
    claims = read(LEGACY / 'audit/claims.json')
    require(len(claims['claims']) == 13, 'published gate cardinality')
    merged = copy.deepcopy(claims)
    merged['claims'] += current
    baselines = read(LEGACY / 'audit/baselines.json')
    baselines['scope']['max_colors'] = 6
    baselines['sources']['bounds_update_2026_10_09'] = {
        'title': 'Current pointwise CPU records, with fixed-prime monotonicity',
        'location': 'reproducibility/bounds_update/results/', 'priority': 'No new literature-wide priority audit.'}
    for claim in current:
        baselines['ordinary_seeds'].append({'claim_id': claim['id']})
        baselines['ring_seeds'].append({'claim_id': claim['id']})
    with tempfile.TemporaryDirectory(prefix='rabung-bounds-closure-') as tmp:
        mod.ROOT = Path(tmp)
        mod.CLAIMS_PATH = mod.ROOT / 'audit/claims.json'
        mod.BASELINES_PATH = mod.ROOT / 'audit/baselines.json'
        mod.CLAIMS_PATH.parent.mkdir()
        mod.CLAIMS_PATH.write_text(json.dumps(merged, sort_keys=True) + '\n')
        mod.BASELINES_PATH.write_text(json.dumps(baselines, sort_keys=True) + '\n')
        result, nodes, winners, _ = mod.closure()
    result['inputs']['claims'] = 'merged published and current registries (canonical JSON)'
    result['inputs']['baselines'] = 'published registry plus current pointwise seeds; max_colors=6 (canonical JSON)'
    result['scope_note'] = 'Finite admitted corpus through length 28; four-colour legacy scope extended to six colours. No fresh literature priority audit.'
    table = ['% Generated by scripts/verify_bounds_update.py --generate.',
             r'\begin{table}[htbp]\centering\small\setlength{\tabcolsep}{4pt}',
             r'\begin{tabular}{rrrrl}\toprule $k$ & direct $p$ & direct strict bound & closed strict bound & method\\\midrule']
    fmt = lambda n: f'{n:,}'.replace(',', r'\,')
    for k in range(17, 29):
        claim = next(x for x in current if x['colors'] == 3 and x['length'] == k)
        win = nodes[winners[('ordinary',3,k)]]
        method = {'direct_rabung': 'direct Rabung', 'bct2018': 'BCT', 'published_seed': 'Berlekamp', 'length_monotonicity': 'Berlekamp + mono.'}[win['method']]
        table.append(f"{k} & ${fmt(claim['prime'])}$ & ${fmt(claim['lower_bound'])}$ & ${fmt(win['lower_bound'])}$ & {method}" + r'\\')
    table += [r'\bottomrule\end{tabular}',
              r'\caption{Updated three-colour bounds, all in the strict convention $W(3,k)>B$. The direct rows use four primitive ternary certificates. Closed values include the two new binary inputs and the admitted published seeds. The unresolved Landman--Robertson statement remains a separate comparator; no unrestricted optimality is asserted.}\label{tab:updated-closed}',
              r'\end{table}', '']
    prim = ['% Generated from retained full CPU executions.',
            r'\begin{table}[htbp]\centering\small',
            r'\begin{tabular}{rrrrrr}\toprule $r$ & minimum $k$ & $p$ & $M$ & $L$ & $K$\\\midrule']
    for claim in current:
        if claim['qualified_at_length'] != claim['length']:
            continue
        row = profiles[(claim['prime'],claim['colors'],claim['length'])]
        prim.append(f"{claim['colors']} & {claim['length']} & ${fmt(claim['prime'])}$ & {row['M']} & {row['L']} & {row['K']}" + r'\\')
    prim += [r'\bottomrule\end{tabular}',r'\caption{Six primitive certificates of the addendum, rechecked by full ordinary CPU scans. Four ternary and two binary triples give fourteen direct inequalities in the reported window.}\label{tab:updated-primitive}',r'\end{table}', '']
    summary = {'schema':1,'new_full_CPU_records':len(profiles),'new_primitive_certificates':6,
               'new_direct_rows':len(current),'claims':current,'published_claims_unchanged':13,
               'published_doi':'10.5281/zenodo.22980965','source_hashes': {n:sha(PACK/n) for n in ['independent_profile.cpp','two_bit_profile.cpp']},
               'limits':['Retained full author-controlled CPU executions; no external replication.',
                         'The default verification does not rescan the full periods.',
                         'No fresh global literature search or full GPU campaign.']}
    return {'CLAIMS_UPDATE.json':json.dumps(summary,indent=2,sort_keys=True)+'\n',
            'CLOSURE_UPDATE.json':json.dumps(result,indent=2,sort_keys=True)+'\n',
            'TABLE_UPDATE.tex':'\n'.join(table), 'PRIMITIVES_UPDATE.tex':'\n'.join(prim)}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--generate', action='store_true')
    args = parser.parse_args()
    files = derive(records())
    for name, content in files.items():
        path = PACK / name
        if args.generate:
            path.write_text(content)
        else:
            require(path.is_file() and path.read_text() == content, 'stale generated ' + name)
    manuscript = (ROOT / 'paper/RABUNG_BOUNDS_EN.tex').read_text()
    for name in ['TABLE_UPDATE.tex', 'PRIMITIVES_UPDATE.tex']:
        require(files[name].strip() in manuscript, 'manuscript table differs: ' + name)
    print(json.dumps({'status':'PASS','full_CPU_records':9,'new_primitive_certificates':6,
                      'new_direct_rows':14,'published_claims_unchanged':13,
                      'full_periods_rescanned_by_this_command':False}))

if __name__ == '__main__':
    main()
