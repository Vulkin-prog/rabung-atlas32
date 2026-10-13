"""Offline editorial audit of imported data and printed tables. No CUDA calls."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import math
import statistics
import sys
import subprocess

if sys.flags.optimize:
    raise SystemExit('Run without -O: inherited finite checks use assertions.')

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / 'reproducibility/gpu_method'

def read(path):
    return json.loads(path.read_text(encoding='utf-8'))

def close(a, b):
    assert math.isclose(a, b, rel_tol=2e-12, abs_tol=2e-12), (a, b)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    # Literal integer colours for the small worked example in the article.
    colours = [pow(x, 50, 101) for x in range(1, 101)]
    runs = []
    start = 0
    for end in range(1, 101):
        if end == 100 or colours[end] != colours[start]:
            runs.append((end-start, start+1))
            start = end
    maximum = max(n for n, _ in runs)
    witness = min(s for n, s in runs if n == maximum)
    assert (maximum, runs[0][0], max(maximum+1, 2*runs[0][0]+2), witness) == (7, 1, 8, 19)
    assert 4*30*4-4*4**2 > 4*101 and 4*4 <= 2*30
    assert all((-25+4*j) % 101 == 4*(19+j) % 101 for j in range(5))
    assert -17*pow(4, -1, 101) % 101 == 21
    provenance = read(DATA / 'PROVENANCE.json')
    for name, record in provenance['files'].items():
        path = ROOT / name
        assert not path.is_symlink() and path.is_file(), name
        content = path.read_bytes()
        assert len(content) == record['bytes'], name
        assert hashlib.sha256(content).hexdigest() == record['sha256'], name

    identity = read(DATA / 'night/HISTORICAL_ENGINE_IDENTITY.json')
    historical_source = (DATA / 'source/engine32/build/original_scan_gpu.cu').read_bytes()
    assert hashlib.sha256(historical_source).hexdigest() == identity['historical_source_sha256']
    prefix = historical_source.split(b'int main(', 1)[0]
    assert hashlib.sha256(prefix).hexdigest() == identity['engine_prefix_sha256']
    assert identity['engine_prefix_identical'] is True

    text = (ROOT / 'paper/GPU_METHOD_EN.tex').read_text(encoding='utf-8')
    medians = read(DATA / 'night/MESURES_FINALES_32_BITS.json')
    cohort = read(DATA / 'night/FINAL_HOLDOUT.json')
    freeze = read(DATA / 'night/FINAL_FREEZE.json')
    assert freeze['new_holdout_not_yet_run'] and freeze['tuning_complete']
    assert cohort['selection_without_timings'] and cohort['used_for_final_validation_only']
    historical_table = []
    for name, size, selection in [('natural128', 128, 'natural32'),
                                   ('stress8', 8, 'stress32'),
                                   ('control15', 15, None)]:
        group = medians['cohorts'][name]
        rows = group['per_prime']
        assert len(rows) == size
        if selection:
            assert [row['p'] for row in rows] == cohort[selection]
        if name == 'stress8':
            assert all((row['p'] - 1) % 2520 == 0 for row in rows)
        for engine in group['sum_medians']:
            close(sum(row[engine] for row in rows), group['sum_medians'][engine])
        old = group['sum_medians']['historical_ref']
        new = group['sum_medians']['historical_best']
        close(old / new, group['speedups']['historical'])
        cells = f"&{size}&{old:.4f}&{new:.4f}&{old/new:.3f}"
        assert cells in text, cells
        historical_table.append({'cohort': name, 'primes': size,
                                 'historical_sum_hot_medians_ms': old,
                                 'selected_sum_hot_medians_ms': new,
                                 'factor': old/new})

    analysis = read(DATA / 'throughput/ANALYSE.json')
    releases = analysis['final_engine_holdout']
    selected = read(DATA / 'throughput/SELECTION_COMPLEMENTAIRE.json')
    assert selected['no_adaptation_to_holdout'] is True
    assert selected['status'] == 'FROZEN_BEFORE_LATE_HOLDOUT'
    counts = Counter()
    throughput_table = []
    windows_checked = 0
    for window in releases['windows']:
        raw = read(DATA / 'throughput' / window['evidence_file'])
        assert raw['status'] == 'PASS' and raw['stage'] == 'holdout'
        matches = [w for w in raw['windows'] if w['lo'] == window['lo']]
        assert len(matches) == 1
        detail = matches[0]
        profiles = detail['profiles_full']
        assert detail['primes'] == window['primes'] == len(profiles) == 1024
        actual = Counter()
        seen = set()
        for row in profiles:
            p = row['p']
            assert p not in seen and window['lo'] <= p < window['hi']
            seen.add(p)
            expected = {r for r in range(2, 10) if (p - 1) % r == 0}
            assert {x['r'] for x in row['profiles']} == expected
            for profile in row['profiles']:
                r = profile['r']; maximum = profile['max']; pre = profile['prefix']
                boundary = 2*pre+2 if ((p-1)//2) % r == 0 else pre+2
                assert profile['K'] == max(maximum+1, boundary)
                assert 0 < pre <= maximum <= (p-1)//r
                assert 0 < profile['witness'] <= (p-1)//2
                actual[r] += 1
        assert sum(actual.values()) == window['profiles'] == detail['profiles']
        assert dict(sorted(actual.items())) == {int(r): n for r, n in window['per_r'].items()}
        counts.update(actual)
        for variant, stats in window['times'].items():
            samples = [s for s in detail['samples'] if s['variant'] == variant]
            repeats = []
            for repeat in range(6):
                batch = sorted((s for s in samples if s['repeat'] == repeat), key=lambda s: s['block'])
                assert [s['block'] for s in batch] == list(range(8))
                assert all(s['primes'] == 128 for s in batch)
                repeats.append(sum(s['wall_seconds'] for s in batch))
            for computed, recorded in zip(repeats, stats['repeat_wall_seconds']):
                close(computed, recorded)
            close(statistics.mean(repeats), stats['mean_seconds'])
            close(statistics.stdev(repeats), stats['stdev_seconds'])
        old = window['times']['baseline']; new = window['times'][window['selected_variant']]
        factor = old['mean_seconds']/new['mean_seconds']
        close(factor, window['speedup'])
        assert all(a>b for a,b in zip(old['repeat_wall_seconds'], new['repeat_wall_seconds']))
        cells = (f"{window['lo']/1e9:.3f}&{window['bits']}&"
                 f"${old['mean_seconds']:.4f}\\pm{old['stdev_seconds']:.4f}$&"
                 f"${new['mean_seconds']:.4f}\\pm{new['stdev_seconds']:.4f}$&{factor:.3f}")
        assert cells in text, cells
        throughput_table.append({'lo': window['lo'], 'hi_exclusive': window['hi'],
                                 'bits': window['bits'], 'primes': 1024,
                                 'profiles': window['profiles'],
                                 'reference_seconds': old['mean_seconds'],
                                 'selected_seconds': new['mean_seconds'], 'factor': factor})
        windows_checked += 1
    assert windows_checked == 6 and releases['primes'] == 6144
    assert sum(counts.values()) == releases['profiles'] == 20512
    assert dict(sorted(counts.items())) == {int(r): n for r,n in releases['per_r'].items()}
    report = {
        'status': 'PASS', 'operation': 'OFFLINE_EDITORIAL_AUDIT',
        'imported_files_verified': len(provenance['files']),
        'historical_engine_prefix_identity_verified': True,
        'tables_checked': ['historical hot medians', 'later repeated batch wall times'],
        'historical_comparison': historical_table, 'later_throughput': throughput_table,
        'final_engine_primes': 6144, 'final_engine_profiles': 20512,
        'per_r': dict(sorted(counts.items())),
        'limitations': ['No CUDA compilation or execution.',
                        'No proof of completeness from witness or hash checks alone.',
                        'Historical per-prime medians are retained; individual hot samples are external.',
                        'Freeze flags are audited as archived declarations, not an independent chronology proof.'],
    }
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    subprocess.run([sys.executable, '-B', str(DATA / 'check_source_bridges.py')], check=True)
    print(json.dumps({k:v for k,v in report.items() if k not in ['historical_comparison','later_throughput','limitations']}))

if __name__ == '__main__':
    main()

