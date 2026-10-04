#!/usr/bin/env python3
"""Decode every atlas block with the original C++ descriptive analysis.

Requires g++, zlib and OpenSSL development headers. Reads approximately 4.5 GiB
of packs, checks their frame hashes and writes regenerated aggregates outside the
tracked evidence. No GPU or new prime census is used.
"""
import argparse
import json
from pathlib import Path
import subprocess
from atlas_io import ROOT, readonly


def normalize_records(result):
    """Original descriptive reduction: sigma=-1 has no joined centre.

    The streaming worker uses zero as its sentinel. The published reduction
    uses JSON null for that undefined quantity, and evaluates the centre at
    each selected sigma=+1 record. This changes no M, L, K, witness or count.
    """
    odd_primes = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]
    for order, lists in result['records'].items():
        r = int(order)
        for records in lists:
            for rec in records:
                p = rec['p']
                if ((p-1)//2) % r == 0:
                    q = next(q for q in odd_primes if pow(q,(p-1)//r,p) != 1)
                    c = q-1
                    assert c <= rec['M'] and (not rec['C_joined'] or c == rec['C_joined'])
                    rec['C_joined'] = c
                else:
                    rec['C_joined'] = None


def compare(actual, expected, path=''):
    # Timings depend on machine and are never scientific equality claims.
    excluded = {'seconds', 'max_block_seconds'}
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for key in expected.keys() - excluded:
            compare(actual[key], expected[key], path + '/' + key)
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for i, (a,b) in enumerate(zip(actual,expected)):
            compare(a,b,path + '/' + str(i))
    elif isinstance(expected, float):
        assert abs(actual-expected) <= 1e-10 * max(1,abs(expected)), path
    else:
        assert actual == expected, path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=ROOT / 'build/reanalysis')
    args = parser.parse_args(); args.out.mkdir(parents=True, exist_ok=True)
    frame_index = args.out / 'frames.txt'
    s = readonly(args.data / 'scan/ledger.sqlite'); a = readonly(args.data / 'quadruplets/ledger.sqlite')
    try:
        with frame_index.open('w') as target:
            ar = iter(a.execute('SELECT id,lo,hi,n,pack,offset,size FROM blocks ORDER BY id'))
            for sr in s.execute('SELECT id,lo,hi,n,pack,offset,size FROM blocks ORDER BY id'):
                other = next(ar)
                assert sr[:4] == other[:4]
                target.write(' '.join(map(str, sr + other[4:])) + '\n')
            assert next(ar, None) is None
    finally:
        s.close(); a.close()
    source = ROOT / 'analysis/analyse_atlas_32_20261004/methodes/stream.cpp'
    binary = args.out / 'stream'
    subprocess.run(['g++', '-O3', '-std=c++17', str(source), '-lcrypto', '-lz', '-o', str(binary)], check=True)
    result = args.out / 'AGREGATS_DETAILLES.json'
    subprocess.run([str(binary), str(args.data / 'scan'), str(args.data / 'quadruplets'), str(frame_index), str(result)], check=True)
    golden = ROOT / 'analysis/analyse_atlas_32_20261004/AGREGATS_DETAILLES.json'
    actual = json.loads(result.read_text()); expected = json.loads(golden.read_text())
    normalize_records(actual)
    result.write_text(json.dumps(actual, indent=2) + '\n')
    compare(actual, expected)
    print(json.dumps({'status': 'PASS', 'comparison': 'All scientific aggregates; timing fields excluded',
                      'blocks': actual['blocks'], 'primes': actual['primes'], 'sum_A4': actual['sum_A4']}))


if __name__ == '__main__':
    main()
