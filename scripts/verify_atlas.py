#!/usr/bin/env python3
"""Verify bounded examples or the complete published pack hashes and indices.

This checks stored data. It does not run a GPU census or re-evaluate every prime.
Use reproduce_analysis.py for complete decoding and descriptive re-analysis.
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import struct
import zlib
from atlas_codec import decode_frame, decode, science_hash
from atlas_io import ROOT, read_json, sha256, readonly, prime_hash


def check_pair(scan, a4, record):
    rows, sm = decode_frame(scan)
    values, am = decode(a4)
    identity = (record['lo'], record['hi'], record['primes'])
    assert (sm['lo'], sm['hi'], sm['count']) == identity
    assert (am['lo'], am['hi'], am['count']) == identity
    assert sm['witnesses'] and len(rows) == len(values)
    assert prime_hash(rows) == record['prime_hash'] == am['prime_hash']
    assert science_hash(rows) == record['science_hash'] == am['profile_science_hash']
    assert am['values_sha256'] == record['values_hash']
    for row, value in zip(rows, values):
        p = row['p']
        assert value == 0 if p == 2 else 0 <= value <= (p - 3) // 2
    return rows, values


def samples(directory):
    manifest = read_json(directory / 'manifest.json')
    count = 0
    for record in manifest['blocks']:
        content = []
        for item in record['frames']:
            path = directory / item['file']
            assert path.stat().st_size == item['bytes'] and sha256(path) == item['sha256']
            content.append(path.read_bytes())
        rows, unused = check_pair(*content, record)
        count += len(rows)
    return {'mode': 'non-contiguous samples', 'blocks': len(manifest['blocks']), 'primes': count}


def full_indices(data, hash_assets=True):
    manifest = read_json(ROOT / 'data/atlas/manifest.json')
    summary = read_json(ROOT / 'data/atlas/summary.json')
    expected_sizes = {v['file']: v['bytes'] for v in manifest['files']}
    if hash_assets:
        for number, item in enumerate(manifest['files'], 1):
            path = data / item['file']
            assert path.stat().st_size == item['bytes'], item['file']
            assert sha256(path) == item['sha256'], item['file']
            if number % 10 == 0:
                print(json.dumps({'checked_assets': number}), flush=True)
    scan = readonly(data / 'scan/ledger.sqlite')
    a4 = readonly(data / 'quadruplets/ledger.sqlite')
    try:
        assert scan.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        assert a4.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        smeta = dict(scan.execute('SELECT key,value FROM meta'))
        ameta = dict(a4.execute('SELECT key,value FROM meta'))
        sc = bytes(32); ac = bytes(32); lo = 1; primes = 0; nblocks = 0; total_a4 = 0
        offsets = {}; histograms = {name: collections.defaultdict(collections.Counter)
                                  for name in ['hist_max', 'hist_prefix', 'hist_K']}
        srows = scan.execute('SELECT id,lo,hi,n,pack,offset,size,prime_hash,science_hash,summary FROM blocks ORDER BY id')
        arows = iter(a4.execute('SELECT id,lo,hi,n,pack,offset,size,prime_hash,profile_hash,values_hash,summary FROM blocks ORDER BY id'))
        for s in srows:
            a = next(arows)
            assert s[:4] == a[:4] and s[0] == nblocks and s[1] == lo and s[2] > lo
            assert 0 < s[3] <= 4096 and s[7:9] == a[7:9]
            for group, row in [('scan', s), ('quadruplets', a)]:
                path = group + '/' + row[4]
                assert row[5] == offsets.get(path, 0) and row[6] > 0
                offsets[path] = row[5] + row[6]
                assert offsets[path] <= expected_sizes[path]
            bound = struct.pack('<QQ', s[1], s[2])
            sc = hashlib.sha256(sc + bytes.fromhex(s[7]) + bytes.fromhex(s[8]) + bound).digest()
            ac = hashlib.sha256(ac + bytes.fromhex(a[7]) + bytes.fromhex(a[8]) + bytes.fromhex(a[9]) + bound).digest()
            ss = json.loads(zlib.decompress(s[9])); aa = json.loads(zlib.decompress(a[10]))
            assert ss['primes'] == s[3]
            for name, orders in histograms.items():
                for r, values in ss[name].items():
                    orders[r].update({key: int(value) for key, value in values.items()})
            total_a4 += aa['sum_A4']; primes += s[3]; nblocks += 1; lo = s[2]
        assert next(arows, None) is None
        assert all(size == expected_sizes[path] for path, size in offsets.items())
        assert lo == 2**32 == int(smeta['next_lo']) == int(ameta['next_lo'])
        assert nblocks == manifest['blocks'] == int(ameta['next_block'])
        assert primes == manifest['primes'] == int(smeta['prime_count']) == int(ameta['prime_count'])
        assert total_a4 == manifest['sum_A4'] == int(ameta['sum_A4'])
        assert sc.hex() == smeta['chain'] == summary['profile_chain_sha256']
        assert ac.hex() == ameta['chain'] == summary['A4_chain_sha256']
        aggregate = json.loads(smeta['aggregate'])
        for name, orders in histograms.items():
            calculated = {r: dict(values) for r, values in orders.items()}
            assert calculated == aggregate[name] == summary['aggregate'][name]
        profiles = sum(sum(v.values()) for v in histograms['hist_max'].values())
        assert profiles == manifest['profiles']
        assert scan.execute('SELECT count(*) FROM audits').fetchone()[0] == summary['GPU_audits']
        assert a4.execute('SELECT count(*) FROM audits').fetchone()[0] == summary['A4_certificate_audits']
        # Decode the endpoints against both independent indexed commitments.
        from atlas_io import frame
        for block in [0, nblocks - 1]:
            s = scan.execute('SELECT lo,hi,n,prime_hash,science_hash FROM blocks WHERE id=?', (block,)).fetchone()
            vh = a4.execute('SELECT values_hash FROM blocks WHERE id=?', (block,)).fetchone()[0]
            check_pair(frame(data, 'scan', block), frame(data, 'quadruplets', block),
                       dict(zip(['lo','hi','primes','prime_hash','science_hash','values_hash'], (*s, vh))))
        return {'mode': 'complete indices', 'asset_hashes_checked': hash_assets,
                'blocks': nblocks, 'primes': primes, 'profiles': profiles,
                'sum_A4': total_a4, 'next_lo': lo, 'profile_chain': sc.hex(), 'A4_chain': ac.hex()}
    finally:
        scan.close(); a4.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, help='Root downloaded by download_data.py; defaults to bounded samples')
    parser.add_argument('--indices-only', action='store_true', help='Skip whole-pack SHA256 checks; verify indices and endpoint frames')
    args = parser.parse_args()
    if args.indices_only and not args.data:
        parser.error('--indices-only requires --data')
    result = full_indices(args.data, not args.indices_only) if args.data else samples(ROOT / 'data/sample')
    print(json.dumps({'status': 'PASS', **result}, indent=2))


if __name__ == '__main__':
    main()
