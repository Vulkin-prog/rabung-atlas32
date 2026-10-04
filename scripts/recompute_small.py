#!/usr/bin/env python3
"""Recompute a bounded CPU subset and compare it with the committed first block.

This is an independent dense profile calculation and an exact A4 calculation.
It never launches a GPU job. The worker sources are unchanged campaign sources.
"""
import argparse
import json
from pathlib import Path
import struct
import subprocess
from atlas_codec import decode_frame, decode
from atlas_io import ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count', type=int, default=64, help='Number of initial primes; 1..4096')
    parser.add_argument('--build-dir', type=Path, default=ROOT / 'build')
    args = parser.parse_args()
    if not 1 <= args.count <= 4096:
        parser.error('The bounded reproduction count is 1..4096')
    rows, unused = decode_frame((ROOT / 'data/sample/scan-00000.vdw').read_bytes())
    golden_a4, unused = decode((ROOT / 'data/sample/quadruplets-00000.a4').read_bytes())
    rows = rows[:args.count]; primes = [v['p'] for v in rows]; generated = []
    requests = []
    for i in range(0, len(primes), 128):
        requests.append(str(i//128) + ' ' + ' '.join(map(str, primes[i:i+128])))
    worker = subprocess.run([str(args.build_dir / 'small_cpu')], input='\n'.join(requests) + '\nQUIT\n',
                            text=True, capture_output=True, check=True, timeout=30)
    for line in worker.stdout.splitlines():
        item = json.loads(line)
        if item['type'] == 'measurement':
            generated.append({'p': item['p'], 'profiles': item['profiles']})
    assert generated == rows, 'Dense profile calculation disagrees with the atlas'
    request = ('0 -1 ' + str(len(primes)) + ' ' + ' '.join(map(str, primes)) + '\nQUIT\n').encode()
    a4 = subprocess.run([str(args.build_dir / 'moment_worker')], input=request, capture_output=True,
                        check=True, timeout=30).stdout
    ready, payload = a4.split(b'\n', 1)
    assert json.loads(ready)['type'] == 'ready'
    header = struct.Struct('<8sQIIQQIIQQ')
    magic, identifier, n, certificate_size, nanos, total, attempts, small, first, last = header.unpack(payload[:header.size])
    assert (magic, identifier, n, certificate_size, first, last) == (b'A4BATCH1', 0, len(primes), 0, primes[0], primes[-1])
    assert len(payload) == header.size + 4*n
    values = list(struct.unpack('<' + 'I'*n, payload[header.size:]))
    assert values == golden_a4[:n] and sum(values) == total
    print(json.dumps({'status': 'PASS', 'primes': n, 'last_prime': last,
                      'profiles': sum(len(v['profiles']) for v in generated), 'sum_A4': total}))


if __name__ == '__main__':
    main()
