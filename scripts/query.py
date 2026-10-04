#!/usr/bin/env python3
"""Print the stored profiles and exact quadratic A4 count at a prime."""
import argparse
import json
from pathlib import Path
from atlas_codec import decode_frame, decode
from atlas_io import ROOT, readonly, frame, read_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('p', type=int)
    parser.add_argument('--data', type=Path, help='Downloaded data root; without this, use the two sample blocks')
    args = parser.parse_args()
    if not 2 <= args.p < 2**32:
        parser.error('The domain is 2 <= p < 2^32')
    if args.data:
        connection = readonly(args.data / 'scan/ledger.sqlite')
        try:
            block = connection.execute('SELECT id FROM blocks WHERE lo<=? AND hi>?', (args.p, args.p)).fetchone()
        finally:
            connection.close()
        if block is None:
            parser.error('No indexed interval contains p')
        s, a = (frame(args.data, group, block[0]) for group in ['scan', 'quadruplets'])
    else:
        directory = ROOT / 'data/sample'
        block = next((v for v in read_json(directory / 'manifest.json')['blocks'] if v['lo'] <= args.p < v['hi']), None)
        if block is None:
            parser.error('p is outside the two sample intervals; download the full data')
        s, a = ((directory / item['file']).read_bytes() for item in block['frames'])
    rows, unused = decode_frame(s); values, unused = decode(a)
    for row, value in zip(rows, values):
        if row['p'] == args.p:
            print(json.dumps({**row, 'A4': value, 'A4_placeholder': args.p == 2}, indent=2)); return
    parser.error('p is not a stored prime')


if __name__ == '__main__':
    main()
