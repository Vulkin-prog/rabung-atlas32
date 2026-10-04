#!/usr/bin/env python3
"""Download pinned release assets with GitHub CLI, checking every SHA256.

The private repository requires gh authentication. No token is stored in this
repository or command output. The same command works after a public release.
Existing valid files are reused; an existing invalid file is preserved and the
download stops. Only completed, checksum-verified assets are installed.
"""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from atlas_io import ROOT, read_json, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dest', type=Path, default=ROOT / 'data/atlas')
    parser.add_argument('--indices-only', action='store_true', help='Download only the two derived SQLite indices')
    parser.add_argument('--list', action='store_true', help='List assets without downloading')
    args = parser.parse_args()
    manifest = read_json(ROOT / 'data/atlas/manifest.json')
    items = [v for v in manifest['files'] if not args.indices_only or v['file'].endswith('.sqlite')]
    if args.list:
        for v in items:
            print(v['asset'], v['bytes'], v['sha256'])
        return
    if shutil.which('gh') is None:
        parser.error('Install GitHub CLI and authenticate with gh auth login')
    args.dest.mkdir(parents=True, exist_ok=True)
    remaining = []
    for item in items:
        rel = Path(item['file'])
        if rel.is_absolute() or '..' in rel.parts or Path(item['asset']).name != item['asset']:
            raise ValueError('Invalid manifest path')
        path = args.dest / rel
        if path.exists():
            if path.stat().st_size != item['bytes'] or sha256(path) != item['sha256']:
                raise ValueError('Existing file has an unexpected checksum: ' + item['file'])
            print('Verified existing ' + item['file'], flush=True)
        else:
            remaining.append(item)
    # Temporary download and final destination share a filesystem; os.replace
    # adds no second copy of an asset. Retain at least 5 GiB free.
    required = sum(v['bytes'] for v in remaining) + 5 * 1024**3
    if shutil.disk_usage(args.dest).free < required:
        raise ValueError('Insufficient free space for assets plus a 5 GiB reserve')
    for number, item in enumerate(remaining, 1):
        with tempfile.TemporaryDirectory(prefix='.atlas-download-', dir=args.dest) as scratch:
            subprocess.run(['gh', 'release', 'download', manifest['release'], '--repo', manifest['repository'],
                            '--pattern', item['asset'], '--dir', scratch], check=True)
            source = Path(scratch) / item['asset']
            if source.stat().st_size != item['bytes'] or sha256(source) != item['sha256']:
                raise ValueError('Downloaded asset failed verification: ' + item['asset'])
            destination = args.dest / item['file']
            destination.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, destination)
        print(f'Installed {number}/{len(remaining)}: {item["file"]}', flush=True)


if __name__ == '__main__':
    main()
