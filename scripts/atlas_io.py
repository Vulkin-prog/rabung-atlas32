"""Shared read-only access to the published atlas. Python standard library only."""
from pathlib import Path
import hashlib
import json
import sqlite3
import struct

ROOT = Path(__file__).resolve().parents[1]


def read_json(path):
    return json.loads(Path(path).read_text())


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro&immutable=1', uri=True)


def prime_hash(rows):
    h = hashlib.sha256(b'vdw-primes-u64le')
    for row in rows:
        h.update(struct.pack('<Q', row['p']))
    return h.hexdigest()


def frame(data, group, block):
    """Read a block using its published index; never open the index for writing."""
    connection = readonly(Path(data) / group / 'ledger.sqlite')
    try:
        row = connection.execute('SELECT pack,offset,size FROM blocks WHERE id=?', (block,)).fetchone()
        if row is None:
            raise ValueError('Unknown block')
        pack, offset, size = row
        if Path(pack).name != pack or offset < 0 or not 0 < size < 128 * 1024**2:
            raise ValueError('Invalid indexed extent')
        with (Path(data) / group / pack).open('rb') as stream:
            stream.seek(offset)
            value = stream.read(size)
        if len(value) != size:
            raise ValueError('Truncated pack')
        return value
    finally:
        connection.close()
