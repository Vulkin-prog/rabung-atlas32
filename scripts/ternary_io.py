"""Read the fixed ternary million in bounded blocks, using only the stdlib.

This reader starts no worker and writes no data. See data/ternary-million/README.md.
"""
from bisect import bisect_left
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import struct
import zlib

DEFAULT_DATA = Path(__file__).resolve().parents[1] / 'data/ternary-million'
ROW = struct.Struct('<34H')
PARENT_ROW = struct.Struct('<HHIBII')
PARENT_HEADER = struct.Struct('<8sIIIII32s32s16s')
MAX_BLOCK_BYTES = 8 * 1024**2


class IntegrityError(ValueError):
    """A stored scientific commitment or format invariant failed."""


def require(condition, message):
    if not condition:
        raise IntegrityError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def decompress_exact(data, limit=MAX_BLOCK_BYTES):
    decoder = zlib.decompressobj()
    try:
        raw = decoder.decompress(data, limit + 1)
    except zlib.error as exc:
        raise IntegrityError('Corrupt zlib payload') from exc
    require(len(raw) <= limit and decoder.eof and not decoder.unused_data
            and not decoder.unconsumed_tail, 'Zlib size or stream boundary mismatch')
    return raw


def integer(data, cursor):
    value = 0
    for shift in range(0, 70, 7):
        require(cursor < len(data), 'Truncated position varint')
        byte = data[cursor]
        cursor += 1
        value |= (byte & 127) << shift
        if byte < 128:
            require((shift == 0 or byte != 0) and value < 2**64,
                    'Noncanonical or overflowing position varint')
            return value, cursor
    raise IntegrityError('Position varint overflow')


@dataclass(frozen=True)
class Profile:
    p: int
    s: int
    index: int
    fields: tuple
    runs: tuple
    parent: tuple

    @property
    def M(self):
        return self.fields[0]

    @property
    def L(self):
        return self.fields[1]

    @property
    def K(self):
        # r=3, p=1 mod 6 implies chi_3(-1)=1 and sigma=1.
        return max(self.M + 1, 2*self.L + 2)

    def kind(self, run):
        start, length, colour = run
        if start == 1:
            return 'initial'
        if start + length - 1 > (self.p - 1)//2:
            return 'central'
        return 'ordinary'

    def ordinary_count(self, threshold):
        require(18 <= threshold <= 27, 'Stored count thresholds are 18..27')
        return self.fields[threshold - 15]

    def canonical_count(self, threshold):
        return (self.ordinary_count(threshold) + int(self.L >= threshold)
                + int(self.fields[2] >= threshold))

    def as_dict(self):
        return {'p': self.p, 'r': 3, 'stratum': self.s, 'index': self.index,
                'M': self.M, 'L': self.L, 'K': self.K,
                'central_length_ge18': self.fields[2],
                'ordinary_counts': {str(k): self.ordinary_count(k) for k in range(18, 28)},
                'canonical_counts': {str(k): self.canonical_count(k) for k in range(18, 28)},
                'ordinary_colour_counts18': list(self.fields[13:16]),
                'first_window_echoes18': dict(zip(('2', '3'), self.fields[16:18])),
                'ordinary_position_bins18': list(self.fields[18:34]),
                'position_minimum_length': 21,
                'canonical_runs21': [
                    {'start': x, 'length': length, 'colour': colour,
                     'kind': self.kind((x, length, colour)),
                     'reflected_start': self.p - x - length + 1}
                    for x, length, colour in self.runs],
                'parent_features': {'cubic_signature': self.parent[2],
                                    'flags': self.parent[3],
                                    'crossing18_mask': self.parent[4],
                                    'crossing21_mask': self.parent[5]}}


def decode_records(fixed, positions, prime_data, parent, meta):
    """Validate all records before returning a block; no partially valid block."""
    n, s, start = meta['n'], meta['s'], meta['start']
    require(len(fixed) == n*ROW.size and len(prime_data) == 8*n
            and len(parent) == n*PARENT_ROW.size, 'Record array shape mismatch')
    profiles = []
    cursor = positions_count = all_runs18 = 0
    for i in range(n):
        fields = ROW.unpack_from(fixed, ROW.size*i)
        old = PARENT_ROW.unpack_from(parent, PARENT_ROW.size*i)
        p, = struct.unpack_from('<Q', prime_data, 8*i)
        M, L, central = fields[:3]
        require(18 <= M < p and 1 <= L <= M and central <= M, 'Invalid M, L or central length')
        require(central == 0 or (central >= 18 and central % 2 == 0), 'Invalid retained central length')
        require(fields[:2] == old[:2], 'M or L differs from the frozen parent profile')
        signature, flags, mask18, mask21 = old[2:]
        require(flags <= 7 and bool(flags & 4) == (s % 2 == 0), 'Parent cube class mismatch')
        require((L >= 2 if s % 2 == 0 else L == 1), 'Initial length contradicts cube class of 2')
        require(signature >> 22 == 0 and all((signature >> (2*j)) & 3 != 3 for j in range(11)),
                'Invalid normalized parent cubic signature')
        require(mask18 >> 30 == 0 and mask21 & ~mask18 == 0, 'Invalid parent crossing masks')
        require(bool(flags & 1) == (central < 18) and bool(flags & 2) == (central < 21),
                'Parent central gate differs from the retained catalogue')
        require(not (mask21 and M < 21), 'Parent rational witness exceeds the maximum')
        require(all(a >= b for a, b in zip(fields[3:12], fields[4:13])), 'Nonmonotone run counts')
        require(sum(fields[13:16]) == sum(fields[18:34]) == fields[3], 'Colour or position-bin totals disagree')
        count, cursor = integer(positions, cursor)
        require(count == fields[6] + int(L >= 21) + int(central >= 21), 'Canonical count at 21 differs')
        runs, ordinary_lengths = [], []
        previous = longest = 0
        for j in range(count):
            delta, cursor = integer(positions, cursor)
            excess, cursor = integer(positions, cursor)
            require(cursor < len(positions), 'Missing position colour')
            colour = positions[cursor]
            cursor += 1
            x, length = previous + delta, excess + 21
            require(delta > 0 and colour < 3 and 1 <= x <= (p-1)//2 and x + length < p,
                    'Invalid canonical run interval or colour')
            require(not runs or x >= runs[-1][0] + runs[-1][1], 'Overlapping maximal runs')
            previous = x
            longest = max(longest, length)
            if x == 1:
                require(length == L and colour == 0, 'Initial run differs from L')
            elif x + length - 1 > (p-1)//2:
                require(2*x + length == p + 1 and length == central, 'Central reflection identity failed')
            else:
                ordinary_lengths.append(length)
            runs.append((x, length, colour))
        for k in range(21, 28):
            require(fields[k-15] == sum(length >= k for length in ordinary_lengths),
                    'Stored positions and ordinary counts disagree')
        require((M < 21) == (central < 21 and fields[6] == 0), 'Maximum and central/ordinary gate disagree')
        require((M < 21 and longest == 0) or longest == M, 'Longest retained position differs from M')
        positions_count += count
        all_runs18 += fields[3] + int(L >= 18) + int(central >= 18)
        profiles.append(Profile(p, s, start+i, fields, tuple(runs), old))
    require(cursor == len(positions), 'Trailing position bytes')
    require(positions_count == meta['retained_canonical_runs21'], 'Block position total differs')
    require(all_runs18 == meta['all_retained_canonical_runs18'], 'Block canonical total at 18 differs')
    return profiles


class Dataset:
    def __init__(self, directory=DEFAULT_DATA):
        self.root = Path(directory).resolve()
        self.manifest = json.loads((self.root / 'manifest.json').read_text())
        require(self.manifest['schema'] == 1 and
                self.manifest['dataset'] == 'fixed-ternary-million-long-runs-v1', 'Unknown dataset schema')
        self.files = {v['path']: v for v in self.manifest['files']}
        require(len(self.files) == len(self.manifest['files']), 'Repeated manifest path')
        self.blocks = {(v['s'], v['start']): v for v in self.manifest['blocks']}
        require(len(self.blocks) == len(self.manifest['blocks']), 'Repeated block identity')
        self._lists = {}
        protocol_data = self.read('PROTOCOLE.json')
        require(sha(protocol_data) == self.manifest['protocol_sha256'], 'Frozen protocol identity differs')
        self.protocol = json.loads(protocol_data)
        require(self.protocol['n'] == self.manifest['profiles'] == 1000000, 'Population size differs')

    def read(self, name):
        require(name in self.files, 'Unlisted scientific file: ' + name)
        path = self.root / name
        require(path.resolve().is_relative_to(self.root) and not path.is_symlink(), 'Unsafe dataset path')
        data = path.read_bytes()
        item = self.files[name]
        require(len(data) == item['bytes'] and sha(data) == item['sha256'], 'File hash or size differs: ' + name)
        return data

    def prime_bytes(self, s):
        require(s in range(4), 'Stratum must be 0..3')
        if s not in self._lists:
            data = self.read(f'primes/stratum-{s}.u64')
            require(len(data) == 250000*8 and sha(data) == self.protocol['parent_list_hashes'][str(s)],
                    'Frozen prime list identity differs')
            numbers = tuple(v[0] for v in struct.iter_unpack('<Q', data))
            band = (2**33, 3*2**32) if s < 2 else (3*2**32, 2**34)
            require(all(p % 6 == 1 and band[0] <= p < band[1] for p in numbers), 'Prime list band or congruence differs')
            require(all(a < b for a, b in zip(numbers, numbers[1:])), 'Prime list is not sorted and distinct')
            self._lists[s] = (data, numbers)
        return self._lists[s][0]

    def read_block(self, entry):
        s, start, n = entry['s'], entry['start'], entry['n']
        require(s in range(4) and start % 4096 == 0 and n == min(4096, 250000-start)
                and 1 <= n <= 4096, 'Invalid bounded block shape')
        require(entry['name'] == f's{s}-{start:06d}', 'Block name differs from its identity')
        meta = json.loads(self.read(entry['metadata']))
        require(meta['status'] == 'PASS_COUNT_CATALOGUE_BLOCK' and
                (meta['s'], meta['start'], meta['n']) == (s, start, n), 'Block metadata identity differs')
        require(meta['protocol_sha256'] == self.manifest['protocol_sha256'], 'Block protocol differs')
        packed = self.read(entry['pack'])
        require(sha(packed) == meta['packed_sha256'] and len(packed) == meta['packed_bytes'], 'Count pack identity differs')
        prime_data = self.prime_bytes(s)[8*start:8*(start+n)]
        prime_hash = sha(b'vdw-compact-primes-u64le-v1' + prime_data)
        require(prime_hash == meta['prime_hash'], 'Block prime commitment differs')
        payload = decompress_exact(packed)
        require(len(payload) >= 4, 'Missing fixed-record size')
        size, = struct.unpack_from('<I', payload)
        require(size == n*ROW.size and len(payload) >= 4+size, 'Fixed-record size differs')
        fixed, positions = payload[4:4+size], payload[4+size:]
        require(sha(fixed) == meta['counts_sha256'] and sha(positions) == meta['positions21_sha256'],
                'Counts or positions commitment differs')
        commitment = (b'vdw-count-catalogue-extension-v1' + bytes.fromhex(meta['protocol_sha256'])
                      + prime_data + fixed + positions)
        require(sha(commitment) == meta['science_sha256'], 'Scientific block commitment differs')
        parent_pack = self.read(entry['parent'])
        require(sha(parent_pack) == meta['parent_frame_sha256'] and len(parent_pack) >= PARENT_HEADER.size,
                'Frozen parent frame identity differs')
        magic, ss, st, nn, zsize, rawsize, ph, sh, tag = PARENT_HEADER.unpack_from(parent_pack)
        require((magic, ss, st, nn) == (b'VDWCSP1\0', s, start, n) and
                len(parent_pack) == PARENT_HEADER.size + zsize and rawsize == n*PARENT_ROW.size,
                'Parent frame header differs')
        require(ph.hex() == prime_hash and sh.hex() == meta['parent_science_sha256'] and
                tag == bytes.fromhex(self.manifest['parent_frozen_sha256'])[:16], 'Parent frame commitments differ')
        parent = decompress_exact(parent_pack[PARENT_HEADER.size:], rawsize)
        require(len(parent) == rawsize and sha(b'vdw-sample-fields-v1' + prime_data + parent) == sh.hex(),
                'Parent exact-field commitment differs')
        return decode_records(fixed, positions, prime_data, parent, meta)

    def iter_profiles(self, stratum=None):
        require(stratum is None or stratum in range(4), 'Stratum must be 0..3')
        for entry in sorted(self.blocks.values(), key=lambda v: (v['s'], v['start'])):
            if stratum is None or entry['s'] == stratum:
                yield from self.read_block(entry)

    def by_index(self, s, index):
        require(0 <= index < 250000, 'Stratum index must be 0..249999')
        start = index//4096*4096
        require((s, start) in self.blocks, 'Unknown stratum or missing block')
        return self.read_block(self.blocks[s, start])[index-start]

    def by_prime(self, p):
        for s in range(4):
            self.prime_bytes(s)
            numbers = self._lists[s][1]
            i = bisect_left(numbers, p)
            if i < len(numbers) and numbers[i] == p:
                return self.by_index(s, i)
        raise KeyError(f'Prime {p} is not in this fixed probability sample')
