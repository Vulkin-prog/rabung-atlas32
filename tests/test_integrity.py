"""Bounded integrity and cross-frame agreement checks, including corruption."""
import sys
from pathlib import Path
import unittest
import zlib
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from atlas_codec import decode_frame, decode
from atlas_io import ROOT
from verify_atlas import samples


class Integrity(unittest.TestCase):
    def test_original_sample_pairs(self):
        result = samples(ROOT / 'data/sample')
        self.assertEqual(result['primes'], 6409)

    def test_damaged_frames_rejected(self):
        for name, decoder in [('scan-00000.vdw', decode_frame), ('quadruplets-00000.a4', decode)]:
            original = (ROOT / 'data/sample' / name).read_bytes()
            for damaged in [original[:-1], original + b'0', bytes([original[0] ^ 1]) + original[1:]]:
                with self.assertRaises((ValueError, zlib.error)):
                    decoder(damaged)


if __name__ == '__main__':
    unittest.main()
