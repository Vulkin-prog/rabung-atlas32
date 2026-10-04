"""Real catalogue fixtures and malformed-input rejection, with no GPU."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
import zlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from ternary_io import Dataset, DEFAULT_DATA, IntegrityError, decompress_exact, integer
from verify_ternary import exact_run_check


class TernaryIntegrity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = Dataset()

    def test_real_graph_witness_and_last_block(self):
        profile = self.data.by_prime(8611114591)
        self.assertEqual((profile.s, profile.index, profile.M, profile.L, profile.K), (0, 1241, 23, 2, 24))
        self.assertEqual(profile.runs, ((2870371521, 23, 2), (3348834145, 21, 2)))
        self.assertEqual(exact_run_check(profile), 2)
        for s in range(4):
            self.assertEqual(self.data.by_index(s, 249999).index, 249999)

    def test_damaged_scientific_pack_is_refused(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            entry = self.data.blocks[0, 0]
            for rel in ('manifest.json', 'PROTOCOLE.json', entry['metadata'], entry['pack']):
                target = directory / rel
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(DEFAULT_DATA / rel, target)
            path = directory / entry['pack']
            damaged = bytearray(path.read_bytes())
            damaged[-1] ^= 1
            path.write_bytes(damaged)
            with self.assertRaisesRegex(IntegrityError, 'File hash or size differs'):
                Dataset(directory).read_block(entry)

    def test_noncanonical_and_truncated_varints_are_refused(self):
        for raw in (b'\x80', b'\x80\x00', b'\xff'*9+b'\x02', b'\x80'*10):
            with self.subTest(raw=raw), self.assertRaises(IntegrityError):
                integer(raw, 0)

    def test_compressed_trailing_truncated_or_excessive_data_is_refused(self):
        packed = zlib.compress(b'x'*128)
        for raw, limit in ((packed + b'extra', 128), (packed[:-1], 128), (packed, 127)):
            with self.subTest(raw=raw[-10:]), self.assertRaises(IntegrityError):
                decompress_exact(raw, limit)

    def test_command_reports_missing_prime_and_unretained_threshold(self):
        command = [sys.executable, '-B', str(DEFAULT_DATA.parents[1] / 'scripts/query_ternary.py')]
        missing = subprocess.run(command + ['7'], text=True, capture_output=True)
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn('not in this fixed probability sample', missing.stderr)
        shorter = subprocess.run(command + ['--export', 'positions', '--minimum', '18'],
                                 text=True, capture_output=True)
        self.assertNotEqual(shorter.returncode, 0)
        self.assertIn('not retained', shorter.stderr)
        queried = subprocess.run(command + ['--stratum', '0', '--index', '1241'],
                                 text=True, capture_output=True, check=True)
        self.assertEqual(json.loads(queried.stdout)['p'], 8611114591)


if __name__ == '__main__':
    unittest.main()
