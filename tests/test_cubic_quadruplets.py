"""Certify the new six-curve oracle, including wrong-curve certificates."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import cubic_quadruplets as cq


class CubicQuadruplets(unittest.TestCase):
    def test_archived_exact_counts(self):
        for p, expected in [(1000033, 36876), (4294967197, 159072704),
                            (4294967311, 159072584)]:
            ec, cc = cq.stored_certificates(p)
            self.assertEqual(cq.evaluate(p, ec, cc)["A4_cubic"], expected)

    def test_bad_moduli_are_refused(self):
        for p in [True, -7, 3, 5, 25, 10**10 + 9]:
            with self.assertRaises(ValueError):
                cq.admissible(p)

    def test_missing_certificate_set_is_refused(self):
        with self.assertRaises(ValueError):
            cq.stored_certificates(31)

    def test_wrong_curve_and_altered_order_are_refused(self):
        ec, cc = cq.stored_certificates(1000033)
        for kind in ["curve", "trace", "order", "swapped"]:
            broken = copy.deepcopy(cc)
            if kind == "curve":
                broken[0]["b"] += 1
            elif kind == "trace":
                broken[0]["trace"] += 1
            elif kind == "order":
                broken[0]["point_order"] += 1
            else:
                broken[0], broken[1] = broken[1], broken[0]
            with self.assertRaises((ValueError, ArithmeticError)):
                cq.evaluate(1000033, ec, broken)


if __name__ == "__main__":
    unittest.main()
