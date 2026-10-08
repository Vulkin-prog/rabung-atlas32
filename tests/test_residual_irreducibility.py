"""Reject reducible fields, mismatched polynomials and altered witnesses."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reproducibility/proofs"))
import residual_irreducibility as ri


class ResidualIrreducibility(unittest.TestCase):
    def test_complete_model_inventory(self):
        root = Path(__file__).resolve().parents[1] / "reproducibility/evidence/new"
        report = ri.verify(root)
        self.assertEqual(report["models"], 100)
        self.assertEqual(report["finite_irreducibility_witnesses"], 200)

    def test_reducible_and_degree_dropping_reductions_are_rejected(self):
        for polynomial, q in [([1, 0, -1], 5), ([5, 1, 1], 5),
                              ([1, 0, 1], 9), ([1, 0, 1], True)]:
            with self.assertRaises(ValueError):
                ri.finite_witness(polynomial, q)

    def test_rational_irreducibility_does_not_suffice_over_gaussian_field(self):
        polynomial = [1, 0, 1]
        model = {"certificate": "example.json", "polynomial": polynomial}
        witness = {**model, "degree": 2, "over_Q": ri.finite_witness(polynomial, 3),
                   "over_Q_i": ri.finite_witness(polynomial, 3)}
        with self.assertRaises(ValueError):
            ri.verify_binding(model, {"polynomial": polynomial}, witness)
        with self.assertRaises(ValueError):
            ri.finite_witness(polynomial, 5)

    def test_polynomial_and_frobenius_tampering_are_rejected(self):
        polynomial = [1, 0, -2]
        model = {"certificate": "example.json", "polynomial": polynomial}
        witness = {**model, "degree": 2, "over_Q": ri.find_witness(polynomial),
                   "over_Q_i": ri.find_witness(polynomial, gaussian=True)}
        ri.verify_binding(model, {"polynomial": polynomial}, witness)
        with self.assertRaises(ValueError):
            ri.verify_binding(model, {"polynomial": [1, 1, -2]}, witness)
        broken = copy.deepcopy(witness)
        broken["over_Q_i"]["frobenius_powers_ascending"][-1][0] += 1
        with self.assertRaises(ValueError):
            ri.verify_binding(model, {"polynomial": polynomial}, broken)


if __name__ == "__main__":
    unittest.main()
