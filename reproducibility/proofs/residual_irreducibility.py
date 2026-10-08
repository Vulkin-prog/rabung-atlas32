"""Bind residual character certificates to irreducible defining polynomials.

The stored primes certify irreducibility over Q and Q(i). At a prime
q = 1 (mod 4), either prime of Z[i] above q has residue field F_q.
Degree-preserving irreducibility there implies irreducibility over Q(i).
All polynomial coefficient lists in this file's JSON interface descend.
"""
import argparse
import json
from pathlib import Path

from exact_fields import irreducible_mod, modpoly, pgcd, ppow, prime


def finite_witness(polynomial, q):
    if (not isinstance(q, int) or isinstance(q, bool) or not prime(q)
            or not polynomial or polynomial[0] % q == 0):
        raise ValueError("The reduction must preserve degree at a prime.")
    ascending = list(reversed(polynomial))
    degree = len(ascending) - 1
    if degree < 1 or not irreducible_mod(ascending, q):
        raise ValueError("The polynomial reduction is reducible.")
    x = modpoly([0, 1], ascending, q)
    powers = [x]
    for _ in range(degree):
        powers.append(ppow(powers[-1], q, ascending, q))
    gcd_degrees = {}
    for ell in range(2, degree + 1):
        if degree % ell == 0 and prime(ell):
            j = degree // ell
            difference = [0] * max(len(powers[j]), len(x))
            for i, coefficient in enumerate(powers[j]):
                difference[i] += coefficient
            for i, coefficient in enumerate(x):
                difference[i] -= coefficient
            gcd_degrees[str(j)] = len(pgcd(ascending, difference, q)) - 1
    return {"prime": q, "frobenius_powers_ascending": powers,
            "proper_subfield_gcd_degrees": gcd_degrees}


def find_witness(polynomial, gaussian=False):
    for q in range(2, 5001):
        if prime(q) and (not gaussian or q % 4 == 1):
            try:
                return finite_witness(polynomial, q)
            except ValueError:
                pass
    raise ValueError("No witness found below 5001.")


def qualified_models(root):
    for degree in range(2, 7):
        data = json.loads((root / f"RESIDUAL_QUALIFICATION_D{degree}.json").read_text())
        for model in data["models"]:
            yield degree, model


def load_witnesses(root):
    data = json.loads((root / "RESIDUAL_IRREDUCIBILITY.json").read_text())
    if data["schema"] != 1 or data["coefficient_order"] != "descending":
        raise ValueError("Unknown irreducibility certificate format.")
    records = data["models"]
    indexed = {entry["certificate"]: entry for entry in records}
    expected = {model["certificate"] for _, model in qualified_models(root)}
    if len(indexed) != len(records) or set(indexed) != expected or len(records) != 100:
        raise ValueError("Irreducibility records must cover exactly the 100 models.")
    return indexed


def verify_binding(model, character_certificate, witness):
    polynomial = model["polynomial"]
    if (character_certificate["polynomial"] != polynomial
            or witness["polynomial"] != polynomial
            or witness["certificate"] != model["certificate"]
            or witness["degree"] != len(polynomial) - 1):
        raise ValueError("Qualification, character and irreducibility polynomials differ.")
    for key in ("over_Q", "over_Q_i"):
        stored = witness[key]
        q = stored["prime"]
        if key == "over_Q_i" and q % 4 != 1:
            raise ValueError("The Gaussian-field witness prime must split in Q(i).")
        if stored != finite_witness(polynomial, q):
            raise ValueError("An irreducibility witness has been altered.")


def verify(root):
    witnesses = load_witnesses(root)
    maxima = {"over_Q": {}, "over_Q_i": {}}
    for degree, model in qualified_models(root):
        character = json.loads((root / model["certificate"]).read_text())
        witness = witnesses[model["certificate"]]
        verify_binding(model, character, witness)
        for key in maxima:
            maxima[key][degree] = max(maxima[key].get(degree, 0), witness[key]["prime"])
    return {"status": "PASS", "models": len(witnesses),
            "finite_irreducibility_witnesses": 2 * len(witnesses),
            "largest_witness_prime_by_degree": maxima}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path,
                        default=Path(__file__).resolve().parent.parent / "evidence/new")
    parser.add_argument("--generate", action="store_true",
                        help="Regenerate the small witness file, then replay it.")
    args = parser.parse_args()
    if args.generate:
        records = [{"certificate": model["certificate"], "degree": degree,
                    "polynomial": model["polynomial"],
                    "over_Q": find_witness(model["polynomial"]),
                    "over_Q_i": find_witness(model["polynomial"], gaussian=True)}
                   for degree, model in qualified_models(args.root)]
        target = args.root / "RESIDUAL_IRREDUCIBILITY.json"
        target.write_text(json.dumps({"schema": 1, "coefficient_order": "descending",
                                      "models": records}, indent=2) + "\n")
    print(json.dumps(verify(args.root), sort_keys=True))


if __name__ == "__main__":
    main()
