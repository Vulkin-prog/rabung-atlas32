# Elliptic integration review, 8 October 2026

This is an internal proof review and a record of finite reproductions for the
companion manuscript. It is not an independent external review. The manuscript
gives an unconditional argument for four moments, using the cited analytic
theorems; independent checking of that argument remains an editorial task
before publication, rather than an additional mathematical hypothesis.

## Proof obligations and changes during integration

| Obligation | Check or clarification in the companion |
| --- | --- |
| Four-position phase | The projective substitution keeps the factor conjugate to the cubic character of (c−d)/(a−b). |
| Genus-two reduction | Both points at infinity and the two elliptic quotients are included; the degree-three isogeny is given explicitly. |
| Consecutive quadruplets | All 27 character terms and the four forbidden zero positions are accounted for. A self-contained three-position Jacobi-sum derivation is included. |
| Related CM twists | Fermat/Jacobi identities give the exact three-twist relation, with separate cube and noncube cases. They are not replaced by independent random variables. |
| Non-CM qualifications | Two ordinary reductions with incompatible endomorphism fields are given for each curve. |
| Symmetric powers and base change | Newton–Thorne supplies the needed cuspidal symmetric powers. Nonzero symmetric-power traces at inert primes exclude every required quadratic self-twist. |
| GL(2) tensor cuspidality | The three curves are nondihedral after base change; their absolute traces at 31 exclude character-twist equivalence. |
| Rankin–Selberg main terms | Every degree-at-most-four pattern is routed to two cuspidal factors, at least one self-contragredient. Distinct absolute traces at 31 also exclude the equal-rank collisions after all required unitary and imaginary norm twists. This strengthens the original routing explanation. |
| Prime powers and bad places | Tempered good elliptic parameters bound higher powers. The fixed bad Euler parameters are handled explicitly before partial summation. |
| Prime ideals versus rational primes | Residue-degree-two places are negligible; each split rational prime has two conjugate degree-one places. Cubic class indicators give weights 1/3 and 2/3. |

The review found no unresolved obstruction to this finite-degree argument.
That conclusion does not replace independent checking of the analytic passage.
The stated conclusions are four moments only. They establish neither a full
limiting distribution, nor centring of the raw count with an o(1) error, nor
critical intensity or a void probability for logarithmic-length maxima.

## Primary analytic references checked

- Newton–Thorne, [Symmetric power functoriality for holomorphic modular forms](https://arxiv.org/abs/2009.07180), Theorem A.
- Ramakrishnan, [Modularity of the Rankin–Selberg L-series, and multiplicity one for SL(2)](https://arxiv.org/abs/math/0007203), Theorem M and the cyclic base-change criterion on printed page 54.
- Gillespie–Ji, [Prime number theorems for Rankin–Selberg L-functions over number fields](https://arxiv.org/abs/0910.3660), Theorem 1.1.
- Chen–Shen, [Effective Sato–Tate distributions for surfaces arising from products of elliptic curves](https://arxiv.org/abs/2309.08848), Section 3.3 and Theorem 3.6, for the CM Hecke character.
- Sutherland, [Lecture 13](https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf) and [Lecture 21](https://math.mit.edu/classes/18.783/2023/LectureNotes21.pdf), for ordinary endomorphism fields and modularity background.

The manuscript states the exact hypotheses used, including fixed objects,
unitary normalization, cuspidality, self-contragredience and pole exclusions.
Its identity proofs and finite certificates do not stand in for these theorems.
The targeted literature check is not a complete bibliographic priority review.

## Finite verification

The retained elliptic package checks 2,625 large trace certificates, the small
direct-count identity cases, and the exact examples in the manuscript. The
moment package checks 207 small primes, 1,242 direct elliptic point counts,
621 CM/Jacobi comparisons, 621 Fermat comparisons, nine large CM trace
certificates, and 315 analytic routing cases for the 70 monomials through
degree four. Routing enumeration checks coverage, not automorphy or a prime
number theorem.

The public count oracle checks six curves. The third non-CM certificate is
transported from the archived fractional-parameter model to the integer model
by (x,y) → (9x,27y), and is then verified on that claimed integer curve. It
does not infer the count from an unchecked fast Jacobi-sum output.

The remaining publication priority is an independent, line-by-line review of
the analytic section and its interfaces with the cited theorems, followed by
an independent review of the algebraic reduction. The PR is a research draft;
merging it is not journal submission or mathematical peer review.
