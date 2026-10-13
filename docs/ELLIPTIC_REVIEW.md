# Elliptic integration review, 8 October 2026

This records internal proof checking and finite reproductions, now supplemented
by the separate priority audit and AI-assisted counter-review supplied on
8 October 2026. The [response](REVIEW_RESPONSE_2026-10-08_FR.md) identifies
those reports and the resulting changes. They do not constitute specialist
human peer review. The manuscript gives an unconditional argument for four
moments using the stated analytic theorems; specialist checking remains an
editorial task before submission, rather than a mathematical hypothesis.

## Proof obligations and changes during integration

| Obligation | Check or clarification in the companion |
| --- | --- |
| Four-position phase and attribution | The identity is an explicit Greene–Lennon corollary. A fractional linear substitution keeps the cubic phase and the projective contribution one; the hypergeometric normalization and trace theorem are stated. |
| Genus-two reduction and attribution | Both points at infinity and elliptic quotients are included. Freeman–Satoh Proposition 4.2 and Guillevic–Vergnaud Section 3.1 are cited for the same family and isogeny. |
| Consecutive quadruplets | All 27 character terms and the four forbidden zero positions are accounted for. A self-contained three-position Jacobi-sum derivation is included. |
| Related CM twists | Fermat/Jacobi identities give the three-twist relation. Explicit maps over K(cuberoot(3)) identify the Hecke twists; Fité–Sutherland is cited for finite-class conditioning. |
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

- Newton–Thorne, [Symmetric power functoriality for holomorphic modular forms, II](https://arxiv.org/abs/2009.07180), Theorem A.
- Ramakrishnan, [Modularity of the Rankin–Selberg L-series, and multiplicity one for SL(2)](https://arxiv.org/abs/math/0007203), Theorem M; Proposition 2.3.1(5), printed page 54; Theorem 4.2.1, pages 93–94, for the regrouping of tensor factors.
- Gillespie–Ji, [A Prime number theorem for Rankin–Selberg L-functions over number fields](https://arxiv.org/abs/0910.3660v1), Theorem 1.1 in v1. Published under the plural title *Prime number theorems*, Science China Mathematics 54 (2011), 35–46, [DOI](https://doi.org/10.1007/s11425-010-4137-x).
- Chen–Shen, [On Effective Sato–Tate Distributions for Surfaces Arising from Products of Elliptic Curves](https://arxiv.org/abs/2309.08848v1), Section 3.3 and Theorem 3.6 in v1, for the CM Hecke character. Published without the initial *On*, Research in the Mathematical Sciences 11 (2024), article 69, [DOI](https://doi.org/10.1007/s40687-024-00470-x).
- Sutherland, [Lecture 13](https://math.mit.edu/classes/18.783/2023/LectureNotes13.pdf) and [Lecture 21](https://math.mit.edu/classes/18.783/2023/LectureNotes21.pdf), for ordinary endomorphism fields and modularity background.

The manuscript states the exact hypotheses used, including fixed objects,
unitary normalization, cuspidality, self-contragredience and pole exclusions.
Its identity proofs and finite certificates do not stand in for these theorems.
The supplied priority audit identified exact antecedents, now cited and
derived at the point of use. Neither it nor this record establishes worldwide
priority for the assembled count or the numerical moment tuples.

## Finite verification

The retained elliptic package checks 2,625 large trace certificates, the small
direct-count identity cases, and the exact examples in the manuscript. The
moment package checks 207 small primes, 1,242 direct elliptic point counts,
621 CM/Jacobi comparisons, 621 Fermat comparisons, nine large CM trace
certificates, and 315 analytic routing cases for the 70 monomials through
degree four. Routing enumeration checks coverage, not automorphy or a prime
number theorem.
The additional hypergeometric check covers 2,360 four-position cases through
prime 97, including 49 cases with A=0 and 1,589 nontrivial cubic phases. It
compares literal sums, the fractional linear identity, the normalized Greene
sum, direct elliptic traces, opposite parameters, and the 3-torsion coordinate
change. The archived geometric certificates are not rewritten.

The public count oracle checks six curves. The third non-CM certificate is
transported from the archived fractional-parameter model to the integer model
by (x,y) → (9x,27y), and is then verified on that claimed integer curve. It
does not infer the count from an unchecked fast Jacobi-sum output.

The remaining publication priority is an independent, line-by-line review of
the analytic section and its interfaces with the cited theorems, followed by
an independent review of the algebraic reduction. The PR is a research draft;
merging it is not journal submission or mathematical peer review.
