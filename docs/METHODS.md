# Scientific contract

The interval is [1,2^32), with every prime p and every r in {2,…,9} dividing
p−1. A profile colours x in {1,…,p−1} by its class in
F_p^×/(F_p^×)^r. Changing a primitive root only relabels colours.

M is the longest consecutive monochromatic run; L is the initial run of the
identity colour, beginning at x = 1. The witness is the least reflected start
among runs attaining M: min(s, p−(s+M−1)). There is no run crossing zero in M.

Write sigma = +1 when r divides (p−1)/2 and sigma = −1 otherwise. The threshold
stored through the codec is

```text
K(p,r) = max(M+1, 2L+2)   if sigma = +1,
K(p,r) = max(M+1, L+2)    if sigma = -1.
```

For a target k and an admissible prime p > k, criterion **(a)** is M < k.
The complete Rabung criterion **(a)+(b)** is K ≤ k, including the boundary
condition. Counts in the target tables must retain this distinction. The
catalogue is a census of this construction within the stated interval, not a
claim that every extremal van der Waerden construction has this form.

For odd p, A4 counts integers x from 1 through p−4 for which x,x+1,x+2,x+3
have the same quadratic character. The p = 2 record has A4 = 0 as a placeholder.
Small primes are counted densely. Larger primes use an exact elliptic
point-count calculation, with independent certificates recorded periodically.
The point-count algorithm and certificate checker are included in `src/census/`.

Production commits had at most 4,096 primes. GPU batches had at most 128 primes
and CPU A4 batches at most 4,096; each had a 30-second deadline. The final
scientific chains and periodic audits are retained. Execution policies can
affect throughput but do not change accepted profile values or certificates.
An engine failure or capacity refusal never licenses accepting a partial batch.

The theoretical examples and density laws are post-census investigations of
selected witnesses. Exact finite computations and algebraic identities have
different roles: rank lower bounds from finite characters are paired with
exact relation upper bounds in number fields. Conditional and all-prime
densities, and neutral and constant phases, are distinguished in the paper.
The remaining profile and bounded search limits are stated explicitly.
