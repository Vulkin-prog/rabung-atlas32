#!/usr/bin/env python3
"""Compare the separate CPU checker with literal power colours and small APs."""
import argparse
import json
import subprocess
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checker', type=Path)
    parser.add_argument('--max-order', type=int, choices=(3, 9), default=9)
    args = parser.parse_args()
    checker = str(args.checker.resolve())
    checked = aps = 0
    for p in range(5, 200):
        if any(p % d == 0 for d in range(2, int(p ** .5) + 1)):
            continue
        for r in range(2, min(args.max_order, p - 1) + 1):
            if (p - 1) % r:
                continue
            colours = [None] + [pow(x, (p - 1) // r, p) for x in range(1, p)]
            runs = []
            start = 1
            for x in range(2, p + 1):
                if x == p or colours[x] != colours[start]:
                    runs.append((x - start, start))
                    start = x
            M = max(v[0] for v in runs)
            w = min(s for length, s in runs if length == M)
            L = runs[0][0]
            for k in sorted({3, min(p - 1, max(3, M)), min(p - 1, max(3, M + 1)), min(p - 1, max(3, 2 * L + 2))}):
                if k >= p:
                    continue
                run = subprocess.run([checker, str(p), str(r), str(k)], capture_output=True, text=True)
                if run.returncode not in (0, 1):
                    raise RuntimeError(run.stderr)
                result = json.loads(run.stdout)
                boundary = all(len({colours[t % p] for t in range(-j, k - j) if t}) > 1 for j in range(k))
                expected = M < k and boundary
                if (result['M'], result['L'], result['w'], result['Bstar'], result['accepted']) != (M, L, w, boundary, expected):
                    raise RuntimeError(f'power-colour mismatch at {(p, r, k)}')
                if (run.returncode == 0) != expected:
                    raise RuntimeError('exit status mismatch')
                if p <= 31 and k <= 5:
                    # A concrete nonconstant assignment at the multiples of p.
                    N = (k - 1) * p + 1
                    alternate = next(v for v in colours[1:] if v != 1)
                    word = [colours[x % p] if x % p else (alternate if x == (k - 1) * p else 1) for x in range(N)]
                    valid = not any(len({word[a + j * d] for j in range(k)}) == 1
                                    for d in range(1, (N - 1) // (k - 1) + 1)
                                    for a in range(N - (k - 1) * d))
                    if valid != expected:
                        raise RuntimeError(f'literal progression mismatch at {(p, r, k)}')
                    aps += 1
                checked += 1
    for values in [('4294967296', '2', '28'), ('25', '2', '4'), ('11', '3', '4'), ('-5', '2', '3')]:
        if subprocess.run([checker, *values], capture_output=True).returncode != 2:
            raise RuntimeError('invalid input not refused')
    print(json.dumps({'status': 'PASS', 'power_colour_profiles': checked, 'literal_AP_checks': aps, 'invalid_inputs': 4}))

if __name__ == '__main__':
    main()
