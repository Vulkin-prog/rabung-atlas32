"""Query or export counts and long-run positions from the fixed ternary million."""
import argparse
import csv
import json
from pathlib import Path
import sys

from ternary_io import Dataset, DEFAULT_DATA, IntegrityError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prime', nargs='?', type=int, help='An exact prime from the fixed sample')
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--stratum', type=int, choices=range(4))
    parser.add_argument('--index', type=int, help='Zero-based index within the selected stratum')
    parser.add_argument('--export', choices=('counts', 'positions'), help='Stream the selected stratum or all strata as CSV to stdout')
    parser.add_argument('--minimum', type=int, default=21, help='Minimum length for position export (at least 21)')
    args = parser.parse_args()
    if args.export and (args.prime is not None or args.index is not None):
        parser.error('Choose one profile or a CSV export')
    if args.prime is not None and (args.index is not None or args.stratum is not None):
        parser.error('Prime lookup does not take a stratum/index')
    if args.index is not None and args.stratum is None:
        parser.error('--index requires --stratum')
    if args.minimum < 21:
        parser.error('Positions below length 21 were not retained for the million')
    if args.minimum != 21 and args.export != 'positions':
        parser.error('--minimum applies to position export')
    if args.prime is None and args.index is None and not args.export:
        parser.error('Specify a prime, --stratum with --index, or --export')
    try:
        dataset = Dataset(args.data)
        if args.export:
            writer = csv.writer(sys.stdout)
            if args.export == 'positions':
                writer.writerow(['stratum', 'index', 'p', 'start', 'length', 'colour', 'kind', 'reflected_start'])
                for profile in dataset.iter_profiles(args.stratum):
                    for run in profile.runs:
                        x, length, colour = run
                        if length >= args.minimum:
                            writer.writerow([profile.s, profile.index, profile.p, x, length,
                                             colour, profile.kind(run), profile.p-x-length+1])
            else:
                writer.writerow(['stratum', 'index', 'p', 'M', 'L', 'K', 'central_length_ge18']
                                + [f'ordinary_N{k}' for k in range(18, 28)]
                                + [f'canonical_N{k}' for k in range(18, 28)]
                                + [f'colour18_{j}' for j in range(3)] + ['echo18_2', 'echo18_3']
                                + [f'position_bin18_{j}' for j in range(16)])
                for profile in dataset.iter_profiles(args.stratum):
                    writer.writerow([profile.s, profile.index, profile.p, profile.M, profile.L, profile.K,
                                     profile.fields[2], *profile.fields[3:13],
                                     *[profile.canonical_count(k) for k in range(18, 28)],
                                     *profile.fields[13:34]])
        else:
            profile = (dataset.by_prime(args.prime) if args.prime is not None
                       else dataset.by_index(args.stratum, args.index))
            print(json.dumps(profile.as_dict(), indent=2))
    except (IntegrityError, KeyError, OSError) as exc:
        parser.exit(1, str(exc) + '\n')


if __name__ == '__main__':
    main()
