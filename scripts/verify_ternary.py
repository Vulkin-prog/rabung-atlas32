"""Verify all one million stored profiles, positions, counts and doubling graphs.

The four fixed indices per block repeat exact full-run CPU colour checks.
This is not a new GPU census and does not recompute runs shorter than 21.
"""
import argparse
from collections import Counter
import csv
import io
import json
from pathlib import Path
import time

from ternary_graph import graph
from ternary_io import Dataset, DEFAULT_DATA, decompress_exact, require, sha


def threshold_table(summary):
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(['stratum', 'profiles', 'threshold', 'ordinary_runs', 'canonical_runs',
                     'zero_ordinary_profiles', 'ordinary_factorial_moment_sum'])
    for s in range(4):
        for k in range(18, 28):
            item = summary['threshold_counts_by_stratum'][str(s)][str(k)]
            writer.writerow([s, 250000, k, *[item[key] for key in
                             ('ordinary_runs', 'canonical_runs', 'zero_ordinary_profiles',
                              'ordinary_factorial_moment_sum')]])
    for k in range(18, 28):
        items = [summary['threshold_counts_by_stratum'][str(s)][str(k)] for s in range(4)]
        writer.writerow(['all', 1000000, k, *[sum(item[key] for item in items) for key in
                         ('ordinary_runs', 'canonical_runs', 'zero_ordinary_profiles',
                          'ordinary_factorial_moment_sum')]])
    return stream.getvalue()


def exact_run_check(profile):
    p, e = profile.p, (profile.p-1)//3
    require((pow(2, e, p) == 1) == (profile.s % 2 == 0), 'Exact cube class check failed')
    root = None
    for x, length, colour in profile.runs:
        numeric = pow(x, e, p)
        if colour == 0:
            require(numeric == 1, 'Neutral colour is not a cube')
        else:
            candidate = numeric if colour == 1 else numeric*numeric % p
            require(candidate != 1 and pow(candidate, 3, p) == 1,
                    'Invalid cubic character root')
            require(root is None or root == candidate, 'Inconsistent phase labels in this profile')
            root = candidate
        require(all(pow(t, e, p) == numeric for t in range(x, x+length)),
                'Stored run is not monochromatic')
        require(x == 1 or pow(x-1, e, p) != numeric, 'Stored run is not left-maximal')
        require(pow(x+length, e, p) != numeric, 'Stored run is not right-maximal')
    return len(profile.runs)


def verify(directory=DEFAULT_DATA):
    began = time.monotonic()
    dataset = Dataset(directory)
    for name in dataset.files:
        dataset.read(name)
    covered = [0]*4
    seen_primes = set()
    # Integer totals permit exact reproduction of counts and factorial moments.
    sums = [[Counter() for k in range(18, 28)] for s in range(4)]
    colours = [[0]*3 for s in range(4)]
    bins = [[0]*16 for s in range(4)]
    echoes = [[0]*2 for s in range(4)]
    lengths = [Counter() for s in range(4)]
    graphs = [Counter() for s in range(4)]
    histograms = [{key: Counter() for key in ('M', 'L', 'K')} for s in range(4)]
    witnesses = []
    independent_checks = profiles_count = canonical18 = canonical21 = a21 = ab21 = 0
    for entry in sorted(dataset.blocks.values(), key=lambda v: (v['s'], v['start'])):
        s, start, n = entry['s'], entry['start'], entry['n']
        require(start == covered[s], 'Block coverage gap or overlap')
        meta = json.loads(dataset.read(entry['metadata']))
        for profile in dataset.read_block(entry):
            require(profile.p not in seen_primes, 'Prime occurs in more than one stratum')
            seen_primes.add(profile.p)
            profiles_count += 1
            canonical18 += profile.canonical_count(18)
            canonical21 += len(profile.runs)
            a21 += profile.M < 21
            ab21 += profile.K <= 21
            for key, value in zip(('M', 'L', 'K'), (profile.M, profile.L, profile.K)):
                histograms[s][key][str(value)] += 1
            for j in range(10):
                number = profile.fields[3+j]
                sums[s][j]['ordinary_runs'] += number
                sums[s][j]['canonical_runs'] += profile.canonical_count(18+j)
                sums[s][j]['zero_ordinary_profiles'] += number == 0
                sums[s][j]['ordinary_factorial_moment_sum'] += number*(number-1)
            for j in range(3):
                colours[s][j] += profile.fields[13+j]
            for j in range(16):
                bins[s][j] += profile.fields[18+j]
            for j in range(2):
                echoes[s][j] += profile.fields[16+j]
            if profile.index-start in (0, 1024, 2048, 3072):
                independent_checks += exact_run_check(profile)
            nodes = [list(run) for run in profile.runs if profile.kind(run) == 'ordinary']
            lengths[s].update(str(run[1]) for run in nodes)
            edges, cycles = graph(profile.p, nodes, s % 2 == 0)
            graphs[s]['primes'] += 1
            graphs[s]['ordinary_runs21'] += len(nodes)
            if edges:
                graphs[s]['primes_with_edge'] += 1
                graphs[s]['edges'] += len(edges)
                graphs[s]['primes_with_cycle'] += bool(cycles)
                graphs[s]['ordinary_runs_in_weak_edge_graph'] += len(
                    {v for edge in edges for v in (edge['source'], edge['target'])})
                for cycle in cycles:
                    graphs[s][f'cycles_length{len(cycle)}'] += 1
                witnesses.append({'s': s, 'index': profile.index, 'p': profile.p,
                                  'M': profile.M, 'L': profile.L,
                                  'central_full_length': profile.fields[2], 'cube2': s % 2 == 0,
                                  'central_eligible21': bool(profile.parent[3] & 2),
                                  'ordinary_nodes21': nodes, 'edges': edges,
                                  'cycles_length1to3': cycles,
                                  'source_science_sha256': meta['science_sha256']})
        covered[s] += n
    require(profiles_count == dataset.manifest['profiles'] and covered == [250000]*4,
            'The full fixed million is not covered')
    expected = dataset.manifest['expected']
    complete = json.loads(dataset.read('COMPLETE.json'))
    require(complete['status'] == 'PASS_FULL_EXISTING_MILLION_COUNT_EXTENSION' and
            complete['n'] == profiles_count and complete['protocol_sha256'] == dataset.manifest['protocol_sha256'],
            'Scientific completion metadata differs')
    require(canonical18 == complete['canonical_runs18'] == expected['canonical_runs18'],
            'Canonical run total at 18 differs')
    require(canonical21 == complete['canonical_runs21'] == expected['canonical_runs21'],
            'Canonical positions total at 21 differs')
    require((a21, ab21) == (expected['a21'], expected['ab21']), 'Criterion (a) or (a)+(b) total differs')
    old_analysis = json.loads(dataset.read('analysis/MILLION_COUNTS_ANALYSIS.json'))
    require(old_analysis['n'] == profiles_count and old_analysis['canonical_positions21'] == canonical21,
            'Historical count audit totals differ')
    require(old_analysis['original_events'] == {'a21': a21, 'ab21': ab21}, 'Historical criteria differ')
    require(old_analysis['protocol_sha256'] == dataset.manifest['protocol_sha256'], 'Historical protocol differs')
    require(old_analysis['complete_sha256'] == complete['original_complete_sha256'], 'Historical completion anchor differs')
    require(independent_checks == old_analysis['repeated_full_run_CPU_checks'], 'Fixed CPU check coverage differs')
    require({str(s): colours[s] for s in range(4)} == old_analysis['colour_counts18_by_stratum'],
            'Historical colour totals differ')
    require({str(s): bins[s] for s in range(4)} == old_analysis['position_bins18_by_stratum'],
            'Historical position-bin totals differ')
    require({str(s): dict(lengths[s]) for s in range(4)} == old_analysis['ordinary_run_length_histograms21'],
            'Historical ordinary-length histograms differ')
    for row in old_analysis['rows']:
        k = row['k']
        for stratum in row['strata']:
            s = stratum['s']
            regenerated = sums[s][k-18]
            require(stratum['n'] == 250000 and stratum['N'] == regenerated['ordinary_runs']
                    and stratum['void'] == regenerated['zero_ordinary_profiles']
                    and stratum['factorial'] == regenerated['ordinary_factorial_moment_sum'],
                    'Historical threshold counts or factorial moments differ')
    historical_blocks = {v['name']: v for v in old_analysis['audited_blocks']}
    require(len(historical_blocks) == len(dataset.blocks), 'Historical audit block count differs')
    for entry in dataset.blocks.values():
        old = historical_blocks[entry['name']]
        require((old['s'], old['start'], old['n']) == (entry['s'], entry['start'], entry['n'])
                and old['metadata_sha256'] == entry['original_metadata_sha256']
                and old['packed_sha256'] == dataset.files[entry['pack']]['sha256'],
                'Historical audited block anchor differs')
    for s in range(4):
        for j, multiplier in enumerate(('2', '3')):
            require(echoes[s][j] == old_analysis['echoes18_by_stratum'][str(s)][multiplier]['observed'],
                    'Historical first-window echo count differs')
    old_graphs = json.loads(dataset.read('analysis/CRITICAL_GRAPHS.json'))
    require(old_graphs['n'] == profiles_count and old_graphs['summary_by_stratum'] ==
            {str(s): dict(graphs[s]) for s in range(4)}, 'Reconstructed graph totals differ')
    require(old_graphs['input_count_audit_sha256'] ==
            dataset.files['analysis/MILLION_COUNTS_ANALYSIS.json']['source_sha256'],
            'Historical graph input anchor differs')
    witness_blob = dataset.read('analysis/CRITICAL_GRAPH_WITNESSES.json.zlib')
    require(sha(witness_blob) == old_graphs['witnesses_sha256'], 'Original graph witness hash differs')
    raw_witnesses = decompress_exact(witness_blob)
    require(sha(raw_witnesses) == old_graphs['raw_witnesses_sha256'], 'Decoded graph witness hash differs')
    require(witnesses == json.loads(raw_witnesses), 'Reconstructed graph positions, phases or cycles differ')
    graph_edges = sum(v['edges'] for v in graphs)
    require(len(witnesses) == old_graphs['witness_prime_records'] == expected['graph_witness_primes']
            and graph_edges == expected['graph_edges'], 'Graph witness or edge total differs')
    summary = {'profiles': profiles_count, 'blocks': len(dataset.blocks),
               'canonical_runs18': canonical18, 'canonical_runs21': canonical21,
               'a21': a21, 'ab21': ab21,
               'threshold_counts_by_stratum': {
                   str(s): {str(k): dict(sums[s][k-18]) for k in range(18, 28)} for s in range(4)},
               'profile_histograms_by_stratum': {str(s): histograms[s] for s in range(4)},
               'ordinary_colour_counts18_by_stratum': {str(s): colours[s] for s in range(4)},
               'ordinary_position_bins18_by_stratum': {str(s): bins[s] for s in range(4)},
               'first_window_echoes18_by_stratum': {str(s): dict(zip(('2', '3'), echoes[s])) for s in range(4)},
               'ordinary_length_histograms21': {str(s): dict(lengths[s]) for s in range(4)},
               'doubling_graphs_by_stratum': {str(s): dict(graphs[s]) for s in range(4)},
               'graph_witness_primes': len(witnesses), 'graph_edges': graph_edges}
    if 'summary.json' in dataset.files:
        require(summary == json.loads(dataset.read('summary.json')), 'Publication integer summary differs')
    if 'threshold_counts.csv' in dataset.files:
        require(threshold_table(summary) == dataset.read('threshold_counts.csv').decode(),
                'Publication threshold table differs')
    return {'status': 'PASS_FULL_FIXED_TERNARY_MILLION',
            'files_verified': len(dataset.files), 'summary': summary,
            'exact_full_run_CPU_checks': independent_checks,
            'all_parent_M_L_match': True, 'all_graphs_reconstructed_and_match': True,
            'seconds': time.monotonic()-began}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--report', type=Path, help='Write the detailed verification result to this new file')
    parser.add_argument('--summary', type=Path, help='Write the regenerated integer summary to this new file')
    parser.add_argument('--table', type=Path, help='Write the regenerated threshold table to this new CSV file')
    args = parser.parse_args()
    for path in (args.report, args.summary, args.table):
        if path is not None:
            require(not path.exists() and not path.resolve().is_relative_to(args.data.resolve()),
                    'Outputs must be new files outside the immutable dataset')
    result = verify(args.data)
    for path, value in ((args.report, result), (args.summary, result['summary']),
                        (args.table, threshold_table(result['summary']))):
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('x') as stream:
                stream.write(value if isinstance(value, str) else json.dumps(value, indent=2) + '\n')
    short = {k: v for k, v in result.items() if k != 'summary'}
    short.update({k: result['summary'][k] for k in ('profiles', 'blocks', 'canonical_runs18',
                                                   'canonical_runs21', 'a21', 'ab21',
                                                   'graph_witness_primes', 'graph_edges')})
    print(json.dumps(short))


if __name__ == '__main__':
    main()
