"""Frozen scientific functions; no GPU, launch or host monitoring side effects."""
import itertools


def graph(p, nodes, cube):
    half = (p-1)//2
    breakpoint = half//2
    edges, adj = [], set()
    for i, (x, length, colour) in enumerate(nodes):
        end = x+length-1
        for j, (y, target_length, target_colour) in enumerate(nodes):
            target_end = y+target_length-1
            # Both monotone branches, retaining parity and an actual integer witness.
            positive = (max(x, (y+1)//2), min(end, breakpoint, target_end//2))
            negative = (max(x, breakpoint+1, (p-target_end+1)//2), min(end, (p-y)//2))
            witness = next((lo for lo, hi in (positive, negative) if lo <= hi), None)
            if witness is None:
                continue
            image = min(2*witness, p-2*witness)
            assert x <= witness <= end and y <= image <= target_end
            delta = (target_colour-colour) % 3
            assert (delta == 0) == cube
            adj.add((i, j))
            edges.append({'source': i, 'target': j, 'source_position': witness,
                          'image_position': image, 'colour_shift': delta})
    cycles = [[i] for i in range(len(nodes)) if (i, i) in adj]
    cycles += [[i, j] for i, j in itertools.combinations(range(len(nodes)), 2)
               if (i, j) in adj and (j, i) in adj]
    for i, j, k in itertools.combinations(range(len(nodes)), 3):
        for order in [(i, j, k), (i, k, j)]:
            if all((order[t], order[(t+1) % 3]) in adj for t in range(3)):
                cycles.append(list(order))
    if edges:
        exponent = (p-1)//3
        character2 = pow(2, exponent, p)
        assert (character2 == 1) == cube
        root = None
        for x, length, c in nodes:
            if c:
                numeric = pow(x, exponent, p)
                root = numeric if c == 1 else numeric*numeric % p
                assert root != 1 and pow(root, 3, p) == 1
                break
        if root is None:
            assert cube
            roots = [1, None, None]
            phase = 0
        else:
            roots = [1, root, root*root % p]
            phase = roots.index(character2)
        for edge in edges:
            assert edge['colour_shift'] == phase
            assert pow(edge['source_position'], exponent, p) == roots[nodes[edge['source']][2]]
            assert pow(edge['image_position'], exponent, p) == roots[nodes[edge['target']][2]]
        assert all(len(cycle)*phase % 3 == 0 for cycle in cycles)
    return edges, cycles
