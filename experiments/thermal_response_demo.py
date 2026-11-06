"""Reproducible exact-rational thermodynamic-response benchmark.

Usage: python experiments/thermal_response_demo.py
No optional scientific dependencies or network access required.
"""
from fractions import Fraction

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.thermal_response import thermal_response_bounds


def first_walks(n, count):
    """Standalone deterministic demonstrator: enumerate at most count SAWs."""
    result = []
    def visit(path):
        if len(result) >= count:
            return
        if len(path) == n:
            result.append(path)
            return
        x, y = path[-1]
        for point in ((x + 1, y), (x, y + 1), (x - 1, y), (x, y - 1)):
            if point not in path:
                visit(path + (point,))
    visit(((0, 0), (1, 0)))
    return tuple(result)


def main():
    for sequence in ('HHHHHH', 'HPHHPHHH', 'HHHHHHHH'):
        q = Fraction(3, 2)
        threshold = 1
        print(f'\n{sequence}; q={q}; event k>={threshold}')
        for depth in (2, len(sequence) - 1, len(sequence)):
            result = thermal_response_bounds(build_frontier(sequence, q, threshold, depth))
            print(f'  depth={depth}: mean=[{result.mean_lower}, {result.mean_upper}]'
                  f'  variance=[{result.variance_lower}, {result.variance_upper}]'
                  f'  event_response=[{result.event_susceptibility_lower},'
                  f' {result.event_susceptibility_upper}]')
        witness_subset = first_walks(len(sequence), 40)
        with_witnesses = thermal_response_bounds(
            build_frontier(sequence, q, threshold, 2),
            witnesses=witness_subset)
        print(f'  depth=2 plus {len(witness_subset)} checked witnesses:'
              f' mean=[{with_witnesses.mean_lower}, {with_witnesses.mean_upper}]'
              f' variance=[{with_witnesses.variance_lower}, {with_witnesses.variance_upper}]'
              f' event_response=[{with_witnesses.event_susceptibility_lower},'
              f' {with_witnesses.event_susceptibility_upper}]')
        spectrum = exact_spectrum(sequence)
        print(f'  independent reference density of states: {spectrum}')


if __name__ == '__main__':
    main()
