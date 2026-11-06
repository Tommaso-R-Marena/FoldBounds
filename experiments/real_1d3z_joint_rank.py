"""Actual 1D3Z ubiquitin NMR structure: exact GNM rank certificates.

The fixture contains only the C-alpha coordinates of residues 2..70,
model 1 from deposited PDB 1D3Z (source below) at exact 0.001 A precision.
The structure and its ten deposited conformers are *experimental
restraint-refined models*, not a validated dynamic trajectory.

Source: https://www.rcsb.org/structure/1D3Z
Reproducible archival mirror:
  encore-similarity/encore, examples/experimental_ubq_ensembles/
  filtered_1D3Z.pdb, git blob a438d36816e15647177e32d2eafa6aa9a0782747.

The untrusted SymPy producer is independently checked by exact integer
B*N=D*I matrix multiplication before any claim is returned.
"""
from __future__ import annotations

import argparse
import csv
import json
from fractions import Fraction
from pathlib import Path

from foldbounds.gnm_joint_rank import (
    classify_contact_boxes, joint_rank_certificate, untrusted_propose_inverse,
)

DEFAULT_FIXTURE = (Path(__file__).resolve().parents[1] /
                   'tests/fixtures/1d3z_model1_ca_milli.csv')


def load_ca(filename: Path):
    with filename.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    if len(rows) != 69 or [int(row['residue']) for row in rows] != list(range(2, 71)):
        raise ValueError('expected 69 deposited C-alpha residues 2 through 70')
    return tuple(tuple(int(row[f'{axis}_milli']) for axis in 'xyz') for row in rows)


def evaluate(eps: Fraction, *, coordinates=None):
    xyz = coordinates if coordinates is not None else load_ca(DEFAULT_FIXTURE)
    certain, possible = classify_contact_boxes(xyz, Fraction(8), eps)
    lo = untrusted_propose_inverse(len(xyz), certain)
    hi = untrusted_propose_inverse(len(xyz), possible)
    result = joint_rank_certificate(len(xyz), certain, possible, lo, hi)
    return {
        'coordinate_error_angstrom': str(eps),
        'residues': len(xyz), 'total_pairs': result.total_pairs,
        'certain_nonlocal_contacts': len(certain) - len(xyz)+1,
        'possible_nonlocal_contacts': len(possible) - len(xyz)+1,
        'separate_diagonal_certified': result.simple_certified,
        'joint_covariance_certified': result.joint_certified,
        'strictly_new_pairs': result.gained,
        'exact_witnesses_verified': 2,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--eps', nargs='+', default=['0.025', '0.05', '0.1', '0.2'],
                        help='hypothetical per-Cartesian-coordinate error in Angstrom')
    args = parser.parse_args()
    outputs = [evaluate(Fraction(eps)) for eps in args.eps]
    print(json.dumps({'source_pdb':'1D3Z','model':1,
                      'source_fixture':'1d3z_model1_ca_milli.csv',
                      'proof_scope':'all independent optional-contact assignments',
                      'baseline':'diagonal Loewner intervals',
                      'candidate':'rank-two joint covariance Loewner intervals',
                      'results':outputs},indent=2))


if __name__=='__main__':
    main()
