"""Reproduce exact anchored cubic map entropy and robust event/no-go bounds."""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from fractions import Fraction

from foldbounds.filled_box_contact_rigidity import (
    extract_nonbackbone_contacts, certify_filled_box_contacts,
    verify_filled_box_certificate, split_proper_and_mirror_orbits,
    certify_pairwise_contact_energy, verify_pairwise_contact_energy,
)
from foldbounds.cubic_ground_thermo import (
    cubic_ground_population, verify_cubic_ground_population,
    CUBE_GROUND_COUNTS, CUBIC_FREE_BLOCK_COUNTS,
)

ROOT = Path(__file__).resolve().parents[1]


def check_cpp_against_python_constants():
    if shutil.which('g++') is None:
        raise RuntimeError('g++ required for --audit-cpp')
    audit = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, arguments in (
            ('cubic_count_independent.cpp', ['3']),
            ('cubic_ground_contact_map_audit.cpp', ['3']),
            ('cubic_short_block_oracle.cpp', ['10']),
        ):
            executable = Path(tmp)/name.replace('.cpp', '')
            subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra',
                            str(ROOT/'experiments'/name),'-o',str(executable)],
                           check=True,timeout=60,capture_output=True,text=True)
            output = subprocess.run([str(executable)]+arguments,check=True,
                                    timeout=60,capture_output=True,text=True).stdout
            audit[name] = output
    total = int(audit['cubic_count_independent.cpp'].strip())
    lines = audit['cubic_ground_contact_map_audit.cpp'].splitlines()
    counts = tuple(int(row.split(',')[1]) for row in
                   audit['cubic_short_block_oracle.cpp'].splitlines())
    if (total != CUBE_GROUND_COUNTS[3] or
        f'first_E_paths={total}' not in lines[0] or
        'distinct_contact_masks=103346' not in lines[0] or
        '8,103346' not in lines[2] or counts != CUBIC_FREE_BLOCK_COUNTS):
        raise AssertionError('independent C++ audits disagree')
    return {'complete_paths':total, 'distinct_contact_maps':103346,
            'multiplicity_per_contact_map':8,
            'exact_short_block_counts':counts}


def serpentine_cube():
    out = []
    for z in range(3):
        ys = range(3) if z % 2 == 0 else range(2,-1,-1)
        for row, y in enumerate(ys):
            xs = range(3) if (z + row) % 2 == 0 else range(2,-1,-1)
            for x in xs:
                out.append((x,y,z))
    assert len(out)==27 and out[:2]==[(0,0,0),(1,0,0)]
    assert len(set(out))==27
    assert all(sum(abs(a-b) for a,b in zip(v,w))==1 for v,w in zip(out,out[1:]))
    return tuple(out)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--audit-cpp',action='store_true',
                        help='compile and independently recompute every 3^3 ground path and contact map')
    args=parser.parse_args()
    path=serpentine_cube()
    contacts=extract_nonbackbone_contacts(path)
    cert=certify_filled_box_contacts(3,3,contacts)
    assert verify_filled_box_certificate(cert)
    assert path in cert.orbit
    proper,mirror=split_proper_and_mirror_orbits(cert)
    pair_coefficients = {(i,j):Fraction(((17*i+11*j)%19)-9,3)
                         for i in range(27) for j in range(i+2,27)}
    energy_check = certify_pairwise_contact_energy(cert,pair_coefficients)
    assert verify_pairwise_contact_energy(cert,pair_coefficients,energy_check)
    bounds=cubic_ground_population(3, Fraction(10),16,
                                    correction_lower=Fraction(1,2),correction_upper=2)
    assert verify_cubic_ground_population(bounds)
    result={
        'model':'anchored all-H 3D cubic lattice (not a physical protein)',
        'n':27,'max_contacts':28,'ground_paths':bounds.ground_paths,
        'ground_contact_maps':bounds.ground_contact_map_count,
        'contact_map_orbit_size':cert.orbit_size,
        'proper_orientation_count':len(proper),
        'reflected_orientation_count':len(mirror),
        'sample_contact_pairs':len(contacts),
        'arbitrary_pairwise_energy_invariant':str(energy_check.energy),
        'block_partition':bounds.total_walk_upper.blocks,
        'total_walk_cap':bounds.total_walk_upper.bound,
        'q':str(bounds.q),'event_at_least':bounds.minimum_contacts,
        'event_probability_lower_robust':str(bounds.event_probability_lower),
        'event_probability_lower_robust_decimal':float(bounds.event_probability_lower),
        'single_exact_ground_geometry_probability_upper_robust':str(bounds.single_fold_probability_upper),
        'single_exact_ground_contact_map_probability_upper_robust':str(bounds.one_ground_contact_map_probability_upper),
        'audit_status':('recomputed independent C++ oracles' if args.audit_cpp
                        else 'fixed audited counts; use --audit-cpp to regenerate'),
    }
    if args.audit_cpp:
        result['independent_cpp_audits']=check_cpp_against_python_constants()
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
