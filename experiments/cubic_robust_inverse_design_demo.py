"""Prove one globally unique 16-mer cubic HP native orbit under energy errors.

Run after installing FoldBounds and with g++ available:
  python experiments/cubic_robust_inverse_design_demo.py --output-dir ./proof_out

By default generates all three cut-tree files. To reuse existing files, use:
  --proofs-dir /path/to/existing_cut_proofs
where files must be named threshold10.cut, threshold11.cut, threshold12.cut.

q is dimensionless Boltzmann weight per QUARTER-energy unit. Output PDBs
are SYNTHETIC lattice traces, NOT experimentally known protein structures.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

from foldbounds.cubic_cut_certificate import certify_native_gap,contact_score
from foldbounds.cubic_robust_inverse_design import (
    _d2,design_robust_unique_orbit,verify_ground_map_cut
)

SEQUENCE='HHHPPHHHHHHHPPHH'
EXCITED_WITNESS=(
    (0,0,0),(1,0,0),(2,0,0),(3,0,0),
    (3,1,0),(2,1,0),(1,1,0),(0,1,0),
    (0,1,1),(1,1,1),(2,1,1),(3,1,1),
    (3,0,1),(2,0,1),(1,0,1),(0,0,1),
)
TARGET_MAP=15_668_042_105_480


def _pdb(sequence,coordinates,kind):
    lines=[
        'HEADER    SYNTHETIC THREE-DIMENSIONAL CUBIC HP REFERENCE',
        'REMARK   1 NOT AN EXPERIMENTAL OR PHYSICAL PROTEIN STRUCTURE',
        'REMARK   2 H=ALA P=GLY PLACEHOLDER RESIDUES FOR PYMOL ONLY',
        f'REMARK   3 RESIDUES {sequence}',
        f'REMARK   4 CONFORMATION {kind}',
        'REMARK   5 FIRST BOND ANCHORED; 3.8 ANGSTROM VISUAL SPACING',
    ]
    for i,(letter,(x,y,z)) in enumerate(zip(sequence,coordinates),1):
        residue='ALA' if letter=='H' else 'GLY'
        lines.append(f'ATOM  {i:5d}  CA  {residue} A{i:4d}    '
                     f'{3.8*x:8.3f}{3.8*y:8.3f}{3.8*z:8.3f}'
                     f'{1.0:6.2f}{0.0:6.2f}           C')
    lines.extend(f'CONECT{i:5d}{i+1:5d}' for i in range(1,len(sequence)))
    lines.extend(['TER','END'])
    return '\n'.join(lines)+'\n'


def compile_checkers(src:Path,build:Path):
    if not shutil.which('g++'):
        raise RuntimeError('g++ required for independent cut-tree checking')
    tools={}
    for name in ('cubic_cut_checker','cubic_ground_map_checker','cubic_cut_producer'):
        exe=build/name
        subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Werror',
                        str(src/(name+'.cpp')),'-o',str(exe)],
                       check=True,capture_output=True,text=True,timeout=60)
        tools[name]=exe
    return tools


def run_demo(output:Path,proofs_dir:Path|None=None)->dict:
    output.mkdir(parents=True,exist_ok=True)
    tools=compile_checkers(Path(__file__).resolve().parent,output)
    proof_dir=proofs_dir or output
    if not proof_dir.is_dir():
        raise ValueError('proof directory does not exist')
    cuts=tuple(proof_dir/f'threshold{t}.cut' for t in (10,11,12))
    if proofs_dir is None:
        for t,file in zip((10,11,12),cuts):
            subprocess.run([str(tools['cubic_cut_producer']),SEQUENCE,str(t),
                            str(file),'geometry'],check=True,
                           capture_output=True,text=True,timeout=180)
    if contact_score(SEQUENCE,EXCITED_WITNESS)!=9:
        raise RuntimeError('incorrect first-excited-level witness')
    print('Replay 3 independent full-tree proofs (slow but exhaustive over all cuts)...',flush=True)
    gap=certify_native_gap(cuts,tools['cubic_cut_checker'],EXCITED_WITNESS,max_seconds=900)
    if (gap.sequence,gap.ground_contacts,gap.ground_conformations,gap.gap)!=(SEQUENCE,11,2032,2):
        raise RuntimeError('unexpected certified 3D excitation gap')
    print('Independently replay ground-map tree and verify ALL 2032 surviving coordinate paths...',flush=True)
    maps=verify_ground_map_cut(cuts[1],tools['cubic_ground_map_checker'],max_seconds=600)
    design=design_robust_unique_orbit(
        gap,maps,selected_map=TARGET_MAP,
        distance_pair=(0,12),squared_distance=2
    )
    assert len(maps.multiplicities)==73
    assert design.target_orbit_size==8
    assert sum(n for _,_,n in design.remaining_high)==2024
    # Crucial contact-only no-go: among *K=11* folds, every map has at
    # least 16 geometries (i.e. at least two distinct 8-fold D4 orbits).
    assert min(count for _,count in maps.multiplicities)==16
    target=design.target_path
    other=next(path for mask,path in maps.entries
               if mask==design.designated_contact_mask
               and _d2(path,0,12)!=design.target_squared_distance)
    assert _d2(other,0,12)==4
    (output/'selected_native_orbit.pdb').write_text(
        _pdb(SEQUENCE,target,'SELECTED NATIVE D4 ORBIT, RESIDUE DISTANCE SQUARED 2'))
    (output/'contact_identical_competitor.pdb').write_text(
        _pdb(SEQUENCE,other,'EXCLUDED CONTACT-IDENTICAL ORBIT, DISTANCE SQUARED 4'))
    (output/'view_in_pymol.pml').write_text('''load selected_native_orbit.pdb, selected_native
load contact_identical_competitor.pdb, excluded_orbit
hide everything
show spheres, name CA
show sticks, name CA
color orange, resn ALA
color cyan, resn GLY
translate [15,0,0], excluded_orbit, camera=0
zoom all
''')
    targets={}
    for error in (0,1):
        for percentage,p in ((50,Fraction(1,2)),(90,Fraction(9,10)),(99,Fraction(99,100))):
            minimum=design.minimum_integer_q(p,radius=error)
            targets[f'error_{error}_percentage_{percentage}']={
                'minimum_integer_q':minimum,
                'certified_lower_at_q':str(design.lower_orbit_probability(Fraction(minimum),radius=error)),
                'certified_lower_at_q_decimal':float(design.lower_orbit_probability(Fraction(minimum),radius=error)),
            }
    manifest={
        'model':'synthetic 3D simple cubic HP model, anchored first bond +x',
        'not_a_real_protein_prediction':True,
        'sequence':SEQUENCE,'all_ground_contact_maps':len(maps.multiplicities),
        'ground_conformations':gap.ground_conformations,
        'exact_uniform_contact_excitation_gap':gap.gap,
        'smallest_contact_map_class':min(n for _,n in maps.multiplicities),
        'selected_contact_mask':TARGET_MAP,
        'selected_residue_distance_pair_zero_based':[0,12],
        'selected_squared_lattice_distance':2,
        'excluded_equal_contact_map_squared_distance':4,
        'selected_geometry_symmetry_orbit_size':8,
        'energy_units':'quarter-energy; preferred H-H contact=8, other H-H contact=4, squared-distance bonus=2',
        'arbitrary_correlated_contact_energy_uncertainty':'each per-contact energy perturbed within ±1 quarter-energy unit',
        'distinct_geometry_contact_only_no_go':
            'any contact-map-only energy selecting K=11 has at least 16 equal-energy conformations',
        'one_extra_distance_restraint_design':'selects a single D4 symmetry orbit among every 3D conformation',
        'native_contact_score':11,'native_designed_energy':90,
        'remaining_exact_high_contact_energy_class_counts':[
            {'missing_target_contacts':a,'distance_restraint_bonus':b,'conformations':cnt}
            for a,b,cnt in design.remaining_high
        ],
        'verified_tree_excited_upper_bins':[
            {'upper_contacts':k,'max_count':c} for k,c in maps.upper_bins
        ],
        'robust_occupancy_targets':targets,
        'verified_cut_sha256':{str(t):sha256(path.read_bytes()).hexdigest()
                               for t,path in zip((10,11,12),cuts)},
    }
    (output/'design_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('VERIFIED selected exact native orbit and *all correlated* contact-error robustness.')
    for key,item in targets.items():
        print(f'  {key}: q>={item["minimum_integer_q"]}, P>={item["certified_lower_at_q_decimal"]:.8%}')
    print('PyMOL witnesses and manifest:',output)
    print('NO REAL PROTEIN CLAIM. No Lean refinement of this C++/Python checking path.')
    return manifest


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True,type=Path)
    p.add_argument('--proofs-dir',type=Path)
    a=p.parse_args()
    run_demo(a.output_dir,a.proofs_dir)
