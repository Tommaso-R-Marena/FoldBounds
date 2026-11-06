"""Reproduce exact model-versus-NMR-restraint audit, without network access."""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from foldbounds.noe_restraint_audit import load_pinned_reference, audit


def main():
    source=Path(__file__).resolve().parents[1]/'data/pinned_unambiguous_noe_2png_6nbn.json'
    reference=load_pinned_reference(source)
    print('One-row atom-specific mapped NMR intervals; unsupported/ambiguous rows EXCLUDED.')
    for entry,records in sorted(reference.items()):
        print(entry, 'mapped',len(records))
        for radius in (0,10,50,100,500):
            r=audit(records,radius)
            print(f'  movement <= {radius/1000:.3f} Å per atom:'
                  f' certainly satisfied {r.certainly_satisfied}, certainly violated {r.certainly_violated},'
                  f' undetermined {r.undecided}; disjoint contradiction lower bound'
                  f' {r.minimum_atoms_requiring_more_than_epsilon_movement} atoms;'
                  f' witness IDs {r.disjoint_witness_ids}')

if __name__=='__main__':
    main()
