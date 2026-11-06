"""Source-pinned, fail-closed NMR restraint repair and holdout demonstration."""
from pathlib import Path
import json
from foldbounds.noe_restraint_audit import load_pinned_reference
from foldbounds.noe_minimax_repair import certificate_from_json,verify_repair
from foldbounds.noe_holdout import five_fold_audit

root=Path(__file__).resolve().parents[1]
fixture=root/'data/pinned_unambiguous_noe_2png_6nbn.json'
wit=json.loads((root/'data/noe_minimax_repair_witnesses.json').read_text())
ref=load_pinned_reference(fixture)
for id in ('2png','6nbn'):
    result=verify_repair(ref[id],certificate_from_json(id,wit['samples'][id]))
    if result is None:raise RuntimeError('ABSTAIN: real NMR repair witness failed')
    folds=five_fold_audit(id,ref[id])
    print(f'{id}: {result.restraints} measured intervals; {result.atoms} atoms; '
          f'{result.moved_atoms} moved; optimal rho in '
          f'({result.strict_lower_radius_ticks / 1e6:.6f}, '
          f'{result.upper_radius_ticks / 1e6:.6f}] A; '
          f'{sum(f.proximity_le_8a for f in folds)} / '
          f'{sum(f.heldout for f in folds)} withheld atom-pairs source-only <=8A; '
          f'{sum(f.entire_interval_guaranteed_by_upper_only for f in folds)} full intervals guaranteed')
