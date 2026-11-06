"""Reproduce conditional moving-anchor real-structure NOE impossibility."""
from pathlib import Path
import json
from foldbounds.noe_scaffold_dual import from_json
from foldbounds.noe_scaffold_branch_bound import default_val69_branch_bound
from foldbounds.noe_moving_scaffold import (
    HeterogeneousScaffoldBudget, UniformMovingScaffoldProof,
    verify_uniform_moving_scaffold,verify_heterogeneous_moving_scaffold,
)

raw=json.loads((Path(__file__).resolve().parents[1]/'data'/'val69_local_noe_dual.json').read_text())
base=default_val69_branch_bound()
for eta,delta in ((0,199),(50,149),(100,99),(150,49)):
    r=verify_uniform_moving_scaffold(raw,UniformMovingScaffoldProof(eta,delta,base))
    assert r is not None
    print(f'CERTIFIED NO-GO: heavy shift <= {eta/1000:.3f} A and local shell change <= {delta/1000:.3f} A')
used=['N','CA','CB','CG2','O']
b=HeterogeneousScaffoldBudget({a:40 for a in used},
    {h:{a:50 for a in used} for h in ('H','HB')})
r=verify_heterogeneous_moving_scaffold(from_json(raw),b)
assert r is not None
print('HETEROGENEOUS CERTIFIED NO-GO with only', ', '.join(r.constrained_anchor_names),'constrained')
print('Exact directional cap:',r.possible_support_mA2,'< required:',r.necessary_support_mA2)
