"""Reproduce a conditional real-6NBN NOE/local-geometry conflict and counterpoint.

No new data are fetched; pinned source is recorded in local JSON.
"""
import json
from pathlib import Path
from foldbounds.noe_scaffold_dual import verify_no_go,verify_local_geometry_primal,with_delta

ROOT=Path(__file__).resolve().parents[1]
raw=json.loads((ROOT/'data/val69_local_noe_dual.json').read_text())
proof=verify_no_go(raw)
primal=verify_local_geometry_primal(raw)
if proof is None or primal is None:
    raise SystemExit('ABSTAIN: certificate failed')
print('Original PDB blob:',raw['source_pdb_git_blob'])
print('Original NEF blob:',raw['source_nef_git_blob'])
print('Processed restraint:',raw['processed_restraint_id'],'upper',raw['processed_restraint_upper_mA']/1000,'Å')
print('Rigorous NO-GO when all 14 local distances change <=',proof.delta_mA/1000,'Å')
print('Required directional dot:',proof.required_mA2,'mÅ²')
print('Exact upper directional dot:',proof.support_mA2,'mÅ²')
print('Strict gap:',proof.strict_margin_mA2,'mÅ²')
print('Constructively feasible local distance cap:',primal.worst_distance_change_ticks/1e6,'Å')
print('14 local distance constraints checked:',primal.checked_local_pairs)
for d in (50,79,100,103,104):
    print('delta',d/1000,'Å; NO-GO proven?',verify_no_go(with_delta(raw,d)) is not None)
print('IMPORTANT: heavy anchors fixed by assumption, model distances not measured restraints.')
