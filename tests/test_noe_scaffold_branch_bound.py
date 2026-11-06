"""Source-only exact box proof for the Val69 local-scaffold conflict."""
import json
from pathlib import Path
from dataclasses import replace
from foldbounds.noe_scaffold_branch_bound import (
    default_val69_branch_bound, verify_branch_bound_no_go)

DATA=Path(__file__).resolve().parents[1]/'data/val69_local_noe_dual.json'

def load():
    return json.loads(DATA.read_text())


def test_source_pinned_stronger_exact_interval_result():
    c=verify_branch_bound_no_go(load(),default_val69_branch_bound())
    assert c is not None
    assert c.delta_mA==199
    assert c.objective_upper_mA2==2635000
    assert c.necessary_lower_mA2==2637307
    assert 0<c.nodes_h<=10000 and 0<c.nodes_hb<=10000
    assert c.eliminated_geometry_h>0 and c.eliminated_geometry_hb>0


def test_forged_box_or_budget_never_passes():
    d=load();p=default_val69_branch_bound()
    assert verify_branch_bound_no_go(d,replace(p,max_nodes_per_endpoint=2)) is None
    assert verify_branch_bound_no_go(d,replace(p,h=replace(p.h,root_mA=((-1,1),)*3))) is None
    assert verify_branch_bound_no_go(d,replace(p,h=replace(p.h,cap_mA2=-1))) is None
    assert verify_branch_bound_no_go(d,replace(p,h=replace(p.h,quartet=('N','N','CG1','CG2')))) is None
    assert verify_branch_bound_no_go(d,replace(p,hb=replace(p.hb,cap_mA2=100))) is None


def test_inconclusive_conditions_do_not_imply_feasibility():
    d=load();p=default_val69_branch_bound()
    assert verify_branch_bound_no_go(d,replace(p,delta_mA=201)) is None
    altered=load()
    altered['processed_restraint_upper_mA']=3500
    assert verify_branch_bound_no_go(altered,p) is None
