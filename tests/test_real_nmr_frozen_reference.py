"""Check frozen published 2KDL/2KDM NMR metrics separately from live retrieval.
The live GitHub CI workflow separately verifies official coordinate SHA and reproduces them.
"""
import json, pathlib, pytest
from foldbounds.ga_gb_nmr import GA95,GB95,MUTATIONS

REF=pathlib.Path(__file__).resolve().parents[1]/'experiments'/'results'/'ga_gb_nmr_verified_20261020.json'

def test_authentic_pdb_reference_20_models_and_exact_mutations():
    o=json.loads(REF.read_text())
    assert o['GA95']['nmr_models']==o['GB95']['nmr_models']==20
    assert o['GA95']['sequence']==GA95 and o['GB95']['sequence']==GB95
    assert [x['position'] for x in o['comparison']['mutations']]==[20,30,45]
    assert o['comparison']['consensus_contacts_shared']==14
    assert o['comparison']['consensus_contacts_only_GA95']==57
    assert o['comparison']['consensus_contacts_only_GB95']==96
    assert o['comparison']['consensus_contact_jaccard']==pytest.approx(14/167)

def test_independent_structural_phenotype_geometry():
    o=json.loads(REF.read_text())
    ga=o['GA95'];gb=o['GB95']
    assert ga['mean_helix_like_i_i3_fraction'] > gb['mean_helix_like_i_i3_fraction']
    assert ga['mean_extended_i_i3_fraction'] < gb['mean_extended_i_i3_fraction']
    assert ga['consensus_contacts_75pct']==71 and gb['consensus_contacts_75pct']==110
    aromatic={tuple(x['sites']):x for x in gb['aromatic_ring_pairs']}
    for pair in [(29,30),(30,52),(45,52)]:
        assert aromatic[pair]['valid_models']==20
        assert aromatic[pair]['fraction_within_7_A']==1.0
    assert aromatic[(29,45)]['fraction_within_7_A']==0.0