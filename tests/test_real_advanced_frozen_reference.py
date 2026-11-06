"""Reproducibility gates for experimental GA95/GB95 contact topology and reporters.
These are retrospective geometric descriptors, not blind fold predictions.
"""
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]/'experiments'/'results'

def test_nmr_model_level_topological_separability():
    d=json.loads((ROOT/'ga_gb_advanced_20261020.json').read_text())
    m=d['contact_map_model_discrimination']
    assert m['cross_fold_pairs']==400
    assert m['within_GA95_pairs']==m['within_GB95_pairs']==190
    assert m['strict_model_level_separation']
    assert m['between_folds']['max']<m['within_GA95']['min']
    assert m['between_folds']['max']<m['within_GB95']['min']
    assert d['coordinate_SHA256']['2KDL']=='fb75cc55eb8b89e0cec5e2be69901b348e60916375b9c410a0573edfca1db591'
    assert d['coordinate_SHA256']['2KDM']=='1c49e5faa25b443285432d935c621151294ed3b636b0106e6b1a7f7ccf2d00ec'

def test_actual_aromatic_ring_plane_geometry_avoids_pi_stack_overclaim():
    d=json.loads((ROOT/'ga_gb_advanced_20261020.json').read_text())
    a={tuple(x['sites']):x for x in d['GB95_aromatic_geometry']}
    assert a[(30,52)]['orientation_class_counts']=={'perpendicular_plane_like':20}
    assert a[(45,52)]['orientation_class_counts']=={'close_other_orientation':20}
    assert a[(29,45)]['orientation_class_counts']=={'separated':20}

def test_recommended_distance_reporter_is_ranked_by_all_models():
    d=json.loads((ROOT/'ga_gb_distance_reporters_20261020.json').read_text())
    assert d['number_of_nonmutated_pairs_with_disjoint_12_to_40_A_NMR_distance_ranges']==349
    first=d['top_ten'][0]
    assert (first['i'],first['j'])==(11,21)
    assert first['aa_i']=='Q' and first['aa_j']=='V'
    assert first['GA95_max_A']<first['GB95_min_A']
    assert first['worst_case_separation_A']>16
    assert first['max_symmetric_CA_to_label_error_A']>8