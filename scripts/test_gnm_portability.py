"""Hermetic and optional frozen source-only cross-protein portability tests."""
import os
from pathlib import Path
import pytest
from gnm_portability_audit import source, contacts, integer_graph, audit


def test_grounded_path_matrix():
    A, idx=integer_graph([(0,1),(1,2)],3,2)
    assert idx=={0:0,1:1}
    assert A==[[1,-1],[-1,2]]


def test_exact_cutoff_and_graph_invariance():
    import numpy as np
    xyz=np.array([[0.,0.,0.],[10.,0.,0.],[10.000000000000002,0.,0.]],dtype=float)
    assert (0,1) in contacts(xyz)
    assert (0,2) not in contacts(xyz)


def test_forged_nested_archive_rejected(tmp_path):
    src=tmp_path/'fake_prospective.zip'
    src.write_bytes(b'not the frozen archive')
    with pytest.raises(ValueError, match='outer ZIP'):
        source(src,'CspB')


def test_cross_protein_source_only_when_archive_provided():
    path=os.environ.get('FOLDBOUNDS_PROSPECTIVE_ZIP')
    if not path:
        pytest.skip('Original signed frozen prospective ZIP not provided')
    rows=[audit(path,name) for name in ('CspB','Crambin')]
    assert [x['contact_edges'] for x in rows]==[543,334]
    assert [x['exact_ties'] for x in rows]==[9,0]
    assert all(x['matrix_identity_checked'] and
               x['frozen_strict_order_violations']==0 and
               x['max_float_difference']<1e-10 and
               not x['target_data_accessed'] for x in rows)