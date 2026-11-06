"""Unit and opt-in original archive regression for source-only portability."""
from __future__ import annotations
import os
from fractions import Fraction
from pathlib import Path
import numpy as np
import pytest
from gnm_coupled_crossprotein import contact_bounds,audit


def test_exact_coordinate_uncertainty_envelopes():
    coords=np.array([[0.,0.,0.],[10.,0.,0.],[9.5,0.,0.]])
    sure,nominal,possible=contact_bounds(coords,Fraction(1,20))
    assert (0,1) in nominal and (0,1) in possible
    assert (0,1) not in sure  # +/-0.05 in both coordinates can cross cutoff
    assert set(sure).issubset(set(nominal))
    assert set(nominal).issubset(set(possible))


def test_source_only_cspb_crambin_when_archive_present():
    path=os.environ.get('FOLDBOUNDS_PROSPECTIVE_ZIP')
    if not path:pytest.skip('Original frozen source archive not provided to CI')
    c=audit(path,'CspB');r=audit(path,'Crambin')
    assert (c['n_uncertain_contacts'],r['n_uncertain_contacts'])==(52,36)
    assert (c['additional_certified_strict_orders'],r['additional_certified_strict_orders'])==(4656,1643)
    assert c['n_nominal_exact_tied_comparisons']==9
    assert r['n_nominal_exact_tied_comparisons']==0
    assert c['target_NMR_data_accessed'] is False and r['target_NMR_data_accessed'] is False