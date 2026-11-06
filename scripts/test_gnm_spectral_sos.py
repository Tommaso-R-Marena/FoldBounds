"""Hermetic exact-integer tests for the source-only SOS spectral checker."""
from fractions import Fraction
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).parent))
from gnm_spectral_sos_independent_replay import verify_spectral_matrix

A=[[2,-1],[-1,2]]
W=[[2,-2],[0,0]]  # one exact rank-one square
P,QDEN,D=1,2,2  # certified lambda=1/2, dyadic proposal denominator=2

def test_exact_rank_one_sos():
    r=verify_spectral_matrix(A,W,P,QDEN,D)
    assert r['min_diagonal_surplus']==4
    assert r['diagonal_square_terms']==2
    for x in range(-8,9):
      for y in range(-8,9):
        quad=2*x*x-2*x*y+2*y*y
        assert Fraction(quad)>=Fraction(1,2)*(x*x+y*y)

def test_untrusted_zero_proposal_is_allowed_when_remainder_psd():
    # The checker is *not* a proximity-to-numerical-Cholesky checker.
    assert verify_spectral_matrix(A,[[0,0],[0,0]],1,2,1)['min_diagonal_surplus']>=0

@pytest.mark.parametrize('Wbad',[
    [[12,-2],[0,0]],
    [[2,-12],[0,0]],
    [[2,-2],[3,0]],
    [[2,-2]],
    [[2,-2],[0,'x']],
    [[2,-2],[0,1<<100]],
])
def test_reject_corrupted_witness(Wbad):
    with pytest.raises(ValueError):verify_spectral_matrix(A,Wbad,P,QDEN,D)

def test_reject_wrong_matrix_shape_or_asymmetry():
    with pytest.raises(ValueError):verify_spectral_matrix([[2,-1],[0,2]],W,P,QDEN,D)
    with pytest.raises(ValueError):verify_spectral_matrix([[2,-1]],W,P,QDEN,D)

def test_reject_invalid_shifts():
    for p,q,den in [(0,2,2),(1,0,2),(2,1,2),(1,2,0),(True,2,2)]:
        with pytest.raises(ValueError):verify_spectral_matrix(A,W,p,q,den)

def test_reject_false_decomposition():
    # Deliberately over-large proposed square makes negative residual diagonal.
    with pytest.raises(ValueError,match='diagonally dominant'):
        verify_spectral_matrix(A,[[100,0],[0,100]],1,2,2)