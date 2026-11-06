import itertools
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import pytest
from gnm_lowrank_witnesses import grounded_incidence,woodbury_inverse,integer_proposal
import gnm_exact_certificates as gc


def test_all_boolean_contact_subsets_exact_residual_checks():
    base=[(0,1),(1,2),(2,3),(3,4)]
    possible=[(0,2),(0,3),(1,3),(1,4)]
    n=5;root=4;Q=1<<50
    A0,idx=gc.grounded_laplacian(base,n,root)
    X=np.linalg.inv(np.asarray(A0,float))
    for mask in range(1<<len(possible)):
        chosen=[possible[j] for j in range(len(possible)) if mask>>j&1]
        A,_=gc.grounded_laplacian(base+chosen,n,root)
        candidate=woodbury_inverse(X,root,n,chosen)
        assert np.allclose(np.asarray(A,float)@candidate,np.eye(n-1),rtol=0,atol=1e-12)
        w=integer_proposal(X,root,n,chosen,Q)
        r=gc.verify_integer_witness(A,w,Q)
        # Guaranteed exact bound on the true resistance of the pair (0,3).
        lo,hi=gc.score_numerators((0,3),w,r,Q,20,idx)
        vec=np.zeros(n-1);vec[idx[0]]+=1;vec[idx[3]]-=1
        actual=float(vec@np.linalg.solve(np.asarray(A,float),vec))
        assert (lo/Q**2)-1e-12<=actual<=(hi/Q**2)+1e-12


def test_invalid_witness_inputs():
    with pytest.raises(ValueError):grounded_incidence(6,5,[])
    with pytest.raises(ValueError):grounded_incidence(0,3,[(1,1)])
    with pytest.raises(ValueError):woodbury_inverse(np.zeros((2,3)),0,3,[])
    with pytest.raises(ValueError):integer_proposal(np.eye(2),0,3,[],0)