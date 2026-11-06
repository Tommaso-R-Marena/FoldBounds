"""Independent exact-rational tests for primal-dual variance certificates."""
from fractions import Fraction as Q
from itertools import product
from random import Random
import pytest
from foldbounds.variance_dual import (
    VarianceDualCertificate, maximize_relaxed_variance, verify_variance_dual,
)


def brute_exact(known, caps):
    """Independent exhaustive extreme models + all pairwise normalized mixes."""
    if not known:
        scores = [lo for lo, _, _ in caps]+[hi for _, hi, _ in caps]
        return Q(max(scores)-min(scores))**2/4
    models = []
    for states in product(range(3), repeat=len(caps)):
        atoms = list(known)
        for state, (lo, hi, cap) in zip(states, caps):
            if state:
                atoms.append((lo if state == 1 else hi, cap))
        z = sum((w for _, w in atoms), Q(0))
        mu = sum((k*w for k, w in atoms), Q(0))/z
        var = sum((w*(k-mu)**2 for k, w in atoms), Q(0))/z
        models.append((mu, var))
    best = Q(0)
    for a, va in models:
        for b, vb in models:
            gap = (a-b)**2
            alpha = max(Q(0), min(Q(1), Q(1, 2)+(va-vb)/(2*gap))) if gap else Q(0)
            best = max(best, alpha*va+(1-alpha)*vb+alpha*(1-alpha)*gap)
    return best


def test_exact_hand_examples():
    examples = (
        (((1, Q(1)),), ()),
        (((0, Q(1)), (3, Q(2))), ()),
        ((), ((0, 3, Q(1)),)),
        ((), ((0, 0, Q(2)),)),
        ((), ((0, 1, Q(1)), (4, 7, Q(2)))),
        (((0, Q(1)),), ((0, 5, Q(1)),)),
        (((1, Q(3)), (7, Q(1))), ((2, 9, Q(7)),)),
    )
    for known, caps in examples:
        proof = maximize_relaxed_variance(known, caps)
        assert verify_variance_dual(known, caps, proof)
        assert proof.upper == brute_exact(known, caps)


def test_random_exact_vs_independent_exhaustive_mixtures():
    rng = Random(20261105)
    for _ in range(160):
        kcount, ccount = rng.randrange(4), rng.randrange(4)
        if not kcount and not ccount:
            ccount = 1
        known = tuple((rng.randrange(7), Q(rng.randrange(1,5),rng.randrange(1,4)))
                      for _ in range(kcount))
        caps = []
        for _ in range(ccount):
            low = rng.randrange(6)
            caps.append((low, low+rng.randrange(4),
                         Q(rng.randrange(1,8),rng.randrange(1,5))))
        caps = tuple(caps)
        proof = maximize_relaxed_variance(known, caps)
        assert verify_variance_dual(known, caps, proof)
        assert proof.upper == brute_exact(known, caps), (known, caps, proof)


def test_rejects_forged_dual_and_invalid_primal():
    known, caps = ((1,Q(1)),), ((0,4,Q(2)),)
    good = maximize_relaxed_variance(known, caps)
    with pytest.raises(ValueError, match="dual"):
        verify_variance_dual(known, caps,
                             VarianceDualCertificate(Q(0), good.center,
                                                     good.allocations, 1))
    with pytest.raises(ValueError, match="primal"):
        verify_variance_dual(known, caps,
                             VarianceDualCertificate(good.upper, good.center,
                                                     ((),), 1))
    with pytest.raises(ValueError, match="exceeds"):
        verify_variance_dual(known, caps,
                             VarianceDualCertificate(good.upper, good.center,
                                                     (((0,Q(10)),),), 1))


def test_strictly_below_range_only_bound():
    known, caps = ((2,Q(99)),), ((0,8,Q(1)),)
    result = maximize_relaxed_variance(known, caps)
    assert Q(0) < result.upper < Q(16)/3


def test_fail_closed_invalid_inputs_and_budget():
    with pytest.raises(ValueError):
        maximize_relaxed_variance(((1,Q(-1)),), ())
    with pytest.raises(ValueError):
        maximize_relaxed_variance((), ((2,1,Q(1)),))
    with pytest.raises(ValueError):
        maximize_relaxed_variance(((2,Q(9)),), ((0,8,Q(1)),),
                                  max_iterations=0)
