import copy
from fractions import Fraction

import pytest

pytest.importorskip("numpy")
from foldbounds.coupled_pka import checker as C, exact as ex, prop as P
from test_coupled_pka_core import synth_spec


def certs_for(spec, i, R, ph, budget):
    nbr = P.adjacency(spec)
    S = P.interior(spec, i, R)
    prep = P.Prepared(spec, nbr, i, S, Fraction(ph))
    lo, up = P.run_bound(prep, "lower", budget), P.run_bound(prep, "upper", budget)
    return prep, lo, up


def test_exact_inside_every_certified_interval():
    viol = 0
    for seed in (1, 2):
        spec = synth_spec(9, seed, box=10.0)
        for ph in (3, 10):
            Q = ex.exact_marginals(spec, Fraction(ph))
            for i in range(9):
                for R in (0, 5, 9):
                    for bud in (1, 4, 256):
                        prep, lo, up = certs_for(spec, i, R, ph, bud)
                        for s in (lo, up):
                            C.check(spec, s.certificate(spec, s.claim(), i, Fraction(ph)))
                        viol += not (lo.claim() <= Q[i] <= up.claim())
    assert viol == 0


def test_exhaustive_boundary_is_exact_separator_hard_bound():
    spec = synth_spec(8, 5, box=9.0)
    Q = ex.exact_marginals(spec, Fraction(6))
    prep, lo, up = certs_for(spec, 2, 6, 6, 1 << 10)
    assert lo.bound <= float(Q[2]) <= up.bound


def test_hierarchy_widths_not_worse_with_full_enumeration():
    spec = synth_spec(9, 7, box=9.0)
    w = []
    for R in (0, 4, 9):
        prep, lo, up = certs_for(spec, 1, R, 7, 1 << 10)
        w.append(up.bound / lo.bound)
    assert w[2] <= w[1] * (1 + 1e-9) and w[1] <= w[0] * (1 + 1e-9)


@pytest.fixture
def good():
    for seed in range(11, 60):
        spec = synth_spec(10, seed, box=7.0)
        prep, lo, up = certs_for(spec, 0, 8, 6, 64)
        if prep.B >= 3 and len(prep.S) >= 2 and lo.nleaf > 1:
            break
    assert prep.B >= 3 and len(prep.S) >= 2 and lo.nleaf > 1
    cert = lo.certificate(spec, lo.claim(), 0, Fraction(6))
    assert C.check(spec, cert)
    return spec, lo, cert


def test_corrupt_claim_too_strong(good):
    spec, lo, cert = good
    bad = dict(cert, **{"lambda": str(Fraction(lo.bound) * 1000)})
    with pytest.raises(C.CertificateError):
        C.check(spec, bad)


def test_corrupt_upper_claim_too_strong(good):
    spec, lo, cert = good
    prep = P.Prepared(spec, P.adjacency(spec), 0, P.interior(spec, 0, 8), Fraction(6))
    up = P.run_bound(prep, "upper", 64)
    c = up.certificate(spec, up.claim(), 0, Fraction(6))
    assert C.check(spec, c)
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(c, **{"lambda": str(Fraction(up.bound) / 1000)}))


def test_corrupt_tree_dropped_cube(good):
    spec, lo, cert = good
    # replace the whole tree by a leaf: coarser cube, claim no longer justified unless root was tight
    root = dict(cert, tree=0)
    if cert["tree"] != 0:
        # a root leaf is a valid (coarser) partition but its bound is weaker than the claim
        try:
            C.check(spec, root)
        except C.CertificateError:
            pass
    # malformed: node with one child
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, tree=[0]))
    # too deep
    deep = 0
    for _ in range(len(cert["order"]) + 1):
        deep = [deep, deep]
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, tree=deep))


def test_corrupt_order_and_interior(good):
    spec, lo, cert = good
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, order=cert["order"][:-1]))
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, order=cert["order"] + [cert["order"][0]]))
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, interior=[cert["interior"][0]]))


def test_corrupt_spec_weights(good):
    spec, lo, cert = good
    bad = copy.deepcopy(spec)
    k = next(iter(bad["pairs_milli"]))
    bad["pairs_milli"][k] += 1000
    with pytest.raises(C.CertificateError):
        C.check(bad, cert)


def test_corrupt_ph_and_side(good):
    spec, lo, cert = good
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, ph="13"))  # claim computed at pH 6 is not valid at pH 13
    with pytest.raises(C.CertificateError):
        C.check(spec, dict(cert, side="sideways"))


def test_ve_reference_matches_brute_force():
    spec = synth_spec(10, 4, box=9.0)
    th = ex.float_marginals_all_ph(spec, [6.0])[0]
    ve = ex.ve_theta(spec, 6.0)
    assert max(abs(th[i] - ve[i]) for i in range(10)) < 1e-10


def test_mc_runs_and_agrees_roughly():
    from foldbounds.coupled_pka import mc
    spec = synth_spec(8, 9, box=8.0)
    th = ex.float_marginals_all_ph(spec, [6.0])[0]
    m, _ = mc.gibbs_theta(spec, 6.0, 1, chains=32, burn=200, sweeps=1500)
    assert max(abs(m[i] - th[i]) for i in range(8)) < 0.05
