import itertools
import random
from fractions import Fraction as F

import pytest

from foldbounds.coupled_mrf import baselines, constructions as cons
from foldbounds.coupled_mrf.certificate import (
    all_local_ratios, annulus_contraction, check_nesting, exhaustive_enclosure,
    local_ratio, shell_hierarchy)
from foldbounds.coupled_mrf.generators import (
    ising_grid, random_sparse_binary, rotamer_potts, rotamer_grid_potts, ising_q)
from foldbounds.coupled_mrf.io import loads_model, model_to_dict, model_sha256
from foldbounds.coupled_mrf.model import Factor, ModelError, as_rational
from foldbounds.coupled_mrf.oracles import brute_force_z, grid_transfer_z, ve_z
from foldbounds.coupled_mrf.separator import (
    Separator, SeparatorError, radius_separator, require_separator, separator_violations)
import json


def exact_ratio(m):
    return ve_z(m, "B") / ve_z(m, "A")


def test_floats_rejected():
    with pytest.raises(ModelError):
        as_rational(0.5)
    with pytest.raises(ModelError):
        as_rational("1e3")
    with pytest.raises(ModelError):
        loads_model('{"format":"foldbounds.coupled_mrf","version":1,"variables":[{"name":"a","domain":2}],'
                    '"shared":[{"scope":["a"],"table":[1,0.5]}]}')


def test_model_json_roundtrip_hash():
    m = ising_grid(3, F(3, 2), F(6, 5), ("field", (1, 1), F(1, 2)))
    m2 = loads_model(json.dumps(model_to_dict(m)))
    assert model_sha256(m) == model_sha256(m2)
    assert brute_force_z(m2, "B") == brute_force_z(m, "B")


@pytest.mark.parametrize("pert", [("field", (2, 1), F(1, 3)), ("edge", ((1, 1), (1, 2)), F(5, 2)), None])
def test_oracles_agree_grid(pert):
    m = ising_grid(4, F(7, 5), F(5, 4), pert)
    for w in "AB":
        z = brute_force_z(m, w)
        assert ve_z(m, w) == z == grid_transfer_z(m, w)


def test_oracles_agree_random_and_rect():
    m = random_sparse_binary(12, 3, seed=3)
    for w in "AB":
        assert brute_force_z(m, w) == ve_z(m, w)
    m = ising_grid(3, F(2), F(1), ("field", (1, 2), F(2)), M=5)
    assert grid_transfer_z(m, "B") == brute_force_z(m, "B")


def test_oracle_known_value():
    # one edge, q=2: Z = 2+1+1+2 = 6; field u on one site: [1,u] -> 2 + 1 + u + 2u
    m = ising_grid(2, F(2), M=1, perturbation=("field", (0, 0), F(3)))
    assert brute_force_z(m, "A") == 6 and brute_force_z(m, "B") == 2 + 1 + 3 + 6


def test_separator_checks():
    m = ising_grid(5, F(3, 2), F(1), ("field", (2, 2), F(2)))
    sep = radius_separator(m, 2)
    require_separator(m, sep)
    # moving a boundary vertex to exterior breaks (S1)
    bad = Separator.of(sep.interior, sep.boundary[1:], set(sep.exterior) | {sep.boundary[0]})
    assert separator_violations(m, bad)
    with pytest.raises(SeparatorError):
        require_separator(m, bad)
    # perturbed variable outside interior+boundary violates (S2)
    far = radius_separator(m, 1, core=[(f"x0_0")])
    assert any("system-specific" in p for p in separator_violations(m, far))
    # not a partition
    assert separator_violations(m, Separator.of(sep.interior, sep.boundary, []))


@pytest.mark.parametrize("seed", range(6))
def test_enclosure_contains_truth_random(seed):
    m = random_sparse_binary(12, 3, seed=seed)
    exact = exact_ratio(m)
    levels = shell_hierarchy(m, [0, 1, 2, 3, 4])
    assert not check_nesting(levels)
    for lev in levels:
        if lev.enclosure is not None:
            assert lev.enclosure.contains(exact), lev.radius
    assert levels[-1].enclosure.lower == levels[-1].enclosure.upper == exact or levels[-1].enclosure.n_states > 1


def test_finite_termination():
    m = random_sparse_binary(10, 3, seed=11)
    exact = exact_ratio(m)
    levels = shell_hierarchy(m, list(range(0, 12)))
    last = [l for l in levels if l.enclosure is not None and l.enclosure.n_states == 1]
    # once the boundary is empty the enclosure is the exact value (Thm 2.5)
    empties = [l for l in levels if l.enclosure is not None and l.enclosure.n_states == 1 and not l.separator.boundary]
    assert empties and all(l.enclosure.lower == l.enclosure.upper == exact for l in empties)


def test_ising_nested_and_contains():
    m = ising_grid(5, F(3, 2), F(6, 5), ("field", (2, 2), F(1, 2)))
    exact = grid_transfer_z(m, "B") / grid_transfer_z(m, "A")
    levels = shell_hierarchy(m, [0, 1, 2, 3])
    assert not check_nesting(levels)
    widths = [l.enclosure.width() for l in levels]
    assert all(l.enclosure.contains(exact) for l in levels)
    assert widths == sorted(widths, reverse=True)


def test_rotamer_mutation_with_domain_change():
    m = rotamer_potts(7, seed=2, rotamers=(2, 3), mutant_rotamers=4, mutated=3)
    exact = exact_ratio(m)
    with pytest.raises(SeparatorError):   # R=0: changed-domain variable on the boundary
        require_separator(m, radius_separator(m, 0))
    for R in (1, 2, 3):
        enc = exhaustive_enclosure(m, radius_separator(m, R))
        assert enc.contains(exact)


def test_local_ratio_matches_tables():
    m = ising_grid(4, F(3, 2), F(1), ("field", (1, 1), F(2)))
    sep = radius_separator(m, 1)
    pa, pb = all_local_ratios(m, sep)
    k = 5
    from foldbounds.coupled_mrf.certificate import decode_state
    st = decode_state(k, [2] * len(sep.boundary))
    assert local_ratio(m, sep, st) == (pa[k], pb[k])


def test_unbounded_when_phi_a_zero():
    # B has a field factor that kills nothing, A has a hard zero on a boundary state
    names = ["i", "b"]
    dom = {"i": 2, "b": 2}
    shared = [Factor(("i", "b"), (F(1), F(0), F(1), F(1)))]       # (i=0,b=1) forbidden
    onlya = [Factor(("i",), (F(1), F(0)))]                        # A forces i=0 => Phi_A(b=1)=0
    from foldbounds.coupled_mrf.model import make_model
    m = make_model(dom, shared=shared, only_a=onlya)
    sep = Separator.of(["i"], ["b"], [])
    enc = exhaustive_enclosure(m, sep)
    assert enc.unbounded and enc.upper is None and enc.unbounded_witness == (1,)


def test_exhaustive_large_boundary_chunked():
    m = ising_grid(7, F(3, 2), F(1), ("field", (3, 3), F(2)))
    sep = radius_separator(m, 3)
    e1 = exhaustive_enclosure(m, sep, chunk_states=2 ** 20)
    e2 = exhaustive_enclosure(m, sep, chunk_states=64)
    assert (e1.lower, e1.upper) == (e2.lower, e2.upper)


def test_contraction_diagnostic():
    m = ising_grid(6, F(3, 2), F(6, 5), ("field", (3, 3), F(1, 2)))
    d = annulus_contraction(m, radius_separator(m, 1), radius_separator(m, 2))
    assert d.identity_ok and d.bound_ok and d.exact
    assert d.dobrushin <= 1
    d2 = annulus_contraction(m, radius_separator(m, 2), radius_separator(m, 3))
    assert d2.identity_ok and d2.bound_ok


# ---------------------------------------------------------------- baselines
def test_tree_bounds_exact_on_trees_and_bracket():
    path = ising_grid(5, F(3, 2), F(1), ("field", (2, 0), F(2)), M=1)
    for w in "AB":
        lb, ub = baselines.tree_bounds(path, w)
        assert lb == ub == brute_force_z(path, w)
    m = ising_grid(4, F(3, 2), F(6, 5), ("field", (2, 2), F(1, 2)))
    for w in "AB":
        lb, ub = baselines.tree_bounds(m, w)
        assert 0 < lb <= brute_force_z(m, w) <= ub
    lo, hi, _ = baselines.separate_ratio_interval(m)
    assert lo <= exact_ratio(m) <= hi


def test_coupled_narrower_than_separate():
    m = ising_grid(6, F(3, 2), F(1), ("field", (3, 3), F(1, 2)))
    lo, hi, _ = baselines.separate_ratio_interval(m)
    enc = exhaustive_enclosure(m, radius_separator(m, 2))
    assert enc.width() < hi - lo


# ------------------------------------------------------------ sharpness
@pytest.mark.parametrize("theta", [F(0), F(1), F(1, 3)])
def test_sharp_exterior_realises_endpoints(theta):
    base = ising_grid(2, F(3, 2), F(6, 5), ("field", (0, 0), F(1, 2)))
    # keep centre (0,0) interior and its two grid neighbours as the boundary
    from foldbounds.coupled_mrf.model import CoupledModel
    keep = ["x0_0", "x0_1", "x1_0"]
    shared = tuple(f for f in base.shared if set(f.scope) <= set(keep))
    only_b = tuple(f for f in base.only_b if set(f.scope) <= set(keep))
    small = CoupledModel(tuple(keep), {v: 2 for v in keep}, {v: 2 for v in keep}, shared, (), only_b, {})
    sep = Separator.of(["x0_0"], ["x0_1", "x1_0"], [])
    enc = exhaustive_enclosure(small, sep)
    model, sep2, lam = cons.sharp_exterior_model(small, sep, theta)
    require_separator(model, sep2)
    assert exact_ratio(model) == lam == brute_force_z(model, "B") / brute_force_z(model, "A")
    if theta == 1:
        assert lam == enc.lower
    if theta == 0:
        assert lam == enc.upper


# ---------------------------------------------------------- soft separator
def test_soft_enclosure_sound_on_ising_and_between_hard_and_exact():
    q = F(3, 2)
    m = ising_grid(5, q, F(1), ("field", (2, 2), F(1, 2)))
    exact = grid_transfer_z(m, "B") / grid_transfer_z(m, "A")
    for R in (1, 2):
        sep = radius_separator(m, R)
        lo, hi = cons.soft_enclosure(m, sep)
        hard = exhaustive_enclosure(m, sep)
        assert lo <= exact <= hi
        assert hard.lower <= lo and hi <= hard.upper


def test_soft_monotone_in_c_and_s0_exact():
    m = ising_grid(5, F(3, 2), F(1), ("field", (2, 2), F(1, 2)))
    sep = radius_separator(m, 1)
    pa, pb = all_local_ratios(m, sep)
    dims = [2] * 4
    w = []
    for c in (F(1), F(11, 10), F(3, 2), F(3), F(50)):
        lo, hi, _, _ = cons.soft_enclosure_from_tables(pa, pb, dims, [c] * 4)
        w.append(hi - lo)
    assert w[0] == 0 and w == sorted(w)
    assert cons.soft_enclosure_from_tables(pa, pb, dims, [F(1)] * 4)[0] == sum(pb) / sum(pa)


def test_soft_hypothesis_violation_detected_and_counterexample():
    m = ising_grid(5, F(3, 2), F(1), ("field", (2, 2), F(1, 2)))
    # a factor on two boundary variables placed in the exterior breaks E_s
    sep = radius_separator(m, 1)
    b = sep.boundary
    extra = Factor((b[0], b[1], "x0_0"), tuple([F(1)] * 8))
    from foldbounds.coupled_mrf.model import CoupledModel
    m2 = CoupledModel(m.variables, m.domains_a, m.domains_b, m.shared + (extra,), m.only_a, m.only_b, {})
    sep2 = Separator.of(sep.interior | {"x0_0"} if False else sep.interior, sep.boundary,
                        sep.exterior)
    # x0_0 is exterior; factor touches two boundary vars and an exterior var
    with pytest.raises(cons.SoftHypothesisError):
        cons.soft_sensitivities(m2, sep2)
    r = cons.two_boundary_counterexample(F(2))
    assert r["violates"]
    assert r["ratio_agree_exterior"] > r["soft_hi"] and r["ratio_disagree_exterior"] < r["soft_lo"]


def test_soft_sound_random_product_exteriors():
    rng = random.Random(5)
    for _ in range(40):
        n = 3
        pa = [F(rng.randint(1, 9)) for _ in range(2 ** n)]
        pb = [F(rng.randint(1, 9)) for _ in range(2 ** n)]
        c = [F(rng.randint(11, 40), 10) for _ in range(n)]
        lo, hi, _, _ = cons.soft_enclosure_from_tables(pa, pb, [2] * n, c)
        for _ in range(30):
            # mixture of exterior states, each a product of fields with ratio <= c_v
            tot_a = tot_b = F(0)
            for _o in range(rng.randint(1, 4)):
                g = F(rng.randint(1, 5))
                f = [[F(1), 1 + (c[v] - 1) * F(rng.randint(0, 10), 10)] for v in range(n)]
                for v in range(n):
                    if rng.random() < .5:
                        f[v].reverse()
                za = zb = F(0)
                for k, t in enumerate(itertools.product((0, 1), repeat=n)):
                    p = F(1)
                    for v in range(n):
                        p *= f[v][t[v]]
                    za += pa[k] * p
                    zb += pb[k] * p
                tot_a += g * za
                tot_b += g * zb
            assert lo <= tot_b / tot_a <= hi
