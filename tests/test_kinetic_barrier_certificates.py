"""Independent reproducibility, adversarial replay, and dynamics checks."""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.kinetic_barriers import (
    BasinCertificate, ProposalLimit, _local_moves, _walks, propose_kinetic_basin,
)
from foldbounds.kinetic_barrier_verifier import (
    KineticCertificateError, _moves_independent, _score_independent,
    _walks_independent, verify_kinetic_basin,
)


@pytest.mark.parametrize("n,count", [(4, 9), (5, 25), (6, 71), (7, 195),
                                     (8, 543), (9, 1479), (10, 4067),
                                     (11, 11025), (12, 30073)])
def test_two_full_state_enumerators_agree(n, count):
    assert set(_walks(n, 40000)) == set(_walks_independent(n, 40000))
    assert len(_walks(n, 40000)) == count


@pytest.mark.parametrize("n", range(4, 9))
def test_two_move_generators_and_inverse_detailed_balance(n):
    all_paths = _walks(n, 40000)
    universe = set(all_paths)
    for p in all_paths:
        produced = set(_local_moves(p))
        independent = _moves_independent(p)
        assert produced == independent
        assert produced <= universe
        for target in independent:
            assert p in _moves_independent(target)
            for sequence in ("H" * n, "HP" * (n // 2) + ("H" if n % 2 else "")):
                kp = _score_independent(p, sequence)
                kt = _score_independent(target, sequence)
                for q in (Fraction(1), Fraction(3, 2), Fraction(10)):
                    forward = q**kp * min(Fraction(1), q**(kt-kp))
                    backward = q**kt * min(Fraction(1), q**(kp-kt))
                    assert forward == backward


def test_strong_all_q_cubic_kinetic_basin_9mer():
    q = Fraction(10000)
    cert = propose_kinetic_basin("H"*9, q, min_barrier=3)
    r = verify_kinetic_basin(cert, q, 10**9)
    assert r.total_states == 1479
    assert r.native_score == 4
    assert r.native_states == 10
    assert r.native_equilibrium_probability > Fraction(998, 1000)
    assert r.basin_states == 5 and r.basin_top_score == 3
    assert r.boundary_edges == 2 and r.boundary_top_score == 0
    assert r.barrier_order == 3
    assert r.stationary_exit_per_step == Fraction(1, 14002800280000)
    assert r.all_q_coefficient == Fraction(1, 14)
    assert r.exit_by_steps_upper < Fraction(1, 10000)
    assert not r.closed_class
    for q_alt in (Fraction(1), Fraction(3, 2), Fraction(2), Fraction(17), Fraction(1000)):
        result = verify_kinetic_basin(cert, q_alt, 23)
        assert result.stationary_exit_per_step <= result.all_q_exit_per_step_bound
        assert result.barrier_order == 3


def test_milliard_update_ultraslow_basin_11mer():
    q = Fraction(10000)
    cert = propose_kinetic_basin("H"*11, q, min_barrier=4)
    r = verify_kinetic_basin(cert, q, 10**9)
    assert r.total_states == 11025
    assert r.native_score == 5 and r.native_states == 96
    assert r.native_equilibrium_probability > Fraction(999, 1000)
    assert r.basin_states == 321
    assert r.basin_top_score == 4 and r.boundary_top_score == 0
    assert r.boundary_edges == 195 and r.barrier_order == 4
    assert r.stationary_exit_per_step == Fraction(13, 696103226403336000)
    assert r.exit_by_steps_upper < Fraction(1, 10**7)
    assert not r.closed_class


def test_frozen_misfold_is_permanent_under_declared_move_set_12mer():
    q = Fraction(10000)
    cert = propose_kinetic_basin("H"*12, q, allow_closed=True, min_barrier=4)
    r = verify_kinetic_basin(cert, q, 10**15)
    assert r.total_states == 30073
    assert r.native_score == 6 and r.native_states == 62
    assert r.native_equilibrium_probability > Fraction(998, 1000)
    assert r.basin_top_score == 5 and r.basin_states == 1
    assert r.boundary_edges == 0 and r.closed_class
    assert r.exit_by_steps_upper == 0
    assert _moves_independent(cert.seed) == set()
    # The frozen state is only frozen under THIS local move set, not all
    # conceivable polymer rearrangements. A nonlocal pivot allows escape.
    pivot_point = cert.seed[1]
    def clockwise(z):
        x, y = z[0]-pivot_point[0], z[1]-pivot_point[1]
        return (pivot_point[0]+y, pivot_point[1]-x)
    pivot = cert.seed[:2] + tuple(clockwise(z) for z in cert.seed[2:])
    assert len(set(pivot)) == 12
    assert pivot != cert.seed


def test_exact_one_step_exit_identity_and_finite_time_union_bound():
    q = Fraction(10)
    cert = propose_kinetic_basin("H"*9, q, min_barrier=3)
    r = verify_kinetic_basin(cert, q, 3)
    S = set(cert.members)
    scores = {p: _score_independent(p, cert.sequence) for p in S}
    dist = {p: q**score / r.basin_mass for p, score in scores.items()}
    exited = Fraction(0)
    for timestep in range(1, 4):
        nxt = {}
        for p, mass in dist.items():
            move_mass = Fraction(0)
            for v in _moves_independent(p):
                kv = _score_independent(v, cert.sequence)
                rate = min(Fraction(1), q**(kv - scores[p])) / r.proposal_slots
                move_mass += rate
                if v in S:
                    nxt[v] = nxt.get(v, Fraction(0)) + mass*rate
                else:
                    exited += mass*rate
            nxt[p] = nxt.get(p, Fraction(0)) + mass*(1-move_mass)
        dist = nxt
        if timestep == 1:
            assert exited == r.stationary_exit_per_step
        assert exited <= timestep * r.stationary_exit_per_step
    assert exited <= r.exit_by_steps_upper


def test_missing_member_and_fake_member_fail_closed():
    cert = propose_kinetic_basin("H"*9, Fraction(100), min_barrier=3)
    examples = [
        replace(cert, members=cert.members[1:]),
        replace(cert, members=cert.members + (cert.members[0],)),
        replace(cert, members=cert.members + (((0,0), (1,0)),)),
        replace(cert, seed=((0,0),(1,0))),
        replace(cert, threshold=4),
        replace(cert, members=()),
    ]
    for candidate in examples:
        with pytest.raises(KineticCertificateError):
            verify_kinetic_basin(candidate, Fraction(100), 100)


def test_adversarial_member_from_other_basin_rejected():
    cert = propose_kinetic_basin("H"*9, Fraction(10), min_barrier=3)
    other = next(p for p in _walks(9, 40000)
                 if p not in cert.members and _score_independent(p, "H"*9) >= cert.threshold)
    with pytest.raises(KineticCertificateError):
        verify_kinetic_basin(replace(cert, members=cert.members + (other,)), Fraction(10), 2)


def test_validation_and_budget_fail_closed():
    with pytest.raises(ProposalLimit):
        propose_kinetic_basin("H"*12, Fraction(3), max_states=50)
    cert = propose_kinetic_basin("H"*9, Fraction(3), min_barrier=3)
    with pytest.raises(KineticCertificateError):
        verify_kinetic_basin(cert, Fraction(3), 1, max_states=50)
    for q in (1, 1.0, Fraction(0), Fraction(1, 2)):
        with pytest.raises((ValueError, KineticCertificateError)):
            verify_kinetic_basin(cert, q, 2)
    for steps in (-1, 1.5, True):
        with pytest.raises(KineticCertificateError):
            verify_kinetic_basin(cert, Fraction(2), steps)
    for seq in ("", "HHX", "H"):
        with pytest.raises(ValueError):
            propose_kinetic_basin(seq, Fraction(2))
    with pytest.raises(ProposalLimit):
        propose_kinetic_basin("PPPPPPPP", Fraction(2))


def test_both_sides_symmetric_and_all_hp_compositions_short():
    # Independent geometry and scores across all sequences for n <= 6.
    for n in range(4, 7):
        paths = _walks_independent(n, 40000)
        for seq in map("".join, product("HP", repeat=n)):
            brute = []
            for p in paths:
                a = sum(seq[i] == seq[j] == "H" and
                        abs(p[i][0]-p[j][0])+abs(p[i][1]-p[j][1]) == 1
                        for i in range(n) for j in range(i+2,n))
                brute.append(a)
                assert a == _score_independent(p, seq)
            assert min(brute) == 0


def test_exact_min_cut_13x_sharper_and_fixed_start_11mer():
    from foldbounds.kinetic_barriers import propose_optimal_kinetic_cut
    q = Fraction(10000)
    level = propose_kinetic_basin("H"*11, q, min_barrier=4)
    baseline = verify_kinetic_basin(level, q, 10**9)
    cert = propose_optimal_kinetic_cut("H"*11, q)
    r = verify_kinetic_basin(cert, q, 10**9)
    assert r.optimal_cut_verified_at_q
    assert r.basin_states == 445
    assert r.seed_score == 4
    assert r.boundary_edges == 15
    assert r.boundary_flux == 15
    assert r.barrier_order == 4
    assert baseline.boundary_flux / r.boundary_flux == 13
    assert r.deterministic_seed_exit_by_steps_upper == Fraction(1, 24000000)
    assert r.deterministic_seed_exit_by_steps_upper < Fraction(1, 1000000)
    assert r.exit_by_steps_upper < baseline.exit_by_steps_upper
    assert r.native_equilibrium_probability == baseline.native_equilibrium_probability


def test_rational_maxflow_min_cut_duality_and_tamper_rejection():
    from foldbounds.kinetic_barriers import propose_optimal_kinetic_cut
    q = Fraction(3, 2)
    cert = propose_optimal_kinetic_cut("H"*9, q)
    r = verify_kinetic_basin(cert, q, 8)
    assert r.optimal_cut_verified_at_q
    assert r.boundary_flux == Fraction(3, 2)
    assert r.deterministic_seed_exit_by_steps_upper <= 1
    # Changing q invalidates the flow proof but the same geometry can still
    # be independently checked as a non-optimal arbitrary cut at any q.
    with pytest.raises(KineticCertificateError):
        verify_kinetic_basin(cert, Fraction(2), 8)
    ordinary = replace(cert, mode="arbitrary", flow_q=None, flow_edges=())
    assert not verify_kinetic_basin(ordinary, Fraction(2), 8).optimal_cut_verified_at_q
    a, b, f = cert.flow_edges[0]
    corruptions = (
        replace(cert, flow_edges=cert.flow_edges[1:]),
        replace(cert, flow_edges=((a, b, f+1),) + cert.flow_edges[1:]),
        replace(cert, flow_edges=((a, b, 10**99),) + cert.flow_edges[1:]),
        replace(cert, flow_edges=((a, a, 1),) + cert.flow_edges[1:]),
        replace(cert, flow_edges=(cert.flow_edges[0],) + cert.flow_edges),
        replace(cert, mode="unknown"),
    )
    for forged in corruptions:
        with pytest.raises(KineticCertificateError):
            verify_kinetic_basin(forged, q, 8)


def test_every_edge_of_min_cut_verified_without_proposal_algorithms():
    from foldbounds.kinetic_barriers import propose_optimal_kinetic_cut
    for q in (Fraction(1), Fraction(3, 2), Fraction(2), Fraction(10), Fraction(10000)):
        c = propose_optimal_kinetic_cut("H"*9, q)
        r = verify_kinetic_basin(c, q, 250)
        assert r.optimal_cut_verified_at_q
        assert r.native_equilibrium_probability >= 0
        assert r.stationary_exit_per_step <= r.all_q_exit_per_step_bound
        assert r.deterministic_seed_exit_by_steps_upper <= min(
            Fraction(1), 250 * r.deterministic_all_q_exit_per_step_bound)
