"""JEPA is an untrusted scheduling policy, NEVER a source of bounds."""
from dataclasses import replace
from fractions import Fraction
import math

import pytest

from foldbounds.adaptive import verify_refinement_run
from foldbounds.jepa_guidance import (
    ACTION_DIM, STATE_DIM, action_features, build_model, candidate_paths,
    collect_checked_decisions, fit_checked_jepa, guided_refine, state_features,
    train_online_jepa_rl, verify_guided_run,
)
from foldbounds.lattice import build_frontier, exact_spectrum


def _oracle(sequence, threshold=1):
    q = Fraction(2)
    spec = exact_spectrum(sequence)
    z = sum((c * q**k for k, c in spec.items()), Fraction())
    e = sum((c * q**k for k, c in spec.items() if k >= threshold), Fraction())
    return z, e / z


@pytest.mark.parametrize("sequence", ["HHHHHH", "HPHHPPHH"])
def test_every_learned_choice_stays_sound_and_replayable(sequence):
    cert = build_frontier(sequence, Fraction(2), 1, 2)
    exact_z, exact_p = _oracle(sequence)

    def nonsense(state, actions):
        return tuple((len(actions) - i) * (-1.0 if state[0] > 0.2 else 1.0)
                     for i in range(len(actions)))

    run = guided_refine(cert, steps=7, scorer=nonsense)
    assert verify_guided_run(run) == run.refinement.bounds
    assert verify_refinement_run(run.refinement) == run.refinement.bounds
    for snapshot in run.refinement.history:
        assert snapshot.partition_lower <= exact_z <= snapshot.partition_upper
        assert snapshot.probability_lower <= exact_p <= snapshot.probability_upper
    for before, after in zip(run.refinement.history, run.refinement.history[1:]):
        assert after.probability_width <= before.probability_width
    assert run.checker_calls >= len(run.refinement.history)


def test_features_dimensions_and_closed_action_set():
    cert = build_frontier("HHHHHH", Fraction(2), 1, 3)
    run = guided_refine(cert, steps=0)
    state = state_features(cert, run.refinement.bounds)
    assert len(state) == STATE_DIM and all(math.isfinite(x) for x in state)
    for path in candidate_paths(cert):
        action = action_features(cert, path)
        assert len(action) == ACTION_DIM and all(math.isfinite(x) for x in action)


def test_malformed_scorer_abstains_before_unsafe_transition():
    cert = build_frontier("HHHHHH", Fraction(2), 1, 3)
    with pytest.raises(ValueError):
        guided_refine(cert, steps=1, scorer=lambda s, a: [math.nan] * len(a))
    with pytest.raises(ValueError):
        guided_refine(cert, steps=1, scorer=lambda s, a: ())


def test_checked_oracle_counterfactuals_and_determinism():
    opts = dict(steps=3, seed=42, lookahead=0, use_matching=True)
    a = collect_checked_decisions(("HHHHHH", "HPHHPPHH"), **opts)
    assert a == collect_checked_decisions(("HHHHHH", "HPHHPPHH"), **opts)
    assert a
    assert any(len(x.candidates) > 1 for x in a)
    for decision in a:
        assert 0 <= decision.chosen < len(decision.candidates)
        for candidate in decision.candidates:
            assert candidate.gain >= 0 and candidate.partition_gain >= 0
            assert candidate.verifier_calls == 1
    cert = build_frontier("HHHHHH", Fraction(2), 1, 2)
    audited = guided_refine(cert, steps=3, collect_counterfactuals=True)
    assert verify_guided_run(audited) == audited.refinement.bounds
    assert audited.checker_calls > len(audited.refinement.history)


def test_tampered_neural_run_is_rejected():
    cert = build_frontier("HHHHHH", Fraction(2), 1, 2)
    run = guided_refine(cert, steps=2, scorer=lambda s, a: [0.0] * len(a))
    wrong = replace(run.refinement, bounds=replace(
        run.refinement.bounds,
        partition_upper=run.refinement.bounds.partition_upper + 1,
    ))
    with pytest.raises(ValueError):
        verify_guided_run(replace(run, refinement=wrong))


def test_offline_jepa_training_smoke():
    pytest.importorskip("torch")
    data = collect_checked_decisions(("HHHHHH", "HPHHPPHH"), steps=2)
    model, losses = fit_checked_jepa(data, epochs=3, seed=1)
    assert len(losses) == 3 and all(math.isfinite(x) for x in losses)
    cert = build_frontier("HPHHPPHH", Fraction(2), 1, 2)
    run = guided_refine(cert, steps=3, scorer=model.score)
    assert verify_guided_run(run) == run.refinement.bounds


def test_online_reinforce_real_checked_rewards_smoke():
    pytest.importorskip("torch")
    model, rewards, calls = train_online_jepa_rl(
        ("HHHHHH", "HPHHPPHH"), episodes=2, steps=2, seed=3,
    )
    assert len(rewards) == 2 and all(math.isfinite(r) for r in rewards)
    assert calls >= 2
    cert = build_frontier("HHHHHH", Fraction(2), 1, 2)
    run = guided_refine(cert, steps=2, scorer=model.score)
    assert verify_guided_run(run) == run.refinement.bounds
