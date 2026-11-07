"""Opt-in JEPA and checker-rewarded RL for *search*, never proof authority.

Actions expand finite-model certificate frontiers. These are computational
interventions, NOT molecular dynamics or physical kinetic transitions.
Only FoldBounds' independent replay verifier authorizes reported intervals.
PyTorch is optional, and no neural output can be a certificate or bound.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import copy
import math
import random
from typing import Callable

from .adaptive import (
    RefinementRun, RefinementSnapshot, _intersect, verify_refinement_run,
)
from .lattice import (
    Certificate, Path, VerifiedBounds, _contacts, _neighbors,
    _possible_hh_pairs, refine_frontier,
)
from .lookahead import verify_lookahead

STATE_DIM = 10
ACTION_DIM = 9


def state_features(cert: Certificate, bounds: VerifiedBounds) -> tuple[float, ...]:
    """Fixed-size model-relative features. Never used to certify a result."""
    n = len(cert.sequence)
    frontier = cert.frontier
    return (
        n / 32.0,
        cert.sequence.count("H") / n,
        math.log1p(float(cert.q - 1)),
        cert.event_at_least / max(1, n),
        sum(x.kind == "bound" for x in frontier) / max(1, n),
        sum(x.kind == "exact" for x in frontier) / max(1, n),
        math.log1p(len(frontier)),
        float(bounds.probability_upper - bounds.probability_lower),
        math.log1p(float(bounds.partition_upper)),
        sum(len(x.path) for x in frontier) / (max(1, len(frontier)) * n),
    )


def action_features(cert: Certificate, path: Path) -> tuple[float, ...]:
    """Untrusted scheduling features without target or exact-enumeration labels."""
    n = len(cert.sequence)
    pairs = _possible_hh_pairs(cert.sequence)
    tip = path[-1]
    xs, ys = zip(*path)
    steps = len(path)
    return (
        steps / n,
        _contacts(path, pairs) / max(1, n),
        len(_neighbors(path)) / 3.0,
        (max(xs) - min(xs)) / n,
        (max(ys) - min(ys)) / n,
        (abs(tip[0]) + abs(tip[1])) / n,
        cert.sequence[:steps].count("H") / steps,
        1.0 if steps < n and cert.sequence[steps] == "H" else 0.0,
        (n - steps) / n,
    )


def candidate_paths(cert: Certificate) -> tuple[Path, ...]:
    """Canonical ordering stabilizes ties and seeded experiments."""
    return tuple(sorted(
        (node.path for node in cert.frontier if node.kind == "bound"),
        key=lambda path: (len(path), path),
    ))


@dataclass(frozen=True)
class CheckedCandidate:
    action: tuple[float, ...]
    successor: tuple[float, ...]
    gain: Fraction
    partition_gain: Fraction
    verifier_calls: int = 1

    @property
    def reward(self) -> float:
        """Verified-width reduction per checker call, used only for learning."""
        return (
            float(self.gain) + 0.01 * math.tanh(float(self.partition_gain))
        ) / max(1, self.verifier_calls)


@dataclass(frozen=True)
class CheckedDecision:
    state: tuple[float, ...]
    candidates: tuple[CheckedCandidate, ...]
    chosen: int
    sequence: str


@dataclass(frozen=True)
class GuidedRun:
    refinement: RefinementRun
    checker_calls: int
    decisions: tuple[CheckedDecision, ...]


def _settings(*, lookahead: int, use_matching: bool) -> dict:
    if type(lookahead) is not int or lookahead < 0 or type(use_matching) is not bool:
        raise ValueError("invalid checker settings")
    return {"lookahead": lookahead, "use_matching": use_matching}


def _outcome(
    cert: Certificate, bounds: VerifiedBounds, path: Path, settings: dict
) -> tuple[Certificate, VerifiedBounds, CheckedCandidate]:
    next_cert = refine_frontier(cert, path)
    checked = verify_lookahead(next_cert, **settings)
    updated = _intersect(bounds, checked)
    p_gain = (
        bounds.probability_upper - bounds.probability_lower
        - (updated.probability_upper - updated.probability_lower)
    )
    z_gain = (
        bounds.partition_upper - bounds.partition_lower
        - (updated.partition_upper - updated.partition_lower)
    )
    if p_gain < 0 or z_gain < 0:
        raise AssertionError("independently checked refinement widened interval")
    return next_cert, updated, CheckedCandidate(
        action_features(cert, path), state_features(next_cert, updated),
        p_gain, z_gain,
    )


def collect_checked_decisions(
    sequences: tuple[str, ...], *, q: Fraction = Fraction(2),
    threshold: int = 1, depth: int = 2, steps: int = 5,
    seed: int = 0, lookahead: int = 0, use_matching: bool = True,
) -> tuple[CheckedDecision, ...]:
    """Exploratory rollouts with independently checked *all-action* targets.

    Producing labels costs one verifier invocation for every candidate action.
    Keep whole sequences disjoint between training and evaluation splits.
    """
    from .lattice import build_frontier
    if type(steps) is not int or steps < 0:
        raise ValueError("steps must be a nonnegative integer")
    rng = random.Random(seed)
    settings = _settings(lookahead=lookahead, use_matching=use_matching)
    examples: list[CheckedDecision] = []
    for sequence in sequences:
        cert = build_frontier(sequence, q, threshold, depth)
        current = verify_lookahead(cert, **settings)
        for _ in range(steps):
            paths = candidate_paths(cert)
            if not paths:
                break
            outcomes = [_outcome(cert, current, path, settings) for path in paths]
            chosen = rng.randrange(len(paths))
            examples.append(CheckedDecision(
                state_features(cert, current),
                tuple(result[2] for result in outcomes), chosen, sequence,
            ))
            cert, current, _ = outcomes[chosen]
    return tuple(examples)


def guided_refine(
    initial: Certificate, *, steps: int,
    scorer: Callable[
        [tuple[float, ...], tuple[tuple[float, ...], ...]],
        tuple[float, ...],
    ] | None = None,
    lookahead: int = 0, use_matching: bool = True,
    collect_counterfactuals: bool = False,
) -> GuidedRun:
    """Schedule refinements, always rechecking the full certificate.

    Bad scorer outputs are rejected before any proposed transition. When no
    scorer is supplied, deterministically expand the first (breadth) node.
    """
    if type(steps) is not int or steps < 0:
        raise ValueError("steps must be a nonnegative integer")
    settings = _settings(lookahead=lookahead, use_matching=use_matching)
    cert = initial
    current = verify_lookahead(cert, **settings)
    checker_calls = 1
    frontiers = [cert]
    history = [RefinementSnapshot(
        0, None, current.bounded_subtrees, current.exact_paths,
        current.partition_lower, current.partition_upper,
        current.probability_lower, current.probability_upper,
    )]
    decisions: list[CheckedDecision] = []
    for step in range(1, steps + 1):
        paths = candidate_paths(cert)
        if not paths:
            break
        state = state_features(cert, current)
        actions = tuple(action_features(cert, path) for path in paths)
        if scorer is None:
            index = 0
        else:
            scores = tuple(scorer(state, actions))
            if len(scores) != len(paths) or any(
                not isinstance(score, (float, int)) or
                not math.isfinite(score) for score in scores
            ):
                raise ValueError("untrusted scorer returned invalid action scores")
            index = max(range(len(paths)), key=lambda i: (scores[i], -i))
        if collect_counterfactuals:
            outcomes = [_outcome(cert, current, path, settings) for path in paths]
            checker_calls += len(paths)
            next_cert, updated, chosen = outcomes[index]
            decisions.append(CheckedDecision(
                state, tuple(result[2] for result in outcomes), index, cert.sequence,
            ))
        else:
            next_cert, updated, chosen = _outcome(
                cert, current, paths[index], settings,
            )
            checker_calls += 1
            decisions.append(CheckedDecision(state, (chosen,), 0, cert.sequence))
        cert, current = next_cert, updated
        frontiers.append(cert)
        history.append(RefinementSnapshot(
            step, paths[index], current.bounded_subtrees, current.exact_paths,
            current.partition_lower, current.partition_upper,
            current.probability_lower, current.probability_upper,
        ))
    # The policy string is metadata; the independent verifier ignores it.
    run = RefinementRun(
        cert, current, tuple(history), "jepa", tuple(frontiers),
        lookahead, use_matching, False, 0, False, 2,
    )
    verify_refinement_run(run)
    return GuidedRun(run, checker_calls + len(frontiers), tuple(decisions))


def verify_guided_run(run: GuidedRun) -> VerifiedBounds:
    """No learned score is used by the independent audit replay."""
    return verify_refinement_run(run.refinement)


def build_model(*, latent_dim: int = 32, hidden_dim: int = 64):
    """Task-conditioned JEPA, EMA teacher, decoder, actor, and critic.

    PyTorch is imported lazily. Reconstructions and checked-gain losses provide
    additional training signals; they are not mathematical guarantees.
    """
    try:
        import torch
        from torch import nn
    except ImportError as exc:
        raise ImportError("Install PyTorch to use JEPA learning") from exc
    if latent_dim < 2 or hidden_dim < 2:
        raise ValueError("latent_dim and hidden_dim must be >=2")

    class JEPA(nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = nn.Sequential(
                nn.Linear(STATE_DIM, hidden_dim), nn.GELU(),
                nn.Linear(hidden_dim, latent_dim),
            )
            self.target_encoder = copy.deepcopy(self.encoder)
            for parameter in self.target_encoder.parameters():
                parameter.requires_grad_(False)
            self.action_encoder = nn.Sequential(
                nn.Linear(ACTION_DIM, hidden_dim), nn.GELU(),
                nn.Linear(hidden_dim, latent_dim),
            )
            self.predictor = nn.Sequential(
                nn.Linear(2 * latent_dim, hidden_dim), nn.GELU(),
                nn.Linear(hidden_dim, latent_dim),
            )
            self.decoder = nn.Linear(latent_dim, STATE_DIM)
            self.actor = nn.Sequential(
                nn.Linear(2 * latent_dim, hidden_dim), nn.GELU(),
                nn.Linear(hidden_dim, 1),
            )
            self.critic = nn.Sequential(
                nn.Linear(2 * latent_dim, hidden_dim), nn.GELU(),
                nn.Linear(hidden_dim, 1),
            )

        def heads(self, states, actions):
            state_z = self.encoder(states)
            action_z = self.action_encoder(actions)
            conditioned = torch.cat((state_z, action_z), dim=-1)
            return (
                self.actor(conditioned).squeeze(-1),
                self.critic(conditioned).squeeze(-1),
                self.predictor(conditioned),
            )

        @torch.no_grad()
        def update_target(self, momentum: float = 0.99):
            for target, online in zip(
                self.target_encoder.parameters(), self.encoder.parameters(),
            ):
                target.lerp_(online, 1.0 - momentum)

        @torch.no_grad()
        def score(self, state, actions) -> tuple[float, ...]:
            self.eval()
            states = torch.tensor([state] * len(actions), dtype=torch.float32)
            choices = torch.tensor(actions, dtype=torch.float32)
            scores, _, _ = self.heads(states, choices)
            return tuple(float(value) for value in scores)

    return JEPA()


def fit_checked_jepa(
    decisions: tuple[CheckedDecision, ...], *, epochs: int = 40,
    lr: float = 1e-3, seed: int = 0, model=None,
):
    """Offline JEPA representation + checker-rewarded counterfactual ranking."""
    import torch
    from torch import nn
    if not decisions or type(epochs) is not int or epochs < 1 or lr <= 0:
        raise ValueError("nonempty decisions, positive epochs and lr required")
    torch.manual_seed(seed)
    model = build_model() if model is None else model
    optimizer = torch.optim.AdamW(
        (p for p in model.parameters() if p.requires_grad), lr=lr,
    )
    states, actions, successors, rewards = [], [], [], []
    groups = []
    for record in decisions:
        start = len(states)
        for candidate in record.candidates:
            states.append(record.state)
            actions.append(candidate.action)
            successors.append(candidate.successor)
            rewards.append(candidate.reward)
        groups.append((start, len(states)))
    states_t = torch.tensor(states, dtype=torch.float32)
    actions_t = torch.tensor(actions, dtype=torch.float32)
    successors_t = torch.tensor(successors, dtype=torch.float32)
    rewards_t = torch.tensor(rewards, dtype=torch.float32)
    losses = []
    for _ in range(epochs):
        model.train()
        actor, critic, prediction = model.heads(states_t, actions_t)
        with torch.no_grad():
            target = model.target_encoder(successors_t)
        prediction_loss = nn.functional.smooth_l1_loss(prediction, target)
        reconstruction_loss = nn.functional.mse_loss(
            model.decoder(prediction), successors_t,
        )
        value_loss = nn.functional.smooth_l1_loss(critic, rewards_t)
        ranking_losses = []
        for start, end in groups:
            if end - start < 2:
                continue
            best = start + int(rewards_t[start:end].argmax())
            ranking_losses.append(nn.functional.cross_entropy(
                actor[start:end][None, :],
                torch.tensor([best - start]),
            ))
        ranking_loss = (
            torch.stack(ranking_losses).mean()
            if ranking_losses else actor.sum() * 0.0
        )
        loss = prediction_loss + reconstruction_loss + value_loss + ranking_loss
        if not torch.isfinite(loss):
            raise ArithmeticError("non-finite JEPA loss")
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(
            (p for p in model.parameters() if p.requires_grad), 2.0,
        )
        optimizer.step()
        model.update_target()
        losses.append(float(loss.detach()))
    return model, tuple(losses)


def train_online_jepa_rl(
    sequences: tuple[str, ...], *, episodes: int = 8,
    steps: int = 5, seed: int = 0,
    q: Fraction = Fraction(2), threshold: int = 1, depth: int = 2,
    lookahead: int = 0, use_matching: bool = True,
    lr: float = 1e-3, model=None,
):
    """Actual checked online REINFORCE/critic + JEPA successor training.

    No counterfactual reward is presumed: only the sampled action is checked.
    Each episode returns its own sum of verifier-derived training rewards.
    """
    from .lattice import build_frontier
    import torch
    from torch import nn
    if (
        not sequences or type(episodes) is not int or episodes < 1 or
        type(steps) is not int or steps < 0 or lr <= 0
    ):
        raise ValueError("nonempty sequences and valid episodes/steps/lr required")
    torch.manual_seed(seed)
    rng = random.Random(seed)
    model = build_model() if model is None else model
    optimizer = torch.optim.AdamW(
        (p for p in model.parameters() if p.requires_grad), lr=lr,
    )
    settings = _settings(lookahead=lookahead, use_matching=use_matching)
    returns = []
    checker_calls = 0
    for _ in range(episodes):
        sequence = rng.choice(sequences)
        cert = build_frontier(sequence, q, threshold, depth)
        bounds = verify_lookahead(cert, **settings)
        checker_calls += 1
        rewards = []
        for _ in range(steps):
            paths = candidate_paths(cert)
            if not paths:
                break
            state = state_features(cert, bounds)
            states_t = torch.tensor(
                [state] * len(paths), dtype=torch.float32,
            )
            actions_t = torch.tensor(
                [action_features(cert, path) for path in paths],
                dtype=torch.float32,
            )
            model.train()
            logits, values, predictions = model.heads(states_t, actions_t)
            distribution = torch.distributions.Categorical(logits=logits)
            index = int(distribution.sample())
            next_cert, next_bounds, outcome = _outcome(
                cert, bounds, paths[index], settings,
            )
            checker_calls += 1
            rewards.append(outcome.reward)
            reward_t = torch.tensor(outcome.reward, dtype=torch.float32)
            advantage = reward_t - values[index]
            actor_loss = (
                -distribution.log_prob(torch.tensor(index))
                * advantage.detach()
            )
            critic_loss = (values[index] - reward_t).square()
            successor = torch.tensor(
                [outcome.successor], dtype=torch.float32,
            )
            with torch.no_grad():
                target = model.target_encoder(successor)[0]
            predictive = nn.functional.smooth_l1_loss(predictions[index], target)
            reconstruction = nn.functional.mse_loss(
                model.decoder(predictions[index]), successor[0],
            )
            loss = (
                actor_loss + 0.5 * critic_loss + predictive + reconstruction
                - 0.01 * distribution.entropy()
            )
            if not torch.isfinite(loss):
                raise ArithmeticError("non-finite RL training loss")
            optimizer.zero_grad()
            loss.backward()
            nn.utils.clip_grad_norm_(
                (p for p in model.parameters() if p.requires_grad), 2.0,
            )
            optimizer.step()
            model.update_target()
            cert, bounds = next_cert, next_bounds
        returns.append(sum(rewards))
    return model, tuple(returns), checker_calls
