"""Experimental time-conditioned geometry-invariant JEPA for MD trajectories.

This is NOT a physical kinetic certificate or validated folding predictor.
Frames must come from genuine time-ordered dynamics, not unordered NMR models.
Generated predictions never enter FoldBounds' rigorous certificate checks.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import random

from .jepa_guidance import STATE_DIM, ACTION_DIM, build_model


@dataclass(frozen=True)
class Frame:
    time_ps: float
    xyz_angstrom: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class Trajectory:
    source_id: str
    frames: tuple[Frame, ...]


def validate_trajectory(trajectory: Trajectory) -> None:
    if not trajectory.source_id or len(trajectory.frames) < 2:
        raise ValueError("trajectory requires source_id and >=2 time-ordered frames")
    n = len(trajectory.frames[0].xyz_angstrom)
    if n < 5:
        raise ValueError("need at least five ordered atoms")
    last_time = float("-inf")
    for frame in trajectory.frames:
        if not math.isfinite(frame.time_ps) or frame.time_ps <= last_time:
            raise ValueError("trajectory time_ps must increase strictly")
        last_time = frame.time_ps
        if len(frame.xyz_angstrom) != n or any(
            len(point) != 3 or not all(math.isfinite(x) for x in point)
            for point in frame.xyz_angstrom
        ):
            raise ValueError("coordinates must be finite fixed-topology 3D triplets")


def observables(frame: Frame) -> tuple[float, ...]:
    """Distance observables invariant to rigid rotation and translation.

    Contact thresholds are merely geometric features, not certified contacts.
    """
    coords = frame.xyz_angstrom
    n = len(coords)
    if n < 5:
        raise ValueError("at least five atoms required")
    if any(len(p) != 3 or not all(math.isfinite(v) for v in p) for p in coords):
        raise ValueError("nonfinite or malformed coordinates")
    centroid = tuple(sum(p[k] for p in coords) / n for k in range(3))
    radius = math.sqrt(sum(
        sum((p[k] - centroid[k]) ** 2 for k in range(3))
        for p in coords
    ) / n)

    def distance(i, j):
        return math.sqrt(sum(
            (coords[i][k] - coords[j][k]) ** 2 for k in range(3)
        ))

    remote = [distance(i, j) for i in range(n) for j in range(i + 4, n)]
    near = []
    for offset in (1, 2, 4):
        ds = [distance(i, i + offset) for i in range(n - offset)]
        near.append(sum(ds) / len(ds) / 10.0)
    occupancy = tuple(
        sum(d < cutoff for d in remote) / max(1, len(remote))
        for cutoff in (4.0, 6.0, 8.0, 12.0)
    )
    result = (
        math.log1p(n) / 4.0, radius / (1.0 + radius),
        *occupancy, *near, sum(remote) / max(1, len(remote)) / 20.0,
    )
    assert len(result) == STATE_DIM
    return result


def trajectory_pairs(
    trajectories: tuple[Trajectory, ...], *, max_stride: int = 2,
):
    """Chronological pairs without mixing replicas; time delta conditions JEPA."""
    if type(max_stride) is not int or max_stride < 1:
        raise ValueError("max_stride must be a positive integer")
    records = []
    for trajectory in trajectories:
        validate_trajectory(trajectory)
        frames = trajectory.frames
        states = tuple(observables(frame) for frame in frames)
        for i in range(len(frames) - 1):
            for stride in range(
                1, min(max_stride, len(frames) - i - 1) + 1
            ):
                delta = frames[i + stride].time_ps - frames[i].time_ps
                action = (
                    math.log1p(delta), stride / max_stride,
                ) + (0.0,) * (ACTION_DIM - 2)
                records.append((
                    trajectory.source_id, states[i],
                    action, states[i + stride],
                ))
    return tuple(records)


def load_npz(path: str, *, source_id: str) -> Trajectory:
    """Import OpenMM/MD exports: time_ps[T], xyz_angstrom[T,N,3].

    Units and source provenance must be verified by the producing simulation.
    """
    try:
        import numpy as np
    except ImportError as exc:
        raise ImportError("numpy required for NPZ trajectory loading") from exc
    with np.load(path, allow_pickle=False) as data:
        times = data["time_ps"]
        xyz = data["xyz_angstrom"]
        if (
            times.ndim != 1 or xyz.ndim != 3 or
            xyz.shape != (len(times), xyz.shape[1], 3)
        ):
            raise ValueError("expected time_ps[T] and xyz_angstrom[T,N,3]")
        frames = tuple(
            Frame(float(time), tuple(tuple(map(float, point)) for point in coords))
            for time, coords in zip(times, xyz)
        )
    trajectory = Trajectory(source_id, frames)
    validate_trajectory(trajectory)
    return trajectory


def fit_trajectory_jepa(
    trajectories: tuple[Trajectory, ...], *, epochs: int = 40,
    seed: int = 0, lr: float = 1e-3,
    max_stride: int = 2, mask_probability: float = 0.15, model=None,
):
    """Train masked, time-conditioned successor embedding with EMA target.

    Evaluate on entire held-out replicas, not random frames from the same run.
    This is a forecasting objective, not a rate or rare-event guarantee.
    """
    import torch
    from torch import nn
    if (
        type(epochs) is not int or epochs < 1 or lr <= 0 or
        not (0 <= mask_probability < 1)
    ):
        raise ValueError("invalid epochs, learning rate or mask probability")
    pairs = trajectory_pairs(trajectories, max_stride=max_stride)
    if not pairs:
        raise ValueError("no genuine chronological pairs")
    torch.manual_seed(seed)
    rng = random.Random(seed)
    model = build_model() if model is None else model
    optimizer = torch.optim.AdamW(
        (p for p in model.parameters() if p.requires_grad), lr=lr,
    )
    source = torch.tensor([r[1] for r in pairs], dtype=torch.float32)
    action = torch.tensor([r[2] for r in pairs], dtype=torch.float32)
    successor = torch.tensor([r[3] for r in pairs], dtype=torch.float32)
    losses = []
    for _ in range(epochs):
        masked = source.clone()
        mask = torch.tensor(
            [
                [rng.random() < mask_probability for _ in range(STATE_DIM)]
                for _ in range(len(pairs))
            ], dtype=torch.bool,
        )
        masked[mask] = 0.0
        _, _, prediction = model.heads(masked, action)
        with torch.no_grad():
            target = model.target_encoder(successor)
        loss = (
            nn.functional.smooth_l1_loss(prediction, target)
            + nn.functional.mse_loss(model.decoder(prediction), successor)
        )
        if not torch.isfinite(loss):
            raise ArithmeticError("nonfinite trajectory JEPA loss")
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(
            (p for p in model.parameters() if p.requires_grad), 2.0,
        )
        optimizer.step()
        model.update_target()
        losses.append(float(loss.detach()))
    return model, tuple(losses)


def evaluate_trajectory_jepa(
    model, trajectories: tuple[Trajectory, ...], *, max_stride: int = 2,
):
    """Held-out next-observable MSE versus no-change persistence baseline."""
    import torch
    pairs = trajectory_pairs(trajectories, max_stride=max_stride)
    if not pairs:
        raise ValueError("no evaluation pairs")
    source = torch.tensor([r[1] for r in pairs], dtype=torch.float32)
    action = torch.tensor([r[2] for r in pairs], dtype=torch.float32)
    target = torch.tensor([r[3] for r in pairs], dtype=torch.float32)
    model.eval()
    with torch.no_grad():
        _, _, prediction = model.heads(source, action)
        decoded = model.decoder(prediction)
        model_mse = float((decoded - target).square().mean())
        baseline_mse = float((source - target).square().mean())
    return {
        "pairs": len(pairs),
        "model_mse": model_mse,
        "persistence_mse": baseline_mse,
    }
