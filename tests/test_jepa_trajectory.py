"""Geometric JEPA respects symmetry and rejects invalid temporal data."""
from dataclasses import replace
import math

import pytest

from foldbounds.jepa_trajectory import (
    Frame, Trajectory, evaluate_trajectory_jepa, fit_trajectory_jepa,
    load_npz, observables, trajectory_pairs, validate_trajectory,
)


def _trajectory(source, shift=0.0):
    frames = tuple(
        Frame(float(t), tuple(
            (i * 1.5 + shift + t * 0.02, math.sin(i) * 0.25,
             math.cos(i) * 0.3)
            for i in range(8)
        ))
        for t in range(5)
    )
    return Trajectory(source, frames)


def test_geometry_invariant_to_rigid_rotation_and_translation():
    frame = _trajectory("source").frames[0]
    transformed = Frame(
        frame.time_ps,
        tuple((-y + 24.0, x - 13.0, z + 8.0)
              for x, y, z in frame.xyz_angstrom),
    )
    for x, y in zip(observables(frame), observables(transformed)):
        assert x == pytest.approx(y, abs=1e-12)


def test_temporal_pairs_respect_replica_boundaries():
    a, b = _trajectory("a"), _trajectory("b", shift=3.0)
    pairs = trajectory_pairs((a, b), max_stride=2)
    assert len(pairs) == 14
    assert {record[0] for record in pairs} == {"a", "b"}
    assert all(len(record[1]) == 10 and len(record[2]) == 9 for record in pairs)


def test_reject_bad_time_and_coordinates():
    original = _trajectory("id")
    validate_trajectory(original)
    with pytest.raises(ValueError):
        validate_trajectory(replace(original, source_id=""))
    with pytest.raises(ValueError):
        validate_trajectory(replace(original, frames=(
            original.frames[0], original.frames[0],
        )))
    corrupted = replace(original.frames[1], xyz_angstrom=(
        *original.frames[1].xyz_angstrom[:-1], (math.nan, 0.0, 0.0),
    ))
    with pytest.raises(ValueError):
        validate_trajectory(replace(original, frames=(
            original.frames[0], corrupted, *original.frames[2:],
        )))


def test_npz_loader_roundtrip(tmp_path):
    np = pytest.importorskip("numpy")
    original = _trajectory("replica-1")
    path = tmp_path / "md.npz"
    np.savez(path,
             time_ps=[f.time_ps for f in original.frames],
             xyz_angstrom=[f.xyz_angstrom for f in original.frames])
    loaded = load_npz(str(path), source_id=original.source_id)
    assert trajectory_pairs((loaded,)) == trajectory_pairs((original,))


def test_trajectory_jepa_smoke_without_physical_claim():
    pytest.importorskip("torch")
    train = _trajectory("train-replica")
    test = _trajectory("heldout-replica", shift=0.4)
    model, loss = fit_trajectory_jepa((train,), epochs=3, seed=4)
    report = evaluate_trajectory_jepa(model, (test,))
    assert len(loss) == 3 and all(math.isfinite(value) for value in loss)
    assert report["pairs"] == 7
    assert math.isfinite(report["model_mse"])
    assert math.isfinite(report["persistence_mse"])
    # It may be worse than persistence. This is NOT a demonstrated predictor.
