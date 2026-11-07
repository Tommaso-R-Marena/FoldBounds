"""Export *time-ordered* MD frames to FoldBounds JEPA NPZ with provenance.

Requires optional mdtraj/numpy. PDB/NMR MODEL collections are NOT chronological
MD. MDTraj uses nanometers; exported xyz_angstrom uses Angstroms.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def export_md(
    trajectory_path: str, topology_path: str, output_path: str, *,
    source_id: str, selection: str = "protein and name CA",
    stride: int = 1, timestep_ps: float | None = None,
    first_time_ps: float = 0.0,
) -> dict:
    import mdtraj as md
    import numpy as np
    if type(stride) is not int or stride < 1 or not source_id:
        raise ValueError("stride must be positive and source ID nonempty")
    if timestep_ps is not None and (
        not math.isfinite(timestep_ps) or timestep_ps <= 0
    ):
        raise ValueError("timestep_ps must be finite positive")
    if not math.isfinite(first_time_ps):
        raise ValueError("first_time_ps must be finite")
    source = Path(trajectory_path)
    topology = Path(topology_path)
    output = Path(output_path)
    if source.suffix.lower() in (".pdb", ".cif", ".mmcif"):
        raise ValueError("PDB/NMR ensembles are not chronological MD trajectories")
    if source.suffix.lower() == ".dcd" and timestep_ps is None:
        raise ValueError("supply original frame timestep_ps explicitly for DCD")
    if output.exists() or output.with_suffix(output.suffix + ".json").exists():
        raise FileExistsError("output and provenance metadata must be new")
    trajectory = md.load(str(source), top=str(topology), stride=stride)
    atom_ids = trajectory.topology.select(selection)
    if len(atom_ids) < 5:
        raise ValueError("select at least five ordered atoms")
    if trajectory.n_frames < 2:
        raise ValueError("need at least two time-ordered MD frames")
    xyz = np.asarray(trajectory.xyz[:, atom_ids, :], dtype=np.float64) * 10.0
    if timestep_ps is None:
        time = np.asarray(trajectory.time, dtype=np.float64)
    else:
        time = first_time_ps + (
            np.arange(trajectory.n_frames, dtype=np.float64)
            * timestep_ps * stride
        )
    if (
        not np.isfinite(time).all() or not np.isfinite(xyz).all()
        or not np.all(np.diff(time) > 0)
    ):
        raise ValueError("nonfinite coordinates or invalid time order")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, time_ps=time, xyz_angstrom=xyz)
    metadata = {
        "source_id": source_id,
        "trajectory_sha256": _sha256(source),
        "topology_sha256": _sha256(topology),
        "topology_selection": selection,
        "selected_atom_count": len(atom_ids),
        "frame_count": trajectory.n_frames,
        "stride": stride,
        "units": {"time": "ps", "xyz": "Angstrom"},
        "time_source": (
            "explicit-frame-spacing" if timestep_ps is not None
            else "mdtraj-trajectory-time"
        ),
        "first_time_ps": float(time[0]),
        "last_time_ps": float(time[-1]),
        "mdtraj_version": md.__version__,
        "prediction_claim": "none",
    }
    output.with_suffix(output.suffix + ".json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trajectory", required=True)
    parser.add_argument("--topology", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--selection", default="protein and name CA")
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--timestep-ps", type=float, default=None)
    parser.add_argument("--first-time-ps", type=float, default=0.0)
    args = parser.parse_args()
    print(json.dumps(export_md(
        args.trajectory, args.topology, args.output,
        source_id=args.source_id,
        selection=args.selection, stride=args.stride,
        timestep_ps=args.timestep_ps, first_time_ps=args.first_time_ps,
    ), indent=2))


if __name__ == "__main__":
    main()
