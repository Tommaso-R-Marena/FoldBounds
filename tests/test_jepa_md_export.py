"""Technical MD file conversion fixture; NOT physical trajectory validation."""
import importlib.util
import math
from pathlib import Path

import pytest


def _exporter():
    path = Path(__file__).resolve().parents[1] / "scripts/export_jepa_md.py"
    spec = importlib.util.spec_from_file_location("jepa_md_export", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_mdtraj_export_units_time_and_checksum(tmp_path):
    np = pytest.importorskip("numpy")
    md = pytest.importorskip("mdtraj")
    top = md.Topology()
    chain = top.add_chain()
    for i in range(8):
        residue = top.add_residue("ALA", chain, resSeq=i + 1)
        top.add_atom("CA", md.element.carbon, residue)
    coords = np.array([
        [[0.38 * i, 0.1 * math.sin(i), 0.05 * math.cos(i + t / 10)]
         for i in range(8)]
        for t in range(5)
    ], dtype=np.float32)
    traj = md.Trajectory(
        coords, top, time=np.array([0, 2, 4, 6, 8], dtype=np.float32),
    )
    md_path, topology_path = tmp_path / "frames.xtc", tmp_path / "top.pdb"
    traj.save_xtc(str(md_path))
    traj[0].save_pdb(str(topology_path))
    output = tmp_path / "md.npz"
    exporter = _exporter()
    info = exporter.export_md(
        str(md_path), str(topology_path), str(output),
        source_id="synthetic-format-test",
    )
    with np.load(output, allow_pickle=False) as data:
        assert data["time_ps"].tolist() == [0, 2, 4, 6, 8]
        assert data["xyz_angstrom"].shape == (5, 8, 3)
        assert np.allclose(data["xyz_angstrom"][0, :, 0],
                           np.arange(8) * 3.8, atol=1e-2)
    assert info["trajectory_sha256"] == exporter._sha256(md_path)
    assert info["units"] == {"time": "ps", "xyz": "Angstrom"}
    assert output.with_suffix(".npz.json").is_file()
    with pytest.raises(FileExistsError):
        exporter.export_md(
            str(md_path), str(topology_path), str(output),
            source_id="synthetic-format-test",
        )


def test_reject_unordered_nmr_models_before_loading(tmp_path):
    exporter = _exporter()
    with pytest.raises(ValueError, match="not chronological"):
        exporter.export_md(
            str(tmp_path / "nmr.pdb"), str(tmp_path / "top.pdb"),
            str(tmp_path / "output.npz"), source_id="bad",
        )
    with pytest.raises(ValueError, match="timestep_ps"):
        exporter.export_md(
            str(tmp_path / "dynamics.dcd"), str(tmp_path / "top.pdb"),
            str(tmp_path / "output.npz"), source_id="bad",
        )
