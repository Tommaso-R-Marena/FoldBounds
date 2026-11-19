#!/usr/bin/env python3
"""Run ONE verification gate claimed by draft PR #32 and record what happened.

This is a measurement/replication harness, not a proof and not a fixer: every
command is executed exactly as written in ``GATES`` from the repository root,
its combined stdout/stderr goes to ``/tmp/baseline_v2/logs/gate_<id>.log`` and a
JSON record goes to ``/tmp/baseline_v2/gate_<id>.json``.  Nothing in the
repository is modified (``PYTHONDONTWRITEBYTECODE=1``, pytest cache plugin
disabled through ``PYTEST_ADDOPTS``, git called with ``--no-optional-locks``).

    python3 replicate_gates.py --gate 5          # run one gate
    python3 replicate_gates.py --gate 1p         # variant on the pristine export
    python3 replicate_gates.py --list
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/workspace")
SCRATCH = Path("/tmp/baseline_v2")
LOGS = SCRATCH / "logs"
PRISTINE = SCRATCH / "pristine"          # `git archive fc3de05 | tar -x`
EXPECTED_HEAD = "fc3de052d67f54e0316d56ef7577c272d4afd02d"
MAX_SECONDS = 1200                       # hard rule: nothing longer than 20 minutes

RG_PATTERN = r'native_decide|sorry|admit|implemented_by|unsafe|@\[extern|axiom '
RG_CMD = f'rg -n "{RG_PATTERN}" FoldBounds scripts/*.lean'

# id -> (title, [(step label, shell command)], cwd, extra env, decisive-line regex)
GATES: dict[str, dict] = {
    "1": dict(title="Lean source hygiene scan",
              steps=[("scan", "python3 scripts/check_lean_sources.py")],
              decisive=r"^Source audit: .*"),
    "1p": dict(title="Lean source hygiene scan (pristine `git archive` of fc3de05)",
               steps=[("scan", "python3 scripts/check_lean_sources.py")],
               cwd=str(PRISTINE), decisive=r"^Source audit: .*"),
    "2": dict(title="Executable certificate accept/reject tests",
              steps=[("lean", "lake env lean FoldBounds/ExecutableTests.lean")],
              decisive=None),
    "3": dict(title="Concrete Lean HP lattice examples",
              steps=[("lean", "lake env lean FoldBounds/LatticeModelTests.lean")],
              decisive=None),
    "4": dict(title="Cross-language (Python vs Lean) lattice semantics differential",
              steps=[("diff", "PYTHONPATH=src python3 scripts/cross_language_diff.py")],
              decisive=r"^Lean/Python differential checks passed: .*"),
    "5": dict(title="Kernel-checked version-1 frontier tests",
              steps=[("lean", "lake env lean FoldBounds/FrontierTests.lean")],
              decisive=None),
    "6": dict(title="Differential Python vs Lean version-1 frontier checker (interpreted runner)",
              steps=[("diff", "PYTHONPATH=src python3 experiments/lean_differential.py "
                              "--lake-root . --runner scripts/foldcheck.lean")],
              decisive=r"^TOTAL PROBLEMS: .*"),
    "7": dict(title="Kernel-checked version-2 frontier tests",
              steps=[("lean", "lake env lean FoldBounds/FrontierNextTests.lean")],
              decisive=None),
    "8": dict(title="Version-2 exact-oracle differential, adversarial catalogue and mutants "
                    "(compiled foldcheck2)",
              steps=[("adv", "PYTHONPATH=src python3 experiments/lean_v2_adversarial.py "
                             "--exe2 .lake/build/bin/foldcheck2 --mutants 100 "
                             "--out /tmp/baseline_v2/v2adv")],
              decisive=r"^(differential|adversarial|mutation): .*"),
    "9": dict(title="Axiom audit of headline theorems + repository gate script",
              steps=[("print-axioms",
                      "lake env lean FoldBounds/AxiomAudit.lean > /tmp/baseline_v2/axioms.log"),
                     ("gate-script",
                      "python3 scripts/check_lean_sources.py /tmp/baseline_v2/axioms.log")],
              decisive=r"^(Source audit|Axiom audit): .*"),
    "10": dict(title="Full Python test suite",
               steps=[("pytest", "python -m pytest -q")],
               env={"PYTEST_ADDOPTS": "-p no:cacheprovider"},
               decisive=r"^(=+ )?.*\b(passed|failed|error|errors)\b.* in [0-9.]+s.*"),
    "10p": dict(title="Full Python test suite (pristine `git archive` of fc3de05, PYTHONPATH=src)",
                steps=[("pytest", "python -m pytest -q")], cwd=str(PRISTINE),
                env={"PYTEST_ADDOPTS": "-p no:cacheprovider", "PYTHONPATH": "src"},
                decisive=r"^(=+ )?.*\b(passed|failed|error|errors)\b.* in [0-9.]+s.*"),
    "11": dict(title="Forbidden-token ripgrep scan (live checkout)",
               steps=[("rg", RG_CMD)], decisive=None),
    "11p": dict(title="Forbidden-token ripgrep scan (pristine `git archive` of fc3de05)",
                steps=[("rg", RG_CMD)], cwd=str(PRISTINE), decisive=None),
}


def base_env() -> dict:
    env = dict(os.environ)
    home = env.get("HOME", "/home/ubuntu")
    env["PATH"] = f"{ROOT}/.venv/bin:{home}/.elan/bin:" + env.get("PATH", "")
    env["VIRTUAL_ENV"] = str(ROOT / ".venv")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["GIT_OPTIONAL_LOCKS"] = "0"
    env.pop("PYTHONPATH", None)          # gates set it themselves where the recipe does
    return env


def sh(cmd: str, cwd: Path = ROOT) -> str:
    return subprocess.run(["bash", "-c", cmd], cwd=cwd, env=base_env(), capture_output=True,
                          text=True).stdout.strip()


def snapshot() -> dict:
    """Read-only description of the shared checkout at this instant."""
    status = sh("git --no-optional-locks status --porcelain --untracked-files=all")
    oleans = sorted((ROOT / ".lake/build/lib/lean/FoldBounds").glob("*.olean"))
    newest = max(oleans, key=lambda p: p.stat().st_mtime) if oleans else None
    bins = {}
    for name in ("foldcheck", "foldcheck2"):
        p = ROOT / ".lake/build/bin" / name
        st = p.stat()
        bins[name] = {"bytes": st.st_size, "mtime": dt.datetime.fromtimestamp(
            st.st_mtime, dt.timezone.utc).isoformat(timespec="seconds")}
    tracked_dirty = sh("git --no-optional-locks diff --name-only HEAD -- FoldBounds FoldBounds.lean "
                       "scripts src tests experiments lakefile.toml lake-manifest.json "
                       "lean-toolchain pyproject.toml")
    return {
        "utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
        "head": sh("git --no-optional-locks rev-parse HEAD"),
        "branch": sh("git --no-optional-locks rev-parse --abbrev-ref HEAD"),
        "git_status_porcelain": status.splitlines(),
        "tracked_gate_inputs_modified": tracked_dirty.splitlines(),
        "lean_files_in_FoldBounds": len(list((ROOT / "FoldBounds").glob("*.lean"))),
        "olean_count": len(oleans),
        "newest_olean": None if newest is None else {
            "file": newest.name, "mtime": dt.datetime.fromtimestamp(
                newest.stat().st_mtime, dt.timezone.utc).isoformat(timespec="seconds")},
        "binaries": bins,
        "loadavg": Path("/proc/loadavg").read_text().split()[:3],
    }


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run_gate(gid: str) -> dict:
    g = GATES[gid]
    LOGS.mkdir(parents=True, exist_ok=True)
    cwd = Path(g.get("cwd", str(ROOT)))
    env = base_env()
    env.update(g.get("env", {}))
    before = snapshot()
    steps = []
    remaining = float(MAX_SECONDS)
    for label, cmd in g["steps"]:
        log = LOGS / f"gate_{gid}_{label}.log"
        budget = max(1, int(remaining))
        argv = ["timeout", "--signal=TERM", "--kill-after=20", str(budget), "bash", "-c", cmd]
        r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
        t0 = time.monotonic()
        with open(log, "wb") as fh:
            proc = subprocess.run(argv, cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT)
        wall = time.monotonic() - t0
        r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        remaining -= wall
        text = log.read_text(encoding="utf-8", errors="replace")
        lines = text.splitlines()
        steps.append({
            "label": label, "command": cmd, "cwd": str(cwd), "timeout_s": budget,
            "exit_code": proc.returncode, "timed_out": proc.returncode in (124, 137),
            "wall_s": round(wall, 3),
            "child_user_s": round(r1.ru_utime - r0.ru_utime, 3),
            "child_sys_s": round(r1.ru_stime - r0.ru_stime, 3),
            "log": str(log), "log_bytes": len(text.encode()), "log_lines": len(lines),
            "log_sha256": sha256(log),
            "lines_with_error": sum(1 for l in lines if re.search(r"\berror\b", l, re.I)),
            "lines_with_warning": sum(1 for l in lines if re.search(r"\bwarning\b", l, re.I)),
            "lines_with_sorry": sum(1 for l in lines if "sorry" in l),
            "decisive_lines": ([l for l in lines if re.match(g["decisive"], l)]
                               if g.get("decisive") else []),
            "last_lines": lines[-5:],
        })
        if proc.returncode != 0 and gid not in ("11", "11p"):
            break                           # like `&&`/CI: later steps of a gate need earlier ones
    after = snapshot()
    rec = {
        "gate": gid, "title": g["title"], "extra_env": g.get("env", {}),
        "exit_code": steps[-1]["exit_code"],      # last executed step (a failing step stops the gate)
        "steps_executed": len(steps), "steps_planned": len(g["steps"]),
        "all_steps_exit_zero": all(s["exit_code"] == 0 for s in steps) and len(steps) == len(g["steps"]),
        "wall_s": round(sum(s["wall_s"] for s in steps), 3),
        "steps": steps, "before": before, "after": after,
        "checkout_changed_during_gate": (
            before["git_status_porcelain"] != after["git_status_porcelain"]
            or before["newest_olean"] != after["newest_olean"]
            or before["binaries"] != after["binaries"] or before["head"] != after["head"]),
        "head_is_expected": before["head"] == EXPECTED_HEAD == after["head"],
    }
    (SCRATCH / f"gate_{gid}.json").write_text(json.dumps(rec, indent=1))
    return rec


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate")
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()
    if args.list or not args.gate:
        for gid, g in GATES.items():
            print(gid, "|", g["title"], "|", " && ".join(c for _, c in g["steps"]))
        return
    rec = run_gate(args.gate)
    print(f"gate {rec['gate']}: {rec['title']}")
    for s in rec["steps"]:
        print(f"  [{s['label']}] exit={s['exit_code']} wall={s['wall_s']}s user={s['child_user_s']}s "
              f"sys={s['child_sys_s']}s log_lines={s['log_lines']} errors={s['lines_with_error']} "
              f"warnings={s['lines_with_warning']}")
        for l in s["decisive_lines"] or s["last_lines"]:
            print("     >", l[:300])
    print(f"  load before/after: {rec['before']['loadavg']} / {rec['after']['loadavg']}; "
          f"tracked inputs modified: {rec['before']['tracked_gate_inputs_modified']}; "
          f"checkout changed during gate: {rec['checkout_changed_during_gate']}")
    sys.exit(0)


if __name__ == "__main__":
    main()
