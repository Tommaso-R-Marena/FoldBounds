#!/usr/bin/env python3
"""Shared helpers of the kernel benchmark harness (standard library only).

* paths of the three prebuilt checkers and of the harness files;
* ``timed_run``: run ONE external command with a hard timeout, precise wall
  time (blocking ``wait4``) and the peak RSS of that child process;
* parsers that normalise the JSON reports of foldkernel / foldcheck /
  foldcheck2 into one shape (eight exact endpoints per (q, t));
* exact-rational helpers (enclosure test, log10 gaps, widths);
* self-avoiding-walk counts (OEIS A001411) and the labelled estimate used
  beyond the table.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import tempfile
import threading
import time
from fractions import Fraction
from pathlib import Path
from typing import Optional

BENCH_DIR = Path(__file__).resolve().parent
REPO = BENCH_DIR.parents[2]
BIN_DIR = REPO / ".lake" / "build" / "bin"
EXE = {
    "kernel": BIN_DIR / "foldkernel",
    "v1": BIN_DIR / "foldcheck",
    "v2": BIN_DIR / "foldcheck2",
}
RESULTS_DIR = BENCH_DIR / "results"
FINDINGS_DIR = BENCH_DIR / "findings"
SCRATCH = Path(os.environ.get("KERNEL_BENCH_SCRATCH", "/tmp/kernel_bench"))
PRODUCER_PATH = REPO / "experiments" / "kernel" / "producer.py"

ENDPOINTS = ("zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi")

# c_k = number of k-step self-avoiding walks on Z^2 (OEIS A001411), k = 0..21.
SAW_COUNTS = (
    1, 4, 12, 36, 100, 284, 780, 2172, 5916, 16268, 44100, 120292, 324932,
    881500, 2374444, 6416596, 17245332, 46466676, 124658732, 335116620,
    897697164, 2408806028,
)


def conformations(n: int) -> tuple[int, bool]:
    """(number of anchored conformations of n residues, is_exact).

    Exact value c_{n-1}/4 inside the table; beyond it the ESTIMATE
    1.177 * 2.638^(n-1) * (n-1)^(11/32) / 4 (flagged ``is_exact = False``).
    """
    k = n - 1
    if 1 <= k < len(SAW_COUNTS):
        return SAW_COUNTS[k] // 4, True
    if k == 0:
        return 1, True
    return int(round(1.177 * 2.638 ** k * k ** (11 / 32) / 4)), False


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frac(s) -> Fraction:
    """Parse ``a/b`` or ``a`` (checker output / manifest strings) exactly."""
    if isinstance(s, Fraction):
        return s
    return Fraction(str(s))


def fstr(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def qkey(q, t) -> str:
    q = Fraction(q)
    return f"{q.numerator}/{q.denominator}|{int(t)}"


def log10_ratio(a: Fraction, b: Fraction) -> Optional[float]:
    """log10(a / b) for positive exact rationals (None when undefined)."""
    if a <= 0 or b <= 0:
        return None
    r = a / b

    def lg(n: int) -> float:
        # math.log10 of huge ints is fine in CPython (arbitrary precision aware)
        return math.log10(n)

    return lg(r.numerator) - lg(r.denominator)


def loadavg() -> list[float]:
    return [round(x, 2) for x in os.getloadavg()]


# ----------------------------------------------------------------------------
# One external run: timeout, wall time, peak RSS
# ----------------------------------------------------------------------------

class _HwmPoller(threading.Thread):
    """Samples VmHWM of a pid with a geometrically growing interval.

    ``ru_maxrss`` returned by wait4 is polluted by the pre-exec image of the
    (much larger) Python parent, so for small children (foldkernel uses a few
    MB) the polled high-water mark is the usable number.  VmHWM is monotone,
    and the interval is at most ~10 % of the elapsed time, so the last sample
    is taken in the final tenth of the run.
    """

    def __init__(self, pid: int):
        super().__init__(daemon=True)
        self.pid = pid
        self.hwm_kb: Optional[int] = None
        self.samples = 0
        self._stop_evt = threading.Event()

    def run(self) -> None:
        path = f"/proc/{self.pid}/status"
        interval = 0.0004
        t0 = time.perf_counter()
        while not self._stop_evt.is_set():
            try:
                with open(path, "rb") as fh:
                    data = fh.read()
                i = data.find(b"VmHWM:")
                if i >= 0:
                    self.hwm_kb = int(data[i + 6:data.index(b"kB", i)].strip())
                    self.samples += 1
            except (OSError, ValueError):
                pass
            elapsed = time.perf_counter() - t0
            interval = min(0.25, max(0.0004, 0.1 * elapsed))
            self._stop_evt.wait(interval)

    def stop(self) -> None:
        self._stop_evt.set()


def timed_run(cmd: list[str], timeout_s: float, cwd: Optional[str] = None,
              env: Optional[dict] = None) -> dict:
    """Run ``cmd`` once.  Never raises on timeout or non-zero exit.

    Returns wall_s, timed_out, returncode, stdout, stderr, peak RSS (KB) and
    how it was obtained, user/sys CPU seconds of the child.
    """
    out_f = tempfile.TemporaryFile()
    err_f = tempfile.TemporaryFile()
    self_rss = _self_maxrss_kb()
    state = {"timed_out": False}
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, stdin=subprocess.DEVNULL, stdout=out_f, stderr=err_f,
                            cwd=cwd, env=env, start_new_session=True)

    def kill() -> None:
        state["timed_out"] = True
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            pass

    timer = threading.Timer(timeout_s, kill)
    timer.daemon = True
    poller = _HwmPoller(proc.pid)
    timer.start()
    poller.start()
    try:
        _, status, ru = os.wait4(proc.pid, 0)
    finally:
        t1 = time.perf_counter()
        timer.cancel()
        poller.stop()
    poller.join(timeout=1.0)
    proc.returncode = os.waitstatus_to_exitcode(status)
    out_f.seek(0)
    err_f.seek(0)
    stdout = out_f.read().decode("utf-8", "replace")
    stderr = err_f.read().decode("utf-8", "replace")
    out_f.close()
    err_f.close()
    wait4_kb = int(ru.ru_maxrss)
    if wait4_kb > self_rss:
        rss_kb, rss_src = wait4_kb, "wait4"          # child grew beyond the parent's image: exact
    elif poller.hwm_kb is not None:
        rss_kb, rss_src = poller.hwm_kb, "VmHWM-poll"
    else:
        rss_kb, rss_src = None, "unavailable"
    return {
        "wall_s": t1 - t0,
        "timed_out": bool(state["timed_out"] and proc.returncode == -signal.SIGKILL),
        "returncode": proc.returncode,
        "stdout": stdout,
        "stderr": stderr[-2000:],
        "rss_kb": rss_kb,
        "rss_source": rss_src,
        "rss_wait4_kb": wait4_kb,
        "rss_polled_kb": poller.hwm_kb,
        "rss_poll_samples": poller.samples,
        "cpu_user_s": ru.ru_utime,
        "cpu_sys_s": ru.ru_stime,
    }


def _self_maxrss_kb() -> int:
    import resource
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


# ----------------------------------------------------------------------------
# Checker reports
# ----------------------------------------------------------------------------

def parse_report(checker: str, stdout: str) -> dict:
    """Normalise one checker JSON line.

    Returns {"accepted": bool|None, "verify_s": float|None, "results":
    {qkey: {endpoint: "num/den"}}, "raw": {...small echo fields...}, "error"}.
    """
    lines = [ln for ln in stdout.strip().splitlines() if ln.strip()]
    if not lines:
        return {"accepted": None, "error": "no output", "results": {}, "verify_s": None, "raw": {}}
    try:
        obj = json.loads(lines[-1])
    except ValueError as exc:
        return {"accepted": None, "error": f"unparsable output: {exc}", "results": {},
                "verify_s": None, "raw": {}}
    rep: dict = {"accepted": None, "error": None, "results": {}, "verify_s": None, "raw": {}}
    if "parse_error" in obj:
        rep["accepted"] = False
        rep["error"] = f"parse_error: {obj['parse_error']}"
        return rep
    rep["accepted"] = bool(obj.get("accepted"))
    if checker == "kernel":
        if "verify_us" in obj:
            rep["verify_s"] = obj["verify_us"] / 1e6
        elif "verify_ms" in obj:
            rep["verify_s"] = obj["verify_ms"] / 1e3
        for key in ("sequence", "n", "cert_bytes", "split_nodes", "bound_leaves", "enum_leaves",
                    "exact_leaves", "exact_count", "exact_histogram", "reason", "pos"):
            if key in obj:
                rep["raw"][key] = obj[key]
        for r in obj.get("results", []):
            rep["results"][qkey(frac(r["q"]), r["threshold"])] = {e: r[e] for e in ENDPOINTS}
    else:
        if "verify_ms" in obj:
            rep["verify_s"] = obj["verify_ms"] / 1e3
        for key in ("sequence", "q", "threshold", "nodes", "version"):
            if key in obj:
                rep["raw"][key] = obj[key]
        if rep["accepted"] and all(e in obj for e in ENDPOINTS):
            rep["results"][qkey(frac(obj["q"]), obj["threshold"])] = {e: obj[e] for e in ENDPOINTS}
    if not rep["accepted"] and rep["error"] is None:
        rep["error"] = "rejected" + (f": {obj.get('reason')}" if obj.get("reason") else "")
    return rep


def encloses(end: dict, ref: dict) -> bool:
    """Do the eight endpoints enclose the exact reference values?"""
    g = {e: frac(end[e]) for e in ENDPOINTS}
    z, za, zn, p = (frac(ref[k]) for k in ("Z", "Z_A", "Z_notA", "P"))
    return (g["zLo"] <= z <= g["zHi"] and g["eLo"] <= za <= g["eHi"]
            and g["nLo"] <= zn <= g["nHi"] and g["pLo"] <= p <= g["pHi"])


def is_exact(end: dict) -> bool:
    return all(end[a] == end[b] for a, b in (("zLo", "zHi"), ("eLo", "eHi"),
                                             ("nLo", "nHi"), ("pLo", "pHi")))


def equals_reference(end: dict, ref: dict) -> bool:
    g = {e: frac(end[e]) for e in ENDPOINTS}
    z, za, zn, p = (frac(ref[k]) for k in ("Z", "Z_A", "Z_notA", "P"))
    return (g["zLo"] == z == g["zHi"] and g["eLo"] == za == g["eHi"]
            and g["nLo"] == zn == g["nHi"] and g["pLo"] == p == g["pHi"])


def accuracy_fields(end: dict, ref: Optional[dict]) -> dict:
    """Width and (when a reference exists) log-gap fields of one enclosure."""
    plo, phi = frac(end["pLo"]), frac(end["pHi"])
    zlo, zhi = frac(end["zLo"]), frac(end["zHi"])
    width = phi - plo
    out = {
        "p_width": fstr(width),
        "p_width_float": float(width),
        "log10_zHi_over_zLo": log10_ratio(zhi, zlo),
    }
    if ref is not None:
        z = frac(ref["Z"])
        out["log10_zHi_over_Z"] = log10_ratio(zhi, z)
        out["log10_Z_over_zLo"] = log10_ratio(z, zlo)
        out["encloses_reference"] = encloses(end, ref)
    else:
        out["log10_zHi_over_Z"] = None
        out["log10_Z_over_zLo"] = None
        out["encloses_reference"] = None
    return out


def load_manifest(name: str) -> dict:
    with open(BENCH_DIR / f"manifest_{name}.json", "r", encoding="ascii") as fh:
        return json.load(fh)


def write_scratch(name: str, data: bytes) -> Path:
    SCRATCH.mkdir(parents=True, exist_ok=True)
    path = SCRATCH / name
    with open(path, "wb") as fh:
        fh.write(data)
    return path


if __name__ == "__main__":
    print(json.dumps({k: str(v) for k, v in EXE.items()}, indent=1))
    print("conformations n=22:", conformations(22), " n=30:", conformations(30))
    print(timed_run(["/bin/true"], 5.0))
    sys.exit(0)
