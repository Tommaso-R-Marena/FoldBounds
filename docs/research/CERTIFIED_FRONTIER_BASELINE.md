# Baseline before the frontier kernel

Recorded so later speed/accuracy claims stay comparable.

## Platform

* Lean 4.28.0, Mathlib v4.28.0, Lake 5.0.0
* Python 3.12, pytest via `.venv`
* Host: 4 CPUs, ~15 GB cgroup, Linux 6.12

## Prior checkers (PR #32 head)

| Binary | Size | Notes |
|---|---|---|
| `foldcheck` | ~222 MB | Mathlib-linked, one `(q,t)` per run |
| `foldcheck2` | ~222 MB | Aristotle v2 leaf families; same size class |

Build cost (cached Mathlib): `lake build foldcheck foldcheck2` ≈ 8 min wall /
25 CPU-min on this host.

## Gaps closed by the kernel

* Compiled checkers were specification-style and slow; linked Mathlib.
* Text parser sat outside every theorem.
* One `(q, threshold)` per check.
* Complement-mass caps ignored the threshold (`q^K` vs `q^min(K,t-1)`).

## Kernel response

* Mathlib-free executable module + `foldkernel` (~5 MB).
* String-level `verify_sound` with parser inside the theorem.
* Multi-query evaluation of one accepted tree.
* Threshold-aware complement cap; grid enumerator proved = list enumerator.
