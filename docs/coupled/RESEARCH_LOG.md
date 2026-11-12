# CEC research log: the coordinator's common record

This file records hypotheses, decisions, proof status, benchmarks,
counterexamples, and negative results for the 2026-11 research program.
Entries are append-only. A superseded entry is marked as such, never deleted.

## 2026-11-12: Program start

### Baseline verification, done by the coordinator before any new work
| Check | Result |
|---|---|
| main @ e373a42, `python -m pytest` | 1260 passed (28.8 s) |
| main, `lake build` (Lean 4.28.0, Mathlib v4.28.0) | success, 8041 jobs |
| main, ExecutableTests / LatticeModelTests / cross_language_diff | pass |
| main, axiom audit | 50 declarations, only propext / Classical.choice / Quot.sound |
| integration PR #91 @ 4eecff2, pytest (declared dependencies only) | 2355 passed, 9 skipped; 6 GNM files skipped because `sympy` is undeclared |
| integration PR #91, sympy-gated tests (sympy installed) | 32 passed |
| integration PR #91, `lake build` | 8137/8140 jobs done; MotionRealCrambinGraph.lean compiled for > 30 min (its headline theorem uses `native_decide`) |

### Early audit findings (the full audit is delegated; see docs/audit/)
1. `MotionRealCrambinGraph.frozen_crambin_contact_network_complete` is proved
   by `native_decide`. It therefore trusts the compiler (Lean.ofReduceBool)
   and is NOT kernel-checked, although PR #40 is titled "kernel-check …".
2. `scripts/check_lean_sources.py` scans only `FoldBounds/*.lean`
   non-recursively, does not forbid `native_decide`, and audits only listed
   declarations.
3. The integration branch imports `sympy` (tests) and `numpy` (experiments),
   but `pyproject.toml` declares neither, so tests can be silently skipped.
4. `docs/PRIOR_ART_BASELINE_GATE.md` omits the closest prior art for
   provable partition-function and binding bounds: OSPREY K*/BBK*/MARK* and
   Toulbar2 guaranteed counting.

### Research-direction decision
Primary direction: **coupled ensemble certificates (CEC)**, i.e. rigorous
bounds on relative free energies Z_B/Z_A from shared decompositions.
Reasons:
- It addresses Challenge B directly, and it gives a precise answer to "when
  does coupling help" (Theorem 1.5) together with a sharp limit on local
  certification (Theorems 2.2 and 2.3).
- Its statements are short enough to machine-check completely in this
  session.
- Exact-oracle validation is feasible for HP lattice models, MRFs, and pKa
  models.
- It connects to real structures through multi-site proton binding.

Rejected or deferred alternatives:
- HP-only partition-function improvements, which the repository already
  explores heavily and which have diminishing returns.
- Kinetics / JEPA (secondary per the brief; audit only).
- A Dobrushin comparison-theorem formalization (high value, but large;
  deferred until the core is done).

### Hypotheses (H) and their status
- H1 (T1 coupled enclosure, spec §1): a numerical counterexample search found
  no counterexample in 39,749 exact random instances. Lean: delegated.
- H2 (T2.7 soft separator, spec §2.4): no counterexample in 90,000 random
  mixture-exterior trials. Lean: optional.
- H3: proton-binding couplings are antiferromagnetic in σ for all site-type
  pairs. Verified symbolically (check_af.py).
- H4: local enclosure width is independent of system size for MRFs, while
  separate-bound width grows. Not yet tested (MRF agent).
- H5: for HP mutations, coupled T1 gives only modest gains over the trivial
  ‖V‖∞ bound at matched compute. Not yet tested (HP agent). We expect this
  to be a partially negative result.
- H6: local pKa certificates converge with radius except in strongly coupled
  clusters. Not yet tested (pKa agent).

### Workstreams (isolated worktrees; local branches, not pushed)
| Role | Worktree / branch | Scope |
|---|---|---|
| B: Lean | /home/user/wt-lean, research/cec-lean-20261112 | Formalize T1, T1′, T2, T4; executable checkers |
| C: HP | /home/user/wt-algo, research/cec-algo-20261112 | T1 solver, independent checker, HP mutation benchmarks |
| C2: MRF | /home/user/wt-mrf, research/cec-mrf-20261112 | Separator certificates, scaling, phase transition |
| D: pKa | /home/user/wt-pka, research/cec-pka-20261112 | Real-structure protonation model, pre-registered validation |
| F: literature | docs/coupled/PRIOR_ART.md | Novelty risk per theorem |
| E-audit | docs/audit/ASSURANCE_AUDIT_20261112.md | Deliverable A |
| E: adversary | after the above | Independent attempts to break the results |

### 2026-11-12: Coordinator prior-building experiments (independent of agents)
- T1 counterexample search: `/tmp/claude-0/pi/cex_t1.py`. 39,749 random exact
  instances, 0 failures. Covers T1.1 (exact test, including g₋(λ₋) = 0
  exactly), T1.3 (threshold = vertex brute force), T1.4, and T1.5 (both
  directions).
- T2.7 soft-separator check: `cex_soft.py`. 90,000 mixture-exterior trials,
  0 violations.
- T2.1 / T2.4 enclosure and nesting on random pairwise MRFs (n ≤ 8, exact
  rationals, ball separators): `cex_nest.py`. 688 checks, 0 failures.
- **H5 prior (HP, partially negative):** `hp_ratio_floor.py`, n = 12, q = 2,
  random sequences, H→P at i ∈ {1, n/2, n−2}. For prefix-order regions,
  even with *exact* region masses (no cap looseness), the coupled width from
  ratio uncertainty alone stays near the trivial log q² ≈ 1.386 until the
  prefix depth m ≈ n−2 for interior residues. The width collapses earlier
  only for residues near the chain start. Interpretation: a residue's contact
  environment in a lattice polymer is fixed only by the global fold, so there
  is no local separator. This is consistent with the polymer-vs-MRF dichotomy
  (spec §4). Prefix-frontier T1 gains for HP are expected to be modest.
- Internal prior art: the integration branch (PR #36,
  `docs/COUNTERFACTUAL_FOLDING_REWEIGHTING.md`) already uses
  `Z_mut/Z_parent = E_parent[R]` and conditional-odds identities, with full
  enumeration of 8–10-residue cubic HP models. CEC differs in seeking
  *non-enumerative* enclosures. It must cite that work and not claim the
  identities.
