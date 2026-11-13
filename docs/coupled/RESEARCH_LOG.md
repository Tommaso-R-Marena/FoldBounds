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

## 2026-11-13: Session resumed after a usage-limit stop
- At about 22:00 UTC all six agents stopped on the account's 5-hour usage
  limit. Partial work was saved on disk. The limit reset at 02:10 UTC.
- Resource decision: run at most four concurrent agents. The Lean agent was
  resumed at the strongest configuration, as the brief asks for difficult
  proofs. The implementation, experiment, and literature agents were
  relaunched on a lower-cost model, per the brief's guidance on routine work.
- The local integration `lake build` was killed by my 90-minute timeout
  (exit 124) while still compiling MotionRealCrambinGraph. GitHub CI on #91
  shows every Lean job cancelled at about 45 minutes. See the audit.
- Deliverable A (the audit) and the prior-art review are committed. Prior-art
  verdict: most CEC mathematics is classical (high novelty risk). The
  defensible parts are the sharpness/separation framing, T2.7, the
  machine-checked kernel, and certified model titration. OSPREY combines
  per-state approximations into K*, per the MARK* abstract (provisional).

### MRF results (agent C2; branch research/cec-mrf-20261112 @ 11228b1; not yet merged)
- H4 (size independence) is SUPPORTED. Ising grid, q = 3/2, centre-field
  perturbation [1, 1/2]:
  - The coupled width is identical to 16 digits for every L ≥ 7 at each
    radius R: 0.5, 0.3351, 0.2052, and 0.1059 for R = 0..3.
  - The separate spanning-tree bound grows like L²: log10(UB/LB) goes from
    5.6 at L = 5 to 5680 at L = 128.
  - The exact value lies inside all 44 checked enclosures for L ≤ 18.
  - Coordinator hand check: R = 1 width = ½·(q⁴−1)/(q⁴+1) = 0.33505. ✓
- T2.8 prediction is SUPPORTED. Width at R = 10 for β/β_c =
  0.5 / 0.8 / 1.0 / 1.25 / 2.0 is 0.0007 / 0.099 / 0.384 / 0.477 / 0.499.
  It decays above T_c and plateaus below. R > 4 uses the FKG all-0/all-1
  extremes, which are valid only for ferromagnetic Ising and are not
  verified by the minimal checker.
- Sharpness (T2.2): constructed exteriors realise the endpoints exactly
  (12/12).
- T2.7 soft separator: 0 violations in 48 exact cases. The hypothesis is
  shown NECESSARY: an exterior factor touching two boundary variables gives an
  exact ratio of 2/3, outside the soft enclosure [4/9, 5/9].
- Negative findings:
  - The soft separator with the true sensitivities adds nothing over the hard
    separator at R+1 for ferromagnetic Ising.
  - The Dobrushin coefficient is nearly vacuous (0.99) at and below T_c.
  - Exhaustive enclosures are limited to 2^18 boundary states, i.e. R ≤ 4 in
    2D.

### Lean results (agent B; branch research/cec-lean-20261112 @ d8af789; not yet merged)
- Coordinator re-verification:
  - `lake build` OK.
  - `CoupledTests.lean` OK (kernel evaluation via `decide +kernel`).
  - Existing gate: 18 modules / 143 declarations, standard axioms.
  - New gate: 93 Coupled declarations, standard axioms.
  - No sorry, admit, or native_decide.
- Formalized:
  - All of §1: T1.1/1.2 iff, T1.3 attainment with g(λ) = 0, T1.4,
    T1.5 iff, T1.6, T1.7(a,b), T1′ soundness.
  - The exact rational checker, with exactness theorems.
  - All of §2.1–2.2: separation (the achievable set is exactly
    (0,∞) × [min ρ, max ρ]), instance optimality, hierarchy, termination,
    contraction, and the exact separator checker.
  - §3: T4.1–4.3.
- **Spec errors found, now corrected in CEC_THEORY.md:**
  - Theorem 2.4 needs the inner support condition (formal counterexample).
  - Theorem 4.2 needs k₁ > 0 everywhere (formal counterexample).
  - Lesson: cex_nest.py used strictly positive weights and could not hit the
    support case. Counterexample searches must include degenerate (zero)
    weights.
- **New audit finding (F3e):** the existing axiom gate ignores continuation
  lines of wrapped `#print axioms` reports. The coordinator reproduced it: a
  forged log hiding `sorryAx` passes. A fix is proposed in the audit, not
  applied, pending owner approval.

### HP results (agent C; branch research/cec-algo-20261112 @ 5303e02; not yet merged)
- Coordinator check: 14 tests pass (14 s). The summary was regenerated from
  the CSVs by summarize.py.
- **H5 is CONFIRMED: a negative or modest result.** T1 is sound (0 failures
  over 71,859 oracle rows; the independent checker accepted 4,785 + 2,160
  certificates and rejected none). On 2D HP, however, the coupled ΔΔF interval
  ≈ the trivial ‖V‖∞ bound at matched frontiers for n ≥ 12 with ≤ 3000
  regions (median coupled/trivial = 0.96–1.00). At n = 30 and 50 it is
  *exactly* trivial in every row.
- Gains over trivial appear only when a large share of the tree is explored:
  n = 10 at ≤ 3000 regions (74% of leaves) gives median width 0.14 kT vs
  trivial 1.10. The coupled bound reaches 0.1 kT at a median 51–59% of the
  tree; the separate bound never does.
- The ~10× coupled/separate advantage is mostly an artifact of the separate
  baseline paying loose 3^r·q^Ucap caps twice. It must NOT be presented as the
  headline. The meaningful baseline is the trivial bound.
- Interpretation: consistent with the polymer no-separator prior. A residue's
  contact environment is fixed only by the global fold. CEC adds a
  machine-checkable combination rule for HP, not a new capability. The
  centred growth order gave inconsistent differences, and no claim is made.

## 2026-11-13: Second usage-limit stop; adversarial review phase 1
- The account limit was hit again at about 03:30 UTC (reset 07:20). The
  max-effort adversarial agent was the main consumer. On resumption its scope
  was narrowed to conclusion-changing checks.
- **Review phase 1: 0 CRITICAL, 3 MAJOR, 8 MINOR. All accepted:**
  - M1, gate forgery: fixed in the new gate; recorded as audit F3f for the
    old gate.
  - M2, MRF scaling: reframed. H4 is a construction-level corollary, not
    evidence. A bulk-field soundness rerun gave 0 violations in 36 cases.
  - M3, §2.3 / Theorem 4.4 support conditions: corrected.
  - The minor findings are fixed; see the ADVERSARIAL_REVIEW.md coordinator
    response.
- **H4 status revised:** from "SUPPORTED" to "holds by construction; not
  evidence of a capability". The scientifically meaningful MRF results are
  the per-radius widths and the T2.8 phase dependence.
- The Lean and MRF branches were merged into the dev branch after phase 1,
  with the fixes applied on top.
