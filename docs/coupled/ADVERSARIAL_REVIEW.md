# Adversarial review of the CEC work (Agent E)

**Status: phase 1 complete (2026-11-13).** This phase covers the Lean formalization
(`/home/user/wt-lean` @ d8af789), the MRF work (`/home/user/wt-mrf` @ 11228b1),
`CEC_THEORY.md`, `PRIOR_ART.md`, `RESEARCH_LOG.md` and the audit. Phase 2 (HP algorithms and pKa)
will extend this file.

The reviewer did not write any of the reviewed work. No worktree was modified. All scripts are
under `/tmp/claude-0/adversary/`. The Lean work used two private copies made with `git archive`
(`lean-copy/`, and `gatecopy/` for the gate demonstration), with the shared `.lake/packages`
symlinked and never updated.

## Summary

- **No CRITICAL finding.** No checker accepted an unsound certificate in any test:
  - the Lean `checkEnclosure`/`checkSeparator` on 716 cross-checked inputs plus kernel spot checks;
  - the MRF independent checker on 50 crafted attacks and 6,241 random degenerate certificates.
- **Three MAJOR findings:**
  - (M1) both Lean axiom gates can be bypassed with a guillemet-quoted axiom name;
  - (M2) the headline MRF scaling benchmark cannot carry the weight put on it;
  - (M3) three unformalized spec items (§2.3 twice, and Theorem 4.4) are false in degenerate
    cases, by the same zero-weight mechanism that broke Theorems 2.4 and 4.2.
- **The reported numbers reproduce.**
  - Lean: the build succeeds, and an independent audit of all 225 Coupled constants by `Name`
    finds only the standard axioms.
  - MRF: every scaling and phase-table number I recomputed matches to 12 digits or better.

## Reproduction (independent)

| What | Command (from `/tmp/claude-0/adversary`) | Result |
|---|---|---|
| Build | `cd lean-copy && lake build` | OK, 3 m 25 s, 8045 jobs, Coupled modules, tests and audit included |
| Gates on own log | `lake env lean FoldBounds/AxiomAudit.lean > ../ax_lean.log`, then both `scripts/check_*` | 143 / 93 declarations, standard axioms only |
| Own axiom audit, all constants | `cd lean-copy && lake env lean Adv/AllAxioms.lean` | 225 constants in `CoupledRatio/Separator/Binding` (176 theorems), **0 non-standard axioms**; the only `partial` is the compiler-generated `ratioBounds._unsafe_rec` |
| Token grep | `grep -rnE 'sorry\|admit\|native_decide\|axiom\|implemented_by\|extern\|unsafe\|opaque' FoldBounds/` | none outside `#print axioms` lines |
| T1 / separator cross-check vs Lean `#eval` | `python3 t1/gen_lean.py; (cd lean-copy && lake env lean Adv/T1Cross.lean > ../t1/lean_out.txt); python3 t1/compare.py` | 412 enclosure inputs (398 accepted at my brute-force endpoints and rejected when 10⁻²⁰⁰ tighter; 9 with empty Adm and 5 malformed give `none`) and 304 separator inputs: **0 discrepancies** |
| Kernel agrees with `#eval` | `cd lean-copy && lake env lean Adv/KernelCross.lean` | 3 large-rational accept/reject pairs by `decide +kernel`; standard axioms |
| T1 degenerate sweep | `python3 t1/spec_sweep.py` | 38,886 instances, many of them degenerate: threshold form, `g(λ)=0`, T1.4 and T1.5 have 0 failures |
| MRF checker attacks | `cd mrf && PYTHONPATH=/home/user/wt-mrf/src python3 -B attacks.py` | 50 attacks, **no unsound acceptance** |
| MRF degenerate random | `... python3 -B random_degenerate.py <seed> <n>`, seeds 1–4 | 1,650 models, 6,241 certificates: producer = my brute force; exact ratio inside every enclosure with Z_A > 0; 0 valid certificates rejected; 0 claims tightened by 10⁻⁶⁰ accepted |
| MRF numbers | `cd mrf && python3 -B ising_check.py` | scaling widths, exact ratio and phase table R = 1..6 reproduced; FKG spot check: 360 random boundaries at R = 5, 6, 0 violations |
| T2.7 degenerate | `python3 soft/t27_degenerate.py` | 28,565 exteriors with zero Φ entries, Φ_A = 0 < Φ_B and ternary variables: 0 violations; necessity counterexample reproduced (2/3 ∉ [4/9, 5/9]) |
| Gate bypass | `cd gatecopy && python3 scripts/check_lean_sources.py ../gate/evil_ax.log; python3 scripts/check_coupled_lean.py ../gate/evil_ax.log` | **both exit 0** on a build that proves `(0:ℝ) = 1` (M1) |
| §2.3 / T4.4 counterexamples | `cd lean-copy && lake env lean Adv/Degenerate.lean` | `some (1, 1)`, but the true value is 2 (M3) |

## Findings

Where to reproduce each finding (paths relative to `/tmp/claude-0/adversary`):

| Finding | Reproduction |
|---|---|
| M1 | the gate-bypass row of the table above (`gatecopy/FoldBounds/AdvEvil.lean`) |
| M2 | `python3 -B mrf/ising_check.py`, and the `exact` column of `wt-mrf/experiments/coupled_mrf/results/scaling_ising_field.csv` |
| M3 | `lean-copy/Adv/Degenerate.lean` |
| m1 | `python3 soft/t28_sign.py` |
| m2 | `python3 t1/spec_sweep.py` |
| m3, m4 | `mrf/attacks.py` |
| m5 | the signature of `affinity_mem_Icc` in `CoupledBinding.lean` |
| m6–m8 | the cited lines |

### CRITICAL

None found. The executable checkers are sound on every input tried. The inputs included:
- zero, negative and malformed data;
- rationals with 30-digit and 300-digit parts;
- claims tightened by 10⁻²⁰⁰ (Lean) or 10⁻⁶⁰ to 10⁻¹⁰⁰ (MRF);
- floats and the tokens NaN, Infinity and `1/0`;
- 5000-digit literals;
- truncated tables and models;
- separators that break S1, S2 or S3;
- reordered or omitted boundary variables;
- witnesses that are wrong or truncated.

### MAJOR

**M1. Both axiom gates accept a proof of `0 = 1` (gate soundness).** In the scratch copy,
`FoldBounds/AdvEvil.lean` contains `axiom «Classical.choice» : False` and
`theorem FoldBounds.CoupledRatio.adv_false : (0:ℝ) = 1`. Lean reports
`[«Classical.choice», propext, Classical.choice, Quot.sound]`.
- `check_lean_sources.py` misses the declaration, because its regex `\baxiom\s+\w+` does not
  match `«`.
- Both gates tokenise the report with `[A-Za-z_][A-Za-z0-9_.]*`, which strips the guillemets
  and reads the name as `Classical.choice`.
- `check_coupled_lean.py` scans source only in `Coupled*.lean`.
- Result: the old gate prints "144 declarations, standard axioms only" and the new gate prints
  "94 declarations, standard axioms only".

This refutes `LEAN_FORMALIZATION.md`'s description of `check_coupled_lean.py` as a "fail-closed
source and axiom audit". It is a second gate flaw, independent of F3e (which I also reproduced).
It is **not exploited** in the real tree, as my `Name`-level audit shows.

Fix: compare axiom `Name`s inside Lean (`Lean.collectAxioms`, checking membership of
``propext``, ``Classical.choice``, ``Quot.sound``), reject `«` in reports, and use `\baxiom\b`
in every scanned module.

**M2. The headline MRF scaling experiment is uninformative as evidence** (benchmark design and
significance).
- (a) The instance has no bulk field, so spin-flip symmetry gives `Z_B/Z_A = 3/4` exactly for
  every L. I checked this exactly at L = 4 and to float precision at L = 5. All 44 "exact"
  values in `scaling_ising_field.csv` are 0.75. The local data has the same symmetry, so every
  enclosure is symmetric about 3/4 (lower + upper = 3/2 for each R). The check "exact value
  inside all 44 enclosures" therefore tests almost nothing.
- (b) Size-independence holds by construction. The certificate is a function only of the
  factors inside the radius-R ball, and these are identical for every L ≥ 2R+1. "Identical to
  16 digits for every L ≥ 7" restates the definition.
  - It is still presented as a result: "H4 SUPPORTED" (`RESEARCH_LOG.md` l.121), "the measured
    size-independence" as a candidate contribution (`CEC_THEORY.md` §5), and "empirical
    demonstration of size-independence" (`PRIOR_ART.md` §11.4).
- (c) The separate baseline is a strawman.
  - Zero-field planar Ising is exactly solvable in polynomial time (Kasteleyn/FKT, Barahona).
    The single-site field reduces to it with one pendant "ghost" spin, because Z(h) = Z_ghost/2
    by symmetry.
  - The agents' own transfer-matrix oracle computed the exact ratio, a separate method with
    width 0, for all L ≤ 18. Those are the sizes at which the table calls the separate bound
    vacuous.
  - The caveat in `MRF_SCALING.md` names Weitz and TRW but not this.
  - In the ordered phase, the hard enclosure of a ferromagnet is the classical ± boundary
    "sandwich". The zero-field rows of experiment (ii) also have exact ratio 3/4 at every β.

  The numbers themselves are correct (reproduced). What fails is the inference drawn from them.

**M3. Three unformalized spec reductions are false in degenerate cases.** The mechanism is the
one that broke Theorems 2.4 and 4.2: a state with zero A-weight but positive B-weight.
- (a) **§2.3, interval exterior** ("apply T1 with α_τ ∈ [Φ_A ℓ, Φ_A u] and exact ratio
  ρ(τ)"). Take Φ_A = (1, 0), Φ_B = (1, 1), Ψ ≡ 1. The true value is Z_B/Z_A = 2. The reduction
  gives a region with α ∈ [0, 0], and the Lean checker accepts [1, 1] for any finite stand-in
  ratio. Required: `Φ_A(τ) = 0 ⇒ Φ_B(τ)·u(τ) = 0`; otherwise the upper end is +∞.
- (b) **Theorem 4.4**, "certified by T1′ duals" through α = p·k₁. Take k₁ = (1, 0),
  k₂ = (1, 1) and p uniform. The true selectivity is S = 2; the reduction is accepted with
  [1, 1]. Required: k₁ > 0 on every state where some p ∈ Π is positive. This is the corrected
  4.2 hypothesis, which was not propagated to 4.4.
- (c) **§2.3, soft separator via T1′.** The same α-encoding forces β_τ = 0 when Φ_A(τ) = 0 < Φ_B(τ).
  The fix is to state the linear-fractional program in Ψ directly.

These items may be used in phase 2 (pKa, selectivity).

### MINOR

- **m1. T2.8 identity has a sign error.** For s < 0, `W_R = (e^{−s} − 1)·osc μ`, but the spec's
  `(1 − e^{−s})·osc μ` is negative. With one spin, q = 3 and s = −1, the width is 0.859 and the
  formula gives −0.859. It needs `|1 − e^{−s}|`, with the oscillation taken over feasible τ
  (T⁺). As `PRIOR_ART` notes, "Theorem 2.8" restates the definition of weak spatial mixing; it
  should be called a corollary.
- **m2. The T1.3 sorted-threshold form must skip thresholds with Σα = 0.** With that rule it
  holds: 0 failures in 38,886 degenerate instances. A naive `0/0 → 0` rule gives a sound but
  non-sharp λ₋ in 1,990 of them. This matters for phase-2 implementations.
- **m3. The MRF checker silently ignores unknown model-level and factor-level keys, which the
  producer rejects.** A typo key `only_b` holding a B-changing factor leaves the hash unchanged,
  and the old certificate is accepted. `MRF_SCALING.md` says the checker "rejects unknown keys";
  that is true only of certificate keys.
- **m4. The MRF checker returns `ok=True` when Φ_A ≡ 0, for any claim that has no `lower`.** An
  example is `{upper: "0", tight: true}`. The acceptance is vacuous, since Z_A = 0 for every
  exterior, but the `ok` is misleading.
- **m5. Lean Theorem 4.1 relies on a junk value when Z_L = 0.** `affinity_mem_Icc` and
  `affinity_sharp`, and the k₂ side of 4.2, hold because Lean sets K_L = 0 when Z_L = 0; the spec
  value is +∞. No positivity hypothesis is added and the case is not flagged, contrary to the
  doc's own stated policy.
- **m6. `CEC_THEORY.md` overstates what is formalized.** Its formalization status calls §2.1–2.2
  machine-checked. It omits two unformalized parts that `LEAN_FORMALIZATION.md` does list:
  realisability by pairwise exteriors (T2.2) and multi-layer products (T2.6).
- **m7. Overclaim sentences:**
  - `MRF_SCALING.md` l.106: "vacuous at the sizes where the coupled bound is still sharp". This
    conflates instance-optimal with tight; the width is 0.106 at R = 3.
  - `CEC_THEORY.md` l.334: "This is why mutational stability ΔΔG … is not a local quantity".
    This is untested, yet stated as fact.
  - `CEC_THEORY.md` l.393: "Allostery is the obstruction …". This is a categorical sentence
    inside a paragraph labelled as a hypothesis.
  - `RESEARCH_LOG.md` l.121: "H4 SUPPORTED" (see M2).
- **m8. The contraction result is reported more broadly than computed.** `MRF_SCALING.md` says δ
  was computed "for inner R = 1, 2". `contraction.csv` has δ only for the 1→2 annulus; the
  inner R = 2 rows are blank. "The bound held in every computed case" therefore means 3 cases.

### NOTE

- **n1. Lean statements match the doc.** For every headline item checked:
  - the quantifiers range over all admissible (α, β), over all nonnegative exteriors on any
    finite O, and over interiors with the given marginals;
  - the hypotheses are satisfiable (the tests instantiate them);
  - no soundness theorem refers to a different function from the one that is evaluated.

  Separator soundness holds relative to the listed boundary states. The completeness of that
  enumeration is not verified, as documented.
- **n2. The FKG shortcut is stated accurately and never generalised** (`MRF_SCALING.md` §2 and
  §4, `RESEARCH_LOG`). Its "nonnegative field" condition is stronger than needed. The R = 5, 6
  values were reproduced independently.
- **n3. T2.7 survives a degenerate sweep.** Its proof needs only Φ_A ≢ 0, which is stronger than
  stated. The hypothesis is necessary.
- **n4. T1′ completeness needs only P ≠ ∅.** Boundedness is unnecessary, because λ ≤ Q bounds the
  primal below by 0.
- **n5. Duplicate JSON keys are not detected by the checker.** It takes parsed dicts, so
  duplicates are resolved by the caller's parser (Python keeps the last value). A text entry
  point that rejects duplicates is advisable.
- **n6. The MRF experiments have no Monte Carlo baseline.** The question of a fair Monte Carlo
  baseline therefore does not arise in phase 1.

## Verdict table

| Claim | Correct | Original | Scalable | Significant |
|---|---|---|---|---|
| T1 exact tests and rational checker (Lean) | supported | unsupported (Dinkelbach, Charnes–Cooper) | supported (linear-time check) | partially (assurance kernel) |
| T1.3 threshold form (unformalized) | partially (needs Σα > 0 rule) | unsupported | supported | partially |
| T2.1–2.6, 2.2′ separator theorems (Lean) | supported (2.4 as corrected) | partially (framing of 2.2, 2.2′, 2.3); unsupported for the rest | partially (exponential in the boundary) | partially |
| Lean `checkSeparator` and MRF checker soundness | supported | partially | partially (≤ 2¹⁸ states) | partially |
| T2.7 soft separator | supported (numerical; unformalized) | partially | unsupported (2^‖∂‖ to 3^‖∂‖ candidates) | unsupported (adds nothing for ferromagnetic Ising, by the agents' own data) |
| T2.8 certifiability ⇔ WSM | partially (sign error, m1) | unsupported (definition of WSM) | — | partially (interpretive) |
| §2.3 reductions; T4.4 | **refuted as stated** (M3) | unsupported | supported (LP) | partially |
| T4.1–4.3 (Lean) | supported (k₁ > 0; caveat m5) | partially | — | partially (model-relative) |
| H4 size-independence (MRF) | supported (trivially) | unsupported | supported | **unsupported** (M2) |
| MRF phase / T2.8 numbers | supported (reproduced) | unsupported (classical sandwich and WSM) | partially (FKG holds only for ferromagnets) | partially |
| Lean assurance claims | axioms: supported; "fail-closed gate": **refuted** (M1) | — | — | supported |

## Notes for phase 2

- **T1 producers.** Check the support condition: R_k = +∞ cannot be expressed in the checker.
  Check the Σα = 0 handling of the threshold form (m2) and the degenerate cases of M3 if pKa or
  selectivity uses §2.3 or 4.4.
- **The HP log (2026-11-13)** already concedes that the coupled-versus-separate gap is an
  artifact of the baseline. Verify it with an independent enumerator for n ≤ 13 on at least 30
  rows.

---

## Coordinator response to phase 1 (2026-11-13)

Every finding was accepted. Dispositions:

| Finding | Disposition | Where |
|---|---|---|
| M1: forged `«Classical.choice»` axiom passes both gates | **Fixed in the new CEC gate.** Axiom entries are parsed strictly; the reviewer's evil log and a synthetic wrapped forgery are now rejected. **Recorded for the existing gate** (audit F3f), with a fix proposed; owner approval needed. | `scripts/check_coupled_lean.py`; audit F3f |
| M2: MRF scaling headline (symmetric exact value; size independence by construction; strawman baseline) | **Accepted and reframed.** MRF_SCALING.md now states that size independence follows from the construction, that the tree baseline is a strawman for zero-field planar Ising, and that the meaningful content is width-at-radius, phase dependence, and checked soundness. New run `run_scaling_field.py` with a bulk field (W = 6/5): the exact ratio varies with L, and is inside all 36 enclosures, 0 violations. | MRF_SCALING.md; `results/scaling_ising_bulkfield.csv` |
| M3: §2.3 reductions and Theorem 4.4 fail with zero A-weight | **Accepted.** The spec now requires the support condition (§2.3) and k₁ > 0 (Theorem 4.4), each with the counterexample. They remain unformalized. | CEC_THEORY.md |
| m: T2.8 sign for s < 0 | Fixed (absolute value) | CEC_THEORY.md |
| m: T1.3 threshold rule must skip Σα = 0 | Fixed in the spec. Phase 2 checks the implementations. | CEC_THEORY.md |
| m: MRF checker ignored unknown keys; accepted an upper claim when Φ_A ≡ 0 | **Fixed.** Unknown model and factor keys are rejected, and no lower or upper claim is accepted on empty support. Two regression tests added; 58 MRF tests pass. | `coupled_mrf/checker.py`; tests |
| m: Theorem 4.1 with Z_L = 0 holds only via x/0 = 0 | Flagged in LEAN_FORMALIZATION.md | LEAN_FORMALIZATION.md |
| m: formalization status incomplete (T2.2 pairwise realisability, T2.6 multi-layer) | Fixed | CEC_THEORY.md |
| m: overclaims (MRF "still sharp"; ΔΔG "This is why"; "Allostery is the obstruction") | Softened to "informative", "we conjecture", and "untested hypothesis" | MRF_SCALING.md; CEC_THEORY.md |
| m: Dobrushin coefficient for a single annulus only | Fixed wording | MRF_SCALING.md |
