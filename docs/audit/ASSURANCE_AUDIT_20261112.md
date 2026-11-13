# FoldBounds research and assurance audit (Deliverable A)

**Auditor:** research coordinator. An independent audit agent was started and
then stopped by a usage limit before it wrote anything. Every finding below
was reproduced by the coordinator with the commands shown.
**Dates:** 2026-11-12/13.
**Scope:** `main` @ e373a42; integration PR #91 @ 4eecff2; the JEPA stack
#92–#95 (tip @ 1a30e28); selected source PRs.
**Toolchain:** Lean 4.28.0, Mathlib v4.28.0 (pinned), Python 3.13.18, on a
4-core, 15 GB container.

## Executive summary

1. **`main` is sound and green, within its stated scope.**
   - All 1260 Python tests pass.
   - `lake build` succeeds.
   - The executable-checker tests, the lattice tests, and the Python–Lean
     differential tests pass.
   - The 50 audited declarations use only `propext`, `Classical.choice`, and
     `Quot.sound`.
   - Its documentation describes the Python-to-Lean gap correctly.
2. **The integration branch's (#91) Lean library has never been verified as a
   whole.**
   - On GitHub, **every Lean job on #91 was cancelled at the ~45-minute
     timeout**: 17 Lean jobs among the 56 check runs inspected, while every
     Python job passed.
   - Locally, `lake build` had not finished after **90 minutes** (exit 124);
     `FoldBounds/MotionRealCrambinGraph.lean` was still compiling.
   - Formal claims made *about the integrated tree* are therefore unverified.
     Several component PRs passed their own Lean CI separately; for example,
     Aristotle v2 (#32) `lean-proof` succeeded on 2026-10-20.
3. **`native_decide` is used in a headline library theorem whose PR claims
   kernel checking.**
   - `MotionRealCrambinGraph.frozen_crambin_contact_network_complete` is proved
     by `native_decide`. Its trust base therefore includes the Lean compiler
     and the `Lean.ofReduceBool` axiom.
   - PR #40 is titled "Experiment: kernel-check full 46-residue Crambin contact
     completeness". That is a factual overstatement: the result is
     compiler-checked, not kernel-checked.
   - About 25 integration *test* files also use `native_decide`. That is
     acceptable for tests, but those files must not be cited as proofs.
4. **The proof-hygiene gate has blind spots.** `scripts/check_lean_sources.py`:
   - (a) globs `FoldBounds/*.lean` non-recursively;
   - (b) does not forbid `native_decide`, so `Lean.ofReduceBool` is caught only
     for declarations that happen to be listed in `AxiomAudit.lean`;
   - (c) audits only listed declarations;
   - (d) cannot see Lean files outside `FoldBounds/`, e.g. `scripts/foldcheck*.lean`
     on #91.
5. **Undeclared dependencies silently skip tests.**
   - On #91, six GNM test files `importorskip("sympy")`, but `pyproject.toml`
     declares only `pytest`. In a fresh install: 2355 passed, 9 skipped.
     Installing `sympy` adds 32 passing tests.
   - `numpy` is used by experiments without being declared.
   - The JEPA stack similarly skips its `torch`/`mdtraj` tests, which is
     documented as optional.
6. **Overclaiming is concentrated in PR titles; the documents themselves are
   mostly careful.** Most research documents carry explicit status disclaimers
   ("not a predictor of experimental protein structure", "not Lean-kernel
   checked"). Some titles overstate generality or formal status:
   - #40: "kernel-check" (see 3).
   - #70: "1,200-residue nonenumerative kinetic no-go". This is one
     *constructed* all-H 4×h Hamiltonian-snake family, checked for sizes up to
     n = 1200, not a statement about arbitrary sequences.
   - #45: "50–1000-residue … certificates". These are conditional,
     model-relative, strong-coupling bounds.
7. **The strongest formal asset is not on `main`.** Aristotle v2's
   `FrontierCheck.check_sound` (#26/#27/#32) is a kernel-proved soundness
   theorem for a non-enumerative HP frontier-certificate checker. It closes
   part of gap #4 in `FORMAL_ASSURANCE_BOUNDARY.md`. It exists only on draft
   branches.
8. **The prior-art gate omits the closest prior art.** It omits provable
   partition-function and binding-constant bounds in protein design: OSPREY
   K*/BBK*/MARK* (Donald lab) and Toulbar2 guaranteed weighted counting
   (Viricel et al.). Any FoldBounds claim about "rigorous partition-function
   bounds" or "certified binding" must be positioned against them. See
   `docs/coupled/PRIOR_ART.md` (in progress).
9. **Negative results are preserved and not reinterpreted, in the cases
   checked.**
   - #91's own status file states that the registered weighted-ANM RNase A vs
     2AAS evaluation **FAILED**, and that later GNM work does not reverse it.
   - JEPA #95 contains an explicit "kinetic mixing falsification: stationary
     lag-one is not enough".
   - The sparse-distance work (#77) records a failed NMR-transfer control.
   - The witness experiment logs on `main` record failures.
10. **The JEPA stack (#92–#95) is green and appropriately scoped.**
    - Its 76 tests pass locally (6 optional-dependency skips), and GitHub CI on
      #95 is green, including `lean-proof` and `cpu-jepa`.
    - Its Lean file proves only the elementary facts that detailed balance
      implies stationarity and that mixtures of reversible kernels are closed.
    - The documents state correctly that this is not a refinement of the
      Python decoder, not a mixing theorem, and not kinetics.

## Reproduction log (coordinator)

| Check | Command | Result |
|---|---|---|
| main pytest | `python -m pytest` | 1260 passed (28.8 s) |
| main Lean | `lake build` | success (8041 jobs, 1.5 min with the cache) |
| main checks | `lake env lean FoldBounds/{ExecutableTests,LatticeModelTests}.lean`; `PYTHONPATH=src python3 scripts/cross_language_diff.py` | pass |
| main axioms | `lake env lean FoldBounds/AxiomAudit.lean > ax.log; python3 scripts/check_lean_sources.py ax.log` | 50 declarations, standard axioms only |
| #91 pytest (declared dependencies) | fresh venv, `pip install -e '.[test]'`, `pytest -q` | 2355 passed, 9 skipped (342 s) |
| #91 pytest (+ sympy) | the 6 GNM files | 32 passed |
| #91 Lean | `timeout 5400 lake build` | **not finished in 90 min (exit 124)**; last module MotionRealCrambinGraph |
| #91 GitHub CI | check runs of PR #91 | all Lean jobs **cancelled** at ~45 min; Python jobs green |
| #32 GitHub CI | check runs of PR #32 | `lean-proof` success (33 min) |
| JEPA tip pytest | `PYTHONPATH=src pytest -k jepa` | 76 passed, 6 skipped (torch/mdtraj) |
| JEPA Lean | `lake env lean FoldBounds/JepaKineticStationarity.lean` | compiles; 4 theorems, no sorry or native_decide |
| #95 GitHub CI | check runs of PR #95 | all green |

The integration axiom audit (296 `#print axioms` lines) is being re-run
locally. AxiomAudit does not import the Crambin module, so it can complete.
The result is recorded at the end of this file when available.

## Assurance inventory (six categories)

The categories are: **(1)** kernel-checked theorem; **(2)** independently
verified computational certificate; **(3)** software-tested algorithm with a
mathematical justification; **(4)** empirically supported scientific result;
**(5)** conditional or speculative hypothesis; **(6)** failed, refuted, or
incompletely reproduced.

The inventory is grouped by research line. A row is only as strong as its
weakest link: when a Python producer or checker sits between the Lean theorem
and the reported number, the number is (2) or (3), not (1).

| Research line (representative PRs) | Category | Scope and assumptions | Notes |
|---|---|---|---|
| Abstract finite thermodynamics: Z > 0, probability bounds, block bounds, coupled ratio, Gibbs–Bogoliubov, symmetry (main) | **1** | Finite C, real energies, model-relative | Elementary but correct; axioms audited |
| Executable rational witness/cap checker `Executable.check_sound` (main) | **1** | Enumerates the whole finite universe to check caps | Sound but not scalable |
| Lean HP lattice model: adjacency, injectivity, `checkHP_coupled_sound` (main) | **1** | Anchored 2D SAWs; rational q ≥ 1 | No bijection proof on main (it is on #91: `conformationEquivSAW`) |
| Python frontier, transfer, isoperimetric, monotone, witness, adaptive bounds (main) | **3** | 2D HP, rational q | Python checkers are not refined to Lean; exact-oracle tests up to small n |
| OEIS A001411 count agreement (main) | **2** | n ≤ 10 rooted walks | Independent published data |
| Aristotle v2 non-enumerative frontier checker `FrontierCheck.check_sound`, `checkFlat_sound`, compiled `foldcheck` (#26/#27/#32) | **1** for the checker theorem (on its own branch); **2** for accepted certificates | Trie/flat certificates; text parser unverified | Strongest formal asset; **not on main**; unverified in integrated form |
| HP exact atlases, design phase diagrams, Sturm-sequence all-q certificates, temperature transport, thermal response (#42–#57) | **2/3** | Exact rational enumeration or certified polynomial arithmetic for small n | Exact but enumerative; "all-q" claims rely on certified root isolation in Python |
| Strong-coupling witnesses, 50–1000 residues (#45) | **2** (conditional) | Needs a rigorous SAW-count bound C_n and a common reweighting interval | Conditional on stated bounds; toy model |
| 3D cubic HP: unique-native folds, cut certificates, inverse design, parity no-go, tile families (#53, #60–#66, #69, #71, #73) | **2/3**; some **1** components (CubicParityGeometry, LocalRewireConfidence) | Cubic HP toy model; PDB files are *synthetic* (marked REMARK) | The Lean components have not been verified in the integrated build |
| Infinite-size square contact density (#49) | **1/2** (SquareDensityTail.lean) | Homopolymer; abstains for some parameters | Doc honestly records the abstentions |
| Protein motion: GNM/ANM rational Hessians, inverses, rank and coupling certificates on real coordinates (#39–#43, #50, #55, #58, #62, #72, #76, #78, #79, #84, #86, #87) | **1** (checker soundness, own branches) + **2** (certificates); `MotionRealCrambinGraph` is **compiler-trusted** | Certify properties of an *elastic network model* given fixed PDB coordinates | No evidence of biological predictive validity; the registered RNase A ANM test **FAILED** (6) |
| NMR NOE restraint feasibility, minimax repair, scaffold duals (#80, #82, #88, #90) | **2/3** | Exact feasibility relative to the restraint lists and scaffold assumptions | Negative holdouts reported |
| GA/GB fold-switching fingerprints and error-correcting codes (#38, #74, #81, #85) | **3/4/5** | Small NMR ensembles (GA95/GB95, GA98/GB98) | Transfer claims are empirical and small-n; Lean only for combinatorial optimality |
| S669 assay regret, conformal and motif-blocked validation (#75, #83, #89) | **4** (empirical) + **2** (certified regret given data) | Real S669 measurements; statistical assumptions stated | About decision support, not folding physics |
| Exact kinetic barriers and pivot plateaus (#67, #70) | **2/3** | Specified Markov chains on toy models; #70 is one constructed family | The #70 title overstates generality |
| JEPA proof-gated search, flux oracle, max-flow cuts, local Hastings (#92–#95) | **1** (elementary stationarity) + **3** (Python) | Stationarity is not kinetics; flux is not mixing | Mixing falsification preserved (**6**, by design) |
| Weighted-ANM RNase A vs 2AAS registered test | **6** | Registered prediction | FAILED; not reversed by later work |

## Formal-verification findings (detail)

**F1. The integration Lean build does not complete.**
- *Evidence:* GitHub check runs on #91 (all Lean jobs `cancelled`, 45 min) and
  a local `timeout 5400 lake build` → exit 124.
- *Impact:* no theorem on #91 can be described as "verified in the integrated
  tree".
- *Recommendations:* (i) move `MotionRealCrambinGraph` out of the default
  `FoldBounds.+` glob, or replace its `native_decide` with a kernel-friendly
  certificate; (ii) gate integration CI on a target that excludes heavy
  compiler-trusted computations; (iii) re-run the per-PR Lean CI on the
  integrated tree.

**F2. `native_decide` in a library module.**
- `FoldBounds/MotionRealCrambinGraph.lean:401–402`:
  `theorem frozen_crambin_contact_network_complete : checkContactNetwork edges coords 100 = true := by native_decide`.
- `frozen_crambin_no_missing_contact` depends on it.
- Recommendation: rename the PR claim to "compiler-checked". Either list it in
  AxiomAudit, where `Lean.ofReduceBool` would then be reported and rejected by
  the checker, or re-prove it with kernel `decide` over a smaller certificate.

**F3. Gate blind spots.**
- Proposed hardening, a strengthening only:
  - (a) a recursive glob;
  - (b) reject `native_decide` outside files named `*Tests.lean`, plus an
    explicit allow-list;
  - (c) have every `theorem` in non-test library modules appear in
    AxiomAudit, generated automatically;
  - (d) include `scripts/*.lean`.
- This changes a CI gate, so it requires owner approval. It has not been
  implemented. The new CEC Lean modules ship with their own stricter checker,
  `scripts/check_coupled_lean.py`.

**F4. Python-to-Lean refinement (unchanged gap).** After Aristotle v2 the
remaining obligations are:
1. Python frontier producers and their many specialised bounds (transfer,
   isoperimetric, monotone, terminal-box) are still unverified. Only
   certificates that `foldcheck` re-checks inherit (2).
2. The text parser (`FrontierCLI`) is unverified.
3. `q = exp(βε)`: there is no real-exponential bridge.
4. No physical validity.

## Integration risks

- **File-level synthesis, not merge.** #91 combines tree snapshots of 61 PRs,
  so per-PR history and review context are not preserved in the integrated
  commits. A conflicting edit to a shared file could have been dropped
  silently. The status file says shared files were reconciled by hand.
- **28 workflows.** CI cost is high, and the Lean jobs time out (F1).
- **Undeclared optional dependencies** (sympy, numpy, torch, mdtraj), which
  cause silent test skips.
- **Two parallel Lean audit lists** (main's 50 versus #91's 296) diverge. On
  integration, the axiom gate has not run successfully in CI.

## Highest-priority opportunities

1. Land Aristotle v2's verified frontier checker on `main` on its own, with a
   CI-feasible build. It is the single largest assurance gain.
2. Close the parser and producer gaps, routing every reported HP number
   through `foldcheck`.
3. Fix F1–F3, then re-verify the integrated tree.
4. Position all "rigorous bound" claims against OSPREY/Toulbar2.
5. New research (this program): bounds on *relative* free energies from
   shared decompositions. See `docs/coupled/`.

## Integration axiom audit result (2026-11-13)

`cd /home/user/wt-integ && lake env lean FoldBounds/AxiomAudit.lean`
completed in 77 s, using the module `.olean` files from the interrupted full
build. AxiomAudit does not import `MotionRealCrambinGraph`.
`python3 scripts/check_lean_sources.py <log>` then reports:

    Source audit: 113 Lean modules, no forbidden proof bypasses
    Axiom audit: 296 declarations, standard axioms only

There are 0 occurrences of `Lean.ofReduceBool` or `sorryAx`. Interpretation:
**all 296 audited integration declarations are kernel-checked with standard
axioms.** The `native_decide` Crambin theorems are simply **absent from the
audit list**, which is finding F2. The audit's green status therefore says
nothing about them.

Caveat: this local check relied on `.olean` files from an incomplete
`lake build`. A clean CI run of the integrated tree (F1) is still required
before integration-level formal claims are made.
