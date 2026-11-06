# Aristotle verification run report (2026-10-16)

## Toolchain

* Lean `leanprover/lean4:v4.28.0`
* Mathlib `v4.28.0`, commit `8f9d9cff6bd728b17a24e163c9402775d9e6a365`
* Python 3.11.14. `pytest` is not installed, and the environment has no
  network access.
* Machine: 8 cores, 64 GB RAM, Linux.

This run used a wrapper Lake project whose `lean_lib FoldBounds` has
`srcDir = "FoldBounds-main"`. Inside the repository the equivalent commands
run from the repository root. The repository's own `lakefile.toml` gains a
`foldcheck` executable target.

## Commands executed and results

Every command listed here was actually run, and each result is the one
observed.

| Command | Result |
|---|---|
| `lake build` (whole library, all modules including `FrontierTests`) | **pass**: "Build completed successfully (8053 jobs)" |
| `rg -n "sorry\|admit\|axiom \|implemented_by\|unsafe" FoldBounds/ scripts/*.lean`, excluding `AxiomAudit.lean` | **no matches** |
| `lake env lean FoldBounds/AxiomAudit.lean` | **pass**: 117 `#print axioms` lines, log in `experiments/results/axiom_audit.log` |
| `python3 scripts/check_lean_sources.py experiments/results/axiom_audit.log` | **pass**: "Source audit: 26 Lean modules, no forbidden proof bypasses" and "Axiom audit: 117 declarations, standard axioms only" |
| `lake env lean FoldBounds/ExecutableTests.lean` | **pass** (exit 0) |
| `lake env lean FoldBounds/LatticeModelTests.lean` | **pass** (exit 0) |
| `lake env lean FoldBounds/FrontierTests.lean` | **pass** (exit 0); all examples are `decide +kernel` |
| `lake build foldcheck` | **pass**; the first build compiled native code for the imported Mathlib modules (about 6 minutes) |
| `python3 experiments/lean_differential.py --lake-root . --runner scripts/foldcheck.lean` | **pass**: "TOTAL PROBLEMS: 0", log in `experiments/results/lean_differential.log` |
| `python3 experiments/lean_frontier_benchmark.py --lake-root . --exe .lake/build/bin/foldcheck --exact-max 16` | **pass**: "rows: 78 problems: 0", 17 minutes; results in `experiments/results/lean_frontier_benchmark.{log,json}` |
| `python3 scripts/run_tests_offline.py` (offline stand-in for pytest) | **pass**: "passed=632 failed=0" |
| `python -m pytest` | **not run**: pytest is not installed |

The axiom output records that every audited declaration depends on a subset
of `[propext, Classical.choice, Quot.sound]`. Some depend on fewer, for
example `conformation_nodup` depends on no axioms at all. No
`Lean.ofReduceBool` appears, because no `native_decide` was used.

## Differential and adversarial testing (testing, not proof)

`experiments/lean_differential.py`:

1. **Semantics.** For 18 sequences (all-H, all-P, alternating, and mixed,
   with lengths 2–8), the Lean decoder and scorer
   (`foldcheck enumerate`) produced exactly the same set of site paths and
   contact counts as an independent brute-force Python enumerator. This
   includes reflected conformations, since nothing is quotiented.
2. **Certificates.** 2160 Lean checks: 18 sequences × `q ∈ {1, 3/2, 2}` ×
   thresholds `{0, 1, 2, 4}` × several depths × modes B, T, and TW. Every
   accepted interval contained the exactly computed `Z`, `Z_A`, `Z_¬A`, and
   `P`. In B mode the Lean output was identical to Python's
   `verify_certificate` for all eight endpoints.
3. **Adversarial.** 51 mutated certificates were all **rejected**. The
   mutations were: dropped nodes (omitted branches), duplicated nodes,
   bound→exact and bound→dead relabelling, hidden descendants,
   self-intersecting witnesses, exact→dead relabelling, reversal steps,
   invalid geometry, `q < 1`, duplicate witnesses, and short witnesses.

The kernel-checked tests in `FoldBounds/FrontierTests.lean` cover the
following accepted and rejected cases, with exact expected values:

* full trie, root basic, root transfer, transfer with thresholds 0 and 3,
  basic and transfer with witnesses, rational `q = 3/2`, a mixed flat
  frontier;
* rejections: exact leaf at a nonterminal prefix, dead leaf with valid
  extensions, omitted branch, invalid geometry, node too deep, duplicate
  witness, self-intersecting witness, incomplete witness, `q < 1`, a
  sequence that is too short, flat duplicate, flat overlap, flat missing
  branch.

`hhhh_probability_exact` derives `P = 4/11` from an accepted certificate
via `check_sound`.

## What was not done or remains unresolved

* No verified parser. `FrontierCLI` and `scripts/foldcheck.lean` are
  trusted I/O.
* No formal complexity theorem; there is an informal analysis in the
  certificate specification §7.
* No refinement proof that the Python code computes the Lean-defined
  quantities. Only differential tests support that.
* The Python matching, isoperimetric, bounding-box, and lookahead caps and
  the adaptive refinement are not formalised.
* Higher-memory transfer relaxations and partially directed or
  contact-weighted lower-bound families are unproved research directions.
* Probability intervals remain nearly trivial for the benchmarked
  thresholds at frontier depths up to 8 (see the research report §C).
* No literature comparison was possible offline, so novelty is not
  assessed against published work.
* The Python suite was run with an offline stand-in for `pytest`, not with
  pytest itself.
