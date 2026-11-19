# Certified frontier kernel v1 — report

## Claim

A Mathlib-free Lean 4 checker for `FOLDBOUNDS-KERNEL 1` certificates
(`FoldBounds.Kernel.verify` / `foldkernel`) that:

1. accepts a byte-level certificate tree and one or more `(q, t)` queries;
2. returns exact-rational enclosures of `Z`, `Z_A`, `Z_notA`, `P(A_t)` for the
   Lean HP model whenever it accepts;
3. is proved sound about its **actual** acceptance path (including the grid
   enumerator `enumFast`);
4. matches an independent Python reference/oracle on a large corpus, and
5. is substantially faster and smaller than the Mathlib-linked `foldcheck` /
   `foldcheck2` checkers without losing accuracy on matched-budget workloads.

## What is proved

* `verify_sound` / `verifyBytes_sound` / `certifiesExact_sound` /
  `certifiesInterval_sound` (module `FrontierKernelSound.lean`).
* `enumFast_eq`: on every self-avoiding represented prefix,
  `enumFast = enumGo` (grid packing is not a trusted optimisation).
* Strict Lean-side axiom audit (`FrontierKernelAudit.lean`): every declaration
  of the kernel modules depends only on `propext`, `Classical.choice`,
  `Quot.sound`. Gate-attack script rejects planted `sorry` / `native_decide`.

## Empirical accuracy

* Differential vs Python reference: **0** disagreements on 757 certificates /
  1514 queries / 12112 endpoints (`differential_summary.json`).
* Adversarial mutants: **0** disagreements on 3069 mutants (binary vs
  reference).
* Fixture corpus: 30 accept + 182 reject files, all matching expected verdicts.
* Bench S1: kernel root-`X` exact on 25/25 frozen reference queries.

## Efficiency (no accuracy loss)

From `experiments/kernel/bench/results/SUMMARY.md` (quick/full harness):

* Binary size: `foldkernel` ~5.0 MB vs `foldcheck`/`foldcheck2` ~222 MB (~43×).
* Root-`X` wall: foldcheck2 / kernel ≈ 25–160× on all-H / mixed `n=8..12`;
  pure-Python `hpenum` / kernel ≈ 4–8×.
* Multi-query: one kernel check with 4 queries vs 4 foldcheck2 runs ≈ 2000×
  wall (kernel checks the tree once).
* Non-enumerative `B` leaves: coverage up to ~1.5×10¹² conformations per
  checker step on long sequences.

Where an alternative leaf family is occasionally tighter on `p`-width, it is
orders of magnitude slower per leaf; under matched wall budgets the kernel
remains competitive or wins (see SUMMARY §S4–S5).

## Trust boundary

Proved: string/byte acceptance ⇒ model enclosure for echoed sequence and
queries. Not proved: OS file I/O, JSON printer, Lean compiler/runtime, user
intent for echoed parameters. Producer is untrusted.
