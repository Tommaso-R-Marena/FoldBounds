# Coupled-ratio (T1) certificates for HP point mutations: algorithms and benchmark

Status: research code, model-relative. No novelty is claimed. The identities are classical
(Dinkelbach-style linear-fractional optimisation, LP duality, `Z_mut/Z_parent = E_parent[ratio]`;
the integration branch PR #91 already uses the latter with full enumeration). The question
tested here is only whether T1 (`docs/coupled/CEC_THEORY.md` §1, instantiation §4) gives useful
certified ΔΔF intervals for the 2D HP model without full enumeration.

**Headline (negative/modest).** The coupled T1 interval is always sound and always at least as
tight as the trivial `‖V‖∞` bound, but in practice it equals the trivial bound until a large
fraction of the tree is explored. At n ≥ 12 with up to 4000 regions it is within 0.1% of trivial
(median); at n = 30, 50 it is exactly trivial in every one of 2160 rows. Only at n = 10, where the
frontier reaches 74% of the leaves, does it shrink substantially. It beats "separate bounds" by a
large factor, but that is mostly because separate bounds inherit loose subtree mass caps (see §6).

## 1. Algorithm

Mutation at residue `i` (H↔P). `A` is the reference sequence, `B` the mutant, `Δ(x) = k_B(x) − k_A(x) = s·c_i(x)`,
where `c_i(x)` is the number of H sites non-consecutively adjacent to `i` and `s = −1` for H→P, `+1` for P→H.
`Z_B/Z_A = Σβ/Σα` with `α_T = Σ_{x∈T} q^{k_A(x)}`, `β_T = Σ_{x∈T} q^{k_B(x)}`.
ΔΔF = −ln(Z_B/Z_A) in kT (ε = 1, q = e^{βε}); interval width = `ln(hi/lo)`.

A *frontier* is a set of prefixes partitioning the conformation tree. For each region `T`:

* complete walk: exact `α, β`; dead end (no self-avoiding extension): `α = 0`;
* otherwise `U_T = 3^r q^{known + Ucap}` (r = remaining sites). `Ucap` charges each future contact to
  its later-placed endpoint `j`, which has at most 2 non-bonded lattice neighbours (3 if `j` is a chain end),
  and only opposite-parity H partners can touch (checkerboard). `L_T = q^{k_A(w)}` for a *witness*
  completion `w` supplied by the producer (checked), else 0. (The existing frontier code has `L = 0` for
  subtrees; witnesses are an addition that makes the separate baseline finite and T1 less degenerate.)
* ratio bounds `[q^{Δmin}, q^{Δmax}]` from `c_i ∈ [cmin, cmax]`: `cmin` = placed H contacts of `i`;
  `cmax = cmin + min(#future H partners of right parity, free lattice cells around i not needed by its unplaced chain neighbours)`
  (for unplaced `i`: `min(#parity partners, 2 or 3)`).

Then T1 (exact, `coupled/ratio.py`): `λ₋ = min Σ rα/Σα` over the box by a threshold sort, certified by
`g₋(λ₋) = Σ (r_k − λ)(L_k if r_k ≥ λ else U_k) ≥ 0` (= 0 at the optimum); `λ₊` symmetric. Baselines:

* **trivial**: one region (the root) with its ratio bounds, i.e. `|ΔΔF| ≤ cmax·ln q` (the `‖V‖∞` bound; at most `2 ln q` for interior sites).
* **separate**: independent frontier bounds on `Z_A` and `Z_B` (each with its own caps and witnesses),
  `[ΣL^B/ΣU^A, ΣU^B/ΣL^A]`. (The spec's `S₋ = ΣrL/ΣU`, `S₊ = ΣRU/ΣL` are implemented in `ratio.py` and tested to be dominated by T1, Thm 1.4, but are weaker than this baseline and are not used in tables.)
* **coupled**: T1 on the same frontier.

**Refinement.** BFS-depth frontiers (all prefixes with m sites) and adaptive splitting (priority queue):
`span` score `U·R − L·r`, or `ratio` score `U·(R − r)`. Scores are untrusted heuristics; they were fixed before
looking at results and all three configurations are reported (no selection).

**Centred growth order** (optional, `order="centred"`): the walk grows from the mutated residue as `i, i+1, i−1, i+2, …`.
This is a different encoding but gives the same ratio, because `Z` is a sum over conformations modulo
translation and rotation, the weights `q^{k}` are invariant under both, and fixing the first bond to `R`
picks exactly one representative per rotation orbit for any growth order. The checker verifies centred
certificates (the order is part of the certificate, validated as a contiguous-segment growth order), and
tests confirm exact agreement with the linear-order oracle.

## 2. Certificate format and trust boundary

```
{"seqA","seqB","q":"p/q","i","order":[...](optional),
 "frontier":[{"w":"RUL...","wit":"optional full completion"}...],
 "claim":{"lo":"p/q","hi":"p/q"}}   # bounds on Z_B/Z_A
```

`coupled/checker.py` imports only `fractions` and `json`; it shares no code with `ratio.py`, the producer or
`foldbounds.lattice`. It re-validates every word (self-avoidance, first step R), proves the frontier is a
partition (no duplicate, no ancestor/descendant pair, every valid child of every proper prefix present;
dead ends must be listed), re-derives all region data, and accepts iff `g₋(lo) ≥ 0` and `g₊(hi) ≤ 0` and
the admissible set is nonempty. It fails closed (any exception ⇒ reject). Trusted: the checker's ~150 lines,
Python `Fraction` arithmetic, and the mathematical statements of T1 and the cap argument above. Untrusted:
`ratio.py`, `hp_mutation.py` (all search, scoring, witness finding, orders). The checker is *not* formally
verified (no Lean/Coq proof); soundness rests on the cap argument and on 0 failures in the tests below.
Claims are strings/ints only, floats are rejected.

## 3. Files

`src/foldbounds/coupled/{ratio,checker,hp_mutation,hp_oracle}.py`; tests `tests/test_coupled_{ratio,hp}.py`
(14 tests: 2500 random exact instances vs vertex enumeration, edge cases, T1′ dual certificates, oracle counts
vs A001411(n−1)/4 and vs `foldbounds.lattice.exact_spectrum`, soundness along refinement paths for both orders,
checker/producer tightness agreement, 25 corrupted-certificate checks, an n=30 check);
`experiments/coupled/{run_oracle,run_large,summarize}.py`, `experiments/coupled/results/*.csv`, `SUMMARY.md`.

`ratio.py` review: the previously uncommitted file was reviewed and found correct (tie handling at group level,
`g±`, dual certificate identities); its tests passed unchanged apart from two over-strict assertions in my own
tests (`g₊(hi+ε) < 0` is only `≤ 0`).

## 4. Setup of the experiments

Oracle: n ∈ {10,12,14}, 8 random sequences per n (seed `20261113+n`, p(H)=0.5), **all** n single-site mutations,
q ∈ {3/2, 2, 3}: 71,859 rows (row = one frontier). Configurations: `linear-span`, `linear-ratio`, `centred-ratio`
(adaptive, geometric checkpoints up to 3000/4000/4000 regions) and BFS depths. Exact value from the independent
enumerator. Every row checks exact ∈ [trivial], [separate], [coupled]. The independent checker accepted
4,785 certificates (max 162 KB). Deep run: n=12, two sequences, up to 20000 regions (67% of leaves).
Large n (30, 50): 8 sequences × 5 random sites × 3 q × 3 configs, checkpoints 100/500/1500 regions, 2160 certificates all checked.

## 5. Results (see `experiments/coupled/results/SUMMARY.md` for all tables)

Soundness failures: **0** of 71,859 oracle rows (+ 0 of the deep run); checker rejections of producer certificates: **0**.

Median widths (kT), matched frontier (final frontier within a budget), rows with nonzero trivial width, coupled / separate / trivial:

| n | config | ≤100 regions | ≤1000 | ≤3000 |
|---|---|---|---|---|
| 10 | linear-span | 1.10 / 11.3 / 1.10 | 0.67 / 4.75 / 1.10 | 0.14 / 1.45 / 1.10 |
| 10 | centred-ratio | 1.08 / 11.1 / 1.10 | 0.69 / 6.18 / 1.10 | 0.53 / 4.39 / 1.10 |
| 12 | linear-span | 1.29 / 17.9 / 1.30 | 1.21 / 11.5 / 1.30 | 1.09 / 8.17 / 1.30 |
| 12 | centred-ratio | 1.29 / 18.1 / 1.30 | 1.15 / 12.3 / 1.30 | 1.09 / 9.56 / 1.30 |
| 14 | any | 1.386 / 24 / 1.386 | 1.386 / 18 / 1.386 | 1.386 / 15 / 1.386 |
| 30, 50 | any (≤1500) | 1.386 / 70–150 / 1.386 | same | same |

Pooled median ratios: **coupled/separate = 0.063 (≤100), 0.086 (≤1000), 0.097 (≤3000)**;
**coupled/trivial = 1.000, 0.999, 0.989** (0.999 for n ≥ 12 at ≤3000 regions).

Exploration fraction (trie nodes / all prefix nodes; target reached within the path, otherwise censored):

* n = 10: coupled reaches 0.5 kT in 142/180 `linear-span` rows at median 35% of the tree (1418 regions); separate in 29/180 at 74%.
  0.1 kT: coupled 111/216 at 59%; separate 0/216.
* n = 12 (4000-region cap, 13% of leaves): coupled reaches 0.5 kT in 32/253 rows; deep run (67% of leaves): coupled 48/65 at 34% (centred-ratio: 21/65 at 23%),
  0.1 kT in 38/72 at 51%; separate 5/65 (0.5) and 0/72 (0.1).
* n = 14: neither reaches 0.5 kT within 4000 regions (1.8% of leaves).

## 6. Negative findings and interpretation

1. **Coupling does not remove the ratio floor.** Until the prefix covers the mutated residue's contacts,
   regions have `r = q^{−cmax}`, `R = q^{−cmin}` with `cmin` small, so `λ₊ ≈ max R` and `λ₋ ≈ min r`: the
   coupled interval equals the trivial one. This reproduces the coordinator's prior (`hp_ratio_floor.py`, even with exact masses). Here masses are also only capped, so it is worse.
2. **Gain appears only when much of the tree is explored.** Growth with depth is exponential (≈2.6× nodes per site), so the explored-fraction requirement is a constant fraction of the tree; absolute cost grows like the tree.
3. **The large coupled/separate advantage is mostly a statement about separate bounds**, which pay for `3^r q^{Ucap}` caps (and an own cap for `B`) twice; widths of 15–150 kT are not competitive with the trivial 1.4 kT. The meaningful comparison is coupled vs trivial, which is ≈1.
4. **Growth order**: the centred order helps slightly early at n = 10/12 (coupled width at ≤1000 regions: 0.69 vs 0.71 for `linear-ratio`, 1.15 vs 1.22 at n = 12) and in the deep run reaches 0.5 kT with less of the tree (23% vs 34%), but is worse at the largest budgets for n = 10 and no better at n ≥ 14. Differences are small and not consistent; no claim.
5. The `ratio` score is not better than the neutral `span` score for T1 in these runs (`linear-span` reaches 0.5 kT more often).
6. Larger-n behaviour is unobserved: widths are exactly trivial at n = 30, 50 for ≤ 1500 regions. These runs only show that the certificates are verifiable at that size (≤ 113 KB each).

## 7. Limitations

Model-relative 2D HP only; one point mutation per certificate; sequences are random (not designed or native-like), p(H)=0.5; 8 sequences per n;
checker not formally verified; `L` bounds use a heuristic greedy witness (valid but arbitrary); `3^r` count cap and the parity/2-neighbour contact cap are loose and
are the dominant cost for both separate and coupled mass bounds; a tighter count cap (SAW growth constants) or exact subtree transfer-matrix style bounds were not tried; exploration fractions are measured on checkpoints (×1.3 granularity), so they are upper estimates of the minimum.

## 8. Reproduction

```
PYTHONPATH=src python3 -m pytest tests/test_coupled_*.py
PYTHONPATH=src python3 experiments/coupled/run_oracle.py 10 experiments/coupled/results/oracle_n10.csv 3000 8
PYTHONPATH=src python3 experiments/coupled/run_oracle.py 12 experiments/coupled/results/oracle_n12.csv 4000 8
PYTHONPATH=src python3 experiments/coupled/run_oracle.py 14 experiments/coupled/results/oracle_n14.csv 4000 8
PYTHONPATH=src python3 experiments/coupled/run_oracle.py 12 experiments/coupled/results/deep_n12.csv 20000 2 linear-span,centred-ratio
PYTHONPATH=src python3 experiments/coupled/run_large.py 30 experiments/coupled/results/large_n30.csv
PYTHONPATH=src python3 experiments/coupled/run_large.py 50 experiments/coupled/results/large_n50.csv
PYTHONPATH=src python3 experiments/coupled/summarize.py
```
Runtimes: n=10 ≈ 5 min, n=12 ≈ 14 min, n=14 ≈ 17 min, deep ≈ 12 min, large ≈ 4 and 9 min (single process each).
