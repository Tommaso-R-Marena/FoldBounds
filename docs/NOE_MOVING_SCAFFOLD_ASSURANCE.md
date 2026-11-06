# Moving-anchor chemical geometry extension (conditional, exact)

Stacked additively on FoldBounds draft PR #88 (which itself is stacked on #82 and #80). No existing files changed. This extension addresses the precise limitation that PR #88's seven heavy atoms were held fixed. It is **not** an experimental protein-folding or calibrated structure-refinement discovery.

## General mathematical implication

Let H⁰ and A⁰ⱼ be frozen source model coordinates. Consider any new hydrogen and heavy anchor H', A'ⱼ. If the heavy anchor moves by at most ηⱼ and the H–A'ⱼ distance differs from the source H⁰–A⁰ⱼ distance by at most δ_Hj, the reverse triangle inequality gives

    | |H' - A⁰ⱼ| - |H⁰ - A⁰ⱼ| | <= ηⱼ + δ_Hj.

This inequality holds for arbitrary simultaneous movements and correlations. Under homogeneous budgets η,δ, **all 14** source-anchored hydrogen-to-heavy distances change by at most η+δ. Replaying PR #88's exact rational box certificate at 199 mÅ rules out satisfaction of the processed Val69 H–HB upper NOE when

    0 <= η + δ <= 199 mÅ = 0.199 Å.

Examples: an allowed 0.050 Å heavy-anchor shift still excludes restoration if every associated hydrogen-to-moving-heavy distance changes by at most 0.149 Å; a 0.100 Å anchor shift pairs with 0.099 Å. These are **stipulated scenario limits**, not measured protein variability. If the source scaffold changes more substantially, no conclusion is made.

## Heterogeneous / partially unconstrained anchor theorem

The new verify_heterogeneous_moving_scaffold uses the LP/Farkas directional multipliers and full exact rational vector identities from PR #88. Each signed slab associated with anchor j receives its **own** conservative squared-distance error cap

    E_Hj = 2 ceil(|H⁰-A⁰ⱼ|) (ηⱼ+δ_Hj) + (ηⱼ+δ_Hj)²

(all in integer milliångströms). A nonnegative exact linear combination bounds the direction needed to close the NOE. If the summed upper directional support is strictly less than

    |v|² - U ceil(|v|),    v=HB⁰-H⁰,

then the NOE is impossible within the declared budgets. A verifier recomputes the normal-vector identities and weighted exact support, accepts only valid nonnegative integer/rational weights, and rejects missing budgets for **used** anchors.

Concrete certificate: ηⱼ=40 mÅ, δ_Hj=δ_HBj=50 mÅ only for heavy anchors N, CA, CB, CG2, O. The anchors C and CG1 have **no imposed motion or local-shell constraints at all**; the dual does not use them, so their behavior is irrelevant to this particular incompatibility proof. The dual proves the NOE cannot be restored under the five modeled anchor budgets despite allowing the other two anchors arbitrary motion.

## Reproduction and limits

    PYTHONPATH=src python -m pytest tests/test_noe_scaffold_dual.py tests/test_noe_scaffold_branch_bound.py tests/test_noe_moving_scaffold.py -q

The source-pinned raw 6NBN NEF/PDB evidence, exact checked fixed-scaffold proof, and constructive feasible 0.211263 Å point are inherited from PR #88. The new theorem is conditional on scenario budgets; **no raw experimental error bars** exist for these original-structure-derived shell constraints. This does not include the other NOEs, chemical energy, steric clash constraints or NMR ensemble averaging. As with the parent, Python proof-checking is not a Lean kernel theorem. Generalized triangle inequalities and exact rational LP duality are established mathematics.
