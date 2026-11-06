# FoldBounds: independent real NMR distance-restraint audit with exact uncertainty certificates

## Scientific question and integrity boundary

Can the structural claims and uncertainty assumptions of a deposited molecular model be crosschecked against its **actual deposited NMR distance-restraint intervals**? The predecessor sparse-geometry pilot extracted all its distances from the same coordinate model; those distances were not independent experimental input. This work uses converted NEF restraint intervals from the NMR Exchange Format project's PDBStat test archive paired with the corresponding AMBER-prepared PDB model snapshots.

**This is a real model-versus-experimental-restraint comparison, not a blind structure-prediction result.** Deposited NMR models were themselves originally built using NMR evidence; processed NEF restraints may be derived or converted from original data. A full raw-spectrum validation, re-refinement, or unbiased biological generalization is not claimed. Stereo-specific ambiguous NOEs and effective r^-6 groups require distinct evaluation and are excluded.

Original publicly available sources (each SHA is a Git **blob** digest, *not* a SHA256 digest):

- `NMRExchangeFormat/NEF/data_1_1/PDBStat_developers_it3/Case_Amber/2png/2png.nef` — `b2b8a878bd61b278bfe76013f63dc510fd7834d9`; matching `2png.ambpdb.pdb` — `fb7e0c40aefc358cf3fb6eeb90e83526eaa5ae86`.
- `NMRExchangeFormat/NEF/data_1_1/PDBStat_developers_it3/Case_Amber/6nbn/6nbn.nef` — `da9769dde30d07a48ad973f1351be4e85e3931cd`; matching `6nbn.ambpdb.pdb` — `5b13f49fda78a0d98e018902fab353ae73ce9bf4`.

The source-pinned JSON includes only one-row `restraint_id` groups, chain A, explicit unambiguous atom names, matching residue identities and available exact coordinates in the matched PDB model. All other restraints are **excluded**, not classified as satisfied. This intentionally discards ambiguity groups, pseudoatoms, missing hydrogens, and unsupported assignment semantics. It cannot make statements about the full original restraint set.

## Mathematical certificate: exact distances and sharp perturbation classification

All PDB coordinates and restraint limits are converted to **integer milliångströms**. For two selected atoms at `p,q`, define `s=||p-q||²`, and let the experimentally processed interval be `[L,U]`. With atom-centered uncertainty balls `||δp||≤ε`, `||δq||≤ε`, the pairwise set of attainable distances is **exactly**

    [max(0, sqrt(s)-2 ε), sqrt(s)+2 ε].

Proof: the triangle inequality gives necessity; moving the endpoints along their connecting axis attains the interval extrema. Each pair can therefore be unconditionally *satisfied*, unconditionally *violated*, or *undecidable* under arbitrary allowed displacements. No floating-point approximation is needed:

- certainly satisfied iff `s≥(L+2ε)²` and `U≥2ε` and `s≤(U-2ε)²`;
- certainly violated iff `s>(U+2ε)²` or (`L>2ε` and `s<(L-2ε)²`);
- otherwise undecided.

The statement is sharp **per pair**, but classifications are not a joint feasibility solver for the protein. These are geometric perturbations of coordinate-model atoms; the ε-balls are **hypothetical**, not independently measured experimental uncertainties.

## Disjoint persistent-violation certificate

Construct a graph whose vertices are coordinate-model atom identifiers `(residue, atom_name)`, and whose edges are certainly violated single-valued restraints at ε. Every such edge must have **at least one endpoint moved by more than ε** to become satisfied. A matching of size `m` consists of disjoint edges, so every correction of all these restraints needs at least `m` independently changed atom positions, regardless of how other atoms move. The producer selects a greedy matching, while an independent verifier confirms **each edge's exact persistent violation and vertex disjointness**. `m` is a certified lower bound, not necessarily the maximum matching or a minimum repair.

This theorem is a direct matching/triangle-inequality consequence—not claimed as novel mathematics—and is strictly conditional on the converted NEF/PDB atom mapping and the specified geometric perturbation model.

## Registered observations (never tune on NMR model target)

Only **304 of 1,732** unique 2PNG restraint IDs and **553 of 1,694** unique 6NBN IDs are usable under the deliberately strict atom-specific policy. Long-range pairs (|residue separation|≥5) in this subset number only **7 and 33**, respectively, so the sample provides limited independent global-fold information.

| Data | ε (Å) | Certainly satisfied | Certainly violated | Undecided | Disjoint persistent violations |
|---|---:|---:|---:|---:|---:|
| 2PNG | 0.00 | 299 | 5 | 0 | 5 |
| 2PNG | 0.01 | 299 | 4 | 1 | 4 |
| 2PNG | 0.05 | 294 | 0 | 10 | 0 |
| 2PNG | 0.10 | 289 | 0 | 15 | 0 |
| 6NBN | 0.00 | 548 | 5 | 0 | 4 |
| 6NBN | 0.01 | 548 | 4 | 1 | 4 |
| 6NBN | 0.05 | 548 | 1 | 4 | 1 |
| 6NBN | 0.10 | 539 | 1 | 13 | 1 |

The persistent 6NBN within-residue restraint for `VAL69 H` versus `VAL69 HB` exceeds its listed upper limit by approximately **0.737 Å**. This is a *source-model-vs-converted-restraint discrepancy*. It may reflect preparation, assignment, or modelling details, not a false NMR experimental observation or an independently established biological anomaly. At ε=0.10 Å per atom, the discrepancy remains certified, but full chemical remapping and original restraint validation would be required before drawing a physical conclusion.

## Reproduction and falsification

```
python -m pytest tests/test_noe_restraint_audit.py -q
python experiments/compare_noe_restraints.py
python scripts/replay_original_noe_sources.py --fixture data/pinned_unambiguous_noe_2png_6nbn.json --entry 2png --nef /path/to/2png.nef --pdb /path/to/2png.ambpdb.pdb
python scripts/replay_original_noe_sources.py --fixture data/pinned_unambiguous_noe_2png_6nbn.json --entry 6nbn --nef /path/to/6nbn.nef --pdb /path/to/6nbn.ambpdb.pdb
```

The replay script **fails closed** if either original Git blob digest differs, then independently parses the NEF and PDB input and checks every included and excluded count against the frozen fixture. Hosted CI tests the frozen fixtures offline. It does not automatically retrieve raw NMR data or authenticate raw peak measurements.

## What this does not establish

No sequence-to-structure prediction, chemical structure refinement, unique 3D geometry from NOEs, protein dynamics, or folding thermodynamics. Absence of a flagged violation is not proof of validity: the stringent parser intentionally drops most groups and long-range information. Unlike artificial anchor distances, these processed NMR restraint intervals are external to the coordinate snapshot, but they are **not an unseen validation set** for the original NMR structure. No Lean kernel verification of this module is claimed.


## Additional negative result: selected-restraint rigidity is insufficient

Using one vertex per **mapped atom identifier** and one edge per included NOE atom pair, we verify these data do **not** uniquely determine three-dimensional atom positions on their own. This graph is **not the protein's full chemical geometry graph**. Treating every interval as an exact equality is already stronger than the original inequality constraints.

For N vertices in 3D with M distinct distance edges, the bar-and-joint rigidity Jacobian has at most M independent rows. Thus its linearized nullspace has dimension at least 3N−M. Even after removing up to six global rigid-body motions, at least max(0,3N−M−6) **internal infinitesimal zero modes** remain. This is a *first-order* result; individual nullspace directions need not integrate to finite flexes at singular configurations.

However, if the restraint graph has C connected components, the components can be independently **translated** while preserving *all* included distances exactly. Modulo one global translation, there are at least 3(C−1) continuous relative translation parameters, yielding a constructive **global non-uniqueness** witness. The separate checker constructs such an embedding and directly checks integer squared-distance invariance.

| Selected data | N mapped atom vertices | M distinct measured pairs | Components C | Independent cycles | Linearized internal zero-mode lower bound | Exact component translation parameters | Long-range (≥5 residues) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2PNG | 170 | 304 | 9 | 143 | **200** | **24** | 7 |
| 6NBN | 381 | 553 | 11 | 183 | **584** | **30** | 33 |

The observed subgraphs contain many cycles; they are **not** mostly trees. Crucially, the no-go concerns only the carefully selected **unambiguous atom-specific NOE subset**, not the entire protein structure determination. Chemical covalent bonds, angle constraints, chirality, excluded ambiguous NOEs, steric constraints and additional experiments can couple these components and remove the demonstrated translation freedom. No biological claim of intrinsic protein nonuniqueness follows.

Reproduce with `python -m pytest tests/test_noe_identifiability.py -q`.


## Positive source-only consequence: 469 new conditional atom-pair proximity bounds

A separate path-certificate module uses only converted NEF *upper limits*, never PDB
coordinates, to derive further upper bounds. For any restraint path
`x=v0, v1, ..., vk=y`, the Euclidean triangle inequality establishes

    distance(x,y) <= sum_{j=0}^{k-1} upper_limit(v_j,v_{j+1}).

The producer performs bounded Dijkstra search over the included exact-atom
restraint graph. A **separate checker** replays each proposed path by verifying
its ordered atom identities, source restraint IDs, source upper limits,
integer sum, <=8000 milliangstrom cutoff and absence of a directly measured
pair. Consequently a fault in shortest-path optimization cannot silently
produce an invalid certificate.

On the frozen source-pinned atom-specific subset:

| Source | New previously unmeasured atom pairs conditionally within 8 Å | Of those, sequence separation >= 5 |
|---|---:|---:|
| 2PNG | **110** | **1** |
| 6NBN | **359** | **39** |
| Combined | **469** | **40** |

Independent JavaScript all-pairs Floyd-Warshall crosschecked these counts;
the new Python implementation uses Dijkstra plus checked path witnesses.
No lower-bound distance or precise coordinate is inferred. In both source
examples, no **directly measured** lower limit exceeded the all-path derived
upper limit; this does **not** prove joint Euclidean feasibility of the
original intervals.

**Crucial condition:** these are mathematical implications if the selected
processed NEF upper bounds are simultaneously valid. They are *not* new
measurements, independent predictions, or proven physical contacts. Since
converted restraints and prepared coordinate models have shared provenance,
a held-out experimental test remains necessary. The additional long-range
pairs are atom-specific (not Cα contacts), and no ambiguous NOE restraint is
silently reused. The triangle inequality and shortest paths are established
mathematics; originality or discovery priority is not claimed.

Run `python -m pytest tests/test_noe_metric_closure.py -q`.


## Sharpness theorem: complete upper-only information, not just sufficient paths

The metric-closure result is actually **optimal if only source upper limits
are admitted**. For source atom `s`, let `phi(v)` be the shortest-path
distance from `s` in the restraint upper-bound graph (and set disconnected
components to an arbitrary common constant). Every observed upper-bound edge
`(u,v)` satisfies the dual Lipschitz inequality

    |phi(u) - phi(v)| <= upper(u,v).

Therefore the **explicit three-dimensional coordinate assignment**
`x_v=(phi(v),0,0)` satisfies **every selected upper-distance restraint**.
For the target atom `t`, its Euclidean distance from `s` under this
assignment is **exactly** the shortest path sum `phi(t)`.

The path witness proves **universal upper bound** `d(s,t) <= phi(t)`;
the feasible one-dimensional embedding proves **attainment**. Both sides
are independently checked using integers: path identities/weights (primal)
and every edge's Lipschitz inequality plus terminal potential separation
(dual). This is equivalent to elementary shortest-path LP duality, not
claimed as an original mathematical theorem. It rigorously establishes
that none of the 469 inferred upper bounds can be made tighter using
**only the included upper-limit graph**.

Critical qualification: the attaining one-dimensional embedding can
violate lower-distance restraints, steric exclusions, chemical bonds and
chirality. Hence the result does **not** imply optimality of the bounds
under the full converted NEF intervals or physical molecular geometry.
To improve them, FoldBounds must introduce genuinely additional data or
independently verified constraints, rather than retuning the same
upper-only graph.

Reproduce using the full `tests/test_noe_metric_closure.py` suite; the
test independently checks a primal/dual optimality pair for each of
the **469** data-derived previously unmeasured proximity bounds.
