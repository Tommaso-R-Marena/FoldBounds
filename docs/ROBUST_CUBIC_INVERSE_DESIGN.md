# Contact-topology no-go and a proof-gated unique-orbit 3D HP inverse design

**Scope:** additive, stacked research on #66 (3D cut proofs), #64 (3D native-orbit geometry) and #61. All results concern an **anchored 3D simple-cubic H/P lattice** with a deliberately designed synthetic energy function. They do not establish a new biological folding law, a new protein-design algorithm priority, or a Lean-refined implementation. Original lattice inverse design and native-contact potentials are substantial prior art; see [Lau & Dill (1992)](https://pubmed.ncbi.nlm.nih.gov/1570343/), [Micheletti et al. (1998)](https://doi.org/10.1103/PhysRevLett.80.2237), and [Schiemann, Bachmann & Janke (2005)](https://doi.org/10.1063/1.1814941).

## Checked starting facts: 16-mer, billions of possible 3D conformations

For the sequence `HHHPPHHHHHHHPPHH` (16 sites, first oriented bond +x), the independently validated 3D binary cut-tree certificates from #66 give the exact conformational tails:

| Event | Exact count of conformations satisfying the event |
|---|---:|
| K>=10 | 2,032 |
| K>=11 | 2,032 |
| K>=12 | 0 |

Therefore the entire K=11 sector is **2,032 conformations**, none has 10 contacts, and all remaining conformations satisfy K<=9. An explicit K=9 witness exists. Importantly this is achieved by independently replaying compact **cut proofs**, not enumerating every one of the billions of complete 16-step walks.

A *new*, separately compiled C++ proof checker `cubic_ground_map_checker.cpp` extends the #66 checker with an exact full **contact-map and coordinate census** of the surviving K>=11 leaves. It reads every cut bit, independently recomputes the legal-branch coverage and geometric score upper bound before accepting every cut, and recomputes contact masks from each complete leaf. A Python consumer independently rescans **every reported coordinate path** and checks mask identities, legal self-avoidance, all C++ summaries, and all counts. Any malformed proof, invalid cut, duplicate/missing path, or inconsistent map/score causes abstention.

The audited 2,032 K=11 conformations occupy precisely **73 contact-map classes**. There are **19** map classes of size 16 and **54** classes of size 32. Each class consists entirely of 8-member orbits under the transverse D4 lattice symmetries that preserve the first oriented bond. Consequently, **every** realized K=11 contact map contains at least **two distinct spatial-symmetry orbits**.

### Exact contact-only no-go (limited to the K=11 sector)

Any energy function `F(contact_map)`—even with arbitrary heterogeneous, nonuniform, real-valued energy per contact—assigns **identical energy** to conformations with the same complete H-H contact map. Thus **any K=11 global maximum** of a contact-map-only Hamiltonian is necessarily degenerate over **at least 16 anchored conformations, i.e. at least two distinct D4 orbits**. This claim is specifically about global maxima **inside the K=11 sector**; it does *not* claim that contact-only energy cannot select a unique lower-contact fold under some other Hamiltonian.

This is an exact model-specific counterexample to the sufficiency of *contact topology alone* to uniquely specify a folded geometry.

## Sparse inverse design: one extra geometric observable is sufficient

Among 73 maps, select the size-16 contact mask (decimal) **15668042105480**. It has exactly **two D4 orbits of eight conformations each**. Their pairwise C-alpha-like **squared lattice distance** between **residue indices 0 and 12** (1-based residues 1 and 13) is, respectively, 2 and 4. Both values are constant on the relevant symmetry orbit.

Construct a new dimensionless *score* in **integer quarter-energy units** (higher means more stabilizing):

```
S(C) = 8 * #(H-H contacts in favored contact mask)
     + 4 * #(other H-H contacts)
     + 2 * 1[ squared_distance(residues 0,12) == 2 ].
```

Only the 11 selected contact-pair labels and **one two-residue distance restraint** are specified. The Boltzmann weight of an entire conformation is `q**S(C)` for a dimensionless `q>=1` defined per quarter-energy unit. This is not a physical Kelvin temperature or an experimental energy calibration.

The selected eightfold D4 orbit has **S=90**, since it realizes all 11 favored H-H contacts plus the bonus. The other eightfold orbit with the **identical** contact mask has **S=88** (it fails the distance restraint). For every other K=11 map, at least one favored H-H contact is missing and, because both maps have eleven contacts, at least one unfavored contact replaces it. Such conformations have score <=86 even if they receive the distance bonus (indeed the observed chosen sector has at least two replacements, making them <=82). All K<=9 conformations score <=74. Hence **S=90 is achieved by exactly one eightfold D4 orbit** among *every* possible anchored 3D conformation, with an exact energetic gap of at least 2 quarter-energy units.

No complete 16-mer geometric enumeration was used to establish that last statement: the K=11 sector is completely checked from one cut-tree, while the other two independently checked certificates exclude K=10 or K>=12.

## Adversarial robustness theorem: all correlated contact uncertainties

Let every pair-labelled H-H contact energy have an **arbitrary globally shared perturbation** `eta_e` satisfying `-1 <= eta_e <= 1` in the same quarter-energy units. Perturbations may be correlated/adversarial and need not be sampled, random, independent, or identically distributed. The distance bonus is kept fixed at 2.

For any different K=11 map, let `a>=1` be the number of favored contacts it lacks. Exactly `a` unfavored contacts replace them; these contacts have distinct labels. Shared perturbations on retained contacts **cancel** in an energy difference. Therefore the selected orbit beats every such competitor by at least

```
(8-1)*a - (4+1)*a + 2*(1 - other_distance_bonus)
  = 2*a + 2*(1 - other_distance_bonus) > 0.
```

For conformations sharing the target map but failing the restraint, contact energies including all errors cancel *exactly*, leaving a fixed score gap of 2. For K<=9, with at most k contacts, the target beats them by at least `(8-1)*(11-k) >= 14` quarter-energy units: the worst case is all `k` contacts lying in the favored map, and the selected distance bonus is never weaker than any competitor's bonus. Thus **the same one D4 orbit is the unique global ground-state family for EVERY vector of perturbations in [-1,1]**. This theorem follows from finite contact-set differences and the verified K=11/9 score-gap facts, not random testing.

## Rigorous native-state occupancy for every real q>=1

The target orbit comprises exactly eight geometries, all with identical favored contacts and the same bonus, so all eight have the same energy even after arbitrary pair-specific perturbations.

For every other checked high-contact geometry, record `(a,b)` with `a=#missing favored contacts` and `b=1` iff it satisfies the distance restraint. Its relative Boltzmann weight compared to one selected geometry is **at most**

```
q**(-(2*a + 2*(1-b)))
```

for **every correlated perturbation vector** with componentwise error <=1. These terms are summed across the **2024** other high-contact conformations. The checked K>=11 cut-tree verifier also reports an exact score-resolved *upper bound* on the unvisited lower-contact Boltzmann mass, in disjoint buckets `(U,C_U)` meaning at most `C_U` completions with K<=U. Because the separate gap certificate proves K<=9 for all these states, their relative weight to the selected state is at most

```
q**(-7*(11-min(U,9))).
```

Alternatively, at most `5**14 - 2032` total other low-contact conformations exist, each with relative weight at most `q**(-14)`. Take the **minimum** of these two valid upper bounds. With `R_high` the summed exact high-sector relative majorants and `R_low` that minimum, the certified native-population bound is

```
P_selected_orbit(q) >= 1 / (1 + R_high(q) + R_low(q))
```

for **every real q>=1 and every correlated contact error vector eta_e in [-1,1]**. This ratio-bound calculation is an exact rational function at rational q. Every positive exponent in the denominator's inverse powers makes the bound monotone nondecreasing with q. Only the *checked* cut-tree output may supply upper bins; no untrusted producer estimates are used as guarantees.

The full independent three-cut replay and coordinate verifier (locally ~66.6 seconds) found:

| Desired guaranteed selected-orbit occupancy | Smallest integer q (nominal errors 0) | Smallest integer q (every |eta_e|<=1) |
|---|---:|---:|
| 50% | 4 | 5 |
| 90% | 5 | **6** |
| 99% | 10 | **11** |

At `q=6`, the worst-case guaranteed native occupancy is **>96.49%** using the checked score-resolved upper bins. This `q=6` is **not directly comparable** with the earlier uniform model's `q=2816`: the new model uses nonuniform favored-contact energies and *quarter-energy scaling*, explicitly changing the Hamiltonian, so the claim is a model design result, not a 469-fold improvement under identical physics.

## Reproducible verification and safeguards

The implementation adds these files without editing existing verifier/CI sources:

- `experiments/cubic_ground_map_checker.cpp`: separately compiled full cut-tree verifier and all high-contact coordinate emitter. Contains no producer code.
- `src/foldbounds/cubic_robust_inverse_design.py`: independent Python coordinate rescoring, D4 orbit reconstruction, fail-closed design construction, adversarial margin theorem and exact occupancy bounds.
- `experiments/cubic_robust_inverse_design_demo.py`: regenerates or reads three cut certificates, **independently rechecks all of them**, writes JSON with SHA256 and exact thresholds, and produces two synthetic PyMOL PDB witnesses showing identical H-H contacts but different geometry.
- `tests/test_cubic_robust_inverse_design.py`: always-on complete small-chain proof checking, complete **13,565-walk** independent cubic n=8 exhaustive energy/probability comparison with randomized *globally shared* contact perturbations, malicious bitstream rejection, work-budget checks, plus an opt-in full n=16 three-cut independent replay (`FOLDBOUNDS_FULL_ROBUST_16=1`).

Reproduction from the stacked branch (requires g++ and installed FoldBounds Python package):

```bash
pytest -q tests/test_cubic_robust_inverse_design.py
python experiments/cubic_robust_inverse_design_demo.py \
  --output-dir ./verified_cubic_inverse_design
```

For an already saved, untrusted set of the three `.cut` files, use the `--proofs-dir` flag. Files must be named `threshold10.cut`, `threshold11.cut`, `threshold12.cut`; **none are accepted without independent checker replay**. The PyMOL files are explicitly labeled **synthetic lattice curves**, not physical protein structures.

**Assurance boundary:** The mathematical argument is exact in the given finite cubic model. The cut-tree and new coordinate checker are independently written/compiled C++ and Python reference programs, **not Lean-kernel checked or refinement-verified**. A green Lean CI audit of unchanged modules does not formalize the new algorithm. Independent review remains necessary before any scientific novelty or generalization claim.

## Exact *sharp* robustness radius, with constructive adversarial counterexample

The strict uniqueness guarantee can be strengthened beyond the tested ±1 interval, to the **entire real uncertainty-radius range `0 <= r < 2`** quarter-energy units, independently of any statistical assumptions on the perturbations. For another K=11 map with `a` favored contact losses and `a` unfavored replacements, the selected orbit's energy advantage is at least `(4-2r)*a + 2*(1-bonus)`; it is strictly positive for `r<2`. For conformations with K<=9 the advantage is at least `(8-r)*(11-K) >= 6*2 = 12` whenever `r<=2`. Identical-map other-orbit conformations trail by exactly 2 under **every** contact perturbation vector, because all contact terms cancel.

**The threshold r=2 is sharp for this Hamiltonian and selected fold.** The independently verified ground-state census includes a different K=11 map, `16217814695944`, that satisfies the same squared-distance restraint and differs by precisely **two** favored contact labels and **two** other contact labels. Apply an adversarial perturbation of `−2` on the two favored-only edges (bit positions **7, 9**) and `+2` on the two competitor-only edges (bit positions **24, 39**), with zero perturbation on every other pair label. The selected and alternative conformations then have **exactly the same score of 86 quarter-energy units**. Independently checking the energy of all 2,032 high-contact geometries finds no one exceeding 86; the universal ≤9-contact gap shows every unvisited lower-contact conformation lies at least 12 below the selected one, even at error radius 2. Therefore this is a **genuine global tie** between *different* D4 orbits—not a spurious pairwise crossing among irrelevant configurations.

The added `cubic_contact_radius_sharpness.py` independently checks the full K=11 leaf census, constructs a concrete adversarial pairwise-error vector, validates both tied high-contact structures, checks **every** high-state score under that vector, and certifies the universal low-contact separation. The exhaustive small-chain test additionally finds the global tie by independently visiting all 13,565 complete eight-residue cubic walks. For larger errors `r>2`, scaling the same adversarial assignment beyond magnitude 2 makes the competing orbit strictly outrank the original selected one.

**Interpretation:** we have an **optimal, precisely quantified robustness interval** for this particular engineered cubic-lattice Hamiltonian. The lower occupancy guarantees in the table above were computed at uncertainty radius 0 and 1; they must **not** be extended numerically to r=2, where uniqueness fails.
