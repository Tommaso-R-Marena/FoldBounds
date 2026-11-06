# Exact, adversarially robust mutation-driven compactness transitions in a 3D lattice

**Research status:** mathematical/computational proof-of-concept, **not** a predictor of experimentally observed proteins, not an atomistic energy model, and **not yet fully Lean-verified**.

## Falsifiable central result

A fixed, mutation-independent conformational event—**at least two nonconsecutive backbone-neighbor contacts** on the three-dimensional cubic lattice—can undergo a certified minority-to-majority switch after **one H/P substitution** over a continuous box of Hamiltonian parameters, even with arbitrary per-state energy errors.

- Parent: `PHPHPHPH` (binary mask 170, 8 residues).
- Mutant: `HHPHPHPH` (binary mask 171, P1H).
- All **13,565** anchored 3D self-avoiding conformations of length 8 are included, with one initial bond fixed along +x; no later conformation is omitted.
- Hamiltonian weights are **exactly** `q**(number of nonconsecutive H-H contacts) * b**(number of changes in bond direction)`. `q` is a dimensionless HH energy preference and `b` a dimensionless bending preference; neither is experimentally fitted.
- Uncertainty box: `q in [8,12]`, `b in [3/4,5/4]`, and for **every individual conformation independently**, true/nominal weight multiplier in `[5/6,6/5]`.

The exact rational Bernstein enclosure and sharp adversarial reweighting formula certify:

| | Model population of compact class | Corrected population under worst allowed errors |
|---|---|---|
| Parent | `652185/3730687` to `9160825/40721119` | at most `2198598/7458647` (~29.48%) |
| Mutant | `11359539/16531846` to `20250625/23950326` | at least `31554275/52243503` (~60.40%) |

Thus **every admissible parameter/error realization** has `P_parent(compact) < 0.5` and `P_mutant(compact) > 0.5`. Even the worst-case mutant compactness exceeds the worst-case original compactness by an absolute **`120489737357131/389665846920441` ≈ 30.92 percentage points**.

At fixed `q=8`, `b=1`, the exact parent compact population is `552/2713` (~20.35%) and the mutant probability is `16152/22423` (~72.03%).

These percentages are not GA95/GB95 measurements, native-fold populations, or evidence that an eight-residue peptide is stable in solution.

## The mathematical argument

For every self-avoiding conformation `c`, let `E(c)` be the *geometry-only* event indicator, `K_s(c)` the HH contact count under sequence `s`, and `R(c)` the bond-turn count. The partition function and compact-state mass for any mutation sequence `s` are

`Z_s(q,b) = sum_c q**K_s(c) * b**R(c)`;
`A_s(q,b) = sum_{c:E(c)} q**K_s(c) * b**R(c)`.

Both are exact polynomials with nonnegative integer coefficients. The set of `c` never changes with sequence; only their energy weights change. There is no requirement that population be monotone in `q` or `b`.

For a closed rational rectangle, transform the polynomial pair `(A_s, Z_s-A_s)` into tensor-product Bernstein bases. The basis functions are nonnegative and partition unity. Every denominator Bernstein coefficient is positive. Thus at each `(q,b)` the probability is a convex weighted average of the coefficient ratios `A_ij/Z_ij`. Its range is enclosed by the minimum and maximum of those exact ratios. Subdividing the box may refine the interval but is not necessary for the certificate above: the initial Bernstein convex-hull bound equals the attained corner values in this chosen example. **This is a mathematical argument, not a formally verified proof of the Python coefficient-generation implementation.**

With independent unknown conformation-specific multipliers `w_true(c)/w_nom(c)` in `[l,u]`, a nominal event population `p` has the sharp, attainable bounds

`l*p / (l*p+u*(1-p)) <= P_true(event) <= u*p / (u*p+l*(1-p))`.

The sharp extremes assign the low multiplier to all event states and the high multiplier to all nonevent states (or vice versa). No assumptions of independent errors, Gaussianity or unbiased sampling are needed; a common energy shift cancels exactly. The bounds become biologically informative **only if the physical energy-model error range is independently justified**.

These abstract weighted-ratio and energy-error inequalities already have separate Lean formalizations in the preceding [robust ensemble research](https://github.com/Tommaso-R-Marena/FoldBounds/pull/34); the new exact histogram/fingerprint/Bernstein implementation has independent Python oracles, not a newly added end-to-end Lean theorem.

## Lossless geometry compression

The full 13,565-conformation ensemble reduces to **171 distinct triples** of `(complete nonconsecutive contact-edge set, bend count, geometric compactness label)`, each with exact integer multiplicity. For every H/P sequence `s`, the contact energy `K_s` is determined entirely by its edge set; hence grouping preserves **exactly** each term of `Z_s` and `A_s` for any positive `(q,b)`. This is a sufficient statistic for this coarse-grained model, not an approximate structural fingerprint.

The 256 binary sequences have 1,024 distinct one-residue mutation edges. Checking the exact uncertainty box on every edge yields **36 robust minority-to-majority edges**, representing **18 orbits after sequence-reversal symmetry**. All robust single-site switches here are P-to-H. The counted cases include only strict certified switches; ambiguous cases are not counted as failures of the underlying true population transition.

## An emergent, contact-forbidden two-site interaction

Because the cubic lattice is bipartite, backbone residues whose indices differ by an even number **cannot be lattice neighbors** in any self-avoiding walk. A pair of such mutations therefore cannot gain direct HH contact energy between each other. Nonetheless the ensemble can induce nonadditive responses.

At exactly `q=8,b=1`, all `C(8,2)*2**6 = 1792` two-site mutational cycles were compared using **exact rational odds**. Of the **768** cycles with equal-parity mutated sites, **504** have nonzero odds epistasis despite direct lattice contact between the sites being geometrically impossible. In **64** of these cycles, the parent and both single mutants have compact-state population below 50% but the double mutant exceeds 50%. **30** such double-only crossovers persist under adversarial per-state multipliers in `[10/11,11/10]` at this one fixed `(q,b)` point; we did **not** find those double-only cases uniformly robust over the wider contact/bending parameter box.

A concrete instance is `HPPPPPPP` with P4H and P6H. At `q=8,b=1`, exact populations for parent, position-4 single, position-6 single and double mutant are respectively:

- `552/2713` (20.35%), `8640/26333` (32.81%), `7660/19081` (40.14%), and `24908/43217` (57.63%).

Neither mutation individually crosses 50%; together they do. This is **ensemble-mediated epistasis in a toy Hamiltonian**, not experimental evidence of protein allostery.

## Independent falsification tests

- An entirely separate depth-first cubic-lattice enumerator agrees on exact `(HH contacts, bends, compactness)` histograms for both spotlight sequences, including all 13,565 states.
- An independently calculated dense 9×9 rational parameter grid stays within the exact Bernstein interval for each sequence; equality at the worst-case corner values shows the main result's Bernstein interval is attained.
- Reweighting the exact event/complement masses at each corner with adversarial factors independently matches the symbolic lower/upper transforms; extrema are constructively attainable.
- All 1,024 mutation edges are visited and their exact one-step status is checked; the full 256-sequence scan is reproducible.
- Exact sequence-reversal invariance is checked for every binary 8-mer.
- Input type validation, impossible switches, 400%-range adversarial model uncertainty, malformed mutation certificates and forged landscape counts are rejected.
- Test count in the initial implementation: **25 targeted tests passed** (including parameterized cases). CI must still independently run on the final branch head before merge readiness.

## Experimental translational target: GA95/GB95, not yet analyzed atomistically

Published protein experiments on GA95/GB95 report switching between a 3-alpha and a 4-beta-plus-alpha topology across a small number of mutations; PDB 2KDL and 2KDM contain solution-NMR coordinate ensembles. This result is *not* a prediction of either structure and the current runtime could not retrieve the coordinates. The relevant 2009 and 2025 studies describe the actual mutational mechanism and aromatic packing interactions:

- https://doi.org/10.1073/pnas.0906408106
- https://doi.org/10.1039/D4SC04951A
- https://www.rcsb.org/structure/2KDL
- https://www.rcsb.org/structure/2KDM

A genuine real-protein folding result would require: coordinate/model inputs; calibrated atomistic or statistically valid ensemble energy bounds; preservation of all important conformational basins; prospective mutation predictions blinded to measured fold populations; and successful independent NMR/thermodynamic validation against strong baseline tools. None of that has been completed here. Generic uncertainty mathematics, Bernstein bounds and lattice-polymer transfer methods are not claimed to be novel.

## Reproduce

```bash
PYTHONPATH=src python -m pytest tests/test_mutation_switch.py -q
PYTHONPATH=src python experiments/certified_mutation_screen.py > /tmp/mutation_screen.json
python -m json.tool /tmp/mutation_screen.json >/dev/null
```

`experiments/certified_mutation_screen_results.json` contains the canonical exact-rational output. The certificate checker recomputes actual histograms and coefficient bounds; it does **not** trust a JSON value. Python 3.11/3.12/3.13 GitHub CI should compare recomputed parsed JSON to the checked-in reference.

## Independent longer-chain replication: all 308,981 ten-residue cubic conformations

The same mutation-driven robust compactness switch was also demonstrated on **all 308,981 anchored length-10 cubic self-avoiding walks**, with a *stricter mutation-independent event of at least three nonconsecutive geometric contacts*.

- Parent sequence: `PHPHPHPHPH` (alternating P/H, no possible HH contact because of lattice parity).
- Mutant sequence: `HHPHPHPHPH` (one P1H substitution).
- Exact fingerprint grouping of the full 3D state space: **1,852 unique contact-graph/bend/compactness fingerprints** retaining every geometric contact edge; all 308,981 original states are represented with exact multiplicity.
- Parameter-and-error rectangle is **identical** to the eight-mer: `q∈[8,12]`, `b∈[3/4,5/4]`, independent conformational weight corrections in `[5/6,6/5]`.
- Nominal *parent* compact-state probability lies in `[96869817/882992602, 3166978825/21138324794]` (~10.97%–14.98%).
- Nominal *mutant* compact-state probability lies in `[2983192704/5014904995, 295886625000/375469692019]` (~59.49%–78.80%).
- **Worst-case original upper:** `4560449508/22531795477` (~20.24%).
- **Worst-case mutant lower:** `2071661600/4103373891` (~50.49%).
- Strict separation across *every* parameter and admissible independent error correction: `27965026026403587572/92456381277673691007`, or at least **30.25 percentage points**.
- At nominal q=8,b=1, the original population is exactly `40932/308981` (~13.25%) and the mutant is `1264000/1997213` (~63.29%).

Unlike the original compressed-certificate calculation, `experiments/independent_n10_oracle.py` uses a **separately coded depth-first lattice enumerator** and directly scores all 308,981 complete states. It obtains the same 40,932 compact conformations and the same exact weighted mutant probability at q=8,b=1. This independently validates the reference-case 10-mer histogram and lossless fingerprint totals; the entire continuous parameter box is certified through the rational Bernstein bound, which is also checked against explicit corner values.

The longer-chain replication still does **not** show that the compactness class corresponds to a true molecular fold, or that the parameter/error budget matches a realistic solvent force field. The 10-mer enumeration is exponential in chain length, although its lossless contact-graph grouping makes many sequence changes cheaper to score.

Reproduce:

```bash
PYTHONPATH=src python -m pytest -q tests/test_mutation_switch_10.py
PYTHONPATH=src python experiments/certified_mutation_10_benchmark.py > /tmp/10mer.json
PYTHONPATH=src python experiments/independent_n10_oracle.py
```

## A falsifiable energy-margin diagnostic: unique critical contact weight

The preceding uniform parameter-box certificates show robust majority switching, but they do not indicate how quickly it fails when the contact preference weakens. At the **worst permitted bend weight b=3/4**, the mutant's robust lower probability exceeds 50% precisely when

`factor_low * E_mutant(q,b) - factor_high * C_mutant(q,b) > 0`,

where `factor_low=5/6`, `factor_high=6/5` and E/C are the exact compact-event/complement mass polynomials in q. Since the complete state histogram is finite, this inequality has **exact rational coefficients**. The coefficients have exactly one sign variation, so Descartes' rule proves at most one positive real root; verified opposite endpoint signs and rational bisection establish the existence and a unique positive crossing.

- Eight-residue P1H mutant `HHPHPHPH`: critical q in `[1711769203061771/281474976710656, 427942300765443/70368744177664]`, approximately **6.081425862667**.
- Ten-residue P1H mutant `HHPHPHPHPH`: critical q in `[1112769948121909/140737488355328, 2225539896243819/281474976710656]`, approximately **7.906706032102**.

Both brackets are generated and verified using exact rational coefficients and 48 successive bisections in `src/foldbounds/mutation_root.py` and `experiments/certified_mutation_roots.py`. These are **strictly model-conditional** mathematical stability thresholds, *not measured folding transition energies*. The ten-mer starts its robust box at q=8, just ~0.0933 above the crossover: this should be interpreted as a narrow energy-model margin, not strong real-world robustness. For the specified examples the independent Bernstein parameter-box calculations attain their lower endpoint at q=q_lower, b=3/4, allowing the critical fixed-b root to diagnose failure as q is reduced.

The existence, uniqueness, and bisection are proved mathematically by elementary polynomial sign arguments and independently tested in Python; the concrete root checker is **not a new Lean kernel proof**. The additional eight root-specific tests include malicious-coefficient and incorrect-bracket rejections.