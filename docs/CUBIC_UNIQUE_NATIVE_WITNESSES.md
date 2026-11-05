# Cubic HP unique native-orbit reference: three exact sequences, two reversal families

**Stacked research:** This additive branch builds on the independently audited 3D cubic contact-map and exact all-q continuum checker in #61. It adds only new files, without touching original Lean, physical protein motion, RNase A blind tests, CI, or any other research branch.

## Scientific prior-art boundary (essential)

Exact exhaustive enumeration of cubic three-dimensional HP chains **up to 19 monomers** was already reported in 2005, including designing sequences with nondegenerate ground states:

- Schiemann, Bachmann and Janke, *Exact sequence analysis for three-dimensional hydrophobic-polar lattice proteins*, *J. Chem. Phys.* (2005), DOI [10.1063/1.1814941](https://doi.org/10.1063/1.1814941); PubMed PMID 15836241.
- Schiemann, Bachmann and Janke, *Exact enumeration of three-dimensional lattice proteins*, *Computer Physics Communications* 166 (2005) 8–16, DOI [10.1016/j.cpc.2004.09.011](https://doi.org/10.1016/j.cpc.2004.09.011).
- Micheletti, Seno, Maritan and Banavar, *Protein Design in a Lattice Model of Hydrophobic and Polar Amino Acids*, *Phys. Rev. Lett.* 80 (1998) 2237, DOI [10.1103/PhysRevLett.80.2237](https://doi.org/10.1103/PhysRevLett.80.2237).

This work does **not** claim invention of the HP ground-state design problem, longer enumeration records, new protein physics, or verified experimental protein folding. Its value is a rigorously auditable small-case **FoldBounds test oracle**, independent geometry implementations, explicit structural witnesses, and exact parameter-occupancy tests.

## The geometry and native-state definition

For length n=12, enumerate all **6,989,025** 3D simple-cubic self-avoiding walks with first bond fixed +x. The independent C++17 full-walk/leaf-based contact-mask auditor from #61 exactly matches the reduced Python atlas of **8,738** contact maps. The planar square-lattice atlas is separately audited over all 30,073 walks. Both sources share the same sequence-labelled potential contact-pair bit positions, and the 2D map presence/absence check is exact.

For each of **4,096** H/P sequences, define K as the number of nonconsecutive lattice-neighbor H-H contacts, and let `M=max K` across all conformations. A *native ground conformation* is any conformation with K=M when contact attraction is favorable (`q>1`). The collection of all native conformations can contain multiple spatial symmetry orbits.

The symmetry group fixing the first oriented bond is D4 in the transverse plane. Its orbit sizes are:

- 1 for the straight walk;
- 4 for a planar nonstraight walk;
- 8 for a genuinely nonplanar walk.

A sequence is **unique up to spatial symmetry** if its entire ground-state set is exactly one such orbit. A 3D contact-map weight of 8 is not *by itself* sufficient: it could be two distinct planar weight-4 orbits with the same contact mask. The checker therefore also requires the 3D ground contact map to be **absent from the complete independently audited 2D atlas**. Under this absence and complete 3D enumeration, a ground map with multiplicity 8 necessarily represents exactly one nonplanar D4 orbit.

### Exhaustive result: exactly three nontrivial designing sequences

Among all 4,096 length-12 HP sequences, precisely three have a nonzero-contact native state in **one** cubic D4 orbit:

| Sequence | H residues | Max H-H contacts | Ground conformations | Ground contact-map bitmask |
|---|---:|---:|---:|---:|
| `HHHPHPHHPHPH` | 8 | 7 | 8 | 3,299,410 |
| `HHPHPHHPHPHH` | 8 | 7 | 8 | 26,289,413 |
| `HPHPHHPHPHHH` | 8 | 7 | 8 | 8,131,232 |

The first and third sequences are reversals of one another. The second is a palindrome. Thus there are **two** reversal-inequivalent sequence families. All three native contact maps have **no planar realizations**, and an independently repeated cubic canonical search finds **exactly one nonplanar orbit representative** for each mask.

All three representative native conformations form a **Hamiltonian path through every site of a 2×2×3 cubic lattice prism**: 12 occupied sites, no overlaps, adjacent consecutive chain beads. This gives explicit test fixtures for 3D geometry and contact scoring. The proof obligation is combinatorial and model-relative: no claim of a real aqueous protein native state.

## Exact native occupancy and why ground-state uniqueness is not foldability

Even when an HP sequence has a unique *energy-minimum symmetry orbit*, it does not follow that this fold occupies most of the thermal ensemble.

For a sequence with unique native orbit weight 8 and maximal contact score M=7, its exact probability at dimensionless contact weight q>=1 is

    P_native(q) = 8*q^7 / Z(q)
                = 8 / (8 + sum_{k=0}^6 N_k*q^(k-7)).

All terms in the denominator's sum decrease strictly with q>0, provided a lower-contact state exists. Consequently **P_native(q) increases strictly for all q>0** and has a unique crossing of any target probability strictly between P_native(1) and 1. This continuum uniqueness is an algebraic property, not sampling.

All three sequences have `P_native(1)=8/6,989,025≈1.145×10^-6`. At q=2, native occupancy is only about 0.0065–0.0073%, despite the ground state's uniqueness.

Full contact-score spectra for the two reversal-inequivalent sequence families, as exact conformation counts in order k=0 through 7:

- `HHHPHPHHPHPH` (and its reverse): **[3,956,309, 2,001,292, 721,652, 239,472, 59,228, 10,068, 996, 8]**.
- `HHPHPHHPHPHH`: **[3,262,721, 2,379,240, 971,780, 297,184, 66,920, 10,224, 948, 8]**.

The exact rational checker brackets the **unique** q at which native occupancy reaches a prescribed target:

| Family representative | P_native=50% | P_native=90% |
|---|---:|---:|
| `HHHPHPHHPHPH` | q≈134.29436 | q≈1130.57071 |
| `HHPHPHHPHPHH` | q≈128.93309 | q≈1077.24241 |

These very large dimensionless weights emphasize the toy model's limitations. They are **not temperatures in Kelvin**, realistic folding energy estimates, or experimental kinetic folding yields.

## PyMOL-compatible structural witnesses

`recover_cubic_native_witnesses` independently re-enumerates the D4-canonical cubic walks and returns an explicit *integer-coordinate path* for each of the three ground maps. It crosschecks orbit weights and H-H contact score against the classified certificates.

`lattice_witness_pdb` converts the exact lattice coordinates to a **synthetic C-alpha-only PDB display** with 3.8 Å spacing:

- the H/P chain is rendered as **placeholder ALA/GLY** solely for visualization;
- explicit `CONECT` bonds join consecutive beads;
- prominent `REMARK` lines identify the structure as **synthetic**, not an experimentally solved or predicted physical protein.

To reproduce the full experiment and write three PyMOL-viewable witness structures:

```bash
python -m pytest tests/test_cubic_native_design.py
python experiments/cubic_native_witness_demo.py --output-dir ./cubic_hp_witnesses
```

The demo independently audits the full 3D contact-map histogram with the standalone C++17 scorer when g++ is installed, or uses the slower unsymmetrized Python auditor; independently audits the planar atlas; classifies every 12-mer; independently re-enumerates coordinate witnesses; writes three synthetic PDBs plus a CSV manifest containing **exact rational occupancy brackets**.

For PyMOL, load any generated PDB and display its CA-only trace as spheres, sticks, or bonds. Do not use AlphaFold, DSSP or structural-alignment metrics on these synthetic lattice paths as if they were physical amino-acid structures.

## Trust boundary and next research hurdle

A **reference oracle** is valuable because larger proof-carrying finite-memory/lattice thermodynamic algorithms can be tested on problems with known exact answer, including native-state degeneracy and occupancy. The next consequential step would be reproducing these native-state results *without* enumerating the full cubic conformational space, then moving to richer energy models and independent held-out experimental protein structures.

Neither the Python enumeration nor the C++ audit is Lean-kernel-verified. Independent agreement is strong software validation, not a formal refinement theorem. The underlying problem and many of its conclusions have substantial published prior art, so this contribution should be assessed as infrastructure and method validation, **not** as a newly discovered law of protein folding.
