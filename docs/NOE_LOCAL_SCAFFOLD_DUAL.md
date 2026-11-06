# Exact Val69 NOE / fixed-scaffold directional dual certificate

This is a stand-alone **conditional, mathematical** research result stacked on draft PR #82. It adds new files only, changes no existing checker, Lean proof, CI, or parallel-agent branch, and is not merged. It is NOT a claim of experimental error, chemical impossibility in natural proteins, or an experimentally validated folding predictor.

## Real-data provenance

The original source files in NMRExchangeFormat/NEF are 6nbn.nef (Git blob da9769dde30d07a48ad973f1351be4e85e3931cd) and 6nbn.ambpdb.pdb (blob 5b13f49fda78a0d98e018902fab353ae73ce9bf4), under data_1_1/PDBStat_developers_it3/Case_Amber/6nbn. The matching selected restraint in base PR #80's original fixture has Git blob 96daf4e88743df4e623cdf5495f28b8e2113fc32.

Processed NOE restraint ID 470 maps uniquely to the hydrogen pair VAL69 H and VAL69 HB, with an upper limit of **2.841 Å**. The source coordinate snapshot's exact squared distance is **12,805,246 milliångström squared**, approximately 3.578442 Å separation.

We *hypothetically* hold seven nearby Val69 heavy atoms (N, CA, CB, CG1, CG2, C, O) at their source coordinates. Both hydrogens may move anywhere, provided each of their fourteen distances to these fixed heavy atoms changes by no more than a common tolerance δ. This is a model-defined neighborhood preservation constraint, including nonbonded local distances, NOT a measured chemistry uncertainty bound.

## Exact result

The source-pinned rational dual and explicit integer-coordinate primal prove:

**0.103 Å < minimum feasible δ ≤ 0.211263 Å**.

**Dual lower:** For old hydrogen x, displacement u, and fixed heavy anchor a_i, write D_i = ||x+u-a_i||² − ||x-a_i||². If the anchor distance changes by at most δ, then |D_i| ≤ E_i = 2 ceil(||x-a_i||) δ + δ², with integer-safe ceiling-square-root evaluation. Subtraction cancels the quadratic term, yielding a linear slab inequality for any two anchors and sign s:

s (a_i-a_j)·u ≤ (E_i+E_j)/2.

For each hydrogen, the untrusted producer supplies *three* signed anchor pairs and nonnegative rational dual weights. The trusted verifier recomputes their exact vector identity to the original H-to-HB contact direction v (opposite directions for the two hydrogen displacements), and their exact weighted error cap.

Closing the original H-HB separation to NOE upper limit U would require v·(u_H−u_HB) ≥ ||v|| (||v||−U) ≥ ||v||² − U ceil(||v||). The checked dual bounds it from above.

At δ=103 milliångströms:

- Minimum necessary directional closure: **2,637,307 mÅ²**.
- Maximum allowed by the exact dual:
  **459873142279194887257392 / 174477508968339989 mÅ²**, approximately 2,635,715.887 mÅ².
- The strictly positive contradiction margin is approximately **1,591.113 mÅ²**.

Every inequality uses integers and exact Fractions; no floating-point tolerance decides verification. The certificate remains valid for every real δ≤0.103 Å by monotonicity. It does NOT prove feasibility at 0.104 Å, where this particular dual is inconclusive.

**Primal upper:** An independently checked untrusted proposal moves H by (-260401,-344471,+999736) fine ticks and HB by (-341416,-167092,-21826) fine ticks, where one tick = 0.000001 Å. The verifier confirms the H-HB upper NOE limit through exact squared-distance arithmetic and proves that each of the fourteen anchor-distance changes is ≤211263 ticks by outward integer square-root intervals. This establishes feasibility for the stipulated model at δ=0.211263 Å. It does NOT establish a chemically valid repair; atom bonds, torsions, sterics, other excluded NOEs and ensemble averaging are not enforced.

## Reproduction

From the stacked repository tree run Python tests and a demo:

    python -m pytest tests/test_noe_scaffold_dual.py -q
    python experiments/val69_scaffold_demo.py

To independently replay original public source hashes and data after downloading the pinned upstream files:

    python scripts/audit_val69_raw_sources.py 6nbn.nef 6nbn.ambpdb.pdb

The source re-extractor checks Git blob hashes, nine source atom coordinates and restraint 470's precise selection and limits. The hosted integration test additionally bridges to the complete source-pinned PR #80 restraint fixture. New Python algorithms are NOT Lean-kernel-verified, regardless of green Lean infrastructure CI.

## Scientific interpretation

This is a precise conditional incompatibility between a processed NMR upper restraint and a chosen local source-coordinate preservation model. Because local geometry and fixed anchor assumptions are hypothetical, **there is no inference that the original experiment or structure is erroneous**. Different local rearrangements, additional chemistry and ensemble models can change feasibility. Farkas/LP duality, distance geometry and shortest-path restraint methods are established mathematics, not claimed mathematical invention priority. Require independent NMR/structural biology review before merge.
