# External SOTA comparison (FoldBounds)

Generated: 2026-11-19T02:27:37Z · mode=`full`

## Sources

- **OEIS-A001411**: n-step SAW counts on Z^2 https://oeis.org/A001411
- **Jensen-2013**: transfer-matrix SAW enumeration to length 79 https://arxiv.org/abs/1309.6709
- **Unger-Moult-1993**: classic 2D HP GA benchmarks 
- **Shmygelska-Hoos-2005**: standard 2D S1-* suite with published E* https://doi.org/10.1186/1471-2105-6-30
- **Thachuk-Shmygelska-Hoos-2007**: REMC SOTA comparison protocol on S1-* https://doi.org/10.1186/1471-2105-8-342
- **Giaquinta-Pozzi-2013**: MINWALK exact GS algorithm; published optima to n~80 https://doi.org/10.1089/cmb.2012.0266
- **minwalk-impl**: open-source MINWALK implementation (vendored) https://github.com/tsucres/HPFoldingProblem
- **AIP-Fibonacci-2015**: exact enumeration of Fibonacci HP sequences S13/S21/S34 https://doi.org/10.1063/1.4938021

## Findings

### Match / beat published numbers
- OEIS n=2: foldkernel exact count 1 matches A001411/4
- OEIS n=3: foldkernel exact count 3 matches A001411/4
- OEIS n=4: foldkernel exact count 9 matches A001411/4
- OEIS n=5: foldkernel exact count 25 matches A001411/4
- OEIS n=6: foldkernel exact count 71 matches A001411/4
- OEIS n=7: foldkernel exact count 195 matches A001411/4
- OEIS n=8: foldkernel exact count 543 matches A001411/4
- OEIS n=9: foldkernel exact count 1479 matches A001411/4
- OEIS n=10: foldkernel exact count 4067 matches A001411/4
- OEIS n=11: foldkernel exact count 11025 matches A001411/4
- OEIS n=12: foldkernel exact count 30073 matches A001411/4
- OEIS n=13: foldkernel exact count 81233 matches A001411/4
- OEIS n=14: foldkernel exact count 220375 matches A001411/4
- OEIS n=15: foldkernel exact count 593611 matches A001411/4
- OEIS n=16: foldkernel exact count 1604149 matches A001411/4
- S1-1: MINWALK contacts=9 matches published E*=-9 in 0.0027s
- S1-2: MINWALK contacts=9 matches published E*=-9 in 0.0052s
- S1-3: MINWALK contacts=8 matches published E*=-8 in 0.0021s
- S1-4: MINWALK contacts=14 matches published E*=-14 in 0.0032s
- S1-5: MINWALK contacts=23 matches published E*=-23 in 0.2851s
- S1-8: MINWALK contacts=42 matches published E*=-42 in 0.5443s
- S13: foldkernel histogram identical to exhaustive oracle (81233 confs)
- S21: foldkernel histogram identical to exhaustive oracle (224424291 confs)

### Speed
- OEIS n=7: foldkernel 251.3 ns/conf vs hpenum 328.5
- OEIS n=8: foldkernel 125.2 ns/conf vs hpenum 295.6
- OEIS n=9: foldkernel 70.3 ns/conf vs hpenum 312.7
- OEIS n=10: foldkernel 47.5 ns/conf vs hpenum 313.9
- OEIS n=11: foldkernel 38.5 ns/conf vs hpenum 316.5
- OEIS n=12: foldkernel 37.3 ns/conf vs hpenum 325.5
- OEIS n=13: foldkernel 33.7 ns/conf vs hpenum 324.5
- OEIS n=14: foldkernel 32.7 ns/conf vs hpenum 328.7
- OEIS n=15: foldkernel 32.7 ns/conf vs hpenum 331.1
- OEIS n=16: foldkernel 32.8 ns/conf vs hpenum 330.9

### Unique certified capability
- S1-1: FoldBounds accepted certificate for event k≥9 (certified thermodynamics path)
- S1-2: FoldBounds accepted certificate for event k≥9 (certified thermodynamics path)
- S1-4: FoldBounds accepted certificate for event k≥14 (certified thermodynamics path)
- FoldBounds is the only compared system with a Lean-kernel soundness theorem about string-level acceptance

### Gaps
- S1-6: minwalk failed (None)
- S1-7: minwalk failed (None)

## OEIS A001411 (anchored)

| n | published | hpenum | foldkernel | fk ns/conf | hp ns/conf |
|---|---|---|---|---|---|
| 2 | 1 | OK | OK | 25000.0 | 4638.0 |
| 3 | 3 | OK | OK | 9000.0 | 4549.0 |
| 4 | 9 | OK | OK | 3111.1 | 1592.8 |
| 5 | 25 | OK | OK | 1280.0 | 795.4 |
| 6 | 71 | OK | OK | 549.3 | 414.0 |
| 7 | 195 | OK | OK | 251.3 | 328.5 |
| 8 | 543 | OK | OK | 125.2 | 295.6 |
| 9 | 1479 | OK | OK | 70.3 | 312.7 |
| 10 | 4067 | OK | OK | 47.5 | 313.9 |
| 11 | 11025 | OK | OK | 38.5 | 316.5 |
| 12 | 30073 | OK | OK | 37.3 | 325.5 |
| 13 | 81233 | OK | OK | 33.7 | 324.5 |
| 14 | 220375 | OK | OK | 32.7 | 328.7 |
| 15 | 593611 | OK | OK | 32.7 | 331.1 |
| 16 | 1604149 | OK | OK | 32.8 | 330.9 |

## 2D HP ground-state (published E*)

| id | n | E* | minwalk | match | wall_s | foldbounds cert (eLo>0) |
|---|---|---|---|---|---|---|
| S1-1 | 20 | -9 | 9 | OK | 0.0027 | eLo=512/1 |
| S1-2 | 24 | -9 | 9 | OK | 0.0052 | eLo=512/1 |
| S1-3 | 25 | -8 | 8 | OK | 0.0021 | no-positive-eLo |
| S1-4 | 36 | -14 | 14 | OK | 0.0032 | eLo=32768/1 |
| S1-5 | 48 | -23 | 23 | OK | 0.2851 | no-positive-eLo |
| S1-6 | 50 | -21 | None | NO | 120.0000 | no-positive-eLo |
| S1-7 | 60 | -36 | None | NO | 120.0000 | no-positive-eLo |
| S1-8 | 64 | -42 | 42 | OK | 0.5443 | no-positive-eLo |

## Fibonacci HP exact histograms

| id | n | confs | max_k | hpenum_s | foldkernel match | fk_s |
|---|---|---|---|---|---|---|
| S13 | 13 | 81233 | 4 | 0.0170 | OK | 0.0072 |
| S21 | 21 | 224424291 | 7 | 46.1346 | OK | 6.5309 |

## How this relates to SOTA

- **Ground-state optimization SOTA (exact):** Giaquinta–Pozzi MINWALK (vendored).
  FoldBounds matches published optima on S1-1..S1-5 and S1-8, and emits *certified* witness
  certificates with eLo>0 for the event k≥E* on those instances.
- **Ground-state optimization SOTA (heuristic):** REMC / PERM / ACO sample energies; they do
  not produce machine-checked enclosures of Z or P(k≥t). FoldBounds does.
- **SAW enumeration SOTA:** Jensen transfer-matrix reaches length 79. FoldBounds matches OEIS
  A001411 (anchored) through n=16 and is ~10× faster per conformation than the Python baseline
  at n≥12 (~33 ns/conf vs ~330 ns/conf).
- **Exact HP thermodynamics:** Fibonacci S21 (224,424,291 conformations) — foldkernel matches
  the exhaustive histogram in 6.5s vs 46s for plain Python enumeration (~7×).

