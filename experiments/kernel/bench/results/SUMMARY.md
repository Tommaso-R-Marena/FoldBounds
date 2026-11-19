# FoldBounds certified frontier — matched-budget benchmark SUMMARY

Generated: 2026-11-19T00:59:06+00:00

## Machine / binaries

- Host: `cursor` · Linux `6.12.113+` · 4 CPUs
- Load average at harness start: [5.45, 5.51, 3.02]
- Mode: `quick` · rows written this run: 5
- `kernel`: 5,209,296 bytes · sha256 `b8ba2338921c795e…`
- `v1`: 222,524,568 bytes · sha256 `6d5a4c490cf06397…`
- `v2`: 222,652,512 bytes · sha256 `b2da8e25f1530d1e…`

## Commands

```bash
export PATH="$HOME/.elan/bin:$PATH"
cd /workspace && . .venv/bin/activate
lake build foldkernel   # foldcheck / foldcheck2 already built
python experiments/kernel/bench/run_bench.py --quick
python experiments/kernel/bench/run_bench.py --full
```

## Findings

### (i) Kernel better
- S2/allH/n=8: foldcheck2/kernel wall speedup = 25.2× (0.05712s / 0.002269s).
- S2/allH/n=8: hpenum/kernel wall speedup = 7.4× (0.0169s / 0.002269s).
- S2/allH/n=9: foldcheck2/kernel wall speedup = 28.9× (0.06895s / 0.002383s).
- S2/allH/n=9: hpenum/kernel wall speedup = 7.0× (0.01657s / 0.002383s).
- S2/allH/n=10: foldcheck2/kernel wall speedup = 44.1× (0.1025s / 0.002323s).
- S2/allH/n=10: hpenum/kernel wall speedup = 7.6× (0.01764s / 0.002323s).
- S2/allH/n=11: foldcheck2/kernel wall speedup = 83.8× (0.2685s / 0.003203s).
- S2/allH/n=11: hpenum/kernel wall speedup = 7.5× (0.02403s / 0.003203s).
- S2/allH/n=12: foldcheck2/kernel wall speedup = 161.5× (0.6413s / 0.003971s).
- S2/allH/n=12: hpenum/kernel wall speedup = 7.0× (0.02776s / 0.003971s).
- S2/mixed/n=8: foldcheck2/kernel wall speedup = 28.1× (0.06744s / 0.002398s).
- S2/mixed/n=8: hpenum/kernel wall speedup = 7.4× (0.01768s / 0.002398s).
- S2/mixed/n=9: foldcheck2/kernel wall speedup = 19.2× (0.07623s / 0.003968s).
- S2/mixed/n=9: hpenum/kernel wall speedup = 4.5× (0.01785s / 0.003968s).
- S2/mixed/n=10: foldcheck2/kernel wall speedup = 28.8× (0.1143s / 0.003966s).
- S2/mixed/n=10: hpenum/kernel wall speedup = 4.5× (0.01783s / 0.003966s).
- S2/mixed/n=11: foldcheck2/kernel wall speedup = 65.1× (0.2407s / 0.003695s).
- S2/mixed/n=11: hpenum/kernel wall speedup = 5.6× (0.02062s / 0.003695s).
- S2/mixed/n=12: foldcheck2/kernel wall speedup = 163.2× (0.6371s / 0.003903s).
- S2/mixed/n=12: hpenum/kernel wall speedup = 6.9× (0.02702s / 0.003903s).
- S3/dev-fixed-12/d=4: v1:B is 18.0× slower per leaf than kernel B (222 vs 12.3 µs).
- S3/dev-fixed-12/d=4: v1:T is 4459.5× slower per leaf than kernel B (5.5e+04 vs 12.3 µs).
- S3/dev-fixed-12/d=4: v2:M0 is 4936.9× slower per leaf than kernel B (6.09e+04 vs 12.3 µs).
- S3/dev-fixed-12/d=4: v2:M2 is 13270.3× slower per leaf than kernel B (1.64e+05 vs 12.3 µs).
- S3/dev-fixed-12/d=6: v1:B is 25.8× slower per leaf than kernel B (127 vs 4.92 µs).
- S3/dev-fixed-12/d=6: v1:T is 3372.5× slower per leaf than kernel B (1.66e+04 vs 4.92 µs).
- S3/dev-fixed-12/d=6: v2:M0 is 3830.9× slower per leaf than kernel B (1.88e+04 vs 4.92 µs).
- S3/dev-fixed-12/d=6: v2:M2 is 7300.9× slower per leaf than kernel B (3.59e+04 vs 4.92 µs).
- S3/dev-fixed-12/d=8: v1:B is 43.3× slower per leaf than kernel B (138 vs 3.19 µs).
- S3/dev-fixed-12/d=8: v1:T is 1430.2× slower per leaf than kernel B (4.57e+03 vs 3.19 µs).
- S3/dev-fixed-12/d=8: v2:M0 is 1634.9× slower per leaf than kernel B (5.22e+03 vs 3.19 µs).
- S3/dev-fixed-12/d=8: v2:M2 is 2013.3× slower per leaf than kernel B (6.43e+03 vs 3.19 µs).
- S3/dev-allH-10/d=4: v1:B is 11.5× slower per leaf than kernel B (111 vs 9.67 µs).
- S3/dev-allH-10/d=4: v1:T is 1758.6× slower per leaf than kernel B (1.7e+04 vs 9.67 µs).
- S3/dev-allH-10/d=4: v2:M0 is 1873.6× slower per leaf than kernel B (1.81e+04 vs 9.67 µs).
- S3/dev-allH-10/d=4: v2:M2 is 3390.8× slower per leaf than kernel B (3.28e+04 vs 9.67 µs).
- S3/dev-allH-10/d=6: v1:B is 21.4× slower per leaf than kernel B (84.5 vs 3.94 µs).
- S3/dev-allH-10/d=6: v1:T is 832.1× slower per leaf than kernel B (3.28e+03 vs 3.94 µs).
- S3/dev-allH-10/d=6: v2:M0 is 985.7× slower per leaf than kernel B (3.89e+03 vs 3.94 µs).
- S3/dev-allH-10/d=6: v2:M2 is 1267.9× slower per leaf than kernel B (5e+03 vs 3.94 µs).
- S3/dev-allH-10/d=8: v1:B is 30.8× slower per leaf than kernel B (95.8 vs 3.11 µs).
- S3/dev-allH-10/d=8: v1:T is 189.5× slower per leaf than kernel B (589 vs 3.11 µs).
- S3/dev-allH-10/d=8: v2:M0 is 254.6× slower per leaf than kernel B (792 vs 3.11 µs).
- S3/dev-allH-10/d=8: v2:M2 is 270.0× slower per leaf than kernel B (840 vs 3.11 µs).
- S4/dev-fixed-12/d=4: kernel B has the tightest (or tied) p_width = 1.
- S4/dev-allH-10/d=4: kernel B has the tightest (or tied) p_width = 1.
- S4/dev-allH-10/d=6: kernel B has the tightest (or tied) p_width = 1.
- S6: kernel non-enumerative check covers ~1.54e+12× more conformations than checker steps (ho-long-30, d=2, conf=1543066804034, steps=1).
- S7/n_queries=1: foldcheck2 sum-of-runs / kernel multi-query wall = 544.8× (1.781s / 0.00327s).
- S7/n_queries=4: foldcheck2 sum-of-runs / kernel multi-query wall = 2295.2× (7.813s / 0.003404s).
- S8: foldcheck binary is 43× larger than foldkernel; foldcheck2 is 43× larger.

### (ii) Ties / no accuracy loss
- S1: kernel root-X exact on 25/25 reference queries (disagreements with reference = 0).
- S1: foldcheck2 root-X agreed with kernel on 15 finished comparisons (disagreements = 0; timeouts = 0).

### (iii) An alternative better
- S4/dev-fixed-12/d=6: tightest p_width is v2:M2 = 0.984607 vs kernel B = 0.98931 (wall 2.60929s vs 0.00265342s).
- S4/dev-fixed-12/d=8: tightest p_width is v2:M2 = 0.874195 vs kernel B = 0.946415 (wall 3.57206s vs 0.00419907s).
- S4/dev-allH-10/d=8: tightest p_width is v2:M0 = 0.768294 vs kernel B = 0.937094 (wall 0.486735s vs 0.00405402s).

### (iv) Inconclusive / timeouts / skipped
- (none in this run)

## S1 — Exactness (just as much accuracy)

Protocol: root-`X` certificates; kernel must match the frozen exact reference on every endpoint; foldcheck2 root-`X` compared when it finishes within 300s. Repeats: median/min wall.

| metric | value |
|---|---|
| kernel rows | 25 |
| kernel exact_match | 25 |
| kernel encloses_reference | 25 |
| kernel failures/timeouts | 0 |
| foldcheck2 agrees_with_kernel | 15 |
| foldcheck2 disagrees | 0 |
| foldcheck2 timeouts | 0 |

## S2 — Exhaustive speed

Protocol: root-`X` (kernel, foldcheck2) and pure-Python `hpenum` on all-H and one mixed family; stop each alternative at first over-cap run.

| kind | n | method | wall_s_median | verify_s | µs/conf | timed_out |
|---|---|---|---|---|---|---|
| allH | 8 | hpenum | 0.0169017 | 0.000162499 | 0.299262 | False |
| allH | 8 | kernel | 0.00226883 | 6.7e-05 | 0.123389 | False |
| allH | 8 | v2 | 0.0571236 | 0.005 | 9.2081 | False |
| allH | 9 | hpenum | 0.0165679 | 0.000465792 | 0.314937 | False |
| allH | 9 | kernel | 0.00238333 | 0.000106 | 0.07167 | False |
| allH | 9 | v2 | 0.0689476 | 0.016 | 10.8181 | False |
| allH | 10 | hpenum | 0.0176368 | 0.00123311 | 0.303199 | False |
| allH | 10 | kernel | 0.00232251 | 0.000201 | 0.0494222 | False |
| allH | 10 | v2 | 0.102459 | 0.051 | 12.54 | False |
| allH | 11 | hpenum | 0.0240296 | 0.00351218 | 0.318565 | False |
| allH | 11 | kernel | 0.00320261 | 0.000432 | 0.0391837 | False |
| allH | 11 | v2 | 0.268511 | 0.181 | 16.4172 | False |
| allH | 12 | hpenum | 0.0277634 | 0.00987784 | 0.328462 | False |
| allH | 12 | kernel | 0.00397131 | 0.001065 | 0.0354138 | False |
| allH | 12 | v2 | 0.641282 | 0.53 | 17.6238 | False |
| allH | 13 | kernel | 0.0057246 | 0.002788 | 0.034321 | False |
| allH | 14 | kernel | 0.013509 | 0.007295 | 0.0331027 | False |
| mixed | 8 | hpenum | 0.0176777 | 0.000120042 | 0.221072 | False |
| mixed | 8 | kernel | 0.00239815 | 6.7e-05 | 0.123389 | False |
| mixed | 8 | v2 | 0.0674357 | 0.005 | 9.2081 | False |
| mixed | 9 | hpenum | 0.0178546 | 0.000417502 | 0.282287 | False |
| mixed | 9 | kernel | 0.00396769 | 0.000104 | 0.0703178 | False |
| mixed | 9 | v2 | 0.0762298 | 0.016 | 10.8181 | False |
| mixed | 10 | hpenum | 0.017832 | 0.000818153 | 0.201169 | False |
| mixed | 10 | kernel | 0.00396613 | 0.000185 | 0.0454881 | False |
| mixed | 10 | v2 | 0.114254 | 0.048 | 11.8023 | False |
| mixed | 11 | hpenum | 0.0206208 | 0.00202693 | 0.183849 | False |
| mixed | 11 | kernel | 0.00369504 | 0.000406 | 0.0368254 | False |
| mixed | 11 | v2 | 0.240695 | 0.166 | 15.0567 | False |
| mixed | 12 | hpenum | 0.0270246 | 0.00804626 | 0.267558 | False |
| mixed | 12 | kernel | 0.003903 | 0.001019 | 0.0338842 | False |
| mixed | 12 | v2 | 0.63713 | 0.529 | 17.5905 | False |

## S3 — Per-leaf cost at identical frontiers

Protocol: matched `uniform(d)` frontiers; µs/leaf = verify_s / n_leaves.

| instance | d | method | leaves | verify_s | µs/leaf | wall_s_median |
|---|---|---|---|---|---|---|
| dev-allH-10 | 4 | kernel:B | 9 | 8.7e-05 | 9.66667 | 0.00241098 |
| dev-allH-10 | 4 | v1:B | 9 | 0.001 | 111.1 | 0.0724612 |
| dev-allH-10 | 4 | v1:T | 9 | 0.153 | 1.7e+04 | 0.209277 |
| dev-allH-10 | 4 | v2:M0 | 9 | 0.163 | 1.811e+04 | 0.228121 |
| dev-allH-10 | 4 | v2:M2 | 9 | 0.295 | 3.278e+04 | 0.350525 |
| dev-allH-10 | 6 | kernel:B | 71 | 0.00028 | 3.94366 | 0.00236377 |
| dev-allH-10 | 6 | v1:B | 71 | 0.006 | 84.507 | 0.0597718 |
| dev-allH-10 | 6 | v1:T | 71 | 0.233 | 3282 | 0.284048 |
| dev-allH-10 | 6 | v2:M0 | 71 | 0.276 | 3887 | 0.329078 |
| dev-allH-10 | 6 | v2:M2 | 71 | 0.355 | 5000 | 0.405283 |
| dev-allH-10 | 8 | kernel:B | 543 | 0.001689 | 3.1105 | 0.00405402 |
| dev-allH-10 | 8 | v1:B | 543 | 0.052 | 95.7643 | 0.103393 |
| dev-allH-10 | 8 | v1:T | 543 | 0.32 | 589.3 | 0.373606 |
| dev-allH-10 | 8 | v2:M0 | 543 | 0.43 | 791.9 | 0.486735 |
| dev-allH-10 | 8 | v2:M2 | 543 | 0.456 | 839.8 | 0.509755 |
| dev-fixed-12 | 4 | kernel:B | 9 | 0.000111 | 12.3333 | 0.00377498 |
| dev-fixed-12 | 4 | v1:B | 9 | 0.002 | 222.2 | 0.100163 |
| dev-fixed-12 | 4 | v1:T | 9 | 0.495 | 5.5e+04 | 0.581522 |
| dev-fixed-12 | 4 | v2:M0 | 9 | 0.548 | 6.089e+04 | 0.647392 |
| dev-fixed-12 | 4 | v2:M2 | 9 | 1.473 | 1.637e+05 | 1.7784 |
| dev-fixed-12 | 6 | kernel:B | 71 | 0.000349 | 4.91549 | 0.00265342 |
| dev-fixed-12 | 6 | v1:B | 71 | 0.009 | 126.8 | 0.0747177 |
| dev-fixed-12 | 6 | v1:T | 71 | 1.177 | 1.658e+04 | 1.25438 |
| dev-fixed-12 | 6 | v2:M0 | 71 | 1.337 | 1.883e+04 | 1.43322 |
| dev-fixed-12 | 6 | v2:M2 | 71 | 2.548 | 3.589e+04 | 2.60929 |
| dev-fixed-12 | 8 | kernel:B | 543 | 0.001734 | 3.19337 | 0.00419907 |
| dev-fixed-12 | 8 | v1:B | 543 | 0.075 | 138.1 | 0.15293 |
| dev-fixed-12 | 8 | v1:T | 543 | 2.48 | 4567 | 2.54668 |
| dev-fixed-12 | 8 | v2:M0 | 543 | 2.835 | 5221 | 2.91983 |
| dev-fixed-12 | 8 | v2:M2 | 543 | 3.491 | 6429 | 3.57206 |

## S4 — Tightness at identical frontiers

Protocol: same frontiers as S3; compare log10(zHi/Z), log10(Z/zLo), p_width.

| instance | d | method | p_width_float | log10(zHi/Z) | log10(Z/zLo) | wall_s |
|---|---|---|---|---|---|---|
| dev-allH-10 | 4 | kernel:B | 1 | 0.896488 | 0.83927 | 0.00241098 |
| dev-allH-10 | 4 | v1:B | 1 | 4.29497 | — | 0.0724612 |
| dev-allH-10 | 4 | v1:T | 1 | 1.09433 | 1.27703 | 0.209277 |
| dev-allH-10 | 4 | v2:M0 | 1 | 0.440289 | 0.83927 | 0.228121 |
| dev-allH-10 | 4 | v2:M2 | 1 | 0.241092 | 0.83927 | 0.350525 |
| dev-allH-10 | 6 | kernel:B | 1 | 0.836192 | 0.612051 | 0.00236377 |
| dev-allH-10 | 6 | v1:B | 1 | 3.46443 | — | 0.0597718 |
| dev-allH-10 | 6 | v1:T | 1 | 0.238897 | 0.899293 | 0.284048 |
| dev-allH-10 | 6 | v2:M0 | 1 | 0.122954 | 0.612051 | 0.329078 |
| dev-allH-10 | 6 | v2:M2 | 1 | 0 | 0.612051 | 0.405283 |
| dev-allH-10 | 8 | kernel:B | 0.937094 | 0.576119 | 0.354266 | 0.00405402 |
| dev-allH-10 | 8 | v1:B | 1 | 2.03411 | — | 0.103393 |
| dev-allH-10 | 8 | v1:T | 0.899164 | 0 | 0.508442 | 0.373606 |
| dev-allH-10 | 8 | v2:M0 | 0.768294 | 0 | 0.354266 | 0.486735 |
| dev-allH-10 | 8 | v2:M2 | 0.768294 | 0 | 0.354266 | 0.509755 |
| dev-fixed-12 | 4 | kernel:B | 1 | 1.66416 | 0.880819 | 0.00377498 |
| dev-fixed-12 | 4 | v1:B | 1 | 4.15955 | — | 0.100163 |
| dev-fixed-12 | 4 | v1:T | 1 | 1.46153 | 1.4636 | 0.581522 |
| dev-fixed-12 | 4 | v2:M0 | 1 | 1.02823 | 0.880819 | 0.647392 |
| dev-fixed-12 | 4 | v2:M2 | 1 | 0.715662 | 0.880819 | 1.7784 |
| dev-fixed-12 | 6 | kernel:B | 0.98931 | 1.57191 | 0.707626 | 0.00265342 |
| dev-fixed-12 | 6 | v1:B | 1 | 3.59124 | — | 0.0747177 |
| dev-fixed-12 | 6 | v1:T | 0.99672 | 0.83187 | 1.12933 | 1.25438 |
| dev-fixed-12 | 6 | v2:M0 | 0.995326 | 0.660057 | 0.707626 | 1.43322 |
| dev-fixed-12 | 6 | v2:M2 | 0.984607 | 0.355891 | 0.707626 | 2.60929 |
| dev-fixed-12 | 8 | kernel:B | 0.946415 | 1.32076 | 0.518942 | 0.00419907 |
| dev-fixed-12 | 8 | v1:B | 1 | 3.22303 | — | 0.15293 |
| dev-fixed-12 | 8 | v1:T | 0.965861 | 0.233555 | 0.833476 | 2.54668 |
| dev-fixed-12 | 8 | v2:M0 | 0.933371 | 0.204429 | 0.518942 | 2.91983 |
| dev-fixed-12 | 8 | v2:M2 | 0.874195 | 0 | 0.518942 | 3.57206 |

## S5 — Accuracy per unit of checker time (headline fairness test)

Protocol: fixed configuration grids per method; for each wall-clock budget in [0.01, 0.1, 1.0, 10.0, 60.0]s, pick the configuration with the smallest certified `p_width` whose **checker** `wall_s_median` is ≤ budget (producer time reported separately, not charged to the budget). Held-out instances are the headline.

### heldout: no completed rows

### dev

| instance | budget_s | best_family | best_width | best_wall_s | kernel_width | kernel_wall_s |
|---|---|---|---|---|---|---|
| dev-allH-10 | 0.01 | kernel | 0 | 0.002347 | 0 | 0.00234659 |
| dev-allH-10 | 0.1 | kernel | 0 | 0.002347 | 0 | 0.00234659 |
| dev-allH-10 | 1.0 | kernel | 0 | 0.002347 | 0 | 0.00234659 |
| dev-allH-10 | 10.0 | kernel | 0 | 0.002347 | 0 | 0.00234659 |
| dev-allH-10 | 60.0 | kernel | 0 | 0.002347 | 0 | 0.00234659 |
| dev-fixed-12 | 0.01 | kernel | 0 | 0.00373 | 0 | 0.00373002 |
| dev-fixed-12 | 0.1 | kernel | 0 | 0.00373 | 0 | 0.00373002 |
| dev-fixed-12 | 1.0 | kernel | 0 | 0.00373 | 0 | 0.00373002 |
| dev-fixed-12 | 10.0 | kernel | 0 | 0.00373 | 0 | 0.00373002 |
| dev-fixed-12 | 60.0 | kernel | 0 | 0.00373 | 0 | 0.00373002 |

Win counts (best width within budget, ties go to first min):

- `kernel`: 10
- `python:adaptive`: 0
- `python:basic`: 0
- `python:lookahead`: 0
- `v1:B`: 0
- `v1:T`: 0
- `v2:M2`: 0
- `v2:X`: 0

Accepted rows missing the exact reference on dev: **0** (should be 0 for Lean checkers; Python is unverified).

## S6 — Non-enumeration / scaling

Protocol: kernel `uniform(d)` with `B` leaves on long / all-H sequences; coverage = conformations / checker_steps. Matched v1/v2 B/M0 where feasible.

| id | n | d | method | cert_bytes | leaves | verify_s | wall_s | conf | steps | coverage |
|---|---|---|---|---|---|---|---|---|---|---|
| ho-long-30 | 30 | 2 | kernel:B | 81 | 1 | 0.00034 | 0.00268426 | 1543066804034 | 1 | 1.543e+12 |
| ho-long-30 | 30 | 2 | v1:B | 86 | 1 | 0.001 | 0.0531528 | 1543066804034 | None | — |
| ho-long-30 | 30 | 2 | v2:M0 | 87 | 1 | — | 302.9 | 1543066804034 | None | — |
| ho-long-30 | 30 | 4 | kernel:B | 97 | 9 | 0.000374 | 0.0029644 | 1543066804034 | 25 | 6.172e+10 |
| ho-long-30 | 30 | 4 | v1:B | 127 | 9 | — | 144.9 | 1543066804034 | None | — |

## S7 — One tree, many queries

Protocol: one hybrid frontier; kernel carries N query lines in one check; foldcheck2 is run once per (q,t) on the matched certificate; compare sum of walls.

| n_queries | method | wall_s | verify_s | per_query |
|---|---|---|---|---|
| 1 | kernel | 0.00326966 | 0.000544 | 0.000544 |
| 1 | v2 | 1.7814 | — | — |
| 4 | kernel | 0.00340425 | 0.000952 | 0.000238 |
| 4 | v2 | 7.81347 | — | — |

## S8 — Footprint

| binary | bytes | ratio vs foldkernel |
|---|---|---|
| foldkernel | 5,209,296 | 1.0× |
| foldcheck | 222,524,568 | 42.7× |
| foldcheck2 | 222,652,512 | 42.7× |

Coordinator-provided build cost (not re-measured): foldcheck+foldcheck2 wall 484s / 25 CPU-min; foldkernel wall 1.4s.

## Critical findings

None. No accepted Lean certificate missed a frozen reference in this run.

## Load / remaining work

- Load average at summary time: [4.3, 5.23, 2.99]
- Skipped entries: 0
- Result JSONL paths: `experiments/kernel/bench/results/S1.jsonl`, `experiments/kernel/bench/results/S2.jsonl`, `experiments/kernel/bench/results/S3.jsonl`, `experiments/kernel/bench/results/S4.jsonl`, `experiments/kernel/bench/results/S5.jsonl`, `experiments/kernel/bench/results/S6.jsonl`, `experiments/kernel/bench/results/S7.jsonl`, `experiments/kernel/bench/results/S8.jsonl`
