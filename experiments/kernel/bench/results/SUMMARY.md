# FoldBounds certified frontier — matched-budget benchmark SUMMARY

Generated: 2026-11-19T02:11:22+00:00

## Machine / binaries

- Host: `cursor` · Linux `6.12.113+` · 4 CPUs
- Load average at harness start: [1.31, 1.78, 2.09]
- Mode: `full` · rows in JSONL: 1190 · wall 3545.7s (started 2026-11-19T01:12:01+00:00, finished 2026-11-19T02:11:08+00:00)
- `kernel`: 5,209,296 bytes · sha256 `b8ba2338921c795e…`
- `v1`: 222,524,568 bytes · sha256 `6d5a4c490cf06397…`
- `v2`: 222,652,512 bytes · sha256 `b2da8e25f1530d1e…`

## Related accuracy evidence

- Independent differential suite (Python reference vs foldkernel): `experiments/kernel/results/differential_summary.json` — `disagreements: 0` across 757 certificates / 12112 endpoints.

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
- S2: foldcheck2/kernel wall speedup over 14 matched (kind,n) pairs: min 20.9× · median 60.4× · max 501.6× (max at mixed/n=14: 4.735s / 0.00944s).
- S2: hpenum/kernel wall speedup over 18 matched (kind,n) pairs: min 4.2× · median 7.0× · max 9.7× (max at allH/n=16: 0.5489s / 0.0563s).
- S3: in 76 matched (instance,depth,alt) comparisons the alternative was slower per leaf than kernel B; largest 162860.2× (v2:M4 at ho-rand-n13-p0.8/d=8: 5.28e+05 vs 3.24 µs).
- S4: kernel B is tightest or tied on p_width in 10 matched frontiers.
- S5 held-out: kernel family wins 25/25 budgeted (instance,budget) slots on certified p_width.
- S6: kernel non-enumerative check covers ~6.97e+19× more conformations than checker steps (ho-long-48, d=2, conf=69729298658567331840, steps=1).
- S7/n_queries=1: foldcheck2 sum-of-runs / kernel multi-query wall = 567.1× (2.272s / 0.004006s).
- S7/n_queries=4: foldcheck2 sum-of-runs / kernel multi-query wall = 1473.3× (9.654s / 0.006553s).
- S7/n_queries=12: foldcheck2 sum-of-runs / kernel multi-query wall = 2398.7× (29.92s / 0.01247s).
- S7/n_queries=16: foldcheck2 sum-of-runs / kernel multi-query wall = 2638.8× (40.36s / 0.0153s).
- S8: foldcheck binary is 43× larger than foldkernel; foldcheck2 is 43× larger.

### (ii) Ties / no accuracy loss
- S1: kernel root-X exact on 284/284 reference queries (disagreements with reference = 0).
- S1: foldcheck2 root-X agreed with kernel on 13 finished comparisons (disagreements = 0; timeouts = 0).
- S5 held-out: every accepted Lean row with a reference encloses it (accuracy loss vs oracle = 0; disagreements = 0).

### (iii) An alternative better
- S4: an alternative leaf kind is strictly tighter on p_width in 6 matched frontiers (transfer/memory DP can beat kernel B per leaf). Example: dev-fixed-12/d=6: v2:M4 = 0.969167 vs kernel B = 0.98931 (wall 3.47467s vs 0.00282466s).

### (iv) Inconclusive / timeouts / skipped
- (none in this run)

## S1 — Exactness (just as much accuracy)

Protocol: root-`X` certificates; kernel must match the frozen exact reference on every endpoint; foldcheck2 root-`X` compared once per instance (q=2, mid threshold, n≤14) with a 60s hard cap. Repeats: median/min wall.

| metric | value |
|---|---|
| kernel rows | 284 |
| kernel exact_match | 284 |
| kernel encloses_reference | 284 |
| kernel failures/timeouts | 0 |
| foldcheck2 agrees_with_kernel | 13 |
| foldcheck2 disagrees | 0 |
| foldcheck2 timeouts | 0 |

## S2 — Exhaustive speed

Protocol: root-`X` (kernel, foldcheck2) and pure-Python `hpenum` on all-H and one mixed family; stop each alternative at first over-cap run.

| kind | n | method | wall_s_median | verify_s | µs/conf | timed_out |
|---|---|---|---|---|---|---|
| allH | 8 | hpenum | 0.0162543 | 0.00015811 | 0.291179 | False |
| allH | 8 | kernel | 0.00237766 | 7e-05 | 0.128913 | False |
| allH | 8 | v2 | 0.0593981 | 0.005 | 9.2081 | False |
| allH | 9 | hpenum | 0.016517 | 0.0004492 | 0.303719 | False |
| allH | 9 | kernel | 0.00306772 | 0.00011 | 0.0743746 | False |
| allH | 9 | v2 | 0.0641873 | 0.017 | 11.4943 | False |
| allH | 10 | hpenum | 0.017575 | 0.00123261 | 0.303075 | False |
| allH | 10 | kernel | 0.00245584 | 0.000201 | 0.0494222 | False |
| allH | 10 | v2 | 0.106225 | 0.051 | 12.54 | False |
| allH | 11 | hpenum | 0.0192477 | 0.00361369 | 0.327772 | False |
| allH | 11 | kernel | 0.00459262 | 0.000438 | 0.0397279 | False |
| allH | 11 | v2 | 0.21255 | 0.162 | 14.6939 | False |
| allH | 12 | hpenum | 0.0262571 | 0.00984075 | 0.327229 | False |
| allH | 12 | kernel | 0.00505274 | 0.001088 | 0.0361786 | False |
| allH | 12 | v2 | 0.566783 | 0.513 | 17.0585 | False |
| allH | 13 | hpenum | 0.0441904 | 0.0268073 | 0.330005 | False |
| allH | 13 | kernel | 0.00623467 | 0.002727 | 0.0335701 | False |
| allH | 13 | v2 | 1.66783 | 1.608 | 19.7949 | False |
| allH | 14 | hpenum | 0.0897616 | 0.0724372 | 0.3287 | False |
| allH | 14 | kernel | 0.0118748 | 0.007285 | 0.0330573 | False |
| allH | 14 | v2 | 5.05446 | 4.976 | 22.5797 | False |
| allH | 15 | hpenum | 0.21269 | 0.195988 | 0.330163 | False |
| allH | 15 | kernel | 0.02226 | 0.0195 | 0.0328498 | False |
| allH | 16 | hpenum | 0.548867 | 0.531019 | 0.331029 | False |
| allH | 16 | kernel | 0.0562984 | 0.052386 | 0.0326566 | False |
| allH | 17 | kernel | 0.143571 | 0.140396 | 0.0325644 | False |
| allH | 18 | kernel | 0.383165 | 0.378908 | 0.0326176 | False |
| allH | 19 | kernel | 1.02251 | 1.01921 | 0.0327042 | False |
| allH | 20 | kernel | 2.76098 | 2.75706 | 0.0329086 | False |
| allH | 21 | kernel | 7.51577 | 7.51199 | 0.0334723 | False |
| mixed | 8 | hpenum | 0.0164462 | 0.000142049 | 0.2616 | False |
| mixed | 8 | kernel | 0.00297416 | 8.6e-05 | 0.158379 | False |
| mixed | 8 | v2 | 0.0721835 | 0.005 | 9.2081 | False |
| mixed | 9 | hpenum | 0.0166366 | 0.000421309 | 0.284861 | False |
| mixed | 9 | kernel | 0.0022686 | 0.000101 | 0.0682894 | False |
| mixed | 9 | v2 | 0.0721305 | 0.015 | 10.142 | False |
| mixed | 10 | hpenum | 0.0169812 | 0.00085526 | 0.210293 | False |
| mixed | 10 | kernel | 0.00279252 | 0.000186 | 0.045734 | False |
| mixed | 10 | v2 | 0.142363 | 0.049 | 12.0482 | False |
| mixed | 11 | hpenum | 0.018472 | 0.00196986 | 0.178672 | False |
| mixed | 11 | kernel | 0.00289101 | 0.000402 | 0.0364626 | False |
| mixed | 11 | v2 | 0.202037 | 0.154 | 13.9683 | False |
| mixed | 12 | hpenum | 0.0240722 | 0.00810714 | 0.269582 | False |
| mixed | 12 | kernel | 0.00362856 | 0.001061 | 0.0352808 | False |
| mixed | 12 | v2 | 0.532945 | 0.485 | 16.1274 | False |
| mixed | 13 | hpenum | 0.0331761 | 0.0164487 | 0.202488 | False |
| mixed | 13 | kernel | 0.00503295 | 0.002424 | 0.0298401 | False |
| mixed | 13 | v2 | 1.54488 | 1.496 | 18.4162 | False |
| mixed | 14 | hpenum | 0.0791212 | 0.0615694 | 0.279384 | False |
| mixed | 14 | kernel | 0.00943996 | 0.006841 | 0.0310425 | False |
| mixed | 14 | v2 | 4.73489 | 4.686 | 21.2638 | False |
| mixed | 15 | hpenum | 0.201764 | 0.185071 | 0.311772 | False |
| mixed | 15 | kernel | 0.0213166 | 0.018526 | 0.031209 | False |
| mixed | 16 | hpenum | 0.364246 | 0.347923 | 0.216889 | False |
| mixed | 16 | kernel | 0.0505605 | 0.047883 | 0.0298495 | False |
| mixed | 17 | kernel | 0.129796 | 0.126662 | 0.0293788 | False |

## S3 — Per-leaf cost at identical frontiers

Protocol: matched `uniform(d)` frontiers; µs/leaf = verify_s / n_leaves.

| instance | d | method | leaves | verify_s | µs/leaf | wall_s_median |
|---|---|---|---|---|---|---|
| dev-fixed-12 | 4 | kernel:B | 9 | 0.000102 | 11.3333 | 0.00233086 |
| dev-fixed-12 | 4 | v1:B | 9 | 0.002 | 222.2 | 0.0533026 |
| dev-fixed-12 | 4 | v1:T | 9 | 0.434 | 4.822e+04 | 0.481652 |
| dev-fixed-12 | 4 | v2:M0 | 9 | 0.462 | 5.133e+04 | 0.512189 |
| dev-fixed-12 | 4 | v2:M2 | 9 | 1.425 | 1.583e+05 | 1.49855 |
| dev-fixed-12 | 4 | v2:M4 | 9 | 2.429 | 2.699e+05 | 2.47753 |
| dev-fixed-12 | 6 | kernel:B | 71 | 0.000289 | 4.07042 | 0.00282466 |
| dev-fixed-12 | 6 | v1:B | 71 | 0.009 | 126.8 | 0.0569288 |
| dev-fixed-12 | 6 | v1:T | 71 | 1.171 | 1.649e+04 | 1.2192 |
| dev-fixed-12 | 6 | v2:M0 | 71 | 1.32 | 1.859e+04 | 1.37364 |
| dev-fixed-12 | 6 | v2:M2 | 71 | 2.391 | 3.368e+04 | 2.43796 |
| dev-fixed-12 | 6 | v2:M4 | 71 | 3.422 | 4.82e+04 | 3.47467 |
| dev-fixed-12 | 8 | kernel:B | 543 | 0.001712 | 3.15285 | 0.00402594 |
| dev-fixed-12 | 8 | v1:B | 543 | 0.069 | 127.1 | 0.118396 |
| dev-fixed-12 | 8 | v1:T | 543 | 2.171 | 3998 | 2.21879 |
| dev-fixed-12 | 8 | v2:M0 | 543 | 2.696 | 4965 | 2.74726 |
| dev-fixed-12 | 8 | v2:M2 | 543 | 3.413 | 6285 | 3.46297 |
| dev-fixed-12 | 8 | v2:M4 | 543 | 3.955 | 7284 | 4.00401 |
| dev-fixed-12 | 10 | kernel:B | 4067 | 0.012474 | 3.06713 | 0.0151664 |
| dev-fixed-12 | 10 | v1:B | 4067 | 0.522 | 128.4 | 0.576783 |
| dev-fixed-12 | 10 | v1:T | 4067 | 3.004 | 738.6 | 3.05405 |
| dev-fixed-12 | 10 | v2:M0 | 4067 | 4.16 | 1023 | 4.21609 |
| dev-fixed-12 | 10 | v2:M2 | 4067 | 4.485 | 1103 | 4.53927 |
| dev-fixed-12 | 10 | v2:M4 | 4067 | 4.487 | 1103 | 4.54358 |
| dev-rand-n14-p0.7 | 4 | kernel:B | 9 | 0.000118 | 13.1111 | 0.00322334 |
| dev-rand-n14-p0.7 | 4 | v1:B | 9 | 0.002 | 222.2 | 0.0509526 |
| dev-rand-n14-p0.7 | 4 | v1:T | 9 | 1.15 | 1.278e+05 | 1.2138 |
| dev-rand-n14-p0.7 | 4 | v2:M0 | 9 | 1.196 | 1.329e+05 | 1.26394 |
| dev-rand-n14-p0.7 | 4 | v2:M2 | 9 | 8.239 | 9.154e+05 | 8.29102 |
| dev-rand-n14-p0.7 | 4 | v2:M4 | 9 | 14.447 | 1.605e+06 | 14.5105 |
| dev-rand-n14-p0.7 | 6 | kernel:B | 71 | 0.000304 | 4.28169 | 0.00257886 |
| dev-rand-n14-p0.7 | 6 | v1:B | 71 | 0.013 | 183.1 | 0.0798768 |
| dev-rand-n14-p0.7 | 6 | v1:T | 71 | 5.446 | 7.67e+04 | 5.51108 |
| dev-rand-n14-p0.7 | 6 | v2:M0 | 71 | 5.814 | 8.189e+04 | 6.02628 |
| dev-rand-n14-p0.7 | 6 | v2:M2 | 71 | 44.691 | 6.295e+05 | 44.7584 |
| dev-rand-n14-p0.7 | 6 | v2:M4 | 71 | 24.657 | 3.473e+05 | 33.3773 |
| dev-rand-n14-p0.7 | 8 | kernel:B | 543 | 0.001786 | 3.28913 | 0.00413254 |
| dev-rand-n14-p0.7 | 8 | v1:B | 543 | 0.093 | 171.3 | 0.14819 |
| dev-rand-n14-p0.7 | 8 | v1:T | 543 | 13.061 | 2.405e+04 | 13.1103 |
| dev-rand-n14-p0.7 | 8 | v2:M0 | 543 | 15.486 | 2.852e+04 | 15.6312 |
| dev-rand-n14-p0.7 | 8 | v2:M2 | 543 | 29.556 | 5.443e+04 | 29.6216 |
| dev-rand-n14-p0.7 | 8 | v2:M4 | 543 | 34.905 | 6.428e+04 | 34.9558 |
| dev-rand-n14-p0.7 | 10 | kernel:B | 4067 | 0.013084 | 3.21711 | 0.0156898 |
| dev-rand-n14-p0.7 | 10 | v1:B | 4067 | 0.712 | 175.1 | 0.77808 |
| ho-rand-n13-p0.6 | 4 | kernel:B | 9 | 0.000107 | 11.8889 | 0.00327273 |
| ho-rand-n13-p0.6 | 4 | v1:B | 9 | 0.001 | 111.1 | 0.0504017 |
| ho-rand-n13-p0.6 | 4 | v1:T | 9 | 0.752 | 8.356e+04 | 0.800176 |
| ho-rand-n13-p0.6 | 4 | v2:M0 | 9 | 0.801 | 8.9e+04 | 0.850867 |
| ho-rand-n13-p0.6 | 4 | v2:M2 | 9 | 3.946 | 4.384e+05 | 3.99608 |
| ho-rand-n13-p0.6 | 4 | v2:M4 | 9 | 6.229 | 6.921e+05 | 6.2822 |
| ho-rand-n13-p0.6 | 6 | kernel:B | 71 | 0.000296 | 4.16901 | 0.00263933 |
| ho-rand-n13-p0.6 | 6 | v1:B | 71 | 0.01 | 140.8 | 0.0580233 |
| ho-rand-n13-p0.6 | 6 | v1:T | 71 | 2.642 | 3.721e+04 | 2.68889 |
| ho-rand-n13-p0.6 | 6 | v2:M0 | 71 | 2.93 | 4.127e+04 | 2.98173 |
| ho-rand-n13-p0.6 | 6 | v2:M2 | 71 | 6.87 | 9.676e+04 | 6.91966 |
| ho-rand-n13-p0.6 | 6 | v2:M4 | 71 | 9.256 | 1.304e+05 | 9.30868 |
| ho-rand-n13-p0.6 | 8 | kernel:B | 543 | 0.001735 | 3.19521 | 0.00410375 |
| ho-rand-n13-p0.6 | 8 | v1:B | 543 | 0.081 | 149.2 | 0.138603 |
| ho-rand-n13-p0.6 | 8 | v1:T | 543 | 5.445 | 1.003e+04 | 5.49527 |
| ho-rand-n13-p0.6 | 8 | v2:M0 | 543 | 6.496 | 1.196e+04 | 6.56303 |
| ho-rand-n13-p0.6 | 8 | v2:M2 | 543 | 9.542 | 1.757e+04 | 9.60969 |
| ho-rand-n13-p0.6 | 8 | v2:M4 | 543 | 12.289 | 2.263e+04 | 12.3449 |
| ho-rand-n13-p0.6 | 10 | kernel:B | 4067 | 0.012723 | 3.12835 | 0.0156098 |
| ho-rand-n13-p0.6 | 10 | v1:B | 4067 | 0.606 | 149 | 0.659325 |
| ho-rand-n13-p0.6 | 10 | v1:T | 4067 | 7.633 | 1877 | 7.68453 |
| ho-rand-n13-p0.6 | 10 | v2:M0 | 4067 | 10.381 | 2552 | 10.4405 |
| ho-rand-n13-p0.6 | 10 | v2:M2 | 4067 | 12.658 | 3112 | 12.7158 |
| ho-rand-n13-p0.6 | 10 | v2:M4 | 4067 | 12.991 | 3194 | 13.0506 |
| ho-rand-n13-p0.8 | 4 | kernel:B | 9 | 0.000109 | 12.1111 | 0.00234732 |
| ho-rand-n13-p0.8 | 4 | v1:B | 9 | 0.002 | 222.2 | 0.0531196 |
| ho-rand-n13-p0.8 | 4 | v1:T | 9 | 0.793 | 8.811e+04 | 0.903074 |
| ho-rand-n13-p0.8 | 4 | v2:M0 | 9 | 0.84 | 9.333e+04 | 0.932374 |
| ho-rand-n13-p0.8 | 4 | v2:M2 | 9 | 4.431 | 4.923e+05 | 4.49314 |
| ho-rand-n13-p0.8 | 4 | v2:M4 | 9 | 7.04 | 7.822e+05 | 7.10282 |
| ho-rand-n13-p0.8 | 6 | kernel:B | 71 | 0.0003 | 4.22535 | 0.00259868 |
| ho-rand-n13-p0.8 | 6 | v1:B | 71 | 0.011 | 154.9 | 0.0678123 |
| ho-rand-n13-p0.8 | 6 | v1:T | 71 | 2.918 | 4.11e+04 | 2.9713 |
| ho-rand-n13-p0.8 | 6 | v2:M0 | 71 | 3.207 | 4.517e+04 | 3.27311 |
| ho-rand-n13-p0.8 | 6 | v2:M2 | 71 | 7.896 | 1.112e+05 | 7.96448 |
| ho-rand-n13-p0.8 | 6 | v2:M4 | 71 | 10.096 | 1.422e+05 | 10.1525 |
| ho-rand-n13-p0.8 | 8 | kernel:B | 543 | 0.00176 | 3.24125 | 0.00418662 |
| ho-rand-n13-p0.8 | 8 | v1:B | 543 | 0.082 | 151 | 0.129184 |
| ho-rand-n13-p0.8 | 8 | v1:T | 543 | 6.151 | 1.133e+04 | 6.19743 |
| ho-rand-n13-p0.8 | 8 | v2:M0 | 543 | 7.301 | 1.345e+04 | 7.35685 |
| ho-rand-n13-p0.8 | 8 | v2:M2 | 543 | 11.047 | 2.034e+04 | 11.1086 |
| ho-rand-n13-p0.8 | 8 | v2:M4 | 543 | 286.6 | 5.279e+05 | 286.7 |
| ho-rand-n13-p0.8 | 10 | kernel:B | 4067 | 0.012849 | 3.15933 | 0.015375 |
| ho-rand-n13-p0.8 | 10 | v1:B | 4067 | 0.624 | 153.4 | 0.685521 |
| ho-rand-n13-p0.8 | 10 | v1:T | 4067 | 8.536 | 2099 | 8.58866 |
| ho-rand-n13-p0.8 | 10 | v2:M0 | 4067 | 11.273 | 2772 | 11.4344 |
| ho-rand-n13-p0.8 | 10 | v2:M2 | 4067 | 14 | 3442 | 14.0603 |
| ho-rand-n13-p0.8 | 10 | v2:M4 | 4067 | 13.361 | 3285 | 13.4222 |

## S4 — Tightness at identical frontiers

Protocol: same frontiers as S3; compare log10(zHi/Z), log10(Z/zLo), p_width.

| instance | d | method | p_width_float | log10(zHi/Z) | log10(Z/zLo) | wall_s |
|---|---|---|---|---|---|---|
| dev-fixed-12 | 4 | kernel:B | 1 | 1.66416 | 0.880819 | 0.00233086 |
| dev-fixed-12 | 4 | v1:B | 1 | 4.15955 | — | 0.0533026 |
| dev-fixed-12 | 4 | v1:T | 1 | 1.46153 | 1.4636 | 0.481652 |
| dev-fixed-12 | 4 | v2:M0 | 1 | 1.02823 | 0.880819 | 0.512189 |
| dev-fixed-12 | 4 | v2:M2 | 1 | 0.715662 | 0.880819 | 1.49855 |
| dev-fixed-12 | 4 | v2:M4 | 1 | 0.0820962 | 0.880819 | 2.47753 |
| dev-fixed-12 | 6 | kernel:B | 0.98931 | 1.57191 | 0.707626 | 0.00282466 |
| dev-fixed-12 | 6 | v1:B | 1 | 3.59124 | — | 0.0569288 |
| dev-fixed-12 | 6 | v1:T | 0.99672 | 0.83187 | 1.12933 | 1.2192 |
| dev-fixed-12 | 6 | v2:M0 | 0.995326 | 0.660057 | 0.707626 | 1.37364 |
| dev-fixed-12 | 6 | v2:M2 | 0.984607 | 0.355891 | 0.707626 | 2.43796 |
| dev-fixed-12 | 6 | v2:M4 | 0.969167 | 0 | 0.707626 | 3.47467 |
| dev-fixed-12 | 8 | kernel:B | 0.946415 | 1.32076 | 0.518942 | 0.00402594 |
| dev-fixed-12 | 8 | v1:B | 1 | 3.22303 | — | 0.118396 |
| dev-fixed-12 | 8 | v1:T | 0.965861 | 0.233555 | 0.833476 | 2.21879 |
| dev-fixed-12 | 8 | v2:M0 | 0.933371 | 0.204429 | 0.518942 | 2.74726 |
| dev-fixed-12 | 8 | v2:M2 | 0.874195 | 0 | 0.518942 | 3.46297 |
| dev-fixed-12 | 8 | v2:M4 | 0.874195 | 0 | 0.518942 | 4.00401 |
| dev-fixed-12 | 10 | kernel:B | 0.804414 | 0.943716 | 0.315026 | 0.0151664 |
| dev-fixed-12 | 10 | v1:B | 1 | 1.79211 | — | 0.576783 |
| dev-fixed-12 | 10 | v1:T | 0.752965 | 0 | 0.476946 | 3.05405 |
| dev-fixed-12 | 10 | v2:M0 | 0.691116 | 0 | 0.315026 | 4.21609 |
| dev-fixed-12 | 10 | v2:M2 | 0.691116 | 0 | 0.315026 | 4.53927 |
| dev-fixed-12 | 10 | v2:M4 | 0.691116 | 0 | 0.315026 | 4.54358 |
| dev-rand-n14-p0.7 | 4 | kernel:B | 1 | 1.94048 | 1.10736 | 0.00322334 |
| dev-rand-n14-p0.7 | 4 | v1:B | 1 | 6.24205 | — | 0.0509526 |
| dev-rand-n14-p0.7 | 4 | v1:T | 1 | 3.16579 | 1.84049 | 1.2138 |
| dev-rand-n14-p0.7 | 4 | v2:M0 | 1 | 1.39926 | 1.10736 | 1.26394 |
| dev-rand-n14-p0.7 | 4 | v2:M2 | 1 | 1.29711 | 1.10736 | 8.29102 |
| dev-rand-n14-p0.7 | 4 | v2:M4 | 1 | 0.876121 | 1.10736 | 14.5105 |
| dev-rand-n14-p0.7 | 6 | kernel:B | 1 | 1.88019 | 0.954818 | 0.00257886 |
| dev-rand-n14-p0.7 | 6 | v1:B | 1 | 5.90443 | — | 0.0798768 |
| dev-rand-n14-p0.7 | 6 | v1:T | 1 | 1.99698 | 1.55452 | 5.51108 |
| dev-rand-n14-p0.7 | 6 | v2:M0 | 1 | 1.2113 | 0.954818 | 6.02628 |
| dev-rand-n14-p0.7 | 6 | v2:M2 | 1 | 0.910193 | 0.954818 | 44.7584 |
| dev-rand-n14-p0.7 | 6 | v2:M4 | 1 | 0.317115 | 0.954818 | 33.3773 |
| dev-rand-n14-p0.7 | 8 | kernel:B | 1 | 1.74263 | 0.795532 | 0.00413254 |
| dev-rand-n14-p0.7 | 8 | v1:B | 1 | 5.56227 | — | 0.14819 |
| dev-rand-n14-p0.7 | 8 | v1:T | 1 | 1.0708 | 1.2567 | 13.1103 |
| dev-rand-n14-p0.7 | 8 | v2:M0 | 1 | 0.819599 | 0.795532 | 15.6312 |
| dev-rand-n14-p0.7 | 8 | v2:M2 | 1 | 0.378921 | 0.795532 | 29.6216 |
| dev-rand-n14-p0.7 | 8 | v2:M4 | 1 | 0 | 0.795532 | 34.9558 |
| dev-rand-n14-p0.7 | 10 | kernel:B | 0.992796 | 1.4655 | 0.606747 | 0.0156898 |
| dev-rand-n14-p0.7 | 10 | v1:B | 1 | 4.14236 | — | 0.77808 |
| ho-rand-n13-p0.6 | 4 | kernel:B | 1 | 1.80629 | 0.857909 | 0.00327273 |
| ho-rand-n13-p0.6 | 4 | v1:B | 1 | 3.9135 | — | 0.0504017 |
| ho-rand-n13-p0.6 | 4 | v1:T | 1 | 1.7546 | 1.58471 | 0.800176 |
| ho-rand-n13-p0.6 | 4 | v2:M0 | 1 | 1.31525 | 0.857909 | 0.850867 |
| ho-rand-n13-p0.6 | 4 | v2:M2 | 1 | 0.951258 | 0.857909 | 3.99608 |
| ho-rand-n13-p0.6 | 4 | v2:M4 | 1 | 0.33247 | 0.857909 | 6.2822 |
| ho-rand-n13-p0.6 | 6 | kernel:B | 1 | 1.74599 | 0.748525 | 0.00263933 |
| ho-rand-n13-p0.6 | 6 | v1:B | 1 | 3.62306 | — | 0.0580233 |
| ho-rand-n13-p0.6 | 6 | v1:T | 1 | 0.880403 | 1.26595 | 2.68889 |
| ho-rand-n13-p0.6 | 6 | v2:M0 | 1 | 0.772367 | 0.748525 | 2.98173 |
| ho-rand-n13-p0.6 | 6 | v2:M2 | 1 | 0.412836 | 0.748525 | 6.91966 |
| ho-rand-n13-p0.6 | 6 | v2:M4 | 1 | 0.0216488 | 0.748525 | 9.30868 |
| ho-rand-n13-p0.6 | 8 | kernel:B | 0.988567 | 1.60844 | 0.555124 | 0.00410375 |
| ho-rand-n13-p0.6 | 8 | v1:B | 1 | 2.7652 | — | 0.138603 |
| ho-rand-n13-p0.6 | 8 | v1:T | 0.9862 | 0.282882 | 0.909929 | 5.49527 |
| ho-rand-n13-p0.6 | 8 | v2:M0 | 0.985814 | 0.268763 | 0.555124 | 6.56303 |
| ho-rand-n13-p0.6 | 8 | v2:M2 | 0.971961 | 0.0524717 | 0.555124 | 9.60969 |
| ho-rand-n13-p0.6 | 8 | v2:M4 | 0.970639 | 0 | 0.555124 | 12.3449 |
| ho-rand-n13-p0.6 | 10 | kernel:B | 0.866471 | 0.710458 | 0.31916 | 0.0156098 |
| ho-rand-n13-p0.6 | 10 | v1:B | 1 | 1.31252 | — | 0.659325 |
| ho-rand-n13-p0.6 | 10 | v1:T | 0.832819 | 0 | 0.538075 | 7.68453 |
| ho-rand-n13-p0.6 | 10 | v2:M0 | 0.771836 | 0 | 0.31916 | 10.4405 |
| ho-rand-n13-p0.6 | 10 | v2:M2 | 0.771836 | 0 | 0.31916 | 12.7158 |
| ho-rand-n13-p0.6 | 10 | v2:M4 | 0.771836 | 0 | 0.31916 | 13.0506 |
| ho-rand-n13-p0.8 | 4 | kernel:B | 1 | 1.60482 | 1.05938 | 0.00234732 |
| ho-rand-n13-p0.8 | 4 | v1:B | 1 | 6.12027 | — | 0.0531196 |
| ho-rand-n13-p0.8 | 4 | v1:T | 1 | 3.11703 | 1.78618 | 0.903074 |
| ho-rand-n13-p0.8 | 4 | v2:M0 | 1 | 1.11378 | 1.05938 | 0.932374 |
| ho-rand-n13-p0.8 | 4 | v2:M2 | 1 | 1.01842 | 1.05938 | 4.49314 |
| ho-rand-n13-p0.8 | 4 | v2:M4 | 1 | 0.526501 | 1.05938 | 7.10282 |
| ho-rand-n13-p0.8 | 6 | kernel:B | 1 | 1.54453 | 0.914462 | 0.00259868 |
| ho-rand-n13-p0.8 | 6 | v1:B | 1 | 5.80838 | — | 0.0678123 |
| ho-rand-n13-p0.8 | 6 | v1:T | 1 | 1.69072 | 1.46742 | 2.9713 |
| ho-rand-n13-p0.8 | 6 | v2:M0 | 1 | 0.944279 | 0.914462 | 3.27311 |
| ho-rand-n13-p0.8 | 6 | v2:M2 | 1 | 0.64796 | 0.914462 | 7.96448 |
| ho-rand-n13-p0.8 | 6 | v2:M4 | 1 | 0.0581294 | 0.914462 | 10.1525 |
| ho-rand-n13-p0.8 | 8 | kernel:B | 1 | 1.40697 | 0.739838 | 0.00418662 |
| ho-rand-n13-p0.8 | 8 | v1:B | 1 | 4.95379 | — | 0.129184 |
| ho-rand-n13-p0.8 | 8 | v1:T | 1 | 0.542735 | 1.11012 | 6.19743 |
| ho-rand-n13-p0.8 | 8 | v2:M0 | 1 | 0.454776 | 0.739838 | 7.35685 |
| ho-rand-n13-p0.8 | 8 | v2:M2 | 1 | 0.10958 | 0.739838 | 11.1086 |
| ho-rand-n13-p0.8 | 8 | v2:M4 | 1 | 0 | 0.739838 | 286.7 |
| ho-rand-n13-p0.8 | 10 | kernel:B | 0.981026 | 1.12984 | 0.5123 | 0.015375 |
| ho-rand-n13-p0.8 | 10 | v1:B | 1 | 3.5045 | — | 0.685521 |
| ho-rand-n13-p0.8 | 10 | v1:T | 0.966539 | 0 | 0.738054 | 8.58866 |
| ho-rand-n13-p0.8 | 10 | v2:M0 | 0.946904 | 0 | 0.5123 | 11.4344 |
| ho-rand-n13-p0.8 | 10 | v2:M2 | 0.946904 | 0 | 0.5123 | 14.0603 |
| ho-rand-n13-p0.8 | 10 | v2:M4 | 0.946904 | 0 | 0.5123 | 13.4222 |

## S5 — Accuracy per unit of checker time (headline fairness test)

Protocol: fixed configuration grids per method; for each wall-clock budget in [0.01, 0.1, 1.0, 10.0, 60.0]s, pick the configuration with the smallest certified `p_width` whose **checker** `wall_s_median` is ≤ budget (producer time reported separately, not charged to the budget). Held-out instances are the headline.

### heldout

| instance | budget_s | best_family | best_width | best_wall_s | kernel_width | kernel_wall_s |
|---|---|---|---|---|---|---|
| ho-rand-n13-p0.6 | 0.01 | kernel | 0 | 0.004837 | 0 | 0.00483734 |
| ho-rand-n13-p0.6 | 0.1 | kernel | 0 | 0.004837 | 0 | 0.00483734 |
| ho-rand-n13-p0.6 | 1.0 | kernel | 0 | 0.004837 | 0 | 0.00483734 |
| ho-rand-n13-p0.6 | 10.0 | kernel | 0 | 0.004837 | 0 | 0.00483734 |
| ho-rand-n13-p0.6 | 60.0 | kernel | 0 | 0.004837 | 0 | 0.00483734 |
| ho-rand-n13-p0.8 | 0.01 | kernel | 0 | 0.004987 | 0 | 0.00498655 |
| ho-rand-n13-p0.8 | 0.1 | kernel | 0 | 0.004987 | 0 | 0.00498655 |
| ho-rand-n13-p0.8 | 1.0 | kernel | 0 | 0.004987 | 0 | 0.00498655 |
| ho-rand-n13-p0.8 | 10.0 | kernel | 0 | 0.004987 | 0 | 0.00498655 |
| ho-rand-n13-p0.8 | 60.0 | kernel | 0 | 0.004987 | 0 | 0.00498655 |
| ho-rand-n15-p0.6 | 0.01 | kernel | 0.824199 | 0.007808 | 0.824199 | 0.00780843 |
| ho-rand-n15-p0.6 | 0.1 | kernel | 0 | 0.02153 | 0 | 0.0215327 |
| ho-rand-n15-p0.6 | 1.0 | kernel | 0 | 0.02153 | 0 | 0.0215327 |
| ho-rand-n15-p0.6 | 10.0 | kernel | 0 | 0.02153 | 0 | 0.0215327 |
| ho-rand-n15-p0.6 | 60.0 | kernel | 0 | 0.02153 | 0 | 0.0215327 |
| ho-rand-n15-p0.8 | 0.01 | kernel | 0.858153 | 0.009748 | 0.858153 | 0.00974792 |
| ho-rand-n15-p0.8 | 0.1 | kernel | 0 | 0.02152 | 0 | 0.0215213 |
| ho-rand-n15-p0.8 | 1.0 | kernel | 0 | 0.02152 | 0 | 0.0215213 |
| ho-rand-n15-p0.8 | 10.0 | kernel | 0 | 0.02152 | 0 | 0.0215213 |
| ho-rand-n15-p0.8 | 60.0 | kernel | 0 | 0.02152 | 0 | 0.0215213 |
| ho-rand-n17-p0.6 | 0.01 | kernel | 0.96618 | 0.008945 | 0.96618 | 0.00894518 |
| ho-rand-n17-p0.6 | 0.1 | kernel | 0.68116 | 0.06719 | 0.68116 | 0.0671938 |
| ho-rand-n17-p0.6 | 1.0 | kernel | 0 | 0.1342 | 0 | 0.134177 |
| ho-rand-n17-p0.6 | 10.0 | kernel | 0 | 0.1342 | 0 | 0.134177 |
| ho-rand-n17-p0.6 | 60.0 | kernel | 0 | 0.1342 | 0 | 0.134177 |

Win counts (best width within budget, ties go to first min):

- `kernel`: 25
- `python:adaptive`: 0
- `python:basic`: 0
- `python:lookahead`: 0
- `v1:B`: 0
- `v1:T`: 0
- `v2:M0`: 0
- `v2:M2`: 0
- `v2:X`: 0

Accepted Lean rows missing the exact reference on heldout: **0** (should be 0).

### dev

| instance | budget_s | best_family | best_width | best_wall_s | kernel_width | kernel_wall_s |
|---|---|---|---|---|---|---|
| dev-fixed-12 | 0.01 | kernel | 0 | 0.003479 | 0 | 0.00347887 |
| dev-fixed-12 | 0.1 | kernel | 0 | 0.003479 | 0 | 0.00347887 |
| dev-fixed-12 | 1.0 | kernel | 0 | 0.003479 | 0 | 0.00347887 |
| dev-fixed-12 | 10.0 | kernel | 0 | 0.003479 | 0 | 0.00347887 |
| dev-fixed-12 | 60.0 | kernel | 0 | 0.003479 | 0 | 0.00347887 |
| dev-rand-n14-p0.7 | 0.01 | kernel | 0 | 0.009879 | 0 | 0.00987874 |
| dev-rand-n14-p0.7 | 0.1 | kernel | 0 | 0.009879 | 0 | 0.00987874 |
| dev-rand-n14-p0.7 | 1.0 | kernel | 0 | 0.009879 | 0 | 0.00987874 |
| dev-rand-n14-p0.7 | 10.0 | kernel | 0 | 0.009879 | 0 | 0.00987874 |
| dev-rand-n14-p0.7 | 60.0 | kernel | 0 | 0.009879 | 0 | 0.00987874 |
| dev-rand-n16-p0.7 | 0.01 | kernel | 0.958922 | 0.006252 | 0.958922 | 0.00625189 |
| dev-rand-n16-p0.7 | 0.1 | kernel | 0 | 0.05292 | 0 | 0.0529209 |
| dev-rand-n16-p0.7 | 1.0 | kernel | 0 | 0.05292 | 0 | 0.0529209 |
| dev-rand-n16-p0.7 | 10.0 | kernel | 0 | 0.05292 | 0 | 0.0529209 |
| dev-rand-n16-p0.7 | 60.0 | kernel | 0 | 0.05292 | 0 | 0.0529209 |

Win counts (best width within budget, ties go to first min):

- `kernel`: 15
- `python:adaptive`: 0
- `python:basic`: 0
- `python:lookahead`: 0
- `v1:B`: 0
- `v1:T`: 0
- `v2:M0`: 0
- `v2:M2`: 0
- `v2:X`: 0

Accepted Lean rows missing the exact reference on dev: **0** (should be 0).

## S6 — Non-enumeration / scaling

Protocol: kernel `uniform(d)` with `B` leaves on long / all-H sequences; coverage = conformations / checker_steps. Matched v1/v2 B/M0 where feasible.

| id | n | d | method | cert_bytes | leaves | verify_s | wall_s | conf | steps | coverage |
|---|---|---|---|---|---|---|---|---|---|---|
| allH-18 | 18 | 2 | kernel:B | 69 | 1 | 0.000128 | 0.00242072 | 11616669 | 1 | 1.162e+07 |
| allH-18 | 18 | 2 | v1:B | 74 | 1 | 0 | 0.0483133 | 11616669 | None | — |
| allH-18 | 18 | 2 | v2:M0 | 75 | 1 | 1.137 | 1.18653 | 11616669 | None | — |
| allH-18 | 18 | 4 | kernel:B | 85 | 9 | 0.000163 | 0.00222187 | 11616669 | 25 | 4.647e+05 |
| allH-18 | 18 | 4 | v1:B | 115 | 9 | 0.003 | 0.0502901 | 11616669 | None | — |
| allH-18 | 18 | 4 | v2:M0 | 124 | 9 | 7.972 | 8.021 | 11616669 | None | — |
| allH-18 | 18 | 6 | kernel:B | 221 | 71 | 0.000389 | 0.00284926 | 11616669 | 223 | 5.209e+04 |
| allH-18 | 18 | 6 | v1:B | 567 | 71 | 0.02 | 0.067818 | 11616669 | None | — |
| allH-18 | 18 | 6 | v2:M0 | 638 | 71 | — | 30.0036 | 11616669 | None | — |
| allH-18 | 18 | 8 | kernel:B | 1285 | 543 | 0.002019 | 0.00446488 | 11616669 | 1759 | 6604 |
| allH-18 | 18 | 10 | kernel:B | 9373 | 4067 | 0.014798 | 0.0171491 | 11616669 | 13371 | 868.8 |
| allH-18 | 18 | 12 | kernel:B | 69741 | 30073 | 0.108022 | 0.111246 | 11616669 | 99745 | 116.5 |
| allH-19 | 19 | 2 | kernel:B | 70 | 1 | 0.000139 | 0.00238872 | 31164683 | 1 | 3.116e+07 |
| allH-19 | 19 | 2 | v1:B | 75 | 1 | 0 | 0.0485692 | 31164683 | None | — |
| allH-19 | 19 | 2 | v2:M0 | 76 | 1 | 1.507 | 1.55953 | 31164683 | None | — |
| allH-19 | 19 | 4 | kernel:B | 86 | 9 | 0.000176 | 0.00262281 | 31164683 | 25 | 1.247e+06 |
| allH-19 | 19 | 4 | v1:B | 116 | 9 | 0.003 | 0.0502208 | 31164683 | None | — |
| allH-19 | 19 | 4 | v2:M0 | 125 | 9 | 10.921 | 10.9722 | 31164683 | None | — |
| allH-19 | 19 | 6 | kernel:B | 222 | 71 | 0.000413 | 0.00281342 | 31164683 | 223 | 1.398e+05 |
| allH-19 | 19 | 6 | v1:B | 568 | 71 | 0.022 | 0.0707011 | 31164683 | None | — |
| allH-19 | 19 | 6 | v2:M0 | 639 | 71 | — | 30.0037 | 31164683 | None | — |
| allH-19 | 19 | 8 | kernel:B | 1286 | 543 | 0.002036 | 0.00459754 | 31164683 | 1759 | 1.772e+04 |
| allH-19 | 19 | 10 | kernel:B | 9374 | 4067 | 0.014809 | 0.0172873 | 31164683 | 13371 | 2331 |
| allH-19 | 19 | 12 | kernel:B | 69742 | 30073 | 0.110286 | 0.113049 | 31164683 | 99745 | 312.4 |
| allH-20 | 20 | 2 | kernel:B | 71 | 1 | 0.000149 | 0.00226553 | 83779155 | 1 | 8.378e+07 |
| allH-20 | 20 | 2 | v1:B | 76 | 1 | 0 | 0.0478093 | 83779155 | None | — |
| allH-20 | 20 | 2 | v2:M0 | 77 | 1 | 2.138 | 2.19141 | 83779155 | None | — |
| allH-20 | 20 | 4 | kernel:B | 87 | 9 | 0.000187 | 0.00251569 | 83779155 | 25 | 3.351e+06 |
| allH-20 | 20 | 4 | v1:B | 117 | 9 | 0.004 | 0.0508271 | 83779155 | None | — |
| allH-20 | 20 | 4 | v2:M0 | 126 | 9 | 14.981 | 15.0313 | 83779155 | None | — |
| allH-20 | 20 | 6 | kernel:B | 223 | 71 | 0.000437 | 0.00300844 | 83779155 | 223 | 3.757e+05 |
| allH-20 | 20 | 6 | v1:B | 569 | 71 | 0.024 | 0.0747903 | 83779155 | None | — |
| allH-20 | 20 | 6 | v2:M0 | 640 | 71 | — | 30.0035 | 83779155 | None | — |
| allH-20 | 20 | 8 | kernel:B | 1287 | 543 | 0.002255 | 0.00472579 | 83779155 | 1759 | 4.763e+04 |
| allH-20 | 20 | 10 | kernel:B | 9375 | 4067 | 0.01543 | 0.0180694 | 83779155 | 13371 | 6266 |
| allH-20 | 20 | 12 | kernel:B | 69743 | 30073 | 0.11241 | 0.118608 | 83779155 | 99745 | 839.9 |
| allH-21 | 21 | 2 | kernel:B | 72 | 1 | 0.000165 | 0.00228981 | 224424291 | 1 | 2.244e+08 |
| allH-21 | 21 | 2 | v1:B | 77 | 1 | 0 | 0.0518496 | 224424291 | None | — |
| allH-21 | 21 | 2 | v2:M0 | 78 | 1 | 2.884 | 2.93582 | 224424291 | None | — |
| allH-21 | 21 | 4 | kernel:B | 88 | 9 | 0.000217 | 0.00242113 | 224424291 | 25 | 8.977e+06 |
| allH-21 | 21 | 4 | v1:B | 118 | 9 | 0.004 | 0.0506497 | 224424291 | None | — |
| allH-21 | 21 | 4 | v2:M0 | 127 | 9 | 19.705 | 19.7571 | 224424291 | None | — |
| allH-21 | 21 | 6 | kernel:B | 224 | 71 | 0.000451 | 0.00279759 | 224424291 | 223 | 1.006e+06 |
| allH-21 | 21 | 6 | v1:B | 570 | 71 | 0.026 | 0.0742526 | 224424291 | None | — |
| allH-21 | 21 | 6 | v2:M0 | 641 | 71 | — | 30.0036 | 224424291 | None | — |
| allH-21 | 21 | 8 | kernel:B | 1288 | 543 | 0.002333 | 0.00485452 | 224424291 | 1759 | 1.276e+05 |
| allH-21 | 21 | 10 | kernel:B | 9376 | 4067 | 0.01548 | 0.0181425 | 224424291 | 13371 | 1.678e+04 |
| allH-21 | 21 | 12 | kernel:B | 69744 | 30073 | 0.114612 | 0.118287 | 224424291 | 99745 | 2250 |
| allH-22 | 22 | 2 | kernel:B | 73 | 1 | 0.000182 | 0.00242715 | 602201507 | 1 | 6.022e+08 |
| allH-22 | 22 | 2 | v1:B | 78 | 1 | 0 | 0.0479986 | 602201507 | None | — |
| allH-22 | 22 | 2 | v2:M0 | 79 | 1 | 3.379 | 3.43269 | 602201507 | None | — |
| allH-22 | 22 | 4 | kernel:B | 89 | 9 | 0.000219 | 0.00260305 | 602201507 | 25 | 2.409e+07 |
| allH-22 | 22 | 4 | v1:B | 119 | 9 | 0.004 | 0.0528121 | 602201507 | None | — |
| allH-22 | 22 | 4 | v2:M0 | 128 | 9 | 25.625 | 25.6789 | 602201507 | None | — |
| allH-22 | 22 | 6 | kernel:B | 225 | 71 | 0.000471 | 0.00298059 | 602201507 | 223 | 2.7e+06 |
| allH-22 | 22 | 6 | v1:B | 571 | 71 | 0.029 | 0.0771239 | 602201507 | None | — |
| allH-22 | 22 | 6 | v2:M0 | 642 | 71 | — | 30.0037 | 602201507 | None | — |
| allH-22 | 22 | 8 | kernel:B | 1289 | 543 | 0.002363 | 0.00486684 | 602201507 | 1759 | 3.424e+05 |
| allH-22 | 22 | 10 | kernel:B | 9377 | 4067 | 0.017101 | 0.0195698 | 602201507 | 13371 | 4.504e+04 |
| allH-22 | 22 | 12 | kernel:B | 69745 | 30073 | 0.115628 | 0.11863 | 602201507 | 99745 | 6037 |
| ho-long-30 | 30 | 2 | kernel:B | 81 | 1 | 0.000328 | 0.0035633 | 1543066804034 | 1 | 1.543e+12 |
| ho-long-30 | 30 | 4 | kernel:B | 97 | 9 | 0.00037 | 0.00259338 | 1543066804034 | 25 | 6.172e+10 |
| ho-long-30 | 30 | 6 | kernel:B | 233 | 71 | 0.000677 | 0.00296986 | 1543066804034 | 223 | 6.92e+09 |
| ho-long-30 | 30 | 8 | kernel:B | 1297 | 543 | 0.002999 | 0.00515887 | 1543066804034 | 1759 | 8.772e+08 |
| ho-long-30 | 30 | 10 | kernel:B | 9385 | 4067 | 0.020207 | 0.0227732 | 1543066804034 | 13371 | 1.154e+08 |
| ho-long-30 | 30 | 12 | kernel:B | 69753 | 30073 | 0.136624 | 0.139601 | 1543066804034 | 99745 | 1.547e+07 |
| ho-long-48 | 48 | 2 | kernel:B | 99 | 1 | 0.000983 | 0.0032765 | 69729298658567331840 | 1 | 6.973e+19 |
| ho-long-48 | 48 | 4 | kernel:B | 115 | 9 | 0.001006 | 0.00305024 | 69729298658567331840 | 25 | 2.789e+18 |
| ho-long-48 | 48 | 6 | kernel:B | 251 | 71 | 0.001418 | 0.00337138 | 69729298658567331840 | 223 | 3.127e+17 |
| ho-long-48 | 48 | 8 | kernel:B | 1315 | 543 | 0.004272 | 0.00643801 | 69729298658567331840 | 1759 | 3.964e+16 |
| ho-long-48 | 48 | 10 | kernel:B | 9403 | 4067 | 0.026194 | 0.0286983 | 69729298658567331840 | 13371 | 5.215e+15 |
| ho-long-48 | 48 | 12 | kernel:B | 69771 | 30073 | 0.187983 | 0.191183 | 69729298658567331840 | 99745 | 6.991e+14 |

## S7 — One tree, many queries

Protocol: one hybrid frontier; kernel carries N query lines in one check; foldcheck2 is run once per (q,t) on the matched certificate; compare sum of walls.

| n_queries | method | wall_s | verify_s | per_query |
|---|---|---|---|---|
| 1 | kernel | 0.0040057 | 0.001901 | 0.001901 |
| 1 | v2 | 2.27173 | — | — |
| 4 | kernel | 0.00655271 | 0.004215 | 0.00105375 |
| 4 | v2 | 9.65379 | — | — |
| 12 | kernel | 0.0124745 | 0.010025 | 0.000835417 |
| 12 | v2 | 29.9222 | — | — |
| 16 | kernel | 0.0152963 | 0.01278 | 0.00079875 |
| 16 | v2 | 40.3638 | — | — |

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

- Load average at summary time: [1.31, 1.78, 2.09]
- Skipped entries: 0
- Result JSONL paths: `experiments/kernel/bench/results/S1.jsonl`, `experiments/kernel/bench/results/S2.jsonl`, `experiments/kernel/bench/results/S3.jsonl`, `experiments/kernel/bench/results/S4.jsonl`, `experiments/kernel/bench/results/S5.jsonl`, `experiments/kernel/bench/results/S6.jsonl`, `experiments/kernel/bench/results/S7.jsonl`, `experiments/kernel/bench/results/S8.jsonl`
