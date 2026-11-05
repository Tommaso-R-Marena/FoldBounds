# Proof-carrying *shape-ranked* HP thermodynamic upper envelopes

**Experimental, Python/model-relative; not a Lean kernel-checked refinement, new physics, or publication-priority claim.** This is an opt-in complement to the existing partially directed spectral lower witnesses and horizon upper certificates.

## The idea: impose a thermodynamic entropy–geometry tradeoff

The old upper mass gave **every possible completion** the same maximum H-H contact count. That ignores the elementary fact that elongated walks cannot achieve as many contacts as compact walks. This method counts relaxed future trajectories by bounding-box shape and uses the shape-conditioned geometric contact cap. Invalid trajectories may revisit sufficiently old sites, so counting them only enlarges the ensemble and cannot exclude true self-avoiding walks.

A finite-memory transfer recurrence tracks endpoint, last 1–6 directed steps, and exact box extrema. It forbids returning to prefix sites and to recent suffix sites. Each true self-avoiding extension is therefore represented **exactly once** among the relaxed walks. After n vertices, reject relaxed trajectories with box area `width*height < n`: no true n-site SAW can fit in fewer than n grid cells. The isoperimetric identity implies every actual score is at most `n+1-width-height`, intersected with the independent H-sequence/parity bound.

### Sharp worst-case mass for the count-plus-shape relaxation

Suppose a score-cap class `k` contains `M_k` relaxed candidates, each real weighted contribution at rational q>=1 bounded by q^k. Separately, the rooted short-block SAW count certificate gives a bound `C` on the total number of actual conformations extending this prefix. The actual conformations form a subset of the relaxed candidates of cardinality at most C. Therefore the maximum weighted sum permitted by **both** relaxations is obtained by selecting the `min(C,sum M_k)` candidates with largest weights: fill classes in descending k and stop after C. The integer selection and exact rational sum are computed explicitly. This is an exact sharp upper bound *for the stated relaxed selection problem*, not for the true HP model.

No independence assumption between geometric shape and contact energy is used; the two are coupled for each relaxed trajectory. This also gives a conservative event upper bound by ignoring classes whose cap is below the event threshold, and a complement upper bound (zero when the fixed prefix already exceeds the threshold).

## Measured ablation (anchored 2D all-H, q=2)

Comparing to the identical short-block count bound multiplied by the identical prefix box cap, with no independent lower mass credited:

| Residues | Upper mass improvement | Transfer memory |
|---|---:|---:|
| 12 | 6.67x | 4 sites |
| 16 | 10.24x | 4 sites |
| 20 | 13.96x | 4 sites |
| 24 | 15.23x | 4 sites |

At n=20 the factor rises to 30.14x for q=3. At q=1 the improvement is only 1.20x, showing that this is primarily an **energy-conditioned** improvement, not merely another count shortcut.

The present Python transfer grows rapidly with n; n=24 took around 15 seconds including two full replays in the local environment, and n=28–30 with four-site memory was not completed under the experiment time budget. Large-number improvements are comparative bounds, **not speedups, absolute accuracy gains, or real-protein predictions**. Many intervals remain too wide for useful long-chain claims.

## Assurance

- Integer transition dynamic programming over a safe superset; no floating point or probabilistic sampling.
- Trusted checks repeat the transfer histogram, recheck fixed prefixes and suffix memory, the area filter, exact integer rooted short-block caps, and rational descending-class optimization.
- Independent exhaustive true-SAW oracle tests across all-H, all-P, mixed sequences through 10 residues and random prefixes, including threshold-event tests.
- Malformed histograms and false total count caps are rejected. The Python replay shares core transition code with the producer; **full implementation-independent assurance and Lean refinement remain open**.

## Reproduce

```bash
python -m pytest tests/test_shape_rank_upper.py tests/test_shape_rank_thermal_integration.py -q
python experiments/benchmark_shape_rank_upper.py
```

Opt in as `thermal_response_bounds(cert, shape_memory=4)`; combine with `use_partially_directed_spectra=True`, `use_geometric_contact_cap=True`, `use_saw_count_cap=True` or `horizon_steps=6` for independently sound intersections. 
