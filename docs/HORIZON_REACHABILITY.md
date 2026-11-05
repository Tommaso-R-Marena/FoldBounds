# Independently checked prefix-horizon upper certificates

**Status:** exact rational/integer Python research prototype, stacked on the draft FoldBounds directed-spectral PR. This is not a Lean-checked implementation, a new SAW asymptotic theorem, a validated protein structure predictor, or a claim of originality in lattice-polymer branch-and-bound.

## New complement to directed subset witnesses

The existing directed-spectral module supplies valid **lower** Boltzmann mass by compressing exponentially large sublanguages. This module works on the **uncounted complement**: it enumerates a bounded number of genuinely reachable partial walks, rigorously covers every legal branch, and assigns each leaf a conservative Boltzmann **upper** mass.

For a leaf with m fixed vertices, n total residues, contact score floor L and certified cap U, a complete extension's weight at rational q>=1 cannot exceed q^U. The number of extensions cannot exceed the separately checked short-block self-avoiding-walk product A(n-m). Therefore its mass is bounded by `A(n-m)*q**U`. The cap U combines parity-compatible remaining H-H pairs and the prefix-dependent 2D box perimeter inequality. Completed leaves contribute exactly, while proved-impossible leaves contribute zero.

### Finite-component parity obstruction

Every legal continuation starts with an unoccupied neighbor of the tip and stays within its connected vacancy component (it cannot cross a previously occupied vertex). If the component is finite, an r-site suffix requires `ceil(r/2)` sites on the checkerboard parity opposite the tip and `floor(r/2)` on the same parity. Failure of either necessary capacity proves that component has no completion. If *every* neighboring component fails, a prefix is impossible even if an immediate free move exists. The implementation stops flood-filling large components conservatively (accept unknown possibilities) and never assumes that finite search found the full infinite exterior.

### Fail-closed verifier

The producer is not trusted. The independent verifier revalidates anchored self-avoiding geometry, integer coordinate types, leaf completeness and prefix disjointness, every omitted branch, all score ceilings, every exact-leaf mass, every unsupported dead-claim, and each short-block count-cap certificate. Modifying an upper cap downward, omitting a branch, duplicating a leaf, or inflating an exact mass is rejected. A checked envelope is valid at any rational q>=1, not just the q used to select the depth.

Splitting a block can perversely worsen a loose combinatorial bound, so `select_best_horizon_certificate` evaluates and rechecks each candidate depth, then chooses the *smallest* certified partition upper bound. Thus the opt-in path never worsens the baseline upper mass.

## Controlled observed cases

All-H, anchored 2D model, q=2, root prefix, horizon at most 8:

- n=10: coarse/exact bound 4.88x; depth-8 is exactly enumerative and recovers Z.
- n=12: upper mass decreases from 2,342,592 to 704,970 (3.32x improvement), while exact Z is 142,757. The refined bound is still approximately 4.94x the exact Z.
- n=20: upper mass decreases by approximately 1.43x.
- n=30: upper mass decreases by approximately 1.10x.

The improvement is nonuniform and still insufficient to make long-chain strong-coupling thermodynamics tight. The principal new assurance is **two-sided**: an existing compressed valid subset establishes a lower mass and the horizon tree provides an independently checked upper mass. Neither result predicts real proteins.

## Reproduce

```bash
python -m pytest tests/test_horizon_certificate.py tests/test_horizon_thermal_integration.py -q
python experiments/benchmark_horizon_envelopes.py
```

Optional API: `thermal_response_bounds(cert, use_saw_count_cap=True, use_geometric_contact_cap=True, use_partially_directed_spectra=True, horizon_steps=6)`. Increasing `horizon_steps` can cost exponentially more time and memory; it is off by default.

## Open proof obligations

1. Formalize the parity-size obstruction for connected components of the unoccupied lattice after a finite SAW prefix.
2. Prove the leaf-frontier coverage theorem and exact/zero/upper leaf mass recurrence.
3. Prove each rooted block segmentation provides an upper count on actual completions.
4. Prove the contact ceiling and exact rational weighted enclosure.
5. Establish a trusted parser/refinement bridge from Python certificate serialization to Lean-checkable propositions.

Existing Monte Carlo, structural benchmarks, and actual protein folding predictions are outside the result's trust boundary.

## Event-conditioned upper accounting (additive extension)

The integration also retains the horizon checker's **separate event and
complement** mass upper bounds. When independently validated whole-walk
witnesses or a directed-spectrum subset are credited as known mass, their
exact event/complement contributions are removed from the corresponding
upper caps exactly once. The surviving residual class caps are intersected
with the residual total-mass cap before tightening the probability interval.
This prevents the common error of treating a witness as both known mass and
unknown event/complement mass. Existing variance certificates remain sound but
are not yet optimized with these additional class constraints.
