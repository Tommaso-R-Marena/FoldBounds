# Partially directed HP spectra: compressed conformational witnesses

**Status:** exact-integer Python research prototype for the anchored 2D square-lattice HP model; not Lean-kernel verified; not a protein-structure or experimental predictor. Partially directed self-avoiding walk transfer matrices and the discrete perimeter bound are established mathematical ideas. Novelty and publication priority are not claimed.

## Disjoint languages

A canonical path begins at (0,0),(1,0). Count the union of three languages:

1. **E-directed:** moves E,N,S, with only one sign of vertical motion inside any fixed-x column.
2. **N-directed with W:** moves N,E,W, with only one horizontal sign inside any fixed-y row, and **at least one W**.
3. **S-directed with W:** similarly S,E,W, with at least one W.

Every path in each language is self avoiding: its longitudinal coordinate never decreases, and its lateral displacement within a layer has a constant sign. The first family has no W. The second and third cannot overlap when W is present, because their intersection would have no vertical steps but the rooted path begins E and cannot reverse horizontally within its first row. Therefore there is no overlap or inclusion-exclusion correction between these three languages.

Any valid prefix can be imposed as a sequence of *forced first moves*. The returned restricted paths extend that prefix. Different complete frontier subtrees are disjoint by their first divergent lattice step.

## Exact score transfer

For each position index j, the state (d,b,c,seen_W,k) contains:

- d: signed lateral displacement inside the current column/row;
- b: previous layer's first lateral coordinate relative to the current start, or None if no previous layer;
- c: sequence index at the current layer's start;
- seen_W: only used for membership in the N/S family;
- k: exact number of H-H nonconsecutive contacts already formed.

If a lateral step reaches displacement d' inside the previous layer's interval between b and 0, its unique neighboring residue has index

    i = c - 1 - abs(b) + abs(d' - b).

When both residues i,j are H and j-i >= 2, increment k once. A longitudinal step enters a new layer, adjacent only to its backbone predecessor. Nonadjacent layers cannot contact because their longitudinal coordinates differ by more than one. All state multiplicities are exact integers.

Summing terminal counts gives a genuine subset contact polynomial P(x) = sum_k count[k] x^k, which supplies **lower** partition and moment contributions for every rational q >= 1 without enumerating each restricted path. For E-only path counts, the independent known generating function is A(z)=(1+z)/(1-2z-z^2), yielding the recurrence a(0)=1, a(1)=3, a(r)=2a(r-1)+a(r-2).

## Existing frontier / thermal-response interface

The new option

    thermal_response_bounds(cert, use_partially_directed_spectra=True,
                            use_geometric_contact_cap=True,
                            use_saw_count_cap=True)

uses the existing verifier to check a complete frontier before crediting any subset. For every bounded prefix it credits the restricted-family spectrum; individually validated complete-path witnesses already in that same family are subtracted from the spectrum **exactly once**. Other witnesses are credited directly, and terminal exact nodes are never counted twice.

The previous cone-monotone subset and the partially directed subset overlap in general. Simultaneous use of the two spectral flags is rejected instead of unsafely adding their masses.

The optional geometric contact cap invokes the existing prefix-aware integer isoperimetric bound: all H-H contacts are among the n+1-w-h nonbackbone grid contacts, where w,h are final bounding-box dimensions satisfying w*h >= n and including the prefix. Replace the unknown subtree's largest possible score with the minimum of this cap and the parity-based cap, then cap its mass by count_upper*q**score_max. If a prefix already exceeds the geometric maximum for a completed chain, the corresponding completion set is empty and its upper mass is zero. This is optional and composable with the verified short-block SAW counting cap.

The existing PR #51 primal-dual verifier continues to independently check the final relaxed variance upper bound. Python checks do not amount to Lean kernel verification.

## Local results and limitations

At n=30, all-H and q=1, the old monotone subset contains 805,306,365 distinct paths; the new three-language union contains 151,602,131,587, a 188.25-fold increase in certified lower partition mass. The combined lower mean contact bound is around 0.08877, upper mean bound around 19.45317, variance lower bound around 0.00332182, and lower probability of at least one contact around 0.0292. These intervals are still broad.

At n=40, the path count gain rises to about 1,240-fold, but the absolute lower variance remains around 0.00022. At q>1, unknown high-contact mass can still dominate conservative upper bounds, so large *relative* improvements do not imply accurate experimental thermodynamics or a breakthrough in folding physics.

## Verification and reproduction

    python -m pip install -e '.[test]'
    python -m pytest tests/test_partial_directed.py tests/test_directed_integration.py -q
    python experiments/benchmark_partial_directed.py

The independent test enumerator checks all rooted SAWs through n=8, multiple HP strings and prefixes, and full contact spectra. A separate full-ensemble oracle validates exact rational means, variances, event probabilities and susceptibilities across short exhaustive sequences and randomized q/depth/witness/threshold cases. Invalid overlapping modes fail closed.

**Trust boundary:** a proof of Python DP equivalence to Lean's finite lattice model, large-scale experimental calibration, and a novelty literature audit are not supplied by this draft.

### Background prior art

- Redner & Majid (1983), Critical properties of directed self-avoiding walks, J. Phys. A 16 L307, doi:10.1088/0305-4470/16/9/005.
- Brak et al. (2008), Journal of Combinatorial Theory, Series A 115, 623–650, explicit partially directed generating function.
