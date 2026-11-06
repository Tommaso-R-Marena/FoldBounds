# Sharp joint event/complement/partition rational probability intervals

**Status:** Algebraic theorem and tested exact-rational Python computation,
**not Lean-verified**. A compatible optional post-processing step for model-
relative FoldBounds mass intervals. No physical or literature novelty claim.

## Why the usual probability endpoints can be unnecessarily loose

For exact unknown masses `E = Z(event) >=0`, `C = Z(complement) >=0`,
`Z=E+C>0`, suppose independently sound bounds say

    Le <= E <= Ue,   Lc <= C <= Uc,   Lz <= E+C <= Uz,

with strictly positive `Lz`. The standard coupled probability endpoints
`Le/(Le+Uc)` and `Ue/(Ue+Lc)` use only the four box bounds. They may
therefore include impossible pairs violating the certified total upper and
lower partition mass bounds.

The feasible region is the intersection of a closed rectangle with a
diagonal slab. The exact minimum and maximum of `E/(E+C)` over this
finite-dimensional region can be calculated *analytically*, with no floating
point arithmetic or linear programming dependency:

    emin = max(Le, Lz-Uc)
    cmin = max(Lc, Lz-Ue)
    Pmin = emin/(emin + min(Uc, Uz-emin))
    Pmax = min(Ue, Uz-cmin)/(min(Ue, Uz-cmin) + cmin)

The checker additionally sharpens component lower/upper endpoints through
`Emax=min(Ue,Uz-Lc)` and `Cmax=min(Uc,Uz-Le)` and verifies joint
feasibility before returning.

**Proof of sharpness:** For any feasible E,C, `E>=emin`, `C<=min(Uc,Uz-E)`.
The function `E/(E+min(Uc,Uz-E))` is monotone nondecreasing in E on the
admissible range, so the ratio is minimized at E=emin and its maximal
feasible C. The extremal pair is feasible because `emin` enforces the lower
Z slab, while the original bounded rectangle and slab are assumed jointly
feasible. The maximum is the symmetric argument holding C=cmin and
choosing the largest feasible E. `Lz>0` eliminates zero denominator issues.
Thus both exact interval endpoints are *attained in the feasible set*.
This sharpness is only relative to the six provided bounds; the actual HP
ensemble generally obeys additional constraints not encoded here.

## Independent proof-oriented tests

`tests/test_coupled_polytope.py` constructs every vertex by intersection of
two boundary lines among E=Le, E=Ue, C=Lc, C=Uc, E+C=Lz, E+C=Uz;
it rejects infeasible systems and checks sharp rational endpoints against
brute vertex enumeration across thousands of exact rational box/slab examples
and 3000 deterministic random fractional realizations. It additionally
checks component mass validity, rejection of malformed inputs and a strict
improvement regression case.

In the synthetic rational example E in [6,100], C in [8,100], Z in [16,25],
standard `P_lower=6/106≈0.0566`, but the exact joint optimum is `6/25=0.24`.
The sharp upper endpoint is `17/25=0.68` rather than `100/108≈0.926`.
**These are synthetic constraint data, not measured lattice probabilities.**

## FoldBounds integration and trust boundary

The optional `verify_lookahead(...,use_coupled_polytope=True)` computes
previously sound event/complement/total interval values, then applies this
exact algebraic post-processor. No search strategy, certificate generation,
underlying physical model, or transfer/geometry assumptions are changed.
Adaptive replay records the choice and rejects tampering when it changes
the output; the CLI uses `--coupled-polytope`.

For the canonical anchored 2D HP model and rational q>=1, a straight rod
ensures Z>=1, and the existing enhanced checker obtains strictly positive
partition lower mass from genuine ray witnesses. The postprocessor itself
requires a strictly positive certified lower Z and fails closed on invalid
inputs. Neither its Python program equivalence nor the underlying
Python-to-Lean refinement is formalized yet.

**Research limitations:** This is exact rational interval tightening under
known independent constraints, not a claim to a new theorem in fractional
programming. It may achieve no improvement if the total-mass slab does not
exclude any probability-extremizing corner of the event/complement box.
Report those negative cases in matched model benchmarks.