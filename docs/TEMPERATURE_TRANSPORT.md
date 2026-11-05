# Reusing verified thermodynamic bounds across contact weights

**Scope:** The anchored 2D square-lattice H/P model, not experimental protein thermodynamics. The source Python certificate is checked by the existing `verify_lookahead` verifier; the new transport is exact-rational Python arithmetic and is **not** a Lean-kernel-checked implementation. No physical calibration, speedup guarantee, or novelty claim is implied.

## Mathematical guarantee

Let the number of nonconsecutive hydrophobic contacts in conformation `s` be `K(s)` and let its weight be `q^K(s)` for exact `q >= 1`. For a source parameter `q0` and target `q1`, let `r = q1/q0`. Every state weight transforms by precisely `r^K(s)`; this identity requires no assumption about independence of states.

An event defined by `A = {K >= t}` has contact counts in `[max(0,t), M]` when nonempty, whereas its complement has counts in `[0, min(t-1,M)]` when nonempty. `M` is an unconditional safe cap: the total number of odd-chain-separation, nonconsecutive H/H pairs. When `r >= 1`, multiplication by `r^K` is nondecreasing; when `r < 1`, the endpoint factors reverse. Multiplying checked event/complement mass intervals by the respective extreme factors gives valid target-mass intervals. Total mass is additionally intersected with the whole-ensemble factor interval; probabilities are enclosed by cross-coupled ratios of event and complement mass bounds.

Because all event-state `K` values exceed every complement-state `K`, `P_q(K >= t)` is nondecreasing in `q`. Because `K >= 0`, `Z(q)` is also nondecreasing. Thus, endpoint enclosures bound population and partition function **for every real q** in any closed parameter interval `[q_low, q_high]`, not just selected rational grid points. This is mathematically a contact-weight parameter window; it should not be reported as a Kelvin temperature interval without specifying a separate physical mapping.

## API

```python
from fractions import Fraction
from foldbounds.lattice import build_frontier
from foldbounds.temperature_transport import (
    transport_certificate, transport_certificates, certify_temperature_window,
)

cert = build_frontier("HHHH", Fraction(2), 1, depth=4)
point = transport_certificate(cert, Fraction(3))
window = certify_temperature_window((cert,), Fraction(2), Fraction(3))
assert point.probability_upper == Fraction(6, 13)
assert window.guaranteed_majority == "complement"
```

`transport_certificates` accepts multiple certificates for the **same sequence and event threshold**, but potentially different source `q` values and depths. It checks each input independently, transports each to the target, and intersects the resulting intervals. There is no independence assumption. The optional `lookahead` and ordinary `verify_lookahead` flags are forwarded to the source checker. Invalid target weights, empty sets of anchors, mixed problems, and malformed certificates are rejected rather than silently accepted.

The reference four-H chain has seven zero-contact and two one-contact walks. Hence `Z(q) = 7 + 2q` and `P_q(K>=1) = 2q/(7+2q)`. A source certificate exactly checked at `q=2` transports to an exact enclosure at `q=3` because all event states have `K=1` and all complement states `K=0`. The resulting `P=6/13` is strictly less than `1/2` for every `q in [2,3]`. An enclosure spanning `[2,4]` correctly abstains from a minority/majority classification.

## Trust boundary and test strategy

- At each call, the source certificate is rechecked with the repository's existing Python verifier. Externally supplied claimed mass numbers and arbitrary alleged contact caps are not accepted.
- The new arithmetic uses only `fractions.Fraction`; no floating-point integration, independence assumption, or grid sampling is needed for interval guarantees.
- `tests/test_temperature_transport.py` has a separate self-avoiding-walk enumerator and contact scorer (not the production helpers) and checks *all* H/P sequences of length 2 through 6, across forward and reverse parameter changes, five event thresholds, vacuous events, and multi-anchor intersections. It also tests malformed-source rejection and exact same-parameter behavior.
- This is a mathematical lifting of a **model-relative Python certificate**, not a formally proved code refinement into the Lean 4 kernel. The existing main Python CI automatically discovers the new tests in the isolated PR. No preexisting algorithm, Lean source, workflow, or active branch is modified.

## Limitations and possible future improvement

The all-pairs cap `M` can be weak for long chains, especially at large parameter jumps. A separately certified geometry-specific contact cap or a verified subrange contact histogram could sharpen transport without changing the correctness argument. This implementation does not assert physical predictive accuracy, estimate non-lattice protein free energy, or provide a new experimental benchmark.
