# Exact contact-map atlas (independent research utility)

This optional module constructs a **geometry-only** exact histogram of contact graphs for the anchored 2D square-lattice HP model, with the root sites `(0,0)`, `(1,0)` fixed. It does not change the existing lattice verifier, Lean files, CI, or motion research. It is an exponential reference oracle, not a scalable folding algorithm or a biological predictor.

## Why it is useful

A contact map records possible nonconsecutive lattice-neighbor pairs `(i,j)` with `j-i>=3` odd. Each fully self-avoiding walk contributes a bitmask. **The atlas depends only on chain length, not sequence or q**, so one finite geometry enumeration supports unlimited exact H/P mutations and rational contact weights without walking the tree again. For each H/P sequence, intersect its hydrophobic-pair bitmask with each geometry mask, then count set bits to recover its contact score. Grouping walks by equal masks makes subsequent queries much cheaper.

An additional exact symmetry reduction is available: reflection in the x-axis maps each walk with a first vertical step north to one with first vertical step south. It preserves the contact map and pairs all nonstraight walks. The only fixed walk is the all-east straight walk. Accordingly the generator visits north-first representatives and assigns weight two (or one for the straight walk). `audit_atlas` deliberately *does not use this symmetry reduction* and independently recomputes every anchored self-avoiding walk, scoring contacts from completed coordinates. It refuses length above the explicit independent-audit budget rather than issuing an unsupported verification claim.

## API / reproduction

```bash
python -m pytest tests/test_contact_map_atlas.py
python experiments/contact_map_atlas_demo.py
python experiments/benchmark_contact_map_atlas.py
```

```python
from fractions import Fraction
from foldbounds.contact_map_atlas import build_atlas, audit_atlas

a = build_atlas(9, max_nodes=2_000_000)
assert audit_atlas(a)  # full independent re-enumeration
assert sum(a.spectrum('HHHHHHHHH').values()) == a.conformations
z = a.partition('HHPPHHPHH', Fraction(2))
p = a.event_probability('HHPPHHPHH', Fraction(2), threshold=1)
mean, variance = a.contact_moments('HHPPHHPHH', Fraction(2))
```

The exact mean and variance of contact count satisfy `q * Z'(q)/Z(q) = E[K]` and `d E[K]/d(log q) = Var(K) >= 0`. `score_covariance(reference,comparison,q)` computes covariance of both sequence scores **over the same conformations**, weighted by the reference sequence. It is a descriptive finite-model quantity, not protein-level biological covariance. The module checks type/range/structural invariants, rejects float q, and aborts on an enumeration-budget overrun without returning a partial atlas.

## Assurance and limits

- Tests compare **every HP sequence of length 2 through 8** with a separate brute-force walk scorer at three exact rational q values and four event thresholds. Established anchored-walk counts through length 9 are regression checks.
- `audit_atlas` is a second enumeration with different direction ordering and leaf scoring; this is strong **software cross-validation**, but not a Lean kernel proof.
- Constructing an instance directly from arbitrary multiplicities may pass shape checks without being genuine. `audit_atlas` is required before using untrusted data; there is no serialization or integrity claim.
- The search is exponential and can consume substantial time/memory. `max_nodes` provides a fail-closed generation budget; the independent audit has a `max_length` gate. Do not confuse these with thermodynamic bounds for large chains.
- No assertion of algorithmic novelty, experimental folding accuracy, or completeness of real protein ensembles is made.

Potential future connection (separate review required): use atlas-derived exact mutation and covariance tables as fixed test oracles for the experimental crossover and mutation-energy projects without coupling to their production implementations.
