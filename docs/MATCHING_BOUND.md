# A parity-bipartite maximum-flow relaxation for HP lattice contacts

**Research status:** This is an exact, combinatorial, *model-relative*
upper bound and a tested Python implementation. The mathematical argument
has not been formalized in Lean, and publication novelty is unestablished.

## The key observation

In a square lattice each step changes checkerboard color; therefore H-H
contacts can only join amino acid positions with opposite sequence-index
parity. The graph of *potential future contacts* is bipartite.

For a self-avoiding prefix of m vertices:

1. Mark contacts between already placed H residues as fixed.
2. For each H, cap future contacts by its remaining lattice coordination:
   a chain-interior residue has at most 2 non-chain neighbors; an endpoint
   at most 3. For placed H, additionally bound by currently vacant
   neighboring lattice sites and subtract existing H-H contacts.
3. Enumerate potential H-H contact pairs at nonconsecutive odd sequence
   separation that are not already fixed.
4. If the earlier H is placed, remove its edge when no route from the
   current walk tip through still-unoccupied lattice sites can reach a
   vacant neighbor of that H within the available future chain steps.
   We calculate *shortest* distance using ordinary BFS through static
   obstacles; any true self-avoiding extension must obey this necessary
   condition. BFS is deliberately a relaxation (it ignores obstacles that
   would be created by future placements).
5. Give potential edges capacity 1 and H vertices their remaining degree
   capacities. The maximum possible number of future H-H contacts is at
   most the integral maximum flow from even H residues to odd H residues.

**Why sound?** Every complete conformation extending a prefix induces a
set of future H-H contacts. Each such pair appears among candidate edges
(parity and static reachability are necessary). Each uses one available
non-chain neighbor at each of its H endpoints. It is therefore a feasible
integral b-matching, so its cardinality cannot exceed the maximum flow.
Adding fixed contacts yields a valid upper bound on its total score.

The original bound took the minimum of total available candidate edges
and half the total degree capacity. Bipartite b-matching includes both
of those constraints *plus* the distribution of capacity across contact
vertices. It is therefore never weaker than the former score cap,
provided all candidate feasibility tests are applied consistently.

## Integration with probability bounds

The existing certificate checker in lookahead.py accepts the optional
boolean use_matching=True. This substitutes the graph-theoretic score cap
inside the previously documented weighted-mass and lookahead interval
bounds. The certificate checker independently reconstructs the candidate
graph rather than trusting numbers submitted in a certificate.

For reproducibility:

    python -m pip install -e '.[test]'
    python -m pytest
    python -m foldbounds HHHHHHHH --q 2 --depth 2 --enhanced-lookahead 2 --matching

The exhaustive test enumerates *every H/P sequence of length <=8*
and every prefix that appears in a complete anchored 2D self-avoiding
walk, comparing the upper score to the actual maximum completion score
using a separately written enumerator and contact scorer. Other tests
compare the resulting thermodynamic interval enclosures with the exact
partition function and with the older verifier.

## Next mathematical tasks

- Independently formalize parity-bipartite contact graphs in Lean 4.
- Prove placed-residue vacancy and BFS-distance conditions are necessary.
- Prove actual conformational contacts give an integral feasible
  b-matching respecting every vertex capacity.
- Connect an integral flow certificate to a checked residual-network
  optimality argument, or use a directly verified max-flow implementation.
- Establish refinement to the executable Python certificate checker
  without silently assuming that its graph construction is correct.

All existing statistical and mathematical claims remain restricted to the
specified 2D HP lattice toy model; there is no physical-model validation
and no accepted scientific-novelty claim at this stage.
