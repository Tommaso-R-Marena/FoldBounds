# Lean-native anchored HP lattice model

This module defines finite canonical self-avoiding square-lattice walks
as subtype-checked absolute-direction words with an anchored first bond.
It computes hydrophobic-hydrophobic contacts on integer coordinates
and exact rational statistical weights q raised to contact count.

The checkHP function specializes the proved executable finite-rational
certificate checker to the Lean lattice-conformation universe. The
checkHP_coupled_sound theorem states that successful checks enclose
the rational-model event probability, with the same endpoints actually
returned by the checker.

## Trust boundary

- The checker enumerates the full finite word universe to validate
  caps. It is not a scalable certified lattice solver.
- No equivalence to the Python path parser, contact scorer, frontier
  generator, max-flow, isoperimetry, or transfer relaxation is proved.
- This is a combinatorial q model, not experimental molecular physics.
- The geometric representation is not yet proved bijective with every
  canonical lattice path; although decoded words are required to be
  self-avoiding, completeness remains an explicit separate obligation.

## Next formal proof obligations

1. Prove all generated successive sites are cardinal neighbors.
2. Prove bijection between self-avoiding words and anchored site paths.
3. Verify external path serialization against the Lean decoder.
4. Replace full-universe cap inspection with sound family certificates.
5. Benchmark published lattice-polymer methods on aligned models.
## Verified lattice adjacency

The new direction-advance theorem and chain theorem prove that every
conformation decoded from a direction word contains only cardinally
adjacent consecutive vertices. Together with the subtype-enforced
no-repeat property, this establishes both geometric safety conditions
for every Lean-model conformation without assuming Python validity.
This remains distinct from a proof of path-decoder completeness or
cross-language correspondence.
