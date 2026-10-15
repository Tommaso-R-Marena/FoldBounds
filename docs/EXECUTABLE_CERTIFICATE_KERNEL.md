# Lean-native executable rational certificate kernel

## Verified scope

FoldBounds/Executable.lean introduces a small executable and a Lean theorem
about its actual acceptance decision, not merely an abstract mathematical
inequality. An untrusted proposer supplies a subset of witnessed states in
each of two complementary event classes and a per-state rational upper cap.

The checker rejects certificates unless it can **decide** all obligations:

- each witness belongs to the claimed event class;
- rational weights are positive on the finite universe;
- every not-yet-witnessed state's exact model weight is <= its claimed cap.

It computes exact rational lower mass by summing distinct Finset witnesses.
Its upper mass adds the weight cap on the finite unwitnessed tail. The Lean
theorem check_sound states that if check returns some rational endpoints,
both true finite-model event masses belong to the corresponding intervals.

A separate Lean file tests accepted and rejected cases using by decide.
This avoids a user-provided proof hypothesis standing in for a check.

## Important limitations

**This checker visits every element of its explicitly finite universe
to verify every upper-weight obligation.** The exercise establishes a
small trusted, executable specification and a proof of its decision
semantics; it does not establish a non-enumerative algorithm.

It is NOT a formal verification of the Python HP-lattice walk validator,
prefix coverage, maximum flow, isoperimetry, transfer DP, or parser.
It does NOT establish that an arbitrary claimed list of states is the
complete physical conformation universe. It is a rational statistical
weight model; relating q^k to exp(-βE) requires additional assumptions.

There is no novel theorem or protein-physics accuracy claim here.

## Next steps

1. Introduce an independent Lean definition of canonical finite lattice
   conformations and a checked map from parsed paths to that type.
2. Prove complete self-avoiding-walk generation / frontier coverage.
3. Prove contact scores and rational q^k model correctness.
4. Replace fully enumerated per-state upper checks with **verified**
   subtree combinatorial caps, without widening the trusted assumptions.
5. Connect the executable bounds to Aristotle's coupled probability bounds.
6. Compare matched published algorithms at n >= 14 under equal budgets.
