# Soundness boundary and verification argument (Python reference v0.1)

## Physical and combinatorial specification

The sequence is a string over H/P of length n >= 2. A state is a self-avoiding
walk of n vertices on the infinite 2D integer square lattice, with first
vertex (0,0) and second vertex (1,0). Rotations are normalized by the first
edge; reflections are still counted as distinct states. **This is a toy
model, not a predictor of experimental protein structure.**

The contact score k(c) is the number of nonconsecutive H-H residue pairs at
Manhattan distance one. The energy is -epsilon * k(c); the rational parameter
q >= 1 denotes exp(beta * epsilon). All weights q^k, partition sums, and
returned bounds are exact rational numbers. We do not claim to calculate
exp, beta, epsilon, or physical temperature in exact rational form.

For an event A_t = { c | k(c) >= t }, the probability is

    P(A_t) = sum(c in A_t) q^k(c) / sum(c) q^k(c).

## Certificate: a frontier of prefix subtrees

A certificate is a set of prefix nodes, each either a complete state (exact)
or a conservatively bounded subtree (bound). The verifier checks:

1. All prefixes start at the fixed root and are valid self-avoiding walks.
2. All claimed exact paths are terminal.
3. No duplicate or ancestor-shadowed nodes exist.
4. Every valid continuation is covered by a frontier node.
   A dead end has zero completions. The verifier does not trust supplied
   subtree counts, spectra, contact scores, or claimed probabilities.

For a bounded subtree with m fixed vertices, let r = n - m and k0 be the
contacts among placed vertices. Any continuation has a contact score in
[k0, k0 + U], where U counts remaining possible nonconsecutive H-H pairs.
The square lattice is bipartite: pairs with even chain separation cannot
be adjacent. This parity condition reduces U.

A self-avoiding walk has at most three candidate next steps after its first
edge (the predecessor site is occupied). Thus it has at most 3^r completions.
For q >= 1, the subtree's total weight is at most

    3^r * q^(k0 + U).

Its lower weight is conservatively zero. No floating point is needed.

Each exact leaf adds its computed rational weight to the event or its
complement. A bounded subtree contributes certified upper bounds to event
and complement according to the achievable contact range. Write
[EL,EU] for event mass bounds and [NL,NU] for complement mass bounds.
The actual partition function Z is always positive because a straight-line
full-length self-avoiding walk exists. Therefore

    EL / (EL + NU) <= P(A_t) <= EU / (EU + NL),

with safe limiting cases when the relevant known mass is zero. If all
terminal paths are exact, upper and lower bounds coincide.

## What is NOT formally verified

The Python verifier independently checks submitted frontier coverage and
recomputes the relevant weights, but is **not** a Lean-kernel-verified
executable. The argument above is a mathematical proof sketch and pytest
provides software testing, not a formal correctness proof.

Aristotle's Lean 4 work should independently establish the finite-state
theorems, then later prove a refinement theorem relating an executable
certificate checker (including parsing) to those abstract definitions.
Until that bridge exists, do not describe Python outputs as Lean verified.

## Novelty and current limitations

Exact self-avoiding-walk enumeration and elementary bounds on statistical
sums are established ideas. This code currently claims **no new scientific
theorem, state-of-the-art result, experimental validity, or publishable
algorithm**. Novelty requires thorough literature comparison, genuinely
stronger bounds or algorithms, and fair baselines.

The reference enumerator is exponential. Our subtree bounds are deliberately
loose, especially for shallow frontiers. The goal of future research is to
replace these bounds with provably stronger ones and demonstrate practical
speedups without weakening soundness.
