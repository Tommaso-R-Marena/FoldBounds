# Endpoint-conditioned finite-memory transfer upper bounds

**Status:** Optional, exact-rational Python research method. The separate
Aristotle Lean theorem package does not establish correctness of this Python
algorithm, its parser, or the experimental measurements.

## Why endpoint geometry can tighten an upper mass

For a canonical square-lattice self-avoiding walk of n vertices let S be
its occupied sites, and let its bounding box have integer width w and
height h. Its number of *all* noncovalent lattice contacts is bounded by

    n + 1 - w - h,

because every row and column contributes at least two outward perimeter
edges and the backbone accounts for exactly n-1 internal lattice edges.

The fixed prefix imposes a minimum rectangle. An endpoint T = (x, y),
as represented exactly in a finite-memory transfer-DP terminal state,
imposes another compulsory occupied site. For any genuine complete
conformation reaching T, its final bounding box must contain both the
prefix and T. If w0,h0 are the width and height of the smallest box
containing those sites, then the actual full box has w>=w0, h>=h0,
and w*h>=n.

Define

    B(n,w0,h0) =
      max(0, n+1 - min{w+h | w>=w0, h>=h0, w*h>=n}).

Then **every actual self-avoiding completion reaching T** has at most
B(n,w0,h0) total noncovalent contacts, hence at most that many H-H
contacts. The existing prefix-only isoperimetric bound is simply the
case without T.

For an actual completion with k_fixed preexisting H-H contacts,
k_cross contacts connecting prefix to suffix, and k_future between
suffix residues, a separate parity-bipartite matching certificate gives
k_future <= flow_cap. Thus

    k_total <= k_fixed + k_cross +
       min(flow_cap, max(0, B(n,w0,h0)-k_fixed-k_cross)).

This upper score determines a conservative rational Boltzmann
upper weight q^upper_score when q>=1. If the upper score is below
the event threshold, that path contributes zero to the event upper
bound. Complement exclusion still uses the separate guaranteed
lower contact count from visible suffix contacts; this refinement
does not change that lower-bound criterion.

## No invalid relaxation assumptions

The finite-memory DP may count relaxed walks that self-intersect
outside its K-site window, so their geometry may be nonphysical.
These invalid paths correspond to **no actual conformations**. Their
assigned weight is nonnegative and never needed for the validity
of an upper bound. Every genuine conformation remains present
within its correct endpoint class and has an admissible upper
weight based on that endpoint. Therefore summing nonnegative
endpoint-class upper weights covers all genuine model mass.

The checker does not claim that the terminal position determines
the complete path, or that a relaxed path's exact contact score is
known from its truncated history.

## Opt-in and controlled validation

    python -m pytest tests/test_terminal_box_transfer.py
    python experiments/compare_terminal_box.py
    python -m foldbounds HPHHPHHPHHHH --q 3/2 --event-at-least 2 \
      --enhanced-lookahead 0 --matching --transfer \
      --transfer-memory 4 --terminal-box --monotone --turn-family

The original transfer behavior remains the default. The optional
`--terminal-box` mode also flows through adaptive search and
is recorded and replayed by independent audit checks. It cannot
be used without enabling transfer.

The independent tests (1) enumerate complete short trajectories of
the relaxed memory-K model and compare their terminal/cross/local
histograms with DP aggregation, (2) compare exact small conformational
masses to both old and new upper bounds, and (3) verify monotonicity
of the upper mass and adaptive replay. These are Python test oracles,
not Lean formal verification.

**Limits:** Endpoint grouping adds histogram coefficients and costs CPU
time/memory. A tighter mathematical cap does not imply a practical
speedup. This is a known kind of elementary perimeter/isoperimetric
reasoning, not claimed as a novel scientific method or real-protein
structure predictor. Published model-matched solver comparisons
are still necessary before any publication/performance assertion.
