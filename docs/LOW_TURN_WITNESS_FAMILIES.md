# Constructive low-turn conformation witnesses

**Status:** exact rational tested Python research method, not a Lean-verified
executable, a novelty claim, or a real-protein model.

## Why

The existing monotone-quadrant lower-mass family forbids suffix/suffix
nonconsecutive contacts. For event thresholds requiring multiple contacts,
it can produce a zero event lower mass even when valid event conformations
exist. In particular, HPHHPHHPHHHH at q=3/2 and event k>=2 has legitimate
2-contact folds but the initial monotone lower event mass is zero.

## Restricted constructive family

For an arbitrary valid anchored prefix of m residues, generate valid
suffix direction words (first suffix move has no prior suffix direction)
with at most T changes between successive directions, 0<=T<=2. Check
each proposed next site against the entire current path, rejecting every
self-intersection. Every retained completed path is a distinct genuine
canonical lattice self-avoiding walk.

For suffix length r>=1, the number of *unfiltered* direction words
with at most two turns is at most

    4 * [1 + 3(r-1) + 9*binom(r-1,2)].

Thus there are O(r^2) candidates for fixed T=2. With O(r)-sized paths
and a direct full-length contact scorer, the reference implementation
is polynomial-time for this **restricted subset**. This is NOT a
polynomial-time partition algorithm for all self-avoiding walks.

For every accepted completed path, compute k exactly and add the
exact rational q^k to either the event sum or the complement sum.
These are simultaneously sound lower masses since the two event classes
are disjoint. This family's folds may overlap monotone, ray or DFS
families. FoldBounds therefore combines lower masses via separate
maximums, NEVER addition across overlapping families.

## Reproducibility and trust

    python -m pytest tests/test_turn_family.py
    python experiments/compare_turn_family.py
    python -m foldbounds HPHHPHHPHHHH --q 3/2 --event-at-least 2 \
        --enhanced-lookahead 0 --matching --transfer \
        --transfer-memory 4 --monotone --turn-family

The test oracle separately enumerates all genuine short conformations
and computes each path's turn count and contact score from full
coordinate lists. For the mixed sequence, a two-turn family at the
canonical root contributes exact event mass 459/8 at q=3/2, threshold 2.

The production checker treats turn-family outputs as Python-computed
data, NOT as formally proved Lean witness enumeration. Formalizing this
construction requires proving restricted enumeration has no duplicates,
only generates self-avoiding paths, counts each q-weight once, and
connects the checker output to Lean's finite HP model.

## Negative cases and limitations

Restricted low-turn families can be exponentially small compared with
the true fold population and do not dominate monotone families. In
many longer all-H sequences, monotone-family event mass is stronger.
A guaranteed strictly positive event mass may still yield an
almost-vacuous probability lower bound when the complement upper
mass is enormous. Runtime grows with suffix length and independent
frontier nodes. These results are not yet compared with published
lattice-polymer algorithms under matched budgets.
