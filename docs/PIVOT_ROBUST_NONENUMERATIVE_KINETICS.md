# Nonenumerative kinetic/equilibrium separation under local + pivot moves

**Research prototype, with a new independent Python checker.** This is an
explicitly specified 2D square-lattice HP model, not a real protein, a
biological folding-time estimate, a Lean-kernel theorem, or a priority claim.

## Main contribution

The prior finite-state kinetic certificates were exact but needed exhaustive
walk enumeration and featured a local move set with disconnected components.
Here an independently auditable small-*plateau* certificate works without
examining the vast majority of folds and survives the addition of **all seven
nonidentity D4 suffix pivots**, which can radically rearrange conformations.
The enlarged chain was independently found connected for n=5..10 by explicit
full-state graph exploration. Connectivity for arbitrary n is not proved by
our checker; pivot-walk ergodicity is established in the literature.

For even height h >= 10, consider a canonical 4 x h Hamiltonian snake of
n=4h **all-H** lattice residues. The producer proposes its connected
zero-barrier plateau. An independent verifier replays *every* legal distinct
proposal and its contact score, confirms complete closure at contact score
k0=3h-3, and validates a separate 8 x (h/2) Hamiltonian path with score
k1=7h/2-7. The constructive **energy gap** is k1-k0=h/2-4.

This rectangular family is directly checked for several sizes up to n=1200;
we do not claim a Lean proof that the plateau always contains exactly eight
states for *every* h. The one-line geometric score formula is proved by
counting nearest-neighbor lattice edges minus n-1 covalent bonds.

## Dynamics and probability statements

For any valid rooted path of n beads, one update chooses uniformly among
M=11(n-2) slots: for each internal/terminal index i=2..n-1, four immediate
neighbor positions of bead i-1; or, for each pivot index j=1..n-2, one of
the seven nonidentity D4 isometries of the suffix j+1..n-1. Invalid paths,
unchanged paths and rejected Metropolis proposals are self-loops. Every legal
distinct proposal has a reverse slot with precisely equal multiplicity,
because D4 inverse operations occur with the same frequency and local inverse
relocations use the same preceding bead. For rational q >= 1, accept a legal
proposal x->y with min(1, q^(k(y)-k(x))). Hence full-space Gibbs weights
pi(x) proportional to q^k(x) satisfy detailed balance, without assuming mixing.

A proposed plateau S is accepted only if the checker independently verifies:
- all states have the same contact score k0, including a specified seed x;
- every one of the M proposal slots from every state is enumerated, preserving
  multiplicity, and every accepted destination with score >=k0 is already S;
- S is the entire zero-barrier component of x, with a separately validated
  better-score witness outside S;
- every outward proposal has an integer contact penalty d>=1.

Let a_d(s) count outward *proposal slots* losing exactly d contacts from s.
Then the conditional probability of escaping S in the next step, while at s,
is exactly H_q(s)=sum_d a_d(s) q^(-d)/M. For a deterministic x in S,

    P_x[ ever reach {k>k0} by T ] <=
        P_x[exit S by T] <= min(1, T * max_{s in S} H_q(s)).

The proof is a union bound on first exits using a bound on the exit hazard
conditional on *any* surviving state in S. It does **not** assume thermal
stationarity, rapid within-basin mixing, spectral approximations or the
entire conformation enumeration. For **every real q>=1**, not merely tested
rationals, the finite combinatorial witness also implies

    H_q(s) <= B_s/[M q], where B_s=sum_d a_d(s).

The checker computes exact rational finite-q hazards. If instead the initial
state is conditional Gibbs within S (uniform over its equal-score members),
stationarity domination gives a separate upper bound
T * sum_s H_q(s)/|S|. This is not the default fixed-start bound.

The higher-score event A={k>k0} has Gibbs probability at least

    pi(A) >= q^gap / (q^gap + 3^(n-2)).

Reason: the one independently checked witness contributes weight q^(k0+gap)
to A; any complementary path has k<=k0, and there are at most 3^(n-2)
rooted self-avoiding paths. This bound is conservative and independent of
full enumeration. **This is not a certified global maximum or unique native
fold**; A is simply the explicitly specified better-contact event.

## Reproduced results, exact integer/rational evaluation

Choose q=1,000,000,000 and T=1,000,000 proposed updates.

| Beads n | Plateau members | Witness contact improvement | Fixed-start target-hit upper | Equilibrium bad-event upper |
|---:|---:|---:|---:|---:|
| 40 | 8 | 1 | 1.197e-5 | weak (~1) |
| 80 | 8 | 6 | 5.828e-6 | < 1e-16 |
| 200 | 8 | 21 | 2.296e-6 | < 1e-94 |
| 400 | 8 | 46 | 1.143e-6 | < 1e-224 |
| 800 | 8 | 96 | 5.697e-7 | < 1e-483 |
| 1200 | 8 | 146 | 3.795e-7 | < 1e-742 |

The 1200-bead checker took approximately 30 seconds in the original local
Python environment; it enumerated **zero full conformation ensemble states**.
It explicitly inspected all local/pivot proposals from the eight submitted
plateau members. The Boltzmann factor q is dimensionless, not an assumed
physical temperature, and one time step is an attempted artificial move.

## Threat model, limitations, prior work

The producer is entirely untrusted. The independent checker uses separate
geometry, contact, pivot, and local-corner implementations. Its fail-closed
budget refuses partial plateau sets, malformed/missing extra members, invalid
witnesses, inconsistent energies, and resource overruns. No Lean/Python
refinement theorem has been established. A mathematical proof of the generic
Markov first-exit inequality is given above, and pytest validates instances;
these are distinct from kernel verification.

- Kinetic traps in lattice HP protein models and strong move-set dependence
  were documented by Chan & Dill, *J. Chem. Phys.* 100, 9238 (1994),
  DOI 10.1063/1.466677.
- The well-studied pivot algorithm for self-avoiding walks and its ergodicity
  go back decades; see Kennedy, *J. Stat. Phys.* 106, 407 (2002),
  DOI 10.1023/A:1013750203191; Clisby, *J. Stat. Phys.* 140, 349 (2010),
  DOI 10.1007/s10955-010-9994-8.
- Certified kinetic trapping **does not** imply biologically slow folding;
  HP, square lattice, all-H composition, the move rule, and chosen q are
  significant modeling assumptions. At q=1, the claimed huge gap disappears.
- No new theorem originality, actual 3D folding prediction, molecular-dynamics
  validity, or scientific-priority claim is made. The interesting advance is
  proof-carrying scaling plus explicit falsification of simpler move-set traps.

## Reproduce

    PYTHONPATH=src python -m pytest tests/test_pivot_plateau_certificates.py -q
    PYTHONPATH=src python experiments/benchmark_pivot_plateau.py

The test suite checks independent proposal multiplicities and reversibility
on every rooted SAW up to 8 beads, exhaustive pivot-graph connectivity up to
10 beads, all small H/P contact scores, detailed balance for rational q,
large-n proofs, all-q rational examples, and adversarial certificate mutation.
Keep this stacked branch as a **draft**, separate from the earlier kinetic
PR and all other workstreams. Its output is not merged into main.
