# Infinite-size contact-density certificates (2D HP homopolymer)

## What has actually been proved

A deterministic, finite, **exact-integer/rational** checker certifies a sufficient bound for a specified conformational contact event, not a simulation estimate. For the all-H, square-lattice self-avoiding-walk model with Boltzmann contact weight q, it finds an infinite arithmetic progression of square chain sizes on which **the probability that at least a specified fraction of sites participates in the contact-count threshold event** exceeds a requested confidence level. More precisely, the event is K >= rho * n where K is the *number of nonconsecutive H-H nearest-neighbor contacts*, not the number of residues participating in contacts. Each distinct contact uses two sites, so "25% contacts" means K >= n/4, NOT 25% of residues with at least one contact.

The method makes no independence assumption about microscopic corrections: every conformation may be multiplied by any factor in a fixed interval [l,u], with no correlation restrictions, even varying with chain length and q. The resulting guarantee is model-relative and uniform over all such weight assignments.

This is an elementary Peierls-type entropy/energy domination proof. The submultiplicative self-avoiding-walk bound and the strong-coupling intuition are established prior art; **we do not claim a new connective-constant bound, a novel thermodynamic phase transition, or actual protein-fold prediction**. The exact Python checker is NOT a Lean-kernel-checked proof. Independent scrutiny of the proof and implementation remains necessary.

## Construction and exact inequality

Take a square of side m (m divisible by 4), with n=m^2 all-H residues. An explicit serpentine Hamiltonian path visits each site once, starts at (0,0)->(1,0), and therefore is a legal canonical SAW. The square has exactly 2m(m-1) nearest-neighbor edges and the chain uses n-1 of them; the difference is the exact witness contact count

    Kstar(m) = m^2 - 2m + 1 = (m-1)^2.

A universal count of anchored SAWs of length n is bounded by multiplying maximal free-space extension counts. The checker enumerates *only* the small independent 0..16 step free-space trees by DFS, obtaining

    c_14 = 1,604,149
    c_16 = 11,616,669.

Even conditioned on an arbitrary prior self-avoiding prefix, the number of possible k-step extensions is <= c_k, since previous sites are obstacles. With n divisible by 16, the n-2 remaining steps split into n/16-1 blocks of 16 and one block of 14. Because the two **directly verified** integer counts obey c14 <= c16, the entire enormous path count satisfies

    number_of_full_walks(n) <= (c16)^(n/16).

Importantly, no assumption is made that a trapped 14-step walk can always be extended.

For a fixed rational density rho=a/b in (0,1), assume that m is divisible by s=lcm(4,b). This ensures that t=rho*m^2 is an integer. Let the conformational event A be {K>=t}. Every *unobserved non-event* has K<=t-1. The verified square witness belongs to A when Kstar>=t. With each microstate multiplier a_s in [l,u], the ratio of non-event to event partition mass obeys

    N/E <= (u/l) * (c16)^(m^2/16) * q^(t-1-Kstar)
         = F(m)
         = (u/l) * (c16)^(m^2/16) * q^[-(1-rho)m^2+2m-2].

This uses one event-state lower witness and an unconditional upper count for the entire non-event ensemble. Therefore

    P_m(A) >= 1/(1+F(m)).

Everything on the right is **exact rational arithmetic** for rational l,u,rho and integer q.

## The infinite-tail gate (proof by induction, not size sampling)

Define R(m)=F(m+s)/F(m), for sides m divisible by s. By direct algebra,

    R(m+s)/R(m) = (c16)^(s^2/8) * q^[-2(1-rho)s^2].

The sufficient strict inequality

    (c16)^b < q^[16*(b-a)]

is checked with integers, so the right-hand ratio is strictly less than 1 for every eligible m. Hence R(m) decreases geometrically along the progression.

At a selected starting side m0, the checker verifies R(m0)<=1. Consequently **F(m0+ks) is nonincreasing for every integer k>=0**. A single finite exact-integer calculation at m0 therefore supplies the uniform probabilistic lower bound

    P_{m0+ks}(A) >= 1/(1+F(m0))    for all k>=0.

As k tends to infinity, F(m0+ks) tends to zero with an exponentially negative leading term in m^2 (the O(m) boundary term is subleading). Thus the event's lower probability approaches 1 in the stated family, even under arbitrary allowed reweightings. Neither the program nor its claim relies on enumerating a finite grid of all future chain lengths.

## Demonstrations (l=1/2, u=2)

| Contact weight q | Threshold t/n | First certified side m for P>=99% | Chain length n | Witness Kstar | Rigorous lower probability at m |
|---|---:|---:|---:|---:|---:|
| 3 | 1/20 | 100 | 10,000 | 9,801 | > 0.999999999999999999999 |
| 4 | 1/4 | **124** | **15,376** | **15,129** | **> 0.999980** |
| 5 | 1/3 | 60 | 3,600 | 3,481 | > 0.999984 |

For the middle row, the event is K>=3844, and the result applies to **every** m=124+4k. An independent m=128 calculation (n=16384, Kstar=16129, t=4096) gives a bound above 0.999999999. The earliest side is earliest **under this conservative proof strategy and congruence restriction**, not the true model transition.

Notice the bound is not informative at all parameter choices (e.g. q=3, rho=1/4 fails the leading integer inequality). Such failures are **abstentions**, not evidence of absence of a conformational transition.

## Reproduce

    python -m pip install -e '.[test]'
    python -m pytest tests/test_square_density.py
    python experiments/square_density_demo.py

Implementation: src/foldbounds/square_density.py, isolated from any existing solver or Lean proof. All short-block numbers are recomputed, not trusted from submitted files. Tests independently enumerate short-chain walks, recompute contact counts with O(n^2) direct scoring, verify the exact rational ratio recurrence, reject invalid inputs, and check both favorable and unfavorable regimes.

## Formal trust boundary and significance

- Exact arithmetic is used, but the Python implementation itself is **not formally verified by Lean**. The existing Lean CI may pass while this new assertion remains Python-checked only.
- The model is 2D, uses the highly artificial all-H homopolymer, and omits steric side chains, solvent-mediated forces, hydrogen bonding, kinetics, realistic entropic/enthalpic corrections, and finite experimental calibration.
- q is a **dimensionless contact weight**, not a Kelvin temperature. The parameter could be interpreted as exp(beta*epsilon) in this toy Hamiltonian; any mapping to actual molecular epsilon and temperature requires independent physical justification.
- The result establishes an extensive contact count for favorable strong-coupling regimes, not a uniquely determined structure, native-state assignment, or trustworthy prediction of any real protein.
- The bound intentionally ignores all other event conformations and grossly overcounts non-events; it is generally weak at realistic heterogeneous-sequence parameters.
- A relevant prior-art result is the square-lattice connective-constant upper bound (Couronné, arXiv:2211.16146; published 2025). A more recent Oct 2026 preprint by Tejasvi Singh Tomar reports an improved upper bound with Lean-related verification (DOI 10.5281/zenodo.23540764). Neither supports replacing the finite-size c16 bound by mu^n without a separately justified finite-n prefactor.
- A possible next step is to place the finite-sum inequality and the arithmetic-progression induction in Lean, and to upgrade the unstructured all-H argument with verified coarse-grained conformational families and physically grounded energy models.

**Bottom line:** a rigorously checkable *infinite-size scaling certificate for one very simple polymer model*, not a demonstrated new law of real protein folding.
