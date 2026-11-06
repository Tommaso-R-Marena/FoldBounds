# A geometry-only heterogeneous HP concentration certificate

This is a **sound but conservative** exact-rational extension of the
all-H ground-state count in #68. It does not claim that the counted
Hamiltonian squares are low-energy states for any general H/P sequence.

## The contact-loss argument

Fix a genuine self-avoiding walk on the full m×m lattice square. Its
total number of nonconsecutive lattice contacts, before assigning H/P
labels, is exactly K_all=(m-1)^2.

In a square lattice each internal chain vertex spends two of its four
grid incidences on mandatory backbone edges. Thus any *internal* residue
can participate in at most two noncovalent contacts. An endpoint spends
one edge on backbone, and can participate in at most three.

Changing any residue to P deletes from the H-H contact score no more than
the noncovalent degree of that residue. If both ends of one contact are P,
the union bound may count that missing contact twice, which is safely
conservative. Therefore **every** one of the C independently counted
Hamiltonian square walks, for a H/P sequence of length m², has

    K_HH >= L = max(0, (m-1)^2
                     - sum_{i:P} [3 if i is an endpoint else 2]).

This inequality holds for **any arrangement** of the P labels, not just
periodic or randomly sampled sequences.

If t<=L, all C paths are event witnesses for A_t={K_HH>=t}, each with
Boltzmann weight >=q^L for q>=1. Let N be any independently checked upper
count on the *entire* unrestricted canonical SAW ensemble. For arbitrary
independent statewise positive reweighting a_s in [l,u]:

    P_q(A_t) >= l*C*q^L / (l*C*q^L + u*N*q^(t-1)).

The inequality is valid at every q>=1, and a bound evaluated at a
rational q0 automatically holds for **all real q>=q0** because L>t-1.

## Reproducible model-relative examples, n=64

All example P positions are interior chain indices (0-indexed).
The exact m=8 square Hamiltonian count is C=1,511,656,642; the
11-step-block SAW upper count is
N=5,231,553,077,432,424,492,299,000,247.

* P at indices 10 and 20, q=10: L=45 and threshold t=24.
  Conservative event population lower bound **0.9996540389518002**.
  Under correction factors [1/2,2], lower bound **0.9986175905866296**.
* P at indices 10,20,30,40, q=100: L=41 and threshold t=30.
  Lower bound **0.9999965392041905**, or **0.9999861569604857**
  with [1/2,2] correction factors.

This is a concentration claim on a *coarse event* for a synthetic 2D
H/P lattice model at stated strong coupling, NOT a statement about
specific folds, physically meaningful binding affinities, natural protein
sequence design, or a phase transition.

## Verification

A separate small-grid geometry DFS checks the contact-loss lower bound
against all square Hamiltonian paths for selected and exhaustive H/P
sequences. An unrestricted independent small-walk enumerator checks
that the event bounds never exceed the true finite-model Gibbs
probabilities. Invalid certificates, sequences and thresholds fail closed.

    python -m pytest tests/test_heterogeneous_ground.py -q

The certificate verifier is Python, not Lean-refined. The family count
rests on the #68 connected-frontier recurrence; see
[the exact ground-state document](EXACT_GROUND_STATE_DEGENERACY.md).
