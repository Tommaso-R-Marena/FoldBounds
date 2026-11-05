# Local rewiring: a fail-closed structural-confidence ceiling

## Status / scientific claim

This is a **model-relative non-identifiability certificate** for any supplied
*complete, canonical 3D simple-cubic self-avoiding walk* under a fixed positive
**pairwise nearest-neighbor contact-weight matrix**. It **upper bounds the
equilibrium probability of a particular proposed conformation**, without
enumerating the remaining exponentially many conformations.

It is not a guarantee that a natural protein is misfolded, not an
experimental failure of AlphaFold or other contemporary models, not a
general statement about experimentally observable backbone structures,
and not an end-to-end kernel proof of the 3D geometry.

The key distinction from the other FoldBounds workstreams:

- A *high-contact event* can have probability >99%, while any particular
  microscopic structure has extremely low equilibrium probability.
- A *structure-level confidence ceiling* is therefore an essential
  counterweight to misleading claims based on contact-count concentration.
- The method applies to **arbitrary finite residue alphabets** and an
  **explicit positive rational matrix** of contact Boltzmann factors, not
  just the all-H HP toy model. No energy-attraction or q>=1 assumption is
  required.

### Prior art

Highly degenerate HP landscapes have been extensively documented, including
Irbäck & Troein, *Enumerating Designing Sequences in the HP Model* (2002),
https://pmc.ncbi.nlm.nih.gov/articles/PMC3456819/, and Schiemann, Bachmann
& Janke, *Exact enumeration of three-dimensional lattice proteins* (2005),
https://doi.org/10.1016/j.cpc.2004.09.011. A practical distinction between
model equilibrium probability and simulation hit frequency is described in
Cox & Johnston (2006), https://pubmed.ncbi.nlm.nih.gov/16774371/.

Local conformational moves, finite-block transfer and degeneracy are
therefore NOT new discoveries. The implemented contribution is a
**reproducible, conservative, exact-rational, arbitrarily reweighted
confidence ceiling for a user-supplied arbitrary sequence and 3D target**.
It packages the reasoning into a independently testable, fail-closed
certificate. Literature novelty and potential publication worthiness
require independent review.

## The checked model and trust boundary

Input is an untrusted

    StructuralModel(sequence, target, pair_weights, correction_lower, correction_upper)

where target is a complete sequence-indexed 3D lattice SAW starting from
(0,0,0)->(1,0,0), pair_weights is an explicit positive **Fraction** for
each unordered pair of residue types in the sequence (symmetry is enforced
by canonical sorted keys), and every actual microstate weight is its
nominal lattice-contact weight multiplied by an arbitrary factor in
[correction_lower, correction_upper], with 0<lower<=upper.

The checker independently validates:

- Full length, canonical root, exact integer triple coordinates.
- Nearest-neighbor chain steps, global self-avoidance, complete alphabet
  interaction-pair coverage, positive exact rational inputs.
- Every local 2x2x2 cube occupied by **eight consecutive residues**,
  with exact endpoint constraints.
- Every alternate local Hamiltonian path, reconstructed from the small
  complete 144-path motif library already independently tested in the
  cubic family research.
- Both internal and external *noncovalent* contact factors: bridge edges
  remain covalent and are excluded from the energetic score.
- Pairwise spatial **nonadjacency** of all selected cubes. No contact
  edge can connect two selected cubes, ensuring that their relative
  Boltzmann weight changes are truly independent.

An infeasible or malformed input is rejected. A valid target with no
recognized local rewiring gets an uninformative ceiling of 1, **not a
fabricated guarantee of uniqueness**.

## Mathematical certificate

Take a candidate 8-residue interval [i,i+7] occupying all eight sites of
an axis-aligned 2x2x2 cube. Replace its internal seven bonds by any
Hamiltonian path through the same eight sites whose start and end
coordinates match the original interval (and, if i=0, whose first two
coordinates also equal the canonical global root).

Every replacement produces a valid complete walk because:
1. The site set and all exterior coordinates are unchanged.
2. The first/last anchor coordinates of the interval are unchanged,
   preserving the two exterior covalent chain bonds.
3. The interior is again a Hamiltonian unit-step path through exactly
   the same distinct eight vertices.
4. Nonselected site occupancy cannot collide with a replacement because
   the entire interior site set remains identical.

For each candidate patch i, enumerate all such local paths
v∈V_i. Let W(P) be the nominal full-conformation Boltzmann weight:

    W(P) = product_{noncovalent adjacent site pair {j,k}} w[seq[j],seq[k]]

where all w>0 and are exact rationals. Compute exact relative weights

    r_i(v) = W(P with patch i replaced by v) / W(P).

The program **does not evaluate this huge full weight** to get these
ratios: it multiplies only contact factors involving at least one site
in the patch, with every affected lattice edge counted exactly once.
All other factors cancel, and the original path has r_i(original)=1.

Select a subset of patches whose occupied cubes are **pairwise not
adjacent on the lattice**, not just non-overlapping in chain index.
Then changing patch i cannot alter any contact factor involving any
selected patch j. This removes all interaction cross terms, proving

    W(P with all chosen patch variants) / W(P)
       = product_i r_i(v_i).

Different local choices yield different globally indexed coordinate
paths. Summing over the explicitly constructive family gives the
**exact** relative family partition factor

    S = product_i [sum_{v in V_i} r_i(v)]   >= 1.

For the whole (unrestricted) conformational ensemble, every
microstate may have an independently adversarial correction a_s ∈ [l,u].
The target has weight at most u W(P). The other verified family
members collectively have weight at least l W(P)(S−1); any other
states only INCREASE the denominator. Hence

    P(target) <= u / [u + l*(S−1)].

The guarantee is **not statistical** and does not rely on a
monomer-independence assumption: it is exact finite combinatorics.

The checker greedily chooses a set of profitable, nonadjacent cubes;
the greedy schedule is NOT assumed to be optimal, and the resulting
ceiling is sound even if a much stronger selection exists.

## Structural symmetry

The above ceiling concerns an **exact anchored, indexed lattice
microstate**. Many physical discussions identify structures up to
rigid rotations and reflections. The full three-dimensional cubic
point group contains at most **48** signed coordinate permutations,
so an indexed target's congruence class has at most 48 distinct
anchored microstates. Every such transformation preserves the
residue-labeled contact energy.

If S>48, the constructive family's nominal mass **outside** the
target's rigid congruence class is at least (S−48)*W(P), and
the full congruence-class weight is at most 48*u*W(P). Thus

    P(rigid fold class) <= 48*u/[48*u+l*(S−48)].

This intentionally uses the conservative group size 48; the actual
canonical-root residual group is smaller. If S<=48 the checker returns
the trivial fold-class ceiling 1, not an unjustified inference.

## Exact energy-degeneracy subset

For each selected patch the checker also counts local variants with
r_i(v)=1. Any combination of those alternatives has *exactly* the
original model energy. If d_i is that count, the target has at least

    D = product_i d_i

distinct **isoenergetic** conformations in its constructible class.
When all cube residues share a common interaction type, EVERY local
permutation has the same contact energy for ANY fixed positive
species-pair contact matrix, because the cube has the same seven
covalent internal edges, occupies the same sites, and keeps its
external bridge locations.

For an uncorrected Boltzmann model, this yields

    P(target) <= 1/D  for ALL temperatures and all positive q.

With adversarial factors [l,u], the exact degeneracy-only upper is
u/[u+l*(D−1)]. The full exact-rational S typically strengthens it.

### Model resolution requirement

If D distinct conformers have identical baseline energy, and an
additional interaction is permitted to make a favored target
at most R times as heavy as **each** of its D−1 alternatives,
then even with optimal discrimination the favored target's
probability cannot exceed

    1/[1+(D−1)/R].

To achieve a desired target probability p, such a **uniform upper
allowance on target-vs-alternative weight ratios** must be at least

    R_required = (D−1)*p/(1−p).

Equivalently, any additional energy model whose relative target
advantage against *all* D−1 alternatives is bounded by Δ in kBT
units must permit at least Δ >= log(R_required) for that confidence
to be attainable. This is a necessary **information-resolution**
condition, NOT a sufficient design criterion.

## Critical independent checks

Run from the FoldBounds repository root:

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_local_rewire_confidence.py
    python experiments/local_rewire_confidence_demo.py

The test suite:
- Independently rebuilds the 144 ordered 8-vertex Hamiltonian paths.
- Independently calculates FULL sequence-labeled noncovalent contact
  weights for every tested alternate, including synthetic 20-letter
  alphabets and positive arbitrary fractional pair matrices.
- Exhaustively verifies the factorization for two nonadjacent
  2x2x2 cubes by enumerating all their joint alternatives.
- Completely enumerates all 13,565 canonical eight-residue 3D SAWs
  and proves the checker never overclaims the maximal **adversarially
  corrected** equilibrium probability of the tested target.
- Checks pure H chains at several q values, including q<1, and tests
  approximate/invalid inputs and an **abstaining straight-chain case**.

### Limitations

This remains a simple-cubic lattice nearest-neighbor pairwise
Boltzmann model. No steric side-chain geometry, hydrogen-bond
directionality, loop torsion energies, kinetic selection,
solvation model, or experimentally calibrated temperature is
established. A conformational probability ceiling is NOT a
prediction of protein folding rates, experimental structural
heterogeneity, ensemble populations, or a failure of AlphaFold.

The code does not prove the geometry in Lean. The new Lean theorem
proves the conditional rational mass/probability inference only.
The important future step is to generalize to approximate geometric
motifs in physical coordinate ensembles *with a justified energy-error
certificate*, and to benchmark the actual confidence bounds against
blinded experimental structural ensembles. Do not merge this PR before
the CI, independent exact-oracle tests, and peer review are complete.
