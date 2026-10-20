# Composition- and parity-controlled cubic HP folding class switches

**Research model only. Not experimental protein folding or alpha/beta switching.**

## Why this is a stronger control
On an anchored cubic self-avoiding walk, geometric contacts between
nonconsecutive vertices can occur only between opposite **chain index
parities**. Swapping H and P at two positions with the *same* parity
preserves (i) total H count and (ii) both parity-sublattice H counts.
A population switch after this controlled exchange cannot be explained
by merely changing those three gross hydrophobicity counts.

Event A is **geometric compactness**: at least two nonconsecutive contacts
for 8-residue walks, or at least three for 10-residue walks. Every H/P
sequence has the same finite conformation universe. The weight of a walk
is q^(number of HH contacts) * b^(number of bends).

Parameters range over **all** q in [8,12], b in [3/4,5/4].
Every full microstate may independently have its Boltzmann weight multiplied
by any factor in [5/6,6/5], with no statistical independence assumption.
These parameters are mathematically stipulated, not calibrated to proteins.

## Proof method and trust boundary
Enumerate the full anchored 3D cubic walk universe, retaining all contact
pairs and bends. Aggregate exactly into lossless contact-graph fingerprints.
Compute event/complement bivariate polynomials with nonnegative integral
coefficients. Transform both to tensor Bernstein bases on the parameter
rectangle. The weighted convex-ratio identity gives sound exact rational
enclosures of nominal event probability over the entire continuum rectangle.
Worst-case arbitrary statewise reweighting propagates event odds by the
ratio of allowed correction factors. A robust switch holds only when the
parent's corrected upper population is strictly below 1/2 and the mutant's
corrected lower population strictly above 1/2.

The abstract weighted-ratio and uncertainty lemmas are Lean-proved in
the parent research branch. **The concrete cubic enumeration, fingerprint
accounting, and Python Bernstein coefficient conversion are not yet
Lean-kernel proved end-to-end.**

## Results
| n | Parent sequence | Mutant (same number H) | Parent upper | Mutant lower | Distinct walks |
|---|---|---|---|---|---|
|8| PPPHPHHH | HPPHPHPH | 5460756/11934247 | 31554275/52243503 | 13,565 |
|10| PPPHPHHHPH | HPPHPHPHPH | 58746294144/133053135263 | 2071661600/4103373891 | 308,981 |

Both exchange positions **1 and 7 (1-indexed)**, preserving each parity
class H count. The guaranteed population gaps are about **14.64 and
6.33 percentage points**, respectively, but these are restricted to the
mathematical compactness event in this specified model.

Complete graph scan over all 256 n=8 H/P sequences:
768 distinct undirected same-parity composition-preserving exchanges;
12 robust switches, 6 equivalence classes after chain reversal.

Complete graph scan over all 1024 n=10 H/P sequences:
5120 distinct undirected exchanges; 4 robust switches, 2 reversal classes.

These count passing sufficient **certificates**, not all physically
possible true switches; the Bernstein bounds can be conservative.

## How much extra energy error can the certificate tolerate?
For symmetric multiplicative errors [1/s,s], s>=1, let U be the
Bernstein nominal parent upper and L the nominal mutant lower.
The **sufficient-certificate threshold** is

    s² < min((1-U)/U, L/(1-L)).

At n=8 the right-hand side is exactly 19420473/11376575,
so s < about **1.30654**. At n=10 it is exactly
2983192704/2031712291, so s < about **1.21174**.
The chosen s=6/5=1.2 passes, but only narrowly in the larger example.
Beyond this value the current *certificate* fails; this does not
prove a true physical transition fails.

## Reproduction and falsification
The separate test suite constructs a **fresh 3D DFS** without calling the
production geometry enumerator or fingerprint routine, checks every one of
13,565/308,981 full walks, and independently reproduces exact q=8, b=1
event populations for both sequences:

- 8-mer original 2516/7919, mutant 16152/22423.
- 10-mer original 95124/332659, mutant 1264000/1997213.

Additional tests enumerate all 768/5120 mutation edges, check exact rational
error thresholds, reject adversarially altered certificates, and reject
invalid composition/parity changes.

    python -m pytest -q tests/test_composition_swap.py
    PYTHONPATH=src python experiments/composition_swap_benchmark.py

## What is not shown
No experimental structure, mutation phenotype, hydration model, atomistic
free energy, folding kinetics, NMR validation or clinical utility is
inferred. Literature on metamorphic proteins and statistical-mechanical
epistasis is extensive. We do not claim novelty versus that literature.
To demonstrate a real breakthrough requires prospective blind prediction
of experimentally characterized **distinct 3D folds** (rather than this
coarse contact threshold), calibrated energetic uncertainties and rigorous
benchmarking against modern ensemble predictors.
