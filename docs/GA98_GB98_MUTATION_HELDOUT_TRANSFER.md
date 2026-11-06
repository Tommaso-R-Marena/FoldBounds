# Mutation-heldout experimental transfer of FoldBounds' frozen fold-geometry codes

**Design freeze:** GA95/GB95 codes and thresholds from official 2009 RCSB
2KDL/2KDM NMR structures, already committed and independently checked at
\`experiments/results/ga_gb_error_code_verified_20261105.json\`.

**External evaluation:** 2012 deposited experimental NMR structures
2LHC (GA98, three-alpha fold) and 2LHD (GB98, four-beta-plus-alpha
fold), with 20 submitted complete 56-residue models each. GA98 and GB98
differ by a **single amino acid (Leu45 versus Tyr45)**. The
endpoints and classifier thresholds are **never refitted** or altered
using the new models.

Source:
- Original 2009 sequence/fold pair: https://www.rcsb.org/structure/2KDL
  and https://www.rcsb.org/structure/2KDM
- Independently deposited GA98: https://www.rcsb.org/structure/2LHC
- Independently deposited GB98: https://www.rcsb.org/structure/2LHD
- He et al., *Mutational tipping points for switching protein folds
  and functions*, Structure (2012), DOI 10.1016/j.str.2011.11.018.

## Precisely frozen protocol

1. Download original and external RCSB PDB coordinate files in a
   separate GitHub Actions run. The 2KDL/2KDM reference file bytes
   must match their previously published SHA-256 hashes.
2. Independently re-evaluate the source 1/3/5/7/9-channel reporter
   selection against only the reference training structures and
   require exact matching residue pairs, distance bands, g values
   and candidate counts with committed original results.
3. Only after freezing these codes, load the separately deposited
   2LHC/2LHD coordinate files. Recheck PDB IDs, all 20 complete 56-site
   C-alpha traces each, and every amino-acid residue against
   expected sequences. GA98 must equal GA95 I30F, and GB98 must
   equal GB95 A20L; they must differ only by Leu45/Tyr45.
4. Calculate C-alpha pair distances on the exact integer 0.001 Å
   coordinate grid, conservatively enclosing square roots by
   integer floor and ceiling. Reject measurements ambiguous at a
   fixed class-threshold rather than arbitrarily deciding.
5. Count unfiltered GA98/GB98 classification accuracy for all five
   frozen optimized codes. Never omit failed conformers.
6. Enumerate ALL subsets of up to f=(K-1)/2 reporter bits forced
   to vote for the **wrong class** in each test conformer. This is
   a genuinely adversarial stress test; it does not silently
   preserve reporters whose "extreme noise" still votes correctly.
7. Output explicit per-model, per-fold classification counts,
   conservative true-fold signed distance margins, and every
   corruption-trial count in machine-readable JSON.

## What can and cannot be concluded

A successful 2LHC/2LHD transfer would be evidence that FoldBounds'
selected **structural distance fingerprint** survives at least this
one independent *sequence variant pair* and its deposited NMR
conformer spread. It would be much more meaningful than in-sample
GA95/GB95 classification or odd/even-model holdouts.

However, this is an **independent deposited structure test**, not
a genuinely blind prospective test: the experimental GA98/GB98
fold identities and structures were published long before this
computational code was designed. The two structure pairs are
from a closely related designed protein family, with possible
methodology and refinement correlations across depositions.
One successful variant pair cannot establish generalizability
to unrelated fold switches or unseen functional conditions.

No physical dye-labeling precision is provided; pair distances
are C-alpha atom positions, not automatically measurable
fluorescence/labeled-protein observables. No population, kinetics,
Boltzmann weights, folding free energy, or sequence-to-fold
prediction is estimated here. The experiments measure whether
an already known fold *can be recognized* from a frozen set
of structural distance observables.

Conversely, a failed external conformer or a failed worst-case
error pattern is valuable negative evidence. It identifies
which frozen reporter loses its fit under a one-residue
mutation, and prevents transfer of an unjustified robustness
certificate to a new structural ensemble.

The algorithmic matching and majority decoding proofs are
Lean-checked in the parent PR #74; actual PDB byte decoding
and mutation-heldout coordinates remain Python checked and
RCSB-provenance guarded.

## Reproduce

The dedicated CI workflow
\`.github/workflows/ga98-gb98-external-transfer.yml\`
downloads four official structures, validates original file hashes,
runs independent tests, then executes:

    python experiments/run_ga98_gb98_transfer.py \
      --ga95 path/to/2KDL.pdb --gb95 path/to/2KDM.pdb \
      --ga98 path/to/2LHC.pdb --gb98 path/to/2LHD.pdb \
      --reference experiments/results/ga_gb_error_code_verified_20261105.json \
      --out frozen_transfer.json

An authentic new result is reportable only after a successful
real-data CI run. Any subsequent data-dependent summary must
quote the observed successes and failures rather than invent
the outcome in advance.
