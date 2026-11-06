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


## Verified primary external-variant transfer: GA98 vs GB98

The first fully successful independent-mutant GitHub Actions run downloaded
the **exact official RCSB files** and parsed all 20 submitted conformers
each from 2LHC and 2LHD. Source SHA-256s (pinned in code and CI):

    2LHC 91ed5cf27675f68799e60951cca5df7a4ecb2e8d19fd2063c4babb300e9507ac
    2LHD 73131f50fdb7fc709e0c1b3e66daf64236a19d911a4bfe5324a23007847c2e6f

**No residue pair, reference interval, decision threshold or code
weight was refitted using either experimental test structure.**
This is the first genuine *variant-heldout deposited-structure*
evaluation of the already frozen FoldBounds reporter code, although
it is retrospective rather than prospective.

| Frozen reporter pairs K | Arbitrary wrong decisions f | Clean GA98/GB98 NMR model classifications | Correct adversarial subsets across both entries | Smallest strictly enclosed true-fold signed decision margin (Å) |
|---|---:|---:|---:|---:|
| 1 | 0 | 40/40 | 40/40 | 7.510 |
| 3 | 1 | 40/40 | 160/160 | 4.977 |
| 5 | 2 | 40/40 | 640/640 | 4.212 |
| 7 | 3 | 40/40 | 2,560/2,560 | 1.055 |
| 9 | 4 | 40/40 | 10,240/10,240 | 1.055 |

Thus **all 13,640 exhaustive adversarial binary-error patterns**
allowed by the five frozen code lengths were correctly classified on
all 40 external deposited NMR conformers.

The strict positive margins are **post-hoc** *finite reference-geometry*
results. They are not modeled physical dye/localization uncertainties.
For example, the 3-channel model has a conservative 4.977 Å
true-class threshold margin in the worst of the 40 external
conformers, and therefore a strict pairwise distance perturbation
smaller than 4.977 Å cannot reverse any individual correct reporter
vote *on those deposited conformers*. If both endpoint CA sites
each move by at most r, a sufficient geometric condition is
r < 4.977/2 = 2.4885 Å for this observed finite set, by the
reverse triangle inequality. Neither endpoint movement bound nor
label positional accuracy has been experimentally established.

All exact model-by-model results, SHA-256 source hashes, and
individual worst-case wrong-vote pattern counts are retained in
\`experiments/results/ga98_gb98_frozen_transfer_20261105.json\`.
The workflow re-downloads all four structures and fails closed
unless the recomputed file agrees with the retained data fields.

**Interpretation:** successful transfer across a one-amino-acid
fold-switch pair within the very same engineered protein family.
It supports generalization of these *structural readout choices*
beyond the original GA95/GB95 NMR depositions, not sequence-based
protein-fold prediction, equilibrium population inference, or
generalization to unrelated protein families.

## Second external challenge: a distinct reciprocal mutation cascade

The original experimental study additionally resolved 2LHG
(GB98-T25I, **3α fold**, ten submitted conformers) and 2LHE
(GB98-T25I,L20A, **4β+α fold**, twenty conformers). The T25I
mutation induces a fold switch back to a helix-rich structure,
and the subsequent L20A mutation switches back to a beta-rich
structure (He et al., Structure 2012; PMID 22325777).
This is a more challenging stress test because position 25
is **itself one of the originally selected reporter sites** in
some codes. The baseline pair choices and thresholds remain
frozen from GA95/GB95.

All exact results from this follow-up must be read from the new
\`reciprocal-mutant-nmr\` CI job, not assumed to be a success.
The planned test reports every classified structure, wrong reporter
bit, and adversarial stress pattern, with no reoptimization.


## Verified reciprocal second-cascade results: 2LHG and 2LHE

The independent reciprocal-mutant GitHub Actions run **successfully**
re-downloaded, parsed and audited the official follow-up files
and their published fold labels (He et al., Structure 2012,
https://pmc.ncbi.nlm.nih.gov/articles/PMC3278708/):

    2LHG 946bb4d44567de7725aeecd334d7c1cec086737d6445d157cd94976326c84cff
    2LHE cd2fa17ac762106c4f21c464bd7709cb612e3cd96d227ee18fdccd99608e83d1

Each hash is pinned independently in BOTH Python and GitHub Actions.
The same 2009 GA95/GB95 reporter pairs and thresholds were used
without any refitting for these **30 additional deposited NMR models**.
Their experimental labels are independently established: 2LHG
(GB98-T25I) is three-helix GA-like, while 2LHE
(GB98-T25I,L20A) is beta-rich GB-like.

| Frozen code length K | Correct independent NMR conformers | Correct exhaustive wrong-vote subset patterns | Smallest true-fold reporter threshold margin (Å) |
|---|---:|---:|---:|
| 1 | **30/30** | 30/30 | 6.961 |
| 3 | **30/30** | 120/120 | 0.570 |
| 5 | **30/30** | 480/480 | 0.570 |
| 7 | **30/30** | 1,920/1,920 | 0.570 |
| 9 | **30/30** | 7,680/7,680 | 0.570 |

The experimental conformer whose distance approaches closest to
its trained threshold has a margin of only **0.570 Å** in the
three-channel (and longer) codes, substantially smaller than
the corresponding GA98/GB98 margins. This is an important sign
that success on the current deposited conformers does NOT validate
the original 6.474 Å reference-model measurement tolerance on
arbitrary related mutants. Additional drift of the reporter
distance by more than 0.570 Å in the unfavorable direction
could flip the closest individual channel decision.

All 30 models and all **10,230** allowed adversarial wrong-vote
patterns passed for this additional mutation cascade. The
unfiltered per-model outcomes and individual reporter error
diagnostics are committed at
\`experiments/results/ga98_t25i_followup_20261105.json\`.
The dedicated workflow requires exact regeneration of this file
from the four authenticated RCSB sources.

### Total external experimental-structure transfer evidence

Taken together, the **four distinct 2012 PDB entries**
2LHC (20 models), 2LHD (20 models), 2LHG (10 models),
and 2LHE (20 models) provide **70 heldout deposited
conformers**, all recognized by every one of the five
GA95/GB95-derived code lengths with no refitting. Across them,
all **23,870** enumerated allowed wrong-channel corruption
subsets were correctly decoded.

This is a useful, concrete generalization within a closely
related designed-protein family. It does NOT constitute 70
independent biological replicates, a blind structure-prediction
challenge, or a protein-folding mechanism breakthrough.

A meaningful next independent scientific gate is to measure
real state-conditional reporter distributions on new constructs,
verify label-site neutrality, and evaluate fully prospective
variants. A single C-alpha geometric difference by itself is
not a calibrated FRET observable or a thermodynamic population
prediction.
