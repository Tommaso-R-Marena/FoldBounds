# Proof-gated JEPA and JEPA-RL — experimental research branch

**Status (7 November 2026):** opt-in model and tests, not merged, not a
validated new protein-folding algorithm. This branch is additive, based on
`main` at `e373a422`; it intentionally does not rewrite the 61-PR draft
integration tree in PR #91. Reconcile its imports after that integration is
independently reviewed.

## Two explicitly separate models

### A. Certificate-search JEPA / JEPA-RL

A finite HP frontier state `s` consists of a sequence, exact rational
Boltzmann factor `q`, threshold, complete frontier, and independently
verified current event-probability interval. Its feasible actions are **only**
`refine_frontier(cert, bounded_path)`. The model learns

\[
z_t=f_\\theta(s_t),\qquad
\\widehat z_{t+1}=g_\\phi(z_t,a_t),\qquad
z_{t+1}^{\mathrm{target}}=f_{\\bar\\theta}(s_{t+1}).
\]

The target encoder follows an exponential moving average (EMA). Auxiliary
next-state reconstruction reduces collapse risk. Offline training ranks
refinements using **independently verified** improvements. Online RL samples
an action and uses REINFORCE with a learned critic and realized reward

\[
r_t = \\frac{(w_t-w_{t+1}) +
  0.01\\tanh(Z_t^{\\mathrm{width}}-Z_{t+1}^{\\mathrm{width}})}
 {\\text{checker calls charged to the action}}.
\]

Widths and gains are computed with exact `Fraction` arithmetic before
conversion to untrusted floating-point ML labels. The partition-width term
is bounded by tanh to prevent excessive reward scale. Neither the learner
nor the reward equation proves the underlying probability enclosures.

**Boundary:** `refine_frontier`, `verify_lookahead`, and
`verify_refinement_run` remain authoritative. The actor's scores are
advisory and checked for shape and finite values. The full frontier audit
is replayed with no ML. A bad model can waste computation, but cannot
rewrite a certificate or a bound. This is a Python software assurance
claim only; Lean does not yet refine these Python checkers end to end.

### B. Time-conditioned trajectory JEPA

`jepa_trajectory.py` consumes actual chronological, fixed-topology MD
frames, `time_ps[T]` and `xyz_angstrom[T,N,3]`, with replica/source IDs.
It builds 10 translation- and rotation-invariant structural observables
(contact-distance fractions, local distances and radius of gyration),
then uses the same latent-prediction/EMA machinery conditioned on the
**physical time gap** instead of a computational expansion action.
Training includes feature masking and future-observable reconstruction.

**Boundary:** This is representation learning, **not** a certified dynamical
model. It does not calculate reliable transition probabilities, free-energy
landscapes, folding rates, rare-event kinetics, or medical affinities.
A stochastic MD trajectory cannot be certified by one deterministic
latent forecast. Do not train the temporal model on arbitrarily ordered
NMR ensemble models. No experimental or genuine OpenMM/GROMACS trajectories
were used to validate this implementation in this PR.

## Reproducibility

```bash
python -m pip install -e '.[test]'
python -m pip install torch numpy
python -m pytest -q tests/test_jepa_guidance.py tests/test_jepa_trajectory.py
python experiments/benchmark_jepa.py --train --rl --seed 7 --steps 5
```

The PyTorch package is optional; installing the base FoldBounds package
does not import it. `collect_checked_decisions` charges one verifier call
per candidate, and the benchmark reports its training-data-check cost
separately from deployment. Both `guided_refine` and the probability
greedy comparator report final **audited** bounds. Measured wall-clock
includes the model/scorer and a final independent audit replay.

For genuine MD data, export .npz with arrays `time_ps` (ps),
`xyz_angstrom` (Å) from consistent, topology-validated trajectories,
then call:

```python
from foldbounds.jepa_trajectory import (
    load_npz, fit_trajectory_jepa, evaluate_trajectory_jepa,
)
train = load_npz("md_replica_1.npz", source_id="train-1")
heldout = load_npz("md_replica_2.npz", source_id="heldout-2")
model, losses = fit_trajectory_jepa((train,), epochs=40)
print(evaluate_trajectory_jepa(model, (heldout,)))
```

Compare against the **persistence baseline** on wholly held-out replicas;
also compare protein-held-out evaluations for any generalization claim.
A lower training loss is not evidence of improved MD forecasting.
Keep train/test trajectory IDs disjoint, including frames linked by
near-duplicate starting conformations. Record force field, protonation,
thermostat, integration timestep, solvent/ions, equilibration, replica
seed, source hashes, frame spacing, protein/PDB IDs, and leakage checks.

## Remaining scientific gates (NOT completed)

1. Run the actual FoldBounds full pytest and pinned Lean audit on the research
   branch, then check combined import/test compatibility with draft PR #91.
2. Freeze and run a larger prospective held-out HP benchmark with matched
   **verifier calls** and wall-clock (the greedy oracle's scoring cost counts).
   Compare random, breadth, mass, probability-greedy, supervised-only,
   JEPA-only, and JEPA+RL across multiple seeds; publish failed cases.
3. Evaluate MD forecasting on real OpenMM/GROMACS trajectories and compare
   persistence, supervised MLP, contact-only and latent-dynamics models at
   identical data and compute budgets. Split by independent replica and
   ideally protein; do not claim transferable kinetics without that test.
4. Integrate credible external experimental data and uncertainty models before
   investigating closed-loop simulation/experiment selection. Any new bound
   must still be checked by a separate trusted kernel.
5. Consider a formal Lean refinement layer for the discrete action/replay
   protocol. The neural policy itself need not be verified if it remains
   outside the trust boundary.

**No scientific breakthrough or empirical neural advantage is claimed.**

## Export actual MD simulations (optional MDTraj)

The dedicated exporter accepts time-ordered XTC and DCD molecular-dynamics
trajectories and an explicit PDB topology. It records SHA-256 provenance,
replica ID, selection, coordinate/time units and frame counts, and rejects
unordered PDB/NMR ensembles as fake time series. MDTraj coordinates in nm
are multiplied by 10 to obtain the advertised Å array.

```bash
python -m pip install mdtraj numpy
python scripts/export_jepa_md.py \
  --trajectory replica_01.xtc \
  --topology topology.pdb \
  --output replica_01.npz \
  --source-id proteinA-replica01 \
  --selection 'protein and name CA' \
  --stride 10
```

For DCD, pass `--timestep-ps` specifying the **original frame spacing**
in ps. The exporter refuses DCD without this value, rather than inventing
simulation times. Run against *actual independently equilibrated dynamics*,
with replica-level isolation, force-field/protonation details and separate
equilibrium/nonequilibrium validation. The exporter and its synthetic format
test prove only that units and time stamps are transmitted as intended.
