import Mathlib

/-!
# Exact error-and-erasure safety margins for NMR fold fingerprint codes

A K-channel binary repetition code with at most E erased reporters
and T arbitrarily wrong remaining reporters is guaranteed to decode
the correct fold whenever 2*T+E<K, provided all other reporters
are correctly classified.

This is classical coding mathematics, not a new coding-theory theorem.
We formalize the strict quantitative assurance used by FoldBounds'
mutation-invariant pairwise C-alpha geometric reporter decoder.

The mapping of absent experimental residue sites to erasures relies
on the verified pairwise residue-disjointness and must be audited by
the Python source certificate checker. Protein atom coordinates and
external physics are NOT formalized in this module.
-/

namespace FoldBounds.MutationSafeCodes

/-- In an exact K-site binary reporter code, each erased vote counts
against the available majority, while a wrong vote counts twice. -/
theorem majority_under_errors_and_erasures
    {K E T good wrong : ℕ}
    (hPartition : good + wrong + E = K)
    (hErrors : wrong ≤ T)
    (hBudget : 2*T + E < K) :
    wrong < good := by
  omega

/-- A bound on missing residue labels dominates missing reporter
channels once each reporter uses a pair of entirely disjoint sites. -/
theorem site_erasure_bound_implies_channel_budget
    {K T missingSites erasedChannels : ℕ}
    (hErasure : erasedChannels ≤ missingSites)
    (hBudget : 2*T + missingSites < K) :
    2*T + erasedChannels < K := by
  omega

/-- The combined theorem shows that K distinct reporter channels
retain a strict correct-class majority after missing sites and wrong
remaining readouts, under explicitly declared disjointness premise. -/
theorem site_loss_plus_fault_majority
    {K lostSites erasedChannels T good wrong : ℕ}
    (hErasure : erasedChannels ≤ lostSites)
    (hPartition : good + wrong + erasedChannels = K)
    (hErrors : wrong ≤ T)
    (hBudget : 2*T + lostSites < K) :
    wrong < good := by
  exact majority_under_errors_and_erasures hPartition hErrors
    (site_erasure_bound_implies_channel_budget hErasure hBudget)

/-- At equality K=2T+E, the adversary CAN produce a tied vote:
there are T wrong votes, T correct votes, and E erasures.
Therefore the strict budget threshold cannot be relaxed in general. -/
theorem bound_is_sharp_at_tie (T E : ℕ) :
    ∃ good wrong : ℕ,
      good + wrong + E = 2*T + E ∧ good = wrong := by
  exact ⟨T, T, by omega, rfl⟩

/-- If all K original channels vote for the true class, and E
are erased and T of the remainder are flipped, the reconstructed
vote counts must give a strict correct majority under 2T+E<K. -/
theorem exact_adversarial_budget
    {K T E : ℕ} (hBudget : 2*T + E < K) :
    T < K - E - T := by
  omega

/-- A stronger composition: if some channels were ALREADY wrong before
adversarial faults, count them explicitly in the error budget. -/
theorem preexisting_errors_consume_budget
    {K E preexisting newFaults good wrong : ℕ}
    (hPartition : good + wrong + E = K)
    (hErrors : wrong ≤ preexisting + newFaults)
    (hBudget : 2*(preexisting + newFaults) + E < K) :
    wrong < good := by
  exact majority_under_errors_and_erasures hPartition hErrors hBudget

#print axioms FoldBounds.MutationSafeCodes.majority_under_errors_and_erasures
#print axioms FoldBounds.MutationSafeCodes.site_erasure_bound_implies_channel_budget
#print axioms FoldBounds.MutationSafeCodes.site_loss_plus_fault_majority
#print axioms FoldBounds.MutationSafeCodes.bound_is_sharp_at_tie
#print axioms FoldBounds.MutationSafeCodes.exact_adversarial_budget
#print axioms FoldBounds.MutationSafeCodes.preexisting_errors_consume_budget

end FoldBounds.MutationSafeCodes
