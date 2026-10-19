# Aristotle lattice verification report

This report covers the 2026-10 extension of FoldBounds: the representation
theorem, contact semantics, compact frontier certificates, witness lower
bounds, and the transfer relaxation. Each item lists the exact Lean names, the
statement, the prerequisites, and the limitations. Every theorem named here
compiles with Lean 4.28.0 / Mathlib v4.28.0. None contains `sorry` or
`admit`, and none adds an axiom. All of the headline theorems are in
`FoldBounds/AxiomAudit.lean`, which now audits 117 declarations; each depends
on a subset of `propext`, `Classical.choice`, `Quot.sound`. See
`docs/ARISTOTLE_VERIFICATION_RUN_REPORT.md`.

New modules, all under `FoldBounds/`:

| Module | Content |
|---|---|
| `LatticeEquivalence.lean` | direction words ≃ canonical SAWs, sum preservation, nonemptiness |
| `LatticeContacts.lean` | checkerboard parity, Lean/Python contact equivalence, `q^k` monotonicity |
| `Frontier.lean` | completion sets, tree partition, `3^(n-m)` counting bound, contact decomposition and caps |
| `FrontierDegree.lean` | degree lemma (≤ 2 or 3 non-chain contacts per residue), degree-capped contact caps |
| `FrontierFamily.lean` | certified diagonal families of completions: `2^r` lower bounds without enumeration |
| `Transfer.lean` | generic merged transfer DP soundness, lattice relaxation, domination of genuine completions |
| `FrontierCheck.lean` | certificate trie `FTree`, the checker `checkNode`/`check`, per-leaf soundness |
| `FrontierSoundness.lean` | end-to-end `check_sound`, real-valued probability form, compactness witness |
| `FrontierPartition.lean` | leaf partition and unique cover, transfer ≤ degree ≤ basic hierarchy |
| `FrontierFlat.lean` | flat (Python-style) certificates, strict trie builder, `checkFlat_sound` |
| `FrontierCLI.lean` | text parser and JSON printer (**unverified boundary**, no theorems) |
| `FrontierTests.lean` | kernel-checked accepted and rejected examples (`decide +kernel`) |

## 1. Representation theorem (Priority 1, done)

* `LatticeModel.CanonicalSAW steps`: site lists `path` with
  `path.length = steps + 2`, `path.take 2 = [(0,0), (1,0)]`, consecutive
  sites grid-adjacent (`adjacentChain`), and `path.Nodup`.
* `LatticeModel.stepDir`, `decodeSites`, `pathDirections`: the computable
  inverse, which reads each unit step back as one of the four direction
  constructors.
* `LatticeModel.conformationEquivSAW steps : Conformation steps ≃ CanonicalSAW steps`.
  The forward map is `walkFrom ∘ List.ofFn`, and the two round-trip proofs are
  `pathDirections_walkFrom` and `walkFrom_pathDirections`.
  `conformationEquivSAW_apply` shows that the forward map is literally the
  existing decoder.
* `toSAW_injective` and `toSAW_surjective` are the injectivity and
  surjectivity of the decoder.
* `sum_conformation_eq_sum_saw`: `∑ c, w (toSAW c) = ∑ P : CanonicalSAW, w P`
  for any `AddCommMonoid`-valued weight. `card_conformation_eq_card_saw` gives
  equal cardinalities. `massQ_weightQ_eq_sum_saw` specialises this to
  FoldBounds' exact rational contact weights.
* Preserved structure: length, the fixed anchor and first bond, adjacency,
  distinctness, chain order (the list order itself), and orientation (absolute
  directions, nothing quotiented). Multiplicity is preserved because the map
  is a bijection.

**Nonemptiness.** `walkFrom_straight_nodup`, `straight`, and the instance
`conformation_nonempty : Nonempty (Conformation steps)` hold for every
`steps`. `checkHP_coupled_sound_unconditional` is the existing
`checkHP_coupled_sound` with the `[Nonempty]` premise discharged. The
existing theorem statement was left unchanged for compatibility; it now
fires automatically through instance resolution.

## 2. Contact semantics (Priority 2, done)

* `colour_walkFrom`: the checkerboard colour `(x+y) mod 2` of site `i` is
  `(i+1) mod 2`.
* `index_parity_of_adjacent`: lattice-adjacent sites of a decoded walk have
  chain indices of opposite parity. `not_adjacent_of_even_gap` and
  `not_adjacent_of_gap_two` follow: an even index gap, and in particular a
  gap of 2, can never form a contact.
* `candidatePair`: Python's `_possible_hh_pairs` (both H, `i + 3 ≤ j`,
  `(j − i)` odd). `parityContactCount` is Python's `_contacts(path, pairs)`.
  It is also defined on prefixes, counting only pairs with `j < len(path)`,
  as Python does.
* `contactCount_eq_parityContactCount`: for every decoded walk whose length
  matches the sequence, the Lean scorer (`j ≥ i + 2`, both H, adjacent) equals
  the parity scorer. `saw_contactCount_eq_parity` states the same for
  `CanonicalSAW`.
  *Indexing audit:* the Lean scorer already excludes `j = i + 1`, and
  `j = i + 2` is excluded geometrically by parity. Each unordered pair is
  counted once because `i < j` in both scorers. A length mismatch is excluded
  by hypothesis. No discrepancy between the two scorers was found.
* `weightQ_eq_of_contacts_eq`, `weightQ_eq_parity`, `qpow_mono`
  (`1 ≤ q → k₁ ≤ k₂ → q^k₁ ≤ q^k₂`), `qpow_pos`, and `weightQ_mono` cover the
  weights. They concern exact **rational** weights `q^k`. The real Boltzmann
  model `exp(−βE)` from the original library is separate, and the
  identification `q = exp(βε)` is not formalised.

## 3. Conformation tree (Priority 3, done)

Prefix words `p : List (Fin 4)` index tree nodes. `Frontier.ext steps p` is
the finite set of complete valid words extending `p`. It is a semantic object
that the checker never computes.

* `ext_nil`: the root's completion set is the whole conformation space. It is
  identified with `Conformation` through `validWords`, `mem_validWords`, and
  `sum_conformation_eq_sum_validWords`.
* `ext_eq_biUnion` and `ext_pairwiseDisjoint`: for a nonterminal prefix, the
  completion sets of its one-step children are pairwise disjoint and their
  union is `ext p`. `sum_ext_split` is the additive form.
* Blocked prefixes: `ext_eq_empty_of_not_nodup` and
  `ext_eq_empty_of_length_gt`. Terminal prefixes: `ext_eq_singleton`.
* `not_nodup_reverse`: an immediate reversal always collides.

## 4. Non-enumerative counting and contact bounds (done)

* `Frontier.card_ext_le_pow_three : (ext steps p).card ≤ 3 ^ (steps − |p|)`.
  The proof is by induction through the child partition and
  `not_nodup_reverse`. It is about the exact conformation space, not about
  unrestricted words.
* `contacts_decompose`: total contacts = known prefix contacts + cross
  contacts (placed–unplaced) + future–future contacts.
* `contacts_le_kmax`: `k(c) ≤ known(p) + futurePairs(m)`, which is Python's
  `K_max`. `knownContacts_le` gives the lower side.
* `contactsAt_le_slots` (the **degree lemma**): in a complete SAW, residue
  `j` touches at most 2 earlier candidate residues, or 3 if it is the last
  residue. `contacts_le_degree_caps` gives
  `k(c) ≤ known + futureCap ≤ known + futurePairs` (`futureCap_le_futurePairs`)
  and the cross-score refinement `k(c) ≤ known + cross(c) + ffCap`.

## 5. Transfer relaxation (Priority 5, done for one relaxation level)

* Generic part: `Transfer.wsum_dp`. For any finitely branching transition
  system and any comparator used to sort and merge equal states, the merged
  forward DP equals the backward path-sum specification `spec`. This is
  transfer-matrix soundness, and merging is proved exact. `spec_nonneg`
  shows the specification is nonnegative.
* Lattice part: `RState = (site, last dir, cross score, next index)` and
  `relaxedNext`, which takes non-reversing steps that never land on a fixed
  prefix site and may revisit new sites. `crossContacts_snoc` shows that
  appending a site adds exactly its cross contacts.
* `relaxed_sound` and `sum_ext_le_dp`: for any nonnegative `g`,
  `∑_{l ∈ ext p} g(cross(l)) ≤ wsum g (relaxedTable p)`. This is the
  inclusion of genuine completions into relaxed walks, combined with the
  cross-score correspondence.
* `relaxed_le_basic`, `basicZ_degree_le_basic`, and `transfer_total_le_basic`
  form a **hierarchy theorem**: the transfer bound ≤ the degree-capped
  envelope ≤ Python's basic envelope `3^r q^Kmax`. Enlarging the memory
  further (forbidding short loops in the suffix, for example) was **not**
  formalised. A monotone hierarchy beyond these three levels is a research
  direction, not a theorem.

## 5b. Certified families for lower bounds (§6.3, done for one family class)

`FoldBounds/FrontierFamily.lean`:

* `pot a b s = (±x) + (±y)`, `dirX`, `dirY`: the diagonal potential and the
  two steps that increase it by exactly one (`pot_advance_dirX`,
  `pot_advance_dirY`).
* `trace_diag`: along any word over those two steps, the potential strictly
  increases. So the new sites are pairwise distinct and all lie strictly
  above the starting potential.
* `diagOK a b p` and `famOK p`: the locally checkable condition that the
  prefix end maximises the potential over the prefix walk.
* `nodup_family`: under `diagOK`, every family continuation is
  self-avoiding. `card_family` shows the family has exactly `2^r` members,
  and `family_subset_ext` shows they are genuine completions.
* `family_lower`: `2^r · c ≤ ∑_{l ∈ C(p)} F l` whenever `F ≥ 0` and
  `F ≥ c` on `C(p)`. This lower bound requires no enumeration.
* Checker integration: `FrontierCheck.famLow`, `famRaise`, `famLow_le`, and
  `famRaise_sound`. In transfer mode, bound leaves take the maximum of the
  witness lower bound and `2^r q^known`; the event and complement lower
  bounds are raised only when the classification is certain.

This family class is simple, and its count `2^r` is much smaller than the
`≈ 2.64^r` growth of all completions. Partially directed families
(`≈ 2.41^r`) and contact-weighted family DPs were **not** formalised.

## 6. The checker and its soundness (Priority 3/4, done)

* `FrontierCheck.check : List Bool → ℚ → ℕ → FTree → Option FrontierResult`;
  see `docs/FORMALLY_VERIFIED_FRONTIER_CERTIFICATE.md`.
* `exactEnc_sound`, `boundEnc_sound`, `deadOK_ext_empty`, and
  `checkNode_sound` give soundness of each accepted subtree enclosure, by
  structural induction on the trie.
* **`check_sound`** (flagship): for every sequence `seq : Fin (steps+2) → Bool`,
  `q : ℚ`, `thr : ℕ`, and certificate `t`, `check (List.ofFn seq) q thr t =
  some r` implies that `Z`, `Z_A`, `Z_¬A`, and `Z_A/Z` (all about the Lean HP
  model `weightQ`/`contactEvent`) lie in the returned intervals. There are no
  further hypotheses.
* `check_coupled_probability_sound` is the real-valued coupled form, matching
  `checkHP_coupled_sound`. `checkFlat_sound` covers flat certificates.
* Witness lower bounds (Priority 4) are inside `boundEnc_sound`. Verified
  distinct witnesses in one leaf contribute their exact weight sum to the
  lower bounds of total, event, and complement, and tighten the upper bounds.
  `wit_subset` shows the witnesses lie in the leaf's completion set.
  Disjointness across leaves follows from `leaf_cover_unique` and
  `sum_ext_eq_sum_leaves`.
* `check_root_bound_accepted`: a single bound leaf is accepted for every
  admissible model (compactness). `Z_le_basic_root` gives
  `Z ≤ 3^steps · q^{futurePairs}` for every model, with no enumeration.
* `FrontierTests.hhhh_probability_exact`: an exact, kernel-verified
  `P(k ≥ 1) = 4/11` for HHHH at `q = 2`, obtained by applying `check_sound`
  to an accepted certificate.

## 7. Limitations

* The theorems concern the Lean definitions. That the Python code computes
  the same quantities is supported only by differential tests (see the run
  report), not by a refinement proof.
* The text parser and the command-line runner are unverified. Executing
  `check` natively trusts the Lean compiler; only `decide +kernel` examples
  avoid that.
* There is no formal complexity theorem; §7 of the certificate specification
  gives an informal cost analysis.
* The bounds are sound but can still be loose for long chains (see
  `docs/NONENUMERATIVE_BOUNDS_RESEARCH_REPORT.md`).
* Python's isoperimetric, bounding-box, and matching caps and its adaptive
  refinement are **not** formalised. The Lean checker uses different caps
  (degree caps and the transfer DP), and the relative strength of the two
  sets of caps has not been compared systematically.
* Lower bounds come from exact leaves, individually verified witnesses, and
  the single certified diagonal-family class (§5b). That family is
  exponentially smaller than the full completion set, so lower bounds on long
  chains remain loose. Richer families were not formalised.
