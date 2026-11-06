# Formally verified frontier certificates

This document says exactly what a FoldBounds frontier certificate claims,
what the Lean checker checks, and what follows once the checker accepts.
The checker is `FoldBounds.FrontierCheck.check` in
`FoldBounds/FrontierCheck.lean`. The flat-list front end is
`FoldBounds.FrontierCheck.checkFlat` in `FoldBounds/FrontierFlat.lean`.
Soundness is proved by `check_sound` and `checkFlat_sound`, and both are
in the axiom audit.

## 1. The model these certificates are about

* **Sequence.** `S : List Bool` of length `n = steps + 2 ≥ 2`, where `true` means H.
* **Conformations.** `LatticeModel.Conformation steps` is the set of direction
  words `Fin steps → Fin 4` (directions `0:+x 1:-x 2:+y 3:-y`) whose decoded
  walk `walkFrom` starts at `(0,0) → (1,0)` and visits no site twice. By
  `conformationEquivSAW` this set is in explicit bijection with the canonical
  anchored self-avoiding paths `CanonicalSAW steps`. Reflections and rotations
  are **not** quotiented out, which matches the Python state-counting
  convention.
* **Contacts.** `contactCount S W` counts unordered pairs `i < j` with `j ≥ i+2`,
  both residues H, at lattice distance 1. `contactCount_eq_parityContactCount`
  shows this equals the Python-style parity-filtered count (gap ≥ 3 and odd).
* **Weights.** `weightQ = q ^ contactCount`, with exact `q : ℚ` and `q ≥ 1`.
* **Event.** `contactEvent thr` = conformations with at least `thr` contacts.
* **Masses.** `Z = massQ weightQ univ`, `Z_A = massQ weightQ (contactEvent thr)`,
  `Z_¬A = massQ weightQ (contactEvent thr)ᶜ`, and `P(A) = Z_A / Z`.

These are rational combinatorial weights. They become a physical Boltzmann
model only through the separate identification `q = exp(βε)`, which is not
formalised here (see `docs/FORMAL_ASSURANCE_BOUNDARY.md`).

## 2. Certificate structure

### 2.1 Trie form (`FTree`)

```
inductive FTree
  | exact                          -- leaf: a complete conformation
  | bound (transfer : Bool)        -- leaf: a whole subtree, bounded
          (witnesses : List (List (Fin 4)))
  | dead                           -- leaf: a prefix with no valid completion
  | node (t0 t1 t2 t3 : FTree)     -- branch on the next direction
```

The position of a trie node determines its prefix word `p`: the root is
`[]`, and child `d` of `p` is `p ++ [d]`. As a result, **duplicate frontier
entries, overlapping ancestor/descendant entries, hidden descendants of
covered nodes, and ambiguous encodings cannot be represented at all.** Each
prefix has exactly one position in the trie.

### 2.2 Flat form (`FlatCert`, Python-compatible)

This is a list of `(word, kind)` entries. `buildTrie` turns it into an
`FTree` and **rejects**:

* two entries with the same word (duplicate);
* an entry together with an entry for one of its proper extensions
  (ancestor/descendant overlap or hidden descendant);
* a missing child whose prefix is still self-avoiding (omitted branch).
  A missing child whose prefix already collides is filled in as `dead`; its
  deadness is then re-checked by `checkNode` like any other `dead` leaf.

`buildTrie` itself needs no proof. Whatever trie it produces, `check` still
validates it independently.

### 2.3 Text form (unverified parser, `FoldBounds/FrontierCLI.lean`)

```
FOLDBOUNDS-FRONTIER 1
sequence HHPHH...           (H/P only, length >= 2)
q 3/2                       (natural or natural/positive natural)
threshold 2                 (natural)
E 0203                      exact leaf at word 0203
B 02 21 23                  basic bound leaf at 02 with witness suffixes 21, 23
T 3 .                       transfer bound leaf at 3 with witness suffix "" (only if complete)
D 001                       dead leaf
end
```

Words are strings over `0123`; `.` denotes the empty word. Any other token,
an extra space, a missing header, or a missing `end` makes the parser reject
the file. **The parser is an unverified boundary.** It has no theorems. The
runner prints the parsed sequence, `q`, and threshold back out so that they
can be compared with what was intended. All soundness claims are relative to
the *parsed* parameters.

## 3. What the checker verifies

`check S q thr t` first requires `2 ≤ |S|` and `1 ≤ q`. It then runs
`checkNode (|S|-2) S q thr [] t`, which visits each trie node once:

| Node at prefix `p` | Local test (all decided in Lean) | Contribution |
|---|---|---|
| `exact` | `|p| = steps` and `walkFrom p` has no repeated site | the exact weight `q^k(p)`, counted as event or complement |
| `dead` | `walkFrom p` collides, **or** `|p| < steps` and all 4 one-step extensions collide | 0 |
| `bound tr ws` | `|p| ≤ steps`, `walkFrom p` self-avoiding, the witness list `ws` has no duplicates, and every `p ++ w` is a complete self-avoiding word | an `Enclosure` (see §4) |
| `node t0..t3` | `|p| < steps`, `walkFrom p` self-avoiding, and all four children accepted | the sum of the children's enclosures |

So the checker rejects invalid geometry (self-intersection), wrong anchoring
(which cannot occur, because the anchor is built into `walkFrom`), exact
leaves that are not terminal or not self-avoiding, dead leaves that have a
valid extension one step on (they fail the local test), descent below the
chain length, duplicate witnesses, witnesses that are incomplete, too long,
or self-intersecting, `q < 1`, and sequences shorter than 2. Witnesses are
always re-scored by the checker; a certificate cannot state a contact count.

**Every certificate quantity that the conclusion depends on is recomputed or
checked by the checker.** The untrusted producer chooses only *where* to
stop refining and which witnesses to propose.

## 4. Bound leaves: the enclosure of a whole subtree

Let `m = |p| + 2` be the number of placed sites, `known = knownContacts S p`,
and `r = steps - |p|`.

* **Basic mode (`tr = false`, text kind `B`).** Upper bound
  `3^r · q^Kmax`, where `Kmax = known + futurePairs S m` (the number of
  parity-compatible H–H pairs whose larger index is not yet placed). The
  event part is `0` if `thr > Kmax`; the complement part is `0` if
  `known ≥ thr`. This is the same rule as Python's `verify_certificate`, and
  the differential test checks that the numbers agree exactly.
* **Transfer mode (`tr = true`, text kind `T`).** Each upper bound is the
  minimum of
  1. the degree-capped counting envelope `3^r · q^(known + futureCap S m)`, and
  2. a sum over the merged dynamic-programming table `relaxedTable` of relaxed
     walks. These walks never reverse and never step on a fixed prefix site,
     and each table entry records its exact cross-contact score `c` with the
     prefix. Its weight is `q^min(known + c + ffCap, known + futureCap)`.
     The event and complement versions add the obvious indicator tests.
* **Witnesses.** Each verified witness `w` adds its exact weight to the lower
  bound. It also replaces its own envelope term in the upper bound by the
  exact weight: the checker subtracts `envelope(w) − exact(w) ≥ 0`. Witnesses
  are distinct (`Nodup`) and all lie in this leaf's subtree, so the leaf's
  witness sums never count a conformation twice. Witnesses in different
  leaves lie in disjoint subtrees (`leaf_cover_unique`).

* **Certified diagonal family (transfer mode only, no certificate data).**
  The checker computes `famOK p` from the prefix alone. It holds when the
  last prefix site maximises one of the four diagonal potentials `±x ± y`
  over all prefix sites. In that case all `2^r` continuations that use only
  the two potential-increasing steps are distinct valid completions
  (`nodup_family`, `card_family`, `family_subset_ext`). Each has at least
  `known` contacts, so `famLow = 2^r · q^known ≤ Z(p)` (`family_lower`,
  `famLow_le`). The lower bound `zLo` becomes `max(witness sum, famLow)`.
  The event lower bound is raised only if `known ≥ thr`, since then every
  family member is in the event. The complement lower bound is raised only
  if `kcap < thr`, since then every member is in the complement. Witnesses
  and the family may overlap, so the two are combined by **max, never by
  addition** (`famRaise_sound`).

Event and complement upper bounds are computed **separately**. A subtree
that may contain both kinds of conformation contributes to both upper
bounds, so `eHi + nHi ≥ zHi` is allowed. Nothing in the proofs assumes that
these two upper bounds add up to the total.

## 5. Probability endpoints and degenerate cases

`finalize` returns `pLo = eLo / (eLo + nHi)` and `pHi = eHi / (eHi + nLo)`,
using Lean's convention `x / 0 = 0`. The soundness proof uses the
coupled-ratio theorem `FoldBounds.ratio_mem_coupled` together with
`massZ_root_pos`: the straight rod is always a conformation, so `Z > 0` for
every admissible model. Degenerate outputs are therefore still sound:

* `eLo = nHi = 0` gives `pLo = 0`, a trivial lower bound;
* `eHi = nLo = 0` gives `pHi = 0`. This can only happen when the true event
  mass is 0, since `eHi ≥ Z_A ≥ 0`. In that case `P(A) = 0` and the enclosure
  `[0, 0]` is correct.

## 6. What acceptance implies

```
theorem FrontierCheck.check_sound (steps) (seq : Fin (steps+2) → Bool) (q : ℚ) (thr : ℕ)
    (t : FTree) (r : FrontierResult) (h : check (List.ofFn seq) q thr t = some r) :
    Z ∈ [r.zLo, r.zHi] ∧ Z_A ∈ [r.eLo, r.eHi] ∧ Z_¬A ∈ [r.nLo, r.nHi] ∧
    Z_A / Z ∈ [r.pLo, r.pHi]
```

Here `Z`, `Z_A`, `Z_¬A` are the `Executable.massQ` sums over the Lean HP
conformation type defined above. The theorem has **no other hypotheses**:
`q ≥ 1` and `|S| ≥ 2` are checked by `check` itself, and nonemptiness of the
conformation space is proved, not assumed. `check_coupled_probability_sound`
restates the probability claim over `ℝ` in exactly the form of the earlier
`checkHP_coupled_sound`. `checkFlat_sound` gives the same conclusion for flat
certificates.

Further guarantees:

* `leaf_cover_unique`: for an accepted certificate, every valid complete
  conformation lies in the completion set of **exactly one** frontier leaf.
  Hence the frontier partitions the conformation space, and dead leaves
  cover the empty set.
* `sum_ext_eq_sum_leaves`: any additive quantity over the conformation space
  is the sum over the leaves' completion sets.
* `check_root_bound_accepted`: the one-node certificate `bound tr []` is
  accepted for every admissible model. So acceptance never requires
  enumerating the conformation space.
* `transfer_total_le_basic` and `basicZ_degree_le_basic`: the transfer-mode
  total upper bound of a leaf is never larger than the basic-mode bound.

## 7. Cost of checking (analysis, not a kernel-checked complexity theorem)

`checkNode` is structurally recursive and visits each trie node exactly once.
The work at each node is:

* `exact` and `dead`: `O(n²)` for the self-avoidance tests and the contact
  count.
* `bound` in basic mode: `O(n²)` (contact caps are `O(n²)` pair sums).
* `bound` in transfer mode: a dynamic program with `r` layers. In layer `i`
  the states are (site, last direction, cross score, index). The site lies in
  a box of side `2i+1` around the prefix end, so there are at most
  `4·(2i+1)²·(C+1)` distinct states after merging, where `C ≤ 2n` bounds the
  cross score. Each state does `O(n)` work for its cross contacts against the
  prefix, and merging sorts the layer. This gives a polynomial cost per leaf,
  roughly `O(r³·n·C·log)`, which is **independent of the number of
  completions** (that number grows like `2.64^r`).
* witnesses: `O(n²)` each.

The total is therefore polynomial in the certificate size and `n`. It is not
proportional to the number of conformations. This analysis is **not** a
formal Lean complexity theorem. The formal compactness statement is
`check_root_bound_accepted`, together with the fact that `checkNode` is
defined by structural recursion on the submitted trie.

## 8. Explicit trust boundary

Trusted:

* the Lean kernel and the pinned toolchain;
* the definitions of the model (`walkFrom`, `Conformation`, `contactCount`,
  `weightQ`, `contactEvent`, `massQ`);
* the statement of `check_sound` / `checkFlat_sound`;
* for executed checks only, the Lean compiler/interpreter. The proofs are
  about the Lean function, and running compiled code trusts code generation.
  The kernel-checked examples in `FoldBounds/FrontierTests.lean` use
  `decide +kernel`, which does not involve the compiler.

Untrusted, with every claim it makes re-checked: the Python producer,
frontier choice, witness search, and certificate serialiser.

Unverified boundary: the text parser and JSON printer (`FrontierCLI`) and the
command-line runner `scripts/foldcheck.lean`.
