module

public import FoldBounds.FrontierDegree

/-!
# Verified square-lattice geometric contact bounds

For a self-avoiding square-lattice chain `W` of `N ≥ 1` sites, let
`rows` (resp. `cols`) be the number of distinct `y` (resp. `x`)
coordinates it occupies.  We prove the edge-counting (perimeter)
inequality

  `#{(i, j) : i + 1 < j, W i ~ W j} ≤ N + 1 - rows - cols`,

i.e. the number of nonconsecutive lattice contacts is at most
`N + 1 - rows - cols` (`contacts_le_rows_cols`).  Because a chain is
connected, every row between two occupied sites is occupied
(`rows_ge_span`), and `N ≤ rows * cols`.  Hence, for any list `K` of
sites known to lie on the chain (e.g. the fixed prefix, the endpoint, or
retained history), the number of contacts is at most

  `geomCap N K = N + 1 - minPerim N (xSpan K) (ySpan K)`

(`parityContacts_le_geomCap`), where `minPerim N w h` is the least
`c + r` over integer rectangles `c × r` with `c ≥ w`, `r ≥ h`, `r c ≥ N`.

For an all-H chain of 18 residues this gives at most `10` contacts
(versus `19` for the degree cap) at the root, and it is exactly the
classical bound `K ≤ n + 1 - w - h` for the bounding rectangle `w × h`.
-/
@[expose] public section

namespace FoldBounds
namespace Geometry

open LatticeModel Frontier Finset

/-! ### Fibre-maximum lemma -/

/-- If `s` preserves `g` and strictly increases `f`, then in every
`g`-fibre the `f`-maximal site has no `s`-image in `T`.  Hence the sites
whose `s`-image lies in `T`, plus the number of fibres, is at most `|T|`. -/
theorem card_shift_add_image_le (T : Finset Site) (f g : Site → ℤ) (s : Site → Site)
    (hg : ∀ a, g (s a) = g a) (hf : ∀ a, f a < f (s a)) :
    (T.filter (fun a => s a ∈ T)).card + (T.image g).card ≤ T.card := by
  have hsplit := Finset.card_filter_add_card_filter_not (s := T) (fun a => s a ∈ T)
  have himg : T.image g ⊆ (T.filter (fun a => s a ∉ T)).image g := by
    intro y hy
    obtain ⟨a0, ha0, rfl⟩ := Finset.mem_image.1 hy
    have hne : (T.filter (fun b => g b = g a0)).Nonempty := ⟨a0, by simp [ha0]⟩
    obtain ⟨a, ha, hmax⟩ := Finset.exists_max_image _ f hne
    rw [Finset.mem_filter] at ha
    refine Finset.mem_image.2 ⟨a, Finset.mem_filter.2 ⟨ha.1, fun hsa => ?_⟩, ha.2⟩
    have := hmax (s a) (Finset.mem_filter.2 ⟨hsa, by rw [hg, ha.2]⟩)
    exact absurd (hf a) (not_lt.2 this)
  have := (Finset.card_le_card himg).trans Finset.card_image_le
  omega

/-! ### Index pairs and sites -/

/-- Shift of a site. -/
def sh (dx dy : ℤ) (a : Site) : Site := (a.1 + dx, a.2 + dy)

section Chain

variable (W : List Site)

/-- All ordered index pairs below `|W|`. -/
def idxPairs : Finset (ℕ × ℕ) := range W.length ×ˢ range W.length

/-- Ordered index pairs whose sites differ by the shift `(dx, dy)`. -/
def dPairs (dx dy : ℤ) : Finset (ℕ × ℕ) :=
  (idxPairs W).filter (fun x => siteAt W x.2 = sh dx dy (siteAt W x.1))

/-- Ordered lattice-adjacent index pairs. -/
def adjPairs : Finset (ℕ × ℕ) :=
  (idxPairs W).filter (fun x => manhattan1 (siteAt W x.1) (siteAt W x.2))

/-- Nonconsecutive lattice contacts (one orientation). -/
def contactPairs : Finset (ℕ × ℕ) :=
  (idxPairs W).filter (fun x => x.1 + 1 < x.2 ∧ manhattan1 (siteAt W x.1) (siteAt W x.2))

/-- Occupied rows. -/
def rows : Finset ℤ := W.toFinset.image Prod.snd

/-- Occupied columns. -/
def cols : Finset ℤ := W.toFinset.image Prod.fst

variable {W}

theorem siteAt_mem {i : ℕ} (hi : i < W.length) : siteAt W i ∈ W := by
  rw [siteAt_eq_getElem hi]; exact List.getElem_mem hi

theorem siteAt_inj (hn : W.Nodup) {i j : ℕ} (hi : i < W.length) (hj : j < W.length)
    (h : siteAt W i = siteAt W j) : i = j := by
  rw [siteAt_eq_getElem hi, siteAt_eq_getElem hj] at h
  exact (List.Nodup.getElem_inj_iff hn).1 h

theorem card_dPairs_le (hn : W.Nodup) (dx dy : ℤ) :
    (dPairs W dx dy).card ≤ (W.toFinset.filter (fun a => sh dx dy a ∈ W.toFinset)).card := by
  apply Finset.card_le_card_of_injOn (fun x => siteAt W x.1)
  · intro x hx
    simp only [dPairs, idxPairs, Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_product,
      Finset.mem_range] at hx
    simp only [Finset.coe_filter, Set.mem_setOf_eq, List.mem_toFinset]
    exact ⟨siteAt_mem hx.1.1, hx.2 ▸ siteAt_mem hx.1.2⟩
  · intro x hx y hy hxy
    simp only [dPairs, idxPairs, Finset.coe_filter, Set.mem_setOf_eq, Finset.mem_product,
      Finset.mem_range] at hx hy
    have h1 : x.1 = y.1 := siteAt_inj hn hx.1.1 hy.1.1 hxy
    have h2 : x.2 = y.2 := siteAt_inj hn hx.1.2 hy.1.2 (by rw [hx.2, hy.2]; simp only at hxy; rw [hxy])
    exact Prod.ext h1 h2

theorem adjPairs_subset :
    adjPairs W ⊆ dPairs W 1 0 ∪ dPairs W (-1) 0 ∪ dPairs W 0 1 ∪ dPairs W 0 (-1) := by
  intro x hx
  simp only [adjPairs, dPairs, Finset.mem_filter, Finset.mem_union] at hx ⊢
  obtain ⟨hp, hm⟩ := hx
  generalize siteAt W x.1 = a at hm ⊢
  generalize siteAt W x.2 = b at hm ⊢
  obtain ⟨a1, a2⟩ := a
  obtain ⟨b1, b2⟩ := b
  simp only [manhattan1, decide_eq_true_eq, sh, Prod.mk.injEq, hp, true_and] at hm ⊢
  omega

theorem card_adjPairs_le (hn : W.Nodup) :
    (adjPairs W).card + 2 * (rows W).card + 2 * (cols W).card ≤ 4 * W.length := by
  have hT : W.toFinset.card = W.length := List.toFinset_card_of_nodup hn
  have e1 := card_shift_add_image_le W.toFinset Prod.fst Prod.snd (sh 1 0)
    (fun a => by simp [sh]) (fun a => by simp [sh])
  have e2 := card_shift_add_image_le W.toFinset (fun a => -a.1) Prod.snd (sh (-1) 0)
    (fun a => by simp [sh]) (fun a => by simp [sh])
  have e3 := card_shift_add_image_le W.toFinset Prod.snd Prod.fst (sh 0 1)
    (fun a => by simp [sh]) (fun a => by simp [sh])
  have e4 := card_shift_add_image_le W.toFinset (fun a => -a.2) Prod.fst (sh 0 (-1))
    (fun a => by simp [sh]) (fun a => by simp [sh])
  have d1 := card_dPairs_le hn 1 0
  have d2 := card_dPairs_le hn (-1) 0
  have d3 := card_dPairs_le hn 0 1
  have d4 := card_dPairs_le hn 0 (-1)
  have hu0 := Finset.card_le_card (adjPairs_subset (W := W))
  have hu1 := Finset.card_union_le (dPairs W 1 0 ∪ dPairs W (-1) 0 ∪ dPairs W 0 1) (dPairs W 0 (-1))
  have hu2 := Finset.card_union_le (dPairs W 1 0 ∪ dPairs W (-1) 0) (dPairs W 0 1)
  have hu3 := Finset.card_union_le (dPairs W 1 0) (dPairs W (-1) 0)
  simp only [rows, cols] at *
  omega

/-- Backbone pairs `(i, i+1)`. -/
def fwdPairs : Finset (ℕ × ℕ) := (range (W.length - 1)).image (fun i => (i, i + 1))

theorem card_fwdPairs : (fwdPairs (W := W)).card = W.length - 1 := by
  rw [fwdPairs, Finset.card_image_of_injective _ (fun a b h => by simpa using h)]
  simp

theorem adjacent_siteAt (hc : adjacentChain W = true) {i : ℕ} (hi : i + 1 < W.length) :
    manhattan1 (siteAt W i) (siteAt W (i + 1)) = true := by
  rw [siteAt_eq_getElem (by omega), siteAt_eq_getElem hi, ← adjacent_eq_manhattan1]
  exact adjacentChain_getElem hc hi

theorem manhattan1_symm {a b : Site} (h : manhattan1 a b = true) : manhattan1 b a = true := by
  rw [← adjacent_eq_manhattan1] at h ⊢; exact adjacent_symm h

theorem card_contact_lower (hc : adjacentChain W = true) :
    2 * (W.length - 1) + 2 * (contactPairs W).card ≤ (adjPairs W).card := by
  set B1 := fwdPairs (W := W)
  set B2 := B1.image Prod.swap
  set C := contactPairs W
  set Cs := C.image Prod.swap
  have hB2 : B2.card = B1.card := Finset.card_image_of_injective _ Prod.swap_injective
  have hCs : Cs.card = C.card := Finset.card_image_of_injective _ Prod.swap_injective
  have mB1 : ∀ x ∈ B1, x.2 = x.1 + 1 ∧ x.1 + 1 < W.length := by
    intro x hx
    obtain ⟨i, hi, rfl⟩ := Finset.mem_image.1 hx
    simp at hi ⊢; omega
  have mB2 : ∀ x ∈ B2, x.1 = x.2 + 1 ∧ x.2 + 1 < W.length := by
    intro x hx
    obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hx
    have := mB1 y hy; simp; omega
  have mC : ∀ x ∈ C, x.1 + 1 < x.2 ∧ x.2 < W.length ∧
      manhattan1 (siteAt W x.1) (siteAt W x.2) = true := by
    intro x hx
    simp only [C, contactPairs, idxPairs, Finset.mem_filter, Finset.mem_product,
      Finset.mem_range] at hx
    exact ⟨hx.2.1, hx.1.2, hx.2.2⟩
  have mCs : ∀ x ∈ Cs, x.2 + 1 < x.1 ∧ x.1 < W.length ∧
      manhattan1 (siteAt W x.1) (siteAt W x.2) = true := by
    intro x hx
    obtain ⟨y, hy, rfl⟩ := Finset.mem_image.1 hx
    have := mC y hy
    exact ⟨this.1, this.2.1, manhattan1_symm this.2.2⟩
  have hsub : B1 ∪ B2 ∪ C ∪ Cs ⊆ adjPairs W := by
    intro x hx
    simp only [Finset.mem_union] at hx
    simp only [adjPairs, idxPairs, Finset.mem_filter, Finset.mem_product, Finset.mem_range]
    rcases hx with ((h | h) | h) | h
    · obtain ⟨h1, h2⟩ := mB1 x h
      refine ⟨⟨by omega, by omega⟩, ?_⟩
      rw [h1]; exact adjacent_siteAt hc h2
    · obtain ⟨h1, h2⟩ := mB2 x h
      refine ⟨⟨by omega, by omega⟩, ?_⟩
      rw [h1]; exact manhattan1_symm (adjacent_siteAt hc h2)
    · obtain ⟨h1, h2, h3⟩ := mC x h
      exact ⟨⟨by omega, h2⟩, h3⟩
    · obtain ⟨h1, h2, h3⟩ := mCs x h
      exact ⟨⟨h2, by omega⟩, h3⟩
  have d12 : Disjoint B1 B2 := Finset.disjoint_left.2 (fun x h1 h2 => by
    have := mB1 x h1; have := mB2 x h2; omega)
  have d3 : Disjoint (B1 ∪ B2) C := Finset.disjoint_left.2 (fun x h1 h2 => by
    have := mC x h2
    rcases Finset.mem_union.1 h1 with h | h
    · have := mB1 x h; omega
    · have := mB2 x h; omega)
  have d4 : Disjoint (B1 ∪ B2 ∪ C) Cs := Finset.disjoint_left.2 (fun x h1 h2 => by
    have := mCs x h2
    rcases Finset.mem_union.1 h1 with h | h
    · rcases Finset.mem_union.1 h with h | h
      · have := mB1 x h; omega
      · have := mB2 x h; omega
    · have := mC x h; omega)
  have := Finset.card_le_card hsub
  rw [Finset.card_union_of_disjoint d4, Finset.card_union_of_disjoint d3,
    Finset.card_union_of_disjoint d12, hB2, hCs, card_fwdPairs] at this
  omega

/-- **Perimeter inequality.** A self-avoiding lattice chain of `N ≥ 1`
sites occupying `rows` rows and `cols` columns has at most
`N + 1 - rows - cols` nonconsecutive lattice contacts. -/
theorem contacts_le_rows_cols (hn : W.Nodup) (hc : adjacentChain W = true)
    (hpos : 0 < W.length) :
    (contactPairs W).card + (rows W).card + (cols W).card ≤ W.length + 1 := by
  have h1 := card_adjPairs_le hn
  have h2 := card_contact_lower hc
  omega

/-! ### Connectivity: occupied rows and columns form intervals -/

/-- Discrete intermediate value theorem for unit-step integer sequences. -/
theorem discrete_ivt (f : ℕ → ℤ) (N : ℕ) (hstep : ∀ k, k + 1 < N → |f (k + 1) - f k| ≤ 1)
    (y : ℤ) : ∀ d i, i + d < N → (f i ≤ y ∧ y ≤ f (i + d) ∨ f (i + d) ≤ y ∧ y ≤ f i) →
      ∃ k < N, f k = y := by
  intro d
  induction d with
  | zero =>
      intro i hi h
      exact ⟨i, by omega, by simp only [add_zero] at h; omega⟩
  | succ d ih =>
      intro i hi h
      have hs := hstep (i + d) (by omega)
      rw [abs_le] at hs
      by_cases hy : f (i + (d + 1)) = y
      · exact ⟨i + (d + 1), hi, hy⟩
      · apply ih i (by omega)
        have : i + (d + 1) = i + d + 1 := by ring
        rw [this] at h hy
        omega

theorem coord_step (hc : adjacentChain W = true) (π : Site → ℤ)
    (hπ : ∀ a b : Site, manhattan1 a b = true → |π b - π a| ≤ 1) :
    ∀ k, k + 1 < W.length → |π (siteAt W (k + 1)) - π (siteAt W k)| ≤ 1 :=
  fun _ hk => hπ _ _ (adjacent_siteAt hc hk)

theorem exists_index_of_mem {a : Site} (ha : a ∈ W) : ∃ i < W.length, siteAt W i = a := by
  obtain ⟨i, hi, rfl⟩ := List.getElem_of_mem ha
  exact ⟨i, hi, siteAt_eq_getElem hi⟩

/-- Interval property of a coordinate along a connected chain. -/
theorem coord_interval (hc : adjacentChain W = true) (π : Site → ℤ)
    (hπ : ∀ a b : Site, manhattan1 a b = true → |π b - π a| ≤ 1)
    {a b : Site} (ha : a ∈ W) (hb : b ∈ W) :
    Finset.Icc (π a) (π b) ⊆ W.toFinset.image π := by
  intro y hy
  rw [Finset.mem_Icc] at hy
  obtain ⟨i, hi, rfl⟩ := exists_index_of_mem ha
  obtain ⟨j, hj, rfl⟩ := exists_index_of_mem hb
  have key : ∃ k < W.length, π (siteAt W k) = y := by
    rcases le_total i j with hij | hij
    · obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hij
      exact discrete_ivt (fun k => π (siteAt W k)) W.length (coord_step hc π hπ) y d i hj
        (Or.inl hy)
    · obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hij
      exact discrete_ivt (fun k => π (siteAt W k)) W.length (coord_step hc π hπ) y d j hi
        (Or.inr hy)
  obtain ⟨k, hk, rfl⟩ := key
  exact Finset.mem_image.2 ⟨siteAt W k, List.mem_toFinset.2 (siteAt_mem hk), rfl⟩

theorem rows_ge_span (hc : adjacentChain W = true) {a b : Site} (ha : a ∈ W) (hb : b ∈ W) :
    (b.2 - a.2 + 1).toNat ≤ (rows W).card := by
  have := Finset.card_le_card (coord_interval hc Prod.snd (fun a b h => by
    obtain ⟨a1, a2⟩ := a; obtain ⟨b1, b2⟩ := b
    simp only [manhattan1, decide_eq_true_eq] at h; simp only; rw [abs_le]; omega) ha hb)
  rw [Int.card_Icc] at this
  have e : b.2 - a.2 + 1 = b.2 + 1 - a.2 := by ring
  rw [e]; exact this

theorem cols_ge_span (hc : adjacentChain W = true) {a b : Site} (ha : a ∈ W) (hb : b ∈ W) :
    (b.1 - a.1 + 1).toNat ≤ (cols W).card := by
  have := Finset.card_le_card (coord_interval hc Prod.fst (fun a b h => by
    obtain ⟨a1, a2⟩ := a; obtain ⟨b1, b2⟩ := b
    simp only [manhattan1, decide_eq_true_eq] at h; simp only; rw [abs_le]; omega) ha hb)
  rw [Int.card_Icc] at this
  have e : b.1 - a.1 + 1 = b.1 + 1 - a.1 := by ring
  rw [e]; exact this

theorem length_le_rows_mul_cols (hn : W.Nodup) :
    W.length ≤ (rows W).card * (cols W).card := by
  have hT : W.toFinset.card = W.length := List.toFinset_card_of_nodup hn
  have hsub : W.toFinset ⊆ (cols W) ×ˢ (rows W) := by
    intro a ha
    exact Finset.mem_product.2 ⟨Finset.mem_image.2 ⟨a, ha, rfl⟩,
      Finset.mem_image.2 ⟨a, ha, rfl⟩⟩
  have := Finset.card_le_card hsub
  rw [Finset.card_product, hT] at this
  linarith [mul_comm (rows W).card (cols W).card]

theorem cols_le_length : (cols W).card ≤ W.length :=
  Finset.card_image_le.trans (List.toFinset_card_le W)

end Chain

/-! ### Executable spans and the geometric cap -/

/-- Maximum of a list of naturals (0 for the empty list). -/
def lmax (l : List ℕ) : ℕ := l.foldr max 0

theorem le_lmax {l : List ℕ} {x : ℕ} (h : x ∈ l) : x ≤ lmax l := by
  induction l with
  | nil => simp at h
  | cons a l ih =>
      simp only [lmax, List.foldr_cons] at ih ⊢
      rcases List.mem_cons.1 h with rfl | h
      · exact le_max_left _ _
      · exact (ih h).trans (le_max_right _ _)

theorem lmax_le {l : List ℕ} {X : ℕ} (h : ∀ x ∈ l, x ≤ X) : lmax l ≤ X := by
  induction l with
  | nil => simp [lmax]
  | cons a l ih =>
      simp only [lmax, List.foldr_cons] at ih ⊢
      exact max_le (h a (by simp)) (ih (fun x hx => h x (by simp [hx])))

/-- Width (number of columns) spanned by a list of known sites. -/
def xSpan (K : List Site) : ℕ := lmax (K.flatMap (fun a => K.map (fun b => (b.1 - a.1 + 1).toNat)))

/-- Height (number of rows) spanned by a list of known sites. -/
def ySpan (K : List Site) : ℕ := lmax (K.flatMap (fun a => K.map (fun b => (b.2 - a.2 + 1).toNat)))

/-- Least `c + r` over rectangles with `c ≥ w` columns, `r ≥ h` rows and
`r c ≥ N` cells (columns range over `w ≤ c ≤ N`). -/
def minPerim (N w h : ℕ) : ℕ :=
  ((List.range (N + 1)).filter (fun c => w ≤ c)).foldr
    (fun c acc => min acc (c + max h ((N + c - 1) / c))) (2 * N + 2)

theorem minPerim_le (N w h c r : ℕ) (hw : w ≤ c) (hcN : c ≤ N) (hh : h ≤ r) (hrc : N ≤ r * c)
    (hc : 0 < c) : minPerim N w h ≤ c + r := by
  have hmem : c ∈ (List.range (N + 1)).filter (fun c => w ≤ c) := by
    simp [hw]; omega
  have hdiv : (N + c - 1) / c ≤ r := by
    rw [Nat.div_le_iff_le_mul_add_pred hc]
    have : N ≤ c * r := by rw [mul_comm]; exact hrc
    omega
  have key : ∀ l : List ℕ, c ∈ l → l.foldr (fun c acc => min acc (c + max h ((N + c - 1) / c)))
      (2 * N + 2) ≤ c + max h ((N + c - 1) / c) := by
    intro l hl
    induction l with
    | nil => simp at hl
    | cons a l ih =>
        simp only [List.foldr_cons]
        rcases List.mem_cons.1 hl with rfl | hl
        · exact min_le_right _ _
        · exact (min_le_left _ _).trans (ih hl)
  exact (key _ hmem).trans (by omega)

/-- **Geometric contact cap** for a chain of `N` sites that is known to
contain all sites of `K`. -/
def geomCap (N : ℕ) (K : List Site) : ℕ := N + 1 - minPerim N (xSpan K) (ySpan K)

theorem xSpan_le {W K : List Site} (hc : adjacentChain W = true) (hK : ∀ a ∈ K, a ∈ W) :
    xSpan K ≤ (cols W).card :=
  lmax_le (fun x hx => by
    simp only [List.mem_flatMap, List.mem_map] at hx
    obtain ⟨a, ha, b, hb, rfl⟩ := hx
    exact cols_ge_span hc (hK a ha) (hK b hb))

theorem ySpan_le {W K : List Site} (hc : adjacentChain W = true) (hK : ∀ a ∈ K, a ∈ W) :
    ySpan K ≤ (rows W).card :=
  lmax_le (fun x hx => by
    simp only [List.mem_flatMap, List.mem_map] at hx
    obtain ⟨a, ha, b, hb, rfl⟩ := hx
    exact rows_ge_span hc (hK a ha) (hK b hb))

/-- **Geometric cap theorem (chains).** -/
theorem contactPairs_le_geomCap {W K : List Site} (hn : W.Nodup) (hc : adjacentChain W = true)
    (hpos : 0 < W.length) (hK : ∀ a ∈ K, a ∈ W) :
    (contactPairs W).card ≤ geomCap W.length K := by
  have h1 := contacts_le_rows_cols hn hc hpos
  have hcols : 0 < (cols W).card := by
    obtain ⟨a, ha⟩ := List.exists_mem_of_length_pos hpos
    exact Finset.card_pos.2 ⟨a.1, Finset.mem_image.2 ⟨a, List.mem_toFinset.2 ha, rfl⟩⟩
  have hp := minPerim_le W.length (xSpan K) (ySpan K) (cols W).card (rows W).card
    (xSpan_le hc hK) cols_le_length (ySpan_le hc hK) (length_le_rows_mul_cols hn) hcols
  unfold geomCap
  omega

/-- The Python-style H–H contact count is at most the number of
nonconsecutive lattice contacts. -/
theorem parityContactCount_le_contactPairs (S : List Bool) (W : List Site)
    (hS : S.length = W.length) : parityContactCount S W ≤ (contactPairs W).card := by
  rw [parityContactCount_eq_pairSum]
  unfold pairSum
  rw [← Finset.sum_product' (f := fun i j => if (candidatePair S i j && decide (j < W.length) &&
      manhattan1 (siteAt W i) (siteAt W j)) = true then 1 else 0)]
  rw [Finset.sum_boole]
  simp only [Nat.cast_id]
  apply Finset.card_le_card
  intro x hx
  simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_range, Bool.and_eq_true,
    decide_eq_true_eq] at hx
  simp only [contactPairs, idxPairs, Finset.mem_filter, Finset.mem_product, Finset.mem_range]
  have := candidatePair_lt hx.2.1.1
  exact ⟨⟨by omega, by omega⟩, by omega, hx.2.2⟩

/-- **Geometric cap for complete conformations.** Every valid complete
word `l` whose walk contains all sites of `K` has at most
`geomCap (steps + 2) K` H–H contacts. -/
theorem parityContacts_le_geomCap {steps : ℕ} (S : List Bool) (hS : S.length = steps + 2)
    {l : List (Fin 4)} (hl : l ∈ validWords steps) {K : List Site}
    (hK : ∀ a ∈ K, a ∈ walkFrom l) :
    parityContactCount S (walkFrom l) ≤ geomCap (steps + 2) K := by
  rw [mem_validWords] at hl
  have hlen : (walkFrom l).length = steps + 2 := by rw [walkFrom_length, hl.1]
  have := contactPairs_le_geomCap hl.2 (walkFrom_adjacentChain l) (by omega) hK
  rw [hlen] at this
  exact (parityContactCount_le_contactPairs S _ (by omega)).trans this

end Geometry
end FoldBounds
