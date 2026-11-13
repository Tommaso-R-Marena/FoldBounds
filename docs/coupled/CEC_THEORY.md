# Coupled Ensemble Certificates (CEC): mathematical specification

**Status (2026-11-12): working specification, written by the coordinating
investigator.** The proofs below are written by hand. A result is reported as
*established* only after it is (i) machine-checked in Lean 4.28 / Mathlib with
an audited axiom footprint, and (ii) independently reviewed. Until then, every
statement here is a *candidate theorem with a proof sketch*. Novelty is **not**
claimed for any identity in this file; see the prior-art review
(`docs/coupled/PRIOR_ART.md`, pending) before quoting any originality.

## Formalization status (2026-11-13)

Spec §1, §2.1–2.2 and §3 are machine-checked in Lean 4.28 / Mathlib. There
are 93 audited declarations, using only the standard axioms (`propext`,
`Classical.choice`, `Quot.sound`), with no `sorry` or `native_decide`. They
live on branch `research/cec-lean-20261112`. The theorem-by-theorem map,
including exact hypotheses, is in `docs/coupled/LEAN_FORMALIZATION.md`.

The formalization corrected two statements below (Theorems 2.4 and 4.2) and
made several hypotheses precise:
- Theorems 1.1 and 1.2 do not need Adm ≠ ∅.
- `λ₊ ≤ S₊` needs ΣL > 0.
- Theorem 1.6 needs some U_k > 0.
- Theorem 1.5 needs S₋ > 0. Counterexample: regions (L,U,r) = (0,1,1) and
  (1,1,0).
- Theorem 2.2 needs only a nonempty support T⁺.

Not formalized: T2.7, T2.8, §2.3, Theorem 4.4, T1′ completeness, and the
sorted-threshold form of Theorem 1.3. (Theorem 1.3 is formalized as
attainment plus g₋(λ₋) = 0.)

## 0. The question

For two related ensembles A and B (wild type and mutant; ligand 1 and
ligand 2; protonated and deprotonated site) we want rigorous bounds on

    Q = Z_B / Z_A,      ΔF = -β⁻¹ log Q,

that remain informative when the individual partition functions Z_A and
Z_B cannot be bounded usefully, and that can be checked without enumerating
the conformational space.

The organising observation is elementary, and that is deliberate. Separate
bounds `Z_A ∈ [L_A, U_A]` and `Z_B ∈ [L_B, U_B]` throw away the fact that A
and B share almost all of their structure. Section 1 says exactly how much a
shared decomposition recovers. Section 2 says exactly what can and cannot be
certified from local information, with matching (sharp) lower bounds.

Conventions. All weights are nonnegative reals. Executable versions use exact
rationals. A "region" is any index of a decomposition shared by A and B.
`osc(f) = max f − min f`.

---

## 1. Coupled linear-fractional enclosures (T1)

### 1.1 Data and admissible set

A finite index set `K` of regions. For each `k ∈ K`:

* mass bounds `0 ≤ L_k ≤ U_k` for the A-mass `α_k` of region k;
* ratio bounds `0 ≤ r_k ≤ R_k` for the regional ratio, i.e.
  `r_k α_k ≤ β_k ≤ R_k α_k`, where `β_k` is the B-mass of region k.

Admissible set:

    Adm = { (α, β) ∈ ℝ^K_{≥0} × ℝ^K_{≥0} :
            L_k ≤ α_k ≤ U_k,  r_k α_k ≤ β_k ≤ R_k α_k  (all k),  Σ_k α_k > 0 }.

The target quantity is `Q(α, β) = Σ_k β_k / Σ_k α_k`. Note `Σ α > 0` and
`β_k ≤ R_k α_k` encode a support assumption: a region with no A-mass has no
B-mass. When this fails for a region, its ratio is unbounded above and the
upper enclosure must be reported as +∞.

### 1.2 The test functions

    g₋(λ) = Σ_k min((r_k − λ) L_k, (r_k − λ) U_k)
          = Σ_k (r_k − λ) · (L_k if r_k ≥ λ else U_k)

    g₊(λ) = Σ_k max((R_k − λ) L_k, (R_k − λ) U_k)
          = Σ_k (R_k − λ) · (U_k if R_k ≥ λ else L_k)

### 1.3 Theorems

**Theorem 1.1 (exact lower test).** Assume `Adm ≠ ∅`. For every real λ,

    (∀ (α,β) ∈ Adm, λ ≤ Q(α,β))   ⟺   g₋(λ) ≥ 0.

*Proof.* (⇐) For admissible (α,β):
`Σβ − λΣα ≥ Σ (r_k − λ) α_k ≥ g₋(λ) ≥ 0`. The middle step holds because
`(r_k − λ) α_k` is linear in `α_k ∈ [L_k, U_k]`, so it is at least the smaller
endpoint value. Divide by `Σα > 0`.
(⇒) Suppose `g₋(λ) < 0`. Set `α_k = L_k` if `r_k ≥ λ`, else `U_k`, and set
`β_k = r_k α_k`. Then `Σ (r_k − λ) α_k = g₋(λ) < 0`. If `Σα` were 0, every
term would vanish, so `Σα > 0`. Thus (α,β) is admissible and `Q < λ`. ∎

**Theorem 1.2 (exact upper test).** Assume `Adm ≠ ∅`. For every real λ,
`(∀ adm, Q ≤ λ) ⟺ g₊(λ) ≤ 0`. Symmetric proof, with `β_k = R_k α_k` and
`α_k = U_k` if `R_k ≥ λ`, else `L_k`.

**Theorem 1.3 (sharp enclosure, threshold form).** `λ₋ = min_Adm Q` and
`λ₊ = max_Adm Q` exist and are attained with `α_k ∈ {L_k, U_k}`. Sort the
regions by `r_k` ascending. Then `λ₋` is the minimum, over the `|K|+1`
thresholds t, of the ratio obtained by giving `U` to the t smallest-r regions
and `L` to the others. Symmetrically, `λ₊` gives `U` to the largest-R
regions. At the optimum `g₋(λ₋) = 0` and `g₊(λ₊) = 0`. A certificate for a
claimed `λ₋` is therefore the single exact check `g₋(λ₋) ≥ 0`, which is linear
time. No sorting is needed to *check*; only the proposer sorts.
*Proof sketch.* Dinkelbach optimality. At a minimiser α*, each coordinate
minimises `(r_k − λ₋) α_k` over `[L_k, U_k]`. Otherwise moving it would make
`Σ(r − λ₋)α < 0`, contradicting minimality by Theorem 1.1. ∎

**Theorem 1.4 (dominance over separate enclosures).** Let
`S₋ = Σ r_k L_k / Σ U_k` and `S₊ = Σ R_k U_k / Σ L_k`. These are the
endpoints obtained by bounding `Z_B` and `Z_A` separately from the same data.
Then

    λ₋ ≥ max(S₋, min{ r_k : U_k > 0 }),    λ₊ ≤ min(S₊, max{ R_k : U_k > 0 }).

**Theorem 1.5 (exactly when coupling helps).** Suppose `S₋ > 0`. Then
`λ₋ = S₋` if and only if `r_k (U_k − L_k) = 0` for every k. So the coupled
lower endpoint is *strictly* better than the separate one unless every region
with uncertain mass has zero lower ratio. The gain factor `λ₋ / S₋` is
unbounded: with one region, `λ₋ = r₁` while `S₋ = r₁ L₁ / U₁`. Symmetrically,
when `ΣL > 0`, `λ₊ = S₊` iff `R_k (U_k − L_k) = 0` for every k.
*Proof.* For admissible α, `Σ r α / Σ α ≥ Σ r L / Σ α ≥ Σ r L / Σ U`. The
first step is an equality iff `r_k (α_k − L_k) = 0` for all k. The second is
an equality iff `α = U` (because `S₋ > 0` forces `Σ r L > 0`). Both hold for
some admissible α iff `r_k (U_k − L_k) = 0` for all k. ∎

**Theorem 1.6 (uniform-ratio exactness).** If `r_k = R_k = ρ` for every k
with `U_k > 0`, then `λ₋ = λ₊ = ρ`, whatever the mass bounds are.

**Theorem 1.7 (monotone refinement).** (a) Shrinking any interval
`[L_k,U_k]` or `[r_k,R_k]` can only raise λ₋ and lower λ₊. (b) Split region
k into k₁, k₂ with consistent data: `L_{k1}+L_{k2} ≥ L_k`,
`U_{k1}+U_{k2} ≤ U_k`, and `[r_{kj}, R_{kj}] ⊆ [r_k, R_k]`. The refined
enclosure is then contained in the coarse one. *Proof.* Every refined
admissible configuration maps, by summing the split regions, to a coarse
admissible configuration with the same Q. ∎

### 1.4 Polytope generalisation (T1′) and LP dual certificates

Let the A-masses obey general linear constraints
`P = { α ≥ 0 : G α ≤ h }`, for example boxes, total-mass slabs, or
ratio-of-mass constraints `α_j ≤ c α_k` obtained from injections. Ratio bounds
are as before. **Certificate for `λ ≤ Q` on all admissible configurations:**
a vector `y ≥ 0` with

    (r − λ·1) + Gᵀ y ≥ 0   (componentwise)     and     hᵀ y ≤ 0.

*Soundness.* For `α ∈ P`:
`Σ(r_k − λ) α_k ≥ −(Gᵀy)ᵀα = −yᵀ(Gα) ≥ −yᵀh ≥ 0`. ∎
Completeness holds by LP duality whenever P is nonempty and bounded; only
soundness is needed for checking. The upper bound is symmetric.

---

## 2. Separator (Markov-blanket) certificates (T2, T3)

### 2.1 Abstract factorised form

There are finite sets: boundary states `T`, interior states `I_A` and `I_B`,
and exterior states `O`. All functions are nonnegative:

* `a : I_A × T → ℝ≥0`, the interior weight of system A, including every
  factor that touches the interior;
* `b : I_B × T → ℝ≥0`, the same for system B;
* `e : T × O → ℝ≥0`, the exterior weight, **shared** by A and B.

    Z_A = Σ_{i,τ,o} a(i,τ) e(τ,o),       Z_B = Σ_{i,τ,o} b(i,τ) e(τ,o).

This is exactly the situation in a factor graph whose variables are split as
interior ∪ boundary ∪ exterior when (i) no factor touches both the interior
and the exterior, and (ii) every factor that differs between A and B lies
inside interior ∪ boundary. Interior domains may differ between A and B, as
when a mutation changes the rotamer set.

Define `Φ_X(τ) = Σ_i x(i,τ)` for X ∈ {A, B} and `Ψ(τ) = Σ_o e(τ,o)`.

**Lemma 2.0 (factorisation).** `Z_X = Σ_τ Φ_X(τ) Ψ(τ)`.

Support: let `T⁺ = { τ : Φ_A(τ) > 0 }`. If some τ has `Φ_A(τ) = 0 < Φ_B(τ)`,
the sharp upper bound below is +∞. States with `Φ_A = Φ_B = 0` are inert. On
T⁺ define the **local ratio** `ρ(τ) = Φ_B(τ) / Φ_A(τ)`.

### 2.2 Theorems

**Theorem 2.1 (local enclosure).** Suppose `Z_A > 0` and there is no τ with
`Φ_A(τ) = 0 < Φ_B(τ)`. Then

    min_{τ ∈ T⁺} ρ(τ)  ≤  Z_B / Z_A  ≤  max_{τ ∈ T⁺} ρ(τ).

*Proof.* `Z_B = Σ_{τ∈T⁺} ρ(τ) Φ_A(τ) Ψ(τ)`, which is a convex combination of
the values ρ(τ) with weights `Φ_A Ψ / Z_A`. ∎

**Theorem 2.2 (sharpness, instance optimality).** Fix the local data
(a, b), and assume T⁺ = T (all Φ_A > 0). For every λ with
`min ρ ≤ λ ≤ max ρ`, there is an exterior with Z_A > 0 and `Z_B/Z_A = λ`. One
such exterior has a single state and is supported on at most two boundary
states. For every λ in the open interval, a *strictly positive* exterior
achieves λ.
*Proof.* Let τ₋ and τ₊ be the argmin and argmax, and write
`λ = θ ρ(τ₋) + (1−θ) ρ(τ₊)`. Set `Ψ(τ₋) = θ/Φ_A(τ₋)`,
`Ψ(τ₊) = (1−θ)/Φ_A(τ₊)`, and Ψ = 0 elsewhere. Then `Z_A = 1` and `Z_B = λ`.
For a positive exterior, mix in ε mass on all states and re-solve for θ; this
works because the open interval is the image of the open simplex. ∎
*Realisability by pairwise exteriors.* One hidden variable h with domain T,
hard pairwise couplings `[h_v = σ_v]` to each boundary variable v, and a unary
factor ψ(h) produce `Ψ = ψ`. So sharpness holds even among exteriors built
from pairwise factors. Hard 0/1 couplings can be replaced by soft couplings in
the limit.

**Corollary 2.2′ (instance optimality).** Consider any procedure whose output
depends only on the local data (I_A, I_B, T, a, b) and which is sound for
every exterior. Its output interval must contain `[min ρ, max ρ]`. The
separator enclosure is therefore the tightest sound local certificate.

**Theorem 2.3 (separation: absolute vs relative).** Fix the local data. For
every `z > 0` and every `λ ∈ [min ρ, max ρ]`, some exterior gives
`Z_A = z` and `Z_B / Z_A = λ`. Consequently:
(i) no sound local procedure can bound Z_A, and hence F_A, within any finite
interval; (ii) the ratio is confined to `[min ρ, max ρ]`. Separate local
bounds yield only the vacuous ratio interval `[0, +∞]`. The coupled bound is
finite.
*Proof.* Scale the exterior of Theorem 2.2 by z. ∎

**Theorem 2.4 (nested hierarchy, monotone).** *[Corrected 2026-11-13 after
Lean formalization. The inner support condition is required; see
`CoupledSeparator.hierarchy_needs_support`. Counterexample: Φ_A = (1, 0),
Φ_B = (1, 1), K ≡ 1, which gives ρ = 1 on T⁺ but ρ′ = 2.]* Assume the inner
support condition `Φ_A(τ) = 0 ⇒ Φ_B(τ) = 0` for all τ; it then propagates to
the outer separator (`push_supp`). Suppose a larger separator satisfies, for
X ∈ {A, B},

    Φ'_X(τ') = Σ_τ Φ_X(τ) K(τ, τ')

with a **shared** nonnegative kernel K. This holds when the old boundary and
an annulus of variables are absorbed into the interior, and the annulus
factors are identical in A and B. Then for every τ' with `Φ'_A(τ') > 0`,

    ρ'(τ') = Σ_τ κ(τ'; τ) ρ(τ),    κ(τ'; τ) = Φ_A(τ) K(τ,τ') / Φ'_A(τ').

Here κ(τ'; ·) is a probability vector, so `[min ρ', max ρ'] ⊆ [min ρ, max ρ]`.

**Theorem 2.5 (finite termination).** If the separator grows to cover the
whole system (T' a singleton and the exterior trivial), the enclosure
collapses to the exact value `Z_B/Z_A`.

**Theorem 2.6 (a-posteriori contraction).**
`osc(ρ') ≤ δ(κ) · osc(ρ)`, where
`δ(κ) = max_{τ'₁,τ'₂} ½ Σ_τ |κ(τ'₁;τ) − κ(τ'₂;τ)|` is the Dobrushin
ergodicity coefficient of the annulus kernel. The proof uses
`|p·ρ − p'·ρ| ≤ TV(p,p')·osc(ρ)` for probability vectors p, p'. Products of
layer coefficients bound the convergence of the hierarchy. This is a
diagnostic and is not needed for soundness.

### 2.3 Partially known exteriors

* *Interval exterior.* Suppose `Ψ(τ) ∈ [ℓ(τ), u(τ)]`, for example from bounds
  on the exterior's conditional partition functions. Apply T1 with regions τ,
  `α_τ ∈ [Φ_A ℓ, Φ_A u]`, and exact ratio ρ(τ).
* *Bounded exterior coupling (soft separator).* Suppose
  `Ψ(τ)/Ψ(τ') ≤ c(τ,τ')`, for example
  `c = exp(Σ_{v: τ_v ≠ τ'_v} s_v)` when the exterior couples to boundary
  variable v with log-sensitivity at most s_v. The sharp enclosure is then a
  linear-fractional program over a polyhedral cone. Apply T1′ with LP dual
  certificates.

### 2.4 Soft separators: weakly coupled exteriors (T2.7)

Setting: boundary variables v ∈ ∂ with finite domains D_v, and known local
functions `Φ_A, Φ_B : T = Π_v D_v → ℝ_{>0}`. **Exterior class E_s.** The
exterior interacts with the boundary only through factors that each touch at
most ONE boundary variable; factors among exterior variables are arbitrary.
Equivalently, `Ψ(τ) = Σ_o g(o) Π_v f_v(τ_v, o)` with `g ≥ 0` and
`f_v(·,o) > 0`. For every exterior state o,
`max_a f_v(a,o) / min_a f_v(a,o) ≤ e^{s_v}`. For pairwise couplings `ψ_{vu}`
to exterior variables u, one valid choice is
`s_v = Σ_u max_b [log max_a ψ_{vu}(a,b) − log min_a ψ_{vu}(a,b)]`, which is
computable from the factors touching v.

For a field `F = (F_v)_v` with `F_v : D_v → ℝ_{>0}`, let
`Z_X(F) = Σ_τ Φ_X(τ) Π_v F_v(τ_v)` and `R(F) = Z_B(F)/Z_A(F)`. A **vertex
field** has `F_v(a) ∈ {1, e^{s_v}}` for all v and a.

**Theorem 2.7 (soft-separator enclosure, sharp).** For every exterior in E_s,

    min_{F vertex} R(F)  ≤  Z_B/Z_A  ≤  max_{F vertex} R(F),

and both endpoints are attained by exteriors in E_s: a single exterior state
whose couplings realise the optimal vertex field.

*Proof.* (1) `Q = Σ_o w_o R(f(·,o))` with convex weights
`w_o ∝ g(o) Z_A(f(·,o))`, so Q lies between the extreme R(f(·,o)). (2) R is
invariant under rescaling each F_v, so normalise each `f_v(·,o)` to
`[1, e^{s_v}]`; this gives a point of the box `B = Π_{v,a} [1, e^{s_v}]`.
(3) In each single coordinate F_v(a), with the others fixed, `Z_X = c_X +
d_X F_v(a)` with c_X, d_X ≥ 0. So R is linear-fractional in that coordinate,
hence monotone, and its extremes over the interval are at endpoints. Moving
coordinates one at a time gives a vertex at least as extreme as any box
point. (4) Sharpness: the optimal vertex field is a member of E_s. ∎

Remarks. With s = 0 the exterior is irrelevant and Q = Z_B(1)/Z_A(1) is exact.
As s_v → ∞ for all v, the enclosure tends to `[min ρ, max ρ]` (Theorem 2.1).
For binary boundary variables only the effective fields ±s_v matter, giving
2^{|∂|} candidates. The class E_s matters: if exterior factors may touch
several boundary variables, Ψ need only be log-Lipschitz. That cone strictly
contains the mixtures of bounded products, and the sharp bound becomes the
LP of §2.3. (Numerical check by the coordinator: 0 violations in 90,000
random mixture-exterior trials, file /tmp/claude-0/pi/cex_soft.py.)

### 2.5 When do local certificates converge? (T2.8)

Let `V_s = s·1{σ_v = a}` be a local field at site v, and let
`W_R(s) = osc_τ ρ_R(τ)` be the hard-separator width at radius R. Since
`ρ_R(τ) = E^τ_{B_R}[e^{−V_s}] = 1 − (1 − e^{−s}) μ^τ_{B_R}(σ_v = a)`,

    W_R(s) = (1 − e^{−s}) · osc_τ μ^τ_{B_R}(σ_v = a).

**Theorem 2.8 (local certifiability ⇔ weak spatial mixing).** The local
enclosures for every field perturbation at v shrink to zero as R → ∞
iff the marginal at v becomes asymptotically independent of the boundary
condition (weak spatial mixing at v). By Theorem 2.2 this is also a statement
about *every* sound local method. The exact identity above is elementary; the
content is the interpretation.

Consequences, as hypotheses to test numerically:
* 2D Ising, zero field: below T_c the width plateaus at about
  (1 − e^{−s})·m* (phase coexistence). Above T_c it decays exponentially.
  At T_c it decays as a power law.
* Fixed-backbone side-chain or protonation models conditional on the folded
  state: expect fast decay. Global two-state folding is a finite-size analogue
  of phase coexistence. This is why mutational *stability* ΔΔG (folded vs
  unfolded) is not a local quantity, while the folded-state term is.

**Structural fact (proton binding is antiferromagnetic).** With protonation
variables σ ∈ {0,1}, `q_acid = σ − 1`, `q_base = σ`, and Coulomb
`W_ij > 0`, the σ_iσ_j coefficient of `W_ij q_i q_j` is `+W_ij` for every
site-type pair. Multi-site titration is therefore an antiferromagnetic 2-spin
system with arbitrary fields. Consequences: (i) self-avoiding-walk-tree
(Weitz) recursions are monotone with alternating direction, so depth-R
truncation with all-protonated or all-deprotonated free leaves gives rigorous
two-evaluation bounds; (ii) worst-case hardness of approximating
antiferromagnetic 2-spin partition functions beyond uniqueness (Sly–Sun)
applies to the model *class*. (Prior-art status: pending, PRIOR_ART.md.)

---

## 3. Binding and selectivity (T4)

The model is discrete: receptor variables, a ligand with a finite set of
poses or conformations, and an explicit standard-state convention. The
convention is the same pose grid for both ligands, so the standard-state
constants cancel in ratios. **Locality assumption (A-loc):** the ligand
interacts only with the pocket variables, whose joint state is τ ∈ T.

Write `Ψ_R(τ)` for the apo receptor weight with the pocket in state τ, summed
over everything else; then `Z_R = Σ_τ Ψ_R(τ)`. Write `Φ_L(τ)` for the
ligand's summed weight in the pocket in state τ, and `Z_L` for its free
weight. Then

    K_L = Z_RL / (Z_R Z_L) = Σ_τ p_apo(τ) k_L(τ),   with
    k_L(τ) = Φ_L(τ) / Z_L   and   p_apo(τ) = Ψ_R(τ)/Z_R.

This is the conformational-selection identity, which is classical.

**Theorem 4.1 (certified affinity).**
`K_L ∈ [min_τ k_L(τ), max_τ k_L(τ)]`, and this is sharp over receptor
ensembles. The absolute binding constant is itself a ratio, so it *is* locally
certifiable. The receptor's own free energy is not (Theorem 2.3).

**Theorem 4.2 (certified selectivity).** *[Corrected 2026-11-13: this
requires `k₁(τ) > 0` for every pocket state τ. Counterexample
`CoupledBinding.selectivity_needs_positive_k1`: k₁ = (1, 0), k₂ = (1, 1) with
a uniform receptor gives S = 2, although k₂/k₁ ≤ 1 wherever k₁ > 0.]*
`S = K_{L2}/K_{L1} ∈ [min_τ k₂/k₁, max_τ k₂/k₁]`, and this is sharp. If
`min_τ k₂(τ)/k₁(τ) > 1`, then L2 binds more tightly than L1 for **every**
receptor ensemble.

**Theorem 4.3 (non-identifiability).** If
`min_τ k₂/k₁ < 1 < max_τ k₂/k₁`, there exist receptor ensembles that make L2
the tighter binder and others that make L1 the tighter binder. Local
structural information alone then cannot certify the selectivity ordering.

**Theorem 4.4 (partial ensemble knowledge).** Suppose experimental or
computed bounds on apo pocket populations are available:
`p_apo ∈ Π = { p ≥ 0, Σp = 1, A p ≤ c }`. The sharp selectivity interval is
then a linear-fractional program over Π, certified by T1′ duals.

Interpretation (a hypothesis to test, not a theorem). Long-range allosteric
coupling of the pocket is exactly what widens `[min ρ, max ρ]` across larger
separators. Allostery is the obstruction to local certification of
selectivity.

---

## 4. Instantiations

* **2D HP lattice (FoldBounds native).** There is no fixed separator, because
  polymers lack a spatial Markov blanket in sequence space. Use T1 with
  regions given by prefix-frontier subtrees. Each region's mass bounds come
  from the existing frontier machinery. Its ratio bounds come from local
  contact accounting at the mutated residue i: for an H→P mutation,
  `k_B(x) = k_A(x) − c_i(x)`, so `ρ_T ∈ [q^{−max c_i}, q^{−min c_i}]` over
  completions in T. Expected outcome: modest gains, because subtree ratios
  vary strongly and caps are loose. This is a *testbed with exact oracles*,
  not the flagship.
* **Factor graphs / MRFs (rotamer, Potts, and protonation models).** Exact
  separators exist. Prediction to test: the local enclosure width does not
  depend on system size, while separate-bound widths grow extensively.
* **Multi-site proton binding (pKa).** Sites are binary protonation variables
  with pairwise electrostatic couplings, a cutoff, and model intrinsic pKa
  values. The perturbation "site i protonated vs deprotonated" is local, so
  `P(σ_i=1)/P(σ_i=0)` is a separator ratio. This gives certified *model* pKa
  intervals for proteins too large for exact inference. It is **not** a claim
  of experimental pKa accuracy.

## 5. What these results do not establish

* Nothing about real molecules beyond the specified discrete model.
* No efficiency guarantee in general. Computing `min_τ ρ(τ)` exactly can
  require enumerating boundary states, which is exponential in separator
  size. Certified relaxations (separable bounds and branch-and-bound) are
  sound but may be loose.
* No claim that the identities (FEP, the Markov property, conformational
  selection, linear-fractional duality) are new. The candidate contributions
  are the sharpness, separation, and instance-optimality statements; their
  machine-checked form; executable certificates; and the measured
  size-independence.
