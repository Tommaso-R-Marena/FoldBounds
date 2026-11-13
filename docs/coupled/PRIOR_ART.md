# Prior-art review: Coupled Ensemble Certificates (CEC)

**Status (2026-11-13).** Targeted literature review (Agent F) for
`docs/coupled/CEC_THEORY.md`. Citations marked **[verified]** were confirmed
against Crossref, arXiv, or a publisher or author page during this review.
Anything marked **[unverified]** is from memory or secondary reading and must
be checked before it is quoted. The review is targeted, not exhaustive. Not
searched: HP-lattice mutation literature, Coq/Isabelle libraries, any
non-English or pre-1990 statistical-physics sources beyond those named.
Ratings: **high** = essentially known, **medium** = known ingredients with a
possibly new combination or statement, **low** = no prior statement found.

## 0. Summary

Every mathematical ingredient in CEC is classical: the Zwanzig identity,
the DLR/conditional-specification identity, the mediant inequality,
linear-fractional programming, Dobrushin contraction, and the
spatial-mixing characterisation. The two strongest overlaps are:

1. **Leisink and Kappen's Bound Propagation** already bounds a marginal by
   an LP over the unknown distribution of the Markov blanket. This is T2.1
   and T2.2 for marginals.
2. **The repository's own PR #91 documents** already contain
   `Z_mut/Z_parent = E_parent[R]` and a coupled-mass polytope.

The defensible novelty is narrow (Section 6). It is not any identity or
the basic enclosure.

## 1. Free-energy perturbation family (underlies T1, T2.1, HP mutation ratios)

**Established.**
- Zwanzig, R. W. (1954). *High-Temperature Equation of State by a
  Perturbation Method. I. Nonpolar Gases.* J. Chem. Phys. 22, 1420.
  doi:10.1063/1.1740409 [verified]. Gives `Z_B/Z_A = <e^{-βΔV}>_A`.
- Bennett, C. H. (1976). *Efficient estimation of free energy differences
  from Monte Carlo data.* J. Comput. Phys. 22, 245.
  doi:10.1016/0021-9991(76)90078-4 [verified].
- Shirts, M. R. & Chodera, J. D. (2008). *Statistically optimal analysis of
  samples from multiple equilibrium states.* J. Chem. Phys. 129, 124105.
  doi:10.1063/1.2978177 [verified].
- Gibbs–Bogoliubov inequality: `F_B ≤ F_A + <ΔV>_A` (Bogoliubov 1958,
  Feynman 1955) [unverified citations]. Directly implied: since `ΔV` lies in
  `[min ΔV, max ΔV]`, `|ΔF| ≤ ‖ΔV‖∞` is immediate from Zwanzig. Textbook
  material.
- Reweighting on exact enumerations is already used internally: PR #91,
  `COUNTERFACTUAL_FOLDING_REWEIGHTING.md`, uses `Z_mut/Z_parent = E_parent[R]`
  and conditional-odds identities on fully enumerated tiny cubic HP models.

**Distinct.** Nothing in the identities. Applying them conditional on a
boundary state, with an explicit support condition, is a routine
specialisation.

**Risk: HIGH** for the FEP identity, and for the HP mutation-ratio
application at the identity level.

**Do NOT say:** "we introduce/derive a reweighting identity for mutation
ratios"; "novel free-energy bound"; "beyond FEP". Say instead: "we apply
Zwanzig's identity conditionally and bound it."

## 2. Group T1 and T1′: coupled linear-fractional enclosures

**Established.**
- The mediant inequality (`Σβ/Σα` lies between the smallest and largest
  `β_k/α_k`) is classical, often credited to Cauchy [unverified specific
  citation]. It is T2.1 and the trivial case of T1.
- Charnes, A. & Cooper, W. W. (1962). *Programming with linear fractional
  functionals.* Nav. Res. Logist. Q. 9, 181–186. doi:10.1002/nav.3800090303
  [verified]. The transformation turns T1 and T1′ into LPs. Hence
  strong-duality certificates (T1′) are standard LP theory.
- Dinkelbach, W. (1967). *On nonlinear fractional programming.* Management
  Sci. 13(7), 492–498 [unverified; Crossref lookup failed]. The `g(λ) ≥ 0`
  test (T1.1–1.3) is the Dinkelbach parametric test. Schaible's later
  analyses of Dinkelbach's method (e.g. Management Sci. 22 (1976)) [unverified].
  Box constraints with sorting by `r_k` make the extreme-point structure and
  the closed form a routine consequence.
- Interval / monotone-in-coordinate enclosures (T1.4–T1.7) are standard
  interval-arithmetic and monotonicity facts.

**What practitioners do (key question).**
- K*: Lilien, R. H., Stevens, B. W., Anderson, A. C. & Donald, B. R. (2005).
  *A novel ensemble-based scoring and search algorithm for protein redesign…*
  J. Comput. Biol. 12, 740. doi:10.1089/cmb.2005.12.740 [verified].
- Georgiev, I., Lilien, R. H. & Donald, B. R. (2008). *The minimized
  dead-end elimination criterion and its application to protein redesign in
  a hybrid scoring and search algorithm.* J. Comput. Chem. 29, 1527.
  doi:10.1002/jcc.20909 [verified].
- BBK*: Ojewole, A. A., Jou, J. D., Fowler, V. G., Donald, B. R. (2018).
  J. Comput. Biol. doi:10.1089/cmb.2017.0267 [verified].
- MARK*: Jou, J. D., Holt, G. T., Lowegard, A. U., Donald, B. R. (2020).
  *Minimization-aware recursive K*.* J. Comput. Biol.
  doi:10.1089/cmb.2019.0315 [verified]. Conference version: RECOMB 2019,
  LNCS, doi:10.1007/978-3-030-17083-7_7 [verified].
- FRIES/EWAK*: Lowegard, A. U., Frenkel, M. S., Holt, G. T., Jou, J. D.,
  Ojewole, A. A., Donald, B. R. (2020). *Novel, provable algorithms for
  efficient ensemble-based computational protein design and their
  application to the redesign of the c-Raf-RBD:KRas protein–protein
  interface.* PLoS Comput. Biol. 16, e1007447.
  doi:10.1371/journal.pcbi.1007447 [verified].
- COMETS: Hallen, M. A. & Donald, B. R. (2016). J. Comput. Biol.
  doi:10.1089/cmb.2015.0188 [verified].
- BWM*: Jou, J. D., Jain, S., Georgiev, I. S., Donald, B. R. (2016, JCB
  doi:10.1089/cmb.2015.0194; RECOMB 2015, LNCS doi:10.1007/978-3-319-16706-0_16)
  [verified]. The brief's "2018" date is wrong. The paper is 2015/2016.
- LUTE: Hallen, M. A., Jou, J. D., Donald, B. R. (2017). J. Comput. Biol.
  doi:10.1089/cmb.2016.0136 [verified]. Related, not a ratio bound.
- Toulbar2 weighted counting: Viricel, C., Simoncini, D., Barbe, S.,
  Schiex, T. (2016). *Guaranteed Weighted Counting for Affinity Computation:
  Beyond Determinism and Structure.* CP 2016, LNCS.
  doi:10.1007/978-3-319-44953-1_46 [verified]. Earlier: Viricel, Simoncini,
  Allouche, de Givry et al. (2015), *Approximate Counting with Deterministic
  Guarantees for Affinity Computation*, doi:10.1007/978-3-319-18167-7_15
  [verified].

**Answer to the key question.** From the PLoS Comput. Biol. text I read
(the first ~100k characters), K* is defined as a ratio of partition functions,
and each state's (P, L, PL) partition function is approximated separately
to within a user-specified ε. No shared decomposition of the ratio is
described there. How per-state bounds are combined into K* bounds is
deferred to the BBK* paper, which I did not read: **[unverified]**. The
expected combination is the ratio of the separate endpoints. That is
exactly the "separate" baseline `S₋, S₊` of T1.4. Toulbar2 weighted counting
bounds each Z. I found no OSPREY or Toulbar2 paper bounding Z_B/Z_A through
a shared decomposition, but I did not read BBK* or MARK* in full and the
search was not exhaustive. A reader who finds one would collapse the T1
novelty for protein design. Check the full text of BBK* and MARK* first.

**Distinct.** Using the shared structure of A and B (coupled regions with
regional ratio bounds) rather than ratio-of-separate-bounds is not, as far
as I found, done in OSPREY or Toulbar2. But T1 is an elementary
consequence of linear-fractional programming. PR #91's
`COUPLED_MASS_POLYTOPE.md` already treats the two-region
(event/complement) case with a sharp joint interval. T1 generalises it and
shares its elementary status.

**Risk: HIGH** for the mathematics (T1.1–T1.3, T1′). **MEDIUM** for the
statement T1.4/T1.5 ("coupled strictly beats separate unless …,
unbounded gain"): the theorem is elementary but I found no published
statement. **MEDIUM** for use in protein design, pending the BBK*/MARK*
read.

**Do NOT say:** "new class of bounds"; "first coupled bound on partition
function ratios"; "beats OSPREY/K*" (no head-to-head exists); "novel
linear-fractional theorem". Say instead: "Dinkelbach-type exact tests applied
to coupled regional bounds."

## 3. Group T2.0–T2.3: separator enclosure, sharpness, instance optimality, absolute-vs-relative separation

**Established.**
- DLR/conditional-specification identity: a Gibbs measure's marginal is an
  average of conditional specifications over the exterior distribution.
  Georgii, H.-O. (1988; 2nd ed. 2011), *Gibbs Measures and Phase
  Transitions*, de Gruyter. doi:10.1515/9783110850147 [verified] and
  10.1515/9783110250329 [verified]. T2.0/T2.1 are exactly this (chapter
  references [unverified]). The ratio form is a convex-combination bound.
- **Leisink, M. A. R. & Kappen, H. J. (2003).** *Bound Propagation.* J.
  Artif. Intell. Res. 19, 139–154. doi:10.1613/jair.1130 [verified]. Per the
  abstract and paper pages I read: bounds on marginals of a cluster are
  computed by an LP whose objective is the cluster marginal and whose
  constraints describe the allowed distributions of the surrounding Markov
  blanket, iterated to propagate bounds between clusters. Earlier:
  Leisink & Kappen (2001), *A Tighter Bound for Graphical Models*, Neural
  Comput., doi:10.1162/089976601750399344 [verified, title only].
  **This is the most direct overlap.** "A marginal is a convex combination
  of Markov-blanket conditionals, with the bound sharp over an unknown
  exterior distribution on the blanket" is the Bound Propagation
  construction. Whether Leisink–Kappen *state* the sharpness-over-all-exteriors
  claim and the realisability by pairwise exteriors (T2.2) I could not
  verify from the pages read.
- Mooij, J. M. & Kappen, H. J. (2008). *Bounds on marginal probability
  distributions.* NIPS 21, 1105–1112 [verified via NeurIPS page].
  Long version: *Novel Bounds on Marginal Probabilities*, arXiv:0801.3797
  [verified]. Subtree box propagation and a self-avoiding-walk-tree bound,
  for arbitrary discrete factor graphs. Whether sharpness over exteriors is
  stated: **[unverified]**.
- Ihler, A. (2007). *Accuracy bounds for belief propagation.* UAI 2007,
  183–190; arXiv:1206.5277 [verified metadata; content only from secondary
  reading]. Dynamic-range message errors give marginal bounds.
- Wainwright, M. J., Jaakkola, T. S., Willsky, A. S. (2005). *A new class of
  upper bounds on the log partition function.* IEEE Trans. Inf. Theory
  51(7), 2313. doi:10.1109/tit.2005.850091 [verified]. Upper bounds on
  log Z only. Not ratios.
- Dechter, R. & Rish, I. (2003). *Mini-buckets: a general scheme for bounded
  inference.* J. ACM 50(2), 107. doi:10.1145/636865.636866 [verified].
  Bounds on each Z (and marginals) with no ratio coupling.
- Hard separators and nested hierarchies: Dobrushin–Shlosman's finite-volume
  criteria [see T2.8] and conditioning/cutset ideas in general.

**Distinct.**
(i) Stating the same-boundary-state *ratio* `Φ_B/Φ_A` as the certified
quantity (i.e. bounding a mutation/ligand ratio rather than a marginal) is a
small reformulation. Note that a marginal ratio `P(σ=1)/P(σ=0)` is already a
marginal-type quantity covered by Bound Propagation.
(ii) The explicit instance-optimality *corollary* 2.2′ (any sound local
procedure must contain `[min ρ, max ρ]`) and the realisation by
pairwise-factor exteriors: I did not find these stated. The proof of 2.2 is
simple (convex combination of two atoms), so the content is the framing,
not the mathematics.
(iii) T2.3 (no sound local procedure can bound Z_A; only the ratio is
certifiable): it is common knowledge that log Z is extensive and not locally
determined; I found no formal "local bounds cannot control Z but can control
ratios" theorem.

**Risk: HIGH** for T2.0/T2.1. **MEDIUM–HIGH** for T2.2 (sharpness over
all exteriors is natural for LP-based bounds; whether Leisink–Kappen state it
needs a full-text check). **MEDIUM** for T2.3 and Corollary 2.2′ as
formulated.

**Do NOT say:** "we discover that marginals are convex combinations of
conditionals"; "new local bound on marginals"; "first Markov-blanket
certificate". Do not claim that Bound Propagation, box propagation or
mini-buckets are "unsound/loose" without a measured comparison.

## 4. Group T2.4–T2.6: nested hierarchy, finite termination, Dobrushin contraction

**Established.**
- Monotone nested hierarchy of enclosures with larger regions: Bound
  Propagation iterates with growing clusters (Leisink–Kappen 2003).
  Mini-buckets with larger i-bound are monotone (Dechter–Rish 2003).
  Finite termination at the whole system is trivial (exact inference).
- Dobrushin ergodicity coefficient and contraction of total variation:
  Dobrushin, R. L. (1968). *The problem of uniqueness of a Gibbsian random
  field and the problem of phase transitions.* Funct. Anal. Appl. 2, 302.
  doi:10.1007/bf01075682 [verified]. Dobrushin & Shlosman (1985),
  *Constructive criterion for the uniqueness of Gibbs field*, in
  *Statistical Physics and Dynamical Systems*, Birkhäuser,
  doi:10.1007/978-1-4899-6653-7_20 [verified]. Föllmer, H. (1982). *A
  covariance estimate for Gibbs measures.* J. Funct. Anal. 46, 387.
  doi:10.1016/0022-1236(82)90053-2 [verified]. Georgii's book treats
  Dobrushin comparison and uniqueness [chapter refs unverified]. The
  inequality `|p·ρ − p'·ρ| ≤ TV(p,p')·osc(ρ)` is the standard TV bound.

**Distinct.** Only the specific bookkeeping (shared kernel `K`, annulus
absorption, `κ` weights) as applied to ratio oscillation. This is
bookkeeping on classical tools.

**Risk: HIGH** (T2.5 trivial; T2.6 is Dobrushin contraction; T2.4 is the
standard monotonicity of nested conditioning).

**Do NOT say:** "new convergence theory"; "novel contraction"; "we prove
convergence of the hierarchy" without citing Dobrushin. T2.6 must be cited
as an application of Dobrushin's coefficient.

## 5. T2.7: soft separators (weakly coupled exteriors)

**Established.**
- The vertex/box structure follows from coordinatewise monotonicity of a
  linear-fractional function (`R` is linear-fractional in each `F_v(a)`).
  This is the Charnes–Cooper/Dinkelbach structure again.
- Box propagation and dynamic-range bounds treat an unknown environment by
  bounded likelihood ratios: Mooij & Kappen 2008 (arXiv:0801.3797) and
  Ihler 2007 (UAI), above. The log-oscillation sensitivity `s_v` is the
  standard Dobrushin influence/oscillation quantity.
- The restricted class E_s ("each exterior factor touches ≤ 1 boundary
  variable") is a product-field mixture. Mixtures of product measures over
  an arbitrary exterior state are the general form for conditionally
  independent boundary couplings (a mean-field-like closure).

**Distinct.** The sharp vertex-field enclosure with attainment inside E_s,
and the observation that the broader log-Lipschitz cone gives a strictly
larger LP, are not stated in the sources I read. Sharpness is relative to a
purposely chosen exterior class. The coordinator's numerical check (90,000
trials, zero violations) is a test, not evidence of novelty.

**Risk: MEDIUM** (known components; the exact class-relative sharpness
statement not found).

**Do NOT say:** "sharp for weakly coupled exteriors in general" (sharp only
within E_s as defined); "rigorous bound for proteins" (the class E_s is a
modelling assumption); "a general solution to the environment problem".

## 6. T2.8: local certifiability ⇔ weak spatial mixing

**Established.**
- Weak spatial mixing and the identity `ρ_R(τ) = 1 − (1−e^{-s}) μ^τ(σ_v=a)`
  are elementary. WSM and its equivalence with decay of boundary influence
  go back to Dobrushin–Shlosman (1985) and the Martinelli–Olivieri
  programme [unverified specific citation]; see Georgii's book.
- Weitz, D. (2006). *Counting independent sets up to the tree threshold.*
  STOC '06, 140–149. doi:10.1145/1132516.1132538 [verified]. Reduces
  marginals to a self-avoiding-walk tree; correlation decay ⇒ algorithms.
- Bandyopadhyay, A. & Gamarnik, D. (2008). *Counting without sampling:
  asymptotics of the log-partition function for certain statistical physics
  models.* Random Struct. Alg. 33, 452. doi:10.1002/rsa.20236 [verified].
- Li, L., Lu, P., Yin, Y. (2013). *Correlation decay up to uniqueness in spin
  systems.* SODA 2013. doi:10.1137/1.9781611973105.5 [verified].
- Sinclair, A., Srivastava, P., Štefankovič, D., Yin, Y. (2017). *Spatial
  mixing and the connective constant: optimal bounds.* Probab. Theory
  Relat. Fields 168, 153. doi:10.1007/s00440-016-0708-2 [verified; volume
  and pages from memory, [unverified]].
- Two-spin 2D Ising phase behaviour (plateau of boundary influence below
  T_c, exponential decay above, power law at T_c) is textbook (Georgii).

**Distinct.** The only candidate is the *interpretive* link: "every sound
local certificate (T2.2) shrinks iff the boundary-dependence of the
marginal vanishes", i.e. restating WSM as a property of instance-optimal
certificates. The theorem's own text calls the identity elementary.

**Risk: HIGH.** The mathematical statement is the definition of WSM
restated through Theorem 2.2 (itself MEDIUM–HIGH).

**Do NOT say:** "we discover/prove a new characterisation of spatial
mixing"; "new connection between certification and phase transitions".
Say: "Corollary: instance-optimal local ratio certificates converge exactly
when WSM holds, by T2.2."

## 7. Group T4: binding, selectivity, non-identifiability

**Established.**
- Conformational selection with `K_obs = Σ p_i K_i` and population-weighted
  binding: Monod, J., Wyman, J., Changeux, J.-P. (1965). *On the nature of
  allosteric transitions: a plausible model.* J. Mol. Biol. 12, 88.
  doi:10.1016/s0022-2836(65)80285-6 [verified]. Boehr, D. D., Nussinov, R.,
  Wright, P. E. (2009). *The role of dynamic conformational ensembles in
  biomolecular recognition.* Nat. Chem. Biol. 5, 789. doi:10.1038/nchembio.232
  [verified]. Wyman, J. (1964). *Linked functions and reciprocal effects in
  hemoglobin: a second look.* Adv. Protein Chem. 19, 223.
  doi:10.1016/s0065-3233(08)60190-4 [verified]. Binding polynomials / linkage:
  Wyman & Gill, *Binding and Linkage* (1990) [unverified].
- Ratio of such averages with unknown weights `p` is again a
  linear-fractional program (T4.4 is Charnes–Cooper over a polytope).
- OSPREY/COMETS optimise affinity and specificity over sequence (above), but
  over given conformational ensembles and not as certified over unknown
  receptor ensembles.

**Distinct.** T4.1/T4.2 (certified affinity/selectivity over *all* receptor
ensembles given local pocket data) and T4.3 (explicit non-identifiability
when `min k₂/k₁ < 1 < max`) are the model-relative restatement of
T2.1/T2.2 for ligand ratios. I found no paper stating them in this form.
The conformational-selection identity is classical.

**Risk: HIGH** for the identity (K_obs = Σ p K) and T4.1. **MEDIUM** for
T4.2–T4.3 as framed. **MEDIUM/LOW** for T4.4 beyond Charnes–Cooper.

**Do NOT say:** "new model of selectivity"; "predicts selectivity";
"allostery is the obstruction" as a theorem (it is flagged a hypothesis in
the source; keep it so). Never imply experimental affinity or selectivity
accuracy.

## 8. Antiferromagnetic structure of proton binding

**Established.**
- Multi-site titration is a statistical-mechanical model over `2^N`
  protonation microstates with pairwise electrostatic coupling:
  Tanford, C. & Kirkwood, J. G. (1957). *Theory of protein titration
  curves. I.* J. Am. Chem. Soc. 79, 5333. doi:10.1021/ja01577a001 [verified].
  Bashford, D. & Karplus, M. (1990). *pKa's of ionizable groups in
  proteins: atomic detail from a continuum electrostatic model.*
  Biochemistry 29, 10219. doi:10.1021/bi00496a010 [verified]. Beroza, P.,
  Fredkin, D. R., Okamura, M. Y., Feher, G. (1991). *Protonation of
  interacting residues in a protein by a Monte Carlo method.* PNAS 88,
  5804. doi:10.1073/pnas.88.13.5804 [verified]. Gilson, M. K. (1993).
  *Multiple-site titration and molecular modeling: two rapid methods…*
  Proteins 15, 266. doi:10.1002/prot.340150305 [verified]. Georgescu, R. E.,
  Alexov, E. G., Gunner, M. R. (2002). *Combining conformational
  flexibility and continuum electrostatics for calculating pKas in
  proteins.* Biophys. J. 83, 1731. doi:10.1016/s0006-3495(02)73940-4
  [verified; MCCE]. Alexov & Gunner (1997), Biophys. J. 72, 2075.
  doi:10.1016/s0006-3495(97)78851-9 [verified]. Anandakrishnan, Aguilar,
  Onufriev (2012), *H++ 3.0*, Nucleic Acids Res. 40, W537.
  doi:10.1093/nar/gks375 [verified]. Olsson, M. H. M., Søndergaard, C. R.,
  Rostkowski, M., Jensen, J. H. (2011). *PROPKA3.* J. Chem. Theory Comput.
  7, 525. doi:10.1021/ct100578z [verified]. Pahari, S., Sun, L., Alexov, E.
  (2019). *PKAD.* Database, baz024. doi:10.1093/database/baz024 [verified].
  Onufriev & Alexov (2013), Q. Rev. Biophys. 46, 181.
  doi:10.1017/s0033583513000024 [verified].
- The sign structure (like-charge repulsion between sites penalises joint
  protonation of same-type sites) is a standard feature of these models. I
  did not find a source that names it "antiferromagnetic 2-spin system" in
  the Sly–Sun/Weitz sense. **[unverified]** that it has never been stated.
  Interacting-site titration as an Ising-type (lattice-gas) model is
  commonly said in the field. I could not retrieve a specific verified
  citation for the Ising phrase.
- Theory side: Weitz 2006; Li–Lu–Yin 2013; Sinclair, Srivastava, Thurley
  (2014), *Approximation algorithms for two-state anti-ferromagnetic spin
  systems on bounded degree graphs*, J. Stat. Phys. 155, 666.
  doi:10.1007/s10955-014-0947-5 [verified; vol./pages from memory]. Sly, A.
  & Sun, N. (2014). *Counting in two-spin models on d-regular graphs.* Ann.
  Probab. 42, 2383. doi:10.1214/13-aop888 [verified]; FOCS 2012 version
  arXiv:1203.2602 [verified]. Galanis, Štefankovič, Vigoda (2015).
  *Inapproximability for antiferromagnetic spin systems in the tree
  nonuniqueness region.* J. ACM 62(6). doi:10.1145/2785964 [verified].
  Mooij–Kappen's SAW-tree bound (arXiv:0801.3797) already gives rigorous
  marginal bounds for arbitrary factor graphs.
- Searches for rigorous error bounds on multi-site titration/pKa found only
  empirical or statistical error statements. For example, Myers, J.,
  Grothaus, G., Narayanan, S., Onufriev, A. (2006), *A simple clustering
  algorithm can be accurate enough for use in calculations of pKs in
  macromolecules*, Proteins 63, 928, doi:10.1002/prot.20922 [verified],
  reports a 95%-confidence error on a validation set, not a worst-case
  bound. This is a null result from a limited search, not proof of absence.

**Distinct.** Explicitly (a) casting the cutoff-electrostatics titration
model as an antiferromagnetic 2-spin system with arbitrary external fields,
and (b) using Weitz/Li–Lu–Yin truncation with monotone alternating leaf
bounds to give rigorous two-evaluation marginal bounds, as an application to
pKa. I found no paper doing this.

**Caution (does not weaken the above but limits claims).** Sly–Sun and
Galanis et al. hardness is for *families* of regular graphs of growing
size with uniform parameters beyond the uniqueness threshold. Real
titration instances are finite, dense (long-range Coulomb) and
non-uniform. Hardness does not transfer to a specific protein, and
uniqueness thresholds from bounded-degree theory do not apply to dense
graphs without a cutoff. Also, the Weitz tree has the exponential size of
the graph's connective constant, which is poor for dense couplings.

**Risk: HIGH** for "titration is an Ising-like model with repulsive
coupling" (field common knowledge). **MEDIUM** for the formal
antiferromagnetic 2-spin identification plus monotone SAW-tree bounds.
**LOW–MEDIUM** for *rigorous* certified intervals on titration (null
search; not exhaustive).

**Do NOT say:** "pKa is NP-hard / hard to approximate"; "titration
exhibits a computational phase transition"; "uniqueness regime of proteins";
"certified pKa" without "model" and "for the stated cutoff model"; "first to
note proton binding is antiferromagnetic".

## 9. Application summary

| Application | Established | Distinct | Risk |
|---|---|---|---|
| HP mutation ratios | Zwanzig; PR #91 `Z_mut/Z_parent=E_parent[R]` and event/complement polytope; Giaquinta–Pozzi, Lee, Hsieh et al. for exact HP/ISAW thermodynamics (see `PRIOR_ART_BASELINE_GATE.md`) | T1 with prefix-frontier regions and local contact-cap ratio bounds | HIGH for the identities; MEDIUM for the regional T1 enclosure (theory predicts modest gains) |
| MRF separator certificates | DLR/conditional specification; Leisink–Kappen 2003; Mooij–Kappen 2008; Ihler 2007; TRW; mini-bucket | Ratio-certified form; sharp over exteriors with realisability by pairwise exteriors; E_s vertex formula; measured size-independence | HIGH (T2.1), MEDIUM (T2.2/T2.7 packaging) |
| Multi-site titration / pKa | Tanford–Kirkwood; Bashford–Karplus; Beroza 1991; Gilson; MCCE; H++; PROPKA3; PKAD; Weitz/Li–Lu–Yin/Mooij–Kappen SAW bounds | Rigorous model-relative intervals from antiferromagnetic SAW truncation (none found) | LOW–MEDIUM (application); HIGH (structure fact) |
| Ligand selectivity | Conformational selection; MWC; Wyman linkage; Charnes–Cooper; COMETS/BBK* sequence-space affinity/specificity | Certified selectivity over unknown receptor ensembles; explicit non-identifiability | MEDIUM |

## 10. Formal verification

Is Dobrushin uniqueness or any Gibbs-measure theorem formalised?
- **Douglas, M. R., `MarkovSemigroups` (Lean 4).** Per its own page:
  a Dobrushin uniqueness module marked proved, with exponential correlation
  decay, zero sorries in that module and no axioms beyond Lean's built-ins.
  The project as a whole uses 8 postulated textbook axioms and has two
  sorries in a work-in-progress file. <https://reservoir.lean-lang.org/@mrdouglasny/MarkovSemigroups>
  [verified from page text; proofs not independently checked].
- **Douglas, `LGT` (Lean 4, lattice gauge theory).** Page reports a Gibbs
  specification, DLR identity, and Dobrushin condition verification
  <https://reservoir.lean-lang.org/@mrdouglasny/LGT> [search snippet only;
  unverified].
- **Dillies et al., `GibbsMeasure`.** Early-stage: definition of Gibbs
  measures, aiming at existence and uniqueness
  <https://reservoir.lean-lang.org/@YaelDillies/GibbsMeasure> [snippet;
  unverified].
- No Mathlib-level Dobrushin theorem was found. No Coq or Isabelle search
  was done. **[unverified]** that none exists.

Consequence: formalising the finite-state CEC statements (convex-combination
enclosure, Dinkelbach tests, sharpness) is feasible but is *not* formalising
a Gibbs-measure theorem. The formal claim must be scoped to the finite
algebraic statements. Do not say "first formal verification of Gibbs-measure
bounds"; the above projects exist.

## 11. Most defensible contribution

Not the identities or the basic enclosures. In order of defensibility:

1. **An executable, machine-checked certificate kernel** for ratio
   enclosures: exact-rational checkers for T1 (the `g(λ)` test) and T1′
   (LP duals) and the separator enclosure, with a Lean-verified soundness
   proof, in a setting where published provable-design tools certify each Z
   separately. This is an engineering-plus-assurance contribution. Its
   strength depends on the BBK*/MARK* read confirming that those tools do
   not already certify ratios jointly.
2. **A clean impossibility/separation statement** (T2.2 + T2.3): for fixed
   local data, the ratio is confined to exactly `[min ρ, max ρ]`, every
   value is attained, and `Z_A` is wholly unconstrained, even with exteriors
   restricted to pairwise factors. That packages known mechanisms as a
   precise optimality statement. Frame it as a corollary, not a discovery.
3. **Class-relative sharp enclosure T2.7** with attainment, plus the
   observation that the answer changes (to a larger LP) when exterior
   factors may touch several boundary variables.
4. **Empirical demonstration of size-independence** of local-ratio widths
   versus extensive growth of separate-bound widths on models with an exact
   oracle, and a rigorous model-relative titration interval benchmark.

Anything beyond (1)–(4) needs a stronger literature check first.

## 12. Evidence an independent expert would need

- Full-text checks: BBK*/MARK* (how K* bounds are formed; any coupled ratio
  bound), Leisink–Kappen (is sharpness over blanket distributions stated;
  pairwise realisability), Mooij–Kappen long version.
- A head-to-head on identical instances against Bound Propagation / box
  propagation / SAW-tree bounds (for marginal-ratio quantities) and against
  ratio-of-separate-bounds K*-style intervals, reporting widths and runtime.
  No such comparison exists in this repository yet.
- Exact-oracle validation (full enumeration) including adversarial
  exteriors that attain the endpoints, using exact rationals.
- The Lean axiom footprint (`#print axioms`) and an independent review of
  the formal statements against the paper statements.
- For titration: comparison against MCCE/Monte Carlo and PKAD-based
  experimental pKa as *model* validation, separately from the certified
  model-relative intervals, with all cutoff and parameter choices disclosed.
- A documented search log, including the unsearched items listed at the top
  (HP mutation literature, Coq/Isabelle, pre-1990 mixing literature).

## 13. Phrasing rules (global)

Never write: "new identity", "novel bound on partition-function ratios",
"first rigorous pKa", "beats K*/OSPREY/MCCE", "predicts selectivity",
"proves hardness of pKa", "discovers spatial mixing", "certified free
energy" without "model-relative". Always say: "model-relative certificate",
"applies a classical identity conditionally", "sharp over all exteriors
consistent with the local data", and cite Zwanzig, Georgii (DLR),
Leisink–Kappen, Charnes–Cooper/Dinkelbach, Dobrushin, and Weitz where the
corresponding step is used.
