import FoldBounds
import FoldBounds.MotionCoordinates
import FoldBounds.MotionVariableDirections
import FoldBounds.MotionInverseFourth
import FoldBounds.MotionProteinHessian
import FoldBounds.MotionResidual
import FoldBounds.MotionProteinResponse
import FoldBounds.MotionPairOrdering
import FoldBounds.MotionVerifiedInverse
import FoldBounds.MotionExactSolvability
import FoldBounds.MotionRealPairSamples
import FoldBounds.MotionScaledCertificates
import FoldBounds.MotionRoundedOps
import FoldBounds.MotionScaledOrdering
import FoldBounds.MotionBulkChain
import FoldBounds.MotionBulkRealSmoke
import FoldBounds.MotionDyadicDivision
import FoldBounds.MotionSparseResponse
import FoldBounds.MotionContactCompleteness
import FoldBounds.MotionRobustResponse
import FoldBounds.MotionRobustOrdering
import FoldBounds.MotionRobustConcordance
import FoldBounds.MotionGNMExact
import FoldBounds.MotionGNMTopologyRobust
import FoldBounds.MotionGNMAdaptiveCases

/-! Axiom audit: every headline theorem depends only on the standard axioms
`propext`, `Classical.choice`, `Quot.sound`. -/

#print axioms FoldBounds.Z_pos
#print axioms FoldBounds.prob_mem_Icc
#print axioms FoldBounds.Z_mem_Icc_of_blocks
#print axioms FoldBounds.Z_mem_Icc_of_uniform
#print axioms FoldBounds.upper_antitone_enumeration
#print axioms FoldBounds.lower_monotone_enumeration
#print axioms FoldBounds.blockBound_refine_partition
#print axioms FoldBounds.blockBound_sub_le_mul
#print axioms FoldBounds.prob_mem_naive
#print axioms FoldBounds.prob_mem_coupled
#print axioms FoldBounds.prob_mem_of_blocks
#print axioms FoldBounds.coupled_lower_isLeast
#print axioms FoldBounds.coupled_upper_isGreatest
#print axioms FoldBounds.coupled_upper_lt_naive_upper
#print axioms FoldBounds.coupled_width_le
#print axioms FoldBounds.prob_width_of_blocks
#print axioms FoldBounds.coupled_enclosure_mono
#print axioms FoldBounds.prob_enclosure_refine
#print axioms FoldBounds.freeEnergy_mem_Icc
#print axioms FoldBounds.freeEnergy_mem_Icc_of_blocks
#print axioms FoldBounds.freeEnergy_le_energy
#print axioms FoldBounds.le_freeEnergy_of_energy_ge
#print axioms FoldBounds.freeEnergy_diff_mem_Icc
#print axioms FoldBounds.Z_eq_sum_orbits
#print axioms FoldBounds.Z_mem_Icc_of_reduced
#print axioms FoldBounds.Z_coupledEnergy_mem_Icc
#print axioms FoldBounds.freeEnergy_coupledEnergy_mem_Icc
#print axioms FoldBounds.prob_mem_Icc_of_rat
#print axioms FoldBounds.toy_prob_mem
#print axioms FoldBounds.prob_lower_endpoint_attained
#print axioms FoldBounds.prob_upper_endpoint_attained
#print axioms FoldBounds.Z_ge_gibbs_bogoliubov
#print axioms FoldBounds.freeEnergy_le_gibbs_bogoliubov
#print axioms FoldBounds.Z_coupledEnergy_ge_gibbs
#print axioms FoldBounds.Executable.massQ_mem_Icc
#print axioms FoldBounds.Executable.accepts_sound
#print axioms FoldBounds.Executable.check_sound
#print axioms FoldBounds.Executable.accepts_weights_pos
#print axioms FoldBounds.Executable.accepts_coupled_probability_sound
#print axioms FoldBounds.Executable.check_implies_accepts
#print axioms FoldBounds.Executable.check_coupled_probability_sound
#print axioms FoldBounds.LatticeModel.conformation_nodup
#print axioms FoldBounds.LatticeModel.conformation_length
#print axioms FoldBounds.LatticeModel.checkHP_coupled_sound
#print axioms FoldBounds.LatticeModel.advance_adjacent
#print axioms FoldBounds.LatticeModel.trace_adjacentChain
#print axioms FoldBounds.LatticeModel.conformation_adjacentChain
#print axioms FoldBounds.LatticeModel.advance_injective
#print axioms FoldBounds.LatticeModel.trace_injective
#print axioms FoldBounds.LatticeModel.walkFrom_injective

/-! Executable exact-rational contact and coordinate-box certification. -/
#print axioms FoldBounds.MotionContactChecker.accepts_valid
#print axioms FoldBounds.MotionContactChecker.weight_interval_sound
#print axioms FoldBounds.MotionContactChecker.check_valid
#print axioms FoldBounds.MotionContactChecker.check_energy_sound
#print axioms FoldBounds.MotionContactChecker.all_admissible_energy_sound
#print axioms FoldBounds.MotionCoordinates.difference_sound
#print axioms FoldBounds.MotionCoordinates.squared_sound
#print axioms FoldBounds.MotionCoordinates.squaredDistance_sound
#print axioms FoldBounds.MotionCoordinates.contactBoxes_admissible
#print axioms FoldBounds.MotionCoordinates.coordinate_contact_energy_sound

/-! Variable-direction rank-one contact energy certificates. -/
#print axioms FoldBounds.MotionVariableDirections.addBox_sound
#print axioms FoldBounds.MotionVariableDirections.scaleBox_sound
#print axioms FoldBounds.MotionVariableDirections.stretchBox_sound
#print axioms FoldBounds.MotionVariableDirections.energyTerm_sound
#print axioms FoldBounds.MotionVariableDirections.variable_direction_energy_sound

/-! Inverse-fourth physical spring model, collision-safe executable gate,
rank-one PSD energy and regularized coercivity. -/
#print axioms FoldBounds.MotionInverseFourth.squared_lower_nonneg
#print axioms FoldBounds.MotionInverseFourth.squared_ratio_enclosure
#print axioms FoldBounds.MotionInverseFourth.safeNetwork_valid
#print axioms FoldBounds.MotionInverseFourth.weightedEnergyTerm_sound
#print axioms FoldBounds.MotionInverseFourth.safeNetwork_energy_sound
#print axioms FoldBounds.MotionInverseFourth.weightedEnergyTerm_nonneg
#print axioms FoldBounds.MotionInverseFourth.weightedNetworkEnergy_nonneg
#print axioms FoldBounds.MotionInverseFourth.weightedRidgeEnergy_coercive
#print axioms FoldBounds.MotionInverseFourth.weightedNetworkEnergy_interval

#print axioms FoldBounds.MotionInverseFourth.checkNetworkEnergy_sound

/-! Exact rational Hessian, physical Cartesian incidence, and coercivity. -/
#print axioms FoldBounds.MotionHessian.hessianQ_symmetric
#print axioms FoldBounds.MotionHessian.matVecQ_hessian
#print axioms FoldBounds.MotionHessian.dot_matVec_hessian
#print axioms FoldBounds.MotionHessian.hessianQ_coercive
#print axioms FoldBounds.MotionProteinHessian.incidence_dot
#print axioms FoldBounds.MotionProteinHessian.contactWeight_nonneg
#print axioms FoldBounds.MotionProteinHessian.weight_stretch_energy
#print axioms FoldBounds.MotionProteinHessian.flat_dot
#print axioms FoldBounds.MotionProteinHessian.proteinHessian_energy
#print axioms FoldBounds.MotionProteinHessian.proteinHessian_symm
#print axioms FoldBounds.MotionProteinHessian.proteinHessian_coercive

/-! Exact rational Hessian residual bounds and executable Option checker. -/
#print axioms FoldBounds.MotionResidual.dotQ_comm
#print axioms FoldBounds.MotionResidual.dotQ_vsub_left
#print axioms FoldBounds.MotionResidual.dotQ_vsub_right
#print axioms FoldBounds.MotionResidual.matVecQ_vsub
#print axioms FoldBounds.MotionResidual.dotQ_self_nonneg
#print axioms FoldBounds.MotionResidual.dotQ_matVec_symmetric
#print axioms FoldBounds.MotionResidual.dual_gap_identity
#print axioms FoldBounds.MotionResidual.residual_square_ineq
#print axioms FoldBounds.MotionResidual.symmetric_residual_enclosure
#print axioms FoldBounds.MotionResidual.rankOne_residual_enclosure
#print axioms FoldBounds.MotionResidual.checkResponse_sound

/-! Normalized protein pair-response gate, using no square roots. -/
#print axioms FoldBounds.MotionProteinResponse.checkPairScore_sound

/-! Fail-closed comparative predicted-motion rank certificates. -/
#print axioms FoldBounds.MotionPairOrdering.strictlyOrdered_sound
#print axioms FoldBounds.MotionPairOrdering.checkPairOrder_sound

/-! Exact rational right-inverse witness; assumption-free pair score. -/
#print axioms FoldBounds.MotionVerifiedInverse.checkRightInverse_valid
#print axioms FoldBounds.MotionVerifiedInverse.checked_inverse_solves
#print axioms FoldBounds.MotionVerifiedInverse.checkPairScoreWithInverse_sound

/-! Finite positive-ridge rational Hessians have unique exact solutions;
accepted score intervals now have a uniquely existing ground-truth response. -/
#print axioms FoldBounds.MotionExactSolvability.dotQ_self_eq_zero_iff
#print axioms FoldBounds.MotionExactSolvability.coercive_matrix_injective
#print axioms FoldBounds.MotionExactSolvability.coercive_matrix_solvable
#print axioms FoldBounds.MotionExactSolvability.rankOne_solvable
#print axioms FoldBounds.MotionExactSolvability.checkResponse_exists_unique_sound
#print axioms FoldBounds.MotionExactSolvability.checkPairScore_exists_unique_sound

/-! Two frozen crystal-derived exact rational distance checks. -/
#print axioms FoldBounds.MotionRealPairSamples.cspb_frozen_sq_distance
#print axioms FoldBounds.MotionRealPairSamples.crambin_frozen_sq_distance
#print axioms FoldBounds.MotionRealPairSamples.cspb_reference_pair_eligible
#print axioms FoldBounds.MotionRealPairSamples.crambin_reference_pair_eligible

/-! Actual acceptance of external fixed-point intervals by exact model. -/
#print axioms FoldBounds.MotionScaledCertificates.binaryScale_pos
#print axioms FoldBounds.MotionScaledCertificates.checkValue_sound
#print axioms FoldBounds.MotionScaledCertificates.checkEnclosure_sound
#print axioms FoldBounds.MotionScaledCertificates.checkScaledPairScore_exists_unique_sound

/-! Outward 2^96 arithmetic: integer-gated signed multiplication soundness. -/
#print axioms FoldBounds.MotionRoundedOps.scaleInt_cast
#print axioms FoldBounds.MotionRoundedOps.scaleInt_pos
#print axioms FoldBounds.MotionRoundedOps.addBox_sound
#print axioms FoldBounds.MotionRoundedOps.negBox_sound
#print axioms FoldBounds.MotionRoundedOps.checkMul_sound

/-! Verified strict ordering from two accepted scaled integer score enclosures. -/
#print axioms FoldBounds.MotionScaledOrdering.checkScaledOrder_sound

/-! Verifiable bulk ranking arithmetic on true frozen endpoint examples. -/
#print axioms FoldBounds.MotionBulkChain.separated_sound
#print axioms FoldBounds.MotionBulkChain.checkAdjacent_sound
#print axioms FoldBounds.MotionBulkChain.separated_scaled_sound
#print axioms FoldBounds.MotionBulkChain.scores_ordered_of_certificates
#print axioms FoldBounds.MotionBulkRealSmoke.cspb_frozen_rank_arithmetic
#print axioms FoldBounds.MotionBulkRealSmoke.crambin_frozen_rank_arithmetic
#print axioms FoldBounds.MotionBulkRealSmoke.cspb_frozen_strict
#print axioms FoldBounds.MotionBulkRealSmoke.crambin_frozen_strict

/-! Signed exact integer quotient/remainder enclosure. -/
#print axioms FoldBounds.MotionDyadicDivision.checkQuotient_sound

/-! Sparse true-protein Hessian residual and fail-closed score soundness. -/
#print axioms FoldBounds.MotionSparseResponse.sparseRowTerm_eq
#print axioms FoldBounds.MotionSparseResponse.sparseMatVec_eq
#print axioms FoldBounds.MotionSparseResponse.checkSparseScaledPair_sound

/-! Full cutoff graph validation: no missing, spurious or duplicate contacts. -/
#print axioms FoldBounds.MotionContactCompleteness.checkContactNetwork_sound
#print axioms FoldBounds.MotionContactCompleteness.checked_contact_present
#print axioms FoldBounds.MotionContactCompleteness.checked_contact_active
#print axioms FoldBounds.MotionContactCompleteness.checkCompleteSparseScore_sound

/-! Uniform response validity across arbitrary structural uncertainty and contact switches. -/
#print axioms FoldBounds.MotionRobustResponse.uniform_score_sound

/-! Uniform strict rank ordering across variable-coordinate and variable-contact worlds. -/
#print axioms FoldBounds.MotionRobustOrdering.uniform_strict_pair_order

/-! Aggregate strict rank concordance over structural worlds. -/
#print axioms FoldBounds.MotionRobustConcordance.alignedSign_bounds
#print axioms FoldBounds.MotionRobustConcordance.alignedSign_interval
#print axioms FoldBounds.MotionRobustConcordance.sum_alignedSign_interval

/-! Exact scalar GNM resistance and two-extreme-network soundness. -/
#print axioms FoldBounds.MotionGNMExact.checkExactGNM_sound
#print axioms FoldBounds.MotionGNMExact.gnm_residual_enclosure
#print axioms FoldBounds.MotionGNMTopologyRobust.scalarNetworkEnergy_mono
#print axioms FoldBounds.MotionGNMTopologyRobust.scalarResponse_antitone
#print axioms FoldBounds.MotionGNMTopologyRobust.scalarResponse_extremes

/-! Case-split graph dominance and universal protein rank concordance. -/
#print axioms FoldBounds.MotionGNMAdaptiveCases.case_contact_sandwich
#print axioms FoldBounds.MotionGNMAdaptiveCases.case_response_sandwich
#print axioms FoldBounds.MotionGNMAdaptiveCases.all_contact_worlds_concordance
