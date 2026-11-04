import FoldBounds
import FoldBounds.MotionCoordinates
import FoldBounds.MotionVariableDirections
import FoldBounds.MotionInverseFourth
import FoldBounds.MotionProteinHessian
import FoldBounds.MotionResidual
import FoldBounds.MotionProteinResponse
import FoldBounds.MotionPairOrdering

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
