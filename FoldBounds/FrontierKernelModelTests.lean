import FoldBounds.FrontierKernelModelTestsHHHH
import FoldBounds.FrontierKernelModelTestsH8Interval

/-!
# Semantic bridge tests for the frontier kernel (aggregator)

Individual theorems live in separate modules so each `decide +kernel` goal
elaborates in its own Lean process.  A full root-`X` decide for `n = 8`
(543 conformations) while `FrontierKernelSound` is loaded exceeds the
memory cgroup; exactness for `n = 8` is covered by the Mathlib-free
`FrontierKernelTests` suite and by the compiled `foldkernel` binary, and
the semantic bridge for `n = 8` is the certified enclosure below.
-/
