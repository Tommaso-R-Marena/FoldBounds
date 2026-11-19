# Vendored: tsucres/HPFoldingProblem (MINWALK)

Source: https://github.com/tsucres/HPFoldingProblem  
Upstream algorithm: Giaquinta & Pozzi, *J Comput Biol* 20(8):593–609 (2013),
doi:10.1089/cmb.2012.0266.

Vendored for head-to-head comparison against FoldBounds on the community
2D HP ground-state benchmarks. Local fix applied to `src/Grid.hpp`
(`Point` equality compared `rhs.y == rhs.y`; corrected to `lhs.y == rhs.y`).

This code is **not** part of the Lean-verified kernel trust boundary.
