"""Classical lattice-isoperimetric ground-state degeneracy no-go.

For n=L^3 all-H residues on the simple cubic lattice, every complete
self-avoiding walk has <=2*n-3*L^2+1 nonconsecutive nearest-neighbor
contacts, by the discrete Loomis-Whitney / AM-GM edge-isoperimetric
inequality. The LxLxL box attains equality, so any complete Hamiltonian
walk through it is a GLOBAL contact-maximizing state.

Combine that established inequality with the separately checked
local-rewiring certificate: a globally optimal Hamiltonian target can
have exponentially many structurally different, exactly isoenergetic
competitors. This bounds the equilibrium probability of one exact
ground-state microstate for ALL positive contact weights.

The edge-isoperimetric theorem is an external established mathematical
result used here as a soundness premise; this Python wrapper is NOT
a Lean proof of discrete Loomis-Whitney or an experimental protein
folding result.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .local_rewire_confidence import (
    StructuralModel, StructuralConfidence, check_structural_confidence,
    hp_pair_weights, make_macro_target,
)
from .cubic_tiling import OFFSETS


@dataclass(frozen=True)
class CubicGroundStateNoGo:
    box_side: int
    residues: int
    global_contact_upper: int
    target_contacts: int
    globally_optimal_for_attractive_hh: bool
    confidence: StructuralConfidence


def certify_cubic_groundstate_no_go(
    box_side: int, *, q: Fraction = Fraction(4),
    correction_lower: Fraction = Fraction(1),
    correction_upper: Fraction = Fraction(1),
) -> CubicGroundStateNoGo:
    """A geometric optimum and a checked lower bound on isoenergetic rivals.

    The user is not allowed to submit unverified contact counts or
    degeneracy numbers; the complete target is generated then
    independently checked by the structural confidence checker.

    The **known** edge-isoperimetric inequality is relied upon, rather
    than claimed to have been proved here. For q>1 this target is a
    true global ground state in the stated all-H contact Hamiltonian.
    For q<1 it is only a valid high-contact target (not a ground state).
    """
    if type(box_side) is not int or box_side < 4 or box_side % 2:
        raise ValueError("box_side must be an even integer >=4")
    if not isinstance(q, Fraction) or q <= 0:
        raise ValueError("contact q must be a positive exact Fraction")
    path = make_macro_target(box_side // 2, box_side // 2, box_side // 2)
    n = box_side**3
    occupied = set(path)
    if len(path) != n or len(occupied) != n:
        raise AssertionError("full-cube canonical target size failure")
    if occupied != {
        (x, y, z) for x in range(box_side)
        for y in range(box_side) for z in range(box_side)
    }:
        raise AssertionError("Hamiltonian target does not occupy perfect cube")

    # Independently reconstruct all occupied nearest-neighbor lattice
    # edges, then subtract the n-1 known covalent bonds. Every edge is
    # counted exactly once by checking only positive unit directions.
    edges = sum(
        (x + 1, y, z) in occupied
        + ((x, y + 1, z) in occupied)
        + ((x, y, z + 1) in occupied)
        for x, y, z in occupied
    )
    contacts = edges - (n - 1)
    known_global_upper = 2 * n - 3 * box_side**2 + 1
    if contacts != known_global_upper:
        raise AssertionError("box does not saturate the classical isoperimetric bound")

    seq = "H" * n
    result = check_structural_confidence(StructuralModel(
        seq, path, hp_pair_weights(q),
        correction_lower, correction_upper,
    ))
    if (result.exactly_isoenergetic_conformations !=
            result.constructive_conformations):
        raise AssertionError("all-H local rewires must be energy-invariant")
    if result.relative_partition_lower != result.constructive_conformations:
        raise AssertionError("all-H local rewires changed contact energy")
    return CubicGroundStateNoGo(
        box_side, n, known_global_upper, contacts, q > 1, result
    )
