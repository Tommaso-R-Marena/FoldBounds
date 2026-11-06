"""Couple 3D HP contact-event guarantees with exact target-fold confidence ceilings.

Both certificates MUST use the SAME complete sequence, canonical 3D target,
contact fugacity, per-state adversarial correction interval, and ensemble.

Result is a model-relative diagnostic: a high-probability coarse contact
event may coexist with a very low-probability exact conformation. An
event-population certificate is never promoted into fold confidence.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .cubic_tiling import TiledCertificate, TiledResult, check_tiled_family
from .local_rewire_confidence import (
    StructuralModel, StructuralConfidence, check_structural_confidence,
)

Point = tuple[int, int, int]
_NEIGHBORS = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1))


@dataclass(frozen=True)
class ContactVersusFoldResult:
    residues: int
    hydrophobic_residues: int
    polar_residues: int
    q: Fraction
    event_threshold: int
    target_hh_contacts: int
    contact_event_probability_lower: Fraction
    target_microstate_probability_upper: Fraction
    target_rigid_fold_probability_upper: Fraction
    contact_certificate: TiledResult
    structure_certificate: StructuralConfidence

    @property
    def proves_event_majority_without_unique_fold_majority(self) -> bool:
        return (
            self.contact_event_probability_lower > Fraction(1, 2)
            and self.target_microstate_probability_upper < Fraction(1, 2)
        )


def _target_hh_contacts(sequence: str, path: tuple[Point, ...]) -> int:
    coords = {p: i for i, p in enumerate(path)}
    count = 0
    for i, (x, y, z) in enumerate(path):
        if sequence[i] != "H":
            continue
        for dx, dy, dz in _NEIGHBORS:
            j = coords.get((x + dx, y + dy, z + dz))
            if j is not None and j > i + 1 and sequence[j] == "H":
                count += 1
    return count


def certify_contact_vs_fold(
    model: StructuralModel,
    macro_dimensions: Point,
    event_at_least: int,
) -> ContactVersusFoldResult:
    """Fail closed unless both certificate systems refer to one HP model.

    The event certificate quantifies all unrestricted 3D anchored SAWs;
    the constructed Hamiltonian family supplies just an event-mass lower
    bound. The structural checker bounds an exact target's microstate
    probability against all unrestricted SAWs with the SAME model.
    """
    if (not isinstance(model, StructuralModel)
            or not isinstance(model.sequence, str)
            or not model.sequence
            or set(model.sequence) - {"H", "P"}):
        raise ValueError("this coupling requires an H/P sequence")
    if (not isinstance(macro_dimensions, tuple)
            or len(macro_dimensions) != 3
            or any(type(x) is not int or x < 1 for x in macro_dimensions)):
        raise ValueError("macro_dimensions must be three positive integers")
    if type(event_at_least) is not int or event_at_least < 1:
        raise ValueError("event_at_least must be a positive integer")
    # Reject arbitrary pair matrices here: the event certificate assumes
    # only HH contacts are weighted. Full structural checker separately
    # supports arbitrary positive residue-pair matrices.
    weights = model.pair_weights
    if (not isinstance(weights, dict)
            or not isinstance(weights.get(("H", "H")), Fraction)
            or not isinstance(weights.get(("H", "P")), Fraction)
            or not isinstance(weights.get(("P", "P")), Fraction)
            or weights[("H", "P")] != 1
            or weights[("P", "P")] != 1
            or weights[("H", "H")] < 1):
        raise ValueError("coupled HP model requires exact q>=1 and HP=PP=1")

    shape = tuple(2*x for x in macro_dimensions)
    target = model.target
    # The self-avoiding structure checker performs complete coordinate
    # validation before the target is used as a certified result.
    confidence = check_structural_confidence(model)
    expected_sites = {
        (x, y, z)
        for x in range(shape[0])
        for y in range(shape[1])
        for z in range(shape[2])
    }
    if (len(target) != len(expected_sites)
            or set(target) != expected_sites):
        raise ValueError("target must completely occupy the declared cubic box")

    q = weights[("H", "H")]
    event = check_tiled_family(TiledCertificate(
        macro_dimensions, model.sequence, q, event_at_least,
        model.correction_lower, model.correction_upper
    ))
    target_contacts = _target_hh_contacts(model.sequence, target)
    if target_contacts < event_at_least:
        raise ValueError("target is outside the certified event")
    return ContactVersusFoldResult(
        len(model.sequence), model.sequence.count("H"),
        model.sequence.count("P"), q, event_at_least,
        target_contacts, event.probability_lower,
        confidence.microstate_probability_upper,
        confidence.rigid_fold_probability_upper, event, confidence
    )
