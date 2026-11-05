"""Exact structural-confidence ceiling from independent local 3D rewiring.

An *arbitrary* canonical simple-cubic self-avoiding walk is accepted only
after complete geometry validation. Whenever eight consecutive residues
fill exactly one 2x2x2 lattice cube, that fragment can be replaced by
other Hamiltonian paths with its endpoints fixed. Its occupied sites and
outside chain bonds do not change.

For a supplied strictly positive RATIONAL species-pair contact-weight
matrix, calculate the exact nominal Boltzmann weight ratio of each local
rearrangement without enumerating the global ensemble. Select disjoint
rewiring cubes with no lattice adjacency between cubes. The resulting
family has exactly factorizable relative partition mass

    S = product_i sum_{v in allowed_local_paths(i)} W(v)/W(original_i).

Each choice produces a DIFFERENT complete self-avoiding walk, and the
family contains the original path with relative mass exactly one.
If all microstate weights may independently be multiplied by factors in
[lo,hi] (0<lo<=hi), the exact target's equilibrium probability obeys

    P(target) <= hi / (hi + lo*(S-1)).

No Monte Carlo, full-ensemble enumeration, thermodynamic independence,
positive-contact-attraction assumption, or baseline energy calibration
is required. A symmetry-quotiented conservative ceiling uses at most 48
cubic signed-coordinate permutations per fixed anchored target.

This certificate is Python-checked mathematics for a simple-cubic
nearest-neighbor pairwise interaction model. It is NOT an atomistic
folding guarantee or Lean-kernel verification of the 3D geometry.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import product
from typing import Mapping

from .cubic_tiling import OFFSETS, ROOT, macro_snake, tile_paths

Point = tuple[int, int, int]
Path = tuple[Point, ...]
Pair = tuple[str, str]


def hp_pair_weights(q: Fraction) -> dict[Pair, Fraction]:
    """Declared HP model: each nonconsecutive HH contact has weight q>0."""
    if not isinstance(q, Fraction) or q <= 0:
        raise ValueError("HP contact weight must be a positive exact Fraction")
    return {("H", "H"): q,
            ("H", "P"): Fraction(1),
            ("P", "P"): Fraction(1)}


@dataclass(frozen=True)
class StructuralModel:
    """Untrusted proposal; checker recomputes every spatial obligation."""
    sequence: str
    target: Path
    pair_weights: Mapping[Pair, Fraction]
    correction_lower: Fraction = Fraction(1)
    correction_upper: Fraction = Fraction(1)


@dataclass(frozen=True)
class RewirePatch:
    start: int
    cube_origin: Point
    variants: tuple[Path, ...]              # full global coordinates in 8-site cube
    relative_weights: tuple[Fraction, ...]  # corresponding ratios, original = 1
    partition_factor: Fraction
    equal_weight_choices: int

    @property
    def choices(self) -> int:
        return len(self.variants)

    @property
    def occupied(self) -> frozenset[Point]:
        return frozenset(self.variants[0])


@dataclass(frozen=True)
class StructuralConfidence:
    residues: int
    candidate_cubes: int
    selected_cubes: tuple[RewirePatch, ...]
    constructive_conformations: int
    exactly_isoenergetic_conformations: int
    relative_partition_lower: Fraction
    microstate_probability_upper: Fraction
    rigid_fold_probability_upper: Fraction
    correction_lower: Fraction
    correction_upper: Fraction

    @property
    def rules_out_single_structure_majority(self) -> bool:
        return self.microstate_probability_upper < Fraction(1, 2)

    @property
    def rules_out_rigid_fold_majority(self) -> bool:
        return self.rigid_fold_probability_upper < Fraction(1, 2)


def _pair(a: str, b: str) -> Pair:
    return (a, b) if a <= b else (b, a)


def _validate(model: StructuralModel) -> dict[Point, int]:
    if not isinstance(model, StructuralModel):
        raise ValueError("expected a StructuralModel")
    seq, path = model.sequence, model.target
    if not isinstance(seq, str) or not seq or len(seq) < 2:
        raise ValueError("sequence must contain >=2 single-letter residue types")
    if not isinstance(path, tuple) or len(path) != len(seq) or path[:2] != ROOT:
        raise ValueError("target must be complete canonical-rooted 3D walk")
    for point in path:
        if (not isinstance(point, tuple) or len(point) != 3 or
                any(type(x) is not int for x in point)):
            raise ValueError("all 3D coordinates must be integer triples")
    if len(set(path)) != len(path):
        raise ValueError("target is not self-avoiding")
    for left, right in zip(path, path[1:]):
        if sum(abs(x - y) for x, y in zip(left, right)) != 1:
            raise ValueError("target contains a non-unit step")
    if (not isinstance(model.correction_lower, Fraction) or
            not isinstance(model.correction_upper, Fraction) or
            not 0 < model.correction_lower <= model.correction_upper):
        raise ValueError("corrections must be ordered, positive exact Fractions")
    if not isinstance(model.pair_weights, Mapping):
        raise ValueError("pair_weights must be an explicit mapping")
    alphabet = sorted(set(seq))
    expected = {_pair(a, b) for a, b in product(alphabet, repeat=2)}
    if not expected.issubset(set(model.pair_weights)):
        raise ValueError("pair_weights must explicitly cover every unordered pair of sequence types")
    if any(not isinstance(v, Fraction) or v <= 0
           for v in model.pair_weights.values()):
        raise ValueError("pair weights must be positive exact Fractions")
    return {point: i for i, point in enumerate(path)}


def _cube_origin(path: Path) -> Point | None:
    if len(path) != 8 or len(set(path)) != 8:
        return None
    origin = tuple(min(point[i] for point in path) for i in range(3))
    if any(max(point[i] for point in path) != origin[i] + 1
           for i in range(3)):
        return None
    return origin  # eight unique integer points in a 2x2x2 box


def _local_weight(model: StructuralModel, start: int, patch_path: Path,
                  occupied: Mapping[Point, int]) -> Fraction:
    """Product of ALL affected noncovalent pair factors, counted exactly once.

    Only lattice neighbor edges with >=1 endpoint inside the changed cube
    can be affected. Replacing a Hamiltonian path permutes residue indices
    at those eight fixed coordinates; all other site identities are fixed.
    Both internal and external edges are handled, and each covalent
    adjacency is explicitly excluded based on sequence indices.
    """
    assigned = {point: start + i for i, point in enumerate(patch_path)}
    if len(assigned) != 8 or set(assigned) != set(model.target[start:start + 8]):
        raise AssertionError("unverified local cube occupancy")
    seen: set[tuple[Point, Point]] = set()
    mass = Fraction(1)
    for point in assigned:
        for delta in OFFSETS:
            neighbor = tuple(a + b for a, b in zip(point, delta))
            if neighbor not in occupied:
                continue
            edge = (point, neighbor) if point < neighbor else (neighbor, point)
            if edge in seen:
                continue
            seen.add(edge)
            j = assigned.get(point, occupied[point])
            k = assigned.get(neighbor, occupied[neighbor])
            if abs(j - k) <= 1:
                continue
            label = _pair(model.sequence[j], model.sequence[k])
            mass *= model.pair_weights[label]
    return mass


def find_rewire_patches(model: StructuralModel) -> tuple[RewirePatch, ...]:
    """Recalculate every legal anchored 2x2x2 Hamiltonian rewiring window."""
    occupied = _validate(model)
    fragments = tile_paths()  # exact 8-site local motif library
    out: list[RewirePatch] = []
    for start in range(len(model.sequence) - 7):
        if start == 1:  # would change the canonical first-bond endpoint
            continue
        original = model.target[start:start + 8]
        origin = _cube_origin(original)
        if origin is None:
            continue
        local_start = tuple(a - b for a, b in zip(original[0], origin))
        local_end = tuple(a - b for a, b in zip(original[-1], origin))
        allowed_local = [
            p for p in fragments
            if p[0] == local_start and p[-1] == local_end
            and (start != 0 or p[:2] == (
                tuple(a - b for a, b in zip(original[0], origin)),
                tuple(a - b for a, b in zip(original[1], origin))
            ))
        ]
        if len(allowed_local) < 2:
            continue
        variants = tuple(
            tuple(tuple(a + b for a, b in zip(origin, point)) for point in variant)
            for variant in allowed_local
        )
        if len(set(variants)) != len(variants) or original not in variants:
            raise AssertionError("local rewrite library incomplete or duplicated")
        original_weight = _local_weight(model, start, original, occupied)
        if original_weight <= 0:
            raise AssertionError("positive pair-weight model yielded zero weight")
        ratios = tuple(
            _local_weight(model, start, variant, occupied) / original_weight
            for variant in variants
        )
        if not all(r > 0 for r in ratios):
            raise AssertionError("invalid local reweighting")
        assert ratios[variants.index(original)] == 1
        out.append(RewirePatch(
            start, origin, variants, ratios,
            sum(ratios, Fraction(0)), sum(r == 1 for r in ratios)
        ))
    return tuple(out)


def _conflict_neighborhood(patch: RewirePatch) -> set[Point]:
    """Forbidden sites for any other patch: cube plus all lattice neighbors."""
    forbidden = set(patch.occupied)
    for point in patch.occupied:
        for delta in OFFSETS:
            forbidden.add(tuple(a + b for a, b in zip(point, delta)))
    return forbidden


def check_structural_confidence(model: StructuralModel) -> StructuralConfidence:
    """Certified exact-rational upper bound on a proposed structure's probability.

    Greedy selection is an UNTRUSTED optimization heuristic, but the
    checker independently verifies that its selected windows are
    pairwise spatially nonadjacent before using factorization. Every
    chosen patch can be changed independently of every other one.
    """
    patches = find_rewire_patches(model)
    forbidden: set[Point] = set()
    selected: list[RewirePatch] = []
    # Sorting by exact rational local gain, not floating-point heuristics.
    for patch in sorted(patches, key=lambda p: (p.partition_factor, -p.start),
                        reverse=True):
        if patch.occupied & forbidden:
            continue
        if any(patch.occupied & _conflict_neighborhood(other) for other in selected):
            raise AssertionError("nonadjacent patch invariant violated")
        selected.append(patch)
        forbidden.update(_conflict_neighborhood(patch))

    multiplicity = 1
    equal = 1
    factor = Fraction(1)
    for patch in selected:
        multiplicity *= patch.choices
        equal *= patch.equal_weight_choices
        factor *= patch.partition_factor

    lo, hi = model.correction_lower, model.correction_upper
    if factor < 1 or multiplicity < equal or equal < 1:
        raise AssertionError("invalid family partition factor")
    micro = hi / (hi + lo * (factor - 1))
    # At most |O_h|=48 rigid signed-coordinate permutations can map one
    # anchored, indexed cubic walk into a fixed congruence class. Any
    # such target-equivalent walk has the SAME pairwise contact energy.
    # Of the selected family, total nominal mass outside that class
    # is at least (S-48)*nominal(target) when S>48.
    if factor <= 48:
        orbit = Fraction(1)
    else:
        orbit = 48 * hi / (48 * hi + lo * (factor - 48))
    return StructuralConfidence(
        len(model.sequence), len(patches), tuple(sorted(selected, key=lambda p: p.start)),
        multiplicity, equal, factor, micro, orbit, lo, hi
    )


def minimum_selectivity_for_confidence(equal_weight_states: int,
                                       target_confidence: Fraction) -> Fraction:
    """Necessary extra target-vs-each-alternative weight advantage.

    If D distinct states have identical baseline energy and no extra
    model term distinguishes them, the target can have at most 1/D
    probability (or the bounded-correction generalization above).

    To obtain P(target)>=p, a new discriminatory interaction must
    permit a ratio target/each_other of at least (D-1)*p/(1-p)
    in the most favorable case where ALL alternatives are no heavier
    than that relative minimum. This is an INFORMATION-RESOLUTION
    requirement, not a sufficient design prescription.
    """
    if type(equal_weight_states) is not int or equal_weight_states < 1:
        raise ValueError("equal_weight_states must be a positive integer")
    if not isinstance(target_confidence, Fraction) or not 0 < target_confidence < 1:
        raise ValueError("target_confidence must be a Fraction in (0,1)")
    return Fraction(equal_weight_states - 1) * target_confidence / (1 - target_confidence)


def make_macro_target(nx: int, ny: int, nz: int) -> Path:
    """UNTRUSTED deterministic candidate: one Hamiltonian full-box 3D walk.

    The checker independently validates every coordinate. This helper
    chooses endpoints from the exact local 8-site motif enumeration.
    """
    if any(type(x) is not int or x < 1 for x in (nx, ny, nz)) or nx < 2:
        raise ValueError("positive macro dimensions required, nx>=2")
    route = macro_snake(nx, ny, nz)
    motifs = tile_paths()
    all_points: list[Point] = []
    entrance: Point = (0, 0, 0)
    for i, cell in enumerate(route):
        outgoing = None
        if i + 1 < len(route):
            outgoing = tuple(b - a for a, b in zip(cell, route[i + 1]))
            axis = next(axis for axis in range(3) if outgoing[axis])
            face = 1 if outgoing[axis] == 1 else 0
        else:
            axis = face = -1
        viable = [
            p for p in motifs
            if p[0] == entrance and (i != 0 or p[:2] == ROOT)
            and (outgoing is None or p[-1][axis] == face)
        ]
        if not viable:
            raise AssertionError("no local motif for verified macro step")
        chosen = min(viable)
        all_points.extend(
            tuple(2*a + b for a, b in zip(cell, point)) for point in chosen
        )
        if outgoing is not None:
            entrance = tuple(a - b for a, b in zip(chosen[-1], outgoing))
    candidate = tuple(all_points)
    if (candidate[:2] != ROOT or len(candidate) != 8 * nx * ny * nz
            or len(set(candidate)) != len(candidate)
            or any(sum(abs(a - b) for a, b in zip(x, y)) != 1
                   for x, y in zip(candidate, candidate[1:]))):
        raise AssertionError("macro Hamiltonian witness generator failed")
    return candidate
