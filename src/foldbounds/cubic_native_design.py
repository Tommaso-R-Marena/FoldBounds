"""Exact 3D cubic HP native-ground-orbit discovery and structural witnesses.

A native state here means a conformation maximizing H-H contacts in an
anchored cubic HP toy model. A single D4 orbit is a unique geometry modulo
lattice symmetries fixing the first bond; it is NOT a real protein fold.
Geometry completeness requires full independent audit (on by default).
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product

from .contact_map_atlas import ContactMapAtlas, audit_atlas
from .cubic_contact_atlas import CubicContactAtlas, audit_cubic_atlas

Point3 = tuple[int, int, int]
_DIRECTIONS: tuple[Point3, ...] = (
    (1,0,0), (-1,0,0), (0,1,0), (0,-1,0), (0,0,1), (0,0,-1)
)


@dataclass(frozen=True)
class UniqueCubicNativeOrbit:
    sequence: str
    max_contacts: int
    ground_conformations: int
    ground_contact_mask: int
    symmetry_orbit_size: int

    @property
    def hydrophobic_count(self) -> int:
        return self.sequence.count("H")


@dataclass(frozen=True)
class CubicNativeWitness:
    certificate: UniqueCubicNativeOrbit
    coordinates: tuple[Point3, ...]


def find_unique_cubic_native_orbits(
    cubic: CubicContactAtlas, planar: ContactMapAtlas, *,
    require_audit: bool = True, max_audit_leaves: int = 8_000_000,
    max_sequences: int = 4_096,
) -> tuple[UniqueCubicNativeOrbit, ...]:
    """Exhaustively identify sequences with one nontrivial D4 ground orbit.

    An anchored cubic orbit has size 1 (straight), 4 (planar), or 8
    (genuinely nonplanar). For a nonzero-score native state:
      * one contact map, total degeneracy 4 implies one planar orbit;
      * one contact map, total degeneracy 8 AND no planar realization
        implies one nonplanar orbit (not two planar orbits).

    The supplied 2D atlas rules out the latter ambiguity. Both atlases
    must have been independently checked; default does that here.
    When require_audit=False, the caller MUST already have validated the
    geometry with independently implemented full-walk auditors.
    """
    if not isinstance(cubic, CubicContactAtlas) or not isinstance(planar, ContactMapAtlas):
        raise TypeError("expected separately typed cubic and planar contact atlases")
    if cubic.length != planar.length:
        raise ValueError("cubic and planar model lengths must match")
    if type(require_audit) is not bool or type(max_sequences) is not int or max_sequences < 1:
        raise ValueError("invalid audit flag or search budget")
    if (1 << cubic.length) > max_sequences:
        raise ValueError("complete H/P sequence-space search exceeds budget")
    if require_audit:
        if not audit_atlas(planar):
            raise ValueError("planar geometry failed independent walk audit")
        if not audit_cubic_atlas(
            cubic, max_length=cubic.length, max_leaves=max_audit_leaves
        ):
            raise ValueError("cubic geometry failed independent walk audit")

    possible_planar = {mask for mask, _ in planar.multiplicities}
    pairs = cubic.pairs
    results = []
    for residues in range(1 << cubic.length):
        sequence = "".join("H" if residues & (1 << i) else "P"
                           for i in range(cubic.length))
        allowed = sum(1 << bit for bit, (i, j) in enumerate(pairs)
                      if residues & (1 << i) and residues & (1 << j))
        best = -1
        best_mass = 0
        ground_maps: list[tuple[int, int]] = []
        for mask, count in cubic.multiplicities:
            k = (mask & allowed).bit_count()
            if k > best:
                best = k
                best_mass = count
                ground_maps = [(mask, count)]
            elif k == best:
                best_mass += count
                ground_maps.append((mask, count))
        if best == 0 or len(ground_maps) != 1:
            continue
        mask, mass = ground_maps[0]
        if best_mass != mass:
            raise ArithmeticError("ground contact-map mass is inconsistent")
        if mass == 4 and mask in possible_planar:
            orbit = 4
        elif mass == 8 and mask not in possible_planar:
            orbit = 8
        else:
            continue
        results.append(UniqueCubicNativeOrbit(
            sequence, best, best_mass, mask, orbit
        ))
    return tuple(results)


def recover_cubic_native_witnesses(
    cubic: CubicContactAtlas,
    certificates: tuple[UniqueCubicNativeOrbit, ...], *,
    max_nodes: int = 3_000_000,
) -> tuple[CubicNativeWitness, ...]:
    """Independently obtain an explicit coordinate path for each orbit.

    Enumerates canonical paths again, checking exact mask and orbit weight
    at leaves. Never infers a 3D structure from a contact graph alone.
    """
    if not isinstance(cubic, CubicContactAtlas):
        raise TypeError("expected a cubic geometry atlas")
    if not isinstance(certificates, tuple) or any(
        not isinstance(c, UniqueCubicNativeOrbit)
        or len(c.sequence) != cubic.length for c in certificates
    ):
        raise ValueError("certificates must all apply to this cubic length")
    if type(max_nodes) is not int or max_nodes < 1:
        raise ValueError("invalid canonical DFS budget")
    targets = {c.ground_contact_mask for c in certificates}
    if len(targets) != len(certificates):
        raise ValueError("duplicate ground contact maps require separate orbit analysis")
    pair_bits = {p: 1 << k for k, p in enumerate(cubic.pairs)}
    occupied: dict[Point3, int] = {(0,0,0):0, (1,0,0):1}
    path = [(0,0,0), (1,0,0)]
    found: dict[int, list[tuple[int, tuple[Point3, ...]]]] = {
        mask: [] for mask in targets
    }
    nodes = 0

    def visit(x: int, y: int, z: int, sideways: bool, nonplanar: bool, mask: int):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise ValueError("orbit-witness enumeration exceeded budget")
        index = len(path)
        if index == cubic.length:
            if mask in targets:
                orbit_weight = 8 if nonplanar else 4 if sideways else 1
                found[mask].append((orbit_weight, tuple(path)))
            return
        for dx, dy, dz in _DIRECTIONS:
            if (not sideways and (dy < 0 or dz != 0)) or (sideways and not nonplanar and dz < 0):
                continue
            point = (x + dx, y + dy, z + dz)
            if point in occupied:
                continue
            next_mask = mask
            for ax, ay, az in _DIRECTIONS:
                prior = occupied.get((point[0]+ax, point[1]+ay, point[2]+az))
                if prior is not None and prior < index - 1:
                    next_mask |= pair_bits.get((prior, index), 0)
            occupied[point] = index
            path.append(point)
            visit(*point, sideways or dy != 0, nonplanar or dz != 0, next_mask)
            path.pop()
            del occupied[point]

    visit(1,0,0,False,False,0)
    available = dict(cubic.multiplicities)
    witnesses = []
    for c in certificates:
        candidate_orbits = found[c.ground_contact_mask]
        if (len(candidate_orbits) != 1
                or candidate_orbits[0][0] != c.symmetry_orbit_size
                or available[c.ground_contact_mask] != c.ground_conformations):
            raise ArithmeticError("canonical geometry does not prove one ground orbit")
        coords = candidate_orbits[0][1]
        selected = sum(1 << bit for bit, (i,j) in enumerate(cubic.pairs)
                       if c.sequence[i] == c.sequence[j] == "H")
        if (c.ground_contact_mask & selected).bit_count() != c.max_contacts:
            raise ArithmeticError("witness does not realize claimed H-H ground contacts")
        witnesses.append(CubicNativeWitness(c, coords))
    return tuple(witnesses)


def native_ground_probability(
    cubic: CubicContactAtlas, sequence: str, q: Fraction
) -> Fraction:
    """Exact Boltzmann probability of all maximally contacting conformations."""
    if type(q) is not Fraction or q < 1:
        raise ValueError("q must be an exact Fraction >=1")
    spectrum = cubic.spectrum(sequence)
    k = max(spectrum)
    z = sum((count * q**j for j, count in spectrum.items()), Fraction(0))
    return spectrum[k] * q**k / z


def isolate_native_occupancy_weight(
    cubic: CubicContactAtlas, sequence: str, target: Fraction, *,
    width: Fraction = Fraction(1, 1_000_000),
    max_steps: int = 256,
) -> tuple[Fraction, Fraction]:
    """Exact rational q bracket where ground-ensemble probability reaches target.

    P_ground = N_M / (N_M + sum_{k<M} N_k q^(k-M)), hence it is
    strictly increasing for q>0 if any lower-contact states exist.
    This proves global uniqueness of a crossing in (1,infinity),
    unlike an untrusted numerical q-grid.
    """
    if (type(target) is not Fraction or not 0 < target < 1
            or type(width) is not Fraction or width <= 0
            or type(max_steps) is not int or max_steps < 1):
        raise ValueError("target, bracket tolerance or work limit invalid")
    spectrum = cubic.spectrum(sequence)
    maximum = max(spectrum)
    if not any(k < maximum and count for k, count in spectrum.items()):
        raise ValueError("ground probability is constant, not a unique transition")
    def mass(q: Fraction) -> Fraction:
        return native_ground_probability(cubic, sequence, q)
    low = Fraction(1)
    if mass(low) >= target:
        return low, low
    high = Fraction(2)
    for _ in range(max_steps):
        if mass(high) >= target:
            break
        high *= 2
    else:
        raise ValueError("failed to bracket ground-state occupancy")
    for _ in range(max_steps):
        if high - low <= width:
            return low, high
        mid = (low + high) / 2
        v = mass(mid)
        if v == target:
            return mid, mid
        if v < target:
            low = mid
        else:
            high = mid
    raise ValueError("could not isolate target occupancy within budget")


def lattice_witness_pdb(
    witness: CubicNativeWitness, *, spacing_angstrom: Fraction = Fraction(19, 5)
) -> str:
    """PyMOL-compatible synthetic CA-only model, NOT an actual PDB structure."""
    if not isinstance(witness, CubicNativeWitness):
        raise TypeError("expected an explicit native orbit witness")
    if type(spacing_angstrom) is not Fraction or spacing_angstrom <= 0:
        raise ValueError("positive exact lattice visualization spacing required")
    sequence = witness.certificate.sequence
    coordinates = witness.coordinates
    if len(sequence) != len(coordinates) or len(set(coordinates)) != len(coordinates):
        raise ValueError("invalid self-avoiding coordinate witness")
    if any(sum(abs(x-y) for x,y in zip(coordinates[i], coordinates[i+1])) != 1
           for i in range(len(coordinates)-1)):
        raise ValueError("non-neighbor consecutive lattice points")
    lines = [
        "HEADER    SYNTHETIC CUBIC HP LATTICE CONFORMATION",
        "REMARK   1 MODEL IS NOT AN EXPERIMENTAL OR PREDICTED REAL PROTEIN",
        "REMARK   2 ALA REPRESENTS H; GLY REPRESENTS P; CA-ONLY PLACEHOLDERS",
        "REMARK   3 FIRST BOND ANCHORED; SYMMETRY-UNIQUE GROUND ORBIT",
        f"REMARK   4 SEQUENCE {sequence}",
        f"REMARK   5 CONTACTS {witness.certificate.max_contacts}",
    ]
    for i, (r, xyz) in enumerate(zip(sequence, coordinates), 1):
        residue = "ALA" if r == "H" else "GLY"
        x,y,z = (float(v*spacing_angstrom) for v in xyz)
        lines.append(
            f"ATOM  {i:5d}  CA  {residue} A{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}           C"
        )
    for i in range(1, len(sequence)):
        lines.append(f"CONECT{i:5d}{i+1:5d}")
    lines.extend(["TER", "END"])
    return "\n".join(lines) + "\n"
