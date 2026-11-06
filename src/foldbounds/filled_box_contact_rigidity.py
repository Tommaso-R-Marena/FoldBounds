"""Proof-auditable reconstruction of full 2D/3D grid contact maps.

A completely occupied cubic/square grid has its entire (labeled) nearest-
neighbor graph encoded by backbone edges plus its *full* noncovalent contact
map. Recover its unique geometry up to global box automorphisms using graph
metric landmarks, then validate every adjacency and all occupied sites.

The graph-product automorphism claim is a mathematical premise, NOT a
Lean-kernel theorem. No physical protein or experimental prediction is made.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import permutations, product
from math import factorial
from fractions import Fraction
from typing import Iterable


class FilledBoxCertificateError(ValueError):
    """Invalid contact data, incompatible geometry or non-canonical walk."""


Point = tuple[int, ...]
Walk = tuple[Point, ...]
Pair = tuple[int, int]


@dataclass(frozen=True)
class FilledBoxContactCertificate:
    dimension: int
    side: int
    nonbackbone_contacts: tuple[Pair, ...]
    representative: Walk
    orbit: tuple[Walk, ...]

    @property
    def orbit_size(self) -> int:
        return len(self.orbit)


def _validate_shape(side: int, dimension: int) -> int:
    if type(side) is not int or type(dimension) is not int or dimension not in (2, 3):
        raise FilledBoxCertificateError("integer side and dimension 2 or 3 required")
    if side < 2 or side ** dimension > 256:
        raise FilledBoxCertificateError("side >= 2 and at most 256 sites required")
    return side ** dimension


def _normalize_contacts(raw: Iterable[Pair], n: int) -> tuple[Pair, ...]:
    if isinstance(raw, (str, bytes)):
        raise FilledBoxCertificateError("contact pairs required")
    try:
        entries = tuple(raw)
    except TypeError as e:
        raise FilledBoxCertificateError("contacts must be iterable") from e
    parsed = []
    for pair in entries:
        if (not isinstance(pair, tuple) or len(pair) != 2 or
            any(type(x) is not int for x in pair)):
            raise FilledBoxCertificateError("contacts must be pairs of true integers")
        i, j = pair
        if not 0 <= i < j < n or j - i < 2:
            raise FilledBoxCertificateError("contacts must be ordered nonbackbone index pairs")
        if (j - i) % 2 == 0:
            raise FilledBoxCertificateError("grid contacts require odd sequence separation")
        parsed.append((i, j))
    if len(set(parsed)) != len(parsed):
        raise FilledBoxCertificateError("duplicate contact pair")
    return tuple(sorted(parsed))


def _graph(n: int, contacts: tuple[Pair, ...]) -> tuple[frozenset[int], ...]:
    g = [set() for _ in range(n)]
    for i in range(n - 1):
        g[i].add(i + 1)
        g[i + 1].add(i)
    for i, j in contacts:
        g[i].add(j)
        g[j].add(i)
    return tuple(frozenset(neighbors) for neighbors in g)


def _bfs_distances(g: tuple[frozenset[int], ...], source: int) -> tuple[int, ...]:
    dist = [-1] * len(g)
    dist[source] = 0
    q = deque([source])
    while q:
        v = q.popleft()
        for w in g[v]:
            if dist[w] < 0:
                dist[w] = dist[v] + 1
                q.append(w)
    if -1 in dist:
        raise FilledBoxCertificateError("graph is disconnected")
    return tuple(dist)


def _recover_embedding(g: tuple[frozenset[int], ...], m: int, d: int) -> Walk:
    """Reconstruct coordinates using corner-to-axis-end graph distances.

    In a filled grid P_m^d, a corner has d neighbors. Starting along one
    axis, its straight successor is the UNIQUE next neighbor sharing
    exactly one common neighbor with the vertex two steps behind. Turns
    share two common neighbors (a square plaquette). Axis endpoints lie
    m-1 steps from the corner.
    """
    n = len(g)
    expected_degrees = []
    for c in product(range(m), repeat=d):
        expected_degrees.append(sum(1 + (0 < x < m - 1) for x in c))
    if sorted(len(x) for x in g) != sorted(expected_degrees):
        raise FilledBoxCertificateError("incompatible filled-grid degree signature")
    corners = [i for i in range(n) if len(g[i]) == d]
    if len(corners) != 2 ** d:
        raise FilledBoxCertificateError("incorrect number of corners")
    c = min(corners)
    axis_neighbors = sorted(g[c])
    if len(axis_neighbors) != d:
        raise FilledBoxCertificateError("corner must have exactly d outgoing axes")
    opposite_axes = []
    for first in axis_neighbors:
        prev, now = c, first
        for _ in range(2, m):
            choices = [v for v in g[now] if v != prev and len(g[v] & g[prev]) == 1]
            if len(choices) != 1:
                raise FilledBoxCertificateError("no uniquely defined straight axis")
            prev, now = now, choices[0]
        opposite_axes.append(now)
    base = _bfs_distances(g, c)
    tips = [_bfs_distances(g, v) for v in opposite_axes]
    coords = []
    for v in range(n):
        entry = []
        for tip_dist in tips:
            twice = base[v] + m - 1 - tip_dist[v]
            if twice % 2:
                raise FilledBoxCertificateError("nonintegral reconstructed coordinate")
            x = twice // 2
            if not 0 <= x < m:
                raise FilledBoxCertificateError("reconstructed coordinate out of range")
            entry.append(x)
        coords.append(tuple(entry))
    if len(set(coords)) != n:
        raise FilledBoxCertificateError("reconstructed embedding is not bijective")
    occupied = {pos: i for i, pos in enumerate(coords)}
    if len(occupied) != n:
        raise FilledBoxCertificateError("missing grid positions")
    for i, pos in enumerate(coords):
        predicted = set()
        for axis in range(d):
            for delta in (-1, +1):
                u = list(pos)
                u[axis] += delta
                j = occupied.get(tuple(u))
                if j is not None:
                    predicted.add(j)
        if predicted != g[i]:
            raise FilledBoxCertificateError("contact graph is not a complete filled grid")
    return tuple(coords)


def _canonical_orbit(coords: Walk) -> tuple[Walk, ...]:
    """All isometries fixing residue 0 and the directed first +x step."""
    d = len(coords[0])
    diffs = [coords[1][j] - coords[0][j] for j in range(d)]
    axis = [j for j, v in enumerate(diffs) if v]
    if len(axis) != 1 or abs(diffs[axis[0]]) != 1:
        raise FilledBoxCertificateError("first backbone edge is not unit-length")
    first = axis[0]
    remaining = tuple(j for j in range(d) if j != first)
    orbit = []
    for perm in permutations(remaining):
        for signs in product((-1, +1), repeat=d - 1):
            orientation = ((first, diffs[first]),) + tuple(zip(perm, signs))
            path = tuple(tuple(sign * (v[j] - coords[0][j])
                               for j, sign in orientation) for v in coords)
            orbit.append(path)
    if len(set(orbit)) != 2 ** (d - 1) * factorial(d - 1):
        raise FilledBoxCertificateError("unexpected orbit stabilizer")
    return tuple(sorted(orbit))


def extract_nonbackbone_contacts(path: Walk) -> tuple[Pair, ...]:
    """Recompute all contacts using independent coordinate adjacency."""
    if not isinstance(path, tuple) or len(path) < 2:
        raise FilledBoxCertificateError("tuple of visited sites required")
    d = len(path[0])
    if d not in (2, 3) or any(not isinstance(v, tuple) or len(v) != d or
                              any(type(x) is not int for x in v) for v in path):
        raise FilledBoxCertificateError("invalid integer-coordinate path")
    if len(set(path)) != len(path):
        raise FilledBoxCertificateError("not a self-avoiding walk")
    if path[0] != (0,) * d or path[1] != (1,) + (0,) * (d - 1):
        raise FilledBoxCertificateError("path is not anchored with first +x edge")
    if any(sum(abs(a-b) for a, b in zip(v, w)) != 1
           for v, w in zip(path, path[1:])):
        raise FilledBoxCertificateError("non-lattice backbone edge")
    location = {pos: i for i, pos in enumerate(path)}
    found = set()
    for i, v in enumerate(path):
        for axis in range(d):
            for delta in (-1, 1):
                target = list(v)
                target[axis] += delta
                j = location.get(tuple(target))
                if j is not None and j > i + 1:
                    found.add((i, j))
    return tuple(sorted(found))


def certify_filled_box_contacts(side: int, dimension: int,
                                nonbackbone_contacts: Iterable[Pair]
                                ) -> FilledBoxContactCertificate:
    """Reject forged contact maps; reconstruct *every* anchored isometry.

    Correctness assumes the classical automorphism theorem for Cartesian
    products of paths: any isomorphism of a filled d-dimensional cube grid
    to itself is a signed permutation of its coordinate axes. The explicit
    reconstruction and rechecking require no graph-isomorphism oracle.
    """
    n = _validate_shape(side, dimension)
    pairs = _normalize_contacts(nonbackbone_contacts, n)
    g = _graph(n, pairs)
    coords = _recover_embedding(g, side, dimension)
    orbit = _canonical_orbit(coords)
    for path in orbit:
        if extract_nonbackbone_contacts(path) != pairs:
            raise FilledBoxCertificateError("reconstructed contact map mismatch")
        width = [max(v[j] for v in path)-min(v[j] for v in path)+1
                 for j in range(dimension)]
        if width != [side] * dimension:
            raise FilledBoxCertificateError("reconstruction fails to fill cube")
    return FilledBoxContactCertificate(dimension, side, pairs, orbit[0], orbit)


def verify_filled_box_certificate(cert: FilledBoxContactCertificate) -> bool:
    """Fail-closed replay of all geometry, contacts and symmetry multiplicity."""
    if not isinstance(cert, FilledBoxContactCertificate):
        return False
    try:
        expected = certify_filled_box_contacts(cert.side, cert.dimension,
                                               cert.nonbackbone_contacts)
    except (FilledBoxCertificateError, ValueError, TypeError, OverflowError):
        return False
    return expected == cert


def contact_map_confidence_upper(map_orbit_count: int, ground_paths: int, *,
                                 correction_lower: Fraction = Fraction(1),
                                 correction_upper: Fraction = Fraction(1)):
    """Universal Gibbs probability ceiling for a specified all-H ground map.

    The ground maps all have identical contact score, so for q>0 and any
    independent multipliers in [lower,upper], P(map) is at most
      u*g / (u*g + l*(C-g)),
    irrespective of other (non-ground) conformations.
    """
    if (type(map_orbit_count) is not int or type(ground_paths) is not int or
        not 1 <= map_orbit_count <= ground_paths or
        not isinstance(correction_lower, (Fraction, int)) or
        not isinstance(correction_upper, (Fraction, int)) or
        isinstance(correction_lower, bool) or isinstance(correction_upper, bool)):
        raise ValueError("positive integer counts and exact rational factors required")
    l, u = Fraction(correction_lower), Fraction(correction_upper)
    if not 0 < l <= u:
        raise ValueError("correction factors must obey 0<l<=u")
    return u * map_orbit_count / (u * map_orbit_count +
                                  l * (ground_paths - map_orbit_count))


def split_proper_and_mirror_orbits(
        cert: FilledBoxContactCertificate
        ) -> tuple[tuple[Walk, ...], tuple[Walk, ...]]:
    """Separate the 8 anchored cubic realizations into two chiral classes.

    Within the anchored frame first bond is +x. The four proper rotations
    around x preserve an enantiomer's spatial handedness; reflecting the
    transverse yz plane exchanges it. Pairwise distance/contact energies
    assign the *same energy* to the entire 8-fold orbit. No statement of
    chemical chirality in a real protein is inferred from this toy model.
    """
    if not isinstance(cert, FilledBoxContactCertificate) or cert.dimension != 3:
        raise FilledBoxCertificateError("a validated 3D cube certificate is required")
    if not verify_filled_box_certificate(cert):
        raise FilledBoxCertificateError("forged 3D contact certificate")
    p = cert.representative
    def oriented(sign: int) -> tuple[Walk, ...]:
        paths = []
        for step in range(4):
            transformed = []
            for x, y, z in p:
                y *= sign
                if step == 0:
                    yz = (y, z)
                elif step == 1:
                    yz = (-z, y)
                elif step == 2:
                    yz = (-y, -z)
                else:
                    yz = (z, -y)
                transformed.append((x,) + yz)
            paths.append(tuple(transformed))
        return tuple(sorted(paths))
    proper = oriented(+1)
    mirror = oriented(-1)
    if (len(set(proper)) != 4 or len(set(mirror)) != 4 or
        set(proper) & set(mirror) or
        set(proper) | set(mirror) != set(cert.orbit)):
        raise FilledBoxCertificateError("proper/mirror orbit partition mismatch")
    return proper, mirror


@dataclass(frozen=True)
class ContactEnergyInvariance:
    """Checked exact pair-energy equality over one entire filled-box orbit."""
    energy: Fraction
    orbit_size: int
    contact_count: int


def certify_pairwise_contact_energy(
        cert: FilledBoxContactCertificate,
        pair_energies: dict[Pair, int | Fraction]
        ) -> ContactEnergyInvariance:
    """Check arbitrary residue-specific rational nonbackbone contact energy.

    Every term depends only on the index-labeled contact pair, not spatial
    chirality, torsions, or directional interactions. All eight cubic (or two
    square) realizations must then have identical energy. Input coefficients
    may have either sign, as a general coarse-grained pair Hamiltonian.
    """
    if not verify_filled_box_certificate(cert):
        raise FilledBoxCertificateError("invalid filled-box certificate")
    if not isinstance(pair_energies, dict):
        raise FilledBoxCertificateError("energy coefficients must be a dictionary")
    n = len(cert.representative)
    for pair, value in pair_energies.items():
        if (not isinstance(pair, tuple) or len(pair) != 2 or
            any(type(j) is not int for j in pair) or
            not 0 <= pair[0] < pair[1] < n or
            pair[1] - pair[0] < 2 or
            type(value) not in (int, Fraction)):
            raise FilledBoxCertificateError("invalid rational pair-energy specification")
    energies = []
    for path in cert.orbit:
        pairs = extract_nonbackbone_contacts(path)
        energy = sum((Fraction(pair_energies.get(pair, 0)) for pair in pairs),
                     Fraction(0))
        energies.append(energy)
    if len(set(energies)) != 1:
        raise FilledBoxCertificateError("orbit members unexpectedly differ in contact energy")
    return ContactEnergyInvariance(energies[0], len(energies), len(cert.nonbackbone_contacts))


def verify_pairwise_contact_energy(cert: FilledBoxContactCertificate,
                                   pair_energies: dict[Pair, int | Fraction],
                                   energy_certificate: ContactEnergyInvariance) -> bool:
    if not isinstance(energy_certificate, ContactEnergyInvariance):
        return False
    try:
        independently = certify_pairwise_contact_energy(cert, pair_energies)
    except (ValueError, TypeError, OverflowError):
        return False
    return independently == energy_certificate
