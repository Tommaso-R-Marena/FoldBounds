"""Exact vertex/hyperplane zonotope refinement for sparse real protein geometry.

The four-anchor squared-distance error relaxation gives a 3D zonotope.
A convex squared norm reaches its maximum on an error-box vertex; check
all <=256 vertices for an uncertain pair to certify a contact. A linear
supporting hyperplane supplies a safe lower distance for noncontacts.
All certification uses arbitrary-precision integers. This is classical
convex optimization, not a new theorem or physical protein predictor.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isqrt

from .sparse_distance_certificate import (
    SparseDistanceCertificate, reconstruct_exact, _basis4, _sub, _dot,
)


@dataclass(frozen=True)
class ZonotopeContactCertificate:
    all_pairs: int
    guaranteed_contacts: int
    guaranteed_noncontacts: int
    undecided_pairs: int
    nominal_contacts: int
    coarse_undecided_pairs: int
    new_contacts: int
    new_noncontacts: int
    vertex_bounds_checked: int
    hyperplane_bounds_checked: int


def _coefficient_vectors(cof):
    """Columns of the four signed squared-distance errors, scaled by 2 det."""
    return (
        tuple(sum(v[j] for v in cof) for j in range(3)),
        tuple(-v for v in cof[0]),
        tuple(-v for v in cof[1]),
        tuple(-v for v in cof[2]),
    )


def zonotope_contact_bounds(
    cert: SparseDistanceCertificate, *,
    cutoff_mA: int = 8000,
    epsilon_mA: int = 10,
    min_sequence_separation: int = 4,
) -> ZonotopeContactCertificate:
    """Check exact contact/noncontact invariance over a relaxed error zonotope.

    For nonanchor i, each of the four squared-distance errors is enclosed
    by [-eta_i,+eta_i], eta_i = 2 ceil(distance_max)*eps+eps**2. The true
    independent *distance* error ball maps into this overapproximation.

    For a pair (i,j), displacements lie in d + sum_k [-g_k,g_k]/(2 det).
    Checking every 8D error-box corner proves max ||d+error|| <= cutoff
    for contacts, while a signed projection lower bounds minimum pair
    distance for noncontacts. Tightening is attempted only on pairs left
    uncertain by the existing bounding-box certificate, preserving its
    exact results and trust boundary.
    """
    if (any(type(v) is not int for v in
            (cutoff_mA, epsilon_mA, min_sequence_separation))
        or cutoff_mA <= 0 or epsilon_mA < 0
        or min_sequence_separation < 1):
        raise ValueError("invalid error tolerance, cutoff, or separation")

    p = reconstruct_exact(cert)
    _, cof, det, norm = _basis4(cert.anchors)
    ad = abs(det)
    vectors = _coefficient_vectors(cof)
    anchors = set(cert.anchor_indices)

    eta, radii = [], []
    for i, row in enumerate(cert.squared_to_anchors):
        if i in anchors:
            eta.append(0)
            radii.append(0)
        else:
            d = max(isqrt(s) + (isqrt(s) ** 2 < s) for s in row)
            e = 2 * d * epsilon_mA + epsilon_mA ** 2
            eta.append(e)
            radii.append(e * norm)

    total = nominal = good = bad = uncertain = coarse_uncertain = 0
    gained_contact = gained_noncontact = vertices = projections = 0
    threshold = cutoff_mA * ad
    for i in range(len(p)):
        for j in range(i + min_sequence_separation, len(p)):
            total += 1
            d = _sub(p[i], p[j])
            original_squared = _dot(d, d)
            is_contact = original_squared <= cutoff_mA ** 2
            nominal += int(is_contact)
            margin = 2 * (radii[i] + radii[j])
            if (threshold >= margin and
                original_squared * ad ** 2 <= (threshold - margin) ** 2):
                good += 1
                continue
            if original_squared * ad ** 2 > (threshold + margin) ** 2:
                bad += 1
                continue
            coarse_uncertain += 1

            generators = tuple(
                tuple(eta[i] * v for v in direction)
                for direction in vectors if i not in anchors
            )
            if j not in anchors:
                generators += tuple(
                    tuple(eta[j] * v for v in direction)
                    for direction in vectors
                )

            if is_contact:
                vertices += 1
                corners = [tuple(2 * det * x for x in d)]
                for generator in generators:
                    corners = [
                        tuple(x + sign * generator[k] for k, x in enumerate(v))
                        for v in corners for sign in (-1, 1)
                    ]
                if all(
                    _dot(q, q) <= (2 * det * cutoff_mA) ** 2
                    for q in corners
                ):
                    good += 1
                    gained_contact += 1
                else:
                    uncertain += 1
            else:
                projections += 1
                projection_error = sum(abs(_dot(d, g)) for g in generators)
                minimum_projection = 2 * ad * original_squared - projection_error
                if (
                    minimum_projection > 0
                    and minimum_projection ** 2 >
                    (2 * ad * cutoff_mA) ** 2 * original_squared
                ):
                    bad += 1
                    gained_noncontact += 1
                else:
                    uncertain += 1

    if (good + bad + uncertain != total or
        gained_contact + gained_noncontact > coarse_uncertain):
        raise AssertionError("zonotope contact classification inconsistent")
    return ZonotopeContactCertificate(
        total, good, bad, uncertain, nominal, coarse_uncertain,
        gained_contact, gained_noncontact, vertices, projections,
    )
