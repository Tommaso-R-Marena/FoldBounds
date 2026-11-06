"""Independent regression and falsification for the exact zonotope refinement."""
from __future__ import annotations
import itertools
import random
import sys
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "data"))
from experimental_ca_packed import STRUCTURES, decode
from foldbounds.sparse_distance_certificate import (
    make_certificate, robust_contacts, _basis4, _sub, _dot,
)
from foldbounds.sparse_distance_zonotope import (
    _coefficient_vectors, zonotope_contact_bounds,
)


@pytest.mark.parametrize("name,old_undecided,new_undecided,yes,no,additional_yes,additional_no", [
    ("1A8O", 17, 4, 86, 1863, 8, 5),
    ("1LCD", 18, 7, 91, 1030, 4, 7),
    ("2BEG", 2, 0, 1, 252, 1, 1),
    ("7DDO", 229, 66, 1203, 174852, 80, 83),
])
def test_real_pdb_stronger_certificates(
    name, old_undecided, new_undecided, yes, no, additional_yes, additional_no
):
    p = decode(STRUCTURES[name]["packed"])
    original = robust_contacts(make_certificate(p), epsilon_mA=10)
    c = zonotope_contact_bounds(make_certificate(p), epsilon_mA=10)
    assert original.undecided_pairs == old_undecided
    assert c.coarse_undecided_pairs == old_undecided
    assert c.undecided_pairs == new_undecided
    assert (c.guaranteed_contacts, c.guaranteed_noncontacts) == (yes, no)
    assert (c.new_contacts, c.new_noncontacts) == (
        additional_yes, additional_no
    )
    assert c.guaranteed_contacts >= original.guaranteed_contacts
    assert c.guaranteed_noncontacts >= original.guaranteed_noncontacts
    assert c.all_pairs == (
        c.guaranteed_contacts + c.guaranteed_noncontacts + c.undecided_pairs
    )


def test_randomized_enclosures_never_worsen():
    rng = random.Random(1907)
    for _ in range(30):
        p = tuple(
            tuple(rng.randrange(-25000, 25000) for _ in range(3))
            for _ in range(13)
        )
        certificate = make_certificate(p)
        for epsilon in (0, 10, 50):
            prior = robust_contacts(certificate, epsilon_mA=epsilon)
            updated = zonotope_contact_bounds(certificate, epsilon_mA=epsilon)
            assert updated.undecided_pairs <= prior.undecided_pairs
            assert updated.guaranteed_contacts >= prior.guaranteed_contacts
            assert updated.guaranteed_noncontacts >= prior.guaranteed_noncontacts
            assert updated.coarse_undecided_pairs == prior.undecided_pairs


def test_exact_zero_error_and_invalid_settings():
    p = decode(STRUCTURES["2BEG"]["packed"])
    certificate = make_certificate(p)
    for kwargs in (
        {"epsilon_mA": -1},
        {"cutoff_mA": 0},
        {"epsilon_mA": True},
        {"min_sequence_separation": 0},
    ):
        with pytest.raises(ValueError):
            zonotope_contact_bounds(certificate, **kwargs)
    z = zonotope_contact_bounds(certificate, epsilon_mA=0)
    assert z.undecided_pairs == 0
    assert z.guaranteed_contacts == z.nominal_contacts


def test_independent_8d_vertex_generation():
    p = (
        (0, 0, 0), (1000, 0, 0), (0, 1000, 0), (0, 0, 1000),
        (3000, 4000, 5000), (4000, 4000, 4000),
    )
    c = make_certificate(p, (0, 1, 2, 3))
    _, cof, det, _ = _basis4(c.anchors)
    directions = _coefficient_vectors(cof)
    d = _sub(p[4], p[5])
    e4, e5 = 50, 100
    scores = []
    for signs in itertools.product((-1, 1), repeat=8):
        offset = tuple(
            sum(signs[k] * e4 * directions[k][j] for k in range(4))
            + sum(signs[k + 4] * e5 * directions[k][j] for k in range(4))
            for j in range(3)
        )
        q = tuple(2 * det * d[j] + offset[j] for j in range(3))
        scores.append(_dot(q, q))
    assert len(scores) == 256 and max(scores) >= min(scores)
