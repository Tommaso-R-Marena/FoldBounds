import copy
import hashlib
import json
from fractions import Fraction as F

import pytest

from foldbounds.coupled_mrf import checker
from foldbounds.coupled_mrf.certificate import certificate_dict, exhaustive_enclosure
from foldbounds.coupled_mrf.generators import ising_grid, random_sparse_binary, rotamer_potts
from foldbounds.coupled_mrf.io import model_to_dict
from foldbounds.coupled_mrf.model import Factor, make_model
from foldbounds.coupled_mrf.separator import Separator, radius_separator


def make(m, R):
    sep = radius_separator(m, R)
    enc = exhaustive_enclosure(m, sep)
    cert = certificate_dict(m, sep, enc)
    return json.loads(json.dumps(model_to_dict(m))), json.loads(json.dumps(cert))


@pytest.fixture(scope="module")
def good():
    m = ising_grid(5, F(3, 2), F(6, 5), ("field", (2, 2), F(1, 2)))
    return make(m, 2)


def test_sha_matches_hashlib():
    for s in (b"", b"abc", b"a" * 55, b"a" * 56, b"a" * 64, bytes(range(256)) * 3):
        assert checker.sha256_hex(s) == hashlib.sha256(s).hexdigest()


def test_checker_imports_only_fractions_and_json():
    import ast, pathlib
    src = pathlib.Path(checker.__file__).read_text()
    mods = set()
    for n in ast.walk(ast.parse(src)):
        if isinstance(n, ast.Import):
            mods |= {a.name for a in n.names}
        if isinstance(n, ast.ImportFrom):
            mods.add(n.module)
    assert mods <= {"json", "fractions"}


def test_accepts_valid(good):
    md, cert = good
    ok, why = checker.check_certificate(md, cert)
    assert ok, why


@pytest.mark.parametrize("seed", range(4))
def test_accepts_random_models(seed):
    for R in (0, 1, 2):
        ok, why = checker.check_certificate(*make(random_sparse_binary(10, 3, seed), R))
        assert ok, why


def test_accepts_domain_change_and_unbounded():
    m = rotamer_potts(6, seed=1, rotamers=(2, 3), mutant_rotamers=4, mutated=2)
    ok, why = checker.check_certificate(*make(m, 1))
    assert ok, why
    dom = {"i": 2, "b": 2}
    mm = make_model(dom, shared=[Factor(("i", "b"), (F(1), F(0), F(1), F(1)))],
                    only_a=[Factor(("i",), (F(1), F(0)))])
    sep = Separator.of(["i"], ["b"], [])
    enc = exhaustive_enclosure(mm, sep)
    cert = json.loads(json.dumps(certificate_dict(mm, sep, enc)))
    assert cert["claim"]["upper"] == "inf"
    ok, why = checker.check_certificate(model_to_dict(mm), cert)
    assert ok, why


def mutate(good, fn):
    md, cert = copy.deepcopy(good[0]), copy.deepcopy(good[1])
    r = fn(md, cert)
    return (r if r is not None else md), cert


def rejected(good, fn, **kw):
    md, cert = mutate(good, fn)
    ok, why = checker.check_certificate(md, cert, **kw)
    assert not ok, "corrupted certificate accepted"
    return why


def test_reject_lower_too_high(good):
    def f(md, c):
        lo = F(c["claim"]["lower"]); c["claim"]["lower"] = str(lo + F(1, 1000))
    rejected(good, f)


def test_reject_upper_too_low(good):
    def f(md, c):
        up = F(c["claim"]["upper"]); c["claim"]["upper"] = str(up - F(1, 10 ** 6))
    rejected(good, f)


def test_reject_tight_but_loose(good):
    def f(md, c):
        c["claim"]["upper"] = str(F(c["claim"]["upper"]) + 1)
    rejected(good, f)


def test_reject_upper_inf_when_tight(good):
    def f(md, c):
        c["claim"]["upper"] = "inf"
    rejected(good, f)


def test_reject_model_changed(good):
    def f(md, c):
        md["shared"][0]["table"][0] = "7/4"
    rejected(good, f)


def test_reject_hash_wrong(good):
    def f(md, c):
        c["model_sha256"] = "0" * 64
    rejected(good, f)


def test_reject_expected_hash_mismatch(good):
    md, cert = good
    ok, _ = checker.check_certificate(md, cert, expected_sha256="1" * 64)
    assert not ok


def test_reject_separator_not_partition(good):
    def f(md, c):
        c["separator"]["exterior"] = c["separator"]["exterior"][1:]
    rejected(good, f)


def test_reject_factor_touches_interior_and_exterior(good):
    def f(md, c):
        s = c["separator"]
        v = s["boundary"].pop(0)
        s["exterior"].append(v)
    why = rejected(good, f)
    assert any("interior and exterior" in w for w in why)


def test_reject_specific_factor_outside_local(good):
    def f(md, c):
        # shrink to a separator not containing the perturbed variable locally
        c["separator"] = {"interior": ["x0_0"], "boundary": ["x0_1", "x1_0"],
                          "exterior": [v["name"] for v in md["variables"]
                                       if v["name"] not in ("x0_0", "x0_1", "x1_0")]}
    rejected(good, f)


def test_reject_witness_wrong(good):
    def f(md, c):
        w = c["witness"]["lower"]; w[0] = 1 - w[0]
        c["witness"]["lower"] = w
    md, cert = mutate(good, f)
    ok, why = checker.check_certificate(md, cert)
    # a flipped witness is either the same value by symmetry or rejected; make it surely wrong
    cert["witness"]["lower"] = cert["witness"]["upper"]
    ok, why = checker.check_certificate(md, cert)
    assert not ok


def test_reject_float_and_nan_numbers(good):
    def f(md, c):
        c["claim"]["lower"] = 0.5
    rejected(good, f)
    def g(md, c):
        c["claim"]["lower"] = "nan"
    rejected(good, g)
    def h(md, c):
        md["shared"][0]["table"][1] = 0.5
    rejected(good, h)


def test_reject_negative_table_entry(good):
    def f(md, c):
        md["shared"][0]["table"][1] = "-1"
    rejected(good, f)


def test_reject_bad_format_version_unknown_keys(good):
    for key, val in (("format", "x"), ("version", 2), ("evil", 1)):
        def f(md, c, key=key, val=val):
            c[key] = val
        rejected(good, f)


def test_reject_garbage_inputs():
    for md, cert in ((None, None), ({}, {}), ([], []), ("x", {"format": 1})):
        ok, why = checker.check_certificate(md, cert)
        assert not ok


def test_reject_claim_on_empty_support():
    dom = {"i": 2, "b": 2}
    mm = make_model(dom, shared=[Factor(("i", "b"), (F(0),) * 4)])
    sep = Separator.of(["i"], ["b"], [])
    enc = exhaustive_enclosure(mm, sep)
    cert = json.loads(json.dumps(certificate_dict(mm, sep, enc)))
    ok, _ = checker.check_certificate(model_to_dict(mm), cert)
    assert ok  # no claim made
    cert["claim"]["lower"] = "0"; cert["claim"]["upper"] = "inf"
    ok, _ = checker.check_certificate(model_to_dict(mm), cert)
    assert not ok


def test_boundary_limit_fails_closed(good):
    md, cert = good
    ok, why = checker.check_certificate(md, cert, max_boundary_states=2)
    assert not ok


def test_checker_agrees_with_producer_on_loose_valid_claim(good):
    md, cert = copy.deepcopy(good[0]), copy.deepcopy(good[1])
    cert["claim"]["tight"] = False
    cert["claim"]["lower"] = "0"
    cert["claim"]["upper"] = "inf"
    ok, why = checker.check_certificate(md, cert)
    assert ok, why
