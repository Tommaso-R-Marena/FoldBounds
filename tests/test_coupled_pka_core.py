import random
from fractions import Fraction
import itertools

from foldbounds.coupled_pka import model as M, exact as ex, checker as C
from foldbounds.coupled_pka.treewidth import min_fill_width


def synth_sites(n=6, seed=1, box=9.0):
    rnd = random.Random(seed)
    types = ["ASP", "LYS", "GLU", "HIS", "TYR", "ARG", "C-term", "N-term"]
    out = []
    for k in range(n):
        t = types[k % len(types)]
        out.append({"chain": "A", "resnum": k + 1, "icode": "", "rtype": t,
                    "xyz": [rnd.uniform(0, box) for _ in range(3)], "pka": M.MODEL_PKA[t]})
    return out


def synth_spec(n=6, seed=1, box=9.0):
    return M.build_spec(synth_sites(n, seed, box), "synth", rc=12.0)


def test_pow10_rational_definition():
    assert M.pow10_rational(0) == 1
    assert M.pow10_rational(1000) == 10
    assert M.pow10_rational(-2000) == Fraction(1, 100)
    x = M.pow10_rational(500)
    assert abs(float(x) - 10 ** 0.5) < 1e-14
    assert 10 ** 20 % x.denominator == 0  # 15 significant decimal digits


def test_coupling_scale_and_sign():
    assert M.coupling_milli(4.0) == round(1000 * 83.01593 / 16 / M.kt_ln10()) or abs(M.coupling_milli(4.0) - 1000 * 83.01593 / 16 / M.kt_ln10()) < 2
    assert M.coupling_milli(3.0) > M.coupling_milli(6.0) > M.coupling_milli(12.0) > 0


def test_model_is_antiferromagnetic():
    spec = synth_spec()
    mod = M.Model(spec, Fraction(7))
    assert spec["pairs_milli"]
    for key in spec["pairs_milli"]:
        i, j = map(int, key.split(","))
        cross = mod.phi(i, j, 1, 1) * mod.phi(i, j, 0, 0) / (mod.phi(i, j, 1, 0) * mod.phi(i, j, 0, 1))
        assert cross < 1  # coefficient of sigma_i sigma_j in log weight is negative-energy: -w < 0 in weight, i.e. +w in energy


def naive_Q(spec, ph):
    mod = M.Model(spec, Fraction(ph))
    n = mod.n
    Z1 = [Fraction(0)] * n
    Z0 = [Fraction(0)] * n
    for s in itertools.product((0, 1), repeat=n):
        w = Fraction(1)
        for i in range(n):
            w *= mod.u(i, s[i])
            for j in range(i + 1, n):
                w *= mod.phi(i, j, s[i], s[j])
        for i in range(n):
            (Z1 if s[i] else Z0)[i] += w
    return [Z1[i] / Z0[i] for i in range(n)]


def test_exact_enumeration_matches_naive():
    spec = synth_spec(6, 3)
    for ph in (2, 7, 11):
        assert ex.exact_marginals(spec, Fraction(ph)) == naive_Q(spec, ph)


def test_treewidth():
    assert min_fill_width(5, [(i, i + 1) for i in range(4)])[0] == 1
    assert min_fill_width(5, [(i, (i + 1) % 5) for i in range(5)])[0] == 2
    assert min_fill_width(5, [(i, j) for i in range(5) for j in range(i + 1, 5)])[0] == 4


def test_checker_rejects_wrong_spec_hash():
    spec = synth_spec(4, 2)
    cert = {"spec_hash": "0" * 64, "site": 0, "interior": [0], "ph": "7", "side": "lower", "lambda": "0",
            "order": [], "tree": 0}
    try:
        C.check(spec, cert)
        assert False
    except C.CertificateError:
        pass


def test_checker_reports_malformed_numbers_as_certificate_error():
    """Adversarial review P2-m2: a 'nan' lambda must raise CertificateError, not ValueError."""
    import pytest
    from foldbounds.coupled_pka import checker
    spec = {"sites": []}
    cert = {"spec_hash": checker.spec_hash(spec), "ph": "7", "lambda": "nan", "side": "lower"}
    with pytest.raises(checker.CertificateError):
        checker.check(spec, cert)
    with pytest.raises(checker.CertificateError):
        checker.check(spec, {"ph": "7"})  # missing keys
