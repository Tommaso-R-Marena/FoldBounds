"""Independent assurance checks for the GA/GB frozen-JSON replay gate."""
import math
import pytest

from scripts.verify_ga_gb_frozen import compare_json


def test_only_tiny_cross_platform_float_variability_is_allowed():
    frozen = {"source_sha": "abc", "count": 20, "model": "GB95",
              "classes": ["parallel", "other"], "mean_distance_A": 8.12345678912345}
    actual = dict(frozen, mean_distance_A=frozen["mean_distance_A"] + 1e-13)
    errors, rounded, delta = compare_json(actual, frozen)
    assert errors == []
    assert rounded == 1 and 0 < delta < 1e-12


@pytest.mark.parametrize("new", [
    {"source_sha": "tampered"},
    {"count": 19},
    {"model": "GA95"},
    {"classes": ["other", "parallel"]},
    {"classes": ["parallel", "other", "other"]},
    {"mean_distance_A": 8.12345679912345},
    {"unexpected": "field"},
])
def test_wrong_provenance_semantics_counts_or_scientific_value_are_rejected(new):
    frozen = {"source_sha": "abc", "count": 20, "model": "GB95",
              "classes": ["parallel", "other"], "mean_distance_A": 8.12345678912345}
    actual = dict(frozen, **new)
    assert compare_json(actual, frozen)[0]


def test_boolean_vs_integer_and_null_vs_missing_not_conflated():
    assert compare_json({"x": True}, {"x": 1})[0]
    assert compare_json({"x": None}, {})[0]


def test_nonfinite_values_fail_closed():
    assert compare_json({"v": math.nan}, {"v": 1.5})[0]
    assert compare_json({"v": 1.5}, {"v": math.inf})[0]
