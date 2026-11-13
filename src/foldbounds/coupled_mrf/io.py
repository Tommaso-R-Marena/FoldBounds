"""JSON input/output for coupled factor-graph models.

Model file format (``"format": "foldbounds.coupled_mrf"``, ``"version": 1``)
---------------------------------------------------------------------------
::

    {
      "format": "foldbounds.coupled_mrf",
      "version": 1,
      "variables": [
        {"name": "x0", "domain": 2},                      # same domain in A and B
        {"name": "m",  "domain_A": 3, "domain_B": 5}      # domain differs (mutation)
      ],
      "shared": [                                         # identical in A and B
        {"name": "J_x0_x1", "scope": ["x0", "x1"], "table": ["2", "1", "1", "2"]}
      ],
      "only_A": [ {"scope": ["m"], "table": ["1", "1/2", "3"]} ],
      "only_B": [ {"scope": ["m"], "table": ["1", "1", "1", "2", "5/7"]} ],
      "metadata": {"free-form": "ignored by all computations and by the model hash"}
    }

* Numbers are exact rationals written as JSON strings ``"p"`` or ``"p/q"`` (or
  JSON integers).  JSON floats are rejected: a model must be exact.
* ``table`` is row-major over ``scope``: the first scope variable is the most
  significant digit, and the states of a variable with domain size d are
  ``0..d-1``.  Its length must equal the product of the scope's domain sizes
  in the system(s) the factor belongs to.
* Factor names are optional and do not affect semantics.

The *canonical* form (used for hashing and for binding certificates to a model)
lists variables in file order, writes every rational in lowest terms as
``"p"`` or ``"p/q"``, drops names and metadata, and is serialised with
``json.dumps(sort_keys=True, separators=(",", ":"))``.  The SHA-256 of that
string is the model hash.  The independent checker re-implements this
canonicalisation separately.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

from .model import CoupledModel, Factor, ModelError, as_rational

FORMAT = "foldbounds.coupled_mrf"
VERSION = 1


def frac_to_str(x: Fraction) -> str:
    x = Fraction(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def _parse_number(value: Any) -> Fraction:
    if isinstance(value, float):
        raise ModelError(f"JSON float {value!r} is not allowed; write rationals as strings")
    return as_rational(value)


def _factor_to_dict(f: Factor, canonical: bool) -> dict:
    d = {"scope": list(f.scope), "table": [frac_to_str(t) for t in f.table]}
    if f.name and not canonical:
        d["name"] = f.name
    return d


def model_to_dict(model: CoupledModel, canonical: bool = False) -> dict:
    variables = []
    for v in model.variables:
        da, db = model.domains_a[v], model.domains_b[v]
        if da == db:
            variables.append({"name": v, "domain": da})
        else:
            variables.append({"name": v, "domain_A": da, "domain_B": db})
    out = {
        "format": FORMAT,
        "version": VERSION,
        "variables": variables,
        "shared": [_factor_to_dict(f, canonical) for f in model.shared],
        "only_A": [_factor_to_dict(f, canonical) for f in model.only_a],
        "only_B": [_factor_to_dict(f, canonical) for f in model.only_b],
    }
    if not canonical and model.metadata:
        out["metadata"] = dict(model.metadata)
    return out


def _factor_from_dict(d: Any) -> Factor:
    if not isinstance(d, dict):
        raise ModelError("factor must be a JSON object")
    unknown = set(d) - {"scope", "table", "name"}
    if unknown:
        raise ModelError(f"unknown factor keys {sorted(unknown)}")
    scope = d.get("scope")
    table = d.get("table")
    if not isinstance(scope, list) or not all(isinstance(v, str) for v in scope):
        raise ModelError("factor scope must be a list of variable names")
    if not isinstance(table, list):
        raise ModelError("factor table must be a list")
    name = d.get("name", "")
    if not isinstance(name, str):
        raise ModelError("factor name must be a string")
    return Factor(tuple(scope), tuple(_parse_number(t) for t in table), name)


def model_from_dict(d: Any) -> CoupledModel:
    if not isinstance(d, dict):
        raise ModelError("model must be a JSON object")
    if d.get("format") != FORMAT or d.get("version") != VERSION:
        raise ModelError(f"expected format {FORMAT!r} version {VERSION}")
    unknown = set(d) - {"format", "version", "variables", "shared", "only_A", "only_B", "metadata"}
    if unknown:
        raise ModelError(f"unknown model keys {sorted(unknown)}")
    names, da, db = [], {}, {}
    for entry in d.get("variables", []):
        if not isinstance(entry, dict) or not isinstance(entry.get("name"), str):
            raise ModelError("variable entries need a string 'name'")
        name = entry["name"]
        keys = set(entry) - {"name"}
        if keys == {"domain"}:
            a = b = entry["domain"]
        elif keys == {"domain_A", "domain_B"}:
            a, b = entry["domain_A"], entry["domain_B"]
        else:
            raise ModelError(f"variable {name!r}: give 'domain' or both 'domain_A' and 'domain_B'")
        for size in (a, b):
            if isinstance(size, bool) or not isinstance(size, int) or size < 1:
                raise ModelError(f"variable {name!r}: invalid domain size {size!r}")
        names.append(name)
        da[name], db[name] = a, b
    if len(set(names)) != len(names):
        raise ModelError("repeated variable name")
    lists = []
    for key in ("shared", "only_A", "only_B"):
        raw = d.get(key, [])
        if not isinstance(raw, list):
            raise ModelError(f"{key} must be a list")
        lists.append(tuple(_factor_from_dict(f) for f in raw))
    meta = d.get("metadata", {})
    if not isinstance(meta, dict):
        raise ModelError("metadata must be an object")
    return CoupledModel(tuple(names), da, db, lists[0], lists[1], lists[2], meta)


def canonical_model_json(model: CoupledModel) -> str:
    return json.dumps(model_to_dict(model, canonical=True), sort_keys=True, separators=(",", ":"))


def model_sha256(model: CoupledModel) -> str:
    return hashlib.sha256(canonical_model_json(model).encode("utf-8")).hexdigest()


def dump_model(model: CoupledModel, path: str | Path) -> None:
    Path(path).write_text(json.dumps(model_to_dict(model), indent=1) + "\n", encoding="utf-8")


def load_model(path: str | Path) -> CoupledModel:
    def no_float(text: str):
        raise ModelError(f"JSON float {text!r} is not allowed")
    return model_from_dict(json.loads(Path(path).read_text(encoding="utf-8"), parse_float=no_float))


def loads_model(text: str) -> CoupledModel:
    def no_float(t: str):
        raise ModelError(f"JSON float {t!r} is not allowed")
    return model_from_dict(json.loads(text, parse_float=no_float))
