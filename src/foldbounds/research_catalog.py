"""Discovery-only registry of independent FoldBounds research capabilities.

Listing capabilities does NOT import optional heavy dependencies and does not
claim that distinct models form a calibrated end-to-end protein predictor.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass
from importlib.util import find_spec
import argparse
import json


@dataclass(frozen=True)
class Capability:
    name: str
    category: str
    module: str
    assurance: str
    purpose: str


CAPABILITIES = (
    Capability("HP lattice bounds", "core", "foldbounds.lookahead", "Python exact-rational checker / separate Lean certificate theory", "2D HP finite-model probability enclosures"),
    Capability("Contact-map atlas", "thermodynamics", "foldbounds.contact_map_atlas", "model-relative exact enumeration", "Reusable HP contact spectra and mutation queries"),
    Capability("Exact phase diagram", "thermodynamics", "foldbounds.multiway_envelope", "exact model-relative arithmetic", "Sequence-dependent all-q phase envelopes"),
    Capability("3D robust geometry", "geometry", "foldbounds.cubic_robust_inverse_design", "conditional geometric certificates", "Cubic lattice inverse design"),
    Capability("Kinetic plateau", "kinetics", "foldbounds.pivot_plateau_verifier", "independent model-relative checker", "Metastability under stipulated pivot moves"),
    Capability("GNM exact integer network", "structural", "foldbounds.gnm_integer_updates", "exact-integer checker", "Uncertain network response certificates"),
    Capability("GNM exact coordinate boxes", "structural", "foldbounds.gnm_exact_coordinate_boxes", "conditional exact-coordinate model", "Contact uncertainty and invariant ranks"),
    Capability("NMR restraint audit", "structural", "foldbounds.noe_moving_scaffold", "conditional rational geometric bounds", "NOE restraint feasibility under local atom-motion budgets"),
    Capability("Mutation-safe fold fingerprint", "structural", "foldbounds.ga_gb_mutation_safe", "source-frozen held-out NMR structural study", "Fold-family distance reporters"),
    Capability("S669 assay budget", "experimental decision", "foldbounds.anytime_assay_regret", "group-marginal conformal guarantees subject to exchangeability", "Adaptive mutation assay budgeting"),
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Discover experimental FoldBounds research capabilities")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable catalog")
    args = parser.parse_args()
    rows = [{**asdict(c), "module_present": find_spec(c.module) is not None} for c in CAPABILITIES]
    if args.json:
        print(json.dumps(rows, indent=2))
    else:
        for row in rows:
            print(f"{row['category']}: {row['name']} ({row['module']})")
    if not all(row["module_present"] for row in rows):
        parser.exit(2, "Missing integrated research module(s)\n")


if __name__ == "__main__":
    main()
