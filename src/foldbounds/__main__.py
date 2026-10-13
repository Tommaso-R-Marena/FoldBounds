"""Reproducible command-line demonstration of finite-model bounds."""
import argparse
import json
from fractions import Fraction

from .lattice import build_frontier, exact_spectrum, verify_certificate


def main() -> None:
    parser = argparse.ArgumentParser(description="FoldBounds 2D HP-model reference runner")
    parser.add_argument("sequence", help="H/P residue string, at least two residues")
    parser.add_argument("--q", type=Fraction, default=Fraction(2), help="exact q>=1, e.g. 3/2")
    parser.add_argument("--event-at-least", type=int, default=1, help="minimum contact count")
    parser.add_argument("--depth", type=int, default=2, help="certificate frontier depth")
    parser.add_argument("--exact-spectrum", action="store_true", help="exponentially expensive enumeration")
    args = parser.parse_args()
    cert = build_frontier(args.sequence, args.q, args.event_at_least, args.depth)
    bounds = verify_certificate(cert)
    out = {
        "model": "2D square-lattice HP self-avoiding walk, anchored first edge",
        "sequence": args.sequence,
        "q": str(args.q),
        "event_at_least": args.event_at_least,
        "depth": args.depth,
        "partition_interval": [str(bounds.partition_lower), str(bounds.partition_upper)],
        "event_probability_interval": [str(bounds.probability_lower), str(bounds.probability_upper)],
        "exact_paths": bounds.exact_paths,
        "bounded_subtrees": bounds.bounded_subtrees,
        "exact": bounds.exact,
    }
    if args.exact_spectrum:
        out["contact_spectrum"] = exact_spectrum(args.sequence)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
