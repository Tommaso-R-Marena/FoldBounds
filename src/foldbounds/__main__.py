"""Reproducible command-line demonstration of finite-model bounds."""
import argparse
import json
from fractions import Fraction

from .lattice import build_frontier, exact_spectrum, verify_certificate
from .lookahead import verify_lookahead
from .adaptive import adaptive_refine
from .shell_portfolio import optimize_shell_portfolio
from .sparse_monotone import optimize_sparse_monotone
from .lattice import ROOT


def main() -> None:
    parser = argparse.ArgumentParser(description="FoldBounds 2D HP-model reference runner")
    parser.add_argument("sequence", help="H/P residue string, at least two residues")
    parser.add_argument("--q", type=Fraction, default=Fraction(2), help="exact q>=1, e.g. 3/2")
    parser.add_argument("--event-at-least", type=int, default=1, help="minimum contact count")
    parser.add_argument("--depth", type=int, default=2, help="certificate frontier depth")
    parser.add_argument("--enhanced-lookahead", type=int, default=None,
                        help="use witness/geometry bounds with N-step lookahead (N>=0)")
    parser.add_argument("--matching", action="store_true",
                        help="enable max-flow H-H contact cap (requires --enhanced-lookahead)")
    parser.add_argument("--transfer-memory", type=int, default=2,
                        help="recent visited sites kept by transfer DP (2-8)")
    parser.add_argument("--transfer", action="store_true",
                        help="enable prefix-obstacle transfer upper bounds")
    parser.add_argument("--terminal-box", action="store_true",
                        help="condition transfer energy caps on each endpoint's enclosing box")
    parser.add_argument("--monotone", action="store_true",
                        help="count disjoint monotone completion families for lower mass")
    parser.add_argument("--turn-family", action="store_true",
                        help="count checked self-avoiding suffixes with at most two turns")
    parser.add_argument("--shell-core", type=int, default=None,
                        help="exact lower mass from disjoint extremal monotone shells rooted at this core length")
    parser.add_argument("--portfolio-leaves", type=int, default=None,
                        help="max number of leaves in each optimized verified shell-portfolio")
    parser.add_argument("--portfolio-max-depth", type=int, default=None,
                        help="maximum core depth explored by exact portfolio tree-knapsack")
    parser.add_argument("--sparse-leaves", type=int, default=None,
                        help="maximum selected leaves in sparse lower-mass antichain")
    parser.add_argument("--sparse-max-depth", type=int, default=None,
                        help="maximum core depth for sparse monotone-family optimizer")
    parser.add_argument("--sparse-turn-union", action="store_true",
                        help="add exactly counted nonmonotone at-most-two-turn families to sparse certificates")
    parser.add_argument("--coupled-polytope", action="store_true",
                        help="optimize joint event, complement, and total mass intervals exactly")
    parser.add_argument("--witness-nodes", type=int, default=0,
                        help="DFS node budget per bounded subtree for checked lower mass")
    parser.add_argument("--adaptive-steps", type=int, default=None,
                        help="number of independently checked frontier expansions")
    parser.add_argument("--adaptive-policy", choices=("breadth","mass","probability"),
                        default="probability", help="frontier selection strategy")
    parser.add_argument("--exact-spectrum", action="store_true", help="exponentially expensive enumeration")
    args = parser.parse_args()
    if args.matching and args.enhanced_lookahead is None and args.adaptive_steps is None:
        parser.error('--matching requires --enhanced-lookahead')
    if not 2 <= args.transfer_memory <= 8:
        parser.error('--transfer-memory must be between 2 and 8')
    if args.transfer_memory != 2 and not args.transfer:
        parser.error('--transfer-memory above 2 requires --transfer')
    if args.transfer and args.enhanced_lookahead is None and args.adaptive_steps is None:
        parser.error('--transfer requires --enhanced-lookahead or --adaptive-steps')
    if args.terminal_box and not args.transfer:
        parser.error('--terminal-box requires --transfer')
    if args.monotone and args.enhanced_lookahead is None and args.adaptive_steps is None:
        parser.error('--monotone requires --enhanced-lookahead or --adaptive-steps')
    if args.turn_family and args.enhanced_lookahead is None and args.adaptive_steps is None:
        parser.error('--turn-family requires --enhanced-lookahead or --adaptive-steps')
    if args.shell_core is not None and not 2 <= args.shell_core <= len(args.sequence):
        parser.error('--shell-core must lie between 2 and the sequence length')
    if args.shell_core is not None and args.enhanced_lookahead is None and args.adaptive_steps is None:
        parser.error('--shell-core requires --enhanced-lookahead or --adaptive-steps')
    if (args.portfolio_leaves is None) != (args.portfolio_max_depth is None):
        parser.error('--portfolio-leaves and --portfolio-max-depth must be specified together')
    if args.portfolio_leaves is not None:
        if args.portfolio_leaves < 1 or not 2 <= args.portfolio_max_depth <= len(args.sequence):
            parser.error('invalid portfolio leaf budget or maximum depth')
        if args.enhanced_lookahead is None and args.adaptive_steps is None:
            parser.error('portfolio mode requires --enhanced-lookahead or --adaptive-steps')
    if args.sparse_turn_union and args.sparse_leaves is None:
        parser.error('--sparse-turn-union requires --sparse-leaves')
    if (args.sparse_leaves is None) != (args.sparse_max_depth is None):
        parser.error('--sparse-leaves and --sparse-max-depth must be specified together')
    if args.sparse_leaves is not None:
        if args.sparse_leaves < 1 or not 2 <= args.sparse_max_depth <= len(args.sequence):
            parser.error('invalid sparse certificate budget or maximum depth')
        if args.enhanced_lookahead is None and args.adaptive_steps is None:
            parser.error('sparse certificates require --enhanced-lookahead or --adaptive-steps')
    if args.witness_nodes and args.enhanced_lookahead is None and args.adaptive_steps is None:
        parser.error('--witness-nodes requires --enhanced-lookahead or --adaptive-steps')
    cert = build_frontier(args.sequence, args.q, args.event_at_least, args.depth)
    portfolio_plans = ()
    if args.portfolio_leaves is not None:
        portfolio_plans = tuple(optimize_shell_portfolio(
            args.sequence, ROOT, args.q, args.event_at_least,
            leaf_budget=args.portfolio_leaves,
            max_core_length=args.portfolio_max_depth,
            objective=objective) for objective in ('event', 'complement'))
    portfolio_certificates = tuple(p.certificate for p in portfolio_plans)
    sparse_plans=()
    if args.sparse_leaves is not None:
        sparse_plans=tuple(optimize_sparse_monotone(
            args.sequence,ROOT,args.q,args.event_at_least,
            leaf_budget=args.sparse_leaves,max_core_length=args.sparse_max_depth,
            objective=obj,include_nonmonotone_turns=args.sparse_turn_union)
            for obj in ('event','complement'))
    sparse_certificates=tuple(p.certificate for p in sparse_plans)
    run = (adaptive_refine(cert, steps=args.adaptive_steps,
                           lookahead=args.enhanced_lookahead or 0,
                           use_matching=args.matching,use_transfer=args.transfer,
                           witness_nodes=args.witness_nodes,
                           use_monotone=args.monotone,
                           use_turn_family=args.turn_family,
                           terminal_box=args.terminal_box,
                           shell_core=args.shell_core,
                           shell_portfolios=portfolio_certificates,
                           sparse_certificates=sparse_certificates,
                           use_coupled_polytope=args.coupled_polytope,
                           transfer_memory=args.transfer_memory,
                           policy=args.adaptive_policy)
           if args.adaptive_steps is not None else None)
    bounds = (run.bounds if run is not None else
              verify_certificate(cert) if args.enhanced_lookahead is None else
              verify_lookahead(cert,lookahead=args.enhanced_lookahead,
                               use_matching=args.matching,use_transfer=args.transfer,
                               witness_nodes=args.witness_nodes,
                               use_monotone=args.monotone,
                               use_turn_family=args.turn_family,
                               terminal_box=args.terminal_box,
                               shell_core=args.shell_core,
                               shell_portfolios=portfolio_certificates,
                               sparse_certificates=sparse_certificates,
                               use_coupled_polytope=args.coupled_polytope,
                               transfer_memory=args.transfer_memory))
    out = {
        "model": "2D square-lattice HP self-avoiding walk, anchored first edge",
        "method": "baseline" if args.enhanced_lookahead is None else "geometric-witness-lookahead",
        "enhanced_lookahead": args.enhanced_lookahead,
        "matching": args.matching,
        "transfer": args.transfer,
        "transfer_memory": args.transfer_memory,
        "witness_nodes": args.witness_nodes,
        "monotone": args.monotone,
        "turn_family": args.turn_family,
        "shell_core": args.shell_core,
        "portfolio_leaves": args.portfolio_leaves,
        "coupled_polytope": args.coupled_polytope,
        "sparse_leaves": args.sparse_leaves,
        "sparse_turn_union": args.sparse_turn_union,
        "sparse_max_depth": args.sparse_max_depth,
        "sparse_certificates": [dict(objective=p.objective,leaf_count=p.mass.leaf_count,
                                      candidate_prefixes=p.candidate_prefixes,
                                      event_mass=str(p.mass.event),
                                      complement_mass=str(p.mass.complement))
                                for p in sparse_plans],
        "portfolio_max_depth": args.portfolio_max_depth,
        "portfolio_certificates": [dict(objective=p.objective,leaf_count=p.mass.leaf_count,
                                        candidate_prefixes=p.candidate_nodes,
                                        selected_event_mass=str(p.mass.event),
                                        selected_complement_mass=str(p.mass.complement))
                                   for p in portfolio_plans],
        "terminal_box": args.terminal_box,
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
    if run is not None:
        out["adaptive"] = {
            "steps_executed": len(run.history)-1,
            "policy": run.policy,
            "history": [{
                "step": h.step,
                "chosen": h.chosen,
                "bounded_nodes": h.bounded_nodes,
                "partition_interval": [str(h.partition_lower),str(h.partition_upper)],
                "probability_interval": [str(h.probability_lower),str(h.probability_upper)]
            } for h in run.history]
        }
    if args.exact_spectrum:
        out["contact_spectrum"] = exact_spectrum(args.sequence)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
