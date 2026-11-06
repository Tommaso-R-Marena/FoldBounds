"""Reproduce all 17 exact length-16 fixed-composition HP phase diagrams.

The deterministic demonstration enumerates the complete finite lattice
geometry, independently replays the walk atlas, builds exact polynomial
tables, then certifies global optima for all real q>1 by Sturm isolation.

This is a finite toy-model result. Not experimental protein prediction.
"""
from fractions import Fraction

from foldbounds.contact_map_atlas import build_atlas, audit_atlas
from foldbounds.universal_landscape import build_landscape
from foldbounds.multiway_envelope import certify_multiway_envelope


def main():
    atlas = build_atlas(16)
    if not audit_atlas(atlas, max_length=16):
        raise RuntimeError("independent complete-walk audit failed")
    model = build_landscape(
        atlas, require_audit=False, max_length=16,
        max_edge_subsets=3_000_000
    )
    envelopes = tuple(
        certify_multiway_envelope(model, k, width=Fraction(1, 1_000_000))
        for k in range(17)
    )
    total_sequences = sum(c.population for c in envelopes)
    total_switches = sum(len(c.switches) for c in envelopes)
    total_classes = sum(len(set(c.root_free_winners)) for c in envelopes)
    assert (total_sequences, total_switches, total_classes) == (65_536, 25, 42)

    print("FoldBounds: exact all-q 2D HP sequence optimization")
    print("Geometry: n=16, walks=1,604,149, contact maps=31,940")
    print("Full independent geometry replay: PASSED")
    print("k sequences skyline-class final-class pair-roots global-switches winner-classes")
    for c in envelopes:
        print(
            c.h_count, c.population, c.skyline_classes, len(c.final_classes),
            len(c.events), len(c.switches), len(set(c.root_free_winners))
        )
    print(
        f"TOTAL: {total_sequences} sequences, {total_switches} global switches, "
        f"{total_classes} optimal polynomial classes."
    )
    selected = envelopes[9]
    assert selected.population == 11_440 and len(selected.switches) == 4
    print("\nExact nine-H global optimality phase diagram:")
    print("Globally optimal sequence classes by root-free q-region:")
    previous = None
    for winners in selected.root_free_winners:
        if winners != previous:
            print("  " + ", ".join(winners))
            previous = winners
    print("Exactly four global transition roots (q=1+t):")
    for i, event in enumerate(selected.switches, 1):
        assert event.kind == "global_switch"
        print(
            f"  {i}: [{event.t_lower + 1}, {event.t_upper + 1}] "
            f" (~{float((event.t_lower + event.t_upper)/2 + 1):.9f})"
        )
    print(
        "Other pairwise roots are certified irrelevant to the global maximum. "
        "At q=1 all fixed-length sequences have the same partition count."
    )
    print(
        "Trust: exact rational Python, geometry separately re-enumerated; "
        "NOT Lean-refined, NOT experimental protein science."
    )


if __name__ == "__main__":
    main()
