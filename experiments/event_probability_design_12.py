"""Reproduce exact contact-event probability design over all real q>1.

n=12 anchored 2D HP: 4096 sequences, k=0..12, events K>=1,...,5.
Computation uses exact rational/Sturm algebra and a full independent
geometry+landscape audit. Nothing here is an experimental protein predictor.
"""
from fractions import Fraction

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import audit_landscape, build_landscape
from foldbounds.event_probability_envelope import certify_event_probability_envelope


def main():
    atlas = build_atlas(12)
    if not audit_atlas(atlas):
        raise RuntimeError("independent walk re-enumeration failed")
    model = build_landscape(atlas)
    if not audit_landscape(model, atlas):
        raise RuntimeError("independent complete polynomial audit failed")

    cases = [
        certify_event_probability_envelope(model, h_count, threshold)
        for h_count in range(13) for threshold in range(1, 6)
    ]
    switches = sum(len(c.switches) for c in cases)
    classes = sum(len(set(c.root_free_winners)) for c in cases)
    assert len(cases) == 65 and switches == 74 and classes == 136

    print("FoldBounds exact contact-event global design; model 2D HP n=12")
    print("Independent full-walk atlas: 30,073 conformations; 4,096 H/P sequences")
    print("All-sequence polynomial/spectrum independent audit: PASSED")
    print("h_count event_at_least sequences pair_roots switches winner_classes")
    for c in cases:
        print(
            c.h_count, c.event_at_least, c.population, len(c.events),
            len(c.switches), len(set(c.root_free_winners))
        )
    print(f"RESULT: {len(cases)} global design problems; {switches} switches; "
          f"{classes} winning rational probability classes.")

    flagship = next(c for c in cases if c.h_count == 7 and c.event_at_least == 4)
    assert flagship.population == 792
    assert len(flagship.events) == 25 and len(flagship.switches) == 4
    print("\nFlagship: 792 seven-H sequences, maximize P(K>=4) for all q>1.")
    previous = None
    for winners in flagship.root_free_winners:
        if winners != previous:
            print("  global winning class:", ", ".join(winners))
            previous = winners
    print("Exactly four optimum switches with rational q brackets:")
    for i, event in enumerate(flagship.switches, 1):
        print(
            f"  {i}: [{event.t_lower + 1}, {event.t_upper + 1}] "
            f"(~{float((event.t_lower + event.t_upper) / 2 + 1):.9f})"
        )
    print("Important: q is dimensionless, not Kelvin.")
    print("Trust: independent Python exact audits, NOT Lean-refined/real-protein validated.")


if __name__ == "__main__":
    main()
