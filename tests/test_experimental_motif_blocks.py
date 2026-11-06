"""Independent and adversarial tests of motif-blocked physical-data evaluation."""
from decimal import Decimal as D
from fractions import Fraction
from random import Random

import pytest

from foldbounds.experimental_reciprocal_conformal import (
    MODELS, AssuranceError, Measurement,
)
from foldbounds.experimental_motif_blocks import (
    motif_components, motif_blocks, evaluate_motif_blocked, summarize_motif,
)


def examples():
    rng=Random(20261105)
    amino="ACDEFGHIKLMNPQRSTVWY"
    sequences=["".join(rng.choice(amino) for _ in range(30)) for _ in range(25)]
    assert len(set(sequences))==25
    output=[]
    for i,s in enumerate(sequences):
        for j in range(1+i%5):
            estimates={}
            for name in MODELS:
                estimates[name+"_dir"]=D(i+j)/5
                estimates[name+"_inv"]=-D(i+j)/5
            mutant=("A" if s[0]!="A" else "H")+s[1:]
            output.append(Measurement(f"{i}_{j}",s,mutant,D("7"),D("25"),
                                      D(i-j)/3,estimates))
    return output


def test_transitive_component_closure():
    a="AAAAAAAABBBBCCCC"
    b="BBBBCCCCDDDD"
    c="CCCCDDDDFFFF"
    d="XYZXYZXYZXYZXYZ"
    e="NNNNNNNNNNNNNNN"
    connected=motif_components((a,b,c,d,e),k=8)
    assert connected[a]==connected[b]==connected[c]
    assert connected[d]!=connected[a] and connected[e]!=connected[d]
    assert len(set(connected.values()))==3


def test_sequence_aliases_always_share_component():
    seq="ACDEFGHIKLMNPQRSTVWY"
    components=motif_components((seq,seq),k=12)
    assert len(set(components.values()))==1


def test_no_motif_leakage_and_sequence_level_partition_completeness():
    records=examples()
    assignment,groups,reps=motif_blocks(records,k=12)
    assert 5<=len(groups)<=25
    assert sum(map(len,groups.values()))==len(records)
    for f in range(5):
        test={g for g in groups if assignment[g]==f}
        cal={g for g in groups if assignment[g]==(f+1)%5}
        train=set(groups)-test-cal
        assert not test&cal and not test&train and not cal&train
        assert len(test|cal|train)==len(groups)
        assert all(reps[x.wt]==g for g in groups for x in groups[g])


@pytest.mark.parametrize("k",[6,10,12,16,20])
def test_multiple_motif_block_lengths(k):
    records=examples()
    result=evaluate_motif_blocked(records,k=k)
    summary=summarize_motif(result)
    assert summary["mutations_total"]==len(records)
    assert sum(x["test_components"] for x in result)==summary["components_total"]


@pytest.mark.parametrize("bad",[0,1,5,41,12.0,True,"12"])
def test_invalid_k_fails_closed(bad):
    with pytest.raises(AssuranceError):
        motif_components(["ACDEFGHIKLMNPQRSTVWY"],k=bad)


def test_blockwise_guarantee_is_stricter_than_rowwise_negative_control():
    records=examples()
    blocked=evaluate_motif_blocked(records,calibration="component",alpha=Fraction(1,5))
    naive=evaluate_motif_blocked(records,calibration="mutation",alpha=Fraction(1,5))
    for strong,weak in zip(blocked,naive):
        assert strong["test_mutations"]==weak["test_mutations"]
        assert strong["test_components"]==weak["test_components"]
        if strong["cutoff"] is not None and weak["cutoff"] is not None:
            assert D(weak["cutoff"])<=D(strong["cutoff"])
        assert weak["components_covered"]<=strong["components_covered"]


def test_order_independence_and_train_only_physical_model_selection():
    records=examples()
    left=evaluate_motif_blocked(records,strategy="selected_phys")
    right=evaluate_motif_blocked(list(reversed(records)),strategy="selected_phys")
    assert left==right
    assert all(x["predictor"].endswith("_sym") for x in left)


def test_invalid_strategy_and_calibration_rejected():
    with pytest.raises(AssuranceError):
        evaluate_motif_blocked(examples(),strategy="chosen_after_test")
    with pytest.raises(AssuranceError):
        evaluate_motif_blocked(examples(),calibration="individual-with-family-guarantee")


def test_not_claiming_the_cluster_represents_a_biological_family():
    # Some unrelated chains have the same affinity tag; shared motif
    # clustering deliberately merges them, without inferring ancestry.
    tag="HHHHHHHHHHHH"
    a=tag+"ACDEFGHIKLMNPQRSTV"
    b=tag+"VVVVTTTTGGGGAAAARR"
    c="MNPQRSTVWWACDEFGHIK"
    d=motif_components([a,b,c],k=12)
    assert d[a]==d[b]
    assert len(set(d.values()))>=1
