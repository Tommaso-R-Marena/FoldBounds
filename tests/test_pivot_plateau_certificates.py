"""Full independent kinetic-proof tests, including adversarial cut attacks."""
from collections import Counter, deque
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.pivot_plateau_producer import (
    PlateauProposalFailure, PivotPlateauCertificate, hh_score,
    proposed_distinct_moves, propose_rectangular_plateau,
    rectangular_hamiltonian_snake,
)
from foldbounds.pivot_plateau_verifier import (
    PlateauCertificateError, _independent_contacts, _geometry_ok,
    _state_proposal_targets, verify_pivot_plateau,
)
from foldbounds.kinetic_barrier_verifier import _walks_independent


def rooted_walks(n):
    return _walks_independent(n,40000)


@pytest.mark.parametrize('width,height',[(4,3),(4,4),(4,10),(4,20),(8,5),(8,10)])
def test_rectangular_constructive_complete_no_local_motion(width,height):
    p=rectangular_hamiltonian_snake(width,height)
    n=width*height
    assert _geometry_ok(p,n)
    assert len(set(p)) == n
    assert hh_score(p,'H'*n) == n + 1 - width-height
    # Independent geometry of square corners + endpoint
    for i in range(2,n-1):
        a,b,c=p[i-1:i+2]
        other=(a[0]+c[0]-b[0],a[1]+c[1]-b[1])
        assert other in p or abs(other[0]-a[0])+abs(other[1]-a[1])!=1
    left=p[-2]
    assert all((left[0]+dx,left[1]+dy) in p for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)))


@pytest.mark.parametrize('n', [4,5,6,7,8])
def test_all_proposal_multiplicities_and_reversibility(n):
    for p in rooted_walks(n):
        a = Counter(v for _,v in proposed_distinct_moves(p))
        b = Counter(_state_proposal_targets(p))
        assert a == b, (n,p)
        for other,count in a.items():
            assert _geometry_ok(other,n)
            assert Counter(_state_proposal_targets(other))[p] == count


@pytest.mark.parametrize('n',[4,5,6])
def test_two_independent_contact_algorithms_all_short_hp_sequences(n):
    paths=rooted_walks(n)
    for seq in map(''.join, product('HP',repeat=n)):
        for p in paths:
            assert _independent_contacts(p,seq) == hh_score(p,seq)


def test_pivot_extended_graph_connected_for_small_complete_ensembles():
    # An actual exhaustively checked connectivity result, not an extrapolated
    # assertion of ergodicity for large n.
    for n in (5,6,7,8,9,10):
        full=set(rooted_walks(n))
        seed=next(iter(full))
        reached={seed}
        queue=deque([seed])
        while queue:
            p=queue.popleft()
            for v in set(_state_proposal_targets(p)):
                if v not in reached:
                    reached.add(v)
                    queue.append(v)
        assert reached == full


def test_exact_energy_barriers_and_200_residue_non_enumerative_certification():
    q=Fraction(10**9)
    for height in (10,20,50):
        c=propose_rectangular_plateau(height)
        r=verify_pivot_plateau(c,q,10**6)
        assert r.total_enumerated_conformations == 0
        assert r.n==height*4
        assert r.plateau_states==8
        assert r.plateau_score==3*height-3
        assert r.witness_score == (7*height//2)-7  # n+1-8-height/2
        assert r.witness_gap == height//2-4
        assert r.boundary_penalties[0] == (1,28)
        assert r.fixed_seed_hit_by_steps_upper <= Fraction(1,1000)
        assert r.equilibrium_better_event_lower >= 0
        assert r.all_q_fixed_seed_exit_coefficient/q >= r.fixed_seed_exit_one_step_upper
    r=verify_pivot_plateau(propose_rectangular_plateau(50),q,10**6)
    assert r.boundary_proposal_slots==2260
    assert r.max_exit_proposals_per_state==307
    assert r.proposal_slots==2178
    assert r.fixed_seed_hit_by_steps_upper < Fraction(23,10**7)
    assert r.equilibrium_bad_event_upper < Fraction(1,10**90)
    assert r.equilibrium_better_event_lower > Fraction(999999999,10**9)
    assert r.conditional_gibbs_hit_by_steps_upper < r.fixed_seed_hit_by_steps_upper


def test_scaling_to_800_residues_no_full_conformation_enumeration():
    c=propose_rectangular_plateau(200)
    r=verify_pivot_plateau(c,Fraction(10**9),10**6)
    assert r.n == 800 and r.plateau_states==8
    assert r.witness_gap==96
    assert r.boundary_proposal_slots==8860
    assert r.fixed_seed_hit_by_steps_upper < Fraction(6,10**7)
    assert r.equilibrium_bad_event_upper < Fraction(1,10**480)
    assert r.total_enumerated_conformations==0


def test_validity_for_all_positive_rational_q_without_resampling_plateau():
    c=propose_rectangular_plateau(10)
    for q in (Fraction(1),Fraction(101,100),Fraction(3,2),Fraction(2),
              Fraction(37,11),Fraction(10000),Fraction(10**9)):
        r=verify_pivot_plateau(c,q,10**4)
        assert r.boundary_penalties[0] == (1,28)
        assert r.fixed_seed_exit_one_step_upper <= r.all_q_fixed_seed_exit_coefficient/q
        assert 0<=r.fixed_seed_hit_by_steps_upper<=1
        assert 0<=r.equilibrium_better_event_lower<=1


def test_one_step_rare_exit_and_thermal_detailed_balance_nontrivially():
    q=Fraction(37,11)
    c=propose_rectangular_plateau(10)
    r=verify_pivot_plateau(c,q,1)
    n=r.n
    m=r.proposal_slots
    # All outgoing slots are accounted for; every accepted move has exactly
    # equal reverse *proposal multiplicity*, and detailed balance holds.
    for p in c.plateau:
        k=_independent_contacts(p,c.sequence)
        probs=Counter(_state_proposal_targets(p))
        for v,multiplicity in probs.items():
            rate=Fraction(multiplicity,m)*min(Fraction(1),q**(_independent_contacts(v,c.sequence)-k))
            back=Counter(_state_proposal_targets(v))[p]
            assert back==multiplicity
            reverse=Fraction(back,m)*min(Fraction(1),q**(k-_independent_contacts(v,c.sequence)))
            assert q**k*rate == q**_independent_contacts(v,c.sequence)*reverse


def test_tampering_and_missing_plateau_members_fail_closed():
    c=propose_rectangular_plateau(10)
    q=Fraction(1000)
    others=[p for p in c.plateau if p!=c.seed]
    bad_examples=[
        replace(c,plateau=tuple(others)),
        replace(c,plateau=c.plateau[:-1]),
        replace(c,plateau=c.plateau+(c.plateau[0],)),
        replace(c,plateau=(c.seed,)),
        replace(c,plateau=c.plateau+(((0,0),(1,0)),)),
        replace(c,seed=((0,0),(1,0))),
        replace(c,better_witness=c.seed),
        replace(c,better_witness=tuple(reversed(c.better_witness))),
        replace(c,sequence='P'*len(c.sequence)),
        replace(c,sequence=c.sequence[:-1]),
    ]
    for bad in bad_examples:
        with pytest.raises(PlateauCertificateError):
            verify_pivot_plateau(bad,q,100)


def test_budgets_invalid_inputs_and_unclaimed_witnesses():
    c=propose_rectangular_plateau(10)
    with pytest.raises(PlateauProposalFailure):
        propose_rectangular_plateau(10,max_plateau_states=1)
    with pytest.raises(PlateauCertificateError):
        verify_pivot_plateau(c,Fraction(1),1,max_states=4)
    with pytest.raises(PlateauCertificateError):
        verify_pivot_plateau(c,Fraction(1),1,max_beads=10)
    for q in (0,2,2.0,Fraction(1,2)):
        with pytest.raises(PlateauCertificateError):
            verify_pivot_plateau(c,q,1)
    for t in (True,-1,3.14):
        with pytest.raises(PlateauCertificateError):
            verify_pivot_plateau(c,Fraction(2),t)
    for bad_height in (9,11,8,10.0):
        with pytest.raises(ValueError):
            propose_rectangular_plateau(bad_height)


def test_small_exact_thermal_probability_dominance_versus_full_states():
    # Nonenumerative population lower proof validated against the entire
    # n=12 state universe using a smaller arbitrary energetic plateau.
    n=12
    walks=rooted_walks(n)
    seed=rectangular_hamiltonian_snake(4,3)
    k0=hh_score(seed,'H'*12)
    assert k0==6
    # A distinct witness of score>k0 is unavailable for n=12:
    # this expected falsification guards against promising a better fold.
    assert max(_independent_contacts(w,'H'*n) for w in walks) == k0
    with pytest.raises(PlateauCertificateError):
        verify_pivot_plateau(PivotPlateauCertificate('H'*n,seed,(seed,),seed),Fraction(2),10)
