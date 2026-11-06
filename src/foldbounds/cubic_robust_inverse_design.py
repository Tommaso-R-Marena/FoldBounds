"""Proof-gated inverse design of a symmetry-unique cubic HP native orbit.

This is a *finite 3D lattice model*, not protein structure prediction.

The untrusted cut-tree producer proposes pruning. An independently compiled
C++ verifier checks every cut, counts all surviving high-contact leaves, and
emits all of their Cartesian paths. Python independently rescans each path,
reconstructs the H-H contact bitmask, and checks group/orbit identities.
A separate three-threshold verifier certifies the score gap and excludes
all unvisited lower-contact configurations from competing with the design.

The design adds pair-specific attractive contact coefficients plus ONE
Euclidean distance restraint, all invariant under cubic lattice isometries.
The exact robust Boltzmann bound holds for *all correlated*, per-contact
perturbations in the specified uniform interval, not just sampled weights.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path
from hashlib import sha256
import subprocess

from .cubic_cut_certificate import (
    CheckedNativeGap, certify_native_gap, contact_score, inspect_cut_header,
)

Point = tuple[int, int, int]
Walk = tuple[Point, ...]


def cubic_candidate_pairs(length: int) -> tuple[tuple[int, int], ...]:
    if type(length) is not int or not 2 <= length <= 20:
        raise ValueError('invalid chain length')
    return tuple((i,j) for j in range(length) for i in range(j-2)
                 if (j-i)%2 == 1)


def _hh_mask(sequence: str, coordinates: Walk) -> int:
    pairs = cubic_candidate_pairs(len(sequence))
    if len(pairs) > 63:
        raise ValueError('more than 63 candidate pair bits')
    if contact_score(sequence, coordinates) < 0:
        raise AssertionError('unreachable score')
    return sum(1<<bit for bit,(i,j) in enumerate(pairs)
               if sequence[i]==sequence[j]=='H'
               and sum(abs(a-b) for a,b in zip(coordinates[i],coordinates[j]))==1)


def _d2(path: Walk, i: int, j: int) -> int:
    return sum((path[i][k]-path[j][k])**2 for k in range(3))


def _transverse_orbit(path: Walk) -> frozenset[Walk]:
    """The eight lattice symmetries fixing an oriented initial +x bond."""
    return frozenset(tuple((x, sy*(z if swap else y), sz*(y if swap else z))
                     for x,y,z in path)
                     for swap,sy,sz in product((False,True),(-1,1),(-1,1)))


@dataclass(frozen=True)
class CheckedGroundMaps:
    sequence: str
    threshold: int
    cutoff_sha256: str
    high_count: int
    upper_bins: tuple[tuple[int,int], ...]
    # Each row contains an *independently rescored* complete anchored walk.
    entries: tuple[tuple[int,Walk], ...]

    @property
    def multiplicities(self) -> tuple[tuple[int,int], ...]:
        return tuple(sorted(Counter(mask for mask,_ in self.entries).items()))


def verify_ground_map_cut(
    cut_file: str | Path, independent_checker: str | Path,
    *, max_seconds: int = 240,
) -> CheckedGroundMaps:
    """Independently verify a COMPLETE cut tree, then rescore each emitted leaf.

    The checker takes exactly three arguments: cut, output histogram, output
    all-coordinate paths. It must be compiled from cubic_ground_map_checker.cpp.
    It is not sufficient to provide a claimed histogram with no cut tree.
    """
    import tempfile
    cut=Path(cut_file)
    checker=Path(independent_checker)
    header=inspect_cut_header(cut)  # untrusted, for expected output only
    if not checker.is_file():
        raise FileNotFoundError('independent C++ ground-map verifier required')
    if type(max_seconds) is not int or max_seconds < 1:
        raise ValueError('positive time budget required')
    with tempfile.TemporaryDirectory(prefix='fb_verified_ground_') as directory:
        hist=Path(directory)/'maps.tsv'
        coords=Path(directory)/'all_coords.tsv'
        try:
            run=subprocess.run(
                [str(checker.resolve()),str(cut.resolve()),str(hist),str(coords)],
                check=False,capture_output=True,text=True,timeout=max_seconds,
            )
        except subprocess.TimeoutExpired as e:
            raise ValueError('independent proof verifier timed out') from e
        if run.returncode:
            raise ValueError('independent proof verifier rejected cut: '
                             +run.stderr[:350])
        tag=(f'VERIFIED seq={header.sequence} threshold={header.threshold} '
             f'nodes={header.total_nodes} cuts={header.cuts} '
             f'leaves={header.terminal_leaves} '
             f'high={header.qualifying_conformations} mode={header.mode} bins=')
        if not run.stdout.startswith(tag) or run.stdout.count('\n')!=1:
            raise ValueError('independent proof summary/header mismatch')
        bins=[]
        for word in run.stdout[len(tag):].strip().split(','):
            if not word:continue
            pieces=word.split(':')
            if len(pieces)!=2:raise ValueError('malformed upper-score bucket')
            score,count=map(int,pieces)
            if score<0 or score>=header.threshold or count<=0:
                raise ValueError('invalid upper-score bucket')
            bins.append((score,count))
        if bins!=sorted(set(bins)):
            raise ValueError('unsorted or duplicate upper-score buckets')
        map_counts=Counter()
        paths=[]
        n=len(header.sequence)
        for row in coords.read_text().splitlines():
            pieces=list(map(int,row.split()))
            if len(pieces)!=1+3*n:
                raise ValueError('invalid complete-walk coordinates')
            advertised_mask=pieces[0]
            walk=tuple(tuple(pieces[1+3*i:1+3*i+3]) for i in range(n))
            score=contact_score(header.sequence,walk)
            mask=_hh_mask(header.sequence,walk)
            if score<header.threshold or mask!=advertised_mask:
                raise ValueError('coordinate leaf not independently justified')
            paths.append((mask,walk))
            map_counts[mask]+=1
        expected=[]
        for row in hist.read_text().splitlines():
            words=list(map(int,row.split()))
            if len(words)!=2+3*n or words[1]<1:
                raise ValueError('invalid proof map summary')
            mask,count=words[:2]
            witness=tuple(tuple(words[2+3*i:5+3*i]) for i in range(n))
            if _hh_mask(header.sequence,witness)!=mask:
                raise ValueError('map summary witness is inconsistent')
            expected.append((mask,count))
        if tuple(expected)!=tuple(sorted(map_counts.items())):
            raise ValueError('claimed map count differs from independently checked paths')
        if len(paths)!=header.qualifying_conformations or len(set(walk for _,walk in paths))!=len(paths):
            raise ValueError('missing/duplicate high-contact conformations')
        return CheckedGroundMaps(
            header.sequence,header.threshold,sha256(cut.read_bytes()).hexdigest(),
            len(paths),tuple(bins),tuple(paths)
        )


@dataclass(frozen=True)
class RobustUniqueOrbit:
    sequence: str
    ground_contact_count: int
    ground_conformations: int
    gap: int
    designated_contact_mask: int
    target_distance_pair: tuple[int,int]
    target_squared_distance: int
    target_path: Walk
    target_orbit_size: int
    design_energy: int                 # quarter-energy units
    # Every non-designated K=M conformer, grouped by worst-case robust
    # energetic margin from the chosen orbit (radius 0 or 1).
    remaining_high: tuple[tuple[int,int,int], ...]  # missing favored edges,bonus,count
    max_total_conformations: int
    verified_excited_upper_bins: tuple[tuple[int,int], ...]

    def robust_energy_margin(self, missing: int, bonus: int, radius: int) -> int:
        # In the K=M sector: missing favored edges = extra unfavored edges.
        return (4-2*radius)*missing+2*(1-bonus)

    def lower_orbit_probability(self, q: Fraction, *, radius: int=1) -> Fraction:
        """Exact occupancy floor for every correlated contact perturbation.

        Nominal HH energy in quarter-units: 8 for edges in target, 4 for
        others. A D4-invariant squared-distance restraint adds 2. Every
        HH edge coefficient has *fixed but arbitrary* perturbation η_e with
        |η_e|<=radius, for radius 0 or 1. Thus all competitors are bounded
        in ratios relative to the same target orbit with shared errors
        canceling. Full geometry beyond the ground-state sector is NOT
        enumerated. q is the dimensionless contact weight per quarter unit.
        """
        if type(q) is not Fraction or q<1 or type(radius) is not int or radius not in (0,1):
            raise ValueError('q must be exact Fraction >=1; supported radius=0 or 1')
        g=self.target_orbit_size
        high_ratio=sum((Fraction(count,g*q**self.robust_energy_margin(a,b,radius))
                        for a,b,count in self.remaining_high),Fraction(0))
        low_gap=(8-radius)*self.gap
        coarse=Fraction(self.max_total_conformations-self.ground_conformations,
                        g*q**low_gap)
        cuts=sum((Fraction(count,g*q**((8-radius)*(self.ground_contact_count-min(k,self.ground_contact_count-self.gap))))
                  for k,count in self.verified_excited_upper_bins),Fraction(0))
        return Fraction(1)/(1+high_ratio+min(coarse,cuts))

    def minimum_integer_q(self, target: Fraction, *, radius: int=1) -> int:
        if type(target) is not Fraction or not 0<target<1:
            raise ValueError('target must be rational strictly between zero and one')
        lo,hi=1,2
        if self.lower_orbit_probability(Fraction(1),radius=radius)>=target:
            return 1
        while self.lower_orbit_probability(Fraction(hi),radius=radius)<target:
            hi*=2
            if hi>1<<50:
                raise ValueError('occupancy bound threshold exceeds safety limit')
        while hi-lo>1:
            mid=(lo+hi)//2
            if self.lower_orbit_probability(Fraction(mid),radius=radius)>=target:
                hi=mid
            else:
                lo=mid
        return hi


def design_robust_unique_orbit(
    full_gap: CheckedNativeGap,
    checked_maps: CheckedGroundMaps,
    *, selected_map: int | None=None,
    distance_pair: tuple[int,int] | None=None,
    squared_distance: int | None=None,
) -> RobustUniqueOrbit:
    """Select one D4-invariant ground orbit and prove energy separation.

    The cut-tree gap certifies no states in (M-gap, M) or >M. This function
    refuses if the checked full ground-state set cannot isolate exactly one
    full D4 orbit with one two-residue squared-distance indicator.
    """
    s=full_gap.sequence
    if (s!=checked_maps.sequence or full_gap.ground_contacts!=checked_maps.threshold
            or full_gap.ground_conformations!=checked_maps.high_count or full_gap.gap<1):
        raise ValueError('ground-map paths not backed by the validated full excitation gap')
    n=len(s)
    by: dict[int,list[Walk]]=defaultdict(list)
    for mask,path in checked_maps.entries:
        if _hh_mask(s,path)!=mask or mask.bit_count()!=full_gap.ground_contacts:
            raise ValueError('untrusted ground-state map/score')
        by[mask].append(path)
    if sum(map(len,by.values()))!=full_gap.ground_conformations:
        raise ValueError('incomplete full high-contact ground-state census')
    if selected_map is not None and (type(selected_map) is not int or selected_map not in by):
        raise ValueError('requested contact map absent from exact ground sector')
    if distance_pair is not None and (not isinstance(distance_pair,tuple)
            or len(distance_pair)!=2 or any(type(x) is not int for x in distance_pair)
            or not 0<=distance_pair[0]<distance_pair[1]<n):
        raise ValueError('invalid distance restraint pair')
    if squared_distance is not None and (type(squared_distance) is not int or squared_distance<0):
        raise ValueError('invalid squared distance')

    def choose(mask,paths):
        group=set(paths)
        if len(group)!=len(paths):
            raise ValueError('duplicate geometry in ground-state proof')
        remaining=set(paths)
        orbits=[]
        while remaining:
            representative=min(remaining)
            orbit=_transverse_orbit(representative)
            if orbit!=orbit&group or not orbit<=remaining:
                raise ArithmeticError('D4 symmetry does not close ground-map sector')
            orbits.append((representative,orbit))
            remaining-=orbit
        if distance_pair is not None:
            pairs=(distance_pair,)
        else:
            pairs=((i,j) for i in range(n) for j in range(i+2,n))
        for i,j in pairs:
            for representative,orbit in orbits:
                d2=_d2(representative,i,j)
                if squared_distance is not None and squared_distance!=d2:
                    continue
                selected={p for p in group if _d2(p,i,j)==d2}
                if selected==orbit:
                    return representative,orbit,(i,j),d2
        return None

    all_maps=sorted(by,key=lambda key:(len(by[key]),key))
    if selected_map is not None:
        all_maps=[selected_map]
    selection=None
    for mask in all_maps:
        result=choose(mask,by[mask])
        if result is not None:
            selection=(mask,*result)
            break
    if selection is None:
        raise ValueError('no D4-unique orbit can be isolated by one squared-distance restraint')
    mask,rep,orbit,pair,d2=selection
    if len(orbit) not in (4,8):
        raise ValueError('selected D4 orbit is trivial or degenerate')
    class_count=Counter()
    for m,paths in by.items():
        missing=(mask&~m).bit_count()
        extras=(m&~mask).bit_count()
        if missing!=extras:
            raise ArithmeticError('same-score contact-map parity mismatch')
        for p in paths:
            bonus=int(_d2(p,*pair)==d2)
            if m==mask and bonus:
                if p not in orbit:
                    raise ArithmeticError('restraint fails to identify precisely one D4 orbit')
                continue
            margin=(4-2)*missing+2*(1-bonus)  # worst-case at error radius 1
            if margin<=0:
                raise ValueError('no positive robust unique-native separation')
            class_count[missing,bonus]+=1
    if sum(class_count.values())+len(orbit)!=full_gap.ground_conformations:
        raise ArithmeticError('missed a ground-state competitor')
    # Any K<=M-gap state has relative energy gap >=(8-r)*(M-K),
    # valid for all pairwise correlated perturbations |eta_e|<=r<=1.
    if full_gap.gap<1:
        raise ValueError('no non-ground excitation gap')
    for k,amount in checked_maps.upper_bins:
        if type(k) is not int or type(amount) is not int or k<0 or amount<=0:
            raise ValueError('malformed verifier upper bins')
    return RobustUniqueOrbit(
        s,full_gap.ground_contacts,full_gap.ground_conformations,full_gap.gap,
        mask,pair,d2,rep,len(orbit),8*full_gap.ground_contacts+2,
        tuple(sorted((a,b,count) for (a,b),count in class_count.items())),
        full_gap.max_total_conformations,checked_maps.upper_bins,
    )


def validate_complete_robust_design(
    cut_files: tuple[str|Path,str|Path,str|Path],
    gap_checker: str|Path, ground_map_checker: str|Path,
    excited_witness: Walk, *, max_seconds: int=300,
    **design_options,
) -> RobustUniqueOrbit:
    """No trust shortcut: independently replay ALL three original cut trees.

    The ground-map verifier independently replays the middle tree a second
    time and independently rescans every emitted high-contact coordinate.
    """
    gap=certify_native_gap(
        cut_files,gap_checker,excited_witness,max_seconds=max_seconds
    )
    all_ground=verify_ground_map_cut(
        cut_files[1],ground_map_checker,max_seconds=max_seconds
    )
    return design_robust_unique_orbit(gap,all_ground,**design_options)
