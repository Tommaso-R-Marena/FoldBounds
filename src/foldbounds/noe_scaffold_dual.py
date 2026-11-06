"""Exact rational directional NOE / fixed-scaffold incompatibility certificates.

Mathematical result is CONDITIONAL on frozen processed NEF restraint #470,
seven fixed Val69 heavy atoms, and per-H/anchor distance tolerances. The
producer is untrusted; this module independently replays every dual vector
identity and every primal squared-distance comparison with integers/Fraction.
Neither the uncertainty tolerance nor the hydrogen positions are inferred
from independent experiments. No Lean-checked or physical-chemistry claim.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


def _v(x):
    return type(x) in (list, tuple) and len(x) == 3 and all(type(t) is int for t in x)


def _sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def _dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def _ceil(n):
    r = isqrt(n)
    return r if r*r == n else r+1


@dataclass(frozen=True)
class VerifiedNoGo:
    delta_mA: int
    required_mA2: int
    support_mA2: Fraction
    support_h_mA2: Fraction
    support_hb_mA2: Fraction

    @property
    def strict_margin_mA2(self):
        return self.required_mA2-self.support_mA2


@dataclass(frozen=True)
class VerifiedPrimal:
    cap_ticks: int
    worst_distance_change_ticks: int
    checked_local_pairs: int


def parse(raw):
    """Accept only the expected source-pinned structured input; not a proof."""
    if (type(raw) is not dict or raw.get('schema') !=
            'FoldBounds.Val69NOELocalDual.v1'):
        raise ValueError('unsupported certificate')
    return raw


def _validate(raw):
    try:
        if (type(raw) is not dict or raw.get('schema') !=
                'FoldBounds.Val69NOELocalDual.v1' or
                raw.get('sample') != '6nbn' or
                type(raw.get('processed_restraint_id')) is not int or
                raw['processed_restraint_id'] != 470 or
                type(raw.get('processed_restraint_upper_mA')) is not int or
                raw['processed_restraint_upper_mA'] <= 0 or
                type(raw.get('allowed_change_in_each_heavy_to_hydrogen_distance_mA')) is not int or
                raw['allowed_change_in_each_heavy_to_hydrogen_distance_mA'] < 0):
            return False
        a=raw['fixed_heavy_atoms_val69']
        h=raw['target_hydrogens_val69']
        if (type(a) is not dict or len(a) < 4 or len(a)>100 or
                any(type(k) is not str or not _v(v) for k,v in a.items()) or
                len({tuple(v) for v in a.values()}) != len(a) or
                type(h) is not dict or not _v(h['H']) or not _v(h['HB']) or
                type(raw['proof']) is not dict):
            return False
        return True
    except (KeyError,TypeError,ValueError):
        return False


def _upper_support(endpoint, objective, anchors, delta, rows):
    if type(rows) is not list or not 1 <= len(rows) <= 50:
        return None
    direction=[Fraction(0)]*3
    amount=Fraction(0)
    for row in rows:
        if (type(row) is not list or len(row) != 5 or
                type(row[0]) is not int or row[0] not in (-1,1) or
                type(row[1]) is not str or type(row[2]) is not str or
                row[1] == row[2] or row[1] not in anchors or row[2] not in anchors or
                type(row[3]) is not int or row[3] < 0 or
                type(row[4]) is not int or row[4] <= 0):
            return None
        sign, j, k, numerator, denominator = row
        lam=Fraction(numerator,denominator)
        ai,aj=anchors[j],anchors[k]
        diff=_sub(ai,aj)
        for t in range(3):
            direction[t]+=lam*sign*diff[t]
        def err(site):
            displacement=_sub(endpoint,site)
            root=_ceil(_dot(displacement,displacement))
            return 2*root*delta+delta*delta
        amount += lam*Fraction(err(ai)+err(aj),2)
    if tuple(direction) != tuple(objective):
        return None
    return amount


def verify_no_go(raw):
    """Fail-closed replay of a nonnegative rational LP dual.

    If each old-to-new heavy-anchor distance changes by at most delta,
    the squared-distance difference for an anchor has magnitude at most
    E=2 ceil(original_distance) delta + delta**2. Subtracting pairs
    cancels the unknown displacement squared term. Nonnegative signed
    slab multipliers prove a directional displacement bound. An upper
    NOE contact requires a strictly larger directional displacement,
    proving impossibility. Inconclusive returns None, not feasibility.
    """
    if not _validate(raw):
        return None
    try:
        a=raw['fixed_heavy_atoms_val69'];h=raw['target_hydrogens_val69']
        x,y=h['H'],h['HB'];v=_sub(y,x)
        sq=_dot(v,v)
        upper=raw['processed_restraint_upper_mA']
        need=sq-upper*_ceil(sq)
        if need <= 0:
            return None
        delta=raw['allowed_change_in_each_heavy_to_hydrogen_distance_mA']
        cap_h=_upper_support(x,v,a,delta,raw['proof']['H'])
        cap_hb=_upper_support(y,tuple(-z for z in v),a,delta,raw['proof']['HB'])
        if cap_h is None or cap_hb is None or cap_h+cap_hb>=need:
            return None
        return VerifiedNoGo(delta,need,cap_h+cap_hb,cap_h,cap_hb)
    except (KeyError,TypeError,ValueError,ZeroDivisionError,OverflowError):
        return None


def verify_local_geometry_primal(raw):
    """Fail-closed exact integer upper-feasibility witness.

    The constructive positions use fine units of 1e-6 Å. Integer
    floor/ceil square roots enclose every irrational distance, making
    the distance-change comparison conservative.
    """
    if not _validate(raw):
        return None
    try:
        x,y=raw['target_hydrogens_val69']['H'],raw['target_hydrogens_val69']['HB']
        p=raw['constructive_feasible_at_larger_tolerance']
        S=p['scale_ticks_per_milliangstrom']
        cap=p['guaranteed_distance_deviation_cap_ticks']
        if type(S) is not int or S<1 or type(cap) is not int or cap<0:
            return None
        aa=p['moves']['H'];bb=p['moves']['HB']
        if not _v(aa) or not _v(bb):
            return None
        xnew=tuple(S*t+u for t,u in zip(x,aa))
        ynew=tuple(S*t+u for t,u in zip(y,bb))
        dx=_sub(xnew,ynew)
        if _dot(dx,dx) > (S*raw['processed_restraint_upper_mA'])**2:
            return None
        worst=0;count=0
        for old,new in ((x,xnew),(y,ynew)):
            for anchor in raw['fixed_heavy_atoms_val69'].values():
                old_sq=_dot(_sub(old,anchor),_sub(old,anchor))*S*S
                scaled=tuple(S*t for t in anchor)
                new_sq=_dot(_sub(new,scaled),_sub(new,scaled))
                hi0=_ceil(old_sq);lo0=isqrt(old_sq)
                hi1=_ceil(new_sq);lo1=isqrt(new_sq)
                worst=max(worst,hi1-lo0,hi0-lo1)
                count+=1
        if worst>cap:
            return None
        return VerifiedPrimal(cap,worst,count)
    except (KeyError,TypeError,ValueError,OverflowError):
        return None


def with_delta(raw, delta):
    """Copy with a distinct delta; prevents mutating the shared fixture."""
    if type(delta) is not int:
        raise ValueError('integer milliangstrom delta required')
    return {**raw,'allowed_change_in_each_heavy_to_hydrogen_distance_mA':delta}


def check_inherited_restraint(raw, records):
    """Bridge the exact Val69 data to the independent PR #80 source fixture."""
    if not _validate(raw):
        return False
    try:
        from .noe_restraint_audit import AtomRef,NOERestraint,verify_record
        matches=[r for r in records if type(r) is NOERestraint
                 and r.identity == raw['processed_restraint_id']]
        if len(matches)!=1 or not verify_record(matches[0]):
            return False
        r=matches[0]
        ref={AtomRef(69,'H'):tuple(raw['target_hydrogens_val69']['H']),
             AtomRef(69,'HB'):tuple(raw['target_hydrogens_val69']['HB'])}
        return (set((r.left,r.right))==set(ref) and
                r.left_position_mA==ref[r.left] and
                r.right_position_mA==ref[r.right] and
                r.upper_mA==raw['processed_restraint_upper_mA'] and
                r.lower_mA==raw['processed_restraint_lower_mA'])
    except (ImportError,ValueError,TypeError,KeyError,AttributeError):
        return False
