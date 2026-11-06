"""Proof-carrying exact variance maximum for a relaxed rational score ensemble.

Unknown blocks have a total nonnegative rational mass at most cap, distributed
over integral scores in [low, high]. An independently checked dual inequality
and an attaining primal mass allocation prove sharpness within that relaxation.
No Lean refinement, new mathematical inequality, or physical-protein claim.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction

KnownAtom = tuple[int, Fraction]
ScoreCap = tuple[int, int, Fraction]
Allocation = tuple[tuple[int, Fraction], ...]


@dataclass(frozen=True)
class VarianceDualCertificate:
    upper: Fraction
    center: Fraction
    allocations: tuple[Allocation, ...]
    iterations: int


def _validate(known: tuple[KnownAtom, ...], caps: tuple[ScoreCap, ...]) -> None:
    if not isinstance(known, tuple) or not isinstance(caps, tuple):
        raise ValueError("known and caps must be tuples")
    for atom in known:
        if (not isinstance(atom, tuple) or len(atom) != 2 or
                type(atom[0]) is not int or not isinstance(atom[1], Fraction) or
                atom[1] <= 0):
            raise ValueError("known atoms need integral scores and positive rational mass")
    for cap in caps:
        if (not isinstance(cap, tuple) or len(cap) != 3 or
                type(cap[0]) is not int or type(cap[1]) is not int or
                cap[0] < 0 or cap[1] < cap[0] or
                not isinstance(cap[2], Fraction) or cap[2] <= 0):
            raise ValueError("caps require 0 <= low <= high and positive rational cap")
    if not known and not caps:
        raise ValueError("empty possible ensemble")


def _moments(known: tuple[KnownAtom, ...], caps: tuple[ScoreCap, ...],
             allocations: tuple[Allocation, ...]) -> tuple[Fraction, Fraction, Fraction]:
    z = sum((w for _, w in known), Fraction(0))
    m1 = sum((k*w for k, w in known), Fraction(0))
    m2 = sum((k*k*w for k, w in known), Fraction(0))
    for terms in allocations:
        for k, w in terms:
            z += w
            m1 += k*w
            m2 += k*k*w
    return z, m1, m2


def verify_variance_dual(known: tuple[KnownAtom, ...], caps: tuple[ScoreCap, ...],
                         certificate: VarianceDualCertificate) -> bool:
    """Independently check dual universal upper and primal exact attainment."""
    _validate(known, caps)
    if not isinstance(certificate, VarianceDualCertificate):
        raise ValueError("incorrect certificate type")
    c, upper = certificate.center, certificate.upper
    if not isinstance(c, Fraction) or not isinstance(upper, Fraction) or upper < 0:
        raise ValueError("center and upper require exact rational values")
    if not isinstance(certificate.allocations, tuple) or len(certificate.allocations) != len(caps):
        raise ValueError("one allocation tuple per cap required")
    if type(certificate.iterations) is not int or certificate.iterations < 0:
        raise ValueError("nonnegative integer iteration count required")
    residual = sum((w*((k-c)**2-upper) for k, w in known), Fraction(0))
    for low, high, mass in caps:
        residual += mass*max(Fraction(0), max((low-c)**2, (high-c)**2)-upper)
    if residual > 0:
        raise ValueError("dual dominance fails")
    for (low, high, mass), terms in zip(caps, certificate.allocations):
        if not isinstance(terms, tuple):
            raise ValueError("allocation must be a tuple")
        total = Fraction(0)
        for term in terms:
            if (not isinstance(term, tuple) or len(term) != 2 or
                    type(term[0]) is not int or not low <= term[0] <= high or
                    not isinstance(term[1], Fraction) or term[1] < 0):
                raise ValueError("invalid primal allocation")
            total += term[1]
        if total > mass:
            raise ValueError("allocation exceeds block cap")
    z, m1, m2 = _moments(known, caps, certificate.allocations)
    if z <= 0:
        raise ValueError("empty primal distribution")
    if m1 != c*z or m2*z-m1*m1 != upper*z*z:
        raise ValueError("primal fails to attain mean and variance")
    return True


def _oracle(known: tuple[KnownAtom, ...], caps: tuple[ScoreCap, ...],
            center: Fraction) -> tuple[Fraction, tuple[int | None, ...], Fraction, Fraction]:
    """Exact fractional maximizer of centered second moment at fixed center."""
    z0 = sum((w for _, w in known), Fraction(0))
    d0 = sum((w*(k-center)**2 for k, w in known), Fraction(0))
    candidates = []
    for low, high, cap in caps:
        dlo, dhi = (low-center)**2, (high-center)**2
        score = low if dlo >= dhi else high
        candidates.append((max(dlo, dhi), cap, score))
    if not z0:
        i = max(range(len(candidates)), key=lambda j:candidates[j][0])
        choices = tuple(x[2] if j == i else None for j, x in enumerate(candidates))
        k = candidates[i][2]
        return candidates[i][0], choices, Fraction(k), Fraction(k*k)
    mass, numerator = z0, d0
    best = numerator/mass
    for dev, cap, _ in sorted(candidates, reverse=True):
        if dev > best:
            mass += cap
            numerator += cap*dev
            best = numerator/mass
    choices = tuple(score if dev > best else None for dev, _, score in candidates)
    allocation = tuple(((score, caps[i][2]),) if score is not None else ()
                       for i, score in enumerate(choices))
    z, m1, m2 = _moments(known, caps, allocation)
    assert z > 0 and (m2-2*center*m1+center*center*z)/z == best
    return best, choices, m1/z, m2/z


def _restricted_minimum(lines: tuple[tuple[Fraction, Fraction], ...],
                        left: Fraction, right: Fraction) -> tuple[Fraction, Fraction]:
    """Minimize quadratic envelope exactly at vertices and intersections."""
    places = {left, right}
    for mu, _ in lines:
        places.add(min(right, max(left, mu)))
    for i, (mu1, m21) in enumerate(lines):
        for mu2, m22 in lines[i+1:]:
            if mu1 != mu2:
                cut = (m21-m22)/(2*(mu1-mu2))
                if left <= cut <= right:
                    places.add(cut)
    c = min(places, key=lambda x:(x*x+max(m2-2*x*mu for mu, m2 in lines), x))
    value = c*c + max(m2-2*c*mu for mu, m2 in lines)
    return c, value


def maximize_relaxed_variance(known: tuple[KnownAtom, ...], caps: tuple[ScoreCap, ...],
                              *, max_iterations: int = 2048) -> VarianceDualCertificate:
    """Exact rational exchange algorithm; verified independently, fails closed."""
    _validate(known, caps)
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValueError("max_iterations must be positive integer")
    support = [k for k, _ in known]
    for low, high, _ in caps:
        support.extend((low, high))
    left, right = Fraction(min(support)), Fraction(max(support))
    z0 = sum((w for _, w in known), Fraction(0))
    if not z0:
        low_i = next(i for i, (lo, _, _) in enumerate(caps) if lo == left)
        high_i = next(i for i, (_, hi, _) in enumerate(caps) if hi == right)
        c = (left+right)/2
        v = (right-left)**2/4
        groups: list[list[tuple[int, Fraction]]] = [[] for _ in caps]
        mass = (min(caps[low_i][2], caps[high_i][2])/2 if low_i == high_i
                else min(caps[low_i][2], caps[high_i][2]))
        groups[low_i].append((int(left), mass))
        if right != left:
            groups[high_i].append((int(right), mass))
        elif low_i != high_i:
            groups[high_i].append((int(right), mass))
        proof = VarianceDualCertificate(v, c, tuple(tuple(g) for g in groups), 0)
        verify_variance_dual(known, caps, proof)
        return proof
    if not caps:
        z, m1, m2 = _moments(known, caps, ())
        c = m1/z
        proof = VarianceDualCertificate(m2/z-c*c, c, (), 0)
        verify_variance_dual(known, caps, proof)
        return proof
    models: dict[tuple[int | None, ...], tuple[Fraction, Fraction]] = {}
    _, first, mu, m2 = _oracle(known, caps, (left+right)/2)
    models[first] = (mu, m2)
    for iteration in range(1, max_iterations+1):
        c, lower = _restricted_minimum(tuple(models.values()), left, right)
        upper, choice, mu, m2 = _oracle(known, caps, c)
        if upper == lower:
            active = [(p, a) for p, a in models.items()
                      if c*c+a[1]-2*c*a[0] == lower]
            active.append((choice, (mu, m2)))
            under = [(p, a) for p, a in active if a[0] <= c]
            over = [(p, a) for p, a in active if a[0] >= c]
            if not under or not over:
                raise AssertionError("no bracketing active maximizers")
            pleft, (mu_l, _) = min(under, key=lambda x:x[1][0])
            pright, (mu_r, _) = max(over, key=lambda x:x[1][0])
            alpha = Fraction(1) if mu_l == mu_r else (mu_r-c)/(mu_r-mu_l)
            zl = z0+sum((cap for k, (_, _, cap) in zip(pleft, caps)
                         if k is not None), Fraction(0))
            zr = z0+sum((cap for k, (_, _, cap) in zip(pright, caps)
                         if k is not None), Fraction(0))
            norm = alpha/zl+(1-alpha)/zr
            fl, fr = (alpha/zl)/norm, ((1-alpha)/zr)/norm
            groups: list[list[tuple[int, Fraction]]] = [[] for _ in caps]
            for choices, factor in ((pleft, fl), (pright, fr)):
                for i, score in enumerate(choices):
                    if score is not None and factor:
                        groups[i].append((score, caps[i][2]*factor))
            proof = VarianceDualCertificate(upper, c, tuple(tuple(g) for g in groups), iteration)
            verify_variance_dual(known, caps, proof)
            return proof
        if choice in models:
            raise AssertionError("separation oracle produced a repeated invalid cut")
        models[choice] = (mu, m2)
    raise RuntimeError("no sharp variance certificate within iteration budget")
