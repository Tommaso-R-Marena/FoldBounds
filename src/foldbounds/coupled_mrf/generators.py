"""Model generators: Ising grids, random sparse graphs, rotamer-like Potts models.

All generated parameters are exact rationals.  Where a physical parameter is
irrational (e.g. q = exp(2 beta J)), a rational approximation is chosen once
and then *defines* the model; nothing downstream depends on the approximation.

Ising convention: spins in {0, 1}; an edge factor has table [q, 1, 1, q]
(weight q when the two spins agree), i.e. exp(beta J s s') up to the constant
exp(-beta J) per edge, with q = exp(2 beta J).  A field factor has table
[1, w] with w = exp(2 beta h) up to a constant.  Constants cancel in Z_B/Z_A.
"""
from __future__ import annotations

import math
import random
from fractions import Fraction
from typing import Sequence

from .model import CoupledModel, Factor

BETA_C = 0.5 * math.log(1.0 + math.sqrt(2.0))   # 2D square-lattice Ising, J = 1


def rational_approx(x: float, max_den: int = 1000) -> Fraction:
    return Fraction(x).limit_denominator(max_den)


def ising_q(beta_over_beta_c: float, max_den: int = 1000) -> Fraction:
    """Rational q ~ exp(2 beta J) with beta = ratio * beta_c (J = 1).

    exp(2 beta_c) = 1 + sqrt(2) exactly, so q = (1 + sqrt 2)^ratio.
    """
    return rational_approx((1.0 + math.sqrt(2.0)) ** beta_over_beta_c, max_den)


def site(r: int, c: int) -> str:
    return f"x{r}_{c}"


def grid_edges(L: int, M: int | None = None) -> list[tuple[tuple[int, int], tuple[int, int]]]:
    M = L if M is None else M
    edges = []
    for r in range(L):
        for c in range(M):
            if c + 1 < M:
                edges.append(((r, c), (r, c + 1)))
            if r + 1 < L:
                edges.append(((r, c), (r + 1, c)))
    return edges


def edge_factor(u: str, v: str, q: Fraction, name: str = "") -> Factor:
    q = Fraction(q)
    return Factor((u, v), (q, Fraction(1), Fraction(1), q), name)


def field_factor(u: str, w: Fraction, name: str = "") -> Factor:
    return Factor((u,), (Fraction(1), Fraction(w)), name)


def ising_grid(L: int, q: Fraction, w: Fraction = Fraction(1), perturbation: tuple | None = None,
               M: int | None = None, couplings: dict | None = None) -> CoupledModel:
    """L x M open-boundary Ising grid with a local perturbation.

    perturbation:
      ("field", (r, c), u)      -- B has an extra field factor [1, u] at (r, c);
      ("edge", ((r,c),(r2,c2)), qB) -- the edge has coupling q in A and qB in B;
      None                      -- A = B.
    couplings: optional {edge: q_e} overriding q on individual edges (shared).
    """
    M = L if M is None else M
    names = [site(r, c) for r in range(L) for c in range(M)]
    domains = {v: 2 for v in names}
    shared, only_a, only_b = [], [], []
    pert_edge = None
    if perturbation is not None and perturbation[0] == "edge":
        a, b = perturbation[1]
        pert_edge = (tuple(a), tuple(b)) if tuple(a) < tuple(b) else (tuple(b), tuple(a))
    for (a, b) in grid_edges(L, M):
        qe = q if couplings is None else couplings.get((a, b), q)
        f = edge_factor(site(*a), site(*b), qe, f"J{a}{b}")
        if pert_edge == (a, b):
            only_a.append(f)
            only_b.append(edge_factor(site(*a), site(*b), perturbation[2], f"JB{a}{b}"))
        else:
            shared.append(f)
    if Fraction(w) != 1:
        for v in names:
            shared.append(field_factor(v, w, f"h{v}"))
    if perturbation is not None and perturbation[0] == "field":
        r, c = perturbation[1]
        only_b.append(field_factor(site(r, c), perturbation[2], "dh"))
    elif perturbation is not None and perturbation[0] != "edge":
        raise ValueError(f"unknown perturbation {perturbation[0]!r}")
    meta = {"kind": "ising_grid", "L": L, "M": M, "q": str(q), "w": str(w),
            "perturbation": repr(perturbation)}
    return CoupledModel(tuple(names), domains, dict(domains), tuple(shared),
                        tuple(only_a), tuple(only_b), meta)


def center(L: int) -> tuple[int, int]:
    return (L // 2, L // 2)


def center_edge(L: int) -> tuple[tuple[int, int], tuple[int, int]]:
    r, c = center(L)
    return ((r, c - 1), (r, c)) if c >= 1 else ((r, c), (r, c + 1))


# ------------------------------------------------------------ random graphs
def random_regular_edges(n: int, degree: int, rng: random.Random, tries: int = 200
                         ) -> list[tuple[int, int]]:
    """Simple random d-regular graph by the configuration model with rejection."""
    if (n * degree) % 2:
        raise ValueError("n * degree must be even")
    for _ in range(tries):
        stubs = [i for i in range(n) for _ in range(degree)]
        rng.shuffle(stubs)
        edges = set()
        ok = True
        for i in range(0, len(stubs), 2):
            a, b = stubs[i], stubs[i + 1]
            if a == b or (min(a, b), max(a, b)) in edges:
                ok = False
                break
            edges.add((min(a, b), max(a, b)))
        if ok:
            return sorted(edges)
    raise RuntimeError("failed to sample a simple regular graph")


def random_rational(rng: random.Random, lo: Fraction, hi: Fraction, den: int = 16) -> Fraction:
    a = math.ceil(lo * den)
    b = math.floor(hi * den)
    return Fraction(rng.randint(a, b), den)


def random_sparse_binary(n: int, degree: int, seed: int, coupling_range=(Fraction(1, 2), Fraction(3)),
                         field_range=(Fraction(1, 2), Fraction(2)), perturb_scale: Fraction = Fraction(2)
                         ) -> CoupledModel:
    """Random d-regular binary MRF with mixed-sign couplings; B adds a field at node 0."""
    rng = random.Random(seed)
    edges = random_regular_edges(n, degree, rng)
    names = [f"v{i}" for i in range(n)]
    domains = {v: 2 for v in names}
    shared = []
    for (a, b) in edges:
        qe = random_rational(rng, *coupling_range)
        shared.append(edge_factor(names[a], names[b], qe, f"J{a}_{b}"))
    for i in range(n):
        shared.append(field_factor(names[i], random_rational(rng, *field_range), f"h{i}"))
    only_b = [field_factor(names[0], perturb_scale, "dh0")]
    meta = {"kind": "random_regular", "n": n, "degree": degree, "seed": seed}
    return CoupledModel(tuple(names), domains, dict(domains), tuple(shared), (), tuple(only_b), meta)


def _energy_weight(e: int, base: Fraction) -> Fraction:
    return base ** e


def rotamer_potts(n_residues: int, seed: int, rotamers: tuple[int, int] = (5, 12),
                  mutant_rotamers: int | None = None, mutated: int = 0, layout: str = "chain",
                  contact_prob: float = 0.35, base: Fraction = Fraction(3, 2),
                  energy_range: tuple[int, int] = (-2, 2)) -> CoupledModel:
    """Rotamer-like Potts model on a fixed backbone contact graph.

    Residues sit on a chain; residue i contacts i+1 always and i+2, i+3 with
    probability ``contact_prob`` (a crude backbone contact map).  Each residue
    has 5-20 rotamer states.  Pair and self weights are base^(-E) with small
    random integer energies E.  The "mutation" at ``mutated`` replaces that
    residue's rotamer set (new domain size) and all of its factors in B.
    """
    rng = random.Random(seed)
    sizes = [rng.randint(*rotamers) for _ in range(n_residues)]
    if mutant_rotamers is None:
        mutant_rotamers = rng.randint(*rotamers)
    names = [f"r{i}" for i in range(n_residues)]
    da = {names[i]: sizes[i] for i in range(n_residues)}
    db = dict(da)
    db[names[mutated]] = mutant_rotamers
    contacts = set()
    for i in range(n_residues - 1):
        contacts.add((i, i + 1))
        for k in (2, 3):
            if i + k < n_residues and rng.random() < contact_prob:
                contacts.add((i, i + k))
    lo, hi = energy_range

    def table(d1, d2):
        return tuple(_energy_weight(-rng.randint(lo, hi), base) for _ in range(d1 * d2))

    shared, only_a, only_b = [], [], []
    for i in range(n_residues):
        f_a = Factor((names[i],), table(1, da[names[i]]), f"self{i}")
        if i == mutated:
            only_a.append(f_a)
            only_b.append(Factor((names[i],), table(1, db[names[i]]), f"selfB{i}"))
        else:
            shared.append(f_a)
    for (i, j) in sorted(contacts):
        fa = Factor((names[i], names[j]), table(da[names[i]], da[names[j]]), f"pair{i}_{j}")
        if mutated in (i, j):
            only_a.append(fa)
            only_b.append(Factor((names[i], names[j]), table(db[names[i]], db[names[j]]),
                                 f"pairB{i}_{j}"))
        else:
            shared.append(fa)
    meta = {"kind": "rotamer_potts", "n": n_residues, "seed": seed, "mutated": mutated,
            "contacts": len(contacts)}
    return CoupledModel(tuple(names), da, db, tuple(shared), tuple(only_a), tuple(only_b), meta)


def rotamer_grid_potts(rows: int, cols: int, seed: int, rotamers: tuple[int, int] = (5, 10),
                       mutant_rotamers: int | None = None, base: Fraction = Fraction(3, 2),
                       energy_range: tuple[int, int] = (-1, 1), mutated: tuple[int, int] | None = None
                       ) -> CoupledModel:
    """Rotamer-like Potts model on a rows x cols grid contact graph (packing core)."""
    rng = random.Random(seed)
    names = [site(r, c) for r in range(rows) for c in range(cols)]
    sizes = {v: rng.randint(*rotamers) for v in names}
    mut = site(*(mutated if mutated is not None else (rows // 2, cols // 2)))
    db = dict(sizes)
    db[mut] = mutant_rotamers if mutant_rotamers is not None else rng.randint(*rotamers)
    lo, hi = energy_range

    def table(n):
        return tuple(base ** (-rng.randint(lo, hi)) for _ in range(n))

    shared, only_a, only_b = [], [], []
    for v in names:
        f = Factor((v,), table(sizes[v]), f"self_{v}")
        if v == mut:
            only_a.append(f)
            only_b.append(Factor((v,), table(db[v]), f"selfB_{v}"))
        else:
            shared.append(f)
    for (a, b) in grid_edges(rows, cols):
        u, w = site(*a), site(*b)
        f = Factor((u, w), table(sizes[u] * sizes[w]), f"pair_{u}_{w}")
        if mut in (u, w):
            only_a.append(f)
            only_b.append(Factor((u, w), table(db[u] * db[w]), f"pairB_{u}_{w}"))
        else:
            shared.append(f)
    meta = {"kind": "rotamer_grid", "rows": rows, "cols": cols, "seed": seed, "mutated": mut}
    return CoupledModel(tuple(names), sizes, db, tuple(shared), tuple(only_a), tuple(only_b), meta)


def pairwise_binary(model: CoupledModel, which: str, variables: Sequence[str] | None = None):
    """Extract (unary[(phi0, phi1)], edges{(i, j): 2x2 table}) for a binary pairwise system.

    Multiple factors on the same scope are multiplied.  Raises if a factor has
    more than two variables or a variable is not binary.
    """
    view = model.system(which)
    variables = list(model.variables if variables is None else variables)
    index = {v: i for i, v in enumerate(variables)}
    unary = [[Fraction(1), Fraction(1)] for _ in variables]
    edges: dict[tuple[int, int], list[Fraction]] = {}
    const = Fraction(1)
    for v in variables:
        if view.domains[v] != 2:
            raise ValueError("pairwise_binary needs binary variables")
    for f in view.factors:
        if len(f.scope) == 0:
            const *= f.table[0]
        elif len(f.scope) == 1:
            i = index[f.scope[0]]
            unary[i][0] *= f.table[0]
            unary[i][1] *= f.table[1]
        elif len(f.scope) == 2:
            i, j = index[f.scope[0]], index[f.scope[1]]
            t = list(f.table)
            if i > j:
                i, j = j, i
                t = [t[0], t[2], t[1], t[3]]
            cur = edges.setdefault((i, j), [Fraction(1)] * 4)
            for k in range(4):
                cur[k] *= t[k]
        else:
            raise ValueError("pairwise_binary supports unary and pairwise factors only")
    return const, [tuple(u) for u in unary], {e: tuple(t) for e, t in edges.items()}
