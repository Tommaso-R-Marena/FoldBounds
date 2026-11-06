"""Source-restraint-only rigidity obstructions with constructive component motions.

For N observed atom positions and M distinct distance edges, the 3D
rigidity Jacobian has rank <= M and hence nullity >= 3N-M. After removing
at most six rigid motions there are >= max(0,3N-M-6) internal linearized
zero modes. This is a FIRST-ORDER bound, not an integrable flex claim.

For C components, independent translations of components are exact
distance-preserving motions with at least 3(C-1) relative parameters.
Only included atom-pair NOEs are modeled: no chemistry or excluded data.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence

from .noe_restraint_audit import NOERestraint, verify_record


@dataclass(frozen=True)
class NOERigidityNoGo:
    atoms: int
    unique_pair_edges: int
    components: int
    component_sizes: tuple[int, ...]
    independent_cycle_count: int
    max_rigidity_jacobian_rank: int
    minimum_linearized_internal_zero_modes: int
    exact_component_translation_parameters: int
    long_range_restraints: int
    single_residue_restraints: int


def evaluate_rigidity_gap(records: Sequence[NOERestraint], *,
                          long_range_sep: int = 5) -> NOERigidityNoGo:
    if type(long_range_sep) is not int or long_range_sep < 1:
        raise ValueError("invalid long-range sequence separation")
    rows = tuple(records)
    if not rows:
        raise ValueError("at least one restraint required")
    positions, edges, seen = {}, set(), set()
    long = single = 0
    for r in rows:
        if not verify_record(r) or r.identity in seen:
            raise ValueError("invalid/duplicate restraint ID")
        seen.add(r.identity)
        for atom, p in ((r.left,r.left_position_mA),
                        (r.right,r.right_position_mA)):
            if atom in positions and positions[atom] != p:
                raise ValueError("inconsistent coordinates for repeated atom")
            positions[atom] = p
        edge = frozenset((r.left,r.right))
        if len(edge) != 2 or edge in edges:
            raise ValueError("duplicate observed atom pair")
        edges.add(edge)
        long += int(abs(r.left.residue-r.right.residue)>=long_range_sep)
        single += int(r.left.residue == r.right.residue)
    parent = {atom:atom for atom in positions}
    rank = {atom:0 for atom in positions}
    def find(a):
        if parent[a] != a:
            parent[a] = find(parent[a])
        return parent[a]
    for edge in edges:
        u,v = tuple(edge)
        u,v = find(u),find(v)
        if u == v: continue
        if rank[u] < rank[v]:
            u,v=v,u
        parent[v] = u
        if rank[u] == rank[v]:
            rank[u] += 1
    sizes={}
    for atom in positions:
        root=find(atom)
        sizes[root]=sizes.get(root,0)+1
    n,m,c=len(positions),len(edges),len(sizes)
    return NOERigidityNoGo(
        n,m,c,tuple(sorted(sizes.values(),reverse=True)),
        m-n+c,m,max(0,3*n-m-6) if n>=3 else 0,
        3*max(0,c-1),long,single)


def verify_component_translation_witness(records: Sequence[NOERestraint],
                                         data: NOERigidityNoGo) -> bool:
    """Independently replay a nontrivial translated embedding."""
    if type(data) is not NOERigidityNoGo:
        return False
    try:
        expected=evaluate_rigidity_gap(records)
        if expected != data:
            return False
        graph,positions={},{}
        for r in records:
            graph.setdefault(r.left,set()).add(r.right)
            graph.setdefault(r.right,set()).add(r.left)
            positions[r.left]=r.left_position_mA
            positions[r.right]=r.right_position_mA
        groups,visited=[],set()
        for atom in sorted(graph,key=lambda x:(x.residue,x.atom)):
            if atom in visited: continue
            group,stack=set(),[atom]
            while stack:
                v=stack.pop()
                if v in visited: continue
                visited.add(v); group.add(v)
                stack.extend(graph[v]-visited)
            groups.append(group)
        assignment={v:i for i,group in enumerate(groups) for v in group}
        if len(groups)!=data.components or len(assignment)!=data.atoms:
            return False
        for r in records:
            u,v=positions[r.left],positions[r.right]
            a,b=assignment[r.left],assignment[r.right]
            pu=tuple(x+a*(i+1) for i,x in enumerate(u))
            pv=tuple(x+b*(i+1) for i,x in enumerate(v))
            if sum((u[i]-v[i])**2 for i in range(3)) != sum(
                    (pu[i]-pv[i])**2 for i in range(3)):
                return False
        return True
    except (ValueError,TypeError,AttributeError,IndexError,RecursionError):
        return False
