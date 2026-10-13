"""Bipartite b-matching upper bound on HP lattice contacts.

This is a relaxation for a finite toy model, not real-protein prediction.
The Python implementation is tested but NOT yet Lean-verified.
"""
from __future__ import annotations

from collections import deque
from .lattice import DIRECTIONS, Path, _contacts
from .isoperimetric import prefix_contact_cap


def _vacant_neighbor_distance(path: Path, index: int, limit: int) -> int:
    """Shortest walk through vacant sites from the tip to a contact site.

    The occupied tip is the starting vertex; placed sites cannot be revisited.
    Relaxing future self-avoidance cannot overestimate the needed distance.
    """
    occupied = set(path)
    x, y = path[index]
    targets = {(x+dx, y+dy) for dx,dy in DIRECTIONS} - occupied
    if not targets:
        return limit+1
    tip = path[-1]
    visited = {tip}
    queue = deque([(tip, 0)])
    while queue:
        (x,y), dist = queue.popleft()
        if dist >= limit:
            continue
        for dx,dy in DIRECTIONS:
            v=(x+dx,y+dy)
            if v in occupied or v in visited:
                continue
            if v in targets:
                return dist+1
            visited.add(v)
            queue.append((v,dist+1))
    return limit+1


def _maximum_bipartite_bmatching(capacity: dict[int,int],
                                  pairs: tuple[tuple[int,int], ...]) -> int:
    """Integral max flow in a parity-bipartite, unit-edge-capacity graph."""
    left_ids = [i for i,c in capacity.items() if i%2 == 0 and c>0]
    right_ids = [i for i,c in capacity.items() if i%2 == 1 and c>0]
    left = {i:j+1 for j,i in enumerate(left_ids)}
    right = {i:j+len(left)+1 for j,i in enumerate(right_ids)}
    src, sink = 0, len(left)+len(right)+1
    graph: list[list[list[int]]] = [[] for _ in range(sink+1)]

    def edge(u: int, v: int, cap: int) -> None:
        graph[u].append([v,cap,len(graph[v])])
        graph[v].append([u,0,len(graph[u])-1])

    for i,node in left.items():
        edge(src,node,capacity[i])
    for i,node in right.items():
        edge(node,sink,capacity[i])
    for i,j in pairs:
        a,b = (i,j) if i%2 == 0 else (j,i)
        if a in left and b in right:
            edge(left[a],right[b],1)
    total=0
    while True:
        dist=[-1]*len(graph)
        dist[src]=0
        queue=deque([src])
        while queue:
            u=queue.popleft()
            for v,c,_ in graph[u]:
                if c>0 and dist[v]<0:
                    dist[v]=dist[u]+1
                    queue.append(v)
        if dist[sink]<0:
            return total
        pointer=[0]*len(graph)

        def augment(u: int, max_inc: int) -> int:
            if u == sink:
                return max_inc
            while pointer[u]<len(graph[u]):
                v,c,reverse = e = graph[u][pointer[u]]
                if c>0 and dist[v]==dist[u]+1:
                    inc=augment(v,min(c,max_inc))
                    if inc:
                        e[1]-=inc
                        graph[v][reverse][1]+=inc
                        return inc
                pointer[u]+=1
            return 0
        while (inc:=augment(src,len(pairs)+1)):
            total+=inc


def matching_score_cap(sequence: str, path: Path,
                       pairs: tuple[tuple[int,int], ...], *,
                       use_perimeter: bool = True) -> tuple[int,int]:
    """Bounds [realized contacts, max contacts] for all prefix completions.

    Every future H-H contact is an edge joining opposite index parities,
    uses one lattice-contact slot at each endpoint, and must satisfy static
    obstacle reachability if one residue is placed. Hence actual future
    contacts form a feasible b-matching in the relaxed graph.
    """
    if type(use_perimeter) is not bool:
        raise ValueError('use_perimeter must be a bool')
    n,m=len(sequence),len(path)
    known=_contacts(path,pairs)
    occupied=set(path)
    caps: dict[int,int]={}
    for i,residue in enumerate(sequence):
        if residue!='H':
            continue
        bound=3 if i in (0,n-1) else 2
        if i<m:
            x,y=path[i]
            free=sum((x+dx,y+dy) not in occupied for dx,dy in DIRECTIONS)
            existing=sum(1 for a,b in pairs if b<m and (i==a or i==b)
                         and abs(path[a][0]-path[b][0])+abs(path[a][1]-path[b][1])==1)
            bound=min(free,max(0,bound-existing))
        caps[i]=bound
    reach={i:_vacant_neighbor_distance(path,i,n-m) for i in caps if i<m}
    candidates=tuple((i,j) for i,j in pairs if j>=m and
                     (i>=m or reach[i]<=j-m+1))
    cap = known+_maximum_bipartite_bmatching(caps,candidates)
    if use_perimeter:
        cap = min(cap,prefix_contact_cap(path,n))
    return known, cap
