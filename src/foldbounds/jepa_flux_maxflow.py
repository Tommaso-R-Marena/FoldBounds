"""Primal-dual *exact rational* stationary event-flux max-flow certificates.

A complete finite HP pivot-move graph is the trust-domain assumption.
All theorems are within a chosen reversible Monte Carlo model, not
physical folding rates or a new mathematical max-flow theorem.
"""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from fractions import Fraction

from .jepa_reversible_kinetics import Kernel, event_flux, verify_kernel


@dataclass(frozen=True)
class FluxCertificate:
    threshold: int
    flows: tuple[tuple[int,int,Fraction], ...]
    source_side_states: tuple[int, ...]
    total_mass: Fraction
    cut_mass: Fraction


@dataclass(frozen=True)
class ReversibleKernel:
    source: Kernel
    transitions: tuple[tuple[tuple[int,Fraction], ...], ...]
    blend: Fraction
    certificate: FluxCertificate


def _support(kernel: Kernel, threshold: int):
    if type(threshold) is not int:
        raise ValueError("threshold must be an integer")
    a = {i for i,k in enumerate(kernel.scores) if k>=threshold}
    b = set(range(len(kernel.states)))-a
    cross = {(i,j) for i in a for j,p in kernel.transitions[i] if j in b and p>0}
    return a,b,cross


def _maxflow(capacity, neighbors, source, sink):
    """Deterministic Edmonds-Karp on Fraction residual capacities."""
    residual=dict(capacity)
    total=Fraction()
    while True:
        parent={source:None}
        queue=deque([source])
        while queue and sink not in parent:
            u=queue.popleft()
            for v in neighbors[u]:
                if v not in parent and residual[(u,v)]>0:
                    parent[v]=u
                    queue.append(v)
                    if v==sink:
                        break
        if sink not in parent:
            break
        edges=[]
        node=sink
        while node!=source:
            p=parent[node]
            edges.append((p,node))
            node=p
        quantity=min(residual[e] for e in edges)
        for u,v in edges:
            residual[(u,v)]-=quantity
            residual[(v,u)]+=quantity
        total+=quantity
    reachable={source}
    queue=deque([source])
    while queue:
        u=queue.popleft()
        for v in neighbors[u]:
            if v not in reachable and residual[(u,v)]>0:
                reachable.add(v)
                queue.append(v)
    return total,residual,reachable


def propose_maxflow_certificate(kernel: Kernel, threshold: int) -> FluxCertificate:
    """Untrusted exact maxflow producer; verifier never trusts its search."""
    verify_kernel(kernel)
    a,b,cross=_support(kernel,threshold)
    n=len(kernel.states)
    s,t=n,n+1
    graph=[set() for _ in range(n+2)]
    capacity: dict[tuple[int,int],Fraction]={}

    def add_edge(u,v,amount):
        if (u,v) not in capacity:
            capacity[(u,v)]=Fraction()
            capacity[(v,u)]=Fraction()
        capacity[(u,v)]+=amount
        graph[u].add(v)
        graph[v].add(u)

    for i in sorted(a):
        add_edge(s,i,kernel.weights[i])
    # The cut-edge capacity strictly exceeds all available stationary mass.
    for i,j in sorted(cross):
        add_edge(i,j,kernel.partition+1)
    for j in sorted(b):
        add_edge(j,t,kernel.weights[j])
    total,residual,reachable=_maxflow(
        capacity,tuple(tuple(sorted(links)) for links in graph),s,t,
    )
    flows=tuple(sorted(
        (i,j,capacity[(i,j)]-residual[(i,j)])
        for i,j in cross if capacity[(i,j)]>residual[(i,j)]
    ))
    source_side=tuple(sorted(reachable-{s,t}))
    cut=(
        sum((kernel.weights[i] for i in a if i not in reachable),Fraction())
        +sum((kernel.weights[j] for j in b if j in reachable),Fraction())
        +sum((kernel.partition+1 for i,j in cross
              if i in reachable and j not in reachable),Fraction())
    )
    certificate=FluxCertificate(threshold,flows,source_side,total,cut)
    verify_maxflow_certificate(kernel,certificate)
    return certificate


def verify_maxflow_certificate(kernel: Kernel, cert: FluxCertificate) -> None:
    """Independent exact primal feasible-flow and dual-cut equality check.

    A feasible stationary cut-flow is a lower bound on optimum. Any
    source/sink cut gives an upper bound. Their exact equality certifies
    global optimality over reversible kernels with the same edge support.
    The verifier does not run, trust or reconstruct the max-flow search.
    """
    verify_kernel(kernel)
    if not isinstance(cert,FluxCertificate) or type(cert.threshold) is not int:
        raise ValueError("invalid kinetic flux certificate")
    a,b,cross=_support(kernel,cert.threshold)
    flow_a={i:Fraction() for i in a}
    flow_b={j:Fraction() for j in b}
    visited=set()
    value=Fraction()
    for i,j,amount in cert.flows:
        if ((i,j) in visited or (i,j) not in cross or
            type(amount) is not Fraction or amount<=0):
            raise ValueError("duplicate, unsupported or invalid cross-cut flow")
        visited.add((i,j))
        flow_a[i]+=amount
        flow_b[j]+=amount
        value+=amount
    if (any(flow_a[i]>kernel.weights[i] for i in a)
            or any(flow_b[j]>kernel.weights[j] for j in b)):
        raise ValueError("stationary vertex capacity exceeded")
    source_side=set(cert.source_side_states)
    if (len(source_side)!=len(cert.source_side_states)
            or not source_side<=(a|b)):
        raise ValueError("invalid max-flow dual cut vertices")
    cut=(
        sum((kernel.weights[i] for i in a if i not in source_side),Fraction())
        +sum((kernel.weights[j] for j in b if j in source_side),Fraction())
        +sum((kernel.partition+1 for i,j in cross
              if i in source_side and j not in source_side),Fraction())
    )
    if (type(cert.total_mass) is not Fraction
            or type(cert.cut_mass) is not Fraction
            or value!=cert.total_mass or cut!=cert.cut_mass
            or value!=cut):
        raise ValueError("exact feasible flow != exact min-cut certificate")


def construct_reversible_kernel(
    kernel: Kernel, certificate: FluxCertificate, *,
    blend: Fraction=Fraction(1,2),
) -> ReversibleKernel:
    """Combine a fully supported reference MC kernel with optimal cut-flows.

    For 0<blend<1, every positive baseline transition remains positive.
    Conditional event flux improves if the optimal cut flux is higher.
    This does NOT prove faster mixing, shorter MFPT or protein kinetics.
    """
    verify_maxflow_certificate(kernel,certificate)
    if (type(blend) is not Fraction or not 0<blend<1):
        raise ValueError("blend must be exact Fraction strictly between 0 and 1")
    rows=[{} for _ in kernel.states]
    for i,row in enumerate(kernel.transitions):
        for j,p in row:
            if j!=i:
                rows[i][j]=(1-blend)*p
    for i,j,mass in certificate.flows:
        rows[i][j]=rows[i].get(j,Fraction())+blend*mass/kernel.weights[i]
        rows[j][i]=rows[j].get(i,Fraction())+blend*mass/kernel.weights[j]
    for i,row in enumerate(rows):
        out=sum(row.values(),Fraction())
        if out>1:
            raise AssertionError("stationary cut-flow violated row capacity")
        row[i]=1-out
    result=ReversibleKernel(
        kernel,
        tuple(tuple(sorted((j,p) for j,p in row.items() if p)) for row in rows),
        blend,certificate,
    )
    verify_reversible_kernel(result)
    baseline_flux=event_flux(kernel,certificate.threshold)[1]
    improved_flux=event_flux_reversible(result,certificate.threshold)
    if improved_flux<baseline_flux:
        raise AssertionError("checked optimum unexpectedly worsened event crossing")
    return result


def verify_reversible_kernel(candidate: ReversibleKernel) -> None:
    """Recompute every row and exact detailed balance independently."""
    if not isinstance(candidate,ReversibleKernel):
        raise ValueError("invalid reversible Markov result")
    verify_maxflow_certificate(candidate.source,candidate.certificate)
    kernel=candidate.source
    alpha=candidate.blend
    if type(alpha) is not Fraction or not 0<alpha<1:
        raise ValueError("invalid exact mixing parameter")
    rows=[{} for _ in kernel.states]
    for i,row in enumerate(kernel.transitions):
        for j,p in row:
            if i!=j:
                rows[i][j]=(1-alpha)*p
    for i,j,mass in candidate.certificate.flows:
        rows[i][j]=rows[i].get(j,Fraction())+alpha*mass/kernel.weights[i]
        rows[j][i]=rows[j].get(i,Fraction())+alpha*mass/kernel.weights[j]
    for i,row in enumerate(rows):
        out=sum(row.values(),Fraction())
        if out>1:
            raise ValueError("negative self-loop after mixing")
        row[i]=1-out
    expected=tuple(tuple(sorted((j,p) for j,p in row.items() if p)) for row in rows)
    if candidate.transitions!=expected:
        raise ValueError("tampered reversible Markov transition matrix")
    row_maps=[dict(row) for row in candidate.transitions]
    for i,row in enumerate(row_maps):
        if (sum(row.values(),Fraction())!=1
                or any(p<0 for p in row.values())):
            raise ValueError("non-stochastic reversible chain")
        for j,p in row.items():
            if kernel.weights[i]*p!=kernel.weights[j]*row_maps[j].get(i,Fraction()):
                raise ValueError("detailed balance violated by proposed flow")


def event_flux_reversible(candidate: ReversibleKernel,threshold: int) -> Fraction:
    """Exact stationary one-step crossing probability for mixed chain."""
    verify_reversible_kernel(candidate)
    kernel=candidate.source
    return sum((
        kernel.weights[i]*p
        for i,row in enumerate(candidate.transitions)
        if kernel.scores[i]>=threshold
        for j,p in row
        if kernel.scores[j]<threshold
    ),Fraction())/kernel.partition



def stationary_event_lag1(kernel: Kernel, threshold: int) -> Fraction | None:
    """Exact stationary lag-one correlation of the binary event indicator.

    For a stationary chain with event probability p and one-step outward
    event flux Phi, Cov(1_A(X_t),1_A(X_{t+1})) = p(1-p)-Phi.
    Hence rho_1 = 1-Phi/[p(1-p)]. Returns None for constant events.
    It does NOT bound integrated autocorrelation, asymptotic variance,
    mixing time, physical kinetics or sample efficiency after burn-in.
    """
    verify_kernel(kernel)
    p,flux,_=event_flux(kernel,threshold)
    if p==0 or p==1:
        return None
    return 1-flux/(p*(1-p))


def stationary_event_lag1_reversible(
    candidate: ReversibleKernel, threshold: int,
) -> Fraction | None:
    """Exact stationary lag-one event autocorrelation for a checked mixture."""
    verify_reversible_kernel(candidate)
    p,_,_=event_flux(candidate.source,threshold)
    if p==0 or p==1:
        return None
    flux=event_flux_reversible(candidate,threshold)
    return 1-flux/(p*(1-p))
