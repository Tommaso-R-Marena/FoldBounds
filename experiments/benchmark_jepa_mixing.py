"""Numerical long-run mixing falsification for *exact* finite pivot kernels.

The stochastic matrices and equilibrium weights are generated/verified using
exact Fractions. Poisson solves, truncation checks, IAT and ESS are FLOATING
diagnostics, not mathematical certificates or physical folding rates.

Requires optional scipy and numpy; does not run with FoldBounds base install.
"""
from __future__ import annotations

from fractions import Fraction
import json
import math
import time

from foldbounds.jepa_reversible_kinetics import build_kernel,verify_kernel
from foldbounds.jepa_flux_maxflow import (
    propose_maxflow_certificate,construct_reversible_kernel,
    verify_reversible_kernel,
)

CASES = ("HHHHHH", "HPHHPPHH", "HHPPHHHH", "HPHPHPHH", "HPHHPPHHH")


def _matrix(kernel, rows):
    import numpy as np
    from scipy.sparse import csr_matrix
    n=len(kernel.states)
    pi=np.array([float(w/kernel.partition) for w in kernel.weights],dtype=float)
    pi/=pi.sum()
    ijk=[(i,j,float(p)) for i,row in enumerate(rows) for j,p in row]
    P=csr_matrix((
        [p for _,_,p in ijk],
        ([i for i,_,_ in ijk],[j for _,j,_ in ijk]),
    ),shape=(n,n))
    if abs(float(pi @ P.sum(axis=0).A1)-1)>1e-8:
        raise AssertionError("numerical equilibrium mass normalization suspect")
    if max(abs(P.sum(axis=1).A1-1))>1e-10:
        raise AssertionError("numerical transition row not normalized")
    if max(abs(pi @ P - pi))>1e-10:
        raise AssertionError("numerical stationary distribution mismatch")
    return P,pi


def observable_iat(P,pi,observations,*,max_lags:int=3000):
    """Independent Poisson solve and correlation-series cross-validation."""
    import numpy as np
    from scipy.sparse import eye
    from scipy.sparse.linalg import spsolve
    f=np.asarray(observations,dtype=float)
    f-=np.dot(pi,f)
    var=float(np.dot(pi,f*f))
    if var<=1e-15:
        return None
    n=len(f)
    operator=(eye(n,format="csr")-P).tocsr()
    g=np.zeros(n)
    g[1:]=spsolve(operator[1:,1:].tocsc(),f[1:])
    if not np.isfinite(g).all():
        raise ArithmeticError("Poisson equation singular or ill-conditioned")
    residual=float(np.max(np.abs(operator@g-f)))
    if residual>2e-8:
        raise ArithmeticError(f"Poisson equation residual too large: {residual}")
    iat_poisson=2*float(np.dot(pi,f*g))/var-1
    rho1=float(np.dot(pi,f*(P@f)))/var
    evolved=f.copy()
    rho_sum=0.0
    hits=0
    for lag in range(1,max_lags+1):
        evolved=P@evolved
        rho=float(np.dot(pi,f*evolved))/var
        rho_sum+=rho
        hits=hits+1 if abs(rho)<1e-11 else 0
        if hits>=15:
            break
    else:
        raise ArithmeticError("finite-lag correlation series did not converge")
    iat_series=1+2*rho_sum
    if abs(iat_poisson-iat_series)>max(2e-5,abs(iat_poisson)*5e-5):
        raise ArithmeticError("independent correlation sum disagrees with Poisson solve")
    return {
        "tau_int":round(iat_poisson,6),
        "lag_one_rho":round(rho1,6),
        "series_tau_int":round(iat_series,6),
        "last_lag":lag,
        "poisson_residual":residual,
    }


def radius_gyration_sq(path):
    n=len(path)
    cx=sum(x for x,y in path)/n
    cy=sum(y for x,y in path)/n
    return sum((x-cx)**2+(y-cy)**2 for x,y in path)/n


def main():
    for sequence in CASES:
        start=time.perf_counter()
        base=build_kernel(sequence,Fraction(2))
        verify_kernel(base)
        cert=propose_maxflow_certificate(base,1)
        mixture=construct_reversible_kernel(base,cert,blend=Fraction(1,2))
        verify_reversible_kernel(mixture)
        P0,pi0=_matrix(base,base.transitions)
        P1,pi1=_matrix(base,mixture.transitions)
        observables={
            "contact_event":[int(score>=1) for score in base.scores],
            "contact_count":list(base.scores),
            "radius_of_gyration_squared":[radius_gyration_sq(path)
                                          for path in base.states],
        }
        for name,values in observables.items():
            before=observable_iat(P0,pi0,values)
            after=observable_iat(P1,pi1,values)
            if before is None or after is None:
                continue
            print(json.dumps({
                "sequence":sequence,
                "observable":name,
                "states":len(base.states),
                "uniform_tau_int":before["tau_int"],
                "mixture_tau_int":after["tau_int"],
                "effective_sample_size_ratio":round(
                    before["tau_int"]/after["tau_int"],5,
                ),
                "uniform_lag1":before["lag_one_rho"],
                "mixture_lag1":after["lag_one_rho"],
                "direct_series_independent_check":True,
                "float_numerical_diagnostic":True,
                "physical_rate":False,
                "all_comparison_cases_including_regressions_retained":True,
                "elapsed_seconds":round(time.perf_counter()-start,3),
            },sort_keys=True))


if __name__=="__main__":
    main()
