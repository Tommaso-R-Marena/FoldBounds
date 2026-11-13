"""Verify: in protonation variables sigma (1=protonated), every pairwise Coulomb term
W*q_i(sigma_i)*q_j(sigma_j) with W>0 has sigma_i*sigma_j coefficient +W (antiferromagnetic),
for acid-acid, acid-base, base-base; q_acid = sigma-1, q_base = sigma."""
import itertools
q={'acid':lambda s: s-1, 'base':lambda s: s}
W=1.0
for ti,tj in itertools.product(q,q):
    E={(a,b): W*q[ti](a)*q[tj](b) for a in (0,1) for b in (0,1)}
    # interaction coefficient J = E11 - E10 - E01 + E00
    J=E[(1,1)]-E[(1,0)]-E[(0,1)]+E[(0,0)]
    print(ti,tj,"E:",E,"sigma_i*sigma_j coefficient:",J)
