"""Dependency-free independent Fraction validation of outward integer kernel."""
from fractions import Fraction as F
from motion_fixedpoint_core import D, B, Y_SCALE, hessian_boxes, bound_score


def exact_matrix(coords):
    n=len(coords)
    H=[[F(i == j, 5) for j in range(3*n)] for i in range(3*n)]
    for a in range(n):
        for b in range(a+1,n):
            delta=[F(coords[a][k]-coords[b][k], D) for k in range(3)]
            d2=sum((v*v for v in delta), F(0))
            if d2 > 100:
                continue
            assert d2 > 0
            for i in range(3):
                for j in range(3):
                    v=49*delta[i]*delta[j]/(d2*d2)
                    H[3*a+i][3*a+j]+=v
                    H[3*b+i][3*b+j]+=v
                    H[3*a+i][3*b+j]-=v
                    H[3*b+i][3*a+j]-=v
    return H


def check_entries(coords):
    intervals, count=hessian_boxes(coords)
    exact=exact_matrix(coords)
    for i,row in enumerate(exact):
        for j,elem in enumerate(row):
            lo,hi=intervals[i].get(j,(0,0))
            assert F(lo,B)<=elem<=F(hi,B), (i,j,elem,lo,hi)
    return count


def main():
    coords=[[0,0,0],[8*D,0,0]]
    assert check_entries(coords)==1
    rows,_=hessian_boxes(coords)
    # Exact 3D model's two-residue score is 320/277 at τ=1/5.
    y=[F(-1280,277),F(0),F(0),F(1280,277),F(0),F(0)]
    proposal=[round(float(v)*Y_SCALE) for v in y]
    lo,hi=bound_score(rows,(0,1),coords,proposal)
    assert F(lo,B)<=F(320,277)<=F(hi,B)
    assert F(hi-lo,B)<F(1,10**12)

    other=[[0,0,0],[8*D,0,0],[4*D,5*D,0],
           [13*D,-2*D,D//2]]
    assert check_entries(other)>1
    # Check the independent quadratic response for an arbitrary exact witness.
    rows,_=hessian_boxes(other)
    y=[F(i-5,13) for i in range(12)]
    z=[F(0)]*12
    dv=[F(other[0][k]-other[2][k],D) for k in range(3)]
    for k in range(3):
        z[k]=dv[k]
        z[6+k]=-dv[k]
    # This candidate is intentionally inaccurate; sound interval is wide.
    lo,hi=bound_score(rows,(0,2),other,
                      [round(float(a)*Y_SCALE) for a in y])
    H=exact_matrix(other)
    # Exact inverse solved by Gaussian elimination over Fraction.
    n=len(z)
    A=[row[:]+[z[i]] for i,row in enumerate(H)]
    for k in range(n):
        pivot=next(i for i in range(k,n) if A[i][k])
        A[k],A[pivot]=A[pivot],A[k]
        v=A[k][k]
        A[k]=[t/v for t in A[k]]
        for i in range(n):
            if i!=k and A[i][k]:
                c=A[i][k]
                A[i]=[a-c*b for a,b in zip(A[i],A[k])]
    exact_y=[A[i][-1] for i in range(n)]
    d2=sum(v*v for v in dv)
    score=sum(z[i]*exact_y[i] for i in range(n))/d2
    assert F(lo,B)<=score<=F(hi,B)

    for malformed in ([[0,0,0],[0,0,0]], [[0,0],[2,3]]):
        try:
            hessian_boxes(malformed)
        except ValueError:
            pass
        else:
            raise AssertionError("malformed/colliding input accepted")
    try:
        bound_score(rows,(2,2),other,[0]*12)
    except ValueError:
        pass
    else:
        raise AssertionError("zero-length index pair accepted")
    print("PASS: exact Fraction matrix bounds, independent exact solutions, failed-input rejection")


if __name__ == "__main__":
    main()
