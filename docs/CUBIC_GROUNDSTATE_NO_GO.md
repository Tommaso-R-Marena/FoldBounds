# Cubic all-H 3D model: globally optimal folds need not be identifiable

This is an explicit **negative control** for overclaiming high-contact certificates
as unique structure predictions. It uses the classical cubic lattice
edge-isoperimetric theorem and a checked constructive set of independent,
energy-preserving local fold rewires.

It is a mathematical fact about a **simple cubic homopolymer**. It does
not establish a new biological folding law or contradict realistic
sequence-dependent protein energetics. This implementation is
Python-checked; the isoperimetric step is drawn from published
mathematics, NOT newly proved in Lean.

## Global optimum: a classical bound

Let a length-n self-avoiding path occupy an arbitrary set of n distinct
vertices S in Z^3. Let e(S) denote the number of unordered pairs of
occupied sites that are lattice nearest neighbors. Every complete path
has exactly n-1 covalent bonds among those edges, so its count of
noncovalent contacts is

    K = e(S) - (n-1).

For each coordinate direction x,y,z, every nonempty fiber of occupied
sites along that axis contains at most (fiber site count−1) nearest
neighbor edges. If A_x,A_y,A_z are the cardinalities of the three
axis-orthogonal coordinate projections, then

    e(S) <= (n-A_x)+(n-A_y)+(n-A_z)
         = 3n-(A_x+A_y+A_z).

The **discrete Loomis–Whitney inequality** gives

    A_x*A_y*A_z >= n^2.

The arithmetic-geometric mean inequality therefore implies

    A_x+A_y+A_z >= 3*n^(2/3),

and hence the universal 3D cubic lattice contact upper bound

    K <= 2n - 3*n^(2/3) + 1.

For n=L³ this is the exact integer inequality

    K <= 2*L³ - 3*L² + 1.

A completely occupied L×L×L cube has precisely
3L²(L−1) nearest-neighbor edges. Every Hamiltonian path visiting all
L³ vertices thus attains exactly

    K = 3L²(L−1) − (L³−1)
      = 2L³ − 3L² + 1.

So such a path is a **global contact maximizer among ALL length-L³ 3D
self-avoiding conformations**. This is stronger than merely being
the most compact member among a restricted set.

For all-H residues with attractive nearest-neighbor contact factor
q>1, the full-cube Hamiltonian path is consequently a true global
energy minimum (ground state) of the simple lattice HP Hamiltonian.

External mathematical background: see the cubic lattice
edge-isoperimetric problem of
Mainini–Piovano–Schmidt–Stefanelli, *N³/⁴ law in the cubic lattice*
(2019), https://pmc.ncbi.nlm.nih.gov/articles/PMC6733839/ ,
and the general discrete Loomis–Whitney inequality as used by
Strachan–Swanepoel, *Edge Isoperimetry of Lattices* (2026),
https://link.springer.com/article/10.1007/s00026-025-00801-x .

## Globally optimal yet necessarily nonunique

For even L>=4, the checker generates an L³-site Hamiltonian cube walk
using the independent 2×2×2 macro motif library, then **independently**
rechecks every coordinate, occupied site, and nonbonded contact count.
It rejects any mismatch with the sharp global lattice bound above.

The accompanying
\`local_rewire_confidence.py\` checker identifies spatially
nonadjacent, consecutively traversed 2×2×2 subcubes whose interior
can be rearranged while fixing their chain entry/exit coordinates.
Because the chain is all-H and those eight occupied sites remain H,
all rearrangements preserve the total number of HH nonbonded
contacts, regardless of the specific internal Hamiltonian ordering.

If D is the product of the distinct verified choices at the
independently rewirable cubes, then we have constructively witnessed
**at least D globally ground-state conformations** of equal model
energy. No enumeration of the full n-site conformational ensemble
is required.

For an ordinary exact HP Boltzmann ensemble at ANY q>0,

    P(a particular certified Hamiltonian target) <= 1/D.

In particular it holds even in the strongly attractive q→infinity
limit because **all D global maximizers** remain energy-degenerate.
If arbitrary statewise corrections in a positive interval [l,u]
are allowed, the robust ceiling is

    P(target) <= u / [u+l*(D-1)].

If one identifies geometrically congruent rigid folds instead of
index-specific anchored microstates, there are at most 48 cubic
rotations/reflections of any anchored path, hence the conservative
symmetry-quotiented ceiling

    P(rigid fold class) <= 48u/[48u+l*(D-48)],  for D>48.

These are **upper** bounds on a *single structure's confidence*,
not lower bounds on the population of a high-contact event. The
two are compatible: a broad event may capture almost the entire
ensemble while no one of its exponentially many states has much
equilibrium probability.

## Trust boundary and reproducibility

Reproduce:

    python -m pip install -e '.[test]'
    python -m pytest -q tests/test_cubic_groundstate_no_go.py
    python experiments/cubic_groundstate_no_go_demo.py

The independent oracle exhaustively evaluates every one of 13,565
anchored eight-site 3D self-avoiding walks and verifies that the
maximum contact count is exactly **5**, matching the n=2³
isoperimetric bound. For n=4³,6³,8³ it separately computes the
all-pairs contact score of the generated Hamiltonian walk and
checks the exact global bound. It tests arbitrary positive rational
q (including repulsive q<1 where the walk is NOT claimed to be a
ground state) and statewise correction intervals.

**Scientific limitation:** the global bound is a classical theorem
used as an external mathematical premise. The checks do not amount
to a new Lean proof of Loomis–Whitney, nor do they make the all-H
model a biologically realistic protein. Ground-state degeneracy in
HP models is documented prior art; the useful result here is an
independently auditable *no-go confidence certificate* that prevents
a high-contact guarantee from being misreported as a unique-fold
guarantee. A further research step is to apply robust local
reweighting to published heterogeneous 3D lattice benchmark
sequences and, ultimately, to well-calibrated real structural
ensembles with a justified geometric/energy error model.
