# FOLDBOUNDS-KERNEL 1 — certificate format and checker semantics

Status: normative specification of the text format accepted by the Lean
function `FoldBounds.Kernel.verify` (module `FoldBounds/FrontierKernel.lean`,
runner `foldkernel`). The Lean definitions are authoritative; this document
exists so that independent producers, oracles and differential tests can be
written **without importing the checker**.

## 1. Model (unchanged; `FoldBounds.LatticeModel`)

* Sequence `S` over `{H, P}`, length `n >= 2`; residues are indexed `0 .. n-1`.
* A conformation is a self-avoiding walk of `n` sites on `Z^2` whose first
  two sites are `(0,0)` and `(1,0)`. Rotations are fixed by the first bond;
  reflections are **not** identified.
* Directions: `0 = +x`, `1 = -x`, `2 = +y`, `3 = -y`.
* Contact count `k(c)`: the number of pairs `i < j` with `S[i] = S[j] = H`,
  `j >= i + 2`, and sites `i`, `j` at Manhattan distance 1. (On this lattice
  every such pair automatically has `j >= i + 3` and `j - i` odd.)
* Weight `q^k(c)` for an exact rational `q >= 1`; event `A_t = {k >= t}`;
  `Z = sum_c q^k(c)`, `Z_A`, `Z_notA`, `P(A_t) = Z_A / Z`.

## 2. File format

ASCII text, lines terminated by a single `\n`. The file must consist of
exactly these lines, in this order, with no blank lines, no trailing spaces,
no `\r`, and nothing after the final newline:

```
FOLDBOUNDS-KERNEL 1
sequence <S>
query <q> <t>          (one or more query lines)
tree <T>
end
```

* `<S>`: string over `H`/`P`, length `>= 2`.
* `<q>`: canonical nonnegative rational: either a canonical natural number, or
  `<num>/<den>` with canonical naturals, `den >= 2` and `gcd(num, den) = 1`.
  A canonical natural is `0` or a nonempty digit string without a leading `0`.
  The checker rejects the whole file unless every query has `q >= 1`.
* `<t>`: canonical natural number (the event threshold).
* `<T>`: preorder encoding of the certificate tree over the five bytes
  `N E X B -` (section 3), with no separators.

Every deviation is rejected (fail closed). There are no optional fields,
comments or unknown-field escapes, so "duplicate key" ambiguities cannot
arise: each field occurs at one fixed position.

## 3. Certificate tree

The tree is read against the walk built so far. The root corresponds to the
anchored two-site walk `(0,0),(1,0)`. At a node whose walk has `m` placed
sites, `r = n - m` residues remain.

| byte | meaning | accepted when | contribution |
|---|---|---|---|
| `E` | exact leaf (one complete conformation) | `r = 0` | one count in the exact histogram at `k` |
| `X` | exhaustive subtree: the checker enumerates every completion itself | always | exact histogram counts for all completions |
| `B` | bound leaf: the whole subtree is enclosed by a recomputed cap | always | one bound entry `(a, K, cnt, fam)` (section 4) |
| `N c0 c1 c2 c3` | split on the next direction `d = 0,1,2,3` | `r >= 1` | sum of the children |
| `-` | blocked child marker | only as a child `c_d` of `N`, and exactly when the site reached by direction `d` is already occupied | nothing |

For a child `c_d` of `N`: if the next site is occupied the child **must** be
`-`; otherwise it **must** be a node encoding (`E`, `X`, `B` or `N...`) and is
checked with the walk extended by that site. A dead end is `N----`.
The tree must be consumed exactly (no trailing bytes).

Consequences enforced structurally by the encoding plus these checks: every
legal immediate successor is covered exactly once; no duplicate, overlapping,
ancestor/descendant or omitted branch can be represented; the certificate
carries no numbers at all (no claimed scores, masses, caps or dead flags).

## 4. Bound entries (all recomputed by the checker)

For a bound leaf with `m` placed sites, `r = n - m`:

* `a` = number of H–H contacts among the placed residues (as in section 1,
  restricted to indices `< m`).
* `cand(i, j)` := `i + 3 <= j` and `(j - i)` odd and `S[i] = S[j] = H`.
* `slots(j)` := `2` if `j + 1 < n`, else `3`.
* `futureCap(m)` := `sum_{j = m}^{n-1} min(slots(j), #{i < n : cand(i, j)})`.
* `w`, `h` := width and height of the bounding box of the placed sites
  (`max x - min x + 1`, `max y - min y + 1`).
* `minPerim(N, w, h)` := `min( 2N + 2, min_{c = w}^{N} ( c + max(h, ceil(N / c)) ) )`
  where `ceil(N/c) = (N + c - 1) div c`.
* `geomCap` := `(n + 1) - minPerim(n, w, h)` (truncated at 0).
* `K` := `min(a + futureCap(m), geomCap)`.
* `cnt` := `3^r`.
* `fam`: let `e` be the last placed site.
  * if for some axis direction the end site is the **strict** unique extreme
    of the placed sites (`+x`: every other placed site has smaller `x`; `-x`:
    larger `x`; `+y`: smaller `y`; `-y`: larger `y`), `fam := pdCount(r)` where
    `(f,u,d)_0 = (1,1,1)`, `(f,u,d)_{r+1} = (f+u+d, f+u, f+d)`, `pdCount(r) = f_r`;
  * else if the end site maximises one of the four diagonal potentials
    `(±x) + (±y)` over all placed sites (non-strictly), `fam := 2^r`;
  * else `fam := 0`.

Proved meaning (Lean): every completion `c` of the leaf's walk has
`a <= k(c) <= K`; the number of completions is at most `cnt` and at least
`fam`.

## 5. Evaluation of a query `(q, t)`

With exact histogram `ex[k]` (number of exactly counted conformations with
`k` contacts) and bound entries `(a, K, cnt, fam)`:

```
zLo = sum_k ex[k] q^k + sum_entries            fam * q^a
zHi = sum_k ex[k] q^k + sum_entries            cnt * q^K
eLo = sum_{k >= t} ex[k] q^k + sum_{entries: t <= a} fam * q^a
eHi = sum_{k >= t} ex[k] q^k + sum_{entries: t <= K} cnt * q^K
nLo = sum_{k <  t} ex[k] q^k + sum_{entries: K <  t} fam * q^a
nHi = sum_{k <  t} ex[k] q^k + sum_{entries: a <  t} cnt * q^min(K, t-1)
pLo = 0 if eLo = 0 else eLo / (eLo + nHi)
pHi = 0 if eHi = 0 else (1 if nLo = 0 else eHi / (eHi + nLo))
```

The tree check does not depend on `q` or `t`: one accepted tree yields sound
enclosures for **every** rational `q >= 1` and **every** threshold. All
arithmetic is exact (`Nat`, `Int`, `Rat`); no floating point is involved.

The exact histogram has `4n + 4` slots; a contact count outside it makes the
checker reject (this cannot happen for valid input, it is a fail-closed guard).

## 6. Runner output (unverified printer)

`foldkernel check FILE...` prints one JSON object per file:

```
{"file": "...", "accepted": true, "format": "FOLDBOUNDS-KERNEL 1",
 "sequence": "HPH...", "n": 14, "cert_bytes": 123, "split_nodes": 30,
 "bound_leaves": 12, "enum_leaves": 3, "exact_leaves": 5,
 "exact_count": 1234, "exact_histogram": [..],
 "results": [{"q": "3/2", "threshold": 4,
              "zLo": "a/b", "zHi": "a/b", "eLo": "a/b", "eHi": "a/b",
              "nLo": "a/b", "nHi": "a/b", "pLo": "a/b", "pHi": "a/b"}],
 "verify_ms": 12}
```

or `{"file": "...", "accepted": false, "reason": "<category>", ...}`.
Rationals are printed as `num/den` in lowest terms with `den >= 1`. Exit
status: `0` if every file was accepted, `1` if some file was rejected, `2` on
usage or I/O errors.

## 7. Enumeration implementation

`X` leaves (and any exhaustive residual) are counted by `enumFast`, which
is proved equal to the list enumerator `enumGo` on every self-avoiding
prefix (`FrontierKernelSound.enumFast_eq`). Operationally `enumFast` packs
the placed sites into a `(2n+1)²` grid, threads the array linearly, and
stores the contact histogram in the cells after the grid. Contact and
occupancy tests are then constant time per candidate direction. The
histogram length is `4n + 4`.

## 8. Trust boundary

Proved in Lean (`FoldBounds/FrontierKernelSound.lean`): for every input string,
if `verify` accepts and returns a report, then for the echoed sequence and
each echoed `(q, t)` the eight reported rationals enclose the exact `Z`,
`Z_A`, `Z_notA` and `P(A_t)` of the Lean HP model. The proof goes through
the actual acceptance path (`verify` → `runTree` → `enumFast` = `enumGo`).

Not proved: file reading, the JSON printer, the Lean compiler and runtime,
and that the echoed parameters are the ones the user intended (compare the
echo). The certificate producer is untrusted and needs no trust.
