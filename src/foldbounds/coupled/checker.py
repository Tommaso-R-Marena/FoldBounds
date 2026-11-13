"""Independent minimal checker for HP point-mutation T1 certificates.

Imports only ``fractions`` and ``json``. It shares NO code with the producer,
``ratio.py`` or ``foldbounds.lattice``. It fails closed: any malformed input,
any exception, or any unproven claim returns ``(False, reason)``.

Certificate (JSON object)::

  {"seqA": "HPPH...", "seqB": "HPPP...",   # differ at exactly one index
   "q": "3/2",                              # exact rational >= 1
   "i": 3,                                  # the mutated index
   "order": [0,1,2,...],                    # optional growth order (default identity)
   "frontier": [{"w": "RUL...", "wit": "RUL...optional full completion"}, ...],
   "claim": {"lo": "p/q", "hi": "p/q"}}     # claimed bounds on Z_B / Z_A

The growth order pi is a permutation of 0..n-1 such that every prefix of pi is a
contiguous chain segment. Letter t of a word places site pi[t+1] one step
(R=(1,0) U=(0,1) L=(-1,0) D=(0,-1)) from its already placed chain neighbour;
site pi[0] sits at (0,0) and the first letter must be R (this fixes the rotation;
reflections are not quotiented). Z is a sum over conformations modulo translation
and rotation, one representative per C4-orbit; for ANY order the representatives
differ but weights are invariant, so the ratio Z_B/Z_A does not depend on the
order. A word with n-1 letters is a complete walk. Z = sum_x q^{contacts(x)} over
complete SAWs x; contacts = nonconsecutive H-H lattice neighbours.

Per frontier node (prefix with m sites, remaining r = n-m steps):
 * complete (m = n): alpha = q^kA, beta = q^kB exactly, r = R = beta/alpha;
 * dead (m < n, no self-avoiding extension): alpha = 0;
 * otherwise U = 3^r q^(known + Ucap) with Ucap = sum over future H sites j
   of min(#H partners a <= j-3 of opposite chain parity, 2 (3 if j = n-1)),
   L = q^kA(witness) if a valid completion is supplied, else 0.
   Sound: a contact is charged to its later endpoint j; j has 4 lattice
   neighbours, one is j-1, one is j+1 unless j is last; opposite parity is
   forced by the checkerboard. There are at most 3^r completions.
 * ratio: Delta = kB - kA = s*c_i with s = -1 for H->P, +1 for P->H and c_i the
   number of H sites non-consecutively adjacent to i. c_i in [cmin, cmax]
   from placed contacts of i, the free lattice cells around i, and parity.
Then T1 (CEC_THEORY section 1): lo is valid iff g-(lo) >= 0, hi iff g+(hi) <= 0,
provided sum U > 0. Coverage: every valid child of every proper prefix in the
frontier trie must be a frontier word or again a proper prefix; no word is a
proper prefix of another; no duplicates.
"""
import json
from fractions import Fraction

STEP = {"R": (1, 0), "U": (0, 1), "L": (-1, 0), "D": (0, -1)}


class Reject(Exception):
    pass


def _frac(x, what):
    if isinstance(x, bool) or not isinstance(x, (str, int)):
        raise Reject(f"{what}: must be an exact rational string or int")
    try:
        return Fraction(x)
    except Exception:
        raise Reject(f"{what}: unparsable")


def _place(word, order, n):
    """Positions {site: (x, y)} placed by `word`; rejects non-self-avoiding words."""
    if not isinstance(word, str) or not (1 <= len(word) <= n - 1) or word[0] != "R":
        raise Reject("bad word length or first step not R")
    pos = {order[0]: (0, 0)}
    seen = {(0, 0)}
    for t, ch in enumerate(word, start=1):
        if ch not in STEP:
            raise Reject("bad letter")
        s = order[t]
        c = s - 1 if (s - 1) in pos else s + 1       # placed chain neighbour
        d = STEP[ch]
        p = (pos[c][0] + d[0], pos[c][1] + d[1])
        if p in seen:
            raise Reject("self-intersecting word " + word)
        seen.add(p)
        pos[s] = p
    return pos


def _children(pos, order, m, n):
    """Letters that extend a placed prefix of m sites by site order[m] self-avoidingly."""
    s = order[m]
    c = s - 1 if (s - 1) in pos else s + 1
    occ = set(pos.values())
    x, y = pos[c]
    return [ch for ch, (dx, dy) in STEP.items() if (x + dx, y + dy) not in occ]


def _pairs(pos, seq):
    """Nonconsecutive adjacent H-H pairs among placed sites, as frozensets."""
    inv = {p: k for k, p in pos.items()}
    out = set()
    for a, p in pos.items():
        if seq[a] != "H":
            continue
        for dx, dy in STEP.values():
            b = inv.get((p[0] + dx, p[1] + dy))
            if b is not None and abs(b - a) >= 2 and seq[b] == "H":
                out.add(frozenset((a, b)))
    return out


def check(cert):
    """Return (accepted, reason). Never raises."""
    try:
        return _check(cert)
    except Reject as e:
        return False, "reject: " + str(e)
    except Exception as e:  # fail closed
        return False, "reject: malformed (" + type(e).__name__ + ")"


def _check(cert):
    if isinstance(cert, (str, bytes)):
        cert = json.loads(cert)
    if not isinstance(cert, dict):
        raise Reject("certificate must be an object")
    A, B, i = cert["seqA"], cert["seqB"], cert["i"]
    if not (isinstance(A, str) and isinstance(B, str) and len(A) == len(B)
            and set(A) <= {"H", "P"} and set(B) <= {"H", "P"}):
        raise Reject("bad sequences")
    n = len(A)
    if n < 3 or type(i) is not int or not (0 <= i < n):
        raise Reject("bad length or mutated index")
    if [k for k in range(n) if A[k] != B[k]] != [i]:
        raise Reject("sequences must differ exactly at index i")
    q = _frac(cert["q"], "q")
    if q < 1:
        raise Reject("q must be >= 1")
    order = cert.get("order", list(range(n)))
    if not (isinstance(order, list) and sorted(order) == list(range(n))
            and all(type(v) is int for v in order)):
        raise Reject("order is not a permutation")
    for t in range(1, n):            # every prefix must be a contiguous segment
        seg = order[:t]
        if order[t] not in (min(seg) - 1, max(seg) + 1):
            raise Reject("order does not grow a contiguous segment")
    claim = cert["claim"]
    lo, hi = _frac(claim["lo"], "lo"), _frac(claim["hi"], "hi")
    sgn = -1 if A[i] == "H" else 1          # Delta = kB - kA = sgn * c_i
    front = cert["frontier"]
    if not isinstance(front, list) or not front:
        raise Reject("empty frontier")
    hpart = [k for k in range(n) if A[k] == "H" and k != i]   # identical in A and B

    # ---- partition: no duplicates, no overlaps, every valid child present
    words = {}
    for nd in front:
        w = nd["w"]
        if w in words:
            raise Reject("duplicate frontier word")
        words[w] = nd
    trie = {w[:t] for w in words for t in range(1, len(w))}
    if trie & set(words):
        raise Reject("frontier word is a proper prefix of another (overlap)")
    for t in trie:
        pos = _place(t, order, n)
        for ch in _children(pos, order, len(t) + 1, n):
            if t + ch not in trie and t + ch not in words:
                raise Reject("coverage gap below " + t)
    if "R" not in trie and "R" not in words:
        raise Reject("root not covered")

    # ---- per-region data (L, U, r, R)
    regs = []
    for w, nd in words.items():
        pos = _place(w, order, n)
        m = len(w) + 1
        if m == n:
            kA, kB = len(_pairs(pos, A)), len(_pairs(pos, B))
            al, be = q ** kA, q ** kB
            regs.append((al, al, be / al, be / al))
            continue
        if not _children(pos, order, m, n):
            regs.append((Fraction(0),) * 4)           # dead end: no completion
            continue
        known = len(_pairs(pos, A))
        placed = set(pos)
        ucap = 0
        for t in range(m, n):
            j = order[t]
            if A[j] == "H":
                before = placed | set(order[m:t])
                nj = sum(1 for a in before if A[a] == "H" and abs(a - j) >= 3 and (a - j) % 2)
                ucap += min(nj, 3 if j in (0, n - 1) else 2)
        U = Fraction(3) ** (n - m) * q ** (known + ucap)
        L = Fraction(0)
        wit = nd.get("wit")
        if wit is not None:
            if not (isinstance(wit, str) and len(wit) == n - 1 and wit.startswith(w)):
                raise Reject("witness must be a full-length extension of the node word")
            L = q ** len(_pairs(_place(wit, order, n), A))
        # ratio bounds: contacts c_i of site i
        if i in pos:
            x, y = pos[i]
            inv = {p: k for k, p in pos.items()}
            cmin = sum(1 for dx, dy in STEP.values()
                       if (b := inv.get((x + dx, y + dy))) is not None and abs(b - i) >= 2 and A[b] == "H")
            free = sum(1 for dx, dy in STEP.values() if (x + dx, y + dy) not in inv)
            free -= sum(1 for c in (i - 1, i + 1) if 0 <= c < n and c not in pos)
            nf = sum(1 for k in hpart if k not in pos and abs(k - i) >= 3 and (k - i) % 2)
            cmax = cmin + max(0, min(nf, free))
        else:
            cmin = 0
            nf = sum(1 for k in hpart if abs(k - i) >= 3 and (k - i) % 2)
            cmax = min(nf, 3 if i in (0, n - 1) else 2)
        d1, d2 = sgn * cmin, sgn * cmax
        regs.append((L, U, q ** min(d1, d2), q ** max(d1, d2)))

    # ---- T1 exact tests (CEC_THEORY 1.1, 1.2)
    if not any(U > 0 for _, U, _, _ in regs):
        raise Reject("empty admissible set")
    if lo > hi:
        raise Reject("claim lo > hi")
    if sum((r - lo) * (L if r >= lo else U) for L, U, r, R in regs) < 0:
        raise Reject("lower claim not implied: g-(lo) < 0")
    if sum((R - hi) * (U if R >= hi else L) for L, U, r, R in regs) > 0:
        raise Reject("upper claim not implied: g+(hi) > 0")
    return True, "accepted"
