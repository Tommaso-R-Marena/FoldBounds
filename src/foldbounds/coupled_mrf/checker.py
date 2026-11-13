"""Independent minimal checker for coupled-MRF separator certificates.

Imports ONLY ``fractions`` and ``json``; shares no code with the producer.
It re-parses the model, re-derives the canonical hash (pure-Python SHA-256),
re-verifies the separator conditions (S1)-(S3), recomputes Phi_A(tau),
Phi_B(tau) for EVERY boundary state by its own variable elimination in exact
rationals, and compares the certificate's claim with the recomputed
[min rho, max rho].  It fails closed: any malformed input, any exception, or
any mismatch yields ``(False, reasons)``.

Accepted claim: ``lower <= min rho`` and ``upper >= max rho`` (``upper`` may be
the string "inf").  If ``tight`` is true, equality is required.  If some tau has
Phi_A = 0 < Phi_B the upper end must be "inf".  Soundness statement
(CEC_THEORY.md Thm 2.1): for every exterior with Z_A > 0 the true Z_B/Z_A lies
in [lower, upper].
"""
import json
from fractions import Fraction

_FMT = "foldbounds.coupled_mrf.certificate"
_MFMT = "foldbounds.coupled_mrf"


class Reject(Exception):
    pass


def _need(cond, msg):
    if not cond:
        raise Reject(msg)


# ------------------------------------------------------------ pure SHA-256
def _primes(n):
    out, k = [], 2
    while len(out) < n:
        if all(k % p for p in out if p * p <= k):
            out.append(k)
        k += 1
    return out


def _iroot(n, k):
    lo, hi = 0, 1
    while hi ** k <= n:
        hi *= 2
    while lo < hi - 1:
        mid = (lo + hi) // 2
        if mid ** k <= n:
            lo = mid
        else:
            hi = mid
    return lo


def _consts():
    ps = _primes(64)
    # frac(p^(1/3)) * 2^32 = floor(cbrt(p * 2^96)) mod 2^32
    k = [_iroot(p << 96, 3) & 0xFFFFFFFF for p in ps]
    h = [_iroot(p << 64, 2) & 0xFFFFFFFF for p in ps[:8]]
    return k, h


def sha256_hex(data):
    k, h = _consts()
    m = 0xFFFFFFFF
    rr = lambda x, n: ((x >> n) | (x << (32 - n))) & m
    msg = bytearray(data)
    bitlen = len(msg) * 8
    msg.append(0x80)
    while len(msg) % 64 != 56:
        msg.append(0)
    msg += bitlen.to_bytes(8, "big")
    for off in range(0, len(msg), 64):
        w = [int.from_bytes(msg[off + 4 * i: off + 4 * i + 4], "big") for i in range(16)]
        for i in range(16, 64):
            s0 = rr(w[i - 15], 7) ^ rr(w[i - 15], 18) ^ (w[i - 15] >> 3)
            s1 = rr(w[i - 2], 17) ^ rr(w[i - 2], 19) ^ (w[i - 2] >> 10)
            w.append((w[i - 16] + s0 + w[i - 7] + s1) & m)
        a, b, c, d, e, f, g, hh = h
        for i in range(64):
            S1 = rr(e, 6) ^ rr(e, 11) ^ rr(e, 25)
            ch = (e & f) ^ (~e & m & g)
            t1 = (hh + S1 + ch + k[i] + w[i]) & m
            S0 = rr(a, 2) ^ rr(a, 13) ^ rr(a, 22)
            maj = (a & b) ^ (a & c) ^ (b & c)
            t2 = (S0 + maj) & m
            hh, g, f, e, d, c, b, a = g, f, e, (d + t1) & m, c, b, a, (t1 + t2) & m
        h = [(x + y) & m for x, y in zip(h, [a, b, c, d, e, f, g, hh])]
    return "".join(f"{x:08x}" for x in h)


# ----------------------------------------------------------- model parsing
def _num(x):
    _need(not isinstance(x, (bool, float)), "non-exact number")
    if isinstance(x, int):
        return Fraction(x)
    _need(isinstance(x, str), "number must be int or string")
    t = x.strip()
    _need(t != "" and not any(c in t for c in "eEjJnNiI."), "not a rational literal")
    try:
        v = Fraction(t)
    except Exception:
        raise Reject("bad rational literal")
    return v


def _fstr(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def _parse_model(md):
    _need(isinstance(md, dict), "model not an object")
    _need(md.get("format") == _MFMT and md.get("version") == 1, "bad model format")
    # Fail closed on unknown keys (adversarial review m3). A misspelt factor list such as
    # "only_b" would otherwise be ignored silently, while the hash, computed only over the
    # canonical known fields, stayed unchanged.
    _need(set(md) <= {"format", "version", "variables", "shared", "only_A", "only_B", "metadata"},
          "unknown model keys")
    _need(isinstance(md.get("variables"), list), "missing variables list")
    for key in ("shared", "only_A", "only_B"):
        _need(isinstance(md.get(key, []), list), "factor list is not a list")
    names, da, db = [], {}, {}
    for e in md.get("variables", []):
        _need(isinstance(e, dict) and isinstance(e.get("name"), str), "bad variable")
        n = e["name"]
        ks = set(e) - {"name"}
        if ks == {"domain"}:
            a = b = e["domain"]
        else:
            _need(ks == {"domain_A", "domain_B"}, "bad domain keys")
            a, b = e["domain_A"], e["domain_B"]
        for s in (a, b):
            _need(type(s) is int and s >= 1, "bad domain size")
        names.append(n)
        da[n], db[n] = a, b
    _need(len(set(names)) == len(names), "repeated variable")
    out = {}
    for key, doms in (("shared", None), ("only_A", da), ("only_B", db)):
        fl = []
        for f in md.get(key, []):
            _need(isinstance(f, dict), "bad factor")
            _need(set(f) <= {"scope", "table", "name"}, "unknown factor keys")
            sc, tb = f.get("scope"), f.get("table")
            _need(isinstance(sc, list) and isinstance(tb, list), "bad factor")
            _need(all(isinstance(v, str) and v in da for v in sc) and len(set(sc)) == len(sc),
                  "bad factor scope")
            tb = [_num(t) for t in tb]
            _need(all(t >= 0 for t in tb), "negative entry")
            if doms is None:
                _need(all(da[v] == db[v] for v in sc), "shared factor on changed domain")
                dd = da
            else:
                dd = doms
            n = 1
            for v in sc:
                n *= dd[v]
            _need(len(tb) == n, "table size mismatch")
            fl.append((tuple(sc), tb, tuple(dd[v] for v in sc)))
        out[key] = fl
    return names, da, db, out


def _canonical(md, names, da, db, fl):
    def fd(f):
        return {"scope": list(f[0]), "table": [_fstr(t) for t in f[1]]}
    vs = []
    for v in names:
        vs.append({"name": v, "domain": da[v]} if da[v] == db[v]
                  else {"name": v, "domain_A": da[v], "domain_B": db[v]})
    d = {"format": _MFMT, "version": 1, "variables": vs,
         "shared": [fd(f) for f in fl["shared"]],
         "only_A": [fd(f) for f in fl["only_A"]],
         "only_B": [fd(f) for f in fl["only_B"]]}
    return json.dumps(d, sort_keys=True, separators=(",", ":"))


# ------------------------------------------------------------- computation
def _phi(factors, dom, interior, boundary, order):
    """Exact Phi over all boundary states: dict-free list in row-major order."""
    pool = [(sc, tb, ds) for sc, tb, ds in factors]
    # dict tables: scope tuple -> {assignment tuple: value}
    tabs = []
    for sc, tb, ds in factors:
        t = {}
        idx = [0] * len(sc)
        for k in range(len(tb)):
            t[tuple(idx)] = tb[k]
            j = len(sc) - 1
            while j >= 0:
                idx[j] += 1
                if idx[j] < ds[j]:
                    break
                idx[j] = 0
                j -= 1
        tabs.append((sc, t))
    for v in order:
        touch = [x for x in tabs if v in x[0]]
        tabs = [x for x in tabs if v not in x[0]]
        if not touch:
            tabs.append(((), {(): Fraction(dom[v])}))
            continue
        sc = []
        for s, _ in touch:
            for w in s:
                if w != v and w not in sc:
                    sc.append(w)
        full = sc + [v]
        pos = [[full.index(w) for w in s] for s, _ in touch]
        out = {}
        for x in _prod_states([dom[w] for w in full]):
            val = Fraction(1)
            for p, (_, t) in zip(pos, touch):
                val *= t.get(tuple(x[i] for i in p), 0)
                if val == 0:
                    break
            if val:
                key = x[:-1]
                out[key] = out.get(key, 0) + val
        tabs.append((tuple(sc), out))
    res = []
    for tau in _prod_states([dom[b] for b in boundary]):
        asg = dict(zip(boundary, tau))
        val = Fraction(1)
        for s, t in tabs:
            val *= t.get(tuple(asg[w] for w in s), 0)
            if val == 0:
                break
        res.append(val)
    return res


def _prod_states(dims):
    st = [()]
    for d in dims:
        st = [s + (x,) for s in st for x in range(d)]
    return st


def _order(factors, interior, hint):
    inter = set(interior)
    if hint is not None:
        if isinstance(hint, list) and sorted(hint) == sorted(inter) and len(set(hint)) == len(hint):
            return list(hint)
    adj = {v: set() for v in inter}
    for sc, _, _ in factors:
        ins = [v for v in sc if v in inter]
        for a in ins:
            adj[a].update(b for b in ins if b != a)
    order, rem = [], set(inter)
    while rem:
        v = min(rem, key=lambda x: (len(adj[x]), x))
        nb = adj.pop(v)
        for w in nb:
            adj[w].discard(v)
            adj[w] |= nb - {w}
        rem.discard(v)
        order.append(v)
    return order


def _fmtstate(s):
    return None if s is None else list(s)


def check_certificate(model_dict, cert, max_boundary_states=2 ** 18, expected_sha256=None):
    """Return (ok, reasons).  Never raises."""
    try:
        return _check(model_dict, cert, max_boundary_states, expected_sha256)
    except Reject as e:
        return False, [str(e)]
    except BaseException as e:  # fail closed on anything unexpected
        return False, [f"checker error: {type(e).__name__}: {e}"]


def _check(md, cert, max_states, expected):
    _need(isinstance(cert, dict), "certificate not an object")
    allowed = {"format", "version", "model_sha256", "separator", "method", "claim", "witness",
               "hints", "exhaustive", "note"}
    _need(set(cert) <= allowed, "unknown certificate keys")
    _need(cert.get("format") == _FMT and cert.get("version") == 1, "bad certificate format/version")
    names, da, db, fl = _parse_model(md)
    h = sha256_hex(_canonical(md, names, da, db, fl).encode("utf-8"))
    _need(cert.get("model_sha256") == h, "model hash mismatch")
    if expected is not None:
        _need(expected == h, "model hash differs from the expected hash")
    sp = cert.get("separator")
    _need(isinstance(sp, dict) and set(sp) == {"interior", "boundary", "exterior"}, "bad separator")
    I, B, E = sp["interior"], sp["boundary"], sp["exterior"]
    for part in (I, B, E):
        _need(isinstance(part, list) and all(isinstance(v, str) for v in part), "bad separator part")
    allv = set(names)
    _need(len(I) + len(B) + len(E) == len(allv) and set(I) | set(B) | set(E) == allv,
          "separator is not a partition of the variables")
    sI, sB, sE = set(I), set(B), set(E)
    _need(len(sI) == len(I) and len(sB) == len(B) and len(sE) == len(E), "repeated variable in part")
    allf = fl["shared"] + fl["only_A"] + fl["only_B"]
    for sc, _, _ in allf:
        s = set(sc)
        _need(not (s & sI and s & sE), "a factor touches interior and exterior")
    for sc, _, _ in fl["only_A"] + fl["only_B"]:
        _need(set(sc) <= sI | sB, "system-specific factor leaves interior+boundary")
    for v in names:
        _need(da[v] == db[v] or v in sI, "changed-domain variable not interior")
    n = 1
    for v in B:
        n *= da[v]
    _need(n <= max_states, "too many boundary states for the checker")
    loc = sI | sB
    la = [f for f in fl["shared"] + fl["only_A"] if set(f[0]) <= loc]
    lb = [f for f in fl["shared"] + fl["only_B"] if set(f[0]) <= loc]
    hint = (cert.get("hints") or {}).get("elimination_order") if isinstance(cert.get("hints"), dict) else None
    pa = _phi(la, da, sI, B, _order(la, I, hint))
    pb = _phi(lb, db, sI, B, _order(lb, I, hint))
    states = _prod_states([da[v] for v in B])
    lo = hi = None
    lo_s = hi_s = unb_s = None
    for s, a, b in zip(states, pa, pb):
        if a == 0:
            if b != 0 and unb_s is None:
                unb_s = s
            continue
        r = b / a
        if lo is None or r < lo:
            lo, lo_s = r, s
        if hi is None or r > hi:
            hi, hi_s = r, s
    claim = cert.get("claim")
    _need(isinstance(claim, dict) and set(claim) <= {"lower", "upper", "tight"}, "bad claim")
    tight = claim.get("tight", False)
    _need(isinstance(tight, bool), "bad tight flag")
    if lo is None:
        # Z_A = 0 for every exterior: the ratio is undefined. No lower or upper claim may be made
        # (adversarial review m3: previously any upper claim was accepted here).
        _need(claim.get("lower") is None and claim.get("upper") is None,
              "claims a bound but the support is empty")
        return True, ["support empty (Z_A = 0 for every exterior); no ratio claim made"]
    cl = _num(claim.get("lower"))
    _need(cl <= lo, "claimed lower bound exceeds min rho")
    up = claim.get("upper")
    if unb_s is not None:
        _need(up == "inf", "an unbounded boundary state exists but upper is not inf")
        _need(not tight or cl == lo, "tight claim: lower differs from min rho")
    elif up == "inf":
        _need(not tight, "tight claim with infinite upper bound but max rho is finite")
    else:
        cu = _num(up)
        _need(cu >= hi, "claimed upper bound is below max rho")
        if tight:
            _need(cl == lo and cu == hi, "tight claim differs from exact min/max rho")
    if tight and unb_s is None and up != "inf":
        pass
    w = cert.get("witness")
    if isinstance(w, dict) and tight:
        for key, val in (("lower", lo), ("upper", hi)):
            ws = w.get(key)
            if ws is None:
                continue
            _need(isinstance(ws, list) and len(ws) == len(B) and all(
                type(x) is int and 0 <= x < da[v] for x, v in zip(ws, B)), "bad witness state")
            k = 0
            for x, v in zip(ws, B):
                k = k * da[v] + x
            _need(pa[k] > 0 and pb[k] / pa[k] == val, f"{key} witness does not attain the claim")
    return True, ["ok"]
