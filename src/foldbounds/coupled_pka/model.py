"""Rational multi-site protonation model (see experiments/coupled_pka/PROTOCOL.md)."""
from __future__ import annotations

import json
import math
from decimal import Decimal, ROUND_HALF_EVEN, localcontext
from fractions import Fraction
from functools import lru_cache

MODEL_PKA = {"C-term": "3.20", "ASP": "3.80", "GLU": "4.50", "HIS": "6.50", "CYS": "9.00",
             "TYR": "10.00", "LYS": "10.50", "ARG": "12.50", "N-term": "8.00"}
ACID = {"ASP", "GLU", "TYR", "CYS", "C-term"}
CENTRE_ATOMS = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2"), "HIS": ("ND1", "NE2"),
                "CYS": ("SG",), "TYR": ("OH",), "LYS": ("NZ",), "ARG": ("CZ",)}
COULOMB = Fraction("332.0637")  # kcal A / mol / e^2
KB = Fraction("0.001987204")
TEMP = Fraction("298.15")
RC = 12.0


@lru_cache(maxsize=None)
def pow10_rational(e_milli: int) -> Fraction:
    """R15(10^(e/1000)): 60-digit decimal power rounded to 15 significant digits."""
    with localcontext() as ctx:
        ctx.prec = 60
        x = Decimal(10) ** (Decimal(e_milli) / Decimal(1000))
        ctx.prec = 15
        ctx.rounding = ROUND_HALF_EVEN
        x = +x
    return Fraction(x)


def kt_ln10() -> float:
    return float(KB * TEMP) * math.log(10)


def coupling_milli(r: float) -> int:
    """Quantised pair coupling in thousandths of a pK unit (eps(r)=4r)."""
    w = float(COULOMB) / 4.0 / (r * r) / kt_ln10()
    return int(round(1000.0 * w))  # Python round = half-even


# ---------------------------------------------------------------- structure
def parse_pdb(text: str):
    """Return dict chain -> list of residues (resname, resnum, icode, {atom: (x,y,z)}); model 1."""
    chains: dict = {}
    order: dict = {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM"):
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        ch = line[21]
        key = (line[17:20].strip(), int(line[22:26]), line[26].strip())
        atom = line[12:16].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        d = chains.setdefault(ch, {})
        if key not in d:
            d[key] = {}
            order.setdefault(ch, []).append(key)
        d[key].setdefault(atom, xyz)
    return {ch: [(k[0], k[1], k[2], chains[ch][k]) for k in order[ch]] for ch in chains}


def _mean(pts):
    n = len(pts)
    return tuple(sum(p[i] for p in pts) / n for i in range(3))


def build_sites(pdb_text: str, chains: list[str]):
    """Return (sites, dropped). A site = dict(chain,resnum,icode,rtype,xyz,pka)."""
    res = parse_pdb(pdb_text)
    sg = []
    for ch in res:
        for (rn, num, ic, at) in res[ch]:
            if rn == "CYS" and "SG" in at:
                sg.append(((ch, num, ic), at["SG"]))
    ss = set()
    for a in range(len(sg)):
        for b in range(a + 1, len(sg)):
            if math.dist(sg[a][1], sg[b][1]) < 2.5:
                ss.add(sg[a][0]); ss.add(sg[b][0])
    sites, dropped = [], []
    for ch in chains:
        if ch not in res:
            continue
        rl = res[ch]
        for idx, (rn, num, ic, at) in enumerate(rl):
            cands = []
            if rn in CENTRE_ATOMS:
                if rn == "CYS" and (ch, num, ic) in ss:
                    continue
                cands.append((rn, CENTRE_ATOMS[rn]))
            if idx == 0:
                cands.append(("N-term", ("N",)))
            if idx == len(rl) - 1:
                cands.append(("C-term", ("O", "OXT") if "OXT" in at else ("O",)))
            for rtype, names in cands:
                pts = [at[n] for n in names if n in at]
                if len(pts) != len(names):
                    dropped.append((ch, num, ic, rtype, "missing atom"))
                    continue
                sites.append({"chain": ch, "resnum": num, "icode": ic, "rtype": rtype,
                              "xyz": [round(v, 4) for v in _mean(pts)], "pka": MODEL_PKA[rtype]})
    return sites, dropped


def build_spec(sites, name: str = "", rc: float = RC) -> dict:
    """Integer/decimal-string spec: the single trusted input of the checker."""
    n = len(sites)
    out_sites = []
    for s in sites:
        out_sites.append({"id": f"{s['chain']}:{s['rtype']}{s['resnum']}{s['icode']}",
                          "res": [s["chain"], s["resnum"], s["icode"]], "rtype": s["rtype"],
                          "acid": 1 if s["rtype"] in ACID else 0, "pka": s["pka"], "xyz": s["xyz"]})
    pairs = {}
    for i in range(n):
        for j in range(i + 1, n):
            if sites[i]["chain"] == sites[j]["chain"] and sites[i]["resnum"] == sites[j]["resnum"] \
                    and sites[i]["icode"] == sites[j]["icode"]:
                continue
            r = math.dist(sites[i]["xyz"], sites[j]["xyz"])
            if r <= rc:
                k = coupling_milli(r)
                if k > 0:
                    pairs[f"{i},{j}"] = k
    return {"name": name, "T": "298.15", "rc": rc, "sites": out_sites, "pairs_milli": pairs}


def spec_sub(spec: dict, keep: list[int], name: str = "") -> dict:
    """Sub-system on a subset of sites (re-indexed)."""
    pos = {s: k for k, s in enumerate(keep)}
    pairs = {}
    for key, v in spec["pairs_milli"].items():
        i, j = map(int, key.split(","))
        if i in pos and j in pos:
            a, b = sorted((pos[i], pos[j]))
            pairs[f"{a},{b}"] = v
    return {"name": name or spec["name"], "T": spec["T"], "rc": spec["rc"],
            "sites": [spec["sites"][s] for s in keep], "pairs_milli": pairs}


def spec_hash(spec: dict) -> str:
    import hashlib
    return hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()


# ---------------------------------------------------------------- rational MRF
class Model:
    """Rational MRF from a spec at a given pH (Fraction, multiple of 1/1000)."""

    def __init__(self, spec: dict, ph: Fraction):
        self.spec = spec
        self.n = len(spec["sites"])
        self.ph = Fraction(ph)
        assert (self.ph * 1000).denominator == 1
        self.acid = [s["acid"] for s in spec["sites"]]
        self.u1 = []
        for s in spec["sites"]:
            e = (Fraction(s["pka"]) - self.ph) * 1000
            assert e.denominator == 1
            self.u1.append(pow10_rational(int(e)))
        self.nbr = [dict() for _ in range(self.n)]  # j -> k_milli
        for key, v in spec["pairs_milli"].items():
            i, j = map(int, key.split(","))
            self.nbr[i][j] = v
            self.nbr[j][i] = v

    def phi(self, i: int, j: int, si: int, sj: int) -> Fraction:
        k = self.nbr[i].get(j)
        if not k:
            return Fraction(1)
        qi, qj = si - self.acid[i], sj - self.acid[j]
        return pow10_rational(-k * qi * qj)

    def u(self, i: int, s: int) -> Fraction:
        return self.u1[i] if s else Fraction(1)
