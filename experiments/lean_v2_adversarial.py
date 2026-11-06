#!/usr/bin/env python3
"""Differential, adversarial and mutation tests of the version-2 Lean checker.

1. Differential sweep: many small instances (all-H, all-P, alternating, mixed;
   several q, thresholds, depths and leaf modes).  Every certificate must be
   accepted and every accepted interval must enclose the exact values computed
   by the independent Python enumerator ``exact_spectrum``.
2. Adversarial catalogue: hand-written malformed or dishonest certificates
   that must be rejected (parse error or ``accepted: false``).
3. Random mutation testing: valid certificates are mutated (dropped / duplicated
   lines, flipped direction characters, forged kinds, bogus witnesses).  A
   mutant may be rejected, or accepted -- but then it must still enclose the
   exact values (soundness).  Any accepted non-enclosing mutant is a failure.

Tests are evidence, not proof; soundness is the Lean theorem
``FoldBounds.FrontierNext.checkFlat2_sound``.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path as FsPath

HERE = FsPath(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))

from foldbounds.lattice import build_frontier, exact_spectrum  # noqa: E402
from foldbounds.lean_export_v2 import certificate_to_lean_text_v2, raw_lean_text_v2  # noqa: E402
from foldbounds.witnesses import propose_dfs_witnesses  # noqa: E402


def frac(s: str) -> Fraction:
    a, b = s.split("/")
    return Fraction(int(a), int(b))


def run_batch(exe: str, texts: list[str]) -> list[dict]:
    out = []
    with tempfile.TemporaryDirectory() as d:
        files = []
        for i, t in enumerate(texts):
            f = os.path.join(d, f"c{i}.txt")
            FsPath(f).write_text(t)
            files.append(f)
        for i in range(0, len(files), 200):
            res = subprocess.run([exe, "check", *files[i:i + 200]], capture_output=True, text=True,
                                 check=True, timeout=3600)
            out += [json.loads(line) for line in res.stdout.strip().splitlines()]
    return out


def exact(seq: str, q: Fraction, thr: int):
    spec = exact_spectrum(seq)
    z = sum((c * q ** k for k, c in spec.items()), Fraction(0))
    a = sum((c * q ** k for k, c in spec.items() if k >= thr), Fraction(0))
    return z, a, z - a, a / z


def encloses(r: dict, ex) -> bool:
    z, a, n, p = ex
    g = {k: frac(r[k]) for k in ("zLo", "zHi", "eLo", "eHi", "nLo", "nHi", "pLo", "pHi")}
    return (g["zLo"] <= z <= g["zHi"] and g["eLo"] <= a <= g["eHi"] and g["nLo"] <= n <= g["nHi"]
            and g["pLo"] <= p <= g["pHi"] and 0 <= g["zLo"] and 0 <= g["eLo"] and 0 <= g["nLo"])


SEQS = ["HHHH", "HHHHH", "PPPPPP", "HPHPHPH", "HHPPHHPH", "HPHHPHHPH", "HHHHHHHHHH",
        "PHHPPHHPHH", "HPPHPPHHPP"]
QS = [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(5)]
THRS = [0, 1, 2, 3, 5]


def adversarial_catalogue() -> list[tuple[str, str]]:
    s, q = "HHHHHH", Fraction(2)
    good_nodes = [l for l in certificate_to_lean_text_v2(build_frontier(s, q, 2, 3), leaf="M",
                                                          memory=2).splitlines()[4:-1]]
    cat = [
        ("bad header", "FOLDBOUNDS-FRONTIER 3\nsequence HHHHHH\nq 2\nthreshold 2\nX .\nend\n"),
        ("v1 header with X leaf", "FOLDBOUNDS-FRONTIER 1\nsequence HHHHHH\nq 2\nthreshold 2\nX .\nend\n"),
        ("v1 header with M leaf", "FOLDBOUNDS-FRONTIER 1\nsequence HHHHHH\nq 2\nthreshold 2\nM2 .\nend\n"),
        ("q = 0", raw_lean_text_v2(s, Fraction(0), 2, ["X ."])),
        ("q = 1/2", raw_lean_text_v2(s, Fraction(1, 2), 2, ["X ."])),
        ("negative q", "FOLDBOUNDS-FRONTIER 2\nsequence HHHHHH\nq -2\nthreshold 2\nX .\nend\n"),
        ("q denominator 0", "FOLDBOUNDS-FRONTIER 2\nsequence HHHHHH\nq 2/0\nthreshold 2\nX .\nend\n"),
        ("bad sequence char", "FOLDBOUNDS-FRONTIER 2\nsequence HHXHHH\nq 2\nthreshold 2\nX .\nend\n"),
        ("sequence too short", "FOLDBOUNDS-FRONTIER 2\nsequence H\nq 2\nthreshold 0\nX .\nend\n"),
        ("negative threshold", "FOLDBOUNDS-FRONTIER 2\nsequence HHHHHH\nq 2\nthreshold -1\nX .\nend\n"),
        ("missing end", "FOLDBOUNDS-FRONTIER 2\nsequence HHHHHH\nq 2\nthreshold 2\nX .\n"),
        ("unsupported kind", raw_lean_text_v2(s, q, 2, ["Z ."])),
        ("memory without number", raw_lean_text_v2(s, q, 2, ["M ."])),
        ("memory not a number", raw_lean_text_v2(s, q, 2, ["Mx ."])),
        ("negative memory", raw_lean_text_v2(s, q, 2, ["M-1 ."])),
        ("bad direction char", raw_lean_text_v2(s, q, 2, ["M2 0x"])),
        ("exact at nonterminal root", raw_lean_text_v2(s, q, 2, ["E ."])),
        ("dead root", raw_lean_text_v2(s, q, 2, ["D ."])),
        ("duplicate witness", raw_lean_text_v2(s, q, 2, ["M2 . 2113 2113"])),
        ("colliding witness", raw_lean_text_v2(s, q, 2, ["M2 . 2133"])),
        ("short witness", raw_lean_text_v2(s, q, 2, ["M2 . 211"])),
        ("long witness", raw_lean_text_v2(s, q, 2, ["M2 . 21130"])),
        ("witness on exact leaf", raw_lean_text_v2(s, q, 2, ["E 0000 0"])),
        ("enum below chain length", raw_lean_text_v2(s, q, 2, ["X 00000"])),
        ("missing branch", raw_lean_text_v2(s, q, 2, good_nodes[1:])),
        ("duplicated entry", raw_lean_text_v2(s, q, 2, good_nodes + good_nodes[:1])),
        ("ancestor and descendant", raw_lean_text_v2(s, q, 2, ["M2 ."] + good_nodes)),
        ("mem leaf on colliding prefix", raw_lean_text_v2(s, q, 2, good_nodes + ["M0 1"])),
        ("trailing junk on line", raw_lean_text_v2(s, q, 2, ["X . junk"])),
    ]
    return cat


def mutate(text: str, rng: random.Random) -> str:
    lines = text.splitlines()
    head, body, tail = lines[:4], lines[4:-1], lines[-1:]
    op = rng.choice(["drop", "dup", "flip", "kind", "witness", "swapdir"])
    if not body:
        return text
    i = rng.randrange(len(body))
    if op == "drop":
        body = body[:i] + body[i + 1:]
    elif op == "dup":
        body = body + [body[i]]
    elif op == "flip":
        toks = body[i].split(" ")
        w = toks[1]
        if w != ".":
            j = rng.randrange(len(w))
            w = w[:j] + rng.choice("0123") + w[j + 1:]
            toks[1] = w
        body[i] = " ".join(toks)
    elif op == "kind":
        toks = body[i].split(" ")
        toks[0] = rng.choice(["E", "D", "X", "M0", "M3", "B", "T"])
        body[i] = " ".join(toks[:2]) if toks[0] in ("E", "D", "X") else " ".join(toks)
    elif op == "witness":
        toks = body[i].split(" ")
        if toks[0] in ("E", "D", "X"):
            toks[0] = "M1"
        toks.append("".join(rng.choice("0123") for _ in range(rng.randrange(1, 9))))
        body[i] = " ".join(toks)
    else:
        toks = body[i].split(" ")
        if len(toks[1]) > 1:
            w = list(toks[1])
            j = rng.randrange(len(w) - 1)
            w[j], w[j + 1] = w[j + 1], w[j]
            toks[1] = "".join(w)
        body[i] = " ".join(toks)
    return "\n".join(head + body + tail) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe2", default=".lake/build/bin/foldcheck2")
    ap.add_argument("--out", default=str(HERE / "results" / "v2"))
    ap.add_argument("--mutants", type=int, default=600)
    args = ap.parse_args()
    out = FsPath(args.out)
    out.mkdir(parents=True, exist_ok=True)
    report: dict = {}

    # 1. differential sweep
    texts, meta = [], []
    for seq in SEQS:
        n = len(seq)
        for q in QS:
            for thr in THRS:
                for depth in sorted({2, min(n, 3), min(n, (n + 2) // 2)}):
                    cert = build_frontier(seq, q, thr, depth)
                    wit = {nd.path: propose_dfs_witnesses(nd.path, n, node_budget=20)
                           for nd in cert.frontier if nd.kind == "bound"}
                    for mode in ("M0", "M2", "M4", "M3W", "X", "T"):
                        if mode == "X":
                            t = certificate_to_lean_text_v2(cert, enum_remaining=n)
                        elif mode == "T":
                            t = certificate_to_lean_text_v2(cert, leaf="T", witnesses=wit)
                        else:
                            t = certificate_to_lean_text_v2(cert, leaf="M", memory=int(mode[1]),
                                                            witnesses=wit if mode.endswith("W") else None)
                        texts.append(t)
                        meta.append((seq, q, thr, depth, mode))
    results = run_batch(args.exe2, texts)
    cache: dict = {}
    diff_fail = []
    exact_hits = 0
    for (seq, q, thr, depth, mode), r in zip(meta, results):
        key = (seq, q, thr)
        if key not in cache:
            cache[key] = exact(seq, q, thr)
        ok = bool(r.get("accepted")) and encloses(r, cache[key])
        if mode == "X" and ok and frac(r["zLo"]) == frac(r["zHi"]) == cache[key][0]:
            exact_hits += 1
        if not ok:
            diff_fail.append({"case": [seq, str(q), thr, depth, mode], "result": r})
    report["differential"] = {"certificates": len(texts), "failures": diff_fail,
                              "enum_exact": exact_hits,
                              "enum_certificates": sum(1 for m in meta if m[4] == "X")}
    print(f"differential: {len(texts)} certificates, {len(diff_fail)} failures, "
          f"{exact_hits}/{report['differential']['enum_certificates']} X-certificates exact")

    # 2. adversarial catalogue
    cat = adversarial_catalogue()
    res = run_batch(args.exe2, [t for _, t in cat])
    adv = []
    for (name, _), r in zip(cat, res):
        rejected = ("parse_error" in r) or (r.get("accepted") is False)
        adv.append({"name": name, "rejected": rejected, "result": r})
    report["adversarial"] = adv
    n_rej = sum(1 for a in adv if a["rejected"])
    print(f"adversarial: {n_rej}/{len(adv)} rejected")
    for a in adv:
        if not a["rejected"]:
            print("  NOT REJECTED:", a["name"], a["result"])

    # 3. mutation testing
    rng = random.Random(20261019)
    mtexts, mmeta = [], []
    for _ in range(args.mutants):
        seq = rng.choice(SEQS[:8])
        q = rng.choice(QS)
        thr = rng.choice(THRS)
        depth = rng.choice(sorted({2, min(len(seq), 3), min(len(seq), 4)}))
        cert = build_frontier(seq, q, thr, depth)
        base = certificate_to_lean_text_v2(cert, leaf="M", memory=rng.choice([0, 2, 3]))
        mtexts.append(mutate(base, rng))
        mmeta.append((seq, q, thr))
    mres = run_batch(args.exe2, mtexts)
    accepted = rejected = unsound = 0
    bad = []
    for (seq, q, thr), r, t in zip(mmeta, mres, mtexts):
        if r.get("accepted"):
            accepted += 1
            key = (seq, q, thr)
            if key not in cache:
                cache[key] = exact(seq, q, thr)
            if not encloses(r, cache[key]):
                unsound += 1
                bad.append({"text": t, "result": r})
        else:
            rejected += 1
    report["mutation"] = {"mutants": len(mtexts), "accepted_sound": accepted - unsound,
                          "rejected": rejected, "accepted_unsound": unsound, "unsound_cases": bad}
    print(f"mutation: {len(mtexts)} mutants, {rejected} rejected, {accepted - unsound} accepted "
          f"and still enclosing, {unsound} accepted but NOT enclosing")
    (out / "lean_v2_adversarial.json").write_text(json.dumps(report, indent=1, default=str))
    failed = diff_fail or unsound or n_rej != len(adv)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
