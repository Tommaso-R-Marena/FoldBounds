#!/usr/bin/env python3
"""Lightweight REMC-style HP folder for FoldBounds witness search.

Untrusted heuristic: finds a conformation with ``k >= target`` contacts so the
Lean kernel can certify a lower bound. Inspired by the REMC protocol of
Thachuk–Shmygelska–Hoos 2007 (not a bit-exact reimplementation).
"""

from __future__ import annotations

import argparse
import json
import random
import time
from pathlib import Path

DIRS = ((1, 0), (-1, 0), (0, 1), (0, -1))


def contacts(seq: str, walk: list[tuple[int, int]]) -> int:
    occ = {p: i for i, p in enumerate(walk)}
    k = 0
    for i, (x, y) in enumerate(walk):
        if seq[i] != "H":
            continue
        for dx, dy in DIRS:
            j = occ.get((x + dx, y + dy))
            if j is not None and j < i - 1 and seq[j] == "H":
                k += 1
    return k // 2  # each contact counted twice


def grow_greedy(seq: str, rng: random.Random, temp: float = 0.0):
    """Grow an SAW with contact-biased random choices."""
    n = len(seq)
    occ = {(0, 0): 0, (1, 0): 1}
    walk = [(0, 0), (1, 0)]
    score = 0
    for m in range(2, n):
        x, y = walk[-1]
        cands = []
        for dx, dy in DIRS:
            p = (x + dx, y + dy)
            if p in occ:
                continue
            g = 0
            if seq[m] == "H":
                for ddx, ddy in DIRS:
                    j = occ.get((p[0] + ddx, p[1] + ddy))
                    if j is not None and j <= m - 2 and seq[j] == "H":
                        g += 1
            cands.append((g, p))
        if not cands:
            return None, -1
        cands.sort(key=lambda z: -z[0])
        if temp <= 1e-9:
            best = cands[0][0]
            pool = [p for g, p in cands if g == best]
            p = rng.choice(pool)
            g = best
        else:
            # softmax over contact gain
            weights = [pow(2.718281828, g / max(temp, 1e-6)) for g, _ in cands]
            total = sum(weights)
            r = rng.random() * total
            acc = 0.0
            p = cands[-1][1]
            g = cands[-1][0]
            for (gg, pp), w in zip(cands, weights):
                acc += w
                if r <= acc:
                    p, g = pp, gg
                    break
        occ[p] = m
        walk.append(p)
        score += g
    return walk, score


def pull_move(seq: str, walk: list[tuple[int, int]], rng: random.Random):
    """One end-pull / corner-flip style local move; returns new walk or None."""
    n = len(walk)
    i = rng.randrange(1, n - 1)
    # try relocate site i to a free neighbor of i-1 that is also adjacent to i+1 (corner) or
    # a free neighbor of an end.
    # End pull from head:
    if rng.random() < 0.5:
        # move last bead to a free neighbor of previous
        x, y = walk[-2]
        occ = set(walk[:-1])
        opts = [(x + dx, y + dy) for dx, dy in DIRS if (x + dx, y + dy) not in occ]
        if not opts:
            return None
        nw = walk[:-1] + [rng.choice(opts)]
        return nw
    # corner flip: if walk[i-1], walk[i], walk[i+1] form a corner, flip
    ax, ay = walk[i - 1]
    bx, by = walk[i]
    cx, cy = walk[i + 1]
    # proposed new position for i is the other corner of the unit square
    px, py = ax + (cx - bx), ay + (cy - by)
    # equivalently: (ax+cx-bx, ay+cy-by)
    px, py = ax + cx - bx, ay + cy - by
    if (px, py) == (bx, by):
        return None
    occ = set(walk)
    occ.discard((bx, by))
    if (px, py) in occ:
        return None
    # must stay adjacent to both neighbors
    if abs(px - ax) + abs(py - ay) != 1:
        return None
    if abs(px - cx) + abs(py - cy) != 1:
        return None
    nw = list(walk)
    nw[i] = (px, py)
    return nw


def remc_search(seq: str, target: int, seconds: float = 30.0, seed: int = 0):
    rng = random.Random(seed)
    temps = [0.05, 0.2, 0.5, 1.0, 2.0]
    replicas = []
    for t in temps:
        w, s = None, -1
        for _ in range(50):
            w, s = grow_greedy(seq, rng, temp=t)
            if w is not None:
                break
        if w is None:
            w, s = grow_greedy(seq, rng, temp=0.0)
        replicas.append([w, s if w else -1, t])

    best = max(replicas, key=lambda r: r[1])
    if best[1] >= target:
        return best[0], best[1], 0.0

    t0 = time.perf_counter()
    swaps = 0
    steps = 0
    while time.perf_counter() - t0 < seconds:
        for r in replicas:
            if r[0] is None:
                continue
            for _ in range(20):
                steps += 1
                cand = pull_move(seq, r[0], rng)
                if cand is None:
                    continue
                sc = contacts(seq, cand)
                # Metropolis at replica temperature (energy = -contacts)
                dE = -(sc - r[1])
                if dE <= 0 or rng.random() < pow(2.718281828, -dE / max(r[2], 1e-6)):
                    r[0], r[1] = cand, sc
                    if sc > best[1]:
                        best = [list(cand), sc, r[2]]
                    if sc >= target:
                        return cand, sc, time.perf_counter() - t0
        # adjacent replica exchange
        for i in range(len(replicas) - 1):
            a, b = replicas[i], replicas[i + 1]
            if a[0] is None or b[0] is None:
                continue
            # exchange probability on energy difference
            beta_a, beta_b = 1.0 / max(a[2], 1e-6), 1.0 / max(b[2], 1e-6)
            d = (beta_a - beta_b) * (-(a[1] - b[1]))
            if d <= 0 or rng.random() < pow(2.718281828, -d):
                a[0], b[0] = b[0], a[0]
                a[1], b[1] = b[1], a[1]
                swaps += 1
        # occasional fresh greedy seed into hottest replica
        if steps % 5000 == 0:
            w, s = grow_greedy(seq, rng, temp=temps[-1])
            if w is not None and s > replicas[-1][1]:
                replicas[-1][0], replicas[-1][1] = w, s
                if s >= target:
                    return w, s, time.perf_counter() - t0
    return best[0], best[1], time.perf_counter() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--id", action="append", default=[])
    ap.add_argument("--seconds", type=float, default=60.0)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    bench = json.loads((Path(__file__).parent / "benchmarks.json").read_text())
    ids = args.id or ["S1-6", "S1-7"]
    for iid in ids:
        inst = next(x for x in bench["ground_state_2d"] if x["id"] == iid)
        walk, k, wall = remc_search(inst["sequence"], inst["contacts"], seconds=args.seconds, seed=args.seed)
        print(
            json.dumps(
                {
                    "id": iid,
                    "target": inst["contacts"],
                    "found": k,
                    "hit": k >= inst["contacts"],
                    "wall_s": wall,
                    "walk": walk,
                }
            )
        )


if __name__ == "__main__":
    main()
