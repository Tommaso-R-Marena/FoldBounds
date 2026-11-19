#!/usr/bin/env python3
"""Classify every hit of the forbidden-token ripgrep scan.

Runs exactly

    rg -n "native_decide|sorry|admit|implemented_by|unsafe|@\\[extern|axiom " FoldBounds scripts/*.lean

inside ROOT and, for every matched occurrence, decides with a small
comment-aware Lean lexer (nested ``/- -/`` blocks, ``--`` line comments,
string literals) whether the occurrence sits in

* ``comment``            a line/block/doc comment,
* ``string``             a string literal,
* ``print-axioms``       the argument/keyword of a ``#print axioms`` directive,
* ``CODE``               real Lean code (NOT benign: needs human review).

Usage: classify_rg_hits.py ROOT [--json OUT]
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

PATTERN = r'native_decide|sorry|admit|implemented_by|unsafe|@\[extern|axiom '
RG = f'rg -n "{PATTERN}" FoldBounds scripts/*.lean'


def region_map(src: str) -> list[str]:
    """Per character: 'code', 'comment' or 'string'."""
    n = len(src)
    reg = ["code"] * n
    i, depth, in_string = 0, 0, False
    while i < n:
        two = src[i:i + 2]
        if depth:
            if two == "/-":
                depth += 1
                reg[i] = reg[i + 1] = "comment"
                i += 2
            elif two == "-/":
                depth -= 1
                reg[i] = reg[i + 1] = "comment"
                i += 2
            else:
                reg[i] = "comment"
                i += 1
        elif in_string:
            reg[i] = "string"
            if src[i] == "\\" and i + 1 < n:
                reg[i + 1] = "string"
                i += 2
            else:
                in_string = src[i] != '"'
                i += 1
        elif two == "/-":
            depth = 1
            reg[i] = reg[i + 1] = "comment"
            i += 2
        elif two == "--":
            while i < n and src[i] != "\n":
                reg[i] = "comment"
                i += 1
        elif src[i] == "'" and src[i:i + 3] == "'\"'":
            i += 3                               # the character literal '"'
        elif src[i] == '"':
            in_string = True
            reg[i] = "string"
            i += 1
        else:
            i += 1
    return reg


def classify(root: Path) -> dict:
    proc = subprocess.run(["bash", "-c", RG], cwd=root, capture_output=True, text=True)
    hits = []
    cache: dict[str, tuple[str, list[str], list[int]]] = {}
    for line in proc.stdout.splitlines():
        m = re.match(r"^(.*?):(\d+):(.*)$", line)
        if not m:
            hits.append({"raw": line, "class": "UNPARSED"})
            continue
        file, lineno, text = m.group(1), int(m.group(2)), m.group(3)
        if file not in cache:
            src = (root / file).read_text(encoding="utf-8")
            starts = [0]
            for k, ch in enumerate(src):
                if ch == "\n":
                    starts.append(k + 1)
            cache[file] = (src, region_map(src), starts)
        src, reg, starts = cache[file]
        base = starts[lineno - 1]
        end = starts[lineno] - 1 if lineno < len(starts) else len(src)
        physical = src[base:end]
        for mm in re.finditer(PATTERN, physical):
            regions = set(reg[base + mm.start():base + mm.end()])
            if regions == {"comment"}:
                cls = "comment"
            elif regions == {"string"}:
                cls = "string"
            elif physical.lstrip().startswith("#print axioms"):
                cls = "print-axioms"
            else:
                cls = "CODE"
            hits.append({"file": file, "line": lineno, "column": mm.start() + 1,
                         "token": mm.group(0), "class": cls, "benign": cls != "CODE",
                         "text": text.strip()[:240]})
    # Independent cross-check without ripgrep (no ignore rules, no binary detection):
    # every regular file below FoldBounds/ plus scripts/*.lean, scanned with Python `re`.
    scanned = sorted([p for p in (root / "FoldBounds").rglob("*") if p.is_file()]
                     + list((root / "scripts").glob("*.lean")))
    py_lines = set()
    for p in scanned:
        for k, l in enumerate(p.read_text(encoding="utf-8", errors="replace").split("\n"), 1):
            if re.search(PATTERN, l):
                py_lines.add((str(p.relative_to(root)), k))
    rg_lines = {(h["file"], h["line"]) for h in hits if "file" in h}
    files = sorted({h["file"] for h in hits if "file" in h})
    counts: dict[str, int] = {}
    for h in hits:
        counts[h["class"]] = counts.get(h["class"], 0) + 1
    tokens: dict[str, int] = {}
    for h in hits:
        if "token" in h:
            tokens[h["token"]] = tokens.get(h["token"], 0) + 1
    return {"root": str(root), "command": RG, "rg_exit_code": proc.returncode,
            "rg_output_lines": len(proc.stdout.splitlines()), "occurrences": len(hits),
            "by_class": counts, "by_token": tokens, "files_with_hits": files,
            "python_rescan_files": len(scanned),
            "python_rescan_lean_files": sum(1 for p in scanned if p.suffix == ".lean"),
            "python_rescan_matching_lines": len(py_lines),
            "python_rescan_agrees_with_rg": py_lines == rg_lines,
            "python_rescan_only": sorted(py_lines - rg_lines), "rg_only": sorted(rg_lines - py_lines),
            "non_benign": [h for h in hits if not h.get("benign", False)], "hits": hits}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("--json")
    args = ap.parse_args()
    res = classify(Path(args.root))
    if args.json:
        Path(args.json).write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"root {res['root']}: rg exit {res['rg_exit_code']}, {res['rg_output_lines']} matching lines, "
          f"{res['occurrences']} occurrences, by class {res['by_class']}, by token {res['by_token']}")
    print(f"python rescan of {res['python_rescan_files']} files ({res['python_rescan_lean_files']} .lean): "
          f"{res['python_rescan_matching_lines']} matching lines, agrees with rg: "
          f"{res['python_rescan_agrees_with_rg']}")
    for h in res["hits"]:
        if "file" in h:
            print(f"  {h['file']}:{h['line']}:{h['column']} [{h['class']}] {h['token']!r}: {h['text']}")
        else:
            print(f"  UNPARSED {h['raw']}")
    sys.exit(0 if not res["non_benign"] else 1)


if __name__ == "__main__":
    main()
