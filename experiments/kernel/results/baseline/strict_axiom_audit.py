#!/usr/bin/env python3
"""Strict, independent audit of a `#print axioms` log.

The repository gate (``scripts/check_lean_sources.py``) inspects only the first
physical line of each report and extracts axiom names with an ASCII
identifier regex.  This script instead

* parses every report ``'Name' depends on axioms: [...]`` INCLUDING continuation
  lines up to the closing ``]`` (and ``'Name' does not depend on any axioms``);
* splits the bracket content on commas and compares each name by EXACT string
  equality with ``propext``, ``Classical.choice``, ``Quot.sound``;
* requires the log to consist of reports only (any residue is reported);
* compares the ordered list of reported names with EVERY ``#print axioms``
  directive of the audit source (not only ASCII ``FoldBounds.*`` names) after
  comments have been stripped from that source.

Usage: strict_axiom_audit.py LOG AUDIT_SOURCE [--json OUT]
Exit status 0 iff the strict audit is clean.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

STANDARD = ("propext", "Classical.choice", "Quot.sound")
SUSPICIOUS = ("sorryAx", "Lean.ofReduceBool", "Lean.trustCompiler", "Lean.reduceBool",
              "Lean.reduceNat", "Lean.ofReduceNat")
FLAGSHIP = ("FoldBounds.FrontierCheck.check_sound",
            "FoldBounds.FrontierCheck.checkFlat_sound",
            "FoldBounds.FrontierNext.check_sound",
            "FoldBounds.FrontierNextFlat.checkFlat2_sound")

HEAD = re.compile(r"^'(?P<name>.+?)' (?P<kind>depends on axioms: \[|does not depend on any axioms)",
                  re.MULTILINE)


def strip_lean_comments(src: str) -> str:
    """Blank out Lean comments (nested block comments, line comments), keep strings."""
    out = list(src)
    i, n, depth = 0, len(src), 0
    in_string = False
    while i < n:
        two = src[i:i + 2]
        if depth:
            if two == "/-":
                depth += 1
                out[i] = out[i + 1] = " "
                i += 2
            elif two == "-/":
                depth -= 1
                out[i] = out[i + 1] = " "
                i += 2
            else:
                if src[i] != "\n":
                    out[i] = " "
                i += 1
        elif in_string:
            if src[i] == "\\":
                i += 2
            else:
                in_string = src[i] != '"'
                i += 1
        elif two == "/-":
            depth = 1
            out[i] = out[i + 1] = " "
            i += 2
        elif two == "--":
            while i < n and src[i] != "\n":
                out[i] = " "
                i += 1
        elif src[i] == '"':
            in_string = True
            i += 1
        else:
            i += 1
    return "".join(out)


def parse_log(text: str) -> tuple[list[dict], str]:
    reports, residue, pos = [], [], 0
    for m in HEAD.finditer(text):
        if m.start() < pos:
            continue                         # inside a previous report body
        residue.append(text[pos:m.start()])
        if m.group("kind").startswith("does not"):
            end = m.end()
            raw = text[m.start():end]
            reports.append({"name": m.group("name"), "axioms": [], "raw": raw,
                            "physical_lines": 1, "closed": True})
        else:
            close = text.find("]", m.end())
            closed = close != -1
            end = close + 1 if closed else len(text)
            raw = text[m.start():end]
            body = text[m.end():close] if closed else text[m.end():]
            axioms = [a.strip() for a in body.split(",")]
            reports.append({"name": m.group("name"), "axioms": axioms, "raw": raw,
                            "physical_lines": raw.count("\n") + 1, "closed": closed})
        pos = end
    residue.append(text[pos:])
    return reports, "".join(residue)


def audit(log_path: Path, source_path: Path) -> dict:
    text = log_path.read_text(encoding="utf-8")
    reports, residue = parse_log(text)
    source = source_path.read_text(encoding="utf-8")
    code = strip_lean_comments(source)
    directives = [l.strip()[len("#print axioms"):].strip() for l in code.splitlines()
                  if l.strip().startswith("#print axioms")]
    raw_directive_lines = sum(1 for l in source.splitlines() if l.strip().startswith("#print axioms"))
    names = [r["name"] for r in reports]
    by_name = {}
    for r in reports:
        by_name.setdefault(r["name"], []).append(r)
    distinct = sorted({a for r in reports for a in r["axioms"]})
    nonstandard = sorted(a for a in distinct if a not in STANDARD)
    offenders = [{"name": r["name"], "axioms": r["axioms"]} for r in reports
                 if any(a not in STANDARD for a in r["axioms"])]
    guillemets = [r["name"] for r in reports if "«" in r["raw"] or "»" in r["raw"]]
    non_ascii = [r["name"] for r in reports if any(ord(c) > 127 for c in r["raw"])]
    histogram: dict[str, int] = {}
    for r in reports:
        key = "[" + ", ".join(r["axioms"]) + "]" if r["axioms"] else "(no axioms)"
        histogram[key] = histogram.get(key, 0) + 1
    flagship = {}
    for want in FLAGSHIP:
        if want in by_name:
            flagship[want] = {"found": True, "report": by_name[want][0]["raw"]}
        else:
            leaf = want.rsplit(".", 1)[1]
            near = [r["raw"] for r in reports if r["name"].rsplit(".", 1)[-1] == leaf]
            flagship[want] = {"found": False, "report": None,
                              "reports_with_same_final_component": near}
    result = {
        "log": str(log_path), "log_bytes": len(text.encode()), "log_physical_lines": len(text.splitlines()),
        "audit_source": str(source_path),
        "print_axioms_directives_in_source_code": len(directives),
        "print_axioms_lines_in_source_including_comments": raw_directive_lines,
        "directive_names_not_matching_repo_regex": [
            d for d in directives if not re.fullmatch(r"FoldBounds\.[A-Za-z0-9_.]+", d)],
        "audited_declarations": len(reports),
        "distinct_declarations": len(set(names)),
        "duplicate_reports": sorted({n for n in names if names.count(n) > 1}),
        "reports_spanning_multiple_physical_lines": [
            {"name": r["name"], "physical_lines": r["physical_lines"]}
            for r in reports if r["physical_lines"] > 1],
        "unclosed_reports": [r["name"] for r in reports if not r["closed"]],
        "reports_without_axioms": [r["name"] for r in reports if not r["axioms"]],
        "distinct_axiom_names_exact": distinct,
        "axiom_set_histogram": histogram,
        "nonstandard_axiom_names": nonstandard,
        "suspicious_names_present": [s for s in SUSPICIOUS if s in distinct or s in text],
        "declarations_with_nonstandard_axioms": offenders,
        "reports_containing_guillemets": guillemets,
        "reports_containing_non_ascii": non_ascii,
        "unparsed_residue": residue.strip(),
        "reported_names_equal_source_directives_in_order": names == directives,
        "missing_from_log": [d for d in directives if d not in by_name],
        "in_log_but_not_in_source": [n for n in names if n not in directives],
        "flagship_reports": flagship,
    }
    result["strict_audit_clean"] = bool(
        reports and not nonstandard and not offenders and not guillemets and not non_ascii
        and not result["unparsed_residue"] and not result["duplicate_reports"]
        and not result["unclosed_reports"] and result["reported_names_equal_source_directives_in_order"]
        and not result["suspicious_names_present"])
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("audit_source")
    ap.add_argument("--json")
    args = ap.parse_args()
    res = audit(Path(args.log), Path(args.audit_source))
    if args.json:
        Path(args.json).write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print(f"audited declarations: {res['audited_declarations']} "
          f"(distinct {res['distinct_declarations']}; source directives "
          f"{res['print_axioms_directives_in_source_code']})")
    print(f"distinct axiom names (exact match): {res['distinct_axiom_names_exact']}")
    print(f"axiom-set histogram: {res['axiom_set_histogram']}")
    print(f"non-standard axiom names: {res['nonstandard_axiom_names']}")
    print(f"suspicious names present anywhere in log: {res['suspicious_names_present']}")
    print(f"multi-line reports: {len(res['reports_spanning_multiple_physical_lines'])}; "
          f"guillemets: {len(res['reports_containing_guillemets'])}; "
          f"non-ASCII: {len(res['reports_containing_non_ascii'])}; "
          f"residue bytes: {len(res['unparsed_residue'])}")
    print(f"names equal source directives in order: "
          f"{res['reported_names_equal_source_directives_in_order']}; missing: "
          f"{res['missing_from_log']}; extra: {res['in_log_but_not_in_source']}")
    for name, info in res["flagship_reports"].items():
        if info["found"]:
            print(f"FLAGSHIP {info['report']}")
        else:
            print(f"FLAGSHIP '{name}': NO REPORT UNDER THIS NAME; same final component: "
                  f"{info['reports_with_same_final_component']}")
    print("STRICT AUDIT CLEAN" if res["strict_audit_clean"] else "STRICT AUDIT NOT CLEAN")
    sys.exit(0 if res["strict_audit_clean"] else 1)


if __name__ == "__main__":
    main()
