"""Differential and adversarial tests for the ``FOLDBOUNDS-KERNEL 1`` checker.

Run from the repository root::

    python -m pytest experiments/kernel -q

Three independent pieces are compared:

* ``independent_oracle``  - exhaustive enumeration, ground truth of the model;
* ``reference_kernel``    - pure-Python implementation of the specification;
* ``foldkernel``          - the compiled Lean checker (tests are skipped when
  ``.lake/build/bin/foldkernel`` is missing; override with ``FOLDKERNEL_BIN``).

Sections: (a) lemma attacks on the bound entries, (b) reference vs oracle,
(c) binary vs reference, (d) adversarial mutations and permanent fixtures.

Any disagreement between the binary and the reference, any crash, and any
accepted result that fails to enclose the oracle is saved as a reproducer
under ``experiments/kernel/results/findings/`` and fails the test.
"""

from __future__ import annotations

import functools
import hashlib
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import independent_oracle as oracle  # noqa: E402
import mutations  # noqa: E402
import producer  # noqa: E402
import reference_kernel as ref  # noqa: E402

BINARY = Path(os.environ.get("FOLDKERNEL_BIN", str(ROOT / ".lake" / "build" / "bin" / "foldkernel")))
FIXTURES = ROOT / "tests" / "kernel_fixtures"
RESULTS = HERE / "results"
FINDINGS = RESULTS / "findings"
SCRATCH_ROOT = Path(os.environ.get("KERNEL_SCRATCH", "/tmp/kernel_ref"))

Q_SET = [Fraction(1), Fraction(3, 2), Fraction(2), Fraction(3), Fraction(7, 3)]
T_SET = [0, 1, 2, 3, 5, 9]
EXTRA_Q = [Fraction(10, 9), Fraction(22, 7), Fraction(1000001, 1000000), Fraction(97),
           Fraction(2 ** 64 + 1, 2 ** 64), Fraction(5, 4)]
EXTRA_T = [4, 6, 7, 8, 12, 100]
ORACLE_MAX_N = 16
BATCH = 150
COUNTERS = ("n", "cert_bytes", "split_nodes", "bound_leaves", "enum_leaves",
            "exact_leaves", "exact_count")

SUMMARY: dict = {"findings": []}


# ---------------------------------------------------------------------------
# Session plumbing
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def binary_available() -> bool:
    return BINARY.is_file() and os.access(BINARY, os.X_OK)


@pytest.fixture(scope="session", autouse=True)
def _session_summary():
    started = time.time()
    if binary_available():
        SUMMARY["binary"] = {"path": str(BINARY), "sha256_at_start": _sha256(BINARY)}
    yield
    oracle.flush_cache()
    if binary_available() and "binary" in SUMMARY:
        SUMMARY["binary"]["sha256_at_end"] = _sha256(BINARY)
    SUMMARY["wall_seconds"] = round(time.time() - started, 1)
    try:
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / "differential_summary.json").write_text(
            json.dumps(SUMMARY, indent=1, sort_keys=True) + "\n", encoding="ascii")
    except OSError:
        pass


@pytest.fixture(scope="module")
def scratch():
    SCRATCH_ROOT.mkdir(parents=True, exist_ok=True)
    path = Path(tempfile.mkdtemp(prefix="differential-", dir=str(SCRATCH_ROOT)))
    yield path
    if not os.environ.get("KERNEL_KEEP_SCRATCH"):
        shutil.rmtree(path, ignore_errors=True)


def _jsonable(obj):
    if isinstance(obj, Fraction):
        return ref.ratio_str(obj)
    if isinstance(obj, bytes):
        return obj.decode("latin-1")
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    return repr(obj)


def record_finding(kind: str, name: str, data, details: dict) -> str:
    """Save a reproducer (certificate + description); returns its path."""
    FINDINGS.mkdir(parents=True, exist_ok=True)
    payload = data if isinstance(data, bytes) else b""
    digest = hashlib.sha256(payload + name.encode()).hexdigest()[:10]
    stem = f"{kind}-{mutations._slug(name)[:60]}-{digest}"
    if isinstance(data, bytes):
        (FINDINGS / f"{stem}.cert").write_bytes(data)
    (FINDINGS / f"{stem}.json").write_text(
        ref.to_json(_jsonable({"kind": kind, "name": name, **details})) + "\n", encoding="ascii")
    rel = str((FINDINGS / stem).relative_to(ROOT))
    SUMMARY["findings"].append(rel)
    return rel


def _fail_if(problems: list[str], title: str) -> None:
    if problems:
        shown = "\n  ".join(problems[:8])
        more = f"\n  ... and {len(problems) - 8} more" if len(problems) > 8 else ""
        pytest.fail(f"{title}: {len(problems)} problem(s)\n  {shown}{more}", pytrace=False)


# ---------------------------------------------------------------------------
# Runner and comparison helpers
# ---------------------------------------------------------------------------

def _big_int(text: str) -> int:
    negative = text.startswith("-")
    value = ref.nat_from_digits((text[1:] if negative else text).encode("ascii"))
    return -value if negative else value


def parse_ratio(text: str) -> Fraction:
    num, den = text.split("/")
    return Fraction(_big_int(num), _big_int(den))


def strip_zeros(hist) -> list[int]:
    hist = list(hist)
    while hist and hist[-1] == 0:
        hist.pop()
    return hist


def run_binary(paths, timeout: float = 240.0) -> dict[str, dict]:
    """Run ``foldkernel check`` on ``paths``; one report per path.

    A report is the parsed JSON object of the runner.  When the runner produced
    no parsable object for a file (crash, abort, hang) the report is
    ``{"crash": True, ...}`` and the remaining files are re-run.
    """
    reports: dict[str, dict] = {}
    queue = [str(p) for p in paths]
    while queue:
        chunk, queue = queue[:BATCH], queue[BATCH:]
        timed_out = False
        for attempt in range(4):
            try:
                proc = subprocess.run([str(BINARY), "check", *chunk],
                                      capture_output=True, timeout=timeout)
                stdout, stderr, rc = proc.stdout, proc.stderr, proc.returncode
                break
            except subprocess.TimeoutExpired as exc:
                stdout, stderr, rc, timed_out = exc.stdout or b"", exc.stderr or b"", None, True
                break
            except OSError:
                if attempt == 3:
                    pytest.skip("foldkernel binary cannot be executed (being rebuilt?)")
                time.sleep(2.0)
        done = 0
        for name, line in zip(chunk, [ln for ln in stdout.split(b"\n") if ln.strip()]):
            try:
                report = json.loads(line.decode("utf-8"), parse_int=_big_int)
            except (ValueError, UnicodeDecodeError):
                break
            if not isinstance(report, dict) or report.get("file") != name:
                break
            reports[name] = report
            done += 1
        if done < len(chunk):
            reports[chunk[done]] = {
                "crash": True, "returncode": rc, "timeout": timed_out,
                "stderr": stderr[-1500:].decode("utf-8", "replace"),
                "stdout_tail": stdout[-300:].decode("utf-8", "replace"),
            }
            queue = chunk[done + 1:] + queue
        else:
            want = 0 if all(reports[name].get("accepted") is True for name in chunk) else 1
            if rc != want:
                reports[chunk[-1]] = dict(reports[chunk[-1]], exit_status_problem=(
                    f"exit status {rc}, expected {want} for this batch"))
    return reports


def diff_reports(r: dict, b: dict) -> list[str]:
    """Differences between a reference report and a runner report."""
    if b.get("crash"):
        return [f"binary produced no report (returncode={b.get('returncode')}, "
                f"timeout={b.get('timeout')}, stderr={b.get('stderr', '')[-200:]!r})"]
    problems = []
    if "exit_status_problem" in b:
        problems.append(b["exit_status_problem"])
    if b.get("accepted") not in (True, False):
        return problems + [f"binary report without boolean 'accepted': {b}"]
    if r["accepted"] != b["accepted"]:
        return problems + [f"verdict differs: reference accepted={r['accepted']} "
                           f"(reason={r.get('reason')}), binary accepted={b['accepted']} "
                           f"(reason={b.get('reason')})"]
    if not r["accepted"]:
        if not isinstance(b.get("reason"), str) or not b["reason"]:
            problems.append("binary rejection without a reason category")
        return problems
    for key in ("format", "sequence") + COUNTERS:
        if b.get(key) != r[key]:
            problems.append(f"{key}: reference {r[key]!r}, binary {b.get(key)!r}")
    if strip_zeros(b.get("exact_histogram", ["missing"])) != strip_zeros(r["exact_histogram"]):
        problems.append(f"exact_histogram: reference {strip_zeros(r['exact_histogram'])}, "
                        f"binary {b.get('exact_histogram')}")
    rows_r, rows_b = r["results"], b.get("results")
    if not isinstance(rows_b, list) or len(rows_b) != len(rows_r):
        return problems + [f"results: reference has {len(rows_r)} rows, binary {rows_b!r:.200}"]
    for i, (x, y) in enumerate(zip(rows_r, rows_b)):
        try:
            if parse_ratio(x["q"]) != parse_ratio(y["q"]):
                problems.append(f"results[{i}].q: reference {x['q']}, binary {y['q']}")
            if x["threshold"] != y["threshold"]:
                problems.append(f"results[{i}].threshold differs")
            for name in ref.ENDPOINTS:
                if parse_ratio(x[name]) != parse_ratio(y[name]):
                    problems.append(f"results[{i}].{name}: reference {x[name]:.80}, "
                                    f"binary {y[name]:.80}")
                elif x[name] != y[name]:
                    problems.append(f"results[{i}].{name}: equal value but the binary prints "
                                    f"{y[name]:.80} (not lowest terms num/den?)")
        except (KeyError, ValueError, TypeError, ZeroDivisionError) as exc:
            problems.append(f"results[{i}]: malformed binary row ({exc!r})")
    return problems


def enclosure_violations(values: dict, exact: dict, must_be_exact: bool = False) -> list[str]:
    """Violations of ``lo <= exact <= hi`` for the eight endpoints (Fractions)."""
    bad = []
    pairs = (("z", exact["Z"]), ("e", exact["ZA"]), ("n", exact["ZN"]), ("p", exact["P"]))
    for stem, truth in pairs:
        lo, hi = values[stem + "Lo"], values[stem + "Hi"]
        if not lo <= truth <= hi:
            bad.append(f"{stem}Lo={lo} <= {truth} <= {stem}Hi={hi} fails")
        if must_be_exact and not lo == truth == hi:
            bad.append(f"{stem}: no bound leaves, but [{lo}, {hi}] != exact {truth}")
    if not 0 <= values["pLo"] <= values["pHi"] <= 1:
        bad.append(f"0 <= pLo={values['pLo']} <= pHi={values['pHi']} <= 1 fails")
    return bad


def oracle_violations(report: dict) -> list[str]:
    """Check an accepted runner-format report against the oracle (n <= 16)."""
    seq = report["sequence"]
    if len(seq) > ORACLE_MAX_N:
        return []
    bad = []
    exact_only = report["bound_leaves"] == 0
    if exact_only and strip_zeros(report["exact_histogram"]) != oracle.cached_density_of_states(seq) \
            and oracle.cached_density_of_states(seq) != [0]:
        bad.append(f"no bound leaves but histogram {report['exact_histogram']} != oracle "
                   f"{oracle.cached_density_of_states(seq)}")
    for i, row in enumerate(report["results"]):
        q = parse_ratio(row["q"])
        exact = oracle.exact_values(seq, q, row["threshold"])
        values = {name: parse_ratio(row[name]) for name in ref.ENDPOINTS}
        for msg in enclosure_violations(values, exact, must_be_exact=exact_only):
            bad.append(f"results[{i}] (q={row['q']}, t={str(row['threshold'])[:30]}): {msg}")
    return bad


# ---------------------------------------------------------------------------
# Certificate corpus (deterministic)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Case:
    name: str
    mode: str
    seq: str
    tree: str
    queries: tuple

    def data(self) -> bytes:
        return producer.to_text(self.seq, self.queries, self.tree)


def _queries(rng: random.Random) -> tuple:
    out = []
    for _ in range(rng.randint(1, 3)):
        q = rng.choice(Q_SET) if rng.random() < 0.8 else rng.choice(EXTRA_Q)
        t = rng.choice(T_SET) if rng.random() < 0.8 else rng.choice(EXTRA_T)
        out.append((q, t))
    return tuple(out)


def _sequence_pool(rng: random.Random, n: int) -> list[str]:
    if n == 2:
        return ["HH", "HP", "PH", "PP"]
    pool: list[str] = []
    for s in ("H" * n, ("HP" * n)[:n], ("PH" * n)[:n]):
        if s not in pool:
            pool.append(s)
    if n % 3 == 0:
        pool.append("P" * n)
    while len(pool) < 6:
        s = "".join("H" if rng.random() < 0.6 else "P" for _ in range(n))
        if s not in pool:
            pool.append(s)
    return pool


@functools.lru_cache(maxsize=None)
def corpus() -> tuple:
    rng = random.Random(0x5EEDF01D)
    cases: list[Case] = []

    def add(mode: str, seq: str, tree: str) -> None:
        cases.append(Case(f"{len(cases):04d}-{mode}-n{len(seq)}", mode, seq, tree, _queries(rng)))

    for n in range(2, 15):
        pool = _sequence_pool(rng, n)
        for seq in pool:
            add("root", seq, "X" if rng.random() < 0.5 else "B")
            for _ in range(2):
                depth = rng.randint(2, min(n + 1, 9))
                add("uniform", seq, producer.uniform(seq, depth, rng.choice([None, None, 0, 1, 2, 3, 5])))
            add("hybrid", seq, producer.hybrid(seq, rng.randint(2, min(n, 9)), rng.choice([0, 1, 1, 2, 3])))
            for _ in range(2):
                q = rng.choice(Q_SET + EXTRA_Q[:2])
                t = rng.choice([1, 1, 2, 2, 3, 5, 0, 9])
                budget = rng.choice([10, 60, 300, 1500, 8000, 40000])
                add("adaptive", seq, producer.adaptive(seq, q, t, budget))
            add("witness", seq, producer.witness_paths(
                seq, Fraction(2), rng.choice([1, 2, 3, 5]),
                max_witnesses=rng.randint(1, 4), node_limit=3000))
            for _ in range(2):
                add("random", seq, producer.random_tree(
                    seq, rng, p_split=rng.choice([0.5, 0.7, 0.9]),
                    p_enum=rng.choice([0.0, 0.25, 0.6]), max_bytes=1500))
            if n <= 12:
                prefix = producer.random_prefix(n, rng.randint(2, n), rng)
                add("isolate", seq, producer.isolate_prefix(seq, prefix))
        if n <= 10:
            add("full", pool[-1], producer.uniform(pool[-1], n))
    return tuple(cases)


@functools.lru_cache(maxsize=None)
def mutants() -> tuple:
    return tuple(mutations.generate(limit=10))


# ---------------------------------------------------------------------------
# (a) Oracle self-consistency and lemma attacks on the bound entries
# ---------------------------------------------------------------------------

def test_oracle_totals_match_saw_counts():
    for n in range(2, 15):
        for seq in ("H" * n, "P" * n, ("HP" * n)[:n]):
            g = oracle.density_of_states(seq)
            assert sum(g) == oracle.SAW_COUNTS[n - 1] // 4, seq
        assert oracle.density_of_states("P" * n) == [oracle.SAW_COUNTS[n - 1] // 4]
    assert sum(oracle.density_of_states("H" * 16)) == 6416596 // 4
    assert oracle.density_of_states("HHHH") == [7, 2]


def test_oracle_fast_enumerator_matches_naive_scorer_all_sequences_up_to_8():
    info = oracle.self_check(max_n=8, prefix_max_n=6)
    assert info["sequences"] == sum(2 ** n for n in range(2, 9))
    SUMMARY["oracle_self_check"] = info


def _lemma_violations(seq: str, cross_check_a: bool):
    """Compare the reference's bound entry of every prefix with the oracle."""
    n = len(seq)
    bad: list[dict] = []
    stats = {"prefixes": 0, "nonempty": 0, "k_tight": 0, "fam_positive": 0}

    def visit(sites, count, kmin, kmax):
        a, cap, cnt, fam = ref.bound_entry(seq, sites)
        r = n - len(sites)
        stats["prefixes"] += 1
        why = []
        if cnt != 3 ** r:
            why.append("cnt != 3^r")
        if count > cnt:
            why.append("number of completions > cnt")
        if fam > count:
            why.append("fam > number of completions")
        if fam < 0:
            why.append("fam < 0")
        if count:
            stats["nonempty"] += 1
            stats["k_tight"] += 1 if kmax == cap else 0
            if a > kmin:
                why.append("a > min k over completions")
            if kmax > cap:
                why.append("max k over completions > K")
        if fam:
            stats["fam_positive"] += 1
        if cross_check_a and oracle.naive_contacts(seq, sites) != a:
            why.append("a != contacts among the placed residues")
        if why and len(bad) < 25:
            bad.append({"sequence": seq, "prefix": list(sites), "entry": [a, cap, cnt, fam],
                        "completions": count, "kmin": kmin, "kmax": kmax, "violations": why})

    root = oracle.prefix_stats(seq, visit)
    assert root[0] == oracle.SAW_COUNTS[n - 1] // 4
    return bad, stats


def _report_lemma(bad: list[dict], label: str) -> None:
    if bad:
        path = record_finding("lemma-violation", label, None, {"counterexamples": bad})
        pytest.fail(f"bound-entry lemma violated ({label}); first counterexample: {bad[0]}; "
                    f"all saved in {path}.json", pytrace=False)


def test_lemma_attack_all_sequences_up_to_8():
    bad: list[dict] = []
    total = {"sequences": 0, "prefixes": 0, "nonempty": 0, "k_tight": 0, "fam_positive": 0}
    for n in range(2, 9):
        for letters in product("HP", repeat=n):
            b, stats = _lemma_violations("".join(letters), cross_check_a=True)
            bad.extend(b)
            total["sequences"] += 1
            for key, value in stats.items():
                total[key] += value
    SUMMARY["lemma_exhaustive_n_le_8"] = total
    _report_lemma(bad, "all sequences n<=8")
    assert total["sequences"] == 508


def test_lemma_attack_random_sample_9_to_13():
    rng = random.Random(20261118)
    plan = {9: 16, 10: 12, 11: 8, 12: 5, 13: 3}
    bad: list[dict] = []
    total = {"sequences": 0, "prefixes": 0, "nonempty": 0, "k_tight": 0, "fam_positive": 0}
    for n, how_many in plan.items():
        seqs = ["H" * n, ("HP" * n)[:n]]
        while len(seqs) < how_many:
            s = "".join("H" if rng.random() < 0.65 else "P" for _ in range(n))
            if s not in seqs:
                seqs.append(s)
        for seq in seqs:
            b, stats = _lemma_violations(seq, cross_check_a=False)
            bad.extend(b)
            total["sequences"] += 1
            for key, value in stats.items():
                total[key] += value
    SUMMARY["lemma_sample_n_9_to_13"] = total
    _report_lemma(bad, "random sample n=9..13")


def test_incremental_bound_entries_equal_literal_definition():
    """``check_tree`` computes entries incrementally; ``bound_entry`` literally."""
    rng = random.Random(7)
    leaves = 0
    seqs = ["".join(p) for n in range(2, 7) for p in product("HP", repeat=n)]
    seqs += ["".join(rng.choice("HHP") for _ in range(n)) for n in (7, 8, 9, 10, 11) for _ in range(4)]
    for seq in seqs:
        for depth in range(2, len(seq) + 1):
            acc = ref.check_tree(seq, producer.uniform(seq, depth), paranoid=True)
            leaves += acc.bound_leaves
            # depth == n: every leaf is a complete walk, written as E
            assert (acc.bound_leaves == 0) == (depth == len(seq)) or len(seq) == 2
    for _ in range(300):
        n = rng.randint(4, 40)
        seq = "".join(rng.choice("HP") for _ in range(n))
        assert ref._future_cap_table_fast(seq) == list(ref._future_cap_table_literal(seq)), seq
    SUMMARY["incremental_vs_literal_bound_leaves"] = leaves


# ---------------------------------------------------------------------------
# (b) Reference vs oracle
# ---------------------------------------------------------------------------

def test_reference_accepts_specification_example():
    data = (b"FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery 2 1\n"
            b"tree NNE-EE-NEEE-NEE-E\nend\n")
    seq, queries, tree_start = ref.parse(data)
    assert (seq, queries) == ("HHHH", [(Fraction(2), 1)])
    assert data[tree_start:-5] == b"NNE-EE-NEEE-NEE-E"
    out = ref.verify(data)
    assert out == {
        "accepted": True, "format": "FOLDBOUNDS-KERNEL 1", "sequence": "HHHH", "n": 4,
        "cert_bytes": 17, "split_nodes": 4, "bound_leaves": 0, "enum_leaves": 0,
        "exact_leaves": 9, "exact_count": 9, "exact_histogram": [7, 2],
        "results": [{"q": "2/1", "threshold": 1, "zLo": "11/1", "zHi": "11/1",
                     "eLo": "4/1", "eHi": "4/1", "nLo": "7/1", "nHi": "7/1",
                     "pLo": "4/11", "pHi": "4/11"}],
    }
    acc = ref.check_tree("HHHH", b"B")
    assert acc.entries == [(0, 1, 9, 7)] and len(acc.hist) == 2 * 4 + 2
    assert ref.verify(data[:-1])["accepted"] is False


def test_reference_all_exact_certificates_reproduce_oracle():
    rng = random.Random(11)
    seqs = ["".join(p) for n in range(2, 6) for p in product("HP", repeat=n)]
    seqs += ["H" * n for n in range(6, 14)]
    seqs += ["".join(rng.choice("HHP") for _ in range(n)) for n in range(6, 14) for _ in range(3)]
    checked = 0
    for seq in seqs:
        n = len(seq)
        trees = {"X", producer.hybrid(seq, min(n, 5), 0), producer.uniform(seq, min(n, 6), n)}
        if n <= 10:
            trees.add(producer.uniform(seq, n))
        g = oracle.cached_density_of_states(seq)
        for tree in sorted(trees):
            acc = ref.check_tree(seq, tree, paranoid=True)
            assert acc.entries == [] and acc.bound_leaves == 0
            assert len(acc.hist) == 2 * n + 2
            assert strip_zeros(acc.hist) == g, (seq, tree[:40])
            if set(tree) <= set("NE-"):
                assert acc.exact_leaves == sum(g) and acc.enum_leaves == 0
            for q in Q_SET:
                for t in T_SET:
                    bad = enclosure_violations(ref.evaluate(acc, q, t),
                                               oracle.exact_values(seq, q, t), must_be_exact=True)
                    assert not bad, (seq, tree[:40], q, t, bad)
                    checked += 1
    SUMMARY["all_exact_evaluations"] = checked


def test_reference_encloses_oracle_on_corpus():
    problems: list[str] = []
    evaluations = 0
    modes: dict[str, int] = {}
    for case in corpus():
        modes[case.mode] = modes.get(case.mode, 0) + 1
        acc = ref.check_tree(case.seq, case.tree, paranoid=True)
        exact_only = not acc.entries
        for q in Q_SET:
            for t in T_SET:
                values = ref.evaluate(acc, q, t)
                bad = enclosure_violations(values, oracle.exact_values(case.seq, q, t),
                                           must_be_exact=exact_only)
                evaluations += 1
                if bad:
                    path = record_finding(
                        "reference-vs-oracle", f"{case.name}-q{q}-t{t}", case.data(),
                        {"sequence": case.seq, "q": q, "t": t, "violations": bad,
                         "values": values})
                    problems.append(f"{case.name} q={q} t={t}: {bad[0]} [{path}]")
        if len(acc.entries) <= 600:
            q, t = case.queries[0]
            assert ref.evaluate_per_entry(acc, q, t) == ref.evaluate(acc, q, t), case.name
    SUMMARY["corpus"] = {"certificates": len(corpus()), "by_mode": modes,
                         "reference_vs_oracle_evaluations": evaluations}
    _fail_if(problems, "reference kernel fails to enclose the oracle")
    assert len(corpus()) >= 300


def test_command_line_tools(scratch):
    out = scratch / "cli.cert"
    env = dict(os.environ, FOLDBOUNDS_ORACLE_CACHE=str(scratch / "cli_cache.json"))
    for mode, extra in (("uniform", ["--depth", "6", "--enum-r", "1"]),
                        ("hybrid", ["--depth", "6", "--kappa", "1"]),
                        ("adaptive", ["--budget", "2000"]),
                        ("witness", [])):
        proc = subprocess.run(
            [sys.executable, str(HERE / "producer.py"), "--seq", "HPHPPHHPHH", "--query", "2", "4",
             "--query", "3/2", "2", "--mode", mode, *extra, "--out", str(out)],
            capture_output=True, timeout=120)
        assert proc.returncode == 0, proc.stderr
        check = subprocess.run([sys.executable, str(HERE / "reference_kernel.py"), "check", str(out)],
                               capture_output=True, timeout=120)
        assert check.returncode == 0, check.stdout
        report = json.loads(check.stdout)
        assert report["accepted"] is True and report == {"file": str(out), **ref.verify(out.read_bytes())}
        assert not oracle_violations(report)
    proc = subprocess.run([sys.executable, str(HERE / "independent_oracle.py"), "HPHPPHHPHH",
                           "--q", "3/2", "--t", "2"], capture_output=True, timeout=120, env=env)
    assert proc.returncode == 0, proc.stderr
    got = json.loads(proc.stdout)
    assert got["density_of_states"] == oracle.density_of_states("HPHPPHHPHH")
    assert got["Z"] == "41729/8" and got["P"] == "5769/41729"
    bad = subprocess.run([sys.executable, str(HERE / "reference_kernel.py"), "check",
                          str(FIXTURES / "does-not-exist.txt")], capture_output=True, timeout=60)
    assert bad.returncode == 2


# ---------------------------------------------------------------------------
# (c) Binary vs reference
# ---------------------------------------------------------------------------

def _compare_files(named_data: list[tuple[str, bytes]], directory: Path, kind: str):
    """Write files, run both checkers, return (rows, problems).

    ``rows`` are ``(name, data, reference_report, binary_report)``.
    """
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for i, (_, data) in enumerate(named_data):
        path = directory / f"{i:05d}.cert"
        path.write_bytes(data)
        paths.append(path)
    reports = run_binary(paths)
    rows = []
    problems: list[str] = []
    for (name, data), path in zip(named_data, paths):
        r = ref.verify(data)
        b = reports[str(path)]
        rows.append((name, data, r, b))
        diffs = diff_reports(r, b)
        if not diffs and r["accepted"]:
            diffs = [f"binary result does not enclose the oracle: {msg}"
                     for msg in oracle_violations(b)]
        if diffs:
            saved = record_finding(kind, name, data, {
                "differences": diffs, "reference": r,
                "binary": {k: v for k, v in b.items() if k not in ("verify_ms", "verify_us")}})
            problems.append(f"{name}: {diffs[0]} [{saved}.cert]")
    return rows, problems


def test_binary_agrees_with_reference_on_corpus(scratch):
    if not binary_available():
        pytest.skip(f"foldkernel binary not found at {BINARY}")
    cases = corpus()
    rows, problems = _compare_files([(c.name, c.data()) for c in cases],
                                    scratch / "corpus", "binary-vs-reference")
    accepted = sum(1 for _, _, r, _ in rows if r["accepted"])
    SUMMARY["differential"] = {
        "certificates": len(rows),
        "accepted_by_both": sum(1 for _, _, r, b in rows if r["accepted"] and b.get("accepted")),
        "queries_compared": sum(len(r.get("results", [])) for _, _, r, _ in rows),
        "endpoints_compared": 8 * sum(len(r.get("results", [])) for _, _, r, _ in rows),
        "disagreements": len(problems),
    }
    _fail_if(problems, "foldkernel disagrees with the reference kernel on generated certificates")
    assert accepted == len(rows), "every producer output must be a valid certificate"
    assert len(rows) >= 300


def test_binary_agrees_with_reference_on_single_bound_entries(scratch):
    """One bound leaf per certificate, everything else enumerated: any error
    in a single recomputed entry (a, K, cnt, fam) shows up in the sums."""
    if not binary_available():
        pytest.skip(f"foldkernel binary not found at {BINARY}")
    rng = random.Random(424242)
    named = []
    queries = ((Fraction(97), 0), (Fraction(2), 1), (Fraction(3, 2), 2), (Fraction(7, 3), 3))
    for n in range(2, 11):
        for k in range(24 if n > 3 else 4):
            seq = "".join("H" if rng.random() < 0.7 else "P" for _ in range(n))
            prefix = producer.random_prefix(n, rng.randint(2, n), rng)
            tree = producer.isolate_prefix(seq, prefix)
            named.append((f"isolate-n{n}-{k}-{seq}-m{len(prefix)}",
                          producer.to_text(seq, queries, tree)))
    rows, problems = _compare_files(named, scratch / "isolate", "binary-vs-reference")
    for name, _, r, _ in rows:
        assert r["accepted"] and r["bound_leaves"] == 1, name
    SUMMARY["differential_single_entry"] = {"certificates": len(rows),
                                            "disagreements": len(problems)}
    _fail_if(problems, "foldkernel disagrees with the reference kernel on single bound entries")


# ---------------------------------------------------------------------------
# (d) Adversarial mutations
# ---------------------------------------------------------------------------

def test_mutation_catalogue_reference_expectations():
    problems: list[str] = []
    table: dict[str, dict[str, int]] = {}
    for m in mutants():
        r = ref.verify(m.data)
        verdict = "accept" if r["accepted"] else "reject"
        row = table.setdefault(m.category, {"expect": m.expect, "accept": 0, "reject": 0})
        row[verdict] += 1
        label = f"{m.category}/{m.base}/{m.label}"
        if m.expect != "either" and verdict != m.expect:
            saved = record_finding("mutation-expectation", label, m.data,
                                   {"expected": m.expect, "reference": r})
            problems.append(f"{label}: expected {m.expect}, reference says {verdict} [{saved}.cert]")
        elif r["accepted"]:
            bad = oracle_violations(r)
            if bad:
                saved = record_finding("reference-vs-oracle", label, m.data,
                                       {"violations": bad, "reference": r})
                problems.append(f"{label}: {bad[0]} [{saved}.cert]")
    SUMMARY["mutations"] = {"categories": len(table), "mutants": len(mutants()),
                            "by_category": table}
    _fail_if(problems, "mutation catalogue vs reference kernel")
    assert set(table) == {c.name for c in mutations.CATEGORIES}
    assert len(table) >= 40


def test_binary_agrees_with_reference_on_mutations(scratch):
    if not binary_available():
        pytest.skip(f"foldkernel binary not found at {BINARY}")
    ms = mutants()
    rows, problems = _compare_files(
        [(f"{m.category}/{m.base}/{m.label}", m.data) for m in ms],
        scratch / "mutants", "binary-vs-reference-mutation")
    both_reject = sum(1 for _, _, r, b in rows if not r["accepted"] and b.get("accepted") is False)
    both_accept = sum(1 for _, _, r, b in rows if r["accepted"] and b.get("accepted") is True)
    reasons: dict[str, int] = {}
    for _, _, r, b in rows:
        if b.get("accepted") is False:
            reasons[b.get("reason", "?")] = reasons.get(b.get("reason", "?"), 0) + 1
    SUMMARY["mutations_binary"] = {
        "mutants": len(rows), "both_reject": both_reject, "both_accept": both_accept,
        "disagreements_or_crashes": len(problems), "binary_reject_reasons": reasons}
    _fail_if(problems, "foldkernel disagrees with the reference kernel on mutated certificates")
    assert both_reject + both_accept == len(rows)


def test_binary_handles_large_inputs_cleanly(scratch):
    """Megabyte-sized trees, maximal nesting, many queries: no crash, same verdict."""
    if not binary_available():
        pytest.skip(f"foldkernel binary not found at {BINARY}")
    raw = mutations._raw
    n = 400
    chain = mutations.deep_chain(n)
    seq400 = ("HHPHP" * 80)[:n]
    many_queries = tuple((producer.format_q(q), str(t))
                         for q in Q_SET + EXTRA_Q for t in range(0, 40))
    named = [
        ("large/1MB-of-N-short-sequence", raw("HHHH", "N" * 1_000_000), "reject"),
        ("large/1MB-of-E", raw("HPHH", "E" * 1_000_000), "reject"),
        ("large/1MB-of-B", raw("HPHH", "B" * 1_000_000), "reject"),
        ("large/1MB-of-X", raw("HPHH", "X" * 1_000_000), "reject"),
        ("large/1MB-of-blocked-markers", raw("HPHH", "-" * 1_000_000), "reject"),
        ("large/1MB-of-dead-ends", raw("HPHH", "N----" * 200_000), "reject"),
        ("large/chain400-valid", raw(seq400, chain, (("2", "3"), ("3/2", "0"), ("1", "150"))), "accept"),
        ("large/chain400-all-H-valid", raw("H" * n, chain, (("7/3", "100"),)), "accept"),
        ("large/N398-then-tree-ends", raw(seq400, "N" * (n - 2)), "reject"),
        ("large/N398-then-1MB-of-N", raw(seq400, "N" * (n - 2) + "N" * 1_000_000), "reject"),
        ("large/chain400-then-1MB-trailing", raw(seq400, chain + "B" * 1_000_000), "reject"),
        ("large/chain400-one-branch-omitted", raw(seq400, chain[:-2] + "-B"), "reject"),
        ("large/440-queries", raw("HPHHPHHPHH", producer.uniform("HPHHPHHPHH", 6, 1), many_queries),
         "accept"),
        ("large/20000-digit-threshold",
         raw("HPHHPHHPHH", producer.uniform("HPHHPHHPHH", 5), (("2", "7" * 20000),)), "accept"),
        ("large/1MB-sequence-of-garbage", b"FOLDBOUNDS-KERNEL 1\nsequence " + b"h" * 1_000_000
         + b"\nquery 2 1\ntree X\nend\n", "reject"),
        ("large/1MB-query-token", b"FOLDBOUNDS-KERNEL 1\nsequence HHHH\nquery " + b"1" * 500_000
         + b"/" + b"0" * 500_000 + b" 1\ntree X\nend\n", "reject"),
        ("large/1MB-of-NUL", b"\x00" * 1_000_000, "reject"),
        ("large/1MB-of-newlines", b"\n" * 1_000_000, "reject"),
    ]
    rows, problems = _compare_files([(name, data) for name, data, _ in named],
                                    scratch / "large", "binary-vs-reference-large")
    for (name, _, want), (_, _, r, _) in zip(named, rows):
        got = "accept" if r["accepted"] else "reject"
        if got != want:
            problems.append(f"{name}: reference says {got}, expected {want}")
    SUMMARY["large_inputs"] = {"files": len(rows), "problems": len(problems)}
    _fail_if(problems, "large inputs")


# ---------------------------------------------------------------------------
# Permanent fixtures (tests/kernel_fixtures)
# ---------------------------------------------------------------------------

def _manifest() -> dict:
    path = FIXTURES / "MANIFEST.json"
    if not path.is_file():
        pytest.fail(f"{path} is missing; run: python experiments/kernel/mutations.py --write-fixtures",
                    pytrace=False)
    return json.loads(path.read_text(encoding="ascii"))


def test_fixtures_are_what_the_generator_writes():
    plan = mutations.fixture_plan()
    assert _manifest() == mutations.manifest_of(plan), \
        "stale MANIFEST.json; run: python experiments/kernel/mutations.py --write-fixtures"
    for item in plan:
        assert (FIXTURES / item["file"]).read_bytes() == item["data"], item["file"]
    on_disk = sorted(f"{sub}/{p.name}" for sub in ("accept", "reject")
                     for p in (FIXTURES / sub).iterdir())
    assert on_disk == sorted(item["file"] for item in plan)


def test_fixtures_reference_verdicts():
    manifest = _manifest()
    problems = []
    categories = set()
    for entry in manifest["files"]:
        assert entry["expect"] in ("accept", "reject")
        assert entry["file"].startswith(entry["expect"] + "/")
        assert entry.get("category"), entry
        categories.add(entry["category"])
        data = (FIXTURES / entry["file"]).read_bytes()
        r = ref.verify(data)
        if r["accepted"] != (entry["expect"] == "accept"):
            problems.append(f"{entry['file']}: manifest expects {entry['expect']}, "
                            f"reference accepted={r['accepted']} ({r.get('reason')})")
        elif r["accepted"]:
            problems.extend(f"{entry['file']}: {msg}" for msg in oracle_violations(r))
    SUMMARY["fixtures"] = {
        "files": len(manifest["files"]),
        "accept": sum(1 for e in manifest["files"] if e["expect"] == "accept"),
        "reject": sum(1 for e in manifest["files"] if e["expect"] == "reject"),
        "categories": len(categories - {"base"}),
    }
    _fail_if(problems, "fixtures vs reference kernel")
    assert categories - {"base"} == {c.name for c in mutations.CATEGORIES}


def test_fixtures_binary_verdicts():
    if not binary_available():
        pytest.skip(f"foldkernel binary not found at {BINARY}")
    manifest = _manifest()
    paths = [FIXTURES / entry["file"] for entry in manifest["files"]]
    reports = run_binary(paths)
    problems = []
    for entry, path in zip(manifest["files"], paths):
        data = path.read_bytes()
        r = ref.verify(data)
        b = reports[str(path)]
        diffs = diff_reports(r, b)
        if not diffs and b.get("accepted") != (entry["expect"] == "accept"):
            diffs = [f"manifest expects {entry['expect']}, binary accepted={b.get('accepted')}"]
        if not diffs and b.get("accepted"):
            diffs = oracle_violations(b)
        if diffs:
            saved = record_finding("binary-vs-fixture", entry["file"], data,
                                   {"differences": diffs, "reference": r, "binary": b})
            problems.append(f"{entry['file']}: {diffs[0]} [{saved}.cert]")
    _fail_if(problems, "foldkernel vs permanent fixtures")
