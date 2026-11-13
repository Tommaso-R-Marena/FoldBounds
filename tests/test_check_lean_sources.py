"""Regression tests for scripts/check_lean_sources.py (audit findings F3e and F3f).

The gate is run as a subprocess in a temporary copy of `FoldBounds/`, exactly as CI runs it.
Synthetic logs follow Lean's `#print axioms` format. The two exploit reports are REAL Lean
4.28.0 output, captured from a throwaway build (see tests/fixtures/axiom_gate/README.md).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "check_lean_sources.py"
FIXTURE = ROOT / "tests" / "fixtures" / "axiom_gate" / "real_exploit_reports.log"


def audited_names() -> list[str]:
    text = (ROOT / "FoldBounds" / "AxiomAudit.lean").read_text(encoding="utf-8")
    return re.findall(r"^#print axioms (FoldBounds\.[A-Za-z0-9_.]+)\s*$", text, re.MULTILINE)


NAMES = audited_names()


@pytest.fixture
def workdir(tmp_path: Path) -> Path:
    shutil.copytree(ROOT / "FoldBounds", tmp_path / "FoldBounds")
    return tmp_path


def run_gate(workdir: Path, log_text: str | None = None) -> subprocess.CompletedProcess:
    args = [sys.executable, str(SCRIPT)]
    if log_text is not None:
        log = workdir / "axioms.log"
        log.write_text(log_text, encoding="utf-8")
        args.append(str(log))
    return subprocess.run(args, cwd=workdir, capture_output=True, text=True)


def report(name: str, axioms: str = "[propext, Classical.choice, Quot.sound]") -> list[str]:
    return [f"'{name}' depends on axioms: {axioms}"]


def wrapped(name: str) -> list[str]:
    return [f"'{name}' depends on axioms: [propext,", " Classical.choice,", " Quot.sound]"]


def build_log(victim_lines: list[str] | None = None, wrap: bool = False) -> str:
    """Valid reports for every audited name, with the first one optionally replaced."""
    lines: list[str] = []
    for index, name in enumerate(NAMES):
        if index == 0 and victim_lines is not None:
            lines += victim_lines
        else:
            lines += wrapped(name) if wrap else report(name)
    return "\n".join(lines) + "\n"


def captured_report(marker: str) -> list[str]:
    """One real captured report, renamed to the first audited declaration."""
    lines = FIXTURE.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if marker in line)
    end = start + 1
    while end < len(lines) and lines[end][:1].isspace():
        end += 1
    block = lines[start:end]
    block[0] = re.sub(r"^'[^']+'", f"'{NAMES[0]}'", block[0])
    return block


def assert_rejected(result: subprocess.CompletedProcess, needle: str) -> None:
    assert result.returncode != 0, result.stdout + result.stderr
    assert needle in result.stdout + result.stderr, result.stdout + result.stderr


# --- legitimate output keeps passing -------------------------------------------------------

def test_source_audit_alone_passes_on_repository(workdir: Path) -> None:
    result = run_gate(workdir)
    assert result.returncode == 0, result.stderr
    assert "Source audit:" in result.stdout


def test_single_line_reports_pass(workdir: Path) -> None:
    result = run_gate(workdir, build_log())
    assert result.returncode == 0, result.stderr
    assert f"Axiom audit: {len(NAMES)} declarations" in result.stdout


def test_wrapped_reports_with_standard_axioms_pass(workdir: Path) -> None:
    result = run_gate(workdir, build_log(wrap=True))
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("axioms", ["[propext]", "[propext, Quot.sound]", "[Quot.sound, propext]"])
def test_subsets_of_standard_axioms_pass(workdir: Path, axioms: str) -> None:
    result = run_gate(workdir, build_log(report(NAMES[0], axioms)))
    assert result.returncode == 0, result.stderr


def test_no_axiom_report_passes(workdir: Path) -> None:
    result = run_gate(workdir, build_log([f"'{NAMES[0]}' does not depend on any axioms"]))
    assert result.returncode == 0, result.stderr


# --- F3e: wrapped reports are parsed completely ---------------------------------------------

def test_f3e_real_wrapped_report_hiding_sorryax_is_rejected(workdir: Path) -> None:
    block = captured_report("EvilWrap")
    assert any("sorryAx" in line for line in block[1:]), "fixture must hide sorryAx on a continuation line"
    assert_rejected(run_gate(workdir, build_log(block)), "sorryAx")


@pytest.mark.parametrize("hidden", ["sorryAx", "Lean.ofReduceBool", "Lean.trustCompiler"])
def test_f3e_any_axiom_on_a_continuation_line_is_rejected(workdir: Path, hidden: str) -> None:
    block = [f"'{NAMES[0]}' depends on axioms: [propext,", f" {hidden},", " Quot.sound]"]
    assert_rejected(run_gate(workdir, build_log(block)), hidden)


def test_truncated_report_at_end_of_log_is_rejected(workdir: Path) -> None:
    complete = build_log().splitlines()[1:]  # every report except the first ...
    log = "\n".join(complete) + f"\n'{NAMES[0]}' depends on axioms: [propext,\n Classical.choice,\n"
    assert_rejected(run_gate(workdir, log), "Truncated")  # ... which is truncated at EOF


def test_truncated_report_followed_by_another_header_is_rejected(workdir: Path) -> None:
    block = [f"'{NAMES[0]}' depends on axioms: [propext,", " Classical.choice,"]
    assert_rejected(run_gate(workdir, build_log(block)), "Truncated")


def test_trailing_text_after_report_is_rejected(workdir: Path) -> None:
    block = report(NAMES[0], "[propext] sorryAx")
    assert_rejected(run_gate(workdir, build_log(block)), "Trailing text")


def test_unrecognised_report_format_is_rejected(workdir: Path) -> None:
    prefixed = [f"info: FoldBounds/AxiomAudit.lean:5:0: '{NAMES[0]}' depends on axioms: [propext]"]
    assert_rejected(run_gate(workdir, build_log(prefixed)), "Unparseable axiom report")


def test_empty_axiom_entry_is_rejected(workdir: Path) -> None:
    assert_rejected(run_gate(workdir, build_log(report(NAMES[0], "[propext, , Quot.sound]"))), "Empty")


# --- F3f: axiom names must match the authorised set exactly ---------------------------------

def test_f3f_real_guillemet_forged_axiom_is_rejected(workdir: Path) -> None:
    block = captured_report("EvilGuillemet")
    assert "«Classical.choice»" in block[0]
    assert_rejected(run_gate(workdir, build_log(block)), "«Classical.choice»")


@pytest.mark.parametrize(
    "lookalike",
    [
        "«Classical.choice»",
        "«propext»",
        "_private.FoldBounds.Evil.0.propext",
        "Classical.choice'",
        "propext✝",
        "FoldBounds.propext",
        "Classical.choice.extra",
        "Classical",
        "sorryAx",
    ],
)
def test_f3f_lookalike_axiom_names_are_rejected(workdir: Path, lookalike: str) -> None:
    axioms = f"[propext, {lookalike}, Quot.sound]"
    assert_rejected(run_gate(workdir, build_log(report(NAMES[0], axioms))), "Unexpected axioms")


def test_f3f_source_scan_rejects_guillemet_axiom_declaration(workdir: Path) -> None:
    (workdir / "FoldBounds" / "Forged.lean").write_text(
        "import Mathlib\naxiom «Classical.choice» : False\n", encoding="utf-8"
    )
    assert_rejected(run_gate(workdir), "Forged.lean")


@pytest.mark.parametrize(
    "declaration",
    ["axiom bad : False", "private axiom propext : False", "axiom\n  «x» : False", "axiom«x» : False"],
)
def test_source_scan_rejects_axiom_declarations(workdir: Path, declaration: str) -> None:
    (workdir / "FoldBounds" / "Forged.lean").write_text(f"import Mathlib\n{declaration}\n", encoding="utf-8")
    assert_rejected(run_gate(workdir), "Forged.lean")


def test_source_scan_allows_plural_axioms_in_prose(workdir: Path) -> None:
    (workdir / "FoldBounds" / "Prose.lean").write_text(
        "/-! Every headline theorem depends only on the standard axioms. -/\n", encoding="utf-8"
    )
    result = run_gate(workdir)
    assert result.returncode == 0, result.stderr


# --- existing guarantees are kept -------------------------------------------------------------

def test_duplicate_report_is_rejected(workdir: Path) -> None:
    log = build_log() + "\n".join(report(NAMES[0])) + "\n"
    assert_rejected(run_gate(workdir, log), "Duplicate")


def test_missing_report_is_rejected(workdir: Path) -> None:
    lines = build_log().splitlines()[1:]
    assert_rejected(run_gate(workdir, "\n".join(lines) + "\n"), "Missing")


def test_unaudited_extra_report_is_rejected(workdir: Path) -> None:
    log = build_log() + "\n".join(report("FoldBounds.NotAudited.foo")) + "\n"
    assert_rejected(run_gate(workdir, log), "extras")


def test_sorry_in_source_is_still_rejected(workdir: Path) -> None:
    (workdir / "FoldBounds" / "Hole.lean").write_text("theorem t : 0 = 1 := sorry\n", encoding="utf-8")
    assert_rejected(run_gate(workdir), "Hole.lean")
