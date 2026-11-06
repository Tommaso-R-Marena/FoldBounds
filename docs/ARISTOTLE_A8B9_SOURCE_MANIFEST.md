# Upstream Aristotle integration manifest (archival copy)

This preserves the upstream run’s file hashes and instructions **as received**.
Statements below that files were not yet pushed to GitHub describe the
Aristotle environment before ChatGPT integration, not the current PR status.

---

# FoldBounds: integration bundle for run a8b9ef96-c289-471a-b7fc-f8a7cea9be37

This environment had **no authenticated GitHub write access**: no `gh` CLI,
no GitHub token, and outbound connections to github.com were reset. So the
branch `integration/aristotle-frontier-a8b9ef96` and PR #26 were **not**
updated, and no commit SHA exists on GitHub. This bundle is meant for
integrating the results by hand.

## Contents

* `FoldBounds/`: an overlay tree. Paths are relative to the **repository
  root** of `Tommaso-R-Marena/FoldBounds`, using the original layout.
* `foldbounds-a8b9ef96-overlay.tar.gz` / `.zip`: the same overlay tree
  together with this manifest and `rerun_logs/`.
* `rerun_logs/`: outputs from the re-verification run on 2026-10-16
  (see below).

## How to apply (do not overwrite newer `main` work)

Files tagged **NEW** do not exist in the pre-run snapshot. Copy them in
unchanged.

Files tagged **MERGE** existed before this run and were edited by it. The
current `main` (and so the integration branch) may have changed them since.
**Merge them by hand; do not overwrite them.** The edits needed in each are:

| File | What this run changed |
|---|---|
| `FoldBounds.lean` | Adds 11 `public import` lines: `LatticeEquivalence`, `LatticeContacts`, `Frontier`, `FrontierDegree`, `FrontierFamily`, `Transfer`, `FrontierCheck`, `FrontierSoundness`, `FrontierPartition`, `FrontierFlat`, `FrontierCLI`. `FrontierTests` is built through the lakefile glob and is not imported. |
| `FoldBounds/AxiomAudit.lean` | Adds imports of the new modules and new `#print axioms` lines, giving 117 audited declarations in total. Keep any lines `main` has added; the audit script rejects duplicates. |
| `FoldBounds/LatticeModel.lean` | The version used by the verified new modules. Diff it against the branch. If `main` changed it, check that the new modules still build. |
| `lakefile.toml` | Adds the `[[lean_exe]] foldcheck` target (`srcDir = "scripts"`, `root = "foldcheck"`). The lib glob is `FoldBounds.+`. |
| `scripts/check_lean_sources.py` | Source-hygiene and axiom-log audit, run as `python3 scripts/check_lean_sources.py <axiom-log>`. Keep whichever version is stricter. |
| `README.md`, `docs/FORMAL_ASSURANCE_BOUNDARY.md` | Sections on the frontier checker, the trust boundary and integration notes. Merge these sections into the current text. Do not replace descriptions of newer `main` features. |

The Python exporter `src/foldbounds/lean_export.py` and the experiments
depend only on these names from the existing package:
`foldbounds.lattice.{ROOT, Certificate, Node, Path, build_frontier,
exact_spectrum, verify_certificate}` (they use `Certificate.sequence/.q/
.event_at_least/.frontier` and `Node.path/.kind`), and
`foldbounds.witnesses.propose_dfs_witnesses`. If `main` renamed any of
these, adapt the exporter. It is untrusted, so changing it has no effect on
any Lean theorem.

`scripts/run_tests_offline.py` and `scripts/offline_pytest/` are a small
fallback used only because pytest was unavailable here. Under real pytest
they are optional and can be dropped.

## Commands to run after merging (from the repository root)

```
lake build                       # whole library, including FrontierTests
lake build foldcheck             # compiled runner (first build ~6 min)
lake env lean FoldBounds/AxiomAudit.lean > axiom.log
python3 scripts/check_lean_sources.py axiom.log
python -m pytest -q
python3 experiments/lean_differential.py --lake-root . --runner scripts/foldcheck.lean
python3 experiments/lean_frontier_benchmark.py --lake-root . --exe "$PWD/.lake/build/bin/foldcheck" --exact-max 16 --out /tmp/bench.json
```

(`--exe` must be an absolute path, because the script runs it with
`cwd = --lake-root`.)

## Re-verification performed in this environment (2026-10-16)

These checks ran on the pre-`main`-merge snapshot from the original run. They
were **not** run on the merged branch, which was not reachable from here.
Toolchain: Lean `v4.28.0`, Mathlib `v4.28.0`
(`8f9d9cff6bd728b17a24e163c9402775d9e6a365`), Python 3.11.14.

| Check | Result |
|---|---|
| `lake build FoldBounds foldcheck` | pass: "Build completed successfully (16084 jobs)" |
| Search for `sorry` / `admit` / `axiom` / `@[implemented_by]` in the Lean sources | none |
| `lake env lean FoldBounds/AxiomAudit.lean` | 117 declarations. 116 depend only on a subset of `[propext, Classical.choice, Quot.sound]`; `conformation_nodup` depends on no axioms |
| `check_lean_sources.py` on that log | "Source audit: 26 Lean modules, no forbidden proof bypasses" and "Axiom audit: 117 declarations, standard axioms only" |
| `ExecutableTests`, `LatticeModelTests`, `FrontierTests` | each exits 0 (the frontier tests use `decide +kernel`) |
| Python suite (offline runner; pytest is **not installed** and there is no network) | passed=632 failed=0 |
| `lean_differential.py` | 18 sequences with matching semantics; 2160 certificate checks, 0 failures; 51 adversarial certificates, 0 accepted; "TOTAL PROBLEMS: 0" |
| `lean_frontier_benchmark.py --exact-max 16` | "rows: 78 problems: 0"; all 54 rows with exact values enclose them. All bounds, node counts and certificate sizes are **identical** to the committed `experiments/results/lean_frontier_benchmark.log`; only the timings differ |

Differential tests are testing, not proof. They do not establish a verified
correspondence between Python and Lean.

## Trust boundary (unchanged)

The theorems are `FrontierCheck.check_sound`, `check_coupled_probability_sound`,
`checkFlat_sound` and the related results. They state: *if the Lean function
`check` returns `some b` on a certificate value, then `b` encloses the exact
Lean-defined Z, event mass, complement mass and P(A).*

The following are **not** verified:
* the text parser `FrontierCLI` and the runner `scripts/foldcheck.lean`;
* compiled execution, which trusts the Lean compiler; kernel `decide` tests do not;
* the Python proposer and exporter;
* the physical energy model;
* any complexity claim; the cost analysis is informal.

The Python algorithms (matching, isoperimetric, finite-memory transfer,
history-aware pruning, monotone families, two-turn witnesses, endpoint
geometry, adaptive verification) are **not** Lean-verified.

## File list with SHA-256 (paths relative to the repository root)

```
NEW d10a959ab5d65b6229632bdab7724a874b1e25b2ddaa443cdb6612fc82106599  FoldBounds/LatticeEquivalence.lean
NEW 2ff9c461b36126d17c7b6030619114803eef14b27c0d41af027fb1a893e74bd8  FoldBounds/LatticeContacts.lean
NEW 8f0392375241b8e838a50330bb330c335513fdef27e40e5698504eba4894b7f4  FoldBounds/Frontier.lean
NEW c0a84ba08ab19d074d81dabab0b842f606157f7b4d3d5f97e1b72b67cd3a16dc  FoldBounds/FrontierDegree.lean
NEW 3bb63f983dfe9cc7c7443d910eb8e65fe40460f00aa8cadde4fc677f2059566e  FoldBounds/FrontierFamily.lean
NEW c076d192706db0af02e791bce8cfb6579cbc545cf6c7cb7b5ebb2a2960cc8ff4  FoldBounds/Transfer.lean
NEW fdb26e4175653b4ed35dc13d02a9a2b3d11b3a1a122b9fe94ac771b1979d14a1  FoldBounds/FrontierCheck.lean
NEW 2116c862f6cf28a4ac843593a5ef78fc8e9e32bdce6797afae1cb9479f4761fd  FoldBounds/FrontierSoundness.lean
NEW 5a94275940737bf7dd87c8a783e493a25c4541ab921d5911871792ba5240f52e  FoldBounds/FrontierPartition.lean
NEW c117730c8f1343c5a0ee715a78325354bab49034df213aaa58586505645180a9  FoldBounds/FrontierFlat.lean
NEW 0599dc08354c7d48107f2b224bd3e4c348e9dbed0ef7689e6eefda5acb164131  FoldBounds/FrontierCLI.lean
NEW fc6251a9deaf204f05eb250a6120d8e3a31a5fec03f1fce4d5cf434f24e9501c  FoldBounds/FrontierTests.lean
MERGE 5e99a1cfdc84672924d07221139df113c4c9eb9ae3f3f61961f124ab2f89d9f6  FoldBounds/LatticeModel.lean
MERGE f493ceb3524feca75bd418d3d523fe700c7acd9e1954abf3ec92e30f146124ec  FoldBounds.lean
MERGE 82118577ccabe4106e0edb706e6785e9b571f616964be9d90a53c5e196036e6d  FoldBounds/AxiomAudit.lean
MERGE d5084f299e9947bcecf481c434193adff5a9c7dc76f5a519e6db4cd75aebff20  lakefile.toml
NEW 4aa9814a450612c5a2755086419c95216df7673964279fca62951f2fb26ebca0  scripts/foldcheck.lean
MERGE 65950a64babfc80d1b223c9c88bdf0e26cb866a413b3d308b814be412184c357  scripts/check_lean_sources.py
NEW 8b7f058f0c3e1eae4c97827474bf60e81435dfdad229c049b7f824d54f84fffb  scripts/run_tests_offline.py
NEW 91d52c1f956dbf877d13e0f24d439142d18cb7d148f067f3fc2893472a5e2db7  scripts/offline_pytest/pytest.py
NEW 554a81928ec1979956aa40534b7e17357c5435be4f03232d304984a8926d9c11  src/foldbounds/lean_export.py
NEW edb8444443c9284f166cf92f19a017e6f34f2f8b8ae2167110fa74488ca3f4e7  tests/test_lean_export.py
NEW fa659bcd98b2a801bca5398a72ff848282085c014d987c02a1e8fae0223c698b  experiments/lean_differential.py
NEW 742377c57c4552bb09203fa60d81f03b4792cdd8577d0e98f37302a20e35299c  experiments/lean_frontier_benchmark.py
NEW 08335d57a7a8f0dce80e80a7f8c8501be226a9aa6cf5f6e8174de3be824f9317  experiments/results/axiom_audit.log
NEW e604304b9701d9bba5fa47e2dcb5121e8111e5be42545101629b29e159688552  experiments/results/lean_differential.log
NEW 7acc584b7e49224faabd1c0fec6edb729f6abcac4bf77c1cbb70d17fb1171a9e  experiments/results/lean_frontier_benchmark.json
NEW 28dbc240ee60f88a8d50e6b1cc5f4070dba794c1581dd842330b21b12d6ab804  experiments/results/lean_frontier_benchmark.log
NEW 26b7ffd515fffeb4a107fd63c9ab6f2ac4cd45bdfd5a8b970efb6508f83322f6  experiments/results/lean_frontier_benchmark_prefamily_partial.log
NEW 964800aabbcf03a46334f3b514580d27411ae7ee008fd4c8a639be5f95b46060  docs/ARISTOTLE_LATTICE_VERIFICATION_REPORT.md
NEW a0f5ea9bed03a65477c721dcdbaa64b314d011165f5f514070867846043a6257  docs/FORMALLY_VERIFIED_FRONTIER_CERTIFICATE.md
NEW 46b63139f4f79aac8a45651cb3729b4a579addd86ee5abcac234e7e2b798a97a  docs/NONENUMERATIVE_BOUNDS_RESEARCH_REPORT.md
NEW bfd068742f16848b5b88105db268ad3022f5c24d1dfe71b904856270e32f2381  docs/ARISTOTLE_VERIFICATION_RUN_REPORT.md
MERGE 582625c2bb13102b811e70bf3897e47cb907802dc4bee63598378d1ac4ed181f  docs/FORMAL_ASSURANCE_BOUNDARY.md
MERGE 9ed44f012177d397b9e410faf42a10686ea143559a254b8358b77d4388afa8c0  README.md
```
