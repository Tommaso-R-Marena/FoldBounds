# Matched-budget kernel benchmark harness

Fair comparison of `foldkernel` against `foldcheck` (v1), `foldcheck2` (v2),
unverified Python bound methods, and a pure-Python enumerator.

## Run

```bash
export PATH="$HOME/.elan/bin:$PATH"
cd /workspace && . .venv/bin/activate
lake build foldkernel    # foldcheck / foldcheck2 prebuilt on this branch
python experiments/kernel/bench/run_bench.py --quick   # ~10 min
python experiments/kernel/bench/run_bench.py --full    # ~75 min
python experiments/kernel/bench/run_bench.py --summarize-only
```

Results land in `results/*.jsonl` and `results/SUMMARY.md`.
Frozen instances: `manifest_dev.json`, `manifest_heldout.json`.
