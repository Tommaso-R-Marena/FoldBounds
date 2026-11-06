"""Reproduce the S669 real-measurement uncertainty experiment.

Uses pinned third-party published data; the CSV is NOT redistributed here.
One fold tests on a held-out collection of protein sequences, calibrates on a
different collection, and selects its predictor on three other folds.
"""
from __future__ import annotations
import argparse
import json
from fractions import Fraction
from pathlib import Path
from urllib.request import urlopen

from foldbounds.experimental_reciprocal_conformal import (
    SOURCE_REPO, SOURCE_COMMIT, SOURCE_PATH, SOURCE_BLOB_SHA1,
    git_blob_sha1, load_s669, evaluate, summarize,
)

SOURCE_URL=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}"


def main():
    p=argparse.ArgumentParser()
    group=p.add_mutually_exclusive_group(required=True)
    group.add_argument("--csv",type=Path)
    group.add_argument("--download",action="store_true")
    p.add_argument("--seed",type=int,default=20261105)
    p.add_argument("--repeat",type=int,default=1)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    if not 1<=a.repeat<=60:
        p.error("--repeat must be in [1,60]")
    if a.download:
        # Exact commit + source blob are both checked, preventing silent
        # postpublication changes to the input experiment.
        with urlopen(SOURCE_URL,timeout=45) as f:
            raw=f.read(2_000_000)
        if git_blob_sha1(raw)!=SOURCE_BLOB_SHA1:
            raise RuntimeError("upstream S669 hash did not match")
        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile(suffix=".csv",mode="wb") as f:
            f.write(raw);f.flush()
            rows=load_s669(f.name)
    else:
        rows=load_s669(a.csv)
    report={
        "research_scope":"observed real S669 predictions, not a new predictor or physical proof",
        "reference":{"repository":SOURCE_REPO,"commit":SOURCE_COMMIT,
                     "path":SOURCE_PATH,"git_blob_sha1":SOURCE_BLOB_SHA1},
        "rows":len(rows),
        "distinct_protein_sequences":len({x.wt for x in rows}),
        "pH_range":[str(min(x.ph for x in rows)),str(max(x.ph for x in rows))],
        "temperature_range_as_recorded":[str(min(x.temperature for x in rows)),
                                         str(max(x.temperature for x in rows))],
        "assumptions":["new whole protein groups exchangeable with calibration groups",
                       "source ddG convention unchanged",
                       "no leakage of test/calibration labels into predictor selection"],
        "runs":[]
    }
    for seed in range(a.seed,a.seed+a.repeat):
        entry={"seed":seed}
        for level,alpha in [("nominal_90",Fraction(1,10)),
                            ("nominal_80",Fraction(1,5))]:
            entry[level]={
              "fixed_group":summarize(evaluate(rows,seed=seed,alpha=alpha,strategy="fixed")),
              "selected_phys_group":summarize(evaluate(rows,seed=seed,alpha=alpha,strategy="selected_phys")),
              "selected_any_group":summarize(evaluate(rows,seed=seed,alpha=alpha,strategy="selected_any")),
              "naive_mutation_control":summarize(evaluate(rows,seed=seed,alpha=alpha,strategy="fixed",
                                                       calibration_unit="mutation")),
            }
        report["runs"].append(entry)
    report_json=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if a.output:
        a.output.write_text(report_json,encoding="utf-8")
    print(report_json)

if __name__=="__main__":
    main()
