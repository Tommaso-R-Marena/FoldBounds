"""Reproduce the stricter 12-mer-blocked real S669 stability validation."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from urllib.request import urlopen

from foldbounds.experimental_reciprocal_conformal import (
    SOURCE_REPO, SOURCE_COMMIT, SOURCE_PATH, SOURCE_BLOB_SHA1,
    git_blob_sha1, load_s669,
)
from foldbounds.experimental_motif_blocks import (
    motif_components, evaluate_motif_blocked, summarize_motif,
)


def main():
    p=argparse.ArgumentParser()
    source=p.add_mutually_exclusive_group(required=True)
    source.add_argument("--csv",type=Path)
    source.add_argument("--download",action="store_true")
    p.add_argument("--k",type=int,default=12)
    p.add_argument("--seed",type=int,default=20261105)
    p.add_argument("--repeat",type=int,default=60)
    p.add_argument("--output",type=Path)
    args=p.parse_args()
    if not 1<=args.repeat<=60:
        p.error("repeat must be in [1,60]")
    if args.download:
        url=f"https://raw.githubusercontent.com/{SOURCE_REPO}/{SOURCE_COMMIT}/{SOURCE_PATH}"
        with urlopen(url,timeout=40) as response:
            raw=response.read(2_000_000)
        if git_blob_sha1(raw)!=SOURCE_BLOB_SHA1:
            raise RuntimeError("pinned source Git blob mismatch")
        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile(suffix=".csv",mode="wb") as tmp:
            tmp.write(raw)
            tmp.flush()
            records=load_s669(tmp.name)
    else:
        records=load_s669(args.csv)
    components=motif_components((r.wt for r in records),k=args.k)
    result={
        "scope":"shared exact k-mer motif BLOCKS, not inferred homologous families",
        "input_git_blob_sha1":SOURCE_BLOB_SHA1,
        "real_mutants":len(records),
        "wildtype_sequences":len({r.wt for r in records}),
        "motif_length":args.k,
        "motif_components":len(set(components.values())),
        "model_selection":"training blocks only; original published predictions",
        "distribution_assumption":"exchangeable whole motif blocks for rank theorem",
        "splits":[]
    }
    for seed in range(args.seed,args.seed+args.repeat):
        report={"seed":seed}
        for strategy in ("fixed","selected_phys"):
            folds=evaluate_motif_blocked(records,k=args.k,seed=seed,
                     alpha=Fraction(1,10),strategy=strategy)
            report[strategy]=summarize_motif(folds)
        if seed==args.seed:
            folds=evaluate_motif_blocked(records,k=args.k,seed=seed,
                 alpha=Fraction(1,10),calibration="mutation")
            report["naive_row_negative_control"]=summarize_motif(folds)
        result["splits"].append(report)
    output=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(output,encoding="utf-8")
    print(output)

if __name__=="__main__":
    main()
