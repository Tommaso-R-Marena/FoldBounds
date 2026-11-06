"""Run strictly frozen GA95/GB95 -> GA98/GB98 NMR fold-geometry transfer.

No optimization, threshold selection, or conformer filtering on the test
structures. Every conformer is included, and all errors are reported.
"""
import argparse
import json
from pathlib import Path

from ga98_gb98_transfer import run


def main():
    parser=argparse.ArgumentParser()
    for name in ("ga95","gb95","ga98","gb98","reference"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    result=run(args.ga95.read_bytes(),args.gb95.read_bytes(),
               args.ga98.read_bytes(),args.gb98.read_bytes(),
               args.reference)
    serialized=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(serialized)
    print(serialized)


if __name__=="__main__":
    main()
