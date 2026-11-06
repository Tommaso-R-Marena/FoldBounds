"""Run frozen 2009 GA/GB reporter codes on 2012 T25I reciprocal fold switches."""
import argparse
import json
from pathlib import Path
from ga98_t25i_followup import analyze_followup


def main():
    p=argparse.ArgumentParser()
    for name in ("ga95","gb95","t25i","t25i_l20a","reference"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--out",type=Path)
    a=p.parse_args()
    result=analyze_followup(
        a.ga95.read_bytes(),a.gb95.read_bytes(),a.t25i.read_bytes(),
        a.t25i_l20a.read_bytes(),a.reference
    )
    output=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(output)
    print(output)


if __name__=="__main__":
    main()
