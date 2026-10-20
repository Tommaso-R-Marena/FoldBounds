"""Analyze authenticated experimental NMR coordinate files, no network calls.

Download from https://files.rcsb.org/download/{2KDL,2KDM}.pdb separately
and provide their exact byte streams. Do not treat them as equilibrium samples.
"""
import argparse,json
from pathlib import Path
from foldbounds.ga_gb_nmr import analyze_pair
from foldbounds.real_fold_odds import report_known_threshold

def main():
    a=argparse.ArgumentParser()
    a.add_argument('--ga',type=Path,required=True)
    a.add_argument('--gb',type=Path,required=True)
    a.add_argument('--out',type=Path)
    x=a.parse_args()
    result=analyze_pair(x.ga.read_bytes(),x.gb.read_bytes())
    result['experimental_thermodynamic_constraints']=report_known_threshold()
    text=json.dumps(result,indent=2,sort_keys=True)
    if x.out:
        x.out.parent.mkdir(parents=True,exist_ok=True);x.out.write_text(text+'\n')
    print(text)

if __name__=='__main__':main()