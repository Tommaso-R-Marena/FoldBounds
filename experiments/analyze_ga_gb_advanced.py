import argparse,json
from pathlib import Path
from hashlib import sha256
from foldbounds.ga_gb_nmr import GA95,GB95,parse_pdb
from foldbounds.real_aromatic_geometry import analyze_advanced

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--ga',type=Path,required=True);p.add_argument('--gb',type=Path,required=True)
    p.add_argument('--out',type=Path)
    x=p.parse_args()
    b0=x.ga.read_bytes();b1=x.gb.read_bytes()
    ga=parse_pdb(b0,'2KDL',GA95);gb=parse_pdb(b1,'2KDM',GB95)
    summary=analyze_advanced(ga,gb)
    summary['coordinate_SHA256']={'2KDL':sha256(b0).hexdigest(),'2KDM':sha256(b1).hexdigest()}
    summary['method']='Retrospective NMR structural geometry; do not infer fold populations or energy from these models.'
    s=json.dumps(summary,indent=2,sort_keys=True)+'\n'
    if x.out:
        x.out.parent.mkdir(parents=True,exist_ok=True);x.out.write_text(s)
    print(s,end='')
if __name__=='__main__':main()