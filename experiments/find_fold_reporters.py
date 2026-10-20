"""Rank experimentally grounded structural distance reporters; no binding/kinetics claims."""
import argparse,json
from pathlib import Path
from hashlib import sha256
from dataclasses import asdict
from foldbounds.ga_gb_nmr import GA95,GB95,parse_pdb
from foldbounds.real_fold_reporters import reporter_candidates,fret_ca_proxy_envelopes

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--ga',type=Path,required=True);p.add_argument('--gb',type=Path,required=True)
    p.add_argument('--out',type=Path)
    a=p.parse_args()
    ga0=a.ga.read_bytes();gb0=a.gb.read_bytes()
    ga=parse_pdb(ga0,'2KDL',GA95);gb=parse_pdb(gb0,'2KDM',GB95)
    ranked=reporter_candidates(ga,gb)
    result={'number_of_nonmutated_pairs_with_disjoint_12_to_40_A_NMR_distance_ranges':len(ranked),
            'top_ten':[{**asdict(r), 'fret_CA_proxy_R0_25_A':fret_ca_proxy_envelopes(r,25.0,0.0)} for r in ranked[:10]],
            'recommended_followup':['Measure solvent exposure of Q11/V21 and other candidates',
                                    'Model dye attachment and steric perturbation',
                                    'Validate calibrated fluorescence against NMR mutant state populations',
                                    'Do not infer actual state populations from NMR model frequencies'],
            'source_hashes':{'2KDL':sha256(ga0).hexdigest(),'2KDM':sha256(gb0).hexdigest()}}
    s=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if a.out:a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(s)
    print(s,end='')
if __name__=='__main__':main()