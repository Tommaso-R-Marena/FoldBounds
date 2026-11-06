"""Reproduce exact signed/robust sparse-distance certificates on PDB models."""
from __future__ import annotations
import json,sys
from pathlib import Path
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/"data"))
from experimental_ca_packed import STRUCTURES,decode
from foldbounds.sparse_distance_certificate import (
 make_certificate,make_signed_certificate,reconstruct_exact,reconstruct_signed_exact,
 robust_contacts,point_digest)

def main():
    results=[]
    for name in ("1A8O","1LCD","2BEG","7DDO"):
        record=STRUCTURES[name];p=decode(record["packed"])
        four=make_certificate(p)
        three=make_signed_certificate(p)
        if reconstruct_exact(four)!=p or reconstruct_signed_exact(three)!=p:
            raise AssertionError("missing exact coordinate audit")
        result=robust_contacts(four,epsilon_mA=10,cutoff_mA=8000)
        control=robust_contacts(make_certificate(p,(0,1,2,3)),epsilon_mA=10,cutoff_mA=8000)
        results.append({
            "pdb":name,"method":record["method"],"source_blob_sha":record["sha"],
            "ca_sites":len(p),"coordinate_sha256":point_digest(p),
            "anchor_indices":four.anchor_indices,
            "signed_distances":3*len(p)-6,"signed_orientation_bits":len(p)-3,
            "four_anchor_distances":4*len(p)-10,
            "all_pair_distances":len(p)*(len(p)-1)//2,
            "eligible_8A_pairs":result.all_pairs,
            "nominal_contacts":result.nominal_contacts,
            "guaranteed_contacts":result.guaranteed_contacts,
            "guaranteed_noncontacts":result.guaranteed_noncontacts,
            "undecided_pairs":result.undecided_pairs,
            "naive_first_four_undecided":control.undecided_pairs,
            "max_coordinate_error_A":float(result.max_coordinate_error_mA)/1000,
            "assumed_anchor_distance_error_A":0.01,
        })
    print(json.dumps(results,indent=2))
    return results

if __name__=="__main__":
    main()
