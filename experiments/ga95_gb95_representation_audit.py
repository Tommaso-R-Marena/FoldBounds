"""Reproduce RCSB-sequence coarse-graining no-go audit without PDB coordinates."""
import json
from dataclasses import asdict
from foldbounds.ga_gb_indistinguishability import audit_experimental_ga_gb
if __name__=='__main__':
    print(json.dumps(asdict(audit_experimental_ga_gb()),indent=2,sort_keys=True))