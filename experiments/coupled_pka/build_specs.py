"""Download PDB files (if missing), write manifest and model specs (protocol sections 1-3)."""
import json, os, subprocess, sys, datetime
from common import *
from foldbounds.coupled_pka import model as M

os.makedirs(os.path.join(RAW, "pdb"), exist_ok=True)
os.makedirs(os.path.join(DATA, "specs"), exist_ok=True)
wild = pkad_wild()
man = {"retrieved": "2026-11-13", "files": []}
man["files"].append({"what": "PKAD-2 page (wildData array embedded)", "url": "https://compbio.clemson.edu/pkad-2/",
                     "sha256": sha256(os.path.join(RAW, "pkad2", "pkad2_index.html"))})
man["files"].append({"what": "PROPKA 3.5.1 sdist (propka/propka.cfg model_pkas)", "url": "PyPI propka==3.5.1 sdist (pip download --no-binary :all:)",
                     "sha256": sha256(os.path.join(RAW, "propka", "propka-3.5.1.tar.gz"))})
for p in SPEC_PROTEINS:
    f = os.path.join(RAW, "pdb", p + ".pdb")
    if not os.path.exists(f):
        subprocess.check_call(["curl", "-sS", "-f", "-o", f, f"https://files.rcsb.org/download/{p}.pdb"])
    man["files"].append({"what": f"PDB {p}", "url": f"https://files.rcsb.org/download/{p}.pdb", "sha256": sha256(f)})
    chains = chains_for(p, wild)
    sites, dropped = M.build_sites(open(f).read(), chains)
    spec = M.build_spec(sites, p)
    spec["chains"] = chains
    spec["dropped"] = dropped
    json.dump(spec, open(os.path.join(DATA, "specs", p + ".json"), "w"), separators=(",", ":"))
    print(p, chains, len(sites), "sites", len(spec["pairs_milli"]), "pairs", "dropped", dropped, M.spec_hash(spec)[:12])
man["spec_hashes"] = {}
json.dump(man, open(os.path.join(DATA, "MANIFEST.json"), "w"), indent=1)
