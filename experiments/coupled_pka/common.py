import hashlib, json, os, re, sys, datetime
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
RAW = os.path.join(ROOT, "data", "coupled_pka", "raw")
DATA = os.path.join(ROOT, "data", "coupled_pka")
PROTEINS = ["4PTI", "2LZT", "3RN3", "1UBQ", "1PGA", "1STN", "1A2P", "2RN2", "4HHB"]


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def pkad_wild():
    t = open(os.path.join(RAW, "pkad2", "pkad2_index.html")).read()
    return json.loads(re.search(r"var wildData\s*=\s*(\[\[.*?\]\]);", t, re.S).group(1))


def chains_for(pdb, wild):
    if pdb == "4HHB":
        return list("ABCD")
    return sorted({r[3] for r in wild if r[1] == pdb})


def load_spec(pdb):
    return json.load(open(os.path.join(DATA, "specs", pdb + ".json")))
