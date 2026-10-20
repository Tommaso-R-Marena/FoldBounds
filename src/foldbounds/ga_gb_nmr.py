"""Reproducible structural descriptors for REAL NMR structures GA95 (2KDL), GB95 (2KDM).

Important: NMR *models* are conformers selected to represent structural uncertainty;
they are NOT Boltzmann samples. Contact occupancies across these models are
NOT equilibrium fold-state populations or predicted fold preferences.

Pure Python stdlib; refuses wrong or incomplete coordinates instead of fabricating them.
"""
from __future__ import annotations
from dataclasses import dataclass
from hashlib import sha256
from math import dist
from collections import defaultdict,Counter
from pathlib import Path
from typing import Iterable

RES3={"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
GA95="TTYKLILNLKQAKEEAIKELVDAGIAEKYIKLIANAKTVEGVWTLKDEIKTFTVTE"
MUTATIONS={20:("L","A"),30:("I","F"),45:("L","Y")}
def mutate(seq:str,changes:dict[int,tuple[str,str]])->str:
    a=list(seq)
    for i,(old,new) in changes.items():
        if a[i-1]!=old: raise ValueError(f"unexpected reference residue at {i}: {a[i-1]} != {old}")
        a[i-1]=new
    return "".join(a)
GB95=mutate(GA95,MUTATIONS)

@dataclass(frozen=True)
class NMRModel:
    seq: str
    ca: dict[int,tuple[float,float,float]]
    sidechains: dict[int,dict[str,tuple[float,float,float]]]

def parse_pdb(data:bytes, expected_id:str, expected_seq:str, expected_models:int=20)->tuple[NMRModel,...]:
    """Accept only expected real 56-residue chain A, unambiguous atom/residue records."""
    if not isinstance(data,bytes) or len(data)<5000:raise ValueError("PDB byte stream too short")
    text=data.decode("ascii")
    if f"{expected_id.upper()}" not in text[:15000].upper():
        raise ValueError("PDB entry identifier absent from header")
    if "EXPDTA" not in text or "NMR" not in text.upper()[:10000]:
        raise ValueError("expected an NMR structure, not an unrelated coordinate file")
    ca:dict[int,dict[int,tuple[float,float,float]]]=defaultdict(dict)
    side:dict[int,dict[int,dict[str,tuple[float,float,float]]]]=defaultdict(lambda:defaultdict(dict))
    codes:dict[int,dict[int,str]]=defaultdict(dict)
    current=1
    saw_model=False
    for line in text.splitlines():
        if line.startswith("MODEL "):
            saw_model=True
            current=int(line[10:14].strip());continue
        if not line.startswith(("ATOM  ","HETATM")) or len(line)<54:continue
        if not line.startswith("ATOM  "):continue
        if line[21]!="A":continue
        if line[16] not in " A":continue
        res3=line[17:20].strip(); idx_s=line[22:26].strip()
        if res3 not in RES3 or not idx_s:continue
        if line[26] not in " ":raise ValueError("insertion codes not supported")
        idx=int(idx_s)
        if idx<1 or idx>len(expected_seq):raise ValueError("unexpected residue numbering")
        xyz=tuple(float(line[j:j+8]) for j in (30,38,46))
        atom=line[12:16].strip()
        one=RES3[res3]
        if idx in codes[current] and codes[current][idx]!=one:
            raise ValueError("same residue has conflicting identities")
        codes[current][idx]=one
        if atom=="CA":
            if idx in ca[current] and ca[current][idx]!=xyz:raise ValueError("conflicting CA coordinates")
            ca[current][idx]=xyz
        # Sidechain atoms, exclude main chain and hydrogens
        if atom not in {"N","CA","C","O","OXT"} and not atom.startswith("H"):
            side[current][idx][atom]=xyz
    if not saw_model:raise ValueError("explicit NMR MODEL blocks required")
    if len(ca)!=expected_models:raise ValueError(f"expected {expected_models} NMR models, got {len(ca)}")
    result=[]
    for k in sorted(ca):
        if set(ca[k])!=set(range(1,len(expected_seq)+1)):
            raise ValueError(f"model {k}: incomplete C-alpha coordinates: {len(ca[k])}")
        seq="".join(codes[k][i] for i in range(1,len(expected_seq)+1))
        if seq!=expected_seq:
            mismatch=[i for i,(a,b) in enumerate(zip(seq,expected_seq),1) if a!=b]
            raise ValueError(f"model {k}: unexpected sequence, positions {mismatch}")
        result.append(NMRModel(seq,ca[k],dict(side[k])))
    return tuple(result)

def contact_pairs(m:NMRModel, cutoff:float=8.0, min_sep:int=4)->set[tuple[int,int]]:
    if not 0<cutoff<30 or not 1<=min_sep<=30:raise ValueError("invalid contact definition")
    return {(i,j) for i in m.ca for j in m.ca if j>=i+min_sep and dist(m.ca[i],m.ca[j])<cutoff}

def local_ca_geometry(m:NMRModel)->dict[str,float]:
    """Proxy for local secondary structure, explicitly NOT DSSP assignment."""
    n=len(m.seq)
    d3=[dist(m.ca[i],m.ca[i+3]) for i in range(1,n-2)]
    # C-alpha spacing over 3 positions: helix about 5 Å; extended strand about 10 Å
    return {"helix_like_i_to_i3_fraction":sum(4.5<=d<=6.6 for d in d3)/len(d3),
            "extended_i_to_i3_fraction":sum(9.2<=d<=11.6 for d in d3)/len(d3)}

def aromatic_centroid(m:NMRModel,i:int):
    if m.seq[i-1] not in "FYW":return None
    wanted={"F":("CG","CD1","CD2","CE1","CE2","CZ"),
            "Y":("CG","CD1","CD2","CE1","CE2","CZ"),
            "W":("CG","CD1","NE1","CE2","CD2","CE3","CZ2","CZ3","CH2")}[m.seq[i-1]]
    atoms=m.sidechains.get(i,{})
    points=[atoms[k] for k in wanted if k in atoms]
    if len(points)<4:return None
    return tuple(sum(p[ax] for p in points)/len(points) for ax in (0,1,2))

def analyze(models:tuple[NMRModel,...],pdb_id:str,raw_sha256:str)->dict:
    if not models:raise ValueError("no NMR models")
    n=len(models[0].seq)
    paircounts=Counter(p for m in models for p in contact_pairs(m))
    strong={p for p,count in paircounts.items() if count/len(models)>=0.75}
    mutation_sites=sorted(MUTATIONS)
    near={str(i):sum(1 for (a,b) in strong if a==i or b==i) for i in mutation_sites}
    fractions=[local_ca_geometry(m) for m in models]
    arom=[]
    for a,b in ((29,30),(29,45),(29,52),(30,45),(30,52),(45,52)):
        observed=[dist(x,y) for m in models if (x:=aromatic_centroid(m,a)) is not None and (y:=aromatic_centroid(m,b)) is not None]
        arom.append({"sites":[a,b],"valid_models":len(observed),
                     "mean_ring_centroid_distance_A":sum(observed)/len(observed) if observed else None,
                     "fraction_within_7_A":sum(x<7 for x in observed)/len(observed) if observed else None})
    return {"pdb_id":pdb_id,"coordinate_sha256":raw_sha256,"nmr_models":len(models),
            "residues":n,"sequence":models[0].seq,
            "contact_definition":"C-alpha distance <8 A, sequence separation >=4",
            "consensus_contacts_75pct":len(strong),"consensus_contact_pairs":[list(p) for p in sorted(strong)],
            "mutation_site_consensus_contact_degrees":near,
            "mean_helix_like_i_i3_fraction":sum(d["helix_like_i_to_i3_fraction"] for d in fractions)/len(fractions),
            "mean_extended_i_i3_fraction":sum(d["extended_i_to_i3_fraction"] for d in fractions)/len(fractions),
            "aromatic_ring_pairs":arom}

def compare(a:dict,b:dict)->dict:
    ca={tuple(p) for p in a["consensus_contact_pairs"]};cb={tuple(p) for p in b["consensus_contact_pairs"]}
    inter=ca&cb;union=ca|cb
    return {"nmr_model_count_each":[a["nmr_models"],b["nmr_models"]],
            "mutations":[{"position":i,"GA95":GA95[i-1],"GB95":GB95[i-1]} for i in range(1,57) if GA95[i-1]!=GB95[i-1]],
            "consensus_contacts_shared":len(inter),
            "consensus_contacts_only_GA95":len(ca-cb),"consensus_contacts_only_GB95":len(cb-ca),
            "consensus_contact_jaccard":len(inter)/len(union) if union else 1.0,
            "important_note":"NMR deposited models are structural representations, not equilibrium fold-population draws; this is RETROSPECTIVE structure analysis, not predictive scoring."}

def analyze_pair(ga_data:bytes,gb_data:bytes)->dict:
    ga=parse_pdb(ga_data,"2KDL",GA95)
    gb=parse_pdb(gb_data,"2KDM",GB95)
    a=analyze(ga,"2KDL",sha256(ga_data).hexdigest())
    b=analyze(gb,"2KDM",sha256(gb_data).hexdigest())
    return {"GA95":a,"GB95":b,"comparison":compare(a,b)}