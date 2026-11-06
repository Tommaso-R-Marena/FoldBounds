"""Independently verify source PDB bytes and the derived exact CA coordinate data.

Usage: python experiments/audit_public_pdb_sources.py 1A8O /path/1A8O.pdb

The caller supplies original PDB bytes; no silent downloads. Check the
original *Git blob object SHA-1* and re-extract first-model chain-A CA
coordinates using Decimal to avoid float rounding. These are publicly
deposited structural models, not raw experimental measurements.
"""
from __future__ import annotations
import hashlib,sys
from pathlib import Path
from decimal import Decimal
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"data"))
from experimental_ca_packed import STRUCTURES,decode

def extract_ca(original:bytes):
    try:lines=original.decode("ascii").splitlines()
    except UnicodeDecodeError as e:raise ValueError("not ASCII PDB") from e
    model=1;seen=set();sites=[]
    for line in lines:
        if line.startswith("MODEL"):
            try:model=int(line[10:14].strip())
            except ValueError as e:raise ValueError("bad MODEL record") from e
        if (model!=1 or line[:6]!="ATOM  " or
            line[12:16].strip()!="CA" or line[21:22]!="A" or
            line[16:17] not in (" ","A")):continue
        residue=line[22:27].strip()
        if residue in seen:continue
        p=[]
        for at in (30,38,46):
            x=Decimal(line[at:at+8].strip())*1000
            if x!=x.to_integral_value():
                raise ValueError("CA coordinate not 0.001 Angstrom quantized")
            p.append(int(x))
        sites.append(tuple(p));seen.add(residue)
    return tuple(sites)

def audit(name:str,original:bytes):
    name=name.upper()
    if name not in STRUCTURES:raise ValueError("unknown dataset")
    meta=STRUCTURES[name]
    blob=b"blob "+str(len(original)).encode()+b"\0"+original
    sha=hashlib.sha1(blob).hexdigest()
    if sha!=meta["sha"]:raise ValueError("original PDB blob SHA mismatch")
    points=extract_ca(original)
    if points!=decode(meta["packed"]):raise ValueError("source coordinates mismatch")
    return len(points)

if __name__=="__main__":
    if len(sys.argv)!=3:raise SystemExit("usage: audit_public_pdb_sources.py PDB_ID original.pdb")
    print(f"Verified {sys.argv[1]}: {audit(sys.argv[1],Path(sys.argv[2]).read_bytes())} CA coordinates")
