"""Regenerate the pinned, unambiguous atom-pair subset from original NEF/PDB.

Usage (all files are explicitly supplied; NO automatic download):
 python scripts/replay_original_noe_sources.py \
    --fixture data/pinned_unambiguous_noe_2png_6nbn.json \
    --entry 2png --nef /path/to/2png.nef --pdb /path/to/2png.ambpdb.pdb

The original files are under NMRExchangeFormat/NEF, folder
 data_1_1/PDBStat_developers_it3/Case_Amber/<entry>/
The tool checks their original Git blob SHA-1, then reproduces *all* rows,
not just observed violations. This is a format-specific parser, not a general
NEF reader. Non-simple, ambiguous and atom-unmapped rows are excluded and
recorded in the reference manifest. Do not interpret exclusions as matches.
"""
from __future__ import annotations
import argparse
from collections import Counter
from decimal import Decimal
from hashlib import sha1
import json
from pathlib import Path
import re


def blob_sha(contents: bytes) -> str:
    """Git blob identity, not an unauthenticated text-file checksum."""
    return sha1(b'blob '+str(len(contents)).encode('ascii')+b'\0'+contents).hexdigest()


def numeric_mA(s: str) -> int:
    if not re.fullmatch(r'\d+(\.\d{1,3})?', s):
        raise ValueError('unsupported nonexact distance limit')
    milli = Decimal(s)*1000
    if milli != int(milli):
        raise ValueError('noninteger milliangstrom limit')
    return int(milli)


def parse_amber_pdb(contents: bytes) -> dict[tuple[int,str],tuple[str,tuple[int,int,int]]]:
    out = {}
    for line in contents.decode('utf-8').splitlines():
        if not line.startswith('ATOM  '):
            continue
        resid=int(line[22:26]); name=line[12:16].strip()
        key=(resid,name)
        if key in out:
            raise ValueError('duplicate coordinate mapping')
        coords=tuple(int(Decimal(line[j:j+8].strip())*1000) for j in (30,38,46))
        out[key]=(line[17:20].strip(),coords)
    return out


def parse_nef_simple(contents: bytes) -> tuple[list[list[str]],int]:
    src=contents.decode('utf-8')
    start=src.find('save_nef_distance_restraint_list')
    if start<0:
        raise ValueError('no distance-restraint saveframe')
    tail=src[start:]
    header='_nef_distance_restraint.upper_limit'
    pos=tail.find(header)
    if pos<0:
        raise ValueError('wrong or incomplete NEF header')
    rows=[]
    for line in tail[pos+len(header):tail.index('stop_',pos)].splitlines():
        fields=line.strip().split()
        if fields and fields[0].isdigit():
            if len(fields)!=16:
                raise ValueError('noncanonical data row')
            rows.append(fields)
    return rows,len(set(r[1] for r in rows))


def extract_records(nef: bytes,pdb: bytes) -> dict:
    rows,groups=parse_nef_simple(nef)
    atom=parse_amber_pdb(pdb)
    counts=Counter(r[1] for r in rows)
    entries=[]
    for r in rows:
        if counts[r[1]]!=1 or r[2]!='.' or r[3]!='A' or r[7]!='A':
            continue
        if not re.fullmatch(r'[A-Z][A-Z0-9]*',r[6]) or not re.fullmatch(r'[A-Z][A-Z0-9]*',r[10]):
            continue
        p=atom.get((int(r[4]),r[6]))
        q=atom.get((int(r[8]),r[10]))
        if p is None or q is None or p[0]!=r[5] or q[0]!=r[9]:
            continue
        entries.append([int(r[1]),int(r[4]),r[6],*p[1],
                        int(r[8]),r[10],*q[1],numeric_mA(r[14]),numeric_mA(r[15])])
    entries.sort()
    if len(set(r[0] for r in entries))!=len(entries):
        raise ValueError('duplicate included restraint ID')
    return {'total_noe_rows':len(rows),'total_noe_groups':groups,
            'single_row_groups':sum(1 for r in rows if counts[r[1]]==1),
            'explicit_unambiguous_atom_rows':sum(1 for r in rows if counts[r[1]]==1
                                                 and r[2]=='.' and r[3]=='A' and r[7]=='A'
                                                 and re.fullmatch(r'[A-Z][A-Z0-9]*',r[6])
                                                 and re.fullmatch(r'[A-Z][A-Z0-9]*',r[10])),
            'pdb_unique_atoms':len(atom),'records':entries}


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fixture',type=Path,required=True)
    ap.add_argument('--entry',choices=('2png','6nbn'),required=True)
    ap.add_argument('--nef',type=Path,required=True)
    ap.add_argument('--pdb',type=Path,required=True)
    a=ap.parse_args()
    ref=json.loads(a.fixture.read_text(encoding='utf-8'))['samples'][a.entry]
    nb,pb=a.nef.read_bytes(),a.pdb.read_bytes()
    if blob_sha(nb)!=ref['source_nef_sha'] or blob_sha(pb)!=ref['source_pdb_sha']:
        raise SystemExit('ABSTAIN: original NEF or PDB Git blob differs from pinned provenance')
    rebuilt=extract_records(nb,pb)
    for field,value in rebuilt.items():
        if ref[field]!=value:
            raise SystemExit('ABSTAIN: mismatch in '+field)
    print('PASS:',a.entry,'original Git blobs authenticated and all',len(rebuilt['records']),
          'selected restraints independently reproduced; excluded rows were not scored')

if __name__=='__main__':
    main()
