"""Optional OFFLINE audit of original source NEF/PDB bytes.

Fetch these ORIGINAL upstream files separately before running:
NMRExchangeFormat/NEF:data_1_1/PDBStat_developers_it3/Case_Amber/6nbn/6nbn.nef
NMRExchangeFormat/NEF:data_1_1/PDBStat_developers_it3/Case_Amber/6nbn/6nbn.ambpdb.pdb
"""
from pathlib import Path
import argparse,hashlib,json

ROOT=Path(__file__).resolve().parents[1]


def sha(data):
    return hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()


def coord(s):
    v=s.strip()
    sign=-1 if v.startswith('-') else 1
    s=v.lstrip('-')
    a,b=s.split('.')
    if len(b)!=3:
        raise ValueError('expected PDB 0.001 Å coordinates')
    return sign*(int(a)*1000+int(b))


def replay(nef,pdb,fixture):
    if sha(nef)!=fixture['source_nef_git_blob'] or sha(pdb)!=fixture['source_pdb_git_blob']:
        raise ValueError('ABSTAIN: original upstream Git blob hash mismatch')
    observed={}
    for line in pdb.decode('utf-8').splitlines():
        if not line.startswith('ATOM  '):
            continue
        if int(line[22:26])!=69:
            continue
        if line[17:20].strip()!='VAL' or line[16] not in (' ','A'):
            raise ValueError('unexpected Val69 alternate/residue')
        name=line[12:16].strip()
        if name in observed:
            raise ValueError('duplicate Val69 atom')
        observed[name]=tuple(coord(line[i:i+8]) for i in (30,38,46))
    expected={**fixture['fixed_heavy_atoms_val69'],**fixture['target_hydrogens_val69']}
    if any(observed.get(k)!=tuple(v) for k,v in expected.items()):
        raise ValueError('ABSTAIN: source coordinates do not match certified values')
    source=nef.decode('utf-8')
    seg=source.split('save_nef_distance_restraint_list_6nbn.mr',1)
    if len(seg)!=2:
        raise ValueError('NEF distance list missing')
    lines=seg[1].split('_nef_distance_restraint.upper_limit',1)
    if len(lines)!=2:
        raise ValueError('NEF schema unsupported')
    rows=[l.split() for l in lines[1].split('stop_',1)[0].splitlines() if l.strip() and l.lstrip()[0].isdigit()]
    target=[r for r in rows if len(r)==16 and r[1]=='470']
    if len(target)!=1:
        raise ValueError('restraint 470 absent or ambiguous')
    row=target[0]
    if (row[2]!='.' or row[3]!='A' or row[7]!='A' or
        [row[4],row[5],row[6],row[8],row[9],row[10]]
         !=['69','VAL','H','69','VAL','HB'] or
        row[14]!='1.605' or row[15]!='2.841'):
        raise ValueError('restraint 470 assignment/limit mismatch')
    return True


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('nef_file',type=Path)
    ap.add_argument('pdb_file',type=Path)
    args=ap.parse_args()
    source=json.loads((ROOT/'data/val69_local_noe_dual.json').read_text())
    if not replay(args.nef_file.read_bytes(),args.pdb_file.read_bytes(),source):
        raise SystemExit('ABSTAIN')
    print('Source SHA, 9 Val69 atom positions, unique NEF restraint 470: VERIFIED.')
