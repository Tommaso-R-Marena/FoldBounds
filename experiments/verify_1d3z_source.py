"""Audit exact provenance of the small PDB 1D3Z C-alpha fixture.

This is an opt-in networked validation, not used by the exact mathematical
checker. Raw structural data are a public, deposited 1D3Z NMR ensemble.
The source is independently pinned by *Git blob SHA1* of 223801 raw bytes.
"""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path
from urllib.request import urlopen

OWNER = 'encore-similarity'
REPO = 'encore'
COMMIT = '2c99cb3fdf345c7679801378ec63634e1d007a5f'
PATH = 'examples/experimental_ubq_ensembles/filtered_1D3Z.pdb'
BLOB_SHA1 = 'a438d36816e15647177e32d2eafa6aa9a0782747'
URL = f'https://raw.githubusercontent.com/{OWNER}/{REPO}/{COMMIT}/{PATH}'
FIXTURE = (Path(__file__).resolve().parents[1] /
           'tests/fixtures/1d3z_model1_ca_milli.csv')


def extract_ca(raw: bytes):
    models = []
    coords = []
    ids = []
    for line in raw.decode('ascii').splitlines():
        if line.startswith('ENDMDL'):
            models.append((tuple(ids),tuple(coords)))
            ids,coords=[],[]
        if line.startswith('ATOM') and line[12:16].strip() == 'CA' and line[21:22]=='A':
            ids.append(int(line[22:26]))
            coords.append(tuple(round(float(line[a:b])*1000) for a,b in
                                ((30,38),(38,46),(46,54))))
    if ids:
        models.append((tuple(ids),tuple(coords)))
    return tuple(models)


def verify_source_bytes(raw: bytes) -> dict:
    h = hashlib.sha1(b'blob '+str(len(raw)).encode('ascii')+b'\0'+raw).hexdigest()
    if h != BLOB_SHA1:
        raise ValueError('source PDB Git blob hash mismatch')
    models = extract_ca(raw)
    if (len(models)!=10 or any(len(x)!=69 for _,x in models)
            or any(ids!=tuple(range(2,71)) for ids,_ in models)):
        raise ValueError('unexpected NMR model or CA identities')
    with FIXTURE.open(newline='',encoding='utf-8') as f:
        saved=tuple((int(row['residue']), int(row['x_milli']),
                    int(row['y_milli']),int(row['z_milli']))
                    for row in csv.DictReader(f))
    actual=tuple((ids, *xyz) for ids,xyz in zip(*models[0]))
    if saved!=actual:
        raise ValueError('committed fixture not identical to public source model 1')
    return {'git_blob_sha1':h,'nmr_models':len(models),
            'ca_residues':len(saved),'fixture_verified':True}


def main():
    import json
    with urlopen(URL,timeout=60) as r:
        raw=r.read()
    print(json.dumps(verify_source_bytes(raw),indent=2))


if __name__=='__main__':
    main()
