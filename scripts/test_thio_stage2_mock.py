"""Stage2 mock: prove no target access until external SHA validation."""
import sys
import zipfile
import json
import math
from hashlib import sha256
from pathlib import Path
import numpy as np
import pytest
import fb_thioredoxin_stage1 as one
import fb_thioredoxin_stage2 as two

class Atom:
    def __init__(self,i,v):self.resi=str(i);self.coord=list(map(float,v));self.resn='ALA'
class Model:
    def __init__(self,xyz):self.atom=[Atom(i+1,v) for i,v in enumerate(xyz)]
class FakeCmd:
    def __init__(self):
        self.fetched=[];self.objects=set();self.xyz=np.asarray([(3.5*i,3.5*j,2.1*math.sin(i+j))
                         for i in range(8) for j in range(8)],float)
    def get_names(self,who):return list(self.objects)
    def delete(self,name):self.objects.discard(name)
    def fetch(self,pdb,name,type,async_):self.fetched.append(pdb);self.objects.add(name)
    def count_atoms(self,name):return 400 if name in self.objects else 0
    def count_states(self,name):return 20
    def get_model(self,select,state):
        assert 'chain A and name CA' in select
        if '1xoa' in select:
            # All-state target trajectories synthetic, not experimental.
            perturb=(state-10.5)*.006
            xyz=self.xyz + perturb*np.cos(np.arange(len(self.xyz)))[:,None]*np.array([1,1.8,-.7])
        else: xyz=self.xyz
        return Model(xyz)


def test_full_stage2_replay_mock(tmp_path,monkeypatch):
    downloads=tmp_path/'Downloads';downloads.mkdir()
    monkeypatch.setattr(Path,'home',lambda:tmp_path)
    cmd=FakeCmd()
    sys.modules['pymol']=type('FakePyMOL',(),{'cmd':cmd})()
    try:
        frozen=Path(one.fb_thio_freeze())
        digest=sha256(frozen.read_bytes()).hexdigest()
        before=len(cmd.fetched)
        with pytest.raises(ValueError,match='SHA256'):
            two.fb_thio_evaluate('0'*64)
        assert len(cmd.fetched)==before,'Stage2 touched NMR before verifying SHA!'
        output=Path(two.fb_thio_evaluate(digest))
        with zipfile.ZipFile(output) as z:
            assert z.testzip() is None
            meta=json.loads(z.read('FROZEN_EVALUATION.json'))
            assert meta['stage1_archive_sha256']==digest
            assert meta['conformers']==20
            assert meta['primary']=='GANM_half_isotropic'
            assert meta['status'] in ('PASS','FAIL')
            assert len(meta['method_metrics'])==8
            manifest=json.loads(z.read('SHA256SUMS.json'))
            assert set(manifest)==set(z.namelist())-{'SHA256SUMS.json'}
            assert all(sha256(z.read(k)).hexdigest()==v for k,v in manifest.items())
        assert cmd.fetched==['2trx','1xoa']
    finally:
        sys.modules.pop('pymol',None)


if __name__=='__main__':raise SystemExit(pytest.main([__file__,'-q']))