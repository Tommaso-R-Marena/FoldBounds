"""Physics-limit and target-isolation tests; no target PDB coordinates."""
import io
import json
import zipfile
import math
import sys
from pathlib import Path
import numpy as np
import pytest
import hybrid_reference_model as hm
import fb_thioredoxin_stage1 as stage


def mock_coords():
    return np.asarray([(3.5*i,3.5*j,2.1*math.sin(i+j))
                       for i in range(8) for j in range(8)],float)


def test_algebraic_alpha_limits():
    x=mock_coords(); ii=np.array([2,3,4,5]);jj=np.array([18,19,20,30])
    actual,_=hm.frozen_methods(x,ii,jj)
    assert np.max(abs(actual['ANM_weighted_inverse_square']-
                          hm.block_anm(x,ii,jj,0,True)))<1e-13
    assert np.max(abs(actual['ANM_ridge']-hm.block_anm(x,ii,jj,0,False)))<1e-13
    assert np.max(abs(actual['GANM_half_isotropic']-
                          hm.block_anm(x,ii,jj,0.5,True)))<1e-13
    assert np.all(actual['GANM_half_isotropic']>0)
    assert np.all(hm.block_anm(x,ii,jj,1,True)>0)
    assert np.all((actual['GNM_ANM_rank_fusion_75_25']>=0)&
                  (actual['GNM_ANM_rank_fusion_75_25']<=1))


def test_reject_bad_input():
    x=mock_coords()
    x[1]=x[0]
    with pytest.raises(ValueError,match='coincident'):
        hm.contacts(x)
    with pytest.raises(ValueError):
        hm.block_anm(mock_coords(),np.array([1]),np.array([2]),1.5)
    with pytest.raises(ValueError):
        hm.rankdata([0,np.nan])


class FakeAtom:
    def __init__(self,resi,xyz):
        self.resi=str(resi);self.resn='ALA';self.coord=list(map(float,xyz))

class FakeModel:
    def __init__(self,coords):self.atom=[FakeAtom(k+1,r) for k,r in enumerate(coords)]

class FakeCmd:
    def __init__(self):self.fetches=[];self.names=set();self.xyz=mock_coords()
    def get_names(self,who):return list(self.names)
    def delete(self,name):self.names.discard(name)
    def fetch(self,pdb,name,type,async_):
        assert pdb.lower()=='2trx', 'Target leakage: stage1 tried to fetch other PDB'
        self.fetches.append(pdb.lower());self.names.add(name)
    def count_atoms(self,name):return 500 if name in self.names else 0
    def get_model(self,selection,state):
        assert state==1
        assert 'chain A and name CA' in selection
        return FakeModel(self.xyz)


def test_pymol_reference_only_freeze(tmp_path,monkeypatch):
    downloads=tmp_path/'Downloads';downloads.mkdir()
    monkeypatch.setattr(Path,'home',lambda:tmp_path)
    cmd=FakeCmd()
    sys.modules['pymol']=type('FakePyMOL',(),{'cmd':cmd})()
    try:
        result=Path(stage.fb_thio_freeze())
    finally:
        sys.modules.pop('pymol',None)
    assert result.is_file()
    with zipfile.ZipFile(result) as z:
        assert z.testzip() is None
        names=set(z.namelist())
        m=json.loads(z.read('SHA256SUMS.json'))
        assert set(m)==names-{'SHA256SUMS.json'}
        from hashlib import sha256
        assert all(sha256(z.read(k)).hexdigest()==v for k,v in m.items())
        data=np.load(io.BytesIO(z.read('FROZEN_THIO_REFERENCE.npz')),allow_pickle=False)
        pro=json.loads(z.read('FROZEN_PROTOCOL.json'))
        assert pro['sealed_target']=='1XOA NMR all conformers'
        assert pro['primary']=='GANM_half_isotropic'
        assert pro['alpha_isotropic']==0.5
        assert len(data['residues'])==64
        assert len(data['pair_i'])>=100
        assert all('prediction_'+method in data for method in stage.ORIGINAL_BASELINES)
        assert 'prediction_GANM_half_isotropic' in data
        assert 'prediction_GNM_ANM_rank_fusion_75_25' in data
        assert b'def fb_thio_freeze' in z.read('fb_thioredoxin_stage1.py')
    assert cmd.fetches==['2trx']


if __name__=='__main__':
    raise SystemExit(pytest.main([__file__,'-q']))