"""Opt-in exact 64-case source-only geometry+coupled GNM replay."""
import os,json
from pathlib import Path
import pytest
from gnm_geometry_rank_independent_replay import run


def actual_release():
    name=os.environ.get('FOLDBOUNDS_GEOM_RELEASE_DIR')
    if not name:pytest.skip('Pinned frozen reference archives not in hosted CI')
    folder=Path(name)
    return (folder/'FoldBounds_RNaseA_STAGE1_20261104_SEND_BACK.zip',
            folder,folder/'geometry_nogoods.json',folder/'gnm_geometry_rank_verified.json')


def test_exact_original_worst_case():
    s,base,geo,report=actual_release()
    # The original worst-case mask 55 is independently reproduced.
    r=run(s,base,geo,report,55,56,skip_existing_audit=True)
    assert r['status']=='PASS' and r['new_case_order_witnesses']==11870


def test_exact_all_64_protein_cases():
    if os.environ.get('FOLDBOUNDS_GEOM_FULL')!='1':
        pytest.skip('Full 64-case exact-rational replay is an explicit slow gate')
    s,base,geo,report=actual_release()
    r=run(s,base,geo,report)
    assert r['status']=='PASS' and r['new_case_order_witnesses']==986410


def test_wrong_rank_report_rejected(tmp_path):
    s,base,geo,report=actual_release()
    body=json.loads(report.read_text());body['cases'][55]['new_certified_orders']+=1
    out=tmp_path/'fraud.json';out.write_text(json.dumps(body))
    with pytest.raises(ValueError,match='forged case'):
        run(s,base,geo,out,55,56,skip_existing_audit=True)