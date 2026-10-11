"""End-to-end deterministic atlas -> multiple conditions -> original reader."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

exe=Path(sys.argv[1]).resolve();repo=Path(sys.argv[2]).resolve()
sys.path.insert(0,str(repo/'python'))
from rainbow_dataset import PhaseRecord
compose=repo/'hex_atlas/python/compose_atlas.py'

def run(args,expected=0):
    p=subprocess.run([str(a) for a in args],capture_output=True)
    if p.returncode!=expected:
        raise RuntimeError(f'returncode={p.returncode} expected={expected}; stdout={p.stdout!r}; stderr={p.stderr!r}')
    return p

def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

with tempfile.TemporaryDirectory(prefix='ice_atlas_cli_') as tmp:
    root=Path(tmp)
    base=[exe,'--radius-mm','.1','--length-mm','.2','--material','constant','--ior','1.31',
          '--incident-theta-intervals','3','--incident-phi-intervals','2',
          '--theta','12','--phi','24','--area-subdivisions','2']
    paths=[]
    for wave in (550,650):
        a=root/f'atlas{wave}';run([*base,'--wavelength-nm',str(wave),'--out',a]);paths.append(a)
    before={str(a):digest(a/'mass_mm2.npy') for a in paths}
    config={'inclinations_deg':[0,35,80],'orientation':{
        'quadrature':{'azimuth':4,'tilt':4,'spin':3},'components':[
            {'kind':'horizontal_gaussian','sigma_deg':10,'fraction':.75},
            {'kind':'isotropic','fraction':.25}]}}
    cfg=root/'compose.json';cfg.write_text(json.dumps(config),encoding='utf-8')
    cmd=[sys.executable,compose,'--atlas',*paths,'--config',cfg]
    # Composer gets no executable path and cannot call the optical trace backend.
    out=root/'dataset';run([*cmd,'--out',out]);out2=root/'repeat';run([*cmd,'--out',out2])
    plan=json.loads((out/'batch_plan.json').read_text());assert plan['complete']
    assert len(plan['records'])==6
    values=[]
    for r in plan['records']:
        path=out/r['path']
        with PhaseRecord(path,validate=True) as record:
            mass,g=record.mass_and_g();assert abs(mass-1)<1e-12
            assert record.metadata['orientation']['randomness'] is False
            assert record.metadata['atlas']['no_retracing'] is True
            assert abs(record.metadata['transport_audit']['balance_error'])<2e-10
            values.append(g)
        for name in ('phi_cdf.npy','theta_given_phi_cdf.npy','u_edges.npy'):
            assert (path/name).read_bytes()==(out2/r['path']/name).read_bytes()
        assert set(p.name for p in path.iterdir())=={'metadata.json','phi_cdf.npy','theta_given_phi_cdf.npy','u_edges.npy'}
    assert max(values)-min(values)>.01,'Incident condition was lost during composition'
    assert all(digest(a/'mass_mm2.npy')==before[str(a)] for a in paths)
    run([*cmd,'--out',out],2) # no replacement of an existing dataset
    run([*base,'--wavelength-nm','550','--out',paths[0]],2)
    run([*base,'--wavelength-nm','550','--out',root/'fail','--max-internal-hits','1'],2)
    assert not (root/'fail').exists() and (root/'fail.part').exists()
    print('PASS: 2 atlases, 3 incidents, 6 unchanged-reader records; repeat identical; no retrace; failure gates')
