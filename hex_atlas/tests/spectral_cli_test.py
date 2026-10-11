"""End-to-end WB2008 -> deterministic CDF -> actual optional NF reader.

Runs the passed executable: under CUDA CTest this is the CUDA binary, not CPU.
The optional third argument is a spherical-hg-flow checkout; nothing is vendored.
"""
from pathlib import Path
import copy
import json
import os
import subprocess
import sys
import tempfile
import numpy as np

exe=Path(sys.argv[1]).resolve()
repo=Path(sys.argv[2]).resolve()
nf_root=Path(sys.argv[3]).resolve() if len(sys.argv)>3 else None
sys.path[:0]=[str(repo/'python'),str(repo/'hex_atlas/tools')]
from rainbow_dataset import PhaseRecord
from check_nf_compatibility import check_dataset
runner=repo/'hex_atlas/tools/run_spectrum.py'
composer=repo/'hex_atlas/python/compose_atlas.py'


def execute(args, expected=0):
    result=subprocess.run([str(v) for v in args],capture_output=True)
    if result.returncode!=expected:
        raise AssertionError(f'exit={result.returncode}, expected={expected}; stdout={result.stdout!r}; stderr={result.stderr!r}')
    return result


with tempfile.TemporaryDirectory(prefix='ice_spectral_') as temporary:
    root=Path(temporary)
    plan=dict(wavelengths_nm=[450,550,650],precompute=dict(radius_mm=.1,length_mm=.2,
        incident_theta_intervals=3,incident_phi_intervals=2,theta=12,phi=24,
        area_subdivisions=2,batch_rays=8,max_work_mib=32),
        composition=dict(inclinations_deg=[-90,0,35,90],orientation=dict(
            quadrature=dict(azimuth=4,tilt=4,spin=3),components=[
                dict(kind='horizontal_gaussian',sigma_deg=3,fraction=.9),dict(kind='isotropic',fraction=.1)])))
    path=root/'plan.json';path.write_text(json.dumps(plan),encoding='utf-8')
    args=[sys.executable,runner,'--exe',exe,'--plan',path,'--atlas-root',root/'atlases']
    execute([*args,'--out',root/'dry','--dry-run'])
    assert not (root/'atlases').exists() and not (root/'dry').exists()
    execute([*args,'--out',root/'dataset'])
    metas=[];arrays=[];original_times=[]
    for w,n in zip(plan['wavelengths_nm'],[1.3157,1.3110,1.3080]):
        p=root/'atlases'/f'lambda_{w}_nm'
        m=json.loads((p/'metadata.json').read_text())
        assert m['material']['interior_index']==n
        assert m['material']['automatic_ice_dispersion'] is True
        assert m['material']['wavelength_evaluated_nm']==w
        assert m['material']['absorption_model']=='ignored_k_equals_zero'
        assert m['material']['kappa_solver']==0
        arrays.append(np.load(p/'mass_mm2.npy'))
        original_times.append((p/'mass_mm2.npy').stat().st_mtime_ns)
        metas.append(m)
    assert not np.array_equal(arrays[0],arrays[1]) and not np.array_equal(arrays[1],arrays[2])
    index=json.loads((root/'dataset/batch_plan.json').read_text())
    assert index['complete'] and len(index['records'])==12
    cdf=[]
    for row in index['records']:
        with PhaseRecord(root/'dataset'/row['path'],validate=True) as r:
            m=r.metadata
            assert m['coordinate_contract_id']=='rainbow.phase_cdf.coordinates.v1'
            assert m['material']['wavelength_convention']=='vacuum_nm'
            assert m['azimuth_periodic'] is True
            assert m['shape_polar_axis']==[0,-1,0]
            np.testing.assert_allclose(m['hg']['axis_particle_frame'],r.frame[:,2],rtol=0,atol=0)
            np.testing.assert_allclose(row['condition'],[m['wavelength_nm'],-r.frame[1,2]],atol=1e-14)
            cdf.append(np.diff(r.phi_cdf)[:,None]*np.diff(r.theta_cdf,axis=1))
    assert np.max(np.abs(cdf[3]-cdf[4]))>1e-10, 'Wavelength did not change the saved CDF'
    execute([*args,'--out',root/'repeat','--reuse'])
    for w,t in zip(plan['wavelengths_nm'],original_times):
        assert (root/'atlases'/f'lambda_{w}_nm'/'mass_mm2.npy').stat().st_mtime_ns==t
    for row in index['records']:
        for filename in ['phi_cdf.npy','theta_given_phi_cdf.npy','u_edges.npy']:
            assert (root/'repeat'/row['path']/filename).read_bytes()==(root/'dataset'/row['path']/filename).read_bytes()
    if nf_root is not None:
        report=check_dataset(root/'dataset',nf_root)
        assert len(report['records'])==12
        print('UNMODIFIED RainbowReference: 12 records, source blob='+report['reference_module']['git_blob_sha1'])
    changed=copy.deepcopy(plan);changed['precompute']['length_mm']=.3
    path.write_text(json.dumps(changed),encoding='utf-8')
    execute([*args,'--out',root/'reject_cache','--reuse'],2)
    assert not (root/'reject_cache').exists()
    # No implicit promotion of --ior into an ice model; no ignored wavelength.
    basic=[exe,'--radius-mm','.1','--length-mm','.2','--incident-theta-intervals','1',
           '--incident-phi-intervals','1','--theta','12','--phi','24','--area-subdivisions','1']
    execute([*basic,'--wavelength-nm','550','--ior','1.31','--out',root/'bad1'],2)
    execute([*basic,'--wavelength-nm','1101','--out',root/'bad2'],2)
    execute([*basic,'--material','constant','--out',root/'bad3'],2)
    for p in ('bad1','bad2','bad3'):
        assert not (root/p).exists() and not (root/(p+'.part')).exists()
    cfg=copy.deepcopy(plan['composition']);cfg.pop('inclinations_deg');cfg['directions']=[[1,0,1]]
    cfg_path=root/'arbitrary.json';cfg_path.write_text(json.dumps(cfg),encoding='utf-8')
    execute([sys.executable,composer,'--atlas',root/'atlases/lambda_550_nm','--config',cfg_path,'--out',root/'bad_contract'],2)
    assert not (root/'bad_contract').exists()
    print('PASS: dispersive 3-wave/4-incident CDF, no ODF changes, unchanged cached arrays, explicit failure gates')
