#!/usr/bin/env python3
"""Exercise the UNMODIFIED downstream RainbowReference, not a substitute reader.

No PyTorch import or training is required. Module is loaded directly from the
user's spherical-hg-flow checkout; its exact SHA256/Git blob identity is logged.
Inputs are read-only. CDFs are never repaired, rescaled, or given a PDF floor.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
import numpy as np


def checked_reference(nf_root: Path):
    path=nf_root/'src/phaseflow/rainbow.py'
    raw=path.read_bytes()
    name='_downstream_rainbow_reference'
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None:
        raise ImportError(f'Cannot import downstream reference from {path}')
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module.RainbowReference,dict(path=str(path.resolve()),sha256=hashlib.sha256(raw).hexdigest(),
        git_blob_sha1=hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest())


def check_dataset(dataset: Path,nf_root: Path) -> dict:
    Reference,source=checked_reference(nf_root)
    plan=json.loads((dataset/'batch_plan.json').read_text(encoding='utf-8'))
    if plan.get('complete') is not True or not plan.get('records'):
        raise ValueError('Incomplete or empty dataset')
    results=[];root=dataset.resolve()
    for entry in plan['records']:
        path=(root/entry['path']).resolve()
        if root not in path.parents or path.suffix=='.part':
            raise ValueError('Unsafe/incomplete record path')
        with Reference(path) as ref:
            condition=np.array([entry['wavelength_nm'],-entry['ki'][1]])
            if not np.allclose(ref.condition,condition,atol=2e-12,rtol=0):
                raise AssertionError('Downstream conditioning differs')
            mass,g=ref.mass_and_g()
            if abs(mass-1)>1e-12 or abs(g-entry['g'])>5e-12:
                raise AssertionError('Downstream mass/moment differs')
            # Independent checks at every stored cell center, mapped to NF frame.
            u=0.5*(ref.u_edges[:-1]+ref.u_edges[1:])
            ph=-np.pi+2*np.pi*(np.arange(ref.n_phi)+.5)/ref.n_phi
            radius=2*np.sqrt(u*(1-u))
            x=radius[None,:]*np.cos(ph[:,None]);y=radius[None,:]*np.sin(ph[:,None])
            z=np.broadcast_to(1-2*u,x.shape)
            local=np.stack([x,y,z],axis=-1)
            nf=local@ref.source_to_nf.T
            actual=ref.log_prob(nf)
            joint=np.diff(ref.phi_cdf)[:,None]*np.diff(ref.theta_cdf,axis=1)
            density=joint/((4*np.pi/ref.n_phi)*np.diff(ref.u_edges))[None,:]
            positive=density>0
            if not np.array_equal(np.isfinite(actual),positive):
                raise AssertionError('Downstream support mismatch (no floors)')
            error=float(np.max(np.abs(actual[positive]-np.log(density[positive]))))
            if error>2e-11:
                raise AssertionError(f'Downstream density mismatch: {error}')
            # Fixed interior quantiles, not a Monte Carlo transport calculation.
            a=(np.arange(31)+.37)/31;b=(np.arange(29)+.61)/29
            uu,vv=np.meshgrid(a,b,indexing='ij')
            directions,log_pdf=ref.sample(np.stack([uu,vv],axis=-1))
            if not np.isfinite(log_pdf).all() or not np.allclose(np.linalg.norm(directions,axis=-1),1,atol=3e-14,rtol=0):
                raise AssertionError('Invalid NF direction samples')
            if not np.array_equal(ref.log_prob(directions),log_pdf):
                raise AssertionError('NF sample/eval PDF mismatch')
            results.append(dict(path=entry['path'],condition=ref.condition.tolist(),mass=mass,g=g,
                maximum_cell_center_log_pdf_error=error,sampled_directions=int(uu.size),
                source_metadata_repository=ref.metadata.get('source_repository')))
    return dict(complete=True,reference_module=source,records=results,
        tested='reader + condition + exact cell PDF + deterministic CDF inversion; not NF training',
        unchanged_input=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--dataset',type=Path,required=True);p.add_argument('--nf-root',type=Path,required=True)
    p.add_argument('--out',type=Path)
    args=p.parse_args()
    try:
        report=check_dataset(args.dataset,args.nf_root)
        if args.out:
            with args.out.open('x',encoding='utf-8') as f:json.dump(report,f,indent=2,ensure_ascii=True,allow_nan=False)
        print(f"PASS RainbowReference: {len(report['records'])} records; reader blob={report['reference_module']['git_blob_sha1']}")
    except Exception as e:
        print('NF compatibility FAILED: '+ascii(str(e)),file=sys.stderr);raise SystemExit(2)

if __name__=='__main__': main()
