#!/usr/bin/env python3
"""Fixed particle/ODF: spectral atlas generation followed by deterministic CDFs.

No optical simplification, implicit downsampling, RNG, or CPU fallback. Completed
atlases can be reused only with matching arguments and executable/module bytes.
An unfinished .part is never adopted. Run --dry-run before a large job.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))
from compose_atlas import incidents, run as compose, write_json
from atlas_math import orientation_rule, sha256_file

REQUIRED = {'radius_mm', 'length_mm', 'incident_theta_intervals',
            'incident_phi_intervals', 'theta', 'phi', 'area_subdivisions',
            'batch_rays', 'max_work_mib'}
OPTIONAL = {'max_internal_hits', 'tail_tolerance', 'max_unresolved_fraction',
            'max_balance_error'}
INTEGER_KEYS = {'incident_theta_intervals', 'incident_phi_intervals', 'theta',
                'phi', 'area_subdivisions', 'batch_rays', 'max_work_mib', 'max_internal_hits'}


def read_plan(path: Path) -> dict:
    with path.open(encoding='utf-8') as f:
        plan = json.load(f)
    if set(plan) != {'wavelengths_nm', 'precompute', 'composition'}:
        raise ValueError('Plan requires only wavelengths_nm, precompute, composition')
    settings = plan['precompute']
    if not isinstance(settings, dict) or not REQUIRED <= set(settings) <= REQUIRED | OPTIONAL:
        raise ValueError('Missing/unknown precompute parameter')
    for k, v in settings.items():
        if k in INTEGER_KEYS:
            if type(v) is not int or v <= 0:
                raise ValueError(f'{k} must be a positive integer')
        elif type(v) not in (int, float) or not math.isfinite(v) or v < 0:
            raise ValueError(f'{k} must be a finite nonnegative number')
    if min(settings['radius_mm'], settings['length_mm']) <= 0:
        raise ValueError('Particle dimensions must be positive')
    waves = plan['wavelengths_nm']
    if (not isinstance(waves, list) or not waves or
            any(type(w) not in (int, float) or not math.isfinite(w) or not 200 <= w <= 1100 for w in waves) or
            len(set(waves)) != len(waves)):
        raise ValueError('Unique vacuum wavelengths in 200..1100 nm are required')
    if sorted(waves) != waves:
        raise ValueError('Wavelengths must be listed in increasing order (stable w-index)')
    cfg = plan['composition']
    if set(cfg) != {'inclinations_deg', 'orientation'}:
        raise ValueError('This fixed-condition runner requires inclination angles and one orientation model')
    directions = incidents(cfg)
    if len(set(cfg['inclinations_deg'])) != len(directions):
        raise ValueError('Duplicate incident conditions')
    orientation_rule(cfg['orientation'])  # validate before any expensive job
    work = 16 * settings['batch_rays'] * settings['theta'] * settings['phi']
    if work > settings['max_work_mib'] * 1024**2:
        raise ValueError('Per-ray GPU work buffer exceeds max_work_mib; reduce batch_rays explicitly')
    return plan


def name_of(wavelength: float) -> str:
    return 'lambda_' + format(float(wavelength), '.17g').replace('.', 'p') + '_nm'


def verify_atlas(path: Path, wave: float, settings: dict) -> dict:
    if path.name.endswith('.part'):
        raise ValueError('Incomplete atlas')
    m = json.loads((path / 'metadata.json').read_text(encoding='utf-8'))
    if m.get('complete') is not True or m.get('schema') != 'ice.single_particle_atlas.v1':
        raise ValueError('Incomplete/unknown cached atlas')
    mat = m['material']
    if (mat.get('model') != 'ice_ih_warren_brandt_2008_real_only' or
            mat.get('automatic_ice_dispersion') is not True or
            mat.get('wavelength_nm') != wave or mat.get('exterior_index') != 1):
        raise ValueError('Cached atlas is not the requested dispersive ice calculation')
    expected = {'incident_theta_intervals': settings['incident_theta_intervals'],
                'incident_phi_intervals': settings['incident_phi_intervals'],
                'theta_count': settings['theta'], 'phi_count': settings['phi']}
    if any(m.get(k) != v for k, v in expected.items()):
        raise ValueError('Cached angular grid differs from request')
    if (m['geometry']['circumradius_mm'] != settings['radius_mm'] or
            m['geometry']['full_length_mm'] != settings['length_mm'] or
            m['area_quadrature']['subdivisions'] != settings['area_subdivisions']):
        raise ValueError('Cached geometry/area quadrature differs from request')
    for filename in ('mass_mm2.npy', 'u_edges.npy'):
        if not (path / filename).is_file():
            raise ValueError('Missing cached atlas array')
    # compose.Atlas subsequently verifies every node's actual mass/area/tail.
    return m


def run(exe: Path, plan_path: Path, atlas_root: Path, out: Path,
        reuse: bool = False, dry_run: bool = False) -> None:
    plan = read_plan(plan_path)
    exe = exe.resolve()
    if not exe.is_file():
        raise ValueError('Precompute executable not found; build CPU or CUDA target explicitly')
    if out.exists() or out.with_name(out.name + '.part').exists():
        raise ValueError('Use a new CDF output directory; existing output/.part is never replaced')
    s = plan['precompute']
    raw = 8 * (s['incident_theta_intervals'] + 1) * (s['incident_phi_intervals'] + 1) * s['theta'] * s['phi']
    print(f"PLAN: {len(plan['wavelengths_nm'])} wavelengths x {len(plan['composition']['inclinations_deg'])} incidents")
    print(f"RAW ATLAS: {raw / 1024**3:.3f} GiB/wavelength; {raw * len(plan['wavelengths_nm']) / 1024**3:.3f} GiB total (plus headers/metadata/CDFs)")
    print(f"GPU RAY HISTOGRAM: {16*s['batch_rays']*s['theta']*s['phi']/1024**2:.1f} MiB (additional buffers required)")
    args = [item for key, value in sorted(s.items()) for item in ('--' + key.replace('_', '-'), str(value))]
    paths = [atlas_root / name_of(w) for w in plan['wavelengths_nm']]
    if dry_run:
        for w, path in zip(plan['wavelengths_nm'], paths):
            print(subprocess.list2cmdline([str(exe), '--material', 'ice', '--wavelength-nm', str(w), *args, '--out', str(path)]))
        print('DRY RUN: no files were created and no data were generated')
        return
    build = {'executable_sha256': sha256_file(exe)}
    module = exe.parent / 'modules/atlas_trace.fatbin'
    if module.is_file():
        build['adjacent_atlas_module_sha256'] = sha256_file(module)
    atlas_root.mkdir(parents=True, exist_ok=True)
    for index, (wave, path) in enumerate(zip(plan['wavelengths_nm'], paths)):
        receipt = path.with_name(path.name + '.job.json')
        expected = dict(schema='ice.spectral_atlas_job.v1', wavelength_nm=wave,
                        precompute=s, material='ice', build=build)
        if path.exists():
            if not reuse:
                raise FileExistsError('Atlas exists; use --reuse for validated matching jobs: ' + str(path))
            old = json.loads(receipt.read_text(encoding='utf-8'))
            if old.get('request') != expected or old.get('complete') is not True:
                raise ValueError('Cached atlas arguments/build differ; use another atlas-root')
            if old.get('metadata_sha256') != sha256_file(path / 'metadata.json'):
                raise ValueError('Cached metadata changed since the recorded job')
            verify_atlas(path, wave, s)
            print(f'REUSE {index+1}/{len(paths)}: wavelength_nm={wave}')
            continue
        if receipt.exists() or path.with_name(path.name + '.part').exists():
            raise FileExistsError('Unfinished job or stale receipt; inspect it before retrying: ' + str(path))
        log_path = path.with_name(path.name + '.log')
        command = [str(exe), '--material', 'ice', '--wavelength-nm', str(wave), *args, '--out', str(path)]
        started = time.perf_counter()
        with log_path.open('wb') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        if result.returncode != 0:
            raise RuntimeError(f'Precompute failed (exit={result.returncode}); log={log_path}; .part retained, no fallback')
        verify_atlas(path, wave, s)
        write_json(receipt, dict(complete=True, request=expected,
                                metadata_sha256=sha256_file(path / 'metadata.json'),
                                seconds=time.perf_counter()-started))
        print(f'ATLAS {index+1}/{len(paths)}: wavelength_nm={wave}; seconds={time.perf_counter()-started:.3f}', flush=True)
    payload = json.dumps(plan['composition'], sort_keys=True, ensure_ascii=True, allow_nan=False).encode('ascii')
    cfg_path = atlas_root / ('composition_' + hashlib.sha256(payload).hexdigest() + '.json')
    if cfg_path.exists():
        if json.loads(cfg_path.read_text(encoding='utf-8')) != plan['composition']:
            raise ValueError('Composition configuration cache mismatch')
    else:
        write_json(cfg_path, plan['composition'])
    compose(paths, cfg_path, out, output_contract='rainbow')
    print('COMPLETE: ' + ascii(str(out)) + '; check_nf_compatibility.py is the downstream acceptance test')


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--atlas-root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--reuse', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    try:
        run(a.exe, a.plan, a.atlas_root, a.out, a.reuse, a.dry_run)
    except Exception as e:
        print('SPECTRAL GENERATION FAILED: ' + ascii(str(e)), file=sys.stderr)
        raise SystemExit(2)

if __name__ == '__main__':
    main()
