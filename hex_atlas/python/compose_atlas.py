#!/usr/bin/env python3
"""Reuse fixed single-particle atlas(es) to produce multiple rainbow-compatible CDFs.

Composition is CPU/NumPy by design. Both CPU and CUDA optical precomputation
feed this same offline operator. No trace executable is launched, no random
orientation nodes are generated, and no output angular resolution is changed.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys
import time
from typing import Any
import numpy as np
from atlas_math import (Atlas, compose_mass, coupling_for_rule, make_cdf, orientation_rule,
                        sha256_file, unit)


def write_json(path: Path, value: Any) -> None:
    with path.open('x',encoding='utf-8',newline='\n') as f:
        json.dump(value,f,ensure_ascii=True,allow_nan=False,indent=2)
        f.write('\n')


def write_npy(path: Path, value: np.ndarray) -> None:
    with path.open('xb') as f:
        np.save(f,np.ascontiguousarray(value,dtype='<f8'),allow_pickle=False)


def staging_directory(out: Path) -> Path:
    if out.name in ('', '.', '..') or out.name.endswith('.part') or out.exists():
        raise ValueError('Output must be a new named directory, not .part')
    part = out.with_name(out.name+'.part')
    out.parent.mkdir(parents=True,exist_ok=True)
    part.mkdir()  # exclusive: do not remove or adopt existing incomplete runs
    return part


def incidents(config: dict[str,Any]) -> list[tuple[np.ndarray,dict[str,Any]]]:
    have = set(config) & {'inclinations_deg','directions'}
    if len(have)!=1:
        raise ValueError('Specify exactly one of inclinations_deg or directions')
    result = []
    if 'inclinations_deg' in config:
        for value in config['inclinations_deg']:
            a = float(value)
            if not math.isfinite(a) or not -90<=a<=90:
                raise ValueError('Inclination must be in [-90,90] degrees')
            k = np.array([math.cos(math.radians(a)),-math.sin(math.radians(a)),0.])
            if abs(a)==90:
                k = np.array([0., -math.copysign(1,a), 0.])
            result.append((unit(k),{'inclination_deg':a}))
    else:
        for k in config['directions']:
            result.append((unit(k),{}))
    if not result or len(result)>10000:
        raise ValueError('Require 1..10000 incident conditions')
    return result


def metadata(a: Atlas, expected: float, residual: float, report: dict[str,Any],
             f: np.ndarray, plan_hash: str, rule_meta: dict[str,Any], condition: dict[str,Any],
             source_hash: str, azimuth_integration: dict[str,Any], output_contract: str = 'rainbow') -> dict[str,Any]:
    g = report['g']; material = dict(a.meta['material']); geom = a.meta['geometry'].copy()
    geom.update(projected_area_mm2=expected,projected_area_m2=expected*1e-6,
        projected_area_method='deterministic_orientation_quadrature_of_exact_facet_projection')
    total = report['resolved_mass_mm2']
    quality = {k:v for k,v in report.items() if k not in ('g','resolved_mass_mm2')}
    quality.update(source_coverage_certified=False,angular_convergence_certified=False,
        orientation_convergence_certified=False,incident_atlas_convergence_certified=False,
        continuous_wave_solution_certified=False,
        maximum_accepted_l1_mass_error=1e-10,maximum_accepted_lost_probability_mass=1e-12)
    # Material labels from old M3 atlases are not promoted to dispersive ice.
    material.setdefault('wavelength_convention','vacuum_nm')
    common = dict(source_repository='https://github.com/Somatic-Cell/ice-crystal-phase',
        source_base_commit='1cf6758d28473678ed141654f8c0c9b4f21197de',
        azimuth_periodic=True,particle_frame_handedness='right',particle_up_axis=[0,1,0],shape_polar_axis=[0,-1,0],
        phi_range_rad=[-math.pi,math.pi],
        transport_reference_frame='ensemble_reference; not an individual crystal body frame',
        shape_polar_axis_semantics='oriented transport conditioning axis -y; body c-axis is geometry.body_long_axis')
    if output_contract=='rainbow':
        k=f[:,2]
        if abs(k[2])>2e-14 or k[0]<-2e-14:
            raise ValueError('RainbowReference requires canonical ki=(cos(alpha),-sin(alpha),0); no implicit rotation of a non-axisymmetric ensemble')
        common['incident_inclination_degrees']=math.degrees(math.atan2(-k[1],max(0.,k[0])))
    common['nf_condition_axes']=['wavelength_nm','dot([0,-1,0],incident_propagation)']
    common['geometry_and_orientation_fixed_within_dataset']=True
    return dict(**common,schema='rainbow.phase_cdf.numpy.v2',complete=True,
        generator_version='ice.deterministic_atlas.m4.spectral_nf.v1',backend='numpy_deterministic_atlas_convolution',
        source_sha256=source_hash,dtype='<f8',order='C',theta_count=a.nt,phi_count=a.np,
        density_measure='solid_angle_sr',cell_model='constant_density_per_spherical_cell',
        coordinates='u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)',cdf_axis_order=['phi_cell','theta_edge'],
        normalization='integral_p_domega_equals_one',input_polarization='unpolarized',
        direction_convention='physical_propagation',theta_zero='forward',theta_pi='backward',
        coordinate_contract_id='rainbow.phase_cdf.coordinates.v1' if output_contract=='rainbow' else 'ice.phase_cdf.ensemble_coordinates.v1',
        sampling_frame_layout='rows_xyz_columns_e0_e1_k',sampling_frame_map='local_column_to_particle_column' if output_contract=='rainbow' else 'local_column_to_ensemble_reference_column',
        sampling_frame_columns=f.tolist(),query_frame_axis=f[:,2].tolist(),query_frame_e0=f[:,0].tolist(),
        wavelength_nm=material['wavelength_nm'],interior_index=material['interior_index'],exterior_index=material['exterior_index'],
        geometry=geom,material=material,
        cross_sections=dict(model='geometric_nonabsorbing',scattering_mm2=expected,extinction_mm2=expected,absorption_mm2=0,
            scattering_m2=expected*1e-6,extinction_m2=expected*1e-6,absorption_m2=0,wave_scattering_cross_section_computed=False),
        single_scattering_albedo=1,stage='incoherent',diffraction_policy=dict(computed=False,selected_output_includes_diffraction=False),
        inter_path_interference=False,tir_relative_phase_retained=True,
        orientation=dict(representation='deterministic_SO3_number_quadrature',rule=rule_meta,
            model_parameter_policy='explicit_user_parameters_not_aerodynamic_prediction',
            plan_sha256=plan_hash,randomness=False,plan_embedded=False,plan_location='dataset_root',
            azimuth_integration=azimuth_integration),
        atlas=dict(files_sha256=a.hashes,incident_theta_intervals=a.it,incident_phi_intervals=a.ip,
            angular_reconstruction='A_normalized_bilinear_theta_phi_in_incident_spherical_frame',
            symmetry='D6h_including_reflections',transform='phi_shift_or_phi_reflection_then_shift',
            azimuth_remap='exact_overlap_of_piecewise_constant_equal_width_cells',
            analytical_average_of_continuous_optics=False,no_retracing=True),
        transport_audit=dict(escaped_mm2=total,unresolved_mm2=residual,unresolved_fraction=residual/expected,
            balance_error=(total+residual)/expected-1,cdf_normalized_to_resolved_escaped_mass=True,
            residual_estimator='same_covariant_incident_interpolation_as_response'),
        storage_processing=dict(gaussian_sigma_degrees=0,coarsening_tv=0,
            mass_conservative_aggregation=False,g_source=g,g_stored=g),
        hg=dict(g=g,axis_particle_frame=f[:,2].tolist(),method='first_moment_of_saved_cell_pdf',target='saved_cdf',
            cosine_convention='dot(incident_propagation,outgoing_propagation)'),
        condition=condition,quality=quality)


def run(atlas_paths: list[Path], config_path: Path, out: Path, output_contract: str = 'rainbow') -> None:
    with config_path.open(encoding='utf-8') as f:
        config = json.load(f)
    if set(config) - {'orientation','inclinations_deg','directions'} or 'orientation' not in config:
        raise ValueError('Unknown/missing composition configuration keys')
    incident_list = incidents(config)
    atlases: list[Atlas] = []
    part = None
    try:
        for path in atlas_paths:
            atlases.append(Atlas(path))
        if not atlases:
            raise ValueError('At least one atlas is required')
        first = atlases[0]
        geometry = first.meta['geometry']
        rotations,weights,rule_meta=orientation_rule(config['orientation'])
        rule_meta=dict(rule_meta,
            vertical_uniform_azimuth_policy={
                'version':'analytic_uniform_alpha_at_exact_vertical_v1',
                'reference_arrays':'discrete_reference_not_literal_nodes_at_analytic_conditions',
                'eligible_components':['horizontal_gaussian','isotropic'],
                'fixed_components_unchanged':True,'nonvertical_conditions_unchanged':True})
        if output_contract not in ('rainbow','ensemble'):
            raise ValueError('Unknown output contract')
        if output_contract=='rainbow':
            for k,_ in incident_list:
                if abs(k[2])>2e-14 or k[0]<-2e-14:
                    raise ValueError('Noncanonical incident direction is not RainbowReference-compatible; select --output-contract ensemble explicitly')
        for a in atlases[1:]:
            if (a.it,a.ip,a.nt,a.np)!=(first.it,first.ip,first.nt,first.np) or a.meta['geometry']!=geometry or not np.array_equal(a.edges,first.edges) or not np.array_equal(a.area,first.area):
                raise ValueError('Batch atlases must share geometry, angular grids, and areas')
        material_keys=('model','exterior_index','absorption','birefringence')
        for a in atlases[1:]:
            if any(a.meta['material'].get(k)!=first.meta['material'].get(k) for k in material_keys):
                raise ValueError('A spectral batch must use one material model and exterior medium')
            if a.meta['transport']!=first.meta['transport']:
                raise ValueError('A spectral batch must use one transport/audit policy')
        wavelengths = [float(a.meta['material']['wavelength_nm']) for a in atlases]
        if len(set(wavelengths))!=len(wavelengths) or any(not math.isfinite(w) or w<=0 for w in wavelengths):
            raise ValueError('Unique positive wavelength per atlas required')
        if (first.it+1)*(first.ip+1)*2*first.np*8 > 512*1024**2:
            raise ValueError('Coupling buffer would exceed the documented 512 MiB bound')
        part = staging_directory(out)
        write_json(part/'orientation_model.json',rule_meta)
        write_npy(part/'orientation_rotations.npy',rotations)
        write_npy(part/'orientation_weights.npy',weights)
        plan_hashes = {n:sha256_file(part/n) for n in ('orientation_model.json','orientation_rotations.npy','orientation_weights.npy')}
        import hashlib
        plan_hash = hashlib.sha256(json.dumps(plan_hashes,sort_keys=True).encode('ascii')).hexdigest()
        source_files=['compose_atlas.py','atlas_math.py']
        source_hash = hashlib.sha256(''.join(sha256_file(Path(__file__).with_name(n)) for n in source_files).encode('ascii')).hexdigest()
        plan: dict[str,Any] = dict(schema='ice.deterministic_phase_dataset.v1',complete=False,
            nf_condition_order=['wavelength_nm','dot([0,-1,0],incident_propagation)'],
            output_contract=output_contract,per_record_schema='rainbow.phase_cdf.numpy.v2',condition_axes=['incident','wavelength'],
            incident_directions=[k.tolist() for k,_ in incident_list],wavelengths_nm=wavelengths,
            orientation_plan_sha256=plan_hash,orientation_files_sha256=plan_hashes,
            input_config=config,atlases=[a.hashes for a in atlases],source_sha256=source_hash,
            computation='deterministic_offline_mass_convolution_no_raytracing',records=[])
        write_json(part/'generation_plan.json',plan)
        for i,(ki,condition) in enumerate(incident_list):
            start=time.perf_counter()
            c,expected,fr,azimuth_integration = coupling_for_rule(first,rotations,weights,ki,rule_meta)
            for w,a in enumerate(atlases):
                mass,residual = compose_mass(a,c)
                total=math.fsum(mass.flat)
                if abs((total+residual)/expected-1)>2e-10 or residual/expected>a.meta['transport']['max_unresolved_fraction']+2e-15:
                    raise ValueError('Composed mass and unresolved area fail the transport audit')
                phi,theta,report=make_cdf(mass,a.edges)
                relative=Path('records')/f'i{i:04d}'/f'w{w:04d}'
                record=part/relative;record.mkdir(parents=True)
                write_npy(record/'phi_cdf.npy',phi);write_npy(record/'theta_given_phi_cdf.npy',theta);write_npy(record/'u_edges.npy',a.edges)
                m=metadata(a,expected,residual,report,fr,plan_hash,rule_meta,condition,source_hash,azimuth_integration,output_contract)
                write_json(record/'metadata.json',m)
                record_hashes={name:sha256_file(record/name) for name in ('metadata.json','phi_cdf.npy','theta_given_phi_cdf.npy','u_edges.npy')}
                plan['records'].append(dict(incident_index=i,wavelength_index=w,path=relative.as_posix(),
                    condition=[wavelengths[w],-float(ki[1])],ki=ki.tolist(),wavelength_nm=wavelengths[w],g=report['g'],files_sha256=record_hashes,
                    azimuth_integration=azimuth_integration))
                print(f"record i={i} w={w} mass=1 g={report['g']:.17g} A={expected:.17g}")
            print(f"incident_complete i={i} offline_seconds={time.perf_counter()-start:.6f}")
        plan['complete']=True
        write_json(part/'batch_plan.json',plan)
        if out.exists():
            raise FileExistsError('Destination appeared while composing')
        part.rename(out)
    except Exception as e:
        if part is not None and part.exists() and not (part/'failure.json').exists():
            write_json(part/'failure.json',dict(complete=False,error=str(e)))
        raise
    finally:
        for a in atlases:
            a.close()


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--atlas',nargs='+',type=Path,required=True)
    p.add_argument('--config',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--output-contract',choices=['rainbow','ensemble'],default='rainbow',help='rainbow: strict NF canonical incidents; ensemble: explicit legacy arbitrary-vector records')
    args=p.parse_args()
    try:
        run(args.atlas,args.config,args.out,args.output_contract)
    except Exception as e:
        print(f"compose_atlas failed: {e!a}",file=sys.stderr)
        raise SystemExit(2)

if __name__=='__main__':
    main()
