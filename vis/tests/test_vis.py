"""Adapter tests. The analytic backend is a TEST DOUBLE, never a CLI fallback."""
from __future__ import annotations

from pathlib import Path
import csv
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from optics_io import (SELECT_COLUMNS, INTENSITY_COLUMNS, load_optics, make_view,
                       region_mask, write_json)
from analysis_tools import weighted_metrics
from mie_reference import calculate_mie, incident_projection_weights
from compare_mie import run_comparison


def fixture(path: Path, *, nt=12, np_=8, jones=(1,0,0,0), pending=False,
            zero=False, error=False, sphere=True, radius=0.4, wavelength=700.0) -> None:
    meta = {
        "format": "rainbow_patch_optics_v1",
        "quantity": "regular_partial_model_angular_density_NOT_phase_function",
        "optical_complete": "false", "source_coverage_certified": "false",
        "focal_line_phase_applied": "false", "diffraction_applied": "false", "normalized": "false",
        "input_states": "one_coherent_Jones_state", "phasor_convention": "exp(+i*2*pi*q)",
        "interpolation": "transport_then_bilinear_field_and_path_then_propagation_phase",
        "transport": "shortest_great_circle_minimum_rotation_explicit_project_convention",
        "field_components": "s_perpendicular,p_outgoing_cross_s",
        "density_units": "input_field_squared*drop_unit_squared_per_sr",
        "physical_area_factor_mm2": format(radius**2, '.17g'),
        "missing_source_cells": "2", "no_outgoing_source_cells": "4",
        "radius_mm": repr(radius), "wavelength_nm": repr(wavelength),
        "exterior_index": "1", "interior_index": "1.5",
        "incident_grid": "129,129", "incident_grid_half_extent_drop": "1.01",
        "reference_distance_drop": "2", "outgoing_reference_distance_drop": "2",
        "incident_direction": "1,0,0", "incident_basis_x": "0,0,1",
        "incident_field": ','.join(str(x) for x in jones),
        "coefficients": "0,0,0,0,0,0,0,0" if sphere else "0,0,-0.01,0,0,0,0,0",
        "theta_count": str(nt), "phi_count": str(np_),
        "order": "theta_major_phi_minor", "angle_units": "radians",
    }
    i0 = sum(x*x for x in jones)
    with path.open('w', encoding='utf-8', newline='') as stream:
        for key, value in meta.items(): stream.write(f"# {key}={value}\n")
        writer = csv.DictWriter(stream, fieldnames=SELECT_COLUMNS)
        writer.writeheader()
        for i in range(nt):
            t = math.pi*(i+0.5)/nt
            wt = 2*math.pi/np_*(math.cos(math.pi*i/nt)-math.cos(math.pi*(i+1)/nt))
            for k in range(np_):
                phi = -math.pi+2*math.pi*(k+0.5)/np_
                index = i*np_+k
                row = dict.fromkeys(SELECT_COLUMNS, 0)
                row.update(direction_id=index, theta_index=i, phi_index=k, theta_rad=t, phi_rad=phi,
                           solid_angle_sr=wt, wx=math.cos(t), wy=-math.sin(t)*math.sin(phi),
                           wz=math.sin(t)*math.cos(phi), known_hits_complete=1, hit_count=2, evaluated_hits=2)
                for stage in ('incoherent','path'):
                    row[f'regular_partial_{stage}_s'] = i0*(1+0.2*math.cos(t))
                    row[f'regular_partial_{stage}_p'] = i0*0.5
                    row[f'regular_partial_{stage}_total'] = row[f'regular_partial_{stage}_s']+row[f'regular_partial_{stage}_p']
                if zero and index == 0:
                    for name in INTENSITY_COLUMNS: row[name] = 0.0
                if pending and index == 1:
                    row.update(flags=32, known_hits_complete=0, evaluated_hits=1, rejected_hits=1, refinement_hits=1)
                if error and index == 2:
                    row.update(flags=2, known_hits_complete=0, evaluated_hits=0, rejected_hits=2)
                    for name in INTENSITY_COLUMNS: row[name] = float('nan')
                writer.writerow(row)


class AnalyticDipoleBackend:
    """Only for unit-testing adapter conventions; not a replacement for Mie."""
    __version__ = "3.test-double"
    USE_JIT = False
    @staticmethod
    def intensities(m, d, wavelength, mu, *, n_env, norm, n_pole):
        assert norm == 'qsca' and n_pole == 0
        assert d > 0 and wavelength > 0 and m.real > 0 and n_env > 0
        # Integral of (iper+ipar)/2 = 1. pi*a^2 * these has qsca=1.
        return 3/(8*np.pi)*mu*mu, np.full(mu.shape, 3/(8*np.pi))
    @staticmethod
    def efficiencies(m, d, wavelength, *, n_env):
        return 1.0, 1.0, 1.0, 0.0


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.x = self.dir/'sphere_x.csv'; self.y = self.dir/'sphere_y.csv'
        fixture(self.x); fixture(self.y, jones=(0,0,1,0))
    def tearDown(self): self.tmp.cleanup()

    def test_shape_weights_and_units(self):
        g = load_optics(self.x)
        self.assertEqual(g.shape, (12,8))
        self.assertAlmostEqual(g.solid_angle.sum(), 4*np.pi, places=12)
        self.assertTrue(g.is_sphere)
        np.testing.assert_allclose(g.density('path','total'), g.fields['regular_partial_path_total']*0.16)

    def test_shuffled_rows_are_sorted_by_ids(self):
        text = self.x.read_text(); lines=text.splitlines()
        index = next(i for i,s in enumerate(lines) if not s.startswith('#'))
        self.x.write_text('\n'.join(lines[:index+1]+list(reversed(lines[index+1:])))+'\n')
        g=load_optics(self.x)
        self.assertAlmostEqual(g.fields['regular_partial_path_total'][0,0],1.5+0.2*math.cos(np.pi/24))

    def test_duplicate_id_rejected(self):
        s=self.x.read_text().splitlines(); s[-1]=s[-2]; self.x.write_text('\n'.join(s)+'\n')
        with self.assertRaisesRegex(ValueError,'Duplicate or missing'): load_optics(self.x)

    def test_bad_weight_rejected(self):
        self.x.write_text(self.x.read_text().replace('solid_angle_sr','wrong_weight'))
        with self.assertRaisesRegex(ValueError,'Missing optical'): load_optics(self.x)

    def test_non_grid_format_rejected(self):
        self.x.write_text(self.x.read_text().replace('rainbow_patch_optics_v1','rainbow_patch_queries_v1'))
        with self.assertRaisesRegex(ValueError,'Expected'): load_optics(self.x)

    def test_explicit_phase_normalized_input_rejected(self):
        self.x.write_text(self.x.read_text().replace('normalized=false','normalized=true'))
        with self.assertRaisesRegex(ValueError,'Normalized'): load_optics(self.x)

    def test_zero_input_rejected(self):
        fixture(self.x,jones=(0,0,0,0))
        with self.assertRaisesRegex(ValueError,'Jones norm'): load_optics(self.x)

    def test_status_masks_keep_zero_and_pending_separate(self):
        fixture(self.x,pending=True,error=True,zero=True)
        view=make_view(load_optics(self.x))
        self.assertTrue(view.accepted()[0,0])
        self.assertFalse(view.accepted()[0,1])
        self.assertTrue(view.accepted(show_partial=True)[0,1])
        self.assertFalse(view.accepted(show_partial=True)[0,2])
        self.assertEqual(view.values[0,0],0)
        self.assertTrue(np.isnan(view.values[0,2]))

    def test_unpolarized_is_normalized_intensity_average(self):
        fixture(self.y,jones=(0,0,2,0))
        a,b=load_optics(self.x),load_optics(self.y)
        v=make_view(a,orthogonal=b)
        self.assertTrue(v.is_unpolarized)
        np.testing.assert_allclose(v.values,a.density('path','total'))

    def test_nonorthogonal_pair_rejected(self):
        with self.assertRaisesRegex(ValueError,'orthogonal Jones'): make_view(load_optics(self.x),orthogonal=load_optics(self.x))

    def test_mismatched_pair_rejected(self):
        fixture(self.y,radius=0.5,jones=(0,0,1,0))
        with self.assertRaisesRegex(ValueError,'metadata mismatch'): make_view(load_optics(self.x),orthogonal=load_optics(self.y))

    def test_single_linear_polarization_projection(self):
        fixture(self.x,np_=6)
        v=make_view(load_optics(self.x))
        s,p=incident_projection_weights(v)
        np.testing.assert_allclose(s[0],np.sin(v.grid.phi)**2,rtol=0,atol=2e-15)
        np.testing.assert_allclose(s+p,1,rtol=0,atol=2e-15)
        self.assertAlmostEqual(s[0,4],1)  # phi=90 degrees: x input is s
        self.assertAlmostEqual(p[0,4],0)

    def test_circular_polarization_projection(self):
        fixture(self.x,jones=(1,0,0,1))
        s,p=incident_projection_weights(make_view(load_optics(self.x)))
        np.testing.assert_allclose(s,0.5); np.testing.assert_allclose(p,0.5)

    def test_mie_adapter_return_order_and_area(self):
        fixture(self.x,np_=6)
        v=make_view(load_optics(self.x))
        result,info=calculate_mie(v,backend=AnalyticDipoleBackend)
        fs,fp=incident_projection_weights(v)
        expected=3/8*v.grid.radius_mm**2*(fs+fp*np.cos(v.grid.theta[:,None])**2)
        np.testing.assert_allclose(result,expected,rtol=1e-14)
        self.assertAlmostEqual(info['diameter_nm'],800000)
        self.assertAlmostEqual(info['relative_refractive_index'],1.5)

    def test_mie_unpolarized_axisymmetry(self):
        v=make_view(load_optics(self.x),orthogonal=load_optics(self.y))
        result,_=calculate_mie(v,backend=AnalyticDipoleBackend)
        np.testing.assert_allclose(result,np.broadcast_to(result[:,0,None],result.shape))

    def test_non_sphere_mie_rejected(self):
        fixture(self.x,sphere=False)
        with self.assertRaisesRegex(ValueError,'--sphere'): calculate_mie(make_view(load_optics(self.x)),backend=AnalyticDipoleBackend)

    def test_tiny_nonzero_shape_not_silently_spherical(self):
        self.x.write_text(self.x.read_text().replace('coefficients=0,0,0,0,0,0,0,0','coefficients=0,0,1e-15,0,0,0,0,0'))
        self.assertFalse(load_optics(self.x).is_sphere)

    def test_absolute_metrics_do_not_fit_scale(self):
        b=np.array([[1.,2.],[3.,4.]])
        w=np.array([[1.,1.],[4.,4.]])
        m=weighted_metrics(2*b,b,w,np.ones((2,2),dtype=bool))
        self.assertAlmostEqual(m['relative_L1_solid_angle_weighted'],1)
        self.assertAlmostEqual(m['relative_L2_solid_angle_weighted'],1)
        self.assertAlmostEqual(m['integral_ratio'],2)
        self.assertFalse(m['scale_fitted'])

    def test_log_metric_does_not_floor_zero(self):
        m=weighted_metrics(np.array([0.,1.]),np.ones(2),np.ones(2),np.ones(2,dtype=bool))
        self.assertEqual(m['zero_model_positive_mie_samples'],1)
        self.assertEqual(m['log10_metric_sample_count'],1)
        self.assertAlmostEqual(m['relative_L1_solid_angle_weighted'],0.5)

    def test_region_uses_centers_and_full_weights(self):
        g=load_optics(self.x); m=region_mask(g,(120.,150.))
        self.assertEqual(m.sum(),16)
        self.assertAlmostEqual(g.solid_angle[m].sum(),2*np.pi*(np.cos(np.radians(120))-np.cos(np.radians(150))))

    def test_json_is_strict(self):
        write_json(self.dir/'j.json',{'x':float('nan'),'y':np.array([1.,np.inf])})
        self.assertEqual(json.loads((self.dir/'j.json').read_text()),{'x':None,'y':[1.,None]})

    def test_plot_cli_and_linear_npz(self):
        fixture(self.x,pending=True,error=True,zero=True)
        script=Path(__file__).resolve().parents[1]/'plot_optics.py'
        run=subprocess.run([sys.executable,str(script),str(self.x),'--out',str(self.dir/'out')],capture_output=True,text=True)
        self.assertEqual(run.returncode,0,run.stderr)
        self.assertGreater((self.dir/'out/intensity.png').stat().st_size,1000)
        saved=np.load(self.dir/'out/angular_data.npz',allow_pickle=False)
        self.assertEqual(saved['density_mm2_per_sr'][0,0],0)
        self.assertTrue(np.isnan(saved['density_mm2_per_sr'][0,2]))
        self.assertFalse(saved['known_hits_complete'][0,1])

    def test_comparison_cli_requires_explicit_partial_ack(self):
        script=Path(__file__).resolve().parents[1]/'compare_mie.py'
        run=subprocess.run([sys.executable,str(script),str(self.x),'--out',str(self.dir/'out')],capture_output=True,text=True)
        self.assertEqual(run.returncode,1)
        self.assertIn('--allow-partial-model',run.stderr)

    def test_comparison_pipeline_with_test_double(self):
        v=make_view(load_optics(self.x),orthogonal=load_optics(self.y))
        report=run_comparison(v,self.dir/'comparison',theta_range=(120,150),backend=AnalyticDipoleBackend)
        self.assertTrue(report['not_a_validation_certificate'])
        self.assertEqual(report['mie']['miepython_version'],'3.test-double')
        self.assertTrue((self.dir/'comparison/metrics.json').is_file())
        self.assertTrue((self.dir/'comparison/theta_azimuth_mean.png').is_file())


@unittest.skipUnless(importlib.util.find_spec('miepython'), 'Real miepython not installed; adapter test double is NOT Mie validation')
class RealMieTests(unittest.TestCase):
    def test_qsca_normalization_against_efficiency(self):
        import miepython as mie
        mu,w=np.polynomial.legendre.leggauss(256)
        m,d,wavelength=1.5,2.0,3.0
        par,per=mie.intensities(m,d,wavelength,mu,norm='qsca',n_pole=0)
        _,qsca,_,_=mie.efficiencies(m,d,wavelength)
        integral=2*np.pi*np.dot(w,0.5*(par+per))
        self.assertAlmostEqual(integral,float(qsca),delta=1e-10)

    def test_small_sphere_against_independent_dipole_formula(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.csv'
            # size parameter x=0.001; analytic Rayleigh is a deliberate TEST limit.
            wavelength=700.; a=0.001*wavelength/(2*np.pi*1e6)
            fixture(p,nt=24,np_=12,radius=a,wavelength=wavelength)
            view=make_view(load_optics(p)); result,info=calculate_mie(view)
            fs,fp=incident_projection_weights(view)
            k=2*np.pi/(wavelength*1e-6)
            alpha=(1.5**2-1)/(1.5**2+2)
            expected=k**4*a**6*alpha**2*(fs+fp*np.cos(view.grid.theta[:,None])**2)
            np.testing.assert_allclose(result,expected,rtol=3e-6,atol=0)
            self.assertTrue(info['miepython_version'].startswith('3.'))

if __name__ == '__main__': unittest.main()
