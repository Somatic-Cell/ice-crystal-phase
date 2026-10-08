from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO/'tools'))
from phase_analysis import Table, angular_summary, compare_tables, cone_mass, load_table, main


def table(mass, u, *, name='synthetic'):
    w = np.array(mass, dtype=float)
    return Table(Path(name), {'coordinate_contract_id': 'ice.phase_cdf.ensemble_coordinates.v1',
                             'sampling_frame_map': 'local_column_to_ensemble_reference_column',
                             'direction_convention': 'physical_propagation',
                             'theta_zero': 'forward', 'theta_pi': 'backward'},
                 np.eye(3), np.array(u, dtype=float), np.arange(w.shape[0]+1)/w.shape[0], w, {})


def write_record(path, mass, u):
    t = table(mass, u)
    path.mkdir()
    nt = t.mass.shape[1]; np_ = t.mass.shape[0]
    q = t.mass.sum(axis=1)
    a = np.r_[0., q.cumsum()]; a[-1] = 1
    c = np.zeros((np_,nt+1))
    for j in range(np_):
        c[j] = np.r_[0., t.mass[j].cumsum()/q[j]] if q[j] else t.u
    c[:,-1] = 1
    np.save(path/'phi_cdf.npy',a);np.save(path/'theta_given_phi_cdf.npy',c);np.save(path/'u_edges.npy',t.u)
    reconstructed = np.diff(a)[:,None]*np.diff(c,axis=1)
    g = np.sum(reconstructed*(1-t.u[:-1]-t.u[1:]))/reconstructed.sum()
    m = dict(t.metadata, schema='rainbow.phase_cdf.numpy.v2',complete=True,dtype='<f8',order='C',
             theta_count=nt,phi_count=np_,density_measure='solid_angle_sr',
             cell_model='constant_density_per_spherical_cell',coordinates='u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)',
             cdf_axis_order=['phi_cell','theta_edge'],normalization='integral_p_domega_equals_one',
             input_polarization='unpolarized',sampling_frame_columns=np.eye(3).tolist(),
             orientation={'rows':3,'csv_sha256':'test'},
             hg={'g':float(g),'method':'first_moment_of_saved_cell_pdf','target':'saved_cdf',
                 'cosine_convention':'dot(incident_propagation,outgoing_propagation)'})
    (path/'metadata.json').write_text(json.dumps(m),encoding='utf-8')


class AnalysisTests(unittest.TestCase):
    def test_uniform_measure(self):
        u = np.array([0,0.03,0.3,0.7,1.])
        t = table(np.tile(np.diff(u)/3,(3,1)),u)
        np.testing.assert_allclose(t.density_omega,1/(4*np.pi),rtol=2e-15)
        self.assertAlmostEqual(t.g,0,places=15)
        for angle in [0,1,20,90,173,180]:
            self.assertAlmostEqual(cone_mass(t,angle),math.sin(math.radians(angle)/2)**2,places=14)

    def test_forward_cap(self):
        t = table([[0.5,0],[0.5,0]],[0,0.1,1])
        self.assertAlmostEqual(t.g,0.9,places=15)
        self.assertAlmostEqual(cone_mass(t,math.degrees(2*math.asin(math.sqrt(.05)))) ,.5,places=14)

    def test_identity(self):
        t = table([[.1,.2],[.3,.4]],[0,.3,1])
        s = compare_tables(t,t)
        for key in ['tv','l1','l2_per_sqrt_sr','hellinger','jensen_shannon_nats','theta_marginal_tv']:
            self.assertEqual(s[key],0,key)

    def test_disjoint_support(self):
        a = table([[1,0]],[0,.5,1]);b = table([[0,1]],[0,.5,1])
        s = compare_tables(a,b)
        self.assertAlmostEqual(s['tv'],1);self.assertAlmostEqual(s['hellinger'],1)
        self.assertAlmostEqual(s['jensen_shannon_nats'],math.log(2))
        self.assertEqual(s['mass_a_where_b_is_zero'],1)

    def test_isolated_cell_coarsening(self):
        a = table([[1,0],[0,0]],[0,.5,1]);b = table([[1]],[0,1])
        self.assertAlmostEqual(compare_tables(a,b)['tv'],.75)

    def test_common_refinement_preserves_information(self):
        a = table([[.5,0,.5,0]],[0,.25,.5,.75,1])
        b = table([[0,.5,0,.5]],[0,.25,.5,.75,1])
        # Coarsening to two bins would incorrectly report zero; refinement does not.
        self.assertAlmostEqual(compare_tables(a,b)['tv'],1)

    def test_non_nested_unequal_phi(self):
        ua = np.array([0,.1,.5,1]);ub = np.array([0,.2,.6,.8,1])
        a=table(np.tile(np.diff(ua)/3,(3,1)),ua)
        b=table(np.tile(np.diff(ub)/7,(7,1)),ub)
        s=compare_tables(a,b)
        self.assertLess(s['tv'],1e-15)
        self.assertAlmostEqual(s['probability_mass_a_integrated'],1)
        self.assertAlmostEqual(s['probability_mass_b_integrated'],1)

    def test_one_ulp_edge_not_merged(self):
        x=np.nextafter(.5,1.)
        a=table([[1,0]],[0,.5,1]);b=table([[1,0]],[0,x,1])
        s=compare_tables(a,b)
        self.assertEqual(s['common_refinement_shape_phi_theta'],[1,3])
        self.assertGreater(s['tv'],0)
        self.assertAlmostEqual(s['tv']/(x-.5),2,places=12)

    def test_js_tiny_difference_nonnegative(self):
        a=table([[.5,.5]],[0,.5,1]);b=table([[.5+1e-10,.5-1e-10]],[0,.5,1])
        s=compare_tables(a,b)
        self.assertGreater(s['jensen_shannon_nats'],0)
        self.assertLess(s['jensen_shannon_nats'],1e-18)

    def test_frame_refusal(self):
        a=table([[.5,.5]],[0,.5,1]);b=table([[.5,.5]],[0,.5,1]);b.frame=-np.eye(3)
        with self.assertRaises(ValueError):compare_tables(a,b)

    def test_refinement_limit_refusal(self):
        a=table([[.5,.5]],[0,.5,1])
        with self.assertRaises(ValueError):compare_tables(a,a,max_refined_cells=1)

    def test_reader_round_trip_and_no_mutation(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record';write_record(p,[[.25,.75],[0,0]],[0,.2,1])
            before={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in p.iterdir()}
            t=load_table(p);s=angular_summary(t)
            self.assertEqual(s['reader_validation'],'passed')
            self.assertEqual(s['computed_from_saved_arrays']['zero_cells'],2)
            self.assertEqual(before,{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in p.iterdir()})

    def test_reject_bad_cdf(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record';write_record(p,[[.25,.75],[0,0]],[0,.2,1])
            x=np.load(p/'theta_given_phi_cdf.npy');x[0,1]=-1;np.save(p/'theta_given_phi_cdf.npy',x)
            with self.assertRaises(ValueError):load_table(p)

    def test_memory_limit_before_dense_load(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record';write_record(p,[[.25,.75],[0,0]],[0,.2,1])
            with self.assertRaises(ValueError):load_table(p,max_cells=1)

    def test_cli_existing_output_refused(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record';write_record(p,[[.25,.75],[0,0]],[0,.2,1])
            out=Path(root)/'report'
            self.assertEqual(main(['inspect',str(p),'--out',str(out),'--no-plots']),0)
            self.assertEqual(main(['inspect',str(p),'--out',str(out),'--no-plots']),2)

    def test_cli_changed_plan_requires_acknowledgement(self):
        with tempfile.TemporaryDirectory() as root:
            a=Path(root)/'a';b=Path(root)/'b'
            for p in [a,b]:write_record(p,[[.25,.75],[0,0]],[0,.2,1])
            m=json.loads((b/'metadata.json').read_text());m['orientation']['csv_sha256']='changed'
            (b/'metadata.json').write_text(json.dumps(m))
            args=['compare',str(a),str(b),'--out',str(Path(root)/'diff'),'--no-plots']
            self.assertEqual(main(args),2)
            self.assertEqual(main([*args,'--allow-orientation-plan-difference']),0)

    def test_part_refusal(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record.part';write_record(p,[[.5,.5]],[0,.5,1])
            with self.assertRaises(ValueError):load_table(p)

    def test_corrupt_hg_label_refusal(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record';write_record(p,[[.5,.5]],[0,.5,1])
            m=json.loads((p/'metadata.json').read_text());m['hg']['g']=.5
            (p/'metadata.json').write_text(json.dumps(m))
            with self.assertRaises(ValueError):load_table(p)

    def test_nonfinite_metadata_refusal(self):
        with tempfile.TemporaryDirectory() as root:
            p=Path(root)/'record';write_record(p,[[.5,.5]],[0,.5,1])
            m=json.loads((p/'metadata.json').read_text());m['wavelength_nm']=float('inf')
            (p/'metadata.json').write_text(json.dumps(m))
            with self.assertRaises(ValueError):load_table(p)


if __name__=='__main__':unittest.main(verbosity=2)
