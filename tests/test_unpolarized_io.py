"""Schema/adapter tests with synthetic values, not physical Mie validation."""
from __future__ import annotations
import copy
import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vis'))
sys.path.insert(0, str(ROOT / 'vis' / 'tests'))
from optics_io import load_optics, make_view
from mie_reference import incident_projection_weights
from test_vis import fixture

spec = importlib.util.spec_from_file_location('failure_reader', ROOT/'tools/analyze_patch_failures.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def paired_optics(path):
    fixture(path)
    text = path.read_text()
    text = text.replace('rainbow_patch_optics_v1', 'rainbow_patch_optics_v2')
    text = text.replace('regular_partial_model_angular_density_NOT_phase_function', 'model_partial_angular_density_NOT_phase_function')
    text = text.replace('# input_states=one_coherent_Jones_state',
        '# input_states=unpolarized_two_orthogonal_unit_Jones_inputs\n'
        '# incident_polarization=unpolarized\n# incident_total_intensity=1\n# input_coherency=0.5,0,0,0.5')
    text = '\n'.join(line for line in text.splitlines() if not line.startswith('# incident_field='))+'\n'
    path.write_text(text)


def wave_fixture(path):
    """Full synthetic 4x8 angular grid with unit physical area factor."""
    with path.open('w', newline='') as out:
        out.write('# format=rainbow_wave_optics_v2\n# theta_count=4\n# phi_count=8\n# radius_mm=1\n'
                  '# incident_polarization=unpolarized\n# input_states=unpolarized_two_orthogonal_unit_Jones_inputs\n'
                  '# input_coherency=0.5,0,0,0.5\n# incident_total_intensity=1\n')
        names = ['direction_id','theta_index','phi_index','theta_rad','phi_rad','known_hits_complete',
                 'optical_flags','focal_valid','diffraction_valid','incoherent_total','path_total','focal_total','diffraction_total']
        writer=csv.DictWriter(out,fieldnames=names);writer.writeheader()
        for i in range(32):
            row,col=divmod(i,8)
            writer.writerow(dict(zip(names,[i,row,col,np.pi*(row+.5)/4,-np.pi+2*np.pi*(col+.5)/8,1,0,1,1,1,1,1,1])))


class UnpolarizedIOTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.path=self.root/'optics.csv';paired_optics(self.path)
    def tearDown(self):self.tmp.cleanup()
    def test_native_unpolarized_no_second_average(self):
        g=load_optics(self.path);view=make_view(g,'incoherent','total')
        self.assertTrue(g.is_unpolarized and view.is_unpolarized)
        self.assertEqual(g.incident_norm2,1)
        np.testing.assert_allclose(view.values,g.fields['regular_partial_incoherent_total']*.16)
    def test_no_fabricated_single_jones(self):
        with self.assertRaises(ValueError):_ = load_optics(self.path).jones
    def test_no_double_average(self):
        g=load_optics(self.path)
        with self.assertRaises(ValueError):make_view(g,'path','total',g)
    def test_v1_unpolarized_label_rejected(self):
        self.path.write_text(self.path.read_text().replace('rainbow_patch_optics_v2','rainbow_patch_optics_v1'))
        with self.assertRaises(ValueError):load_optics(self.path)
    def test_v2_wrong_coherency_rejected(self):
        self.path.write_text(self.path.read_text().replace('0.5,0,0,0.5','1,0,0,0'))
        with self.assertRaises(ValueError):load_optics(self.path)
    def test_v2_wrong_intensity_rejected(self):
        self.path.write_text(self.path.read_text().replace('incident_total_intensity=1','incident_total_intensity=2'))
        with self.assertRaises(ValueError):load_optics(self.path)
    def test_mie_weights_are_unpolarized(self):
        g=load_optics(self.path);view=make_view(g,'incoherent','total')
        s,p=incident_projection_weights(view)
        np.testing.assert_allclose(s,.5);np.testing.assert_allclose(p,.5)
    def test_native_wave_plot(self):
        p=self.root/'wave.csv';wave_fixture(p)
        result=subprocess.run([sys.executable,str(ROOT/'tools/plot_wave_optics.py'),str(p),'--stage','incoherent','--out',str(self.root/'plot.png')],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertTrue((self.root/'plot.png').is_file())
    def test_native_optical_plot(self):
        result=subprocess.run([sys.executable,str(ROOT/'vis/plot_optics.py'),str(self.path),'--stage','incoherent','--out',str(self.root/'o.png')],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
    def paired_report(self):
        r=json.loads((ROOT/'tests/data/patch_failure_fixture.json').read_text())
        r['format']='rainbow_patch_failure_witness_v2';r['input_polarization']='unpolarized'
        r['config']['input_coherency']=[.5,0,0,.5];r['config']['incident_total_intensity']=1
        r['config'].pop('incident_field_bits',None)
        for p in r['patches']:
            for v in p['vertices']:v['second_input_field_bits']=[0,0,1065353216,0]
        return r
    def test_report_v2(self):
        r=self.paired_report();p=self.root/'report.json';p.write_text(json.dumps(r))
        loaded=audit.load(p);summary,*_=audit.analyze(loaded)
        self.assertEqual(summary['integrity_issues'],[])
    def test_report_missing_column_rejected(self):
        r=self.paired_report();del r['patches'][0]['vertices'][0]['second_input_field_bits']
        p=self.root/'report.json';p.write_text(json.dumps(r))
        with self.assertRaises(ValueError):audit.load(p)

if __name__=='__main__':unittest.main()
