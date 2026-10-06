from __future__ import annotations
import copy
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import struct
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("patch_failure_analysis",ROOT/"tools/analyze_patch_failures.py")
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

def bits(x):
    return struct.unpack(">I",struct.pack(">f",x))[0]

def custom_patch(points):
    # Only fields required by the exact sign audit, not a capture report.
    return {"vertices":[{"direction_bits":[bits(v) for v in p]} for p in points],
            "host_audit":{"hemisphere_axis_bits":[bits(sum(p[k] for p in points)) for k in range(3)]}}

class ExactAuditTests(unittest.TestCase):
    def setUp(self):
        self.report=m.load(ROOT/"tests/data/patch_failure_fixture.json")
    def test_regular_positive(self):
        a=m.exact_patch(self.report['patches'][0])
        self.assertEqual(a['classification'],'regular_positive')
        self.assertEqual(a['signs'],[1,1,1,1])
    def test_regular_negative(self):
        self.assertEqual(m.exact_patch(self.report['patches'][3])['classification'],'regular_negative')
    def test_real_sign_change(self):
        self.assertEqual(m.exact_patch(self.report['patches'][1])['classification'],'jacobian_sign_change')
    def test_exact_zero(self):
        self.assertEqual(m.exact_patch(self.report['patches'][2])['classification'],'jacobian_zero_at_corner')
    def test_small_nonzero_not_clamped(self):
        e=2.0**-120
        p=custom_patch([(1,0,0),(1,e,0),(1,0,e),(1,e,e)])
        a=m.exact_patch(p)
        self.assertEqual(a['signs'],[1,1,1,1])
        self.assertEqual(Fraction(a['jacobians'][0]),Fraction(1,2**240))
    def test_subnormal_bits(self):
        self.assertEqual(m.rational(1),Fraction(1,2**149))
    def test_negative_zero(self):
        self.assertEqual(m.rational(0x80000000),0)
    def test_nan_corner(self):
        p=self.report['patches'][0]
        p['vertices'][0]['direction_bits'][0]=0x7FC00000
        self.assertEqual(m.exact_patch(p)['classification'],'nonfinite_corner')
    def test_invalid_bit_pattern_rejected(self):
        for x in (-1,2**32,True,1.5):
            with self.subTest(x=x), self.assertRaises(ValueError): m.f32(x)
    def test_hemisphere_candidate_failure_is_not_nonexistence(self):
        p=self.report['patches'][0]
        p['host_audit']['hemisphere_axis_bits']=[0,0,0]
        a=m.exact_patch(p)
        self.assertEqual(a['classification'],'original_hemisphere_candidate_not_certified')
        self.assertTrue(a['alternate_exact_sum_hemisphere_certified'])
    def test_enclosures_match_exact(self):
        for p in self.report['patches']:
            self.assertEqual(m.validate_enclosures(m.exact_patch(p),p),[])
    def test_bad_interval_detected(self):
        p=self.report['patches'][0]
        p['host_audit']['jacobian64'][0]=[0.0,0.0]
        self.assertTrue(m.validate_enclosures(m.exact_patch(p),p))
    def test_actual_thin_fp32_inconclusive_fixture(self):
        p=self.report['patches'][4]
        self.assertEqual(p['host_audit']['orientation32'],0)
        self.assertEqual(p['host_audit']['interval_orientation64'],1)
        self.assertEqual(m.exact_patch(p)['classification'],'regular_positive')
        self.assertEqual(m.classify_reason(p,m.exact_patch(p)),'fp32_certificate_inconclusive_fp64_certifies')
    def test_precision_reason(self):
        p=self.report['patches'][0];p['stored_status']=3;p['host_audit']['orientation32']=0
        self.assertEqual(m.classify_reason(p,m.exact_patch(p)),'fp32_certificate_inconclusive_fp64_certifies')
    def test_metric_reason(self):
        p=self.report['patches'][0];p['host_audit']['signed_omega32_bits']=0
        self.assertEqual(m.classify_reason(p,m.exact_patch(p)),'signed_solid_angle_metric_issue')
    def test_does_not_mutate_input(self):
        before=copy.deepcopy(self.report)
        m.analyze(self.report)
        self.assertEqual(self.report,before)
    def test_report_roundtrip(self):
        summary,p,d,h=m.analyze(self.report)
        self.assertEqual(len(p),5);self.assertEqual(len(d),1);self.assertEqual(len(h),1)
        self.assertEqual(summary['integrity_issues'],[])
        self.assertEqual(summary['first_unavailable_stage_counts'],{'optical':1})
        self.assertEqual(h[0]['exact_jacobian_sign_at_stored_uv'],'-1')
    def test_all_recorded_hits_required(self):
        self.report['directions'][0]['hit_count']=2
        with self.assertRaises(ValueError):m.analyze(self.report)
    def test_hit_identity(self):
        self.report['directions'][0]['hits'][0]['patch_id']=123
        with self.assertRaises(ValueError):m.analyze(self.report)
    def test_duplicate_patch(self):
        self.report['patches'].append(self.report['patches'][0])
        with self.assertRaises(ValueError):m.analyze(self.report)
    def test_duplicate_direction(self):
        self.report['directions'].append(self.report['directions'][0])
        with self.assertRaises(ValueError):m.analyze(self.report)
    def test_summary_mismatch(self):
        self.report['summary']['selected_patches']=12
        with self.assertRaises(ValueError):m.analyze(self.report)
    def test_first_stage_focal(self):
        d=self.report['directions'][0];d['optical_flags']=0;d['rejected_hits']=0;d['evaluated_hits']=d['hit_count']
        self.assertEqual(m.first_stage(d),'focal')
    def test_first_stage_diffraction(self):
        d=self.report['directions'][0];d['optical_flags']=d['focal_flags']=0;d['rejected_hits']=0;d['evaluated_hits']=d['hit_count']
        self.assertEqual(m.first_stage(d),'diffraction')
    def test_underresolved_not_invalid(self):
        d=self.report['directions'][0];d['optical_flags']=d['focal_flags']=0;d['rejected_hits']=0;d['evaluated_hits']=d['hit_count'];d['diffraction_flags']=3
        self.assertEqual(m.first_stage(d),'nonfinite_without_explanatory_flag')
    def test_cli(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'audit'
            self.assertEqual(m.main([str(ROOT/'tests/data/patch_failure_fixture.json'),'--out',str(out)]),0)
            self.assertTrue((out/'summary.json').is_file())
            self.assertTrue((out/'hits.csv').is_file())
            self.assertEqual(m.main([str(ROOT/'tests/data/patch_failure_fixture.json'),'--out',str(out)]),1)
    def test_incomplete_report_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'r.json';self.report['report_complete']=False;p.write_text(json.dumps(self.report))
            with self.assertRaises(ValueError):m.load(p)
    def test_json_nonstandard_nan_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'r.json';p.write_text('{"x": NaN}')
            with self.assertRaises(ValueError):m.load(p)
    def test_duplicate_json_key_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'r.json';p.write_text('{"x":0,"x":1}')
            with self.assertRaises(ValueError):m.load(p)
    def test_source_payload_missing_fails(self):
        del self.report['patches'][1]
        self.report['summary']['selected_patches']-=1
        with self.assertRaises(ValueError):m.analyze(self.report)

if __name__=='__main__':unittest.main()
