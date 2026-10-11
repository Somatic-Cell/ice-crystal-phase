"""Polar-incidence symmetry tests; no CUDA and no optical re-tracing."""
from pathlib import Path
import math
import sys
import unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'python'))
from atlas_math import (orientation_rule, coupling, coupling_for_rule, compose_mass,
                        make_cdf, frame)
from test_compose import fixture


def config(na=16):
    return {'quadrature': {'azimuth': na, 'tilt': 4, 'spin': 3}, 'components': [
        {'kind': 'horizontal_gaussian', 'sigma_deg': 10, 'fraction': .9},
        {'kind': 'isotropic', 'fraction': .1}]}


def evaluate(a, cfg, ki, analytic=True):
    r, w, meta = orientation_rule(cfg)
    if analytic:
        c, area, f, info = coupling_for_rule(a, r, w, ki, meta)
    else:
        c, area, f = coupling(a, r, w, ki)
        info = {}
    mass, tail = compose_mass(a, c)
    return mass, tail, area, f, info


class PolarAzimuthTests(unittest.TestCase):
    def test_both_poles_uniform_phi_without_changing_row_mass(self):
        a = fixture()
        a.mass *= .97
        a.residual = .03 * a.area
        for sign in (-1, 1):
            with self.subTest(sign=sign):
                old, old_tail, old_area, _, _ = evaluate(a, config(), [0, sign, 0], False)
                new, tail, area, f, info = evaluate(a, config(), [0, sign, 0])
                expected = np.repeat(old.mean(axis=1, keepdims=True), a.np, axis=1)
                np.testing.assert_allclose(new, expected, rtol=1e-12, atol=2e-17)
                np.testing.assert_allclose(new.sum(axis=1), old.sum(axis=1), rtol=1e-12, atol=2e-17)
                np.testing.assert_allclose(f, frame([0, sign, 0]), rtol=0, atol=0)
                self.assertAlmostEqual(tail, old_tail, delta=1e-16)
                self.assertAlmostEqual(area, old_area, delta=1e-15)
                self.assertAlmostEqual(new.sum() + tail, area, delta=1e-15)
                self.assertEqual(info['analytic_uniform_azimuth_components'], [0, 1])
                self.assertEqual(info['effective_pose_evaluations'], 24)
                self.assertFalse(info['output_post_filter'])

    def test_agrees_with_complete_discrete_cyclic_average(self):
        a = fixture()
        a.mass[:] = 0
        a.mass[:, :, 1, 2] = a.area
        for sign in (-1, 1):
            new, *_ = evaluate(a, config(3), [0, sign, 0])
            reference, *_ = evaluate(a, config(a.np), [0, sign, 0], False)
            np.testing.assert_allclose(new, reference, rtol=2e-12, atol=3e-17)

    def test_independent_of_requested_alpha_node_count(self):
        a = fixture()
        first = evaluate(a, config(1), [0, -1, 0])[0]
        for na in (2, 7, 16, 23):
            new, *_ = evaluate(a, config(na), [0, -1, 0])
            np.testing.assert_array_equal(new, first)

    def test_nonvertical_and_near_vertical_unchanged(self):
        a = fixture()
        cfg = config()
        for ki in ([1, 0, 0], [1, -.7, .2], [1e-13, -1, 0]):
            with self.subTest(ki=ki):
                new, tail, area, f, info = evaluate(a, cfg, ki)
                old, old_tail, old_area, old_f, _ = evaluate(a, cfg, ki, False)
                np.testing.assert_array_equal(new, old)
                np.testing.assert_array_equal(f, old_f)
                self.assertEqual(tail, old_tail)
                self.assertEqual(area, old_area)
                self.assertEqual(info['method'], 'original_fixed_orientation_quadrature')

    def test_fixed_pose_unchanged_even_at_vertical(self):
        a = fixture()
        cfg = config()
        cfg['components'] = [{'kind': 'fixed', 'fraction': 1., 'quaternion_wxyz': [1, 0, 0, 0]}]
        new, *_, info = evaluate(a, cfg, [0, -1, 0])
        old, *_ = evaluate(a, cfg, [0, -1, 0], False)
        np.testing.assert_array_equal(new, old)
        self.assertGreater(float(np.max(np.ptp(new, axis=1))), 0)
        self.assertEqual(info['analytic_uniform_azimuth_components'], [])

    def test_mixed_fixed_component_not_averaged(self):
        a = fixture()
        cfg = config()
        fixed = {'kind': 'fixed', 'fraction': .4, 'quaternion_wxyz': [1, 0, 0, 0]}
        iso = {'kind': 'isotropic', 'fraction': .6}
        cfg['components'] = [fixed, iso]
        mass, tail, area, _, info = evaluate(a, cfg, [0, -1, 0])
        pieces = []
        for comp in (fixed, iso):
            single = {'quadrature': cfg['quadrature'], 'components': [dict(comp, fraction=1.)]}
            pieces.append(evaluate(a, single, [0, -1, 0]))
        np.testing.assert_allclose(mass, .4*pieces[0][0] + .6*pieces[1][0], rtol=2e-12, atol=2e-17)
        self.assertAlmostEqual(area, .4*pieces[0][2] + .6*pieces[1][2])
        self.assertEqual(info['analytic_uniform_azimuth_components'], [1])
        self.assertGreater(float(np.max(np.ptp(mass, axis=1))), 0)

    def test_sparse_distribution_remains_sparse_in_theta_not_phi(self):
        a = fixture()
        a.mass[:] = 0
        a.mass[:, :, 1, 2] = a.area
        mass, *_ = evaluate(a, config(), [0, -1, 0])
        self.assertTrue(np.all(mass[1] > 0))
        np.testing.assert_array_equal(mass[[0, 2, 3]], 0)
        edges = np.linspace(0, 1, a.nt + 1)
        ca, cb, report = make_cdf(mass, edges)
        recovered = (np.diff(ca)[:, None]*np.diff(cb, axis=1)).T
        np.testing.assert_allclose(recovered, mass/mass.sum(), rtol=1e-13, atol=2e-16)
        self.assertLess(report['cdf_l1_mass_error'], 1e-12)
        self.assertEqual(report['cdf_lost_probability_mass'], 0.)

    def test_reject_inconsistent_declared_rule_at_pole(self):
        a = fixture()
        r, w, meta = orientation_rule(config())
        w = w[::-1].copy()
        with self.assertRaisesRegex(ValueError, 'declared'):
            coupling_for_rule(a, r, w, [0, -1, 0], meta)


if __name__ == '__main__':
    unittest.main(verbosity=2)
