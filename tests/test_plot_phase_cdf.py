"""Explicitly invoked debug-tool tests. No point-cloud generation or NF tests."""
from __future__ import annotations

import contextlib
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from plot_phase_cdf import Backend, reconstruct, open_record, draw, plot_limits, main


def fixture(directory: Path, nt: int = 13, np_: int = 17, *, kind: str = "general",
            version: int = 2) -> np.ndarray:
    directory.mkdir(parents=True)
    u = np.sin(np.linspace(0, np.pi / 2, nt + 1)) ** 2
    u[0], u[-1] = 0, 1
    if kind == "nonuniform":
        u = np.linspace(0, 1, nt + 1) ** 3
    a = (4 * np.pi / np_) * np.diff(u)
    if kind == "single":
        mass = np.zeros((np_, nt), dtype=np.float64)
        mass[np_ // 3, nt // 2] = 1
    elif kind == "tiny":
        assert nt == 2 and np_ == 2
        phi = np.array([0, 1e-200, 1], dtype=np.float64)
        theta = np.array([[0, 1e-200, 1], [0, .5, 1]], dtype=np.float64)
        mass = np.diff(phi)[:, None] * np.diff(theta, axis=1)
    else:
        if kind in ("uniform", "nonuniform"):
            density = np.ones((np_, nt))
        else:
            j = np.arange(np_)[:, None]
            i = np.arange(nt)[None, :]
            density = 0.2 + (i + 1) ** 2 + 0.3 * j
            density[(i + 2 * j) % 7 == 0] = 0
            density[0] = 0  # A zero-mass phi column must remain valid.
        mass = density * a[None, :]
        mass /= mass.sum()
    if kind != "tiny":
        q = mass.sum(axis=1)
        phi = np.concatenate(([0.0], np.cumsum(q)))
        phi /= phi[-1]
        conditional = np.divide(mass, q[:, None], out=np.zeros_like(mass), where=q[:, None] > 0)
        for j in np.flatnonzero(q == 0):
            conditional[j] = np.diff(u)  # Defined but never selected.
        theta = np.concatenate((np.zeros((np_, 1)), np.cumsum(conditional, axis=1)), axis=1)
        theta /= theta[:, -1:]
        mass = np.diff(phi)[:, None] * np.diff(theta, axis=1)
    g = float((mass.astype(np.longdouble) * (1 - u[:-1] - u[1:])[None, :]).sum()
              / mass.astype(np.longdouble).sum())
    m = {
        "schema": f"rainbow.phase_cdf.numpy.v{version}", "complete": True,
        "dtype": "<f8", "order": "C", "theta_count": nt, "phi_count": np_,
        "density_measure": "solid_angle_sr", "cell_model": "constant_density_per_spherical_cell",
        "coordinates": "u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)",
        "normalization": "integral_p_domega_equals_one", "input_polarization": "unpolarized",
        "cdf_axis_order": ["phi_cell", "theta_edge"],
        "theta_definition": "angle_between_incident_and_outgoing_propagation_directions",
        "phi_range_rad": [-np.pi, np.pi], "sampling_frame_columns": np.eye(3).tolist(),
        "stage": "synthetic_validation_only", "quality": {},
    }
    if version == 2:
        m["hg"] = {"g": g, "method": "first_moment_of_saved_cell_pdf", "target": "saved_cdf",
                   "cosine_convention": "dot(incident_propagation,outgoing_propagation)"}
    (directory / "metadata.json").write_text(json.dumps(m), encoding="utf-8")
    for name, array in (("phi_cdf.npy", phi), ("theta_given_phi_cdf.npy", theta), ("u_edges.npy", u)):
        np.save(directory / name, np.asarray(array, dtype="<f8"), allow_pickle=False)
    return mass.T / a[:, None]


class PlotTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.directory = self.root / "record"
        self.backend = Backend("cpu")

    def tearDown(self) -> None:
        self.temp.cleanup()

    def run_record(self, **kwargs):
        with open_record(self.directory) as record:
            return reconstruct(record, backend=self.backend, **kwargs)

    def mutate(self, name, fn):
        a = np.load(self.directory / name)
        fn(a)
        np.save(self.directory / name, a, allow_pickle=False)

    def test_uniform_sphere(self):
        fixture(self.directory, kind="uniform")
        data = self.run_record()
        np.testing.assert_allclose(data.log10_pdf, -math.log10(4 * np.pi), rtol=0, atol=2e-13)
        self.assertAlmostEqual(data.summary["g_reconstructed"], 0, places=13)

    def test_nonuniform_u_edges(self):
        fixture(self.directory, nt=9, kind="nonuniform")
        data = self.run_record()
        np.testing.assert_allclose(data.log10_pdf, -math.log10(4 * np.pi), rtol=0, atol=3e-13)
        self.assertFalse(np.allclose(np.diff(data.theta_edges_deg), np.diff(data.theta_edges_deg)[0]))

    def test_reference_pdf_and_g(self):
        reference = fixture(self.directory)
        data = self.run_record(chunk_phi=3)
        np.testing.assert_allclose(10 ** data.log10_pdf, reference, rtol=3e-14, atol=0)
        self.assertAlmostEqual(data.summary["cdf_reconstructed_mass"], 1, places=14)
        self.assertAlmostEqual(data.summary["g_metadata"], data.summary["g_reconstructed"], places=14)

    def test_one_cell_and_zero_mask(self):
        reference = fixture(self.directory, kind="single")
        data = self.run_record()
        self.assertEqual(np.argwhere(np.isfinite(data.log10_pdf)).tolist(), [[13 // 2, 17 // 3]])
        self.assertEqual(data.summary["zero_probability_cells"], reference.size - 1)

    def test_one_by_one(self):
        fixture(self.directory, nt=1, np_=1, kind="uniform")
        data = self.run_record()
        self.assertEqual(data.log10_pdf.shape, (1, 1))
        self.assertAlmostEqual(data.log10_pdf[0, 0], -math.log10(4 * np.pi))

    def test_chunk_invariance(self):
        fixture(self.directory)
        a = self.run_record(chunk_phi=1)
        b = self.run_record(chunk_phi=7)
        np.testing.assert_array_equal(a.log10_pdf, b.log10_pdf)
        self.assertAlmostEqual(a.summary["g_reconstructed"], b.summary["g_reconstructed"], places=14)

    def test_crop_does_not_renormalize(self):
        fixture(self.directory)
        full = self.run_record()
        crop = self.run_record(theta_range=(71, 139), phi_range=(-30, 30))
        i0, i1 = crop.summary["display_theta_indices_half_open"]
        j0, j1 = crop.summary["display_phi_indices_half_open"]
        np.testing.assert_array_equal(crop.log10_pdf, full.log10_pdf[i0:i1, j0:j1])
        self.assertEqual(plot_limits(crop, None, None), plot_limits(full, None, None))

    def test_invalid_phi(self):
        fixture(self.directory)
        self.mutate("phi_cdf.npy", lambda a: a.__setitem__(2, -1))
        with self.assertRaisesRegex(ValueError, "phi CDF"):
            self.run_record()

    def test_decreasing_conditional(self):
        fixture(self.directory)
        self.mutate("theta_given_phi_cdf.npy", lambda a: a.__setitem__((4, 3), .95))
        with self.assertRaisesRegex(ValueError, "decreasing CDF"):
            self.run_record(chunk_phi=2)

    def test_nan_outside_crop_is_rejected(self):
        fixture(self.directory)
        self.mutate("theta_given_phi_cdf.npy", lambda a: a.__setitem__((0, 3), np.nan))
        with self.assertRaisesRegex(ValueError, "nonfinite data"):
            self.run_record(phi_range=(-10, 10))

    def test_invalid_conditional_endpoint(self):
        fixture(self.directory)
        self.mutate("theta_given_phi_cdf.npy", lambda a: a.__setitem__((2, -1), .9999))
        with self.assertRaisesRegex(ValueError, "last endpoint"):
            self.run_record()

    def test_invalid_edges(self):
        fixture(self.directory)
        self.mutate("u_edges.npy", lambda a: a.__setitem__(2, a[1]))
        with self.assertRaisesRegex(ValueError, "u_edges"):
            self.run_record()

    def test_hg_mismatch(self):
        fixture(self.directory)
        path = self.directory / "metadata.json"
        m = json.loads(path.read_text())
        m["hg"]["g"] += .1
        path.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError, "HG label mismatch"):
            self.run_record()

    def test_hg_wrong_convention(self):
        fixture(self.directory)
        path = self.directory / "metadata.json"
        m = json.loads(path.read_text())
        m["hg"]["target"] = "source"
        path.write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError, "convention"):
            self.run_record()

    def test_v1_without_hg(self):
        fixture(self.directory, version=1)
        self.assertIsNone(self.run_record().summary["g_metadata"])

    def test_incomplete_directory(self):
        self.directory = self.root / "record.part"
        fixture(self.directory)
        with self.assertRaises(ValueError):
            self.run_record()

    def test_no_dataset_writes(self):
        fixture(self.directory)
        before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in self.directory.iterdir()}
        self.run_record()
        after = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in self.directory.iterdir()}
        self.assertEqual(before, after)

    def test_wrong_dtype_is_not_promoted(self):
        fixture(self.directory)
        path = self.directory / "phi_cdf.npy"
        np.save(path, np.load(path).astype(np.float32))
        with self.assertRaisesRegex(ValueError, "dtype"):
            self.run_record()

    def test_positive_mass_underflow_not_masked(self):
        fixture(self.directory, nt=2, np_=2, kind="tiny")
        data = self.run_record()
        self.assertEqual(data.summary["positive_mass_product_underflows"], 1)
        self.assertEqual(data.summary["zero_probability_cells"], 0)
        self.assertTrue(math.isfinite(data.log10_pdf[0, 0]))
        self.assertAlmostEqual(data.log10_pdf[0, 0], -400 - math.log10(np.pi), places=12)

    def test_bad_range_and_chunk(self):
        fixture(self.directory)
        with self.assertRaises(ValueError):
            self.run_record(chunk_phi=0)
        with self.assertRaises(ValueError):
            self.run_record(theta_range=(80, 70))
        with self.assertRaises(ValueError):
            self.run_record(phi_range=(170, -170))

    def test_torch_arithmetic_on_cpu(self):
        try:
            import torch
        except ImportError:
            self.skipTest("torch not installed")
        fixture(self.directory)
        numpy_data = self.run_record(chunk_phi=5)
        # Exercise exactly the tensor path, but do NOT claim CUDA validation.
        self.backend.torch, self.backend.device = torch, torch.device("cpu")
        self.backend.name = "torch/cpu test only"
        torch_data = self.run_record(chunk_phi=5)
        np.testing.assert_allclose(torch_data.log10_pdf, numpy_data.log10_pdf, rtol=0, atol=2e-12)
        self.assertAlmostEqual(torch_data.summary["g_reconstructed"], numpy_data.summary["g_reconstructed"], places=13)

    def test_actual_cuda_parity(self):
        try:
            import torch
        except ImportError:
            self.skipTest("torch not installed")
        if not torch.cuda.is_available():
            self.skipTest("CUDA unavailable; actual GPU path not tested")
        fixture(self.directory)
        reference = self.run_record(chunk_phi=3)
        self.backend = Backend("cuda")
        result = self.run_record(chunk_phi=3)
        np.testing.assert_allclose(result.log10_pdf, reference.log10_pdf, rtol=0, atol=2e-12)
        self.assertAlmostEqual(result.summary["g_reconstructed"], reference.summary["g_reconstructed"], places=13)

    def test_pngs_and_summary(self):
        from PIL import Image
        fixture(self.directory, nt=9, np_=11, kind="single")
        path = self.root / "figure.png"
        raw = self.root / "cells.png"
        report = self.root / "report.json"
        with contextlib.redirect_stdout(io.StringIO()):
            code = main([str(self.directory), "--out", str(path), "--cell-image", str(raw),
                         "--report", str(report), "--device", "cpu", "--dpi", "50"])
        self.assertEqual(code, 0)
        self.assertGreater(path.stat().st_size, 1000)
        with Image.open(raw) as image:
            self.assertEqual(image.size, (11, 9))
            rgba = np.asarray(image)
            self.assertEqual(int((rgba[..., 3] > 0).sum()), 1)
            self.assertGreater(rgba[9 // 2, 11 // 3, 3], 0)
        summary = json.loads(report.read_text())
        self.assertFalse(summary["smoothed"])
        self.assertFalse(summary["renormalized"])

    def test_nonuniform_plot(self):
        fixture(self.directory, nt=9, np_=11, kind="nonuniform")
        data = self.run_record()
        draw(data, self.root / "nonuniform.png", log_limits=plot_limits(data, None, None), dpi=50)
        self.assertEqual(data.summary["render_method"], "pcolormesh_flat_no_antialias")

    def test_no_automatic_overwrite(self):
        fixture(self.directory)
        target = self.root / "exists.png"
        target.write_bytes(b"keep this")
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main([str(self.directory), "--out", str(target), "--device", "cpu"]), 1)
        self.assertEqual(target.read_bytes(), b"keep this")

    def test_refuse_outputs_inside_record(self):
        fixture(self.directory)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(main([str(self.directory), "--out", str(self.directory / "plot.png")]), 1)
        self.assertFalse((self.directory / "plot.png").exists())


if __name__ == "__main__":
    unittest.main()
