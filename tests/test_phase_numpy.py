"""Generation-side file/schema/CDF/moment tests; no point sampling or PyTorch."""
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rainbow_dataset import PhaseRecord


def fixture(path: Path, version: str = "v2") -> None:
    path.mkdir()
    edges = np.array([0, .1, .4, 1.], dtype=np.float64)
    phi = np.array([0, .2, 1.], dtype=np.float64)
    theta = np.array([[0, .1, .4, 1.], [0, .3, .3, 1.]], dtype=np.float64)
    g = np.sum(np.diff(phi)[:, None]*np.diff(theta, axis=1)*(1-edges[:-1]-edges[1:]))
    m = dict(schema=f"rainbow.phase_cdf.numpy.{version}", complete=True,
             dtype="<f8", order="C", theta_count=3, phi_count=2,
             density_measure="solid_angle_sr", cell_model="constant_density_per_spherical_cell",
             coordinates="u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)",
             normalization="integral_p_domega_equals_one", input_polarization="unpolarized",
             cdf_axis_order=["phi_cell", "theta_edge"], sampling_frame_columns=np.eye(3).tolist(),
             hg=dict(g=float(g), method="first_moment_of_saved_cell_pdf", target="saved_cdf",
                     cosine_convention="dot(incident_propagation,outgoing_propagation)"))
    for name, a in [("phi_cdf.npy", phi), ("theta_given_phi_cdf.npy", theta), ("u_edges.npy", edges)]:
        np.save(path/name, a, allow_pickle=False)
    (path/"metadata.json").write_text(json.dumps(m), encoding="utf8")


class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.p = Path(self.tmp.name)/"data"
        fixture(self.p)

    def tearDown(self):
        self.tmp.cleanup()

    def test_v2(self):
        with PhaseRecord(self.p) as r:
            mass, g = r.mass_and_g(block_elements=1)
            self.assertAlmostEqual(mass, 1., places=15)
            self.assertAlmostEqual(g, r.metadata["hg"]["g"], places=15)
            self.assertFalse(hasattr(r, "sample"))
            self.assertFalse(hasattr(r, "sample_uniforms"))
            self.assertFalse(hasattr(r, "to_torch"))
        with self.assertRaises(RuntimeError):
            r.mass_and_g()

    def test_v1(self):
        m=json.loads((self.p/"metadata.json").read_text());m["schema"]="rainbow.phase_cdf.numpy.v1";del m["hg"]
        (self.p/"metadata.json").write_text(json.dumps(m))
        with PhaseRecord(self.p) as r:
            self.assertAlmostEqual(r.mass_and_g()[0], 1.)

    def test_label_rejected(self):
        m=json.loads((self.p/"metadata.json").read_text());m["hg"]["g"]+=.01
        (self.p/"metadata.json").write_text(json.dumps(m))
        with self.assertRaisesRegex(ValueError,"HG label"):
            PhaseRecord(self.p)

    def test_monotonic(self):
        a=np.load(self.p/"theta_given_phi_cdf.npy");a[1,2]=.1;np.save(self.p/"theta_given_phi_cdf.npy",a)
        with self.assertRaises(ValueError):PhaseRecord(self.p)

    def test_nonfinite(self):
        a=np.load(self.p/"phi_cdf.npy");a[1]=np.nan;np.save(self.p/"phi_cdf.npy",a)
        with self.assertRaises(ValueError):PhaseRecord(self.p)

    def test_trailing(self):
        with (self.p/"u_edges.npy").open("ab") as f:f.write(b"bad")
        with self.assertRaises(ValueError):PhaseRecord(self.p)

    def test_float32(self):
        a=np.load(self.p/"u_edges.npy");np.save(self.p/"u_edges.npy",a.astype(np.float32))
        with self.assertRaises(ValueError):PhaseRecord(self.p)

    def test_partial(self):
        dst=self.p.with_name("data.part");self.p.rename(dst)
        with self.assertRaises(ValueError):PhaseRecord(dst)

    def test_bad_endpoints(self):
        a=np.load(self.p/"phi_cdf.npy");a[-1]=.9;np.save(self.p/"phi_cdf.npy",a)
        with self.assertRaises(ValueError):PhaseRecord(self.p)

    def test_bad_frame(self):
        m=json.loads((self.p/"metadata.json").read_text());m["sampling_frame_columns"][0][0]=2
        (self.p/"metadata.json").write_text(json.dumps(m))
        with self.assertRaises(ValueError):PhaseRecord(self.p)

if __name__ == "__main__":unittest.main()
