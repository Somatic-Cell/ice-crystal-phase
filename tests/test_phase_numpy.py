from __future__ import annotations
import json
import math
from pathlib import Path
import tempfile
import unittest
import numpy as np
from rainbow_dataset import PhaseRecord

try:
    import torch
except ImportError:
    torch = None


def create_fixture(root: Path, nt: int = 17, np_: int = 13, *, zero_columns: bool = True,
                   uniform: bool = False) -> tuple[np.ndarray, np.ndarray]:
    root.mkdir(parents=True)
    # Independent test fixture, NOT a simulated/validated optical dataset.
    edges = np.sin(np.linspace(0, np.pi/2, nt+1, dtype=np.float64))**2
    edges[0], edges[-1] = 0, 1
    i, j = np.meshgrid(np.arange(nt), np.arange(np_), indexing="ij")
    density = np.ones((nt, np_), dtype=np.float64) if uniform else (
        0.3 + (i+1)**2 + 0.4 * np.sin(2*np.pi*(j+0.5)/np_)**2)
    if not uniform:
        density[(i+2*j) % 5 == 0] = 0
    if zero_columns and np_ > 2:
        density[:, 0] = 0
        density[:, -1] = 0
    weights = density * np.diff(edges)[:, None]
    columns = weights.sum(axis=0)
    marginal = np.concatenate(([0.0], np.cumsum(columns)))
    marginal /= marginal[-1]
    marginal[-1] = 1
    cond = np.empty((np_, nt+1), dtype=np.float64)
    for col in range(np_):
        cond[col] = np.concatenate(([0.0], np.cumsum(weights[:, col]))) / columns[col] if columns[col] else edges
    cond[:, 0], cond[:, -1] = 0, 1
    np.save(root/"phi_cdf.npy", marginal, allow_pickle=False)
    np.save(root/"theta_given_phi_cdf.npy", cond, allow_pickle=False)
    np.save(root/"u_edges.npy", edges, allow_pickle=False)
    metadata = {
        "schema": "rainbow.phase_cdf.numpy.v1", "complete": True,
        "dtype": "<f8", "order": "C", "theta_count": nt, "phi_count": np_,
        "density_measure": "solid_angle_sr", "cell_model": "constant_density_per_spherical_cell",
        "coordinates": "u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)",
        "normalization": "integral_p_domega_equals_one", "input_polarization": "unpolarized",
        "cdf_axis_order": ["phi_cell", "theta_edge"],
        "sampling_frame_columns": [[1,0,0],[0,1,0],[0,0,1]],
        "source": "synthetic_test_fixture_not_optical_ground_truth",
    }
    (root/"metadata.json").write_text(json.dumps(metadata), encoding="utf-8")
    return marginal, cond


class NumpyPhaseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)/"record"
        create_fixture(self.root)
        self.record = PhaseRecord(self.root)

    def tearDown(self):
        # Windows refuses to delete mapped files until the mapping is closed.
        self.record.close()
        self.temp.cleanup()

    def test_memmap_and_bit_exact(self):
        self.assertIsInstance(self.record.theta_cdf, np.memmap)
        for name, a in [("phi_cdf.npy",self.record.phi_cdf),
                        ("theta_given_phi_cdf.npy",self.record.theta_cdf)]:
            expected = np.load(self.root/name,allow_pickle=False)
            np.testing.assert_array_equal(a.view(np.uint64),expected.view(np.uint64))

    def test_sample_pdf_and_zeros(self):
        b = self.record.sample(100_000, rng=np.random.default_rng(9))
        self.assertTrue(np.all(b.pdf_omega > 0))
        self.assertTrue(np.all((b.uv >= 0) & (b.uv < 1)))
        self.assertTrue(np.all((b.cells[:,1] > 0) & (b.cells[:,1] < self.record.np-1)))
        np.testing.assert_allclose(b.pdf_omega, self.record.pdf_omega(b.uv),rtol=5e-15,atol=0)
        np.testing.assert_allclose(np.linalg.norm(self.record.directions(b.uv),axis=1),1,rtol=0,atol=5e-16)

    def test_distribution_frequencies(self):
        b = self.record.sample(300_000,rng=np.random.default_rng(22))
        counts = np.bincount(b.cells[:,1]*self.record.nt+b.cells[:,0],
                             minlength=self.record.nt*self.record.np).reshape(self.record.np,self.record.nt)
        prob = np.diff(self.record.phi_cdf)[:,None]*np.diff(self.record.theta_cdf,axis=1)
        expected = 300_000*prob
        self.assertTrue(np.all(np.abs(counts-expected) < 7*np.sqrt(expected)+10))
        self.assertEqual(int(counts[prob==0].sum()),0)

    def test_edge_uniforms(self):
        # Every marginal edge, including plateaus and exact zero, and the final representable value below one.
        xs = np.unique(np.r_[self.record.phi_cdf[:-1],0,np.nextafter(1.,0.)])
        xs = xs[xs<1]
        xi = np.array([(a,b) for a in xs for b in (0.,np.nextafter(1.,0.),.25,.5)],dtype=np.float64)
        b = self.record.sample_uniforms(xi)
        np.testing.assert_allclose(b.pdf_omega,self.record.pdf_omega(b.uv),rtol=5e-15,atol=0)

    def test_integral(self):
        u = (self.record.u_edges[1:]+self.record.u_edges[:-1])/2
        v = (np.arange(self.record.np)+.5)/self.record.np
        uu,vv=np.meshgrid(u,v,indexing="ij")
        density=self.record.pdf_omega(np.stack((uu.ravel(),vv.ravel()),-1)).reshape(uu.shape)
        integral=(density*np.diff(self.record.u_edges)[:,None]*(4*np.pi/self.record.np)).sum()
        self.assertAlmostEqual(float(integral),1,places=14)

    def test_bad_uniform(self):
        for x in [np.array([[1.,0.]]),np.array([[0.,np.nan]]),np.array([[-1.,0.]])]:
            with self.assertRaises(ValueError):self.record.sample_uniforms(x)
        with self.assertRaises(TypeError):self.record.sample_uniforms(np.zeros((2,2),dtype=np.float32))

    def test_bad_stage_schema(self):
        m=json.loads((self.root/"metadata.json").read_text());m["schema"]="unknown"
        (self.root/"metadata.json").write_text(json.dumps(m))
        with self.assertRaises(ValueError):PhaseRecord(self.root)

    def test_partial(self):
        with self.assertRaises(ValueError):PhaseRecord(self.root.with_name("record.part"))

    def test_empty_batch(self):
        self.assertEqual(self.record.sample(0).uv.shape,(0,2))

    def test_batches(self):
        batches=list(self.record.iter_batches(batch_size=20,batches=3,seed=2))
        self.assertEqual(len(batches),3)
        self.assertFalse(np.array_equal(batches[0].uv,batches[1].uv))

    def test_conditional_validation_chunk_seams(self):
        self.record.validate(block_elements=7)

    def test_external_pdf_outside_domain(self):
        p=self.record.pdf_omega(np.array([[-.1,.3],[.1,1.0],[1.0,.2]],dtype=np.float64))
        np.testing.assert_array_equal(p,0)

    def test_close_is_idempotent(self):
        self.record.close()
        self.record.close()
        with self.assertRaises(RuntimeError):
            self.record.sample(1)

    def test_reject_corrupt_arrays(self):
        # Separate records; never overwrite a memory-mapped live file.
        for case in ("dtype", "decreasing", "truncated", "trailing", "nonfinite"):
            with self.subTest(case=case):
                root=Path(self.temp.name)/case
                create_fixture(root)
                path=root/"theta_given_phi_cdf.npy"
                data=np.load(path,allow_pickle=False)
                if case=="dtype":
                    np.save(path,data.astype(np.float32),allow_pickle=False)
                elif case=="decreasing":
                    data[1,5]=.9;data[1,6]=.1
                    np.save(path,data,allow_pickle=False)
                elif case=="nonfinite":
                    data[1,5]=np.nan
                    np.save(path,data,allow_pickle=False)
                elif case=="truncated":
                    with path.open("r+b") as f:f.truncate(path.stat().st_size-8)
                else:
                    with path.open("ab") as f:f.write(b"extra")
                with self.assertRaises(ValueError):PhaseRecord(root)
                # Constructor-failure cleanup must permit deletion on Windows.
                path.unlink()

    @unittest.skipIf(torch is None,"PyTorch not installed")
    def test_torch_matches_numpy(self):
        s=self.record.to_torch("cpu")
        xi=np.random.default_rng(123).random((12345,2))
        a=self.record.sample_uniforms(xi)
        b=s.sample_uniforms(torch.from_numpy(xi))
        np.testing.assert_array_equal(a.cells,b.cells.numpy())
        np.testing.assert_allclose(a.uv,b.uv.numpy(),rtol=0,atol=2e-16)
        np.testing.assert_allclose(a.pdf_omega,b.pdf_omega.numpy(),rtol=5e-15,atol=0)
        np.testing.assert_allclose(s.pdf_omega(b.uv).numpy(),a.pdf_omega,rtol=5e-15,atol=0)
        self.assertEqual(s.theta_flat.dtype,torch.float64)
        self.assertEqual(s.theta_flat.ndim,1)

    @unittest.skipIf(torch is None,"PyTorch not installed")
    def test_torch_generation(self):
        s=self.record.to_torch("cpu")
        b=s.sample(1000,generator=torch.Generator().manual_seed(3))
        self.assertFalse(b.uv.requires_grad)
        self.assertTrue(bool(torch.all(b.pdf_omega>0)))
        np.testing.assert_allclose(s.pdf_omega(b.uv).numpy(),b.pdf_omega.numpy(),rtol=5e-15,atol=0)
        self.assertEqual(s.sample(0).uv.shape,(0,2))

    @unittest.skipIf(torch is None or not torch.cuda.is_available(),"No CUDA device for PyTorch")
    def test_torch_cuda(self):
        s=self.record.to_torch("cuda")
        xi=np.random.default_rng(123).random((12345,2))
        b=s.sample_uniforms(torch.tensor(xi,device="cuda"))
        a=self.record.sample_uniforms(xi)
        np.testing.assert_array_equal(a.cells,b.cells.cpu().numpy())
        np.testing.assert_allclose(a.uv,b.uv.cpu().numpy(),atol=2e-15,rtol=0)
        np.testing.assert_allclose(a.pdf_omega,b.pdf_omega.cpu().numpy(),rtol=1e-13,atol=0)


class OtherShapes(unittest.TestCase):
    def test_uniform_and_tail_sizes(self):
        for nt,np_ in [(1,1),(3,4),(4,3),(255,3),(256,7),(257,5),(513,11)]:
            with self.subTest(nt=nt,np=np_):
                with tempfile.TemporaryDirectory() as t:
                    path=Path(t)/"record";create_fixture(path,nt,np_,zero_columns=False,uniform=True)
                    r=PhaseRecord(path)
                    b=r.sample(1000,rng=np.random.default_rng(1))
                    np.testing.assert_allclose(b.pdf_omega,1/(4*np.pi),rtol=5e-10,atol=0)
                    np.testing.assert_allclose(b.pdf_omega,r.pdf_omega(b.uv),rtol=5e-15,atol=0)
                    r.close()
    def test_exact_phi_boundaries(self):
        for np_ in (180,720,14400):
            with self.subTest(np=np_):
                with tempfile.TemporaryDirectory() as t:
                    path=Path(t)/"record"
                    create_fixture(path,2,np_,zero_columns=False,uniform=True)
                    with PhaseRecord(path) as r:
                        xi=np.column_stack((np.asarray(r.phi_cdf[:-1]),np.full(np_,.37)))
                        b=r.sample_uniforms(xi)
                        np.testing.assert_allclose(b.pdf_omega,r.pdf_omega(b.uv),rtol=1e-14,atol=0)
                        if torch is not None:
                            ts=r.to_torch("cpu")
                            tb=ts.sample_uniforms(torch.from_numpy(xi))
                            np.testing.assert_allclose(tb.pdf_omega.numpy(),ts.pdf_omega(tb.uv).numpy(),rtol=1e-14,atol=0)


if __name__=="__main__":unittest.main()
