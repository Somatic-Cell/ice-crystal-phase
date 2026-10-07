"""Independent C++ writer -> stock NumPy -> CDF mass interoperability check.

All generated inputs are synthetic test fixtures, not optical ground truth.
"""
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from rainbow_dataset import PhaseRecord


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: check_phase_numpy_cpp.py CPP_TEST_EXECUTABLE")
    executable = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="rainbow_numpy_interop_") as tmp:
        root = Path(tmp) / "record"
        subprocess.run([str(executable), str(root)], check=True)
        # The stock NumPy loader reads the C++ output. No custom binary parser.
        values = np.load(root / "bits.npy", allow_pickle=False)
        expected = np.array([
            0, 0x8000000000000000, 1, 0x3FEFFFFFFFFFFFFF,
            0x3FF0000000000000, 0x7FEFFFFFFFFFFFFF,
            0x7FF0000000000000, 0x7FF8000000001234,
        ], dtype=np.uint64)
        np.testing.assert_array_equal(values.view(np.uint64), expected)
        expected_mass = np.load(root / "expected_mass.npy", allow_pickle=False)
        with PhaseRecord(root) as record:
            actual_mass = np.diff(record.phi_cdf)[:, None]*np.diff(record.theta_cdf, axis=1)
            np.testing.assert_allclose(actual_mass, expected_mass, rtol=5e-12, atol=5e-16)
            mass, g = record.mass_and_g()
            np.testing.assert_allclose(mass, 1.0, rtol=0, atol=2e-15)
            expected_g = np.sum(expected_mass*(1.0-(record.u_edges[:-1]+record.u_edges[1:])))
            np.testing.assert_allclose(g, expected_g, rtol=0, atol=2e-15)
        # Close all mappings before Windows TemporaryDirectory cleanup.
    print("C++ -> NPY -> NumPy: 8 bit patterns and all CDF cell masses/moment passed.")


if __name__ == "__main__":
    main()
