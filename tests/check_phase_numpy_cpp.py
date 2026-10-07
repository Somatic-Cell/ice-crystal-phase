"""Independent C++ writer -> stock NumPy -> sampler interoperability check.

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
        witness = np.load(root / "cpp_samples.npy", allow_pickle=False)
        with PhaseRecord(root) as record:
            actual = record.sample_uniforms(witness[:, :2])
            np.testing.assert_allclose(actual.uv, witness[:, 2:4], rtol=0, atol=3e-16)
            np.testing.assert_allclose(actual.pdf_omega, witness[:, 4], rtol=2e-15, atol=0)
            np.testing.assert_allclose(actual.pdf_omega, record.pdf_omega(actual.uv), rtol=5e-15, atol=0)
        # Close all mappings before Windows TemporaryDirectory cleanup.
    print("C++ -> NPY -> NumPy: 8 bit patterns and 20000 samples/PDFs passed.")


if __name__ == "__main__":
    main()
