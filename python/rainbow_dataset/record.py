from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import numpy as np

SCHEMAS = {"rainbow.phase_cdf.numpy.v1", "rainbow.phase_cdf.numpy.v2"}


class PhaseRecord:
    """Read standard NPYs; never rebuild, renormalize, clamp or transpose the CDF.

    validate=True performs a bounded-memory integrity scan once.
    mmap avoids eagerly allocating all arrays, but random disk access is NOT
    equivalent to GPU memory access. Sampling and learning belong to another project.
    """

    def __init__(self, directory: str | Path, *, validate: bool = True) -> None:
        self.directory = Path(directory)
        if self.directory.name.endswith(".part"):
            raise ValueError("An incomplete .part directory is not a dataset record")
        metadata_path = self.directory / "metadata.json"
        if metadata_path.stat().st_size > 1024 * 1024:
            raise ValueError("Unexpectedly large metadata")
        with metadata_path.open("r", encoding="utf-8") as stream:
            self.metadata: dict[str, Any] = json.load(stream)
        m = self.metadata
        if m.get("schema") not in SCHEMAS or m.get("complete") is not True:
            raise ValueError("Unknown schema or incomplete record")
        for key, expected in {
            "dtype": "<f8", "order": "C", "density_measure": "solid_angle_sr",
            "cell_model": "constant_density_per_spherical_cell",
            "coordinates": "u=(1-cos(theta))/2; v=(phi+pi)/(2*pi)",
            "normalization": "integral_p_domega_equals_one",
            "input_polarization": "unpolarized",
        }.items():
            if m.get(key) != expected:
                raise ValueError(f"Unsupported {key}: {m.get(key)!r}")
        if m.get("cdf_axis_order") != ["phi_cell", "theta_edge"]:
            raise ValueError("Conditional CDF axes are not phi_cell, theta_edge")
        self.nt = self._dimension(m.get("theta_count"))
        self.np = self._dimension(m.get("phi_count"))
        self._mapped_arrays: list[np.memmap] = []
        self._validated = False
        self._closed = False
        try:
            self.phi_cdf = self._load_array("phi_cdf.npy", (self.np + 1,))
            self.theta_cdf = self._load_array("theta_given_phi_cdf.npy", (self.np, self.nt + 1))
            self.u_edges = self._load_array("u_edges.npy", (self.nt + 1,))
            self._theta_flat = self.theta_cdf.reshape(-1)  # view, no matrix replication
            self.frame = np.asarray(m["sampling_frame_columns"], dtype=np.float64)
            if (self.frame.shape != (3, 3) or not np.isfinite(self.frame).all()
                    or not np.allclose(self.frame.T @ self.frame, np.eye(3), rtol=0, atol=2e-12)
                    or not np.isclose(np.linalg.det(self.frame), 1.0, rtol=0, atol=2e-12)):
                raise ValueError("Sampling frame must be a right-handed orthonormal frame")
            if validate:
                self.validate()
        except BaseException:
            # Important on Windows: a rejected record must not keep mapped files open.
            self.close()
            raise

    def __enter__(self) -> PhaseRecord:
        self._require_open()
        return self

    def __exit__(self, *_args: Any) -> None:
        self.close()

    def close(self) -> None:
        if not self._closed:
            # Do not use array views obtained from this record after close().
            for array in self._mapped_arrays:
                array._mmap.close()
            self._closed = True

    def _require_open(self) -> None:
        if self._closed:
            raise RuntimeError("PhaseRecord has been closed")

    @staticmethod
    def _dimension(value: Any) -> int:
        if type(value) is not int or not 1 <= value <= 0x7FFFFFFF:
            raise ValueError("Invalid angular dimension")
        return value

    def _load_array(self, name: str, shape: tuple[int, ...]) -> np.memmap:
        path = self.directory / name
        a = np.load(path, mmap_mode="r", allow_pickle=False)
        if isinstance(a, np.memmap):
            self._mapped_arrays.append(a)
        if (not isinstance(a, np.memmap) or a.dtype.str != "<f8" or a.shape != shape
                or not a.flags.c_contiguous):
            raise ValueError(f"Unexpected shape/dtype/order in {name}")
        if path.stat().st_size != a.offset + a.nbytes:
            raise ValueError(f"Truncated or trailing data in {name}")
        return a

    def validate(self, block_elements: int = 1 << 20) -> None:
        self._require_open()
        if not isinstance(block_elements, int) or block_elements <= 0:
            raise ValueError("block_elements must be positive")
        if (not np.isfinite(self.u_edges).all() or self.u_edges[0] != 0
                or self.u_edges[-1] != 1 or not np.all(np.diff(self.u_edges) > 0)):
            raise ValueError("u_edges are not strictly increasing from 0 to 1")
        p = self.phi_cdf
        if (not np.isfinite(p).all() or p[0] != 0 or p[-1] != 1
                or not np.all(np.diff(p) >= 0)):
            raise ValueError("Marginal CDF is invalid")
        # Flat chunks also support conditional rows larger than a chunk. Compare
        # previous values across chunk seams, excluding only the start of a row.
        stride = self.nt + 1
        for start in range(0, self._theta_flat.size, block_elements):
            end = min(start + block_elements, self._theta_flat.size)
            a = self._theta_flat[start:end]
            if not np.isfinite(a).all() or np.any(a < 0) or np.any(a > 1):
                raise ValueError("Conditional CDF contains nonfinite or out-of-range data")
            positions = np.arange(start, end, dtype=np.int64)
            index_in_row = positions % stride
            if np.any(a[index_in_row == 0] != 0) or np.any(a[index_in_row == self.nt] != 1):
                raise ValueError("Conditional CDF endpoints must be 0 and 1")
            use = index_in_row != 0
            if np.any(a[use] < self._theta_flat[positions[use] - 1]):
                raise ValueError("Conditional CDF decreases; refusing to repair it")
        if self.metadata["schema"].endswith(".v2"):
            hg = self.metadata.get("hg", {})
            if (hg.get("method") != "first_moment_of_saved_cell_pdf" or
                    hg.get("target") != "saved_cdf" or
                    hg.get("cosine_convention") != "dot(incident_propagation,outgoing_propagation)"):
                raise ValueError("Unsupported HG moment convention")
            value = hg.get("g")
            if type(value) not in (int, float) or not np.isfinite(value) or not -1 < value < 1:
                raise ValueError("Invalid HG first moment")
            mass, g = self.mass_and_g(block_elements=block_elements)
            if abs(mass-1.0) > 1e-10 or abs(g-float(value)) > 5e-12:
                raise ValueError("HG label is inconsistent with saved CDF")
        self._validated = True

    def mass_and_g(self, *, block_elements: int = 1 << 20) -> tuple[float, float]:
        """Validation statistic only: exact cell-integral formula, no point generation."""
        self._require_open()
        if not isinstance(block_elements, int) or block_elements <= 0:
            raise ValueError("block_elements must be positive")
        mean_mu = 1.0-(self.u_edges[:-1]+self.u_edges[1:])
        mass = np.longdouble(0)
        axial = np.longdouble(0)
        # Bound memory even when a conditional row is larger than the chunk.
        for j in range(self.np):
            q = np.longdouble(self.phi_cdf[j+1]-self.phi_cdf[j])
            for i in range(0, self.nt, block_elements):
                stop = min(i+block_elements, self.nt)
                w = np.diff(self.theta_cdf[j, i:stop+1]).astype(np.longdouble)*q
                mass += w.sum(dtype=np.longdouble)
                axial += (w*mean_mu[i:stop]).sum(dtype=np.longdouble)
        if not np.isfinite(mass) or not mass > 0 or not np.isfinite(axial):
            raise ValueError("Invalid reconstructed probability mass")
        return float(mass), float(axial/mass)
