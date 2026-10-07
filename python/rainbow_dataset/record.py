from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator
import json
import math
import numpy as np
from numpy.typing import NDArray

SCHEMA = "rainbow.phase_cdf.numpy.v1"


@dataclass(frozen=True)
class SampleBatch:
    """uv[:,0]=(1-cos(theta))/2, uv[:,1]=(phi+pi)/(2*pi)."""
    uv: NDArray[np.float64]
    pdf_omega: NDArray[np.float64]
    cells: NDArray[np.int64]  # [theta_cell, phi_cell]

    @property
    def pdf_uv(self) -> NDArray[np.float64]:
        return self.pdf_omega * (4.0 * np.pi)


class PhaseRecord:
    """Read standard NPYs; never rebuild, renormalize, clamp or transpose the CDF.

    validate=True performs a bounded-memory integrity scan once, not per batch.
    mmap avoids eagerly allocating all arrays, but random disk access is NOT
    equivalent to GPU memory access. Use to_torch() once per resident condition.
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
        if m.get("schema") != SCHEMA or m.get("complete") is not True:
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
        self._validated = True

    @staticmethod
    def _coordinates(values: NDArray[np.float64]) -> NDArray[np.float64]:
        a = np.asarray(values)
        if a.dtype != np.float64 or a.ndim != 2 or a.shape[1] != 2:
            raise TypeError("Expected a float64 array with shape (batch, 2)")
        if not np.isfinite(a).all():
            raise ValueError("Coordinates must be finite")
        return a

    def sample_uniforms(self, xi: NDArray[np.float64]) -> SampleBatch:
        self._require_open()
        xi = self._coordinates(xi)
        if np.any(xi < 0) or np.any(xi >= 1):
            raise ValueError("Uniform coordinates must lie in [0, 1)")
        count = len(xi)
        j = np.searchsorted(self.phi_cdf, xi[:, 0], side="right") - 1
        base = j * (self.nt + 1)
        lo = np.zeros(count, dtype=np.int64)
        hi = np.full(count, self.nt, dtype=np.int64)
        # Element-wise binary search: O(batch) working memory, not O(batch*Ntheta).
        for _ in range(self.nt.bit_length()):
            mid = (lo + hi) // 2
            right = self._theta_flat[base + mid] <= xi[:, 1]
            lo = np.where(right, mid + 1, lo)
            hi = np.where(right, hi, mid)
        i = lo - 1
        a0 = self.phi_cdf[j]
        a = self.phi_cdf[j + 1] - a0
        b0 = self._theta_flat[base + i]
        b = self._theta_flat[base + i + 1] - b0
        if np.any(a <= 0) or np.any(b <= 0):
            raise ValueError("CDF inversion selected an empty interval")
        r0 = (xi[:, 0] - a0) / a
        r1 = (xi[:, 1] - b0) / b
        lower = self.u_edges[i]
        upper = self.u_edges[i + 1]
        du = upper - lower
        u = np.minimum(lower + r1 * du, np.nextafter(upper, lower))
        vlo, vhi = j / self.np, (j + 1) / self.np
        v = np.minimum((j + r0) / self.np, np.nextafter(vhi, vlo))
        pdf = a * b * self.np / (4.0 * np.pi * du)
        return SampleBatch(np.stack((u, v), axis=1), pdf, np.stack((i, j), axis=1))

    def sample(self, count: int, *, rng: np.random.Generator | None = None) -> SampleBatch:
        if type(count) is not int or count < 0:
            raise ValueError("count must be nonnegative")
        rng = np.random.default_rng() if rng is None else rng
        return self.sample_uniforms(rng.random((count, 2), dtype=np.float64))

    def iter_batches(self, *, batch_size: int, batches: int | None = None,
                     seed: int | None = None) -> Iterator[SampleBatch]:
        if type(batch_size) is not int or batch_size < 1:
            raise ValueError("batch_size must be positive")
        if batches is not None and (type(batches) is not int or batches < 0):
            raise ValueError("batches must be nonnegative or None")
        rng = np.random.default_rng(seed)
        index = 0
        while batches is None or index < batches:
            yield self.sample(batch_size, rng=rng)
            index += 1

    def pdf_omega(self, uv: NDArray[np.float64]) -> NDArray[np.float64]:
        self._require_open()
        uv = self._coordinates(uv)
        result = np.zeros(len(uv), dtype=np.float64)
        valid = np.all((uv >= 0) & (uv < 1), axis=1)
        x = uv[valid]
        i = np.searchsorted(self.u_edges, x[:, 0], side="right") - 1
        j = np.minimum(np.floor(x[:, 1] * self.np).astype(np.int64), self.np - 1)
        # Correct quotient rounding at representable phi-cell boundaries.
        j -= x[:, 1] < j / self.np
        j += x[:, 1] >= (j + 1) / self.np
        base = j * (self.nt + 1) + i
        mass = (self.phi_cdf[j + 1] - self.phi_cdf[j]) * (
            self._theta_flat[base + 1] - self._theta_flat[base])
        result[valid] = mass * self.np / (4.0 * np.pi * (self.u_edges[i + 1] - self.u_edges[i]))
        return result

    def pdf_uv(self, uv: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.pdf_omega(uv) * (4.0 * np.pi)

    def directions(self, uv: NDArray[np.float64]) -> NDArray[np.float64]:
        self._require_open()
        uv = self._coordinates(uv)
        if np.any(uv < 0) or np.any(uv >= 1):
            raise ValueError("uv must lie in [0, 1)")
        u, v = uv[:, 0], uv[:, 1]
        mu = 1.0 - 2.0 * u
        radial = 2.0 * np.sqrt(u * (1.0 - u))
        phi = (2.0 * np.pi) * v - np.pi
        local = np.stack((radial * np.cos(phi), radial * np.sin(phi), mu), axis=1)
        return local @ self.frame.T

    def to_torch(self, device: str = "cpu"):
        self._require_open()
        if not self._validated:
            self.validate()
        from .torch_sampler import TorchPhaseSampler
        return TorchPhaseSampler(self, device=device)
