from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING
import math
import torch
if TYPE_CHECKING:
    from .record import PhaseRecord


@dataclass(frozen=True)
class TorchSampleBatch:
    uv: torch.Tensor
    pdf_omega: torch.Tensor
    cells: torch.Tensor

    @property
    def pdf_uv(self) -> torch.Tensor:
        return self.pdf_omega * (4.0 * math.pi)


class TorchPhaseSampler:
    """Float64 CDFs resident on one device. No batch-by-Ntheta CDF replication.

    Sampling is non-differentiable data generation, not a learned transform.
    Keep CDFs AND inversion random numbers in float64. Cast only final training
    coordinates explicitly when the chosen NF uses float32.
    """
    def __init__(self, record: PhaseRecord, *, device: str | torch.device = "cpu") -> None:
        self.device = torch.device(device)
        self.nt, self.np = record.nt, record.np
        self.metadata = record.metadata
        # tensor() copies read-only mmap storage; avoids sharing an unsafe
        # non-writable NumPy buffer. One transfer per condition, never per batch.
        self.phi = torch.tensor(record.phi_cdf, dtype=torch.float64, device=self.device)
        self.theta_flat = torch.tensor(record.theta_cdf, dtype=torch.float64,
                                       device=self.device).reshape(-1)
        self.edges = torch.tensor(record.u_edges, dtype=torch.float64, device=self.device)
        self.frame = torch.tensor(record.frame, dtype=torch.float64, device=self.device)

    def _check(self, x: torch.Tensor) -> None:
        if x.dtype != torch.float64 or x.device != self.phi.device or x.ndim != 2 or x.shape[1] != 2:
            raise TypeError("Expected float64 (batch,2) tensor on the sampler's device")

    @torch.no_grad()
    def sample_uniforms(self, xi: torch.Tensor, *, validate: bool = True) -> TorchSampleBatch:
        self._check(xi)
        if validate and not bool(torch.all(torch.isfinite(xi) & (xi >= 0) & (xi < 1))):
            raise ValueError("Uniform inputs must lie in [0,1)")
        j = torch.searchsorted(self.phi, xi[:, 0].contiguous(), right=True) - 1
        base = j * (self.nt + 1)
        lo = torch.zeros_like(j)
        hi = torch.full_like(j, self.nt)
        for _ in range(self.nt.bit_length()):
            mid = (lo + hi) // 2
            right = self.theta_flat[base + mid] <= xi[:, 1]
            lo = torch.where(right, mid + 1, lo)
            hi = torch.where(right, hi, mid)
        i = lo - 1
        a0 = self.phi[j]
        a = self.phi[j + 1] - a0
        b0 = self.theta_flat[base + i]
        b = self.theta_flat[base + i + 1] - b0
        r0 = (xi[:, 0] - a0) / a
        r1 = (xi[:, 1] - b0) / b
        lower, upper = self.edges[i], self.edges[i + 1]
        du = upper - lower
        u = torch.minimum(lower + r1 * du, torch.nextafter(upper, lower))
        vlo, vhi = j.to(torch.float64) / self.np, (j + 1).to(torch.float64) / self.np
        v = torch.minimum((j + r0) / self.np, torch.nextafter(vhi, vlo))
        pdf = a * b * self.np / ((4.0 * math.pi) * du)
        return TorchSampleBatch(torch.stack((u, v), -1), pdf, torch.stack((i, j), -1))

    @torch.no_grad()
    def sample(self, count: int, *, generator: torch.Generator | None = None) -> TorchSampleBatch:
        if type(count) is not int or count < 0:
            raise ValueError("count must be nonnegative")
        xi = torch.rand((count, 2), dtype=torch.float64, device=self.phi.device, generator=generator)
        # The generator guarantees this interval; no GPU-to-CPU .item() per batch.
        return self.sample_uniforms(xi, validate=False)

    @torch.no_grad()
    def pdf_omega(self, uv: torch.Tensor) -> torch.Tensor:
        self._check(uv)
        valid = torch.all(torch.isfinite(uv) & (uv >= 0) & (uv < 1), dim=-1)
        safe = torch.where(valid[:, None], uv, torch.zeros_like(uv))
        i = torch.searchsorted(self.edges, safe[:, 0].contiguous(), right=True) - 1
        j = torch.clamp(torch.floor(safe[:, 1] * self.np).to(torch.int64), max=self.np - 1)
        j = j - (safe[:, 1] < j.to(torch.float64) / self.np).to(torch.int64)
        j = j + (safe[:, 1] >= (j + 1).to(torch.float64) / self.np).to(torch.int64)
        base = j * (self.nt + 1) + i
        mass = (self.phi[j + 1] - self.phi[j]) * (self.theta_flat[base + 1] - self.theta_flat[base])
        pdf = mass * self.np / ((4.0 * math.pi) * (self.edges[i + 1] - self.edges[i]))
        return torch.where(valid, pdf, torch.zeros_like(pdf))

    @torch.no_grad()
    def pdf_uv(self, uv: torch.Tensor) -> torch.Tensor:
        return self.pdf_omega(uv) * (4.0 * math.pi)

    @torch.no_grad()
    def directions(self, uv: torch.Tensor) -> torch.Tensor:
        self._check(uv)
        u, v = uv[:, 0], uv[:, 1]
        mu = 1.0 - 2.0 * u
        radial = 2.0 * torch.sqrt(u * (1.0 - u))
        phi = (2.0 * math.pi) * v - math.pi
        local = torch.stack((radial * torch.cos(phi), radial * torch.sin(phi), mu), dim=-1)
        return local @ self.frame.T
