"""Deterministic finite-table orientation composition; no random-number calls.

The optical interpolation is an explicit numerical model, not an analytic
Maxwell solution: A-normalized response is bilinear in canonical incident
(theta, phi), with local spherical frames identified. Reorientation is an
exact reflection/shift of that reconstructed piecewise-constant cell density.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
from typing import Any
import numpy as np

PI = math.pi
EPS = np.finfo(np.float64).eps


def positive_int(x: Any, name: str, maximum: int = 1_000_000) -> int:
    if type(x) is not int or not 1 <= x <= maximum:
        raise ValueError(f"{name} must be an integer in [1,{maximum}]")
    return x


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def unit(v: Any) -> np.ndarray:
    a = np.asarray(v, dtype=np.float64)
    if a.shape != (3,) or not np.isfinite(a).all() or np.max(np.abs(a)) == 0:
        raise ValueError('Expected finite nonzero three-vector')
    a = a / np.max(np.abs(a))
    return a / np.linalg.norm(a)


def frame(ki: Any) -> np.ndarray:
    """Right-handed [e0 e1 ki]; match M1's least-aligned-axis convention."""
    k = unit(ki)
    ref = np.eye(3)[int(np.argmin(np.abs(k)))]
    e0 = unit(ref - k*np.dot(k, ref))
    return np.column_stack((e0, unit(np.cross(k, e0)), k))


def ry(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])


def rz(a: float) -> np.ndarray:
    c, s = math.cos(a), math.sin(a)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])


def quaternion_matrix(q: Any) -> np.ndarray:
    q = np.asarray(q, dtype=np.float64)
    if q.shape != (4,) or not np.isfinite(q).all() or abs(float(q@q)-1) > 2e-12:
        raise ValueError('Fixed quaternion must be unit wxyz')
    w, x, y, z = q/np.linalg.norm(q)
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])


def gaussian_tilt_nodes(count: int, sigma_deg: float) -> np.ndarray:
    """Quantile-midpoint quadrature of truncated N(pi/2,sigma^2), d beta.

    Each probability stratum has exactly 1/count mass (erf-CDF definition).
    This is fixed quadrature; inverse CDF is computed by deterministic bisection.
    The density is NOT multiplied by sin(beta). sigma=0 is an explicit delta.
    """
    positive_int(count, 'tilt count')
    sigma = math.radians(float(sigma_deg))
    if not math.isfinite(sigma) or sigma < 0:
        raise ValueError('sigma_deg must be finite and nonnegative')
    if sigma == 0:
        return np.full(count, PI/2)
    bound = PI/(2*math.sqrt(2)*sigma)
    e = math.erf(bound)
    if not e > 0:
        raise ValueError('Gaussian normalization is not representable')
    values = []
    for i in range(count):
        target = 2*((i+.5)/count)-1
        lo, hi = 0., PI
        for _ in range(64):
            mid = (lo+hi)/2
            c = math.erf((mid-PI/2)/(math.sqrt(2)*sigma))/e
            if c < target:
                lo = mid
            else:
                hi = mid
        values.append((lo+hi)/2)
    result = np.asarray(values)
    if count > 1 and np.any(np.diff(result) <= 0):
        raise ValueError('Positive Gaussian width is below angular resolution; use explicit sigma_deg=0 for a delta')
    return result


def orientation_rule(config: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Fixed rotations and NUMBER weights. Each component's fraction is preserved.

    Body +y axis; R=Ry(alpha) Rz(beta) Ry(gamma). alpha in [0,2pi),
    gamma in [0,pi/3) for uniform spin of a regular hexagonal prism. The latter
    is an exact quotient for uniform spin, not an assumption about Parry ice.
    """
    if set(config) != {'quadrature', 'components'}:
        raise ValueError('orientation keys must be quadrature and components')
    q = config['quadrature']
    if set(q) != {'azimuth', 'tilt', 'spin'}:
        raise ValueError('quadrature keys must be azimuth, tilt, spin')
    na, nb, ng = [positive_int(q[k], k, 4096) for k in ('azimuth', 'tilt', 'spin')]
    components = config['components']
    if not isinstance(components, list) or not components:
        raise ValueError('Nonempty orientation components required')
    fractions = [float(c['fraction']) for c in components]
    if any(not math.isfinite(f) or not 0 < f <= 1 for f in fractions) or abs(math.fsum(fractions)-1) > 1e-12:
        raise ValueError('Positive NUMBER fractions must sum to one')
    matrices, weights = [], []
    component_counts = []
    for c, fraction in zip(components, fractions):
        kind = c['kind']
        expected = {'kind', 'fraction'} | ({'sigma_deg'} if kind == 'horizontal_gaussian' else {'quaternion_wxyz'} if kind == 'fixed' else set())
        if set(c) != expected:
            raise ValueError(f'Unknown/missing component keys for {kind}')
        if kind == 'fixed':
            matrices.append(quaternion_matrix(c['quaternion_wxyz'])); weights.append(fraction)
            component_counts.append(1)
            continue
        if kind == 'isotropic':
            betas = np.arccos(1-2*(np.arange(nb)+.5)/nb)
        elif kind == 'horizontal_gaussian':
            betas = gaussian_tilt_nodes(nb, c['sigma_deg'])
        else:
            raise ValueError('Supported models: isotropic, horizontal_gaussian, fixed')
        if na*nb*ng > 2_000_000:
            raise ValueError('Too many explicit quadrature nodes')
        for ia in range(na):
            left = ry(2*PI*(ia+.5)/na)
            for beta in betas:
                tilted = left @ rz(float(beta))
                for ig in range(ng):
                    matrices.append(tilted @ ry((PI/3)*(ig+.5)/ng))
                    weights.append(fraction/(na*nb*ng))
        component_counts.append(na*nb*ng)
    r, w = np.asarray(matrices), np.asarray(weights)
    validate_rotations(r, w)
    return r, w, {'kind': 'deterministic_probability_midpoint_product_quadrature',
                  'config': config, 'component_nodes': component_counts,
                  'rotation': 'Ry(alpha)_Rz(beta)_Ry(gamma)_body_to_ensemble',
                  'horizontal_density_measure': 'd_beta; truncated normal on [0,pi]',
                  'uniform_spin_range': [0, PI/3], 'randomness': False,
                  'cdf_weights_analytic': True, 'optical_integral_analytic': False}


def validate_rotations(r: np.ndarray, w: np.ndarray) -> None:
    if r.ndim != 3 or r.shape[1:] != (3, 3) or w.shape != (len(r),) or len(r) == 0:
        raise ValueError('Malformed orientation quadrature')
    if not np.isfinite(r).all() or not np.isfinite(w).all() or np.any(w <= 0) or abs(math.fsum(w)-1) > 2e-12:
        raise ValueError('Nonfinite / nonpositive / unnormalized orientation weights')
    if not np.allclose(np.swapaxes(r, 1, 2) @ r, np.eye(3), rtol=0, atol=2e-12) or not np.allclose(np.linalg.det(r), 1, rtol=0, atol=2e-12):
        raise ValueError('Orientation nodes must be proper orthogonal rotations')


def canonical_frame(theta: float, phi: float) -> np.ndarray:
    st, ct, sp, cp = math.sin(theta), math.cos(theta), math.sin(phi), math.cos(phi)
    k = np.array([st*cp, ct, st*sp])
    e = np.array([ct*cp, -st, ct*sp])
    return np.column_stack((e, np.cross(k, e), k))


def fold_direction(direction: Any) -> tuple[float, float, np.ndarray]:
    """Return canonical t,p,Q with Q*k(t,p)=s; Q may be a reflection.

    The full D6h group has 24 O(3) operations; they are optical symmetries,
    NOT 24 allowed rotations in SO(3). The pose R remains proper.
    """
    k = unit(direction)
    theta = math.atan2(math.hypot(k[0], k[2]), abs(k[1]))
    az = math.atan2(k[2], k[0]) % (2*PI)
    sector = min(5, int(az/(PI/3)))
    local = az-sector*(PI/3)
    mirror = local > PI/6
    phi = PI/3-local if mirror else local
    # Last ulps only. No finite angular folding/blur tolerance.
    if phi < -64*EPS or phi > PI/6+64*EPS:
        raise ValueError('Canonical azimuth outside fundamental region')
    phi = min(PI/6, max(0., phi))
    angle = (sector+int(mirror))*(PI/3)
    q = ry(-angle) @ np.diag([1., -1. if k[1] < 0 else 1., -1. if mirror else 1.])
    if np.max(np.abs(q @ canonical_frame(theta, phi)[:, 2]-k)) > 2e-12:
        raise ValueError('Symmetry unfolding is inconsistent')
    return theta, phi, q


def interpolation(t: float, p: float, ntheta: int, nphi: int) -> list[tuple[int, int, float]]:
    # Endpoint-inclusive atlas: no nearest-neighbour substitute or extrapolation.
    if not math.isfinite(t) or not math.isfinite(p) or not -64*EPS <= t <= PI/2+64*EPS or not -64*EPS <= p <= PI/6+64*EPS:
        raise ValueError('Incident interpolation outside canonical domain')
    x, y = t/(PI/2)*ntheta, p/(PI/6)*nphi
    x, y = min(float(ntheta), max(0., x)), min(float(nphi), max(0., y))
    i, j = min(ntheta-1, int(x)), min(nphi-1, int(y))
    a, b = x-i, y-j
    return [(i+di, j+dj, u*v) for di, u in ((0, 1-a), (1, a))
            for dj, v in ((0, 1-b), (1, b)) if u*v > 0]


def projected_area(radius: float, length: float, k: Any) -> float:
    k = unit(k)
    normals = np.array([[1,0,0], [.5,0,math.sqrt(3)/2], [-.5,0,math.sqrt(3)/2],
                        [-1,0,0], [-.5,0,-math.sqrt(3)/2], [.5,0,-math.sqrt(3)/2], [0,1,0], [0,-1,0]])
    areas = np.array([radius*length]*6+[3*math.sqrt(3)*radius*radius/2]*2)
    return float(np.dot(areas, np.maximum(0, -normals@k)))


class Atlas:
    def __init__(self, directory: str | Path):
        self.directory = Path(directory)
        if self.directory.name.endswith('.part'):
            raise ValueError('Incomplete atlas directory')
        with (self.directory/'metadata.json').open(encoding='utf-8') as f:
            self.meta = json.load(f)
        m = self.meta
        if m.get('schema') != 'ice.single_particle_atlas.v1' or m.get('complete') is not True:
            raise ValueError('Incomplete or unknown atlas schema')
        if m.get('incident_reconstruction') != 'bilinear_theta_phi_in_spherical_incident_frame' or m.get('symmetry') != 'D6h_O3_full_24_scalar_unpolarized':
            raise ValueError('Unsupported atlas frame/interpolation/symmetry')
        self.it, self.ip = [positive_int(m[k], k, 4096) for k in ('incident_theta_intervals', 'incident_phi_intervals')]
        self.nt, self.np = [positive_int(m[k], k) for k in ('theta_count', 'phi_count')]
        self.mass = np.load(self.directory/'mass_mm2.npy', mmap_mode='r', allow_pickle=False)
        self.edges = np.load(self.directory/'u_edges.npy', allow_pickle=False)
        shape = (self.it+1, self.ip+1, self.nt, self.np)
        if self.mass.shape != shape or self.mass.dtype.str != '<f8' or not self.mass.flags.c_contiguous:
            raise ValueError('Atlas mass layout does not match metadata')
        if self.edges.shape != (self.nt+1,) or self.edges.dtype.str != '<f8' or not np.isfinite(self.edges).all() or self.edges[0] != 0 or self.edges[-1] != 1 or np.any(np.diff(self.edges) <= 0):
            raise ValueError('Invalid atlas u edges')
        if len(m['nodes']) != (self.it+1)*(self.ip+1):
            raise ValueError('Missing incident atlas nodes')
        self.area = np.zeros((self.it+1, self.ip+1))
        self.residual = np.zeros_like(self.area)
        g = m['geometry']; self.radius, self.length = float(g['circumradius_mm']), float(g['full_length_mm'])
        if g['shape'] != 'regular_hexagonal_prism' or not (math.isfinite(self.radius) and math.isfinite(self.length) and self.radius > 0 and self.length > 0):
            raise ValueError('Unsupported geometry')
        for index, n in enumerate(m['nodes']):
            i, j = divmod(index, self.ip+1)
            if (n['i'], n['j']) != (i, j):
                raise ValueError('Node order mismatch')
            f = canonical_frame(i*(PI/2)/self.it, j*(PI/6)/self.ip)
            if not np.allclose(n['ki'], f[:,2], rtol=0, atol=2e-12) or not np.allclose(n['e0'], f[:,0], rtol=0, atol=2e-12):
                raise ValueError('Node frame mismatch')
            cell = self.mass[i,j]
            if not np.isfinite(cell).all() or np.any(cell < 0):
                raise ValueError('Invalid atlas cell mass')
            area, residual, escaped = [float(n[k]) for k in ('projected_area_mm2', 'unresolved_mm2', 'escaped_mm2')]
            if not all(math.isfinite(v) for v in (area, residual, escaped)) or area <= 0 or escaped <= 0 or residual < 0:
                raise ValueError('Invalid node accounting')
            if abs((escaped+residual)/area-1) > 2e-10 or abs(float(np.sum(cell))/escaped-1) > 2e-10 or abs(area/projected_area(self.radius,self.length,n['ki'])-1)>2e-12:
                raise ValueError('Node area / mass / residual audit failed')
            self.area[i,j], self.residual[i,j] = area, residual
        self.hashes = {name: sha256_file(self.directory/name) for name in ('metadata.json','mass_mm2.npy','u_edges.npy')}

    def close(self) -> None:
        if isinstance(self.mass, np.memmap):
            self.mass._mmap.close()


def coupling(atlas: Atlas, rotations: np.ndarray, number_weights: np.ndarray, ki: Any) -> tuple[np.ndarray, float, np.ndarray]:
    """Assemble nonnegative convolution weights without tracing any new rays.

    At each query s, interpolate M_a/A_a and multiply the EXACT projected A(s).
    Division by A_a is NOT division by resolved escaped mass; tails are retained.
    Output weights [node_theta,node_phi,parity,cyclic_shift] multiply RAW M_a.
    """
    validate_rotations(rotations, number_weights)
    output_frame = frame(ki); k = output_frame[:,2]
    w = np.zeros((atlas.it+1, atlas.ip+1, 2, atlas.np), dtype=np.float64)
    expected = []
    for r, nw in zip(rotations, number_weights):
        body = r.T@k
        theta, phi, q = fold_direction(body)
        f = canonical_frame(theta, phi)
        e = r @ q @ f[:,0]
        shift = math.atan2(float(e@output_frame[:,1]), float(e@output_frame[:,0]))/(2*PI)*atlas.np
        nearest = round(shift)
        if abs(shift-nearest) <= 64*EPS*atlas.np:
            shift = float(nearest)  # rounding-sized integer shift, not a blur
        lower = math.floor(shift); fraction = shift-lower
        parity = int(np.linalg.det(q) < 0)
        a = projected_area(atlas.radius, atlas.length, body)
        expected.append(float(nw)*a)
        for i, j, mix in interpolation(theta, phi, atlas.it, atlas.ip):
            c = float(nw)*a*mix/atlas.area[i,j]
            w[i,j,parity,lower%atlas.np] += c*(1-fraction)
            w[i,j,parity,(lower+1)%atlas.np] += c*fraction
    expected_area = math.fsum(expected)
    reconstructed_area = float(np.sum(np.sum(w,axis=(2,3))*atlas.area))
    if abs(reconstructed_area/expected_area-1) > 1e-11 or np.any(w < 0) or not np.isfinite(w).all():
        raise ValueError('Orientation-to-convolution area audit failed')
    return w, expected_area, output_frame



def coupling_for_rule(atlas: Atlas, rotations: np.ndarray,
                      number_weights: np.ndarray, ki: Any,
                      rule_meta: dict[str, Any]
                      ) -> tuple[np.ndarray, float, np.ndarray, dict[str, Any]]:
    """Integrate the uniform world-azimuth angle exactly at vertical incidence.

    For R=Ry(alpha) Rz(beta) Ry(gamma) and ki=+/-y, R.T@ki and
    projected area are independent of alpha. Alpha only cyclically translates
    outgoing phi. Averaging uniformly over alpha is therefore a constant
    circular-convolution kernel, exactly for the declared finite-cell atlas
    reconstruction. Tilt/spin and incident interpolation remain numerical.

    This applies ONLY to declared isotropic/horizontal_gaussian components.
    Fixed poses are not averaged, even when mixed with a continuous component.
    Nonvertical directions (including near-vertical ones) use the original rule.
    The full arrays supplied by orientation_rule are the provenance/reference
    quadrature. At vertical incidence its redundant alpha nodes are replaced
    by an analytic orbit integral; the returned diagnostics document this.
    """
    k = unit(ki)
    vertical = bool(k[0] == 0.0 and k[2] == 0.0)
    config = rule_meta['config']
    components = config['components']
    eligible = [c['kind'] in ('isotropic', 'horizontal_gaussian')
                for c in components]
    if not vertical or not any(eligible):
        w, area, f = coupling(atlas, rotations, number_weights, k)
        return w, area, f, {
            'method': 'original_fixed_orientation_quadrature',
            'exact_vertical_incidence': vertical,
            'analytic_uniform_azimuth_components': [],
            'effective_pose_evaluations': len(rotations),
            'symmetry_axis_ensemble': [0, 1, 0],
            'direction_snapped': False,
            'output_post_filter': False,
        }

    # Never infer invariance from a coincidentally symmetric list of poses.
    # Ensure the model declaration and the supplied provenance arrays agree.
    reference_rotations, reference_weights, _ = orientation_rule(config)
    if (not np.array_equal(rotations, reference_rotations) or
            not np.array_equal(number_weights, reference_weights)):
        raise ValueError('Analytic azimuth requires the declared orientation_rule arrays')
    result = np.zeros((atlas.it+1, atlas.ip+1, 2, atlas.np), dtype=np.float64)
    areas = []
    evaluations = 0
    component_modes = []
    for component, analytic in zip(components, eligible):
        fraction = float(component['fraction'])
        one = dict(component, fraction=1.0)
        quadrature = dict(config['quadrature'])
        if analytic:
            # A single representative alpha suffices; its azimuth is immaterial
            # after integration over the whole uniform orbit.
            quadrature['azimuth'] = 1
        r, nw, _ = orientation_rule({'quadrature': quadrature, 'components': [one]})
        wc, ac, output_frame = coupling(atlas, r, nw, k)
        if analytic:
            # Replace the orbit's discrete delta shifts, not the optical model.
            # For every tilt/spin contribution, integral d alpha/(2*pi) of
            # an azimuth-shifted cell distributes its mass equally to all phi.
            orbit_total = np.sum(wc, axis=-1, keepdims=True, dtype=np.float64)
            wc = np.broadcast_to(orbit_total / atlas.np, wc.shape)
        result += fraction * wc
        areas.append(fraction * ac)
        evaluations += len(r)
        component_modes.append('analytic_uniform_alpha' if analytic else 'fixed_pose')
    area = math.fsum(areas)
    represented = float(np.sum(np.sum(result, axis=(2, 3)) * atlas.area))
    if (not area > 0 or not np.isfinite(result).all() or np.any(result < 0) or
            abs(represented/area-1) > 1e-11):
        raise ValueError('Analytic azimuth coupling failed area audit')
    return result, area, output_frame, {
        'method': 'analytic_uniform_world_azimuth_at_vertical_incidence',
        'exact_vertical_incidence': True,
        'analytic_uniform_azimuth_components': [i for i, x in enumerate(eligible) if x],
        'component_methods': component_modes,
        'reference_pose_count': len(rotations),
        'effective_pose_evaluations': evaluations,
        'symmetry_axis_ensemble': [0, 1, 0],
        'direction_snapped': False,
        'output_post_filter': False,
        'tilt_spin_and_atlas_interpolation_remain_numerical': True,
    }


def compose_mass(atlas: Atlas, weights: np.ndarray) -> tuple[np.ndarray, float]:
    out = np.zeros((atlas.nt,atlas.np),dtype=np.float64)
    # Fixed order, nonnegative direct cyclic convolution; no FFT negatives/floor.
    for i in range(atlas.it+1):
        for j in range(atlas.ip+1):
            for parity in range(2):
                cell = np.asarray(atlas.mass[i,j])
                if parity:
                    cell = cell[:,::-1]  # -phi maps cell j to Nphi-1-j
                for shift, value in enumerate(weights[i,j,parity]):
                    if value > 0:
                        out += float(value)*np.roll(cell,shift,axis=1)
    residual = float(np.sum(np.sum(weights,axis=(2,3))*atlas.residual))
    if not np.isfinite(out).all() or np.any(out < 0):
        raise ValueError('Invalid composed mass')
    return out, residual


def make_cdf(mass: np.ndarray, edges: np.ndarray, max_l1: float = 1e-10, max_lost: float = 1e-12) -> tuple[np.ndarray,np.ndarray,dict[str,Any]]:
    """Ordered FP64 phi marginal + theta|phi CDF, audited against input masses.

    Same external contract as rainbow. Input is integrated mass, not density.
    No floors, repairs, smoothing, conditional-CDF averaging, or hidden coarsening.
    """
    if mass.ndim != 2 or not np.isfinite(mass).all() or np.any(mass < 0) or not float(np.max(mass))>0:
        raise ValueError('CDF requires finite nonnegative nonempty mass')
    nt,np_ = mass.shape
    if edges.shape!=(nt+1,) or not np.isfinite(edges).all() or edges[0]!=0 or edges[-1]!=1 or np.any(np.diff(edges)<=0):
        raise ValueError('Invalid CDF edges')
    raw = mass / float(np.max(mass))
    if np.any((mass>0)&(raw==0)):
        raise ValueError('CDF rescaling underflow')
    b = np.zeros((np_,nt+1),dtype=np.float64)
    b[:,1:] = np.cumsum(raw.T,axis=1)
    cols = b[:,-1].copy()
    a = np.zeros(np_+1,dtype=np.float64); a[1:] = np.cumsum(cols)
    if not a[-1]>0 or not math.isfinite(float(a[-1])):
        raise ValueError('Invalid CDF total')
    a /= a[-1]
    for j in range(np_):
        if cols[j]>0:
            b[j] /= cols[j]
        else:
            b[j] = edges
    if np.any(np.diff(a)<0) or np.any(np.diff(b,axis=1)<0) or a[0]!=0 or a[-1]!=1 or np.any(b[:,0]!=0) or np.any(b[:,-1]!=1):
        raise ValueError('Malformed CDF, refusing repair')
    total = math.fsum(mass.flat)
    target = mass/total
    recovered = np.diff(a)[:,None]*np.diff(b,axis=1)
    if np.any((mass>0)&(target==0)):
        raise ValueError('Normalized mass underflow')
    errors = np.abs(target.T-recovered)
    l1 = math.fsum(errors.flat)
    lost = math.fsum(target.T[(target.T>0)&(recovered==0)])
    g = float(np.sum(recovered*(1-edges[:-1]-edges[1:])[None,:],dtype=np.longdouble)/np.sum(recovered,dtype=np.longdouble))
    cdf_mass = math.fsum(recovered.flat)
    if abs(cdf_mass-1)>1e-10 or l1>max_l1 or lost>max_lost or not -1<g<1:
        raise ValueError(f'CDF audit rejected l1={l1} lost={lost} g={g}')
    return a,b,{'cdf_mass': cdf_mass,'cdf_l1_mass_error':l1,
        'cdf_lost_probability_mass':lost,'cdf_maximum_cell_mass_error':float(np.max(errors)),
        'g':g,'resolved_mass_mm2':total}
