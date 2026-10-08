#!/usr/bin/env python3
"""Read-only validation/plots and exact-cell integration of saved phase CDFs.

Run from the repository or by absolute script path. Never alters the input CDF,
never re-normalizes it, never adds a PDF floor, and never merges unequal u edges.
The comparison target is the saved piecewise-constant density, NOT a Maxwell
solution or an unbinned geometric-optics distribution. See phase_analysis.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
import sys
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'python'))
from rainbow_dataset.record import PhaseRecord

VERSION = 'ice.phase_analysis.v1'
FILES = ('metadata.json', 'phi_cdf.npy', 'theta_given_phi_cdf.npy', 'u_edges.npy')


def checked_int(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return result


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def total(values: np.ndarray) -> float:
    # fsum avoids depending on longdouble being wider than double (MSVC).
    return math.fsum(map(float, np.ravel(values)))


def theta_of_u(u: np.ndarray) -> np.ndarray:
    # Stable at both poles. Only the coordinate visualization uses this inverse.
    return np.where(u <= 0.5, 2*np.arcsin(np.sqrt(u)),
                    np.pi-2*np.arcsin(np.sqrt(1-u)))


@dataclass
class Table:
    directory: Path
    metadata: dict[str, Any]
    frame: np.ndarray
    u: np.ndarray
    v: np.ndarray
    mass: np.ndarray  # [phi_cell, theta_cell], the saved distribution
    hashes: dict[str, str]

    @property
    def du(self) -> np.ndarray:
        return np.diff(self.u)

    @property
    def dv(self) -> np.ndarray:
        # Contract: equal-phi cells have exactly the common factor 1/Nphi.
        # Do not replace it by differences of rounded plotting coordinates.
        return np.full(self.mass.shape[0], 1/self.mass.shape[0])

    @property
    def density_uv(self) -> np.ndarray:
        return self.mass / (self.du[None, :] * self.dv[:, None])

    @property
    def density_omega(self) -> np.ndarray:
        return self.density_uv / (4*np.pi)

    @cached_property
    def theta_mass(self) -> np.ndarray:
        return np.array([math.fsum(map(float, row)) for row in self.mass.T])

    @property
    def g(self) -> float:
        return total(self.mass * (1-self.u[:-1]-self.u[1:])) / total(self.mass)


def load_table(path: str | Path, max_cells: int = 10_000_000) -> Table:
    p = Path(path).resolve()
    if max_cells <= 0:
        raise ValueError('max_cells must be positive')
    if not p.is_dir():
        raise ValueError('Input must be an extracted record directory, not a ZIP.')
    # Require finite JSON numbers and reject duplicate keys.
    def object_pairs(pairs):
        d = {}
        for k, v in pairs:
            if k in d:
                raise ValueError(f'Duplicate metadata key: {k}')
            d[k] = v
        return d
    def invalid_constant(s):
        raise ValueError(f'Nonfinite metadata number: {s}')
    if (p / 'metadata.json').stat().st_size > 1024*1024:
        raise ValueError('Metadata exceeds 1 MiB')
    with (p / 'metadata.json').open(encoding='utf-8') as f:
        m = json.load(f, object_pairs_hook=object_pairs, parse_constant=invalid_constant)
    def finite_json(x):
        if isinstance(x, float) and not math.isfinite(x):
            raise ValueError('Nonfinite metadata number')
        if isinstance(x, dict):
            for value in x.values(): finite_json(value)
        elif isinstance(x, list):
            for value in x: finite_json(value)
    finite_json(m)
    for key, expected in [('direction_convention', 'physical_propagation'),
                          ('theta_zero', 'forward'), ('theta_pi', 'backward')]:
        if m.get(key) != expected:
            raise ValueError(f'Missing or unsupported {key}; no direction convention is inferred.')
    if m.get('coordinate_contract_id') not in {
            'ice.phase_cdf.ensemble_coordinates.v1', 'rainbow.phase_cdf.coordinates.v1'}:
        raise ValueError('Unknown coordinate contract')
    nt, nf = m.get('theta_count'), m.get('phi_count')
    if type(nt) is not int or type(nf) is not int or min(nt, nf) < 1:
        raise ValueError('Invalid dimensions')
    if nt*nf > max_cells:
        raise ValueError(f'{nt*nf} cells exceeds --max-cells={max_cells}; no automatic coarsening.')
    with PhaseRecord(p, validate=True) as r:
        a = np.array(r.phi_cdf, dtype=np.float64)
        u = np.array(r.u_edges, dtype=np.float64)
        mass = np.diff(r.theta_cdf, axis=1) * np.diff(a)[:, None]
        frame = np.array(r.frame, dtype=np.float64)
    if not np.isfinite(mass).all() or (mass < 0).any():
        raise ValueError('Invalid reconstructed cell probability')
    # The current schema fixes v_j=j/Nphi. No implicit phi rotation is allowed.
    v = np.arange(nf+1, dtype=np.float64)/nf
    result = Table(p, m, frame, u, v, mass, {name: sha256(p/name) for name in FILES})
    if not np.isfinite(result.density_omega).all():
        raise ValueError('Density cannot be represented; input not modified.')
    for name, column in [('query_frame_axis', 2), ('query_frame_e0', 0)]:
        if name in m:
            q = np.asarray(m[name], dtype=float)
            if q.shape != (3,) or not np.allclose(q, frame[:, column], rtol=0, atol=2e-12):
                raise ValueError(f'{name} disagrees with sampling frame')
    return result


def cone_mass(t: Table, degrees: float) -> float:
    if not 0 <= degrees <= 180:
        raise ValueError('Cone angle must be in [0,180]')
    edge = math.sin(math.radians(degrees)/2)**2
    # Intersections with cells in u, not linear interpolation in theta.
    overlap = np.maximum(0, np.minimum(t.u[1:], edge)-t.u[:-1])
    return total(t.theta_mass * overlap/t.du)


def angular_summary(t: Table) -> dict[str, Any]:
    m, w = t.metadata, t.mass
    nphi, ntheta = w.shape
    d = t.density_omega
    marginal = t.theta_mass
    theta = np.degrees(theta_of_u(t.u))
    order = np.argsort(-marginal, kind='stable')[:min(10, ntheta)]
    positive = d[d > 0]
    warnings: list[str] = []
    orientation = m.get('orientation', {})
    if orientation.get('rows', 0) <= 2:
        warnings.append('At most two discrete orientations: not a resolved continuous orientation distribution.')
    tv = m.get('storage_processing', {}).get('coarsening_tv')
    if tv is not None and tv > 0.05:
        warnings.append('Logged coarsening TV exceeds the advisory value 0.05; mass conservation does not preserve angular shape.')
    if not m.get('quality', {}).get('angular_convergence_certified', False):
        warnings.append('Angular convergence has not been certified. One record cannot establish ray/orientation/grid convergence.')
    if not orientation.get('plan_embedded', False):
        warnings.append('Orientation plan is not embedded; preserve the original CSV and its hash.')
    if m.get('material', {}).get('automatic_ice_dispersion') is False:
        warnings.append('Index was supplied by the caller; wavelength alone does not prove use of an ice dispersion model.')
    if m.get('histogram', {}).get('forward_axis_mass_mm2', 0) > 0:
        warnings.append('Forward-axis mass is spread over a finite polar cap by storage; this is not physical diffraction.')
    stored_g = m.get('hg', {}).get('g')
    gr = m.get('histogram', {}).get('g_rays')
    gs = m.get('storage_processing', {}).get('g_source')
    audit = m.get('transport_audit', {})
    area = m.get('geometry', {}).get('projected_area_mm2')
    energy_check = None
    if area and 'escaped_mm2' in audit and 'unresolved_mm2' in audit:
        energy_check = (audit['escaped_mm2']+audit['unresolved_mm2'])/area-1
    # Only recomputable values appear in 'computed_from_saved_arrays'.
    return {
        'analysis_version': VERSION,
        'input_directory': str(t.directory),
        'sha256': t.hashes,
        'reader_validation': 'passed',
        'computed_from_saved_arrays': {
            'shape_phi_theta': [nphi, ntheta],
            'probability_mass': total(w),
            'g': t.g,
            'g_minus_metadata': None if stored_g is None else t.g-stored_g,
            'positive_cells': int(np.count_nonzero(w)),
            'zero_cells': int(np.count_nonzero(w == 0)),
            'zero_cell_fraction': float(np.mean(w == 0)),
            'positive_support_solid_angle_sr': total((w > 0)*4*np.pi*t.dv[:, None]*t.du[None, :]),
            'pdf_min_positive_sr_inverse': float(positive.min()),
            'pdf_max_sr_inverse': float(positive.max()),
            'first_polar_band_mass': float(marginal[0]),
            'last_polar_band_mass': float(marginal[-1]),
            'forward_hemisphere_mass': cone_mass(t, 90),
            'cone_probabilities': {str(a): cone_mass(t, a) for a in [1, 2.5, 5, 10, 22, 30, 46, 60, 90, 180]},
            'highest_mass_theta_bands': [{'lower_deg': float(theta[i]), 'upper_deg': float(theta[i+1]),
                                         'mass': float(marginal[i])} for i in order],
        },
        'computed_from_metadata_numbers_not_retraced': {
            'energy_balance_relative': energy_check,
            'g_source_minus_g_rays': None if gr is None or gs is None else gs-gr,
            'g_saved_minus_g_source': None if gs is None else t.g-gs,
            'g_saved_minus_g_rays': None if gr is None else t.g-gr,
        },
        'producer_metadata_not_independently_recomputed': {
            k: m.get(k) for k in ['backend', 'geometry', 'material', 'orientation', 'sampling',
                                 'transport_audit', 'histogram', 'storage_processing', 'quality']},
        'warnings': warnings,
        'scope': 'Saved finite-cell probability distribution only. No claim of physical or sampling convergence.',
    }


def check_frames(a: Table, b: Table) -> float:
    if a.metadata.get('coordinate_contract_id') != b.metadata.get('coordinate_contract_id'):
        raise ValueError('Coordinate contracts differ; comparison would require an explicit transformation.')
    if a.metadata.get('sampling_frame_map') != b.metadata.get('sampling_frame_map'):
        raise ValueError('Coordinate spaces differ')
    error = float(np.max(np.abs(a.frame-b.frame)))
    # Treat only floating-point frame roundoff as equivalent; never rotate/blur.
    if error > 2e-12:
        raise ValueError('Sampling frames differ; no automatic rotation/interpolation is performed.')
    return error


def compare_tables(a: Table, b: Table, max_refined_cells: int = 50_000_000) -> dict[str, Any]:
    frame_error = check_frames(a, b)
    u = np.union1d(a.u, b.u)
    # Common rational phi lattice prevents artificial slivers from j/N rounding.
    na, nb = a.mass.shape[0], b.mass.shape[0]
    denominator = math.lcm(na, nb)
    if denominator > np.iinfo(np.int64).max:
        raise ValueError('Common phi lattice exceeds int64 range')
    sa, sb = denominator//na, denominator//nb
    v = np.union1d(np.arange(na+1, dtype=np.int64)*sa,
                   np.arange(nb+1, dtype=np.int64)*sb)
    count = (len(u)-1)*(len(v)-1)
    if max_refined_cells <= 0 or count > max_refined_cells:
        raise ValueError(f'{count} common-refinement cells exceeds limit; no downsampling or edge snapping.')
    # Index by the LEFT boundary, not a midpoint which can round onto an edge
    # for one-ULP intervals in the common refinement.
    ia = np.searchsorted(a.u, u[:-1], side='right')-1
    ib = np.searchsorted(b.u, u[:-1], side='right')-1
    ja = v[:-1]//sa
    jb = v[:-1]//sb
    da, db = a.density_uv, b.density_uv
    du, dv = np.diff(u), np.diff(v).astype(np.float64)/denominator
    sums = {k: [] for k in ['l1', 'l2', 'h2', 'js', 'a_only', 'b_only', 'mass_a', 'mass_b']}
    # One refined phi row at a time. Density is constant on each intersection.
    for row in range(len(v)-1):
        p, q = da[ja[row], ia], db[jb[row], ib]
        area = du*dv[row]
        diff = p-q
        sums['l1'].append(total(area*np.abs(diff)))
        sums['l2'].append(total(area*diff*diff/(4*np.pi)))
        sums['h2'].append(total(0.5*area*(np.sqrt(p)-np.sqrt(q))**2))
        mixture = 0.5*(p+q)
        delta = np.zeros_like(p)
        use = mixture > 0
        delta[use] = (p[use]-q[use])/(p[use]+q[use])
        # Jensen-Shannon in a non-cancelling form. For tiny differences,
        # F(t)=sum_{n>=1} t^(2n)/(2n*(2n-1)); 5 terms suffice at |t|<0.01.
        ad = np.abs(delta); f = np.zeros_like(ad)
        small = (ad < 0.01)
        x = ad[small]*ad[small]
        f[small] = x*(0.5+x*(1/12+x*(1/30+x*(1/56+x/90))))
        medium = (ad >= 0.01) & (ad < 1)
        x = ad[medium]
        f[medium] = 0.5*((1+x)*np.log1p(x)+(1-x)*np.log1p(-x))
        f[ad == 1] = np.log(2.0)
        sums['js'].append(total(area*mixture*f))
        sums['a_only'].append(total(area*np.where(q == 0, p, 0)))
        sums['b_only'].append(total(area*np.where(p == 0, q, 0)))
        sums['mass_a'].append(total(area*p)); sums['mass_b'].append(total(area*q))
    s = {k: math.fsum(x) for k, x in sums.items()}
    if not all(math.isfinite(x) for x in s.values()):
        raise ValueError('Nonfinite metric; values were not clamped or re-normalized')
    ma, mb = a.theta_mass/a.du, b.theta_mass/b.du
    marginal_tv = 0.5*total(du*np.abs(ma[ia]-mb[ib]))
    return {
        'analysis_version': VERSION,
        'target': 'saved_piecewise_constant_solid_angle_densities',
        'integration': 'union_of_original_u_edges_and_exact_rational_equal_phi_edges; no bin merging, resampling, smoothing, or normalization',
        'common_refinement_shape_phi_theta': [len(v)-1, len(u)-1],
        'frame_max_absolute_difference': frame_error,
        'frame_equivalence_absolute_tolerance': 2e-12,
        'probability_mass_a_integrated': s['mass_a'], 'probability_mass_b_integrated': s['mass_b'],
        'l1': s['l1'], 'tv': 0.5*s['l1'],
        'l2_per_sqrt_sr': math.sqrt(s['l2']),
        'hellinger': math.sqrt(s['h2']), 'jensen_shannon_nats': s['js'],
        'mass_a_where_b_is_zero': s['a_only'], 'mass_b_where_a_is_zero': s['b_only'],
        'theta_marginal_tv': marginal_tv,
        'g_a': a.g, 'g_b': b.g, 'g_a_minus_b': a.g-b.g,
        'input_a_sha256': a.hashes, 'input_b_sha256': b.hashes,
        'scope': 'No confidence interval from a single pair; does not certify physical or statistical convergence.',
    }


def physical_differences(a: Table, b: Table) -> list[str]:
    paths = ['wavelength_nm', 'interior_index', 'exterior_index', 'stage', 'input_polarization',
             'material', 'diffraction_policy', 'inter_path_interference',
             'geometry.shape', 'geometry.circumradius_mm', 'geometry.full_length_mm']
    def at(m, key):
        for part in key.split('.'):
            if not isinstance(m, dict):
                return None
            m = m.get(part)
        return m
    return [p for p in paths if at(a.metadata, p) != at(b.metadata, p)]


def write_json(path: Path, data: Any) -> None:
    with path.open('x', encoding='utf-8', newline='\n') as f:
        json.dump(data, f, indent=2, ensure_ascii=True, allow_nan=False)
        f.write('\n')


def make_out(path: str | Path) -> Path:
    out = Path(path).resolve()
    # Refuse even an existing empty directory: no accidental mixing of runs.
    out.mkdir(parents=True, exist_ok=False)
    return out


def plot_record(t: Table, out: Path, label: str) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    theta = np.degrees(theta_of_u(t.u)); phi = 360*t.v-180
    pdf = t.density_omega.T
    # Every chart is a distinct figure. Default colormap/colors; no interpolation.
    for log in (False, True):
        fig, ax = plt.subplots(figsize=(10, 5.4), layout='constrained')
        values = np.ma.masked_equal(pdf, 0) if log else pdf
        image = ax.pcolormesh(phi, theta, values, shading='flat', norm=LogNorm() if log else None)
        ax.set(xlabel='Azimuth phi [degrees]', ylabel='Scattering polar angle theta [degrees]',
               xlim=(-180, 180), ylim=(180, 0))
        ax.set_title(f'{label}: saved PDF ({"logarithmic" if log else "linear"} scale)\n'
                     f'{pdf.shape[0]} x {pdf.shape[1]} cells; no smoothing'
                     + ('; blank cells = exact zero' if log else ''))
        fig.colorbar(image, ax=ax, label='p(omega) [sr^-1]')
        fig.savefig(out/f'pdf_{"log" if log else "linear"}.png', dpi=180)
        plt.close(fig)
    # Average over azimuth is a density per sr, NOT a theta marginal PDF.
    mean_pdf = t.theta_mass/(4*np.pi*t.du)
    fig, ax = plt.subplots(figsize=(10, 4.5), layout='constrained')
    ax.stairs(mean_pdf, theta, label='Azimuth-averaged PDF')
    ax.set(xlabel='Scattering polar angle theta [degrees]', ylabel='Azimuth-averaged p [sr^-1]', xlim=(0, 180))
    ax.set_title(f'{label}: azimuth average; each step is a stored theta band')
    fig.savefig(out/'azimuth_average.png', dpi=180); plt.close(fig)
    # Exact cell curve for marginal density with respect to theta, shown per deg.
    # p_theta(theta)=m_i*sin(theta)/(2*du_i), not a constant m_i/dtheta.
    fig, ax = plt.subplots(figsize=(10, 4.5), layout='constrained')
    xs, ys = [], []
    for i, probability in enumerate(t.theta_mass):
        x = np.linspace(theta[i], theta[i+1], 16)
        y = probability*np.sin(np.radians(x))/(2*t.du[i])*np.pi/180
        xs.extend([*x, math.nan]); ys.extend([*y, math.nan])
    ax.plot(xs, ys, label='Theta marginal')
    ax.set(xlabel='Scattering polar angle theta [degrees]', ylabel='Probability density per degree', xlim=(0, 180))
    ax.set_title(f'{label}: theta marginal (includes spherical Jacobian)')
    fig.savefig(out/'theta_marginal.png', dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4.5), layout='constrained')
    x = np.linspace(0, 180, 1801)
    ax.plot(x, [cone_mass(t, float(angle)) for angle in x])
    ax.set(xlabel='Forward-cone half-angle [degrees]', ylabel='Probability inside cone', xlim=(0, 180), ylim=(0, 1))
    ax.set_title(f'{label}: integrated forward-cone probability')
    fig.savefig(out/'cone_cdf.png', dpi=180); plt.close(fig)


def plot_comparison(a: Table, b: Table, out: Path) -> None:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 4.5), layout='constrained')
    for t, name in [(a, 'A'), (b, 'B')]:
        ax.stairs(t.theta_mass/(4*np.pi*t.du), np.degrees(theta_of_u(t.u)), label=name)
    ax.set(xlabel='Scattering polar angle theta [degrees]', ylabel='Azimuth-averaged p [sr^-1]', xlim=(0,180))
    ax.set_title('Saved distributions A and B: azimuth-averaged PDFs')
    ax.legend(); fig.savefig(out/'comparison_azimuth_average.png', dpi=180); plt.close(fig)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    one = sub.add_parser('inspect', help='Validate and plot one extracted CDF record')
    one.add_argument('record'); one.add_argument('--label', default='Phase CDF')
    two = sub.add_parser('compare', help='Integrate differences on a common refinement')
    two.add_argument('a'); two.add_argument('b')
    two.add_argument('--allow-orientation-plan-difference', action='store_true',
                     help='Acknowledge different/missing orientation hashes; does not assert the same physical distribution')
    two.add_argument('--allow-physical-differences', action='store_true')
    two.add_argument('--max-refined-cells', type=checked_int, default=50_000_000)
    for q in (one, two):
        q.add_argument('--out', required=True, help='New output directory; existing paths are refused')
        q.add_argument('--no-plots', action='store_true')
        q.add_argument('--max-cells', type=checked_int, default=10_000_000)
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == 'inspect':
            t = load_table(args.record, args.max_cells)
            summary = angular_summary(t)
            out = make_out(args.out)
            write_json(out/'analysis.json', summary)
            profile = {'theta_edges_deg': np.degrees(theta_of_u(t.u)).tolist(),
                       'u_edges': t.u.tolist(), 'theta_band_probability': t.theta_mass.tolist(),
                       'azimuth_average_pdf_sr_inverse': (t.theta_mass/(4*np.pi*t.du)).tolist()}
            write_json(out/'theta_profile.json', profile)
            if not args.no_plots:
                plot_record(t, out, args.label)
            print(f'PASS saved-CDF validation: mass={total(t.mass):.17g}, g={t.g:.17g}')
            for warning in summary['warnings']:
                print('NOTE:', warning)
        else:
            a, b = load_table(args.a, args.max_cells), load_table(args.b, args.max_cells)
            difference = physical_differences(a, b)
            if difference and not args.allow_physical_differences:
                raise ValueError('Physical metadata differ: '+', '.join(difference)+'; explicit --allow-physical-differences required')
            ha = a.metadata.get('orientation', {}).get('csv_sha256')
            hb = b.metadata.get('orientation', {}).get('csv_sha256')
            equal = ha is not None and ha == hb
            if not equal and not args.allow_orientation_plan_difference:
                raise ValueError('Orientation CSV hashes differ or are absent; explicit --allow-orientation-plan-difference required')
            summary = compare_tables(a, b, args.max_refined_cells)
            summary['physical_metadata_differences'] = difference
            summary['orientation_csv_sha256_equal'] = equal
            summary['orientation_plan_difference_acknowledged'] = args.allow_orientation_plan_difference
            out = make_out(args.out)
            write_json(out/'comparison.json', summary)
            if not args.no_plots:
                plot_comparison(a, b, out)
            print(f'PASS saved-density comparison: TV={summary["tv"]:.17g}, delta_g={summary["g_a_minus_b"]:.17g}')
        write_json(out/'analysis_complete.json', {'complete': True, 'analysis_version': VERSION,
                                                'tool_sha256': sha256(Path(__file__)),
                                                'python_version': sys.version, 'numpy_version': np.__version__,
                                                'input_modified': False, 'plots_created': not args.no_plots})
        return 0
    except (OSError, ValueError, MemoryError, ImportError) as e:
        # ASCII escapes prevent Windows CP932 from hiding the real diagnostic.
        print('ERROR: '+ascii(str(e)), file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
