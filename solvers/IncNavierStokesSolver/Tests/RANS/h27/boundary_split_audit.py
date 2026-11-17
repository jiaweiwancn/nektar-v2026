#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only diagnostic of unresolved prescribed MMS k boundary quadrature.

Compare the historical unsplit Gauss16 natural-flux callback action with
Gauss24/32 on faces split at the EXACT MMS cutoff and fixed material heights.
This changes no native solution, tolerance, model, or gate. Optional current
residual reconstruction uses the independent weak operator and the recorded
explicit history loads; it does not repeat their separate history verification.
"""
import argparse
import hashlib
import inspect
import json
import math
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
import coupled_mms_check
import coupled_oracle


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class BoundaryAudit:
    def __init__(self, data, table, deadline):
        self.data, self.deadline = data, deadline
        amplitude = float(data.get('mms_body_amplitude', 0.))
        require(math.isfinite(amplitude), 'nonfinite declared MMS body amplitude')
        if amplitude == 0:
            # Source-v6's frozen helper has only IndependentMMS(table). Calling
            # that shared historical API also preserves the new helper's exact
            # zero-body default, without editing either imported source.
            self.m = coupled_mms_check.IndependentMMS(table)
            self.mms_api = 'legacy-compatible single-argument zero-body constructor'
        else:
            require('body_amplitude' in inspect.signature(coupled_mms_check.IndependentMMS).parameters,
                    'nonzero-body artifact requires a body-aware independent MMS helper; frozen source-v6 supports only body zero')
            self.m = coupled_mms_check.IndependentMMS(table, body_amplitude=amplitude)
            self.mms_api = 'explicit body-aware constructor'
        self.nk = int(data['sizes'][2]); self.elements = data['elements']
        self.free = np.asarray(data['k_free_ids'], dtype=int)
        self.when = float(data['time']); self.faces = []
        require(data['manufactured'] is True and 0 < self.nk <= 4096,
                'diagnostic requires a bounded manufactured fixture (k size <=4096)')
        require(math.isfinite(self.when), 'nonfinite native sample time')
        expected = dict(left=0., right=1., top=.5, d0=self.m.d0, H=self.m.H, h=self.m.h)
        for key, value in expected.items():
            require(abs(data['geometry'][key]-value) <= 1e-12, 'wrong MMS geometry: '+key)
        require(len({e['geometry'] for e in self.elements}) == len(self.elements), 'duplicate geometry')
        visited, essential = set(), set()
        self.mass = np.zeros((self.nk, self.nk)); self.volume = 0.
        a, w = leggauss(8)
        for e in self.elements:
            self.clock()
            xl, xr, rl, rr = e['bounds']; ids = np.asarray(e['k'], dtype=int)
            require(len(ids) == len(e['ks']) == 16 and len(set(ids)) == 16 and
                    np.all((ids >= 0) & (ids < self.nk)) and all(v in (-1., 1.) for v in e['ks']),
                    'invalid actual signed native Q3 k map')
            require(xr > xl and rr > rl and xl >= -1e-12 and xr <= 1+1e-12 and
                    rl >= -1e-12 and rr <= self.m.R+1e-12,
                    'invalid native rectangular cell bounds')
            visited.update(ids.tolist()); self.volume += (xr-xl)*(rr-rl)
            x, r = np.meshgrid((xl+xr)/2+(xr-xl)*a/2, (rl+rr)/2+(rr-rl)*a/2)
            weights = (w[:, None]*w[None, :]*(xr-xl)*(rr-rl)/4).ravel()
            b = coupled_oracle.basis(e, 'k', x.ravel(), r.ravel())[0]
            self.mass[np.ix_(ids, ids)] += b.T@(weights[:, None]*b)
            for side, exterior in enumerate((abs(xl) <= 1e-12, abs(xr-1) <= 1e-12,
                                              abs(rl) <= 1e-12, abs(rr-self.m.R) <= 1e-12)):
                if exterior:
                    self.faces.append((e, side))
                    if side == 0:
                        essential.update(e['k'][j] for j in coupled_oracle.EDGES[side])
        require(visited == set(range(self.nk)) and len(set(self.free)) == len(self.free) and
                set(self.free) == set(range(self.nk))-essential, 'incomplete k map or wrong free mask')
        require(abs(self.volume-self.m.R) <= 1e-12, 'native cells do not cover the MMS area')
        self.riesz = la.cho_factor((self.m.Kin**2*self.mass)[np.ix_(self.free, self.free)], lower=True)

    def clock(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('boundary diagnostic reached its declared elapsed cap')

    def inner(self, a, b):
        return float(a@la.cho_solve(self.riesz, b))/self.volume

    def norm(self, a):
        return math.sqrt(max(0., self.inner(a, a)))

    def face_rule(self, e, side, order, split):
        a, w = leggauss(order); xl, xr, rl, rr = e['bounds']
        if side >= 2:
            return ((xl+xr)/2+(xr-xl)*a/2, np.full(order, rl if side == 2 else rr, dtype=float),
                    w*(xr-xl)/2, [xl, xr])
        xface = xl if side == 0 else xr
        candidates = []
        if split:
            K = float(self.m.functions['K'](xface, self.m.H, self.when)[0])
            candidates = [self.m.cutoff(K), .4*self.m.H, .75*self.m.H, self.m.h, self.m.H]
        # Identify an existing cell edge within arithmetic roundoff; never turn
        # an interior physical cutoff into a mesh-coordinate approximation.
        tolerance = 64*np.finfo(float).eps*max(1., abs(rl), abs(rr))
        cuts = sorted(set([rl, rr]+[v for v in candidates if rl+tolerance < v < rr-tolerance]))
        radii, weights = [], []
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            radii.extend((lo+hi)/2+(hi-lo)*a/2); weights.extend(w*(hi-lo)/2)
        return (np.full(len(radii), xface, dtype=float), np.asarray(radii, dtype=float),
                np.asarray(weights, dtype=float), cuts)

    def rhs_integrand(self, x, r, normal):
        m = self.m
        x, r = np.asarray(x, dtype=float), np.asarray(r, dtype=float)
        require(x.ndim == r.ndim == 1 and x.shape == r.shape and
                np.all(np.isfinite(x)) and np.all(np.isfinite(r)), 'invalid physical face coordinates')
        jets = m.functions['k'](x, r, self.when)
        k, kx, kr = [np.broadcast_to(np.asarray(v, dtype=float), x.shape) for v in jets[:3]]
        K = np.broadcast_to(np.asarray(m.functions['K'](x, r, self.when)[0], dtype=float), x.shape)
        outer = r >= m.H-1e-14*m.H
        # A JSON ground coordinate may be integer 0. Destruction must never
        # inherit that coordinate dtype: integer storage truncates finite Ew
        # to zero. This is floating storage, not a positivity floor.
        E = np.empty(r.shape, dtype=float)
        if np.any(outer):
            E[outer] = m.functions['epsilon'](x[outer], r[outer], self.when)[0]
        if np.any(~outer):
            c = np.asarray([m.cutoff(float(value)) for value in K[~outer]])
            # Ground and inner points use the DEPENDENT destruction, not a
            # dummy/extended transported outer-epsilon field. H keeps H+ owner.
            E[~outer] = K[~outer]**1.5/m.cl/np.maximum(r[~outer], c)
        require(np.all(np.isfinite(k)) and np.all(k > 0) and
                np.all(np.isfinite(E)) and np.all(E > 0),
                'nonpositive/nonfinite independent actual k or dependent/outer destruction before diffusion')
        D = m.nu+m.cd*k*k/(m.sk*E)
        require(np.all(np.isfinite(D)) and np.all(D > 0), 'invalid independent MMS scalar diffusivity')
        # j_out=-D grad(k).n, natural RHS=-integral(phi*j_out).
        return D*(normal[0]*kx+normal[1]*kr)

    def load(self, order, split):
        total = np.zeros(self.nk); contributions = []
        for e, side in self.faces:
            self.clock(); x, r, w, cuts = self.face_rule(e, side, order, split)
            b = coupled_oracle.basis(e, 'k', x, r)[0]
            local = b.T@(w*self.rhs_integrand(x, r, coupled_oracle.NORMALS[side]))
            ids = np.asarray(e['k']); total[ids] += local
            local_full = np.zeros(self.nk); local_full[ids] = local
            contributions.append(dict(geometry=e['geometry'], side=side,
                side_name=('left', 'right', 'ground', 'top')[side], bounds=e['bounds'],
                split_coordinates=cuts, coordinate='r' if side < 2 else 'x',
                essential_reaction_only=side == 0, local_ids=ids.tolist(), local_load=local.tolist(),
                free_dual_norm=self.norm(local_full[self.free])))
        return total, contributions

    def current_residual(self, nx, nr):
        """Independent spatial residual, using the recorded explicit RHS loads."""
        self.clock(); data = self.data; sizes = data['sizes']
        s = coupled_oracle.state(data['coefficients'], sizes)
        h = data['input_history']
        current = coupled_oracle.state(h['current'], sizes)
        previous = coupled_oracle.state(h['previous'], sizes)
        adv = coupled_oracle.state(h['advection'], sizes)
        old_adv = coupled_oracle.state(h['previous_advection'], sizes)
        require(coupled_oracle.state_hash(s) == data['state_sha256'] and
                coupled_oracle.state_hash(current) == data['previous_sha256'] and
                coupled_oracle.state_hash(adv) == data['previous_advection_sha256'],
                'current/history byte identity differs')
        require(data['step'] == h['step']+1 and data['dt'] == h['dt'] > 0 and
                abs(data['time']-data['step']*data['dt']) <= 1e-13,
                'accepted step/history timing differs')
        startup = h['step'] == 0; dt = h['dt']
        history = coupled_oracle.combine(current, (1. if startup else 2.)/dt,
                                         previous, 0. if startup else -.5/dt)
        explicit = coupled_oracle.combine(adv, 1. if startup else 2., old_adv, 0. if startup else -1.)
        o = coupled_oracle.Oracle(data, self.m, nx, nr, self.deadline)
        base = coupled_oracle.combine(o.forcing(s, data['time']), 1., o.mass(history), 1.)
        base = coupled_oracle.combine(base, 1., explicit, -1.)
        residual, norms = o.residual(s, data['time'], (1. if startup else 1.5)/dt, base)
        return residual['k'][self.free], norms


def audit(args, report, deadline):
    packet = json.loads(args.native.read_text()); table = json.loads(args.table.read_text())
    require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True,
            'input is not a registered native run')
    data = packet['state']
    require(data['schema'] == 'WJH27NativeCoupledStateV2' and table['status'] == 'PASS' and
            data['reference_identity'] == sha(args.table) and
            data['model'] == 'FiniteNuIsotropicWallFullWJCollarV2', 'wrong/unbound model or table')
    report.update(native_sha256=sha(args.native), table_sha256=sha(args.table),
                  native_source_identity=data['source_identity'], native_signature=data['signature'],
                  time=data['time'], step=data['step'], native_reported_residuals=data['residuals'])
    original_norm = args.original_k_norm
    if args.oracle:
        previous = json.loads(args.oracle.read_text())
        require(previous['native_sha256'] == report['native_sha256'] and
                previous['table_sha256'] == report['table_sha256'], 'supplied oracle belongs to a different run/table')
        original_norm = float(previous['independent_assemblies'][-1]['current_residual']['k'])
        report['original_oracle'] = dict(path=str(args.oracle), sha256=sha(args.oracle), status=previous['status'])
    report['original_current_k_residual_norm'] = original_norm
    b = BoundaryAudit(data, table, deadline)
    report['independent_MMS_constructor_route'] = b.mms_api
    old, old_faces = b.load(16, False)
    resolved24, faces24 = b.load(24, True)
    resolved32, faces32 = b.load(32, True)
    defect = (old-resolved32)[b.free]; refinement = (resolved24-resolved32)[b.free]
    defect_norm, refinement_norm = b.norm(defect), b.norm(refinement)
    report.update(k_free_ids=b.free.tolist(), domain_area=b.volume, k_reference=b.m.Kin,
        unresolved_Gauss16_load_free=old[b.free].tolist(),
        resolved_split_Gauss24_load_free=resolved24[b.free].tolist(),
        resolved_split_Gauss32_load_free=resolved32[b.free].tolist(),
        predicted_current_residual_defect_free=defect.tolist(),
        boundary_defect_k_dual_norm=defect_norm, resolved24_to32_k_dual_difference=refinement_norm,
        resolved_quadrature_diagnostic=refinement_norm <= 1e-12,
        norms='sqrt(r_free^T [Kin^2 M_k,free]^{-1} r_free / full domain area); no alpha normalization.',
        sign='At an accepted old state, the resolved-equation residual is predicted to contain L16-L32, since A*q=volume_load+L16.',
        faces=dict(unresolved16=old_faces, resolved24=faces24, resolved32=faces32),
        left_face_note='All native NaturalFlux callback faces are reproduced. Left essential rows contribute only reactions; their restriction to k-free tests is zero.')
    if original_norm is not None:
        require(math.isfinite(original_norm) and original_norm >= 0, 'invalid original residual norm')
        report['boundary_defect_to_original_norm_ratio'] = defect_norm/original_norm if original_norm > 0 else None
    report['current_vector_comparison'] = dict(status='NOT_PERFORMED', reason='Enable --compare-current for vector/direction evidence.')
    if args.compare_current:
        residual, norms = b.current_residual(args.nx, args.nr)
        residual_norm = b.norm(residual); remainder = residual-defect
        denominator = residual_norm*defect_norm
        cosine = b.inner(residual, defect)/denominator if denominator > 0 else None
        remainder_norm = b.norm(remainder)
        report['current_vector_comparison'] = dict(status='RECORDED_DIAGNOSTIC_ONLY',
            quadrature=[args.nx, args.nr], current_k_residual_free=residual.tolist(),
            current_k_residual_norm=residual_norm, current_all_block_norms=norms,
            predicted_defect_riesz_cosine=cosine,
            residual_minus_predicted_defect_free=remainder.tolist(),
            residual_minus_predicted_defect_norm=remainder_norm,
            reason='Spatial operators/continuous weak forcing independently reconstructed; recorded BE/BDF2-EX2 history loads are used unchanged, not requalified here.',
            caution='Removing the defect is a causal diagnostic, not a solved corrected system or permission to override any residual gate.')
        report['original_current_k_residual_norm'] = residual_norm if original_norm is None else original_norm
        report['causal_boundary_explanation_consistent'] = bool(
            refinement_norm <= 1e-12 and cosine is not None and cosine >= .99 and remainder_norm <= 1e-10)
        report['causal_diagnostic_criteria'] = dict(resolved24_to32_dual_max=1e-12,
            minimum_riesz_cosine=.99, residual_minus_predicted_defect_dual_max=1e-10,
            note='Diagnostic agreement criteria only; no native residual tolerance or acceptance gate is changed.')
    report['status'] = 'BOUNDARY_SPLIT_DIAGNOSTIC_RECORDED'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--oracle', type=Path)
    parser.add_argument('--original-k-norm', type=float)
    parser.add_argument('--compare-current', action='store_true')
    parser.add_argument('--nx', type=int, default=24)
    parser.add_argument('--nr', type=int, default=32)
    parser.add_argument('--seconds', type=float, default=180.)
    args = parser.parse_args()
    if args.out.exists(): parser.error('immutable diagnostic output already exists')
    started = time.monotonic()
    report = dict(schema='WJH27BoundarySplitAuditV1', status='FAIL_BOUNDARY_SPLIT_AUDIT',
        full_coupled_gate_pass=False, gate_override=False, platform=sys.platform,
        source_sha256=sha(__file__), independent_coupled_source_sha256=sha(coupled_oracle.__file__),
        independent_MMS_source_sha256=sha(coupled_mms_check.__file__),
        independent_basis_source_sha256=sha(Path(__file__).resolve().parents[1]/'h26'/'fixtures.py'),
        native=str(args.native), table=str(args.table),
        scope='Prescribed MMS k boundary-load quadrature diagnosis only. No corrected native state, new tolerance or gate pass.')
    failure = None
    try:
        require(sys.platform.startswith('linux'), 'authoritative boundary diagnosis requires Linux')
        require(0 < args.seconds <= 300 and args.nx >= 20 and args.nr >= 24,
                'invalid diagnostic time/quadrature bounds')
        audit(args, report, started+args.seconds)
    except Exception as error:
        failure = error; report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if failure is not None:
        raise failure
    print(json.dumps({key: report.get(key) for key in ('status', 'boundary_defect_k_dual_norm',
        'resolved24_to32_k_dual_difference', 'causal_boundary_explanation_consistent', 'elapsed_seconds')}, sort_keys=True))


if __name__ == '__main__':
    main()
