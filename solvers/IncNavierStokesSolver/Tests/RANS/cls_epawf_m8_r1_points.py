#!/usr/bin/env python3
"""Later archive-only export of every negative inner covariance evaluation.

The original phase-A reconstruction and ledger remain immutable inputs.  This
helper adds no sampling sets, changes no closure and launches no solver.  It
streams all negative entries in the three already declared inner constitutive
sets, including repeated coordinates when the original finite rules repeat
them.  The full NDJSON remains outside the checkout; the compact receipt is
intended for source control.  Evaluate only on the designated Linux host.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

sys.dont_write_bytecode = True
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
os.environ.update({key: '1' for key in THREAD_KEYS})
import numpy as np

BASELINE = 'ab5c3c6fee19761e8d88692e8e6c331e5889234e'
RECONSTRUCTION_SHA = '0c2700ec408ae12d72bb1442eb978b9ce416be9fd86a85978534a06aa105274e'
LEDGER_SHA = 'ed879358611e31a1d8f90fca30d447bf8f5dccf60c6fddb1d24840f7d5268a26'
TARGETS = ('source_volume_inner', 'physical_sample_K', 'interface_inner')


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source),
                                    *arguments], text=True).strip()


def load_driver(path):
    spec = importlib.util.spec_from_file_location('m8_r1_frozen_archive_reconstruction', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def oriented(vector):
    vector = np.asarray(vector).copy()
    if vector[np.argmax(np.abs(vector))] < 0:
        vector = -vector
    return vector.tolist()


def make_engine(module, raw, source, xml, stream):
    class Export(module.Reconstruction):
        def __init__(self):
            super().__init__(raw, source, xml)
            self.initial_ray_cache = {}
            self.record_count = 0
            self.per_set_totals = {name: 0 for name in TARGETS}
            self.first_exported = self.minimum_exported = None

        def begin_epoch(self, step, at_time):
            self.export_step, self.export_time = step, float(at_time)
            self.calls = {name: 0 for name in TARGETS}
            self.export_counts = {name: 0 for name in TARGETS}
            self.export_minima = {name: None for name in TARGETS}
            self.epoch_minimum = None
            self.epoch_row_start = self.record_count

        def snapshot(self, points, ray, side, eid):
            """Same indexed AWF tensor and physical field rows as phase A."""
            k = self.evaluate(self.d, points, side, eid)[0]
            G = self.velocity(points, side, eid)[1]
            S = G+G.transpose(0, 2, 1)
            p = ray['profile']
            nua = float(p.a)*np.maximum(points[:, 2]-float(p.zv), 0.)
            R = 2*k[:, None, None]/3*np.eye(3)-nua[:, None, None]*S
            eigenvalues, vectors = np.linalg.eigh(R)
            strain_eigenvalues = np.linalg.eigvalsh(S)
            epsilon = ray['lam']*np.asarray(p.fields(points[:, 2]), float)[2]
            return dict(k=k, G=G, S=S, R=R, eigenvalues=eigenvalues, vectors=vectors,
                strain_eigenvalues=strain_eigenvalues, nua=nua, epsilon=epsilon,
                sampled_K=float(p.k), chi=float(p.chi), b=float(ray['b']),
                lambda_value=float(ray['lam']))

        def initial_snapshot(self, points, side, eid):
            current_state, current_rays = self.state, self.rays
            try:
                self.state = np.asarray(self.raw['initial'], float)
                self.rays = self.initial_ray_cache
                # Each original target call belongs to one wall ray.  Refuse
                # to silently assign a different profile to a mixed-ray batch.
                if not np.all(points[:, :2] == points[0, :2]):
                    raise ValueError('The frozen target call unexpectedly contains multiple wall rays')
                ray = self.ray(points[0])
                result = self.snapshot(points, ray, side, eid)
                seed = self.ray(points[0], seed=True)['profile']
                neutral_nua = float(seed.a)*np.maximum(points[:, 2]-float(seed.zv), 0.)
                result['analytic_neutral_shear'] = self.us**2/(self.nu+neutral_nua)
                return result
            finally:
                self.state, self.rays = current_state, current_rays

        def point_record(self, values, index):
            k, K = float(values['k'][index]), values['sampled_K']
            largest = float(values['strain_eigenvalues'][index, -1])
            return dict(k=k, sampled_K=K, k_over_K=k/K,
                nu_A=float(values['nua'][index]), epsilon_profile=float(values['epsilon'][index]),
                G=values['G'][index].tolist(), S=values['S'][index].tolist(), R=values['R'][index].tolist(),
                eigenvalue=float(values['eigenvalues'][index, 0]),
                eigenvalues=values['eigenvalues'][index].tolist(),
                eigenvector=oriented(values['vectors'][index, :, 0]),
                maximum_strain_eigenvalue=largest,
                strain_norm=float(np.linalg.norm(values['S'][index])),
                r_R=float(3*values['nua'][index]*largest/(2*k)) if k > 0 else None,
                r_R_status='finite_positive_local_k' if k > 0 else 'undefined_nonpositive_local_k',
                chi=values['chi'], b=values['b'], lambda_value=values['lambda_value'])

        def sample(self, name, points, inner, constitutive=False, ray=None, side=0, eid=None):
            # This is the complete original calculation, including its finite
            # set summaries and checks.  The added export does not replace it.
            super().sample(name, points, inner, constitutive, ray, side, eid)
            if name not in TARGETS:
                return
            if not inner or not constitutive or ray is None:
                raise ValueError('Frozen target set no longer has its declared constitutive meaning')
            points = np.atleast_2d(points)
            self.calls[name] += 1
            current = self.snapshot(points, ray, side, eid)
            negative = np.flatnonzero(current['eigenvalues'][:, 0] < 0)
            if not len(negative):
                return
            initial = self.initial_snapshot(points, side, eid)
            if self.export_step == 0:
                for key in ('k', 'G', 'S', 'R', 'nua', 'epsilon'):
                    self.near(current[key], initial[key], 'all_initial_samepoint_'+key, 2e-13)
            for index in negative:
                now = self.point_record(current, index)
                before = self.point_record(initial, index)
                before['analytic_neutral_shear'] = float(initial['analytic_neutral_shear'][index])
                point = points[index]
                record = dict(schema=1, record_index=self.record_count,
                    accepted_step=self.export_step, time=self.export_time, set=name,
                    original_sample_call=self.calls[name]-1, original_point_index=int(index),
                    coordinate=point.tolist(), side=int(side), requested_element=eid,
                    velocity_element=int(self.fields[0].element(point, side)['geometry']) if eid is None else int(eid),
                    k_element=int(self.fields[self.d].element(point, side)['geometry']) if eid is None else int(eid),
                    **now, initial_samepoint=before,
                    change_from_initial=dict(k=now['k']-before['k'], sampled_K=now['sampled_K']-before['sampled_K'],
                        k_over_K=now['k_over_K']-before['k_over_K'], nu_A=now['nu_A']-before['nu_A'],
                        strain_norm=now['strain_norm']-before['strain_norm'],
                        covariance_minimum=now['eigenvalue']-before['eigenvalue']),
                    interpretation='An archived original inner AWF covariance failure at a declared finite evaluation entry; no new closure, trajectory or causal attribution.')
                stream.write((json.dumps(record, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode())
                self.record_count += 1
                self.export_counts[name] += 1
                self.per_set_totals[name] += 1
                if self.first_exported is None:
                    self.first_exported = record
                if self.minimum_exported is None or record['eigenvalue'] < self.minimum_exported['eigenvalue']:
                    self.minimum_exported = record
                if self.epoch_minimum is None or record['eigenvalue'] < self.epoch_minimum['eigenvalue']:
                    self.epoch_minimum = record
                minimum = self.export_minima[name]
                self.export_minima[name] = record['eigenvalue'] if minimum is None else min(minimum, record['eigenvalue'])
    return Export()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--ledger', required=True, type=Path)
    parser.add_argument('--reconstruction-driver', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux' or set(os.sched_getaffinity(0)) != {50}:
        raise SystemExit('Run the archive-only supplement on Linux with taskset -c 50')
    args.output.mkdir(parents=True, exist_ok=False)
    started, start_utc = time.monotonic(), stamp()
    report = dict(schema=1, phase='M8-R1 A later all-offending-points supplement', started_utc=start_utc,
        success=False, trajectory_launched=False, solver_or_build_launched=False,
        numerical_workers=1, affinity=sorted(os.sched_getaffinity(0)),
        numerical_threads={key: os.environ[key] for key in THREAD_KEYS}, command=sys.argv,
        helper_sha256=sha(__file__), source=str(args.source), archive=str(args.archive),
        states=[], failures=[], completed_epoch_count=0,
        scope='Every negative covariance evaluation entry in the three existing inner constitutive finite sets, not a new sampling or polynomial-everywhere theorem.',
        chronology='This later archive-only export supplements the original completed phase-A gate; it does not claim to have been collected before the subsequent diagnostic trajectories or rewrite their timing.')
    summary_path, ndjson_path = args.output/'points-summary.json', args.output/'negative-inner-points.ndjson'
    engine = None
    module = None
    try:
        if git(args.source, 'rev-parse', 'HEAD') != BASELINE or git(args.source, 'status', '--porcelain'):
            raise ValueError('Imported numerical oracles require the clean immutable ab5 archive source')
        if sha(args.reconstruction_driver) != RECONSTRUCTION_SHA:
            raise ValueError('The completed phase-A reconstruction driver bytes changed')
        if sha(args.ledger) != LEDGER_SHA:
            raise ValueError('The original complete Linux phase-A ledger bytes changed')
        report['source_commit'] = BASELINE
        report['source_clean_start'] = True
        report['inputs'] = {str(path): dict(bytes=path.stat().st_size, sha256=sha(path)) for path in (
            args.reconstruction_driver, args.ledger, args.archive/'case.evolution.json',
            args.archive/'case.xml', args.archive/'run.log')}
        ledger = json.loads(args.ledger.read_text())
        raw = json.loads((args.archive/'case.evolution.json').read_text())
        xml = (args.archive/'case.xml').read_text()
        report['original_phase_A'] = dict(ledger=str(args.ledger), sha256=LEDGER_SHA,
            driver_sha256=ledger['driver_sha256'], finished_utc=ledger['finished_utc'],
            success=ledger['success'], source_commit=ledger['source_commit'],
            historical_first_negative_inner_event=next(e for e in ledger['events'] if e['type'] == 'first_negative_inner_covariance'))
        if not ledger['success'] or ledger['driver_sha256'] != RECONSTRUCTION_SHA or ledger['source_commit'] != BASELINE:
            raise ValueError('The original phase-A gate is not the required successful frozen receipt')
        if datetime.fromisoformat(start_utc) <= datetime.fromisoformat(ledger['finished_utc']):
            raise ValueError('Supplemental export must be recorded after the original phase-A gate')
        if not (len(raw['accepted_states']) == len(raw['audit_times']) == len(raw['state_audits']) == len(ledger['states']) == 12):
            raise ValueError('Expected every one of the 12 archived accepted epochs')
        module = load_driver(args.reconstruction_driver)
        with ndjson_path.open('xb') as stream:
            engine = make_engine(module, raw, args.source, xml, stream)
            differential = raw['offsets'][raw['dimension']+2]
            for step, (state, at_time, native, old) in enumerate(zip(
                    raw['accepted_states'], raw['audit_times'], raw['state_audits'], ledger['states'])):
                if step == 0:
                    history, alpha = list(state), 0.
                else:
                    history = np.asarray(raw['accepted_states'][step-1], float).copy()
                    alpha = 1/raw['dt'] if step == 1 else 1.5/raw['dt']
                    if step > 1:
                        history[:differential] = (4*np.asarray(raw['accepted_states'][step-1])[:differential]-
                                                np.asarray(raw['accepted_states'][step-2])[:differential])/3
                engine.begin_epoch(step, at_time)
                reconstructed = engine.state_record(state, history, alpha, at_time, native)
                if reconstructed['state_sha256'] != old['state_sha256'] or reconstructed['history_sha256'] != old['history_sha256']:
                    raise ValueError('Original state/history identities changed')
                engine.near(at_time, old['time'], 'original_epoch_time', 0.)
                sets = {}
                for name in TARGETS:
                    current, previous = reconstructed['sets'][name], old['sets'][name]
                    if engine.export_counts[name] != current['negative_covariance_count'] or engine.export_counts[name] != previous['negative_covariance_count']:
                        raise ValueError('Incomplete negative-point export for '+name)
                    if current['count'] != previous['count']:
                        raise ValueError('The original declared finite population changed for '+name)
                    engine.near(current['minimum_covariance'], previous['minimum_covariance'], 'original_minimum_'+name, 1e-11)
                    if engine.export_minima[name] is not None:
                        engine.near(engine.export_minima[name], previous['minimum_covariance'], 'all_exported_minimum_'+name, 1e-11)
                    sets[name] = dict(declared_evaluation_count=current['count'],
                        negative_entries_exported=engine.export_counts[name],
                        original_negative_count=previous['negative_covariance_count'],
                        minimum_covariance=current['minimum_covariance'],
                        minimum_k=current['minimum_k'], invalid_count=current['invalid_count'],
                        exact_original_count_and_complete_export=True)
                report['states'].append(dict(accepted_step=step, time=float(at_time),
                    state_sha256=reconstructed['state_sha256'], history_sha256=reconstructed['history_sha256'],
                    sets=sets, minimum_exported_record=engine.epoch_minimum,
                    ndjson_record_range_half_open=[engine.epoch_row_start, engine.record_count]))
                report['completed_epoch_count'] += 1
                stream.flush()
                module.save(summary_path, report)
                print(json.dumps(dict(accepted_step=step, time=at_time,
                    negative_entries_by_set=engine.export_counts, total_written=engine.record_count)), flush=True)
        report['points_file'] = dict(path=str(ndjson_path), bytes=ndjson_path.stat().st_size,
            sha256=sha(ndjson_path), records=engine.record_count,
            format='UTF-8 compact JSON objects, one per line; exact file-byte SHA256',
            identity='(accepted_step,set,original_sample_call,original_point_index); duplicates in coordinates remain distinct original finite evaluation entries')
        report['negative_entries_by_set'] = engine.per_set_totals
        report['first_exported_record'] = engine.first_exported
        report['global_minimum_record'] = engine.minimum_exported
        report['initial_state_summary'] = report['states'][0]
        report['checks'], report['errors'] = engine.checks, engine.errors
        report['interpretation'] = dict(every_declared_negative_entry_exported=True,
            original_sample_sets_unchanged=True, source_cap_causality_inferred=False,
            original_gate_timing_unchanged=True, physical_model_qualified=False,
            initial_comparison='Each exported current point includes the actual original initial coefficient state, initial physical sampler K/nuA and full tensor at the identical location/side, plus an explicitly separate analytic neutral-shear control.')
        for path, receipt in report['inputs'].items():
            if sha(path) != receipt['sha256']:
                raise ValueError('An immutable input changed: '+path)
        report['source_commit_end'] = git(args.source, 'rev-parse', 'HEAD')
        report['source_clean_end'] = not bool(git(args.source, 'status', '--porcelain'))
        if report['source_commit_end'] != BASELINE or not report['source_clean_end']:
            raise ValueError('Pinned numerical source changed during point export')
        report['imported_source_files'] = {str(Path(m.__file__).resolve()): sha(m.__file__) for m in
            list(sys.modules.values()) if getattr(m, '__file__', None) and
            str(Path(m.__file__).resolve()).startswith(str(args.source.resolve()))}
        report['success'] = True
        report['finished_utc'], report['seconds'] = stamp(), time.monotonic()-started
        module.save(summary_path, report)
        print(json.dumps(dict(success=True, epochs=report['completed_epoch_count'],
            records=engine.record_count, per_set=engine.per_set_totals, output=str(summary_path))), flush=True)
        return 0
    except Exception as error:
        report['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
        report['finished_utc'], report['seconds'] = stamp(), time.monotonic()-started
        if ndjson_path.exists():
            report['partial_points_file'] = dict(path=str(ndjson_path), bytes=ndjson_path.stat().st_size,
                sha256=sha(ndjson_path), complete=False,
                records_written=engine.record_count if engine is not None else None)
        summary_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n')
        print(traceback.format_exc(), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
