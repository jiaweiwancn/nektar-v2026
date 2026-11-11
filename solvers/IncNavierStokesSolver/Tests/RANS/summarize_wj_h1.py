#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only compact H1 publication; inputs remain at their original paths."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numpy.testing import assert_allclose
from scipy.integrate import quad
from validate_wj1 import QE, KAPPA
from wj_h1_cases import NU, US, HW


def load(path):
    # Native JSON may spell signed zero as -0. Preserve its binary64 sign when
    # replaying hashes; other integer metadata retain their natural type.
    return json.loads(Path(path).read_text(), parse_int=lambda s: -0. if s == '-0' else int(s))


def identity(path):
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def unique_run(build, name):
    runs = list((build / ('wj-h1-' + name)).glob('run-*'))
    assert len(runs) == 1, 'ambiguous run; select an immutable build/output epoch'
    assert load(runs[0] / 'result.json')['status'] == 'PASS'
    return runs[0]


def csv_file(path, rows):
    assert rows
    with Path(path).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def timing(path):
    records = [json.loads(s) for s in Path(path).read_text().splitlines() if s.strip()]
    by_id = {(v['pid'], v['id']): v for v in records}
    sums = {}
    for v in records:
        p = v
        while p['parent_id']:
            p = by_id[(p['pid'], p['parent_id'])]
        # Exclude initialization and passive assembly from governing-step cost.
        if p['name'] != 'pimple.sweep':
            continue
        key = v['name']
        entry = sums.setdefault(key, dict(calls=0, inclusive_seconds=0., exclusive_seconds=0.))
        entry['calls'] += 1
        entry['inclusive_seconds'] += v['wall_seconds']
        entry['exclusive_seconds'] += v['exclusive_wall_seconds']
    return dict(scope='descendants of actual pimple.sweep only; inclusive rows overlap', terms=sums)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    args.output.mkdir(parents=True, exist_ok=False)
    heavy, summaries, interfaces, cases, fields = [], {}, [], [], []
    geometry = unique_run(args.build, 'geometry')
    save(args.output / 'geometry.json', load(geometry / 'result.json'))
    # A second oracle for the signed-forcing fixture integrates from its
    # independent inputs, without using the returned J0/J1 as expectations.
    packet = load(geometry / 'geometry.json')
    h, K = .090675, QE * US ** 2  # exact C++ signed-fixture inputs
    a = KAPPA * US
    x = h * np.sqrt(K) / (70 * NU)
    b = 10.7 * (1 - x ** max(.5 - .4 * x ** .7, 1 - .79 * x ** (-.28)))
    zv = NU * b / np.sqrt(K)
    D = lambda r: NU + a * max(r - zv, 0)
    j0 = sum(quad(lambda r: 1 / D(r), l, u, epsabs=1e-12, epsrel=2e-13)[0]
             for l, u in zip([0, h], [h, HW]))
    j1 = quad(lambda r: (r - h) / D(r), h, HW, epsabs=1e-13, epsrel=2e-13)[0]
    signed = []
    for row in packet['signed_forcing']:
        traction = (row['U_H'] - row['F'] * j1) / j0
        production = sum(quad(lambda r: (D(r) - NU) *
            (traction + row['F'] * max(r - h, 0)) ** 2 / D(r) ** 2,
            l, u, epsabs=1e-14, epsrel=2e-13)[0] for l, u in zip([0, h], [h, HW])) / HW
        assert_allclose([row['J0'], row['J1'], row['traction'], row['mean_production']],
                        [j0, j1, traction, production], rtol=2e-12, atol=2e-13)
        signed.append(dict(input_U_H=row['U_H'], input_F=row['F'], independent_traction=traction,
                           independent_mean_production=production, J0=j0, J1=j1))
    save(args.output / 'signed-quadrature.json', signed)
    heavy += [identity(p) for p in geometry.iterdir() if p.is_file()]
    for name in ('hill', 'zero', 'plane'):
        folder = unique_run(args.build, name)
        report_path = folder / 'report.json'
        report, result = load(report_path), load(folder / 'result.json')
        for callback in report['callbacks']:
            for field_name in ('state', 'history'):
                actual = hashlib.sha256(np.asarray(callback[field_name], dtype='<f8').tobytes()).hexdigest()
                assert actual == callback[field_name + '_sha256'], 'callback payload/identity mismatch'
        assert np.array_equal(report['state'], report['callbacks'][-1]['state'])
        for attempt in report['pimple_attempts']:
            snapshot = load(attempt['k_snapshot'])
            assert np.asarray(snapshot['raw'], dtype='<f8').tobytes() == np.asarray(snapshot['corrected'], dtype='<f8').tobytes()
        save(args.output / (name + '-independent.json'), result)
        heavy += [identity(p) for p in folder.rglob('*') if p.is_file()]
        packets = [dict(step=0, time=0, interface=report['startup_interface'])] + report['events']
        for packet in packets:
            for row in packet['interface']['samples']:
                interfaces.append(dict(case=name, step=packet['step'], time=packet['time'],
                    s=row['chart'][0], x=row['physical'][0], z=row['physical'][2],
                    inner_geometry=row['inner_geometry'], outer_geometry=row['outer_geometry'],
                    K=row['K'], epsilon=row['epsilon'], epsilon_target=row['epsilon_target'],
                    cutoff=row['cutoff'], gate=row['cutoff_gate'], stretch=row['metric_stretch'],
                    mass_sum=row['mass_inner'] + row['mass_outer'],
                    traction_x=row['traction_defect'][0], traction_z=row['traction_defect'][2],
                    k_flux_sum=row['k_flux_inner'] + row['k_flux_outer'],
                    mechanical_jump=row['mechanical_jump'], jump_target=row['normal_traction_jump_target'],
                    inner_covariance_min=float(np.linalg.eigvalsh(row['R_inner'])[0]),
                    outer_covariance_min=float(np.linalg.eigvalsh(row['R_outer'])[0]),
                    production_inner=row['production_inner'], production_outer=row['production_outer'],
                    nu_over_nut_H=row['nu_over_nut_H']))
        audits = [report['startup_audit']] + [v['audit'] for v in report['callbacks']]
        summary = dict(report=identity(report_path), design=load(folder / 'design.json'),
            signature=report['signature'], map=report['map'], state_size=report['state_size'],
            field_sizes=report['field_sizes'], volume=report['volume'],
            steps=report['steps'], time=report['time'], dt=report['dt'],
            initialization_work=report['initialization_work'],
            startup_assembly_checks=report['startup_assembly_checks'],
            final_assembly_checks=report['final_assembly_checks'],
            audits=audits, callback_bindings=[{k: v[k] for k in
                ('state_sha256', 'history_sha256', 'time', 'alpha', 'audit_packet')} for v in report['callbacks']],
            events=[dict(step=e['step'], time=e['time'], physical_observations=e['physical_observations'])
                for e in report['events']],
            checkpoints=report['checkpoints'], pimple_attempts=report['pimple_attempts'],
            roughness_budget=report['roughness_budget'], energy_budget=report['energy_budget'],
            seconds=report['seconds'], scope=report['scope'])
        for a in summary['pimple_attempts']:
            # Large inventory branch strings already remain in the original report.
            for domain in ('input_domain', 'flow_domain', 'final_domain'):
                if isinstance(a.get(domain), dict):
                    a[domain].pop('profile_branch_pattern', None)
        summary['timing'] = timing(args.root / 'operations' / ('native-' + name + '-v3') / 'native-timing.jsonl')
        assert summary['timing']['terms']['pimple.sweep']['calls'] == report['steps']
        summary['weak_fields'] = {}
        for phase, path in [('startup', folder / 'report.json.startup-weak.json'), ('final', Path(report['final_weak']))]:
            weak = load(path)
            regions = []
            for inner in (True, False):
                rows = [r for r in weak['samples'] if r['inner'] == inner]
                weights = np.array([r['weight'] for r in rows])
                divergence = np.array([np.trace(r['G']) for r in rows])
                maximum = int(np.argmax(abs(divergence)))
                regions.append(dict(region='inner' if inner else 'outer', measure=float(sum(weights)),
                    divergence_L2=float(np.sqrt(weights @ (divergence ** 2))),
                    maximum_divergence=float(abs(divergence[maximum])), maximum_point=rows[maximum]['physical'],
                    minimum_k=min(r['k'] for r in rows), maximum_k=max(r['k'] for r in rows),
                    minimum_epsilon=min(r['epsilon'] for r in rows), maximum_epsilon=max(r['epsilon'] for r in rows)))
                if name == 'hill':
                    # Same native physical points, downsampled for plotting only.
                    # A prime stride avoids aliasing the 24/32 point rules.
                    for row in rows[::29]:
                        fields.append(dict(phase=phase, region='inner' if inner else 'outer',
                            s=row['chart'][0], rho=row['chart'][2], x=row['physical'][0], z=row['physical'][2],
                            u=row['U'][0], w=row['U'][2], k=row['k'], epsilon=row['epsilon']))
            summary['weak_fields'][phase] = dict(input=identity(path), state_sha256=weak['state_sha256'],
                alpha_mass=weak['alpha_mass'], regions=regions)
            del weak
        summaries[name] = summary
        a = report['pimple_attempts'][0]
        cases.append(dict(case=name, state=report['state_size'], free_velocity=a['velocity_size'],
            pressure_gauge=a['pressure_gauge_size'], flow=a['assembled_flow_size'],
            k=a['k_solve_size'], epsilon=a['epsilon_solve_size'], steps=report['steps'],
            dt=report['dt'], time=report['time'], native_seconds=report['seconds'],
            independent_component_error=max(result['final_independent_loads']['component_max_errors'].values()),
            independent_regional_error=result['final_independent_loads']['regional_max_error'],
            K_change=result['maximum_K_change'], cutoff_change=result['maximum_cutoff_change'],
            epsilon_trace_change=result['maximum_epsilon_trace_change'],
            raw_work_defect=max(a['wj_work_defect_raw_divergence'][0] for a in audits),
            correction_activations=sum(a['bound_activations'] for a in report['pimple_attempts'])))
    save(args.output / 'cases.json', summaries)
    csv_file(args.output / 'cases.csv', cases)
    csv_file(args.output / 'interfaces.csv', interfaces)
    csv_file(args.output / 'hill-field-samples.csv', fields)
    # Include exact command receipts and all pre-scientific failed epochs.
    operations = {}
    for p in sorted((args.root / 'operations').rglob('*')):
        if p.is_file() and p.suffix in ('.json', '.log', '.jsonl', '.py'):
            heavy.append(identity(p))
            if p.name == 'receipt.json' or p.parent == args.root / 'operations' and p.suffix == '.json':
                operations[str(p.relative_to(args.root))] = load(p)
    save(args.output / 'command-receipts.json', operations)
    save(args.output / 'ledger.json', load(args.root / 'checks-accounting.json'))
    for p in (args.root / 'runtime-v3-prelaunch.json', args.root / 'runtime-v3-postbuild.json', args.root / 'source-v1.json',
              args.root / 'source-v2.json', args.root / 'source-v3.json'):
        save(args.output / p.name, load(p)); heavy.append(identity(p))
    save(args.output / 'heavy-index.json', heavy)
    files = [identity(p) for p in args.output.iterdir() if p.is_file()]
    save(args.output / 'manifest.json', dict(schema='WJH1CompactReadOnlyV1', program=identity(__file__),
        root=str(args.root.resolve()), build=str(args.build.resolve()),
        seconds=time.monotonic() - started, compact_files=files,
        verification='Remote files were read and hashed on Linux. Mac transport must verify compact bytes separately.',
        limitations='Controlled coupled prototype; no physical inlet/outlet calibration, convergence or stationarity.'))


if __name__ == '__main__':
    main()
