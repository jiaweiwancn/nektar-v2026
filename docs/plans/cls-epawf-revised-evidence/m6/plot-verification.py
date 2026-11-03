#!/usr/bin/env python3
"""Create four scientific panels from a passed native M6 qualification.

Only recorded numerical errors are plotted. The dt-squared curve is a labelled
reference guide, not a fitted rate. The omitted-grad-chi quantity is the
qualified norm of the missing analytic term, not another solved CFD case.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.update({name: '1' for name in
                   ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')})
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,NullFormatter


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def positive(value, name):
    value = float(value)
    require(np.isfinite(value) and value > 0, 'Log plot requires a positive measured '+name)
    return value


def operator_report(qualification, data):
    """Prefer the actual feet report; locate it through the saved command."""
    candidates = []
    for label in ('operator-feet-d2', 'operator-flat-d2'):
        candidates.append(qualification.parent/label/'probe.json')
        for command in data.get('commands', []):
            if command.get('name') == label:
                arguments = command.get('command', command.get('argv', []))
                for argument in reversed(arguments):
                    path = Path(str(argument))
                    if path.suffix == '.json':
                        candidates.append(path if path.is_absolute() else qualification.parent/label/path)
    receipts = data['sections'].get('operators', [])
    for path in candidates:
        if not path.is_file():
            continue
        raw = json.loads(path.read_text())
        if 'transition' not in raw:
            continue
        samples = raw['transition']['samples']
        actual = np.array([[sample['input'][0], sample['gates'][0], sample['gates'][-1]] for sample in samples])
        for receipt in receipts:
            records = receipt.get('physical_samples', receipt.get('profile_records', []))
            if len(records) != len(samples):
                continue
            qualified = np.array([[record['coordinate'], record['chi'], record['d2']] for record in records])
            if actual.shape == qualified.shape and np.allclose(actual, qualified, atol=1e-10, rtol=1e-10):
                require(raw.get('eta') == 0, 'Operator gate panel must use the unanchored M6 case')
                return raw, path.resolve(), receipt
    raise RuntimeError('Cannot locate native operator samples matching the passed qualification receipt')


def h_inset(axis, base, refined, key, ylabel, omission=False):
    inset = axis.inset_axes([.59, .55, .37, .37])
    h = np.array([1., 1/float(refined['refine'])])
    values = [positive(base['error'][key], key), positive(refined['error'][key], key)]
    inset.semilogy(h, values, 's-', color='#2463a0', markersize=4, linewidth=1.1)
    if omission:
        control = [positive(base['error']['omitted_grad_chi_rms'], 'omission control'),
                   positive(refined['error']['omitted_grad_chi_rms'], 'omission control')]
        inset.semilogy(h, control, ':', color='#946c27', linewidth=1.1)
    inset.set(title='p=3: h refinement', xlabel=r'$h/h_0$')
    inset.set_xticks(h)
    inset.set_xticklabels(['1', f"1/{refined['refine']}"])
    inset.set_ylabel(ylabel, fontsize=7)
    inset.tick_params(labelsize=7)
    inset.xaxis.label.set_size(8)
    inset.title.set_size(8)
    inset.title.set_bbox(dict(facecolor='white',edgecolor='none',pad=1))
    inset.grid(True, which='both', alpha=.18)
    return dict(relative_h=h.tolist(), values=values,
                omitted_grad_chi=[base['error']['omitted_grad_chi_rms'], refined['error']['omitted_grad_chi_rms']] if omission else None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qualification', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    qualification = args.qualification.resolve()
    data = json.loads(qualification.read_text())
    status = data.get('status', data.get('verification_status'))
    require(status == 'passed' and data.get('success', True), 'Qualification has not passed')
    sections = data['sections']
    manufactured = sections.get('manufactured', [])
    temporal = sections.get('temporal', [])
    require(bool(manufactured) and len(temporal) >= 3 and 'h_refinement' in sections,
            'Passed receipt lacks the required p/h/time studies')
    raw, raw_path, operator_receipt = operator_report(qualification, data)
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9,
                         'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    figure, axes = plt.subplots(2, 2, figsize=(12.2, 8.6), constrained_layout=True)
    colors = {'flat': '#2463a0', 'feet': '#a95627'}
    markers = {'flat': 'o', 'feet': '^'}
    epsilon_series, trace_series = [], []
    families = sorted({record['family'] for record in manufactured})
    for family in families:
        records = sorted((record for record in manufactured if record['family'] == family), key=lambda record: record['degree'])
        require(len(records) >= 2, 'Insufficient native p refinement for '+family)
        degrees = [int(record['degree']) for record in records]
        epsilon = [positive(record['error']['epsilon'], 'epsilon RMS') for record in records]
        trace = [positive(record['error']['trace_gradient_rms'], 'trace gradient RMS') for record in records]
        omission = [positive(record['error']['omitted_grad_chi_rms'], 'omitted grad chi RMS') for record in records]
        color, marker = colors.get(family), markers.get(family, 's')
        label = 'flat wall' if family == 'flat' else 'C1 hill feet' if family == 'feet' else family
        style=dict(markerfacecolor='none',markersize=8) if family=='feet' else dict(markersize=5)
        axes[0, 0].semilogy(degrees, epsilon, marker+'-', color=color, label=label,**style)
        axes[0, 1].semilogy(degrees, trace, marker+'-', color=color, label=label+' error',**style)
        axes[0, 1].semilogy(degrees, omission, ':', color=color, linewidth=1.2,
                           label=label+' analytic omission norm')
        epsilon_series.append(dict(family=family, degree=degrees, physical_epsilon_rms=epsilon))
        trace_series.append(dict(family=family, degree=degrees, physical_trace_gradient_rms=trace,
                                 omitted_grad_chi_rms=omission))
    base = next((record for record in manufactured if record['family'] == 'flat' and record['degree'] == 3), None)
    require(base is not None and sections['h_refinement']['refine'] > 1, 'Missing flat p3 h-refinement baseline')
    epsilon_h = h_inset(axes[0, 0], base, sections['h_refinement'], 'epsilon', 'RMS error')
    trace_h = h_inset(axes[0, 1], base, sections['h_refinement'], 'trace_gradient_rms', 'RMS error', True)
    axes[0, 0].set(title='(a) Native MMS: dissipation in physical space',
                   xlabel='Velocity/scalar degree p (base mesh)', ylabel=r'RMS $\epsilon$ error [$\mathrm{m^2\,s^{-3}}$]')
    axes[0, 1].set(title=r'(b) Matching trace: tangential derivative and $\nabla\chi$',
                   xlabel='Velocity/scalar degree p (base mesh)', ylabel=r'RMS $\partial_x\epsilon(H)$ error [$\mathrm{m\,s^{-3}}$]')
    degrees = sorted({int(record['degree']) for record in manufactured})
    for axis in axes[0]:
        axis.set_xticks(degrees)
        axis.set_xlim(min(degrees)-.15, max(degrees)+.15)
        axis.legend(loc='lower left', framealpha=.9)

    temporal = sorted(temporal, key=lambda record: record['dt'])
    dt = np.array([positive(record['dt'], 'time step') for record in temporal])
    mass_error = np.array([positive(record['scaled_scalar_mass_difference_to_dt_0025'], 'scalar mass difference') for record in temporal])
    require(np.all(np.diff(dt) > 0), 'Repeated temporal step')
    guide = mass_error[-1]*(dt/dt[-1])**2
    axes[1, 0].loglog(dt, mass_error, 'o-', color='#2463a0', label='actual native IMEX2 differences')
    axes[1, 0].loglog(dt, guide, '--', color='#555555', label=r'$\Delta t^2$ guide (coarsest point anchor)')
    axes[1, 0].set(title='(c) Fixed mesh and final time: native IMEX2', xlabel=r'Time step $\Delta t$ [s]',
                   ylabel='Scaled scalar mass-norm difference\nto the run with $\\Delta t=0.0025$ s')
    axes[1, 0].set_xticks(dt)
    axes[1, 0].get_xaxis().set_major_formatter(ScalarFormatter())
    axes[1, 0].ticklabel_format(axis='x', style='plain')
    axes[1, 0].xaxis.set_minor_formatter(NullFormatter())
    axes[1, 0].legend(loc='upper left')

    controls = raw['transition']['controls']
    start, end, _, _, baseline = controls[:5]
    samples = sorted(raw['transition']['samples'], key=lambda sample: sample['input'][0])
    ell = np.array([sample['input'][0] for sample in samples])
    gates = np.array([sample['gates'] for sample in samples])
    require(np.all(np.isfinite(gates)) and np.all(gates[:, 0] >= 0) and np.all(gates[:, 0] <= 1), 'Invalid native gates')
    require(np.allclose(gates[:, 0], gates[:, 1]*gates[:, 2]*gates[:, 3], atol=2e-12, rtol=2e-12), 'Native gate product disagrees')
    dense = np.linspace(min(ell.min(), start)-.08, max(ell.max(), baseline)+.08, 800)
    q = np.clip((dense-start)/(end-start), 0, 1)
    spatial = 1-10*q**3+15*q**4-6*q**5
    sample_q = np.clip((ell-start)/(end-start), 0, 1)
    require(np.allclose(gates[:, 1], 1-10*sample_q**3+15*sample_q**4-6*sample_q**5,
                        atol=1e-11, rtol=1e-11), 'Actual operator spatial gate differs from its prescribed mask')
    require(np.all(gates[ell >= end, 0] == 0), 'Actual native operator reactivated downstream')
    axes[1, 1].axvspan(end, baseline, color='#e7be63', alpha=.27, label='trace-projection buffer')
    axes[1, 1].axvspan(baseline, dense[-1], color='#aaaaaa', alpha=.16, label='designated baseline trace faces')
    axes[1, 1].plot(dense, spatial, color='#444444', linewidth=1.3, label=r'prescribed $G_{\mathrm{sp}}$')
    axes[1, 1].plot(ell, gates[:, 0], 'o', color='#2463a0', markersize=4, label=r'qualified native $\chi$ samples')
    support=ell<end
    axes[1, 1].plot(ell[support], gates[support, 2], ':', color='#a95627', linewidth=1., label=r'evaluated $G_{\mathrm{st}}$ on mask support')
    for position in (start, end, baseline):
        axes[1, 1].axvline(position, color='#666666', linestyle=':', linewidth=.7)
    axes[1, 1].set(title='(d) One-way support and baseline trace buffer',
                   xlabel=r'Prescribed physical wall coordinate $\ell$ [m]', ylabel='Weight / gate', ylim=(-.045, 1.08),
                   xlim=(dense[0], dense[-1]))
    axes[1, 1].legend(loc='center right', fontsize=7, framealpha=.9)
    for axis in axes.flat:
        axis.grid(True, which='both', alpha=.2)
    commit = data.get('source_commit', data.get('tested_source_commit'))
    runner_receipt=qualification.parent.parent/'receipt.json'
    if commit is None and runner_receipt.is_file():
        receipt=json.loads(runner_receipt.read_text())
        require(receipt.get('status')=='passed' and receipt.get('finished_utc'), 'Source-pinning runner has not passed')
        require(receipt['source_commit']==receipt['source_end_commit'] and not receipt['source_end_status'], 'Source changed during qualification')
        commit=receipt['source_commit']
    suffix = ' · '+commit[:12] if commit else ''
    figure.suptitle('M6 prescribed transition: short-time native verification'+suffix, fontsize=13)
    files = []
    for extension in ('png', 'pdf'):
        path = args.output/('transition-verification.'+extension)
        figure.savefig(path, dpi=200)
        files.append(dict(path=str(path.resolve()), sha256=digest(path), bytes=path.stat().st_size))
    plt.close(figure)
    reductions = [float(mass_error[j+1]/mass_error[j]) for j in range(len(dt)-1)]
    metadata = dict(schema=1, milestone='M6', qualification=str(qualification),
                    qualification_sha256=digest(qualification), qualification_status=status,
                    source_commit=commit, runner_receipt=str(runner_receipt) if runner_receipt.is_file() else None,
                    runner_receipt_sha256=digest(runner_receipt) if runner_receipt.is_file() else None,
                    generated_utc=datetime.now(timezone.utc).isoformat(),
                    matplotlib_version=matplotlib.__version__, outputs=files,
                    panels=dict(epsilon=dict(series=epsilon_series, h_refinement=epsilon_h),
                                trace_gradient=dict(series=trace_series, h_refinement=trace_h,
                                    omission_meaning='physical RMS norm of the missing analytic (Ec-Ea)*grad(chi) trace term; not another numerical solve'),
                                temporal=dict(dt=dt.tolist(), actual_scaled_scalar_mass_difference=mass_error.tolist(),
                                    finest_reference_dt=.0025, guide=guide.tolist(), guide_power=2,
                                    guide_anchor_dt=float(dt[-1]), guide_anchor_error=float(mass_error[-1]),
                                    successive_coarsening_error_ratios=reductions,
                                    guide_meaning='reference slope anchored to actual coarsest difference; no fitted/extrapolated result'),
                                mask=dict(native_operator_report=str(raw_path), native_operator_sha256=digest(raw_path),
                                    dimension=raw['dimension'], map_kind=int(raw['map_parameters'][0]),
                                    transition_start=float(start), transition_end=float(end), baseline_faces_start=float(baseline),
                                    physical_ell=ell.tolist(), actual_gates=gates.tolist(),
                                    receipt_active_samples=operator_receipt.get('active_samples'),
                                    receipt_inactive_samples=operator_receipt.get('inactive_samples'),
                                    dense_spatial_curve='prescribed analytic mask with the qualified physical endpoints; samples checked against native output')),
                    scope='qualified local-operator and short-time MMS errors; no hill experimental accuracy or long-time stability claim',
                    rendering='Matplotlib Agg; positive measured errors only; no artificial error floors')
    (args.output/'figure-metadata.json').write_text(json.dumps(metadata, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(success=True, figure=str((args.output/'transition-verification.png').resolve()),
                          metadata=str((args.output/'figure-metadata.json').resolve()))))


if __name__ == '__main__':
    main()
