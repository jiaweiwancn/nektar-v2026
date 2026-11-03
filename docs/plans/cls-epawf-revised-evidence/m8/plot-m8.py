#!/usr/bin/env python3
"""Render M8 scientific evidence from immutable Linux JSONs; no solver execution.

The figure deliberately separates a reduced reference-local evolution test,
contradictory necessary fixed-gamma constraints, and zero-time projection
resolution. Optional curation copies small receipts and trims only large 9x9
interval tangent arrays; all complete originals remain under the Linux run root.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                               'MKL_NUM_THREADS', 'BLIS_NUM_THREADS',
                               'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


class Inputs:
    def __init__(self):
        self.files = {}

    def read(self, path):
        path = Path(path).resolve()
        value = json.loads(path.read_text())
        self.files[str(path)] = dict(path=str(path), sha256=digest(path), bytes=path.stat().st_size,
                                    source_commit=value.get('source_commit'))
        return value


def selected_trial(trials, dt):
    matches = [row for row in trials if row.get('operator') == 'p7-n8-xi4' and
               row.get('order') == 2 and row.get('amplitude') == 1e-4 and row.get('dt') == dt]
    if len(matches) != 1:
        raise ValueError('The requested p7/n8/xi4 IMEX2 trial is not unique')
    return matches[0]


def collect(root):
    inputs = Inputs()
    transverse = inputs.read(root/'transverse-001/validation.json')
    plan = inputs.read(root/'transverse-001/experiment-plan.json')
    assert transverse['success'] and transverse['checks'] == 305254
    original_path = Path(transverse['original_control']['path'])
    original = inputs.read(original_path)
    assert digest(original_path) == transverse['original_control']['sha256']
    trial = selected_trial(original['sections']['trials'], plan['native_spatial_time_step'])
    energies = [dict(label='Original CLS', beta=1., gamma_scale=0., trial=trial,
                     source_commit=original['source_commit'], report=str(original_path),
                     scope='Reused original M7 reduced-mode witness, not a new solver execution')]
    for row in transverse['variants']:
        section = row['sections']
        energies.append(dict(beta=section['beta'], gamma_scale=section['gamma_scale'],
                             trial=selected_trial(section['trials'], plan['native_spatial_time_step']),
                             source_commit=transverse['source_commit'], native_report=row['native_report'],
                             native_report_sha256=row['native_report_sha256'],
                             scope=section['qualification']))
    assert [(r['beta'], r['gamma_scale']) for r in energies[1:]] == [(.4, 0.), (.45, 0.), (.4, .25)]
    feasibility = inputs.read(root/'gamma-feasibility-001/validation.json')
    assert feasibility['success'] and feasibility['conclusion'] == 'infeasible_on_declared_box'
    constraints = []
    for row in feasibility['results']:
        assert row['empty_intersection']
        constraints.append(dict(beta=row['beta'], necessary_lower=row['necessary_gamma_lower'],
                                necessary_upper=row['necessary_gamma_upper'],
                                selected_principal_witness=row['selected_principal_witness'],
                                selected_covariance_witness=row['selected_covariance_witness'],
                                scope=row['interpretation']))
    coverage = inputs.read(root/'coverage-001/reference.json')
    p3 = inputs.read(coverage['inputs']['native']['path'])
    assert p3['audit_times'][0] == 0 and p3['eta'] == 1
    # The stored one-step anchored control has exactly identical coefficients.
    # Its aggregate coverage extrema therefore equal its initial projection.
    assert all(state == p3['accepted_states'][0] for state in p3['accepted_states'])
    resolution = [dict(degree=3, time=0., gradient_error=coverage['extrema']['gradient_entry_deviation'],
                       inner_covariance=p3['state_audits'][0]['source_volume_inner'][3],
                       coverage_status=coverage['status'], finite_coverage_pass=coverage['finite_coverage_pass'],
                       native=str(coverage['inputs']['native']['path']),
                       source_evolution_steps=p3['steps'], initial_and_control_coefficients_exactly_identical=True)]
    refinement = inputs.read(root/'reference-refinement-001/receipt.json')
    for case in refinement['cases']:
        raw = inputs.read(case['native_report'])
        sample = inputs.read(root/f'reference-refinement-001/p{case["degree"]}/state-coverage.json')
        assert case['steps'] == 0 and raw['steps'] == 0 and raw['audit_times'] == [0]
        assert raw['robust_parameters'] == p3['robust_parameters']
        resolution.append(dict(degree=case['degree'], time=0.,
            gradient_error=sample['extrema']['gradient_entry_deviation'],
            inner_covariance=raw['state_audits'][0]['source_volume_inner'][3],
            coverage_status=sample['status'], finite_coverage_pass=sample['finite_coverage_pass'],
            native=case['native_report'], native_sha256=case['native_sha256'], source_evolution_steps=0))
    assert [row['degree'] for row in resolution] == [3, 5, 7]
    return dict(schema=1, kind='M8 figure data extracted from completed JSON receipts',
        analysis_only=True, source_files=list(inputs.files.values()),
        source_commits=dict(original_transverse=original['source_commit'],
                            modified_transverse=transverse['source_commit'],
                            fixed_gamma_feasibility=feasibility['source_commit'],
                            reference_projection=refinement['source_commit']),
        transverse=dict(final_time=plan['final_time'], dt=plan['native_spatial_time_step'],
                        operator='p7-n8-xi4', amplitude=1e-4, order=2, trials=energies,
                        scope='Reduced Fourier/native-CG transverse subsystem at the neutral reference; no nonlinear terrain trajectory'),
        necessary_gamma=dict(constraints=constraints,
            box=dict(gradient_entry_radius='.25', k_over_K0=['.5', '1.5'], epsilon_over_epsilon0=['.5', '1.5']),
            scope='Two interval-certified necessary inequalities conflict. No plotted segment is asserted to be sufficient.'),
        initial_projection=dict(records=resolution, robust_parameters=p3['robust_parameters'],
            scope='Initial time only; covariance at native source quadrature and outer certificate coverage on declared finite evaluation sets; no evolved trajectory qualification'))


CURATED = (
    ('proof-kernel-fresh.json', 'kernel-independent-001/validation.json', False),
    ('proof-kernel-fresh-runner.json', 'kernel-independent-001-receipt.json', False),
    ('proof-kernel-replay.json', 'proof-correlated-001/kernel-replay/validation.json', False),
    ('proof-native-replay-provenance.json', 'proof-correlated-001/native-replay-provenance.json', False),
    ('proof-correlated-runner.json', 'proof-correlated-001/receipt.json', False),
    ('proof-local-joint-11aa.json', 'proof-correlated-001/joint/certificate.json', True),
    ('proof-local-joint-aa.json', 'certificate-joint-001/certificate.json', True),
    ('proof-local-joint-aa-runner.json', 'certificate-joint-001/receipt.json', False),
    ('proof-conservative-wide.json', 'proof-correlated-001/wide/certificate.json', True),
    ('proof-fixed-gamma-feasibility.json', 'gamma-feasibility-001/validation.json', False),
    ('proof-fixed-gamma-plan.json', 'gamma-feasibility-001/experiment-plan.json', False),
    ('proof-invariant-controls.json', 'invariant-certificate-001/certificate.json', False),
    ('proof-invariant-controls-runner.json', 'invariant-certificate-001/receipt.json', False),
    ('proof-invariant-screen.json', 'invariant-screen-001/validation.json', False),
    ('proof-invariant-plan.json', 'invariant-screen-001/experiment-plan.json', False),
    ('transverse-validation.json', 'transverse-001/validation.json', False),
    ('transverse-plan.json', 'transverse-001/experiment-plan.json', False),
    ('coverage-reference-p3.json', 'coverage-001/reference.json', False),
    ('coverage-perturbation.json', 'coverage-001/joint.json', False),
    ('coverage-runner.json', 'coverage-001/receipt.json', False),
    ('reference-refinement-runner.json', 'reference-refinement-001/receipt.json', False),
    ('coverage-reference-p5.json', 'reference-refinement-001/p5/state-coverage.json', False),
    ('coverage-reference-p7.json', 'reference-refinement-001/p7/state-coverage.json', False),
    ('reference-p5-energy.json', 'reference-refinement-001/p5/energy-audit.json', False),
    ('reference-p7-energy.json', 'reference-refinement-001/p7/energy-audit.json', False),
    ('runtime-build-runner.json', 'build-002/receipt.json', False),
)


def trim_tangents(value):
    if isinstance(value, dict):
        return {key: trim_tangents(child) for key, child in value.items() if key != 'tangent_intervals'}
    if isinstance(value, list):
        return [trim_tangents(child) for child in value]
    return value


def curate(root, target):
    target.mkdir(parents=True, exist_ok=False)
    manifest = dict(schema=1, kind='Compact M8 proof archive', generated_utc=datetime.now(timezone.utc).isoformat(),
        generator_sha256=digest(__file__), root=str(root), files=[],
        policy='Full copies except declared large 9x9 tangent-array removal. Exact witness states, inequalities, LDL pivots and source/binary/command provenance retained. Complete originals remain external.')
    for name, relative, trim in CURATED:
        source = root/relative
        value = json.loads(source.read_text())
        output = target/name
        if trim:
            value = trim_tangents(value)
            value['curation'] = dict(original_path=str(source), original_sha256=digest(source),
                omitted_key='tangent_intervals', reason='Large 9x9 interval tangent arrays remain in the immutable Linux original; all scalar bounds and interval LDL witnesses are retained.')
            save_json(output, value)
        else:
            output.write_bytes(source.read_bytes())
        manifest['files'].append(dict(name=name, original_path=str(source), original_sha256=digest(source),
            original_bytes=source.stat().st_size, curated_sha256=digest(output), curated_bytes=output.stat().st_size,
            transform='remove tangent_intervals recursively; add explicit curation metadata' if trim else 'byte-for-byte copy'))
    # Preserve the runtime image/binding/environment identity compactly; the
    # complete source/Python-package file manifests stay in the hashed original.
    source = root/'build-002/runtime-before.json'
    raw = json.loads(source.read_text())
    keys = ('schema', 'success', 'source_commit', 'source_root', 'build_root', 'captured_utc', 'clean',
            'comparison_policy', 'effective_environment_policy', 'runtime_artifacts', 'runtime_bindings',
            'runtime_environment', 'runtime_resolved_paths', 'mandatory_originals', 'failures')
    value = {key: raw[key] for key in keys}
    value['curation'] = dict(original_path=str(source), original_sha256=digest(source),
        omitted=['source_inputs', 'python_packages', 'tool_images', 'documentation_tools', 'ldd'],
        source_input_count=len(raw['source_inputs']), runtime_image_count=len(raw['runtime_artifacts']),
        scope='Runtime epoch summary only; not a replacement for the full external fingerprint or final source-equivalence audit')
    output = target/'runtime-epoch-summary.json'
    save_json(output, value)
    manifest['files'].append(dict(name=output.name, original_path=str(source), original_sha256=digest(source),
        original_bytes=source.stat().st_size, curated_sha256=digest(output), curated_bytes=output.stat().st_size,
        transform='declared runtime/source summary; see curation field'))
    save_json(target/'curation-index.json', manifest)
    return manifest


def render(data, output):
    blue, teal, red, orange, gray = '#2563eb', '#0f766e', '#b91c1c', '#b45309', '#64748b'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 12, 'axes.titlesize': 15,
                         'axes.labelsize': 12, 'xtick.labelsize': 11, 'ytick.labelsize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    figure, axes = plt.subplots(1, 3, figsize=(17, 6.8), gridspec_kw={'width_ratios': [1., 1.1, 1.2]})
    figure.subplots_adjust(left=.055, right=.945, bottom=.30, top=.80, wspace=.40)
    figure.suptitle('M8: local verification and failed broad qualification', fontsize=19, y=.965)
    axis = axes[0]
    trials = data['transverse']['trials']
    gains = [row['trial']['gain'] for row in trials]
    labels = ['Original\nCLS']+[rf'$\beta={row["beta"]:g}$'+('\n'+rf'$\Gamma={row["gamma_scale"]:g}$' if row['gamma_scale'] else '') for row in trials[1:]]
    axis.bar(range(len(gains)), gains, color=[red, blue, teal, '#7c3aed'], width=.64, zorder=3)
    axis.set_yscale('log')
    axis.set_ylim(min(gains)*.62, max(gains)*1.62)
    axis.axhline(1, color=gray, linestyle='--', linewidth=1.2)
    axis.set_xticks(range(len(labels)), labels)
    axis.set_ylabel(r'$E(T)/E(0)$')
    axis.set_title('(a) Reduced transverse mode', pad=28)
    axis.text(.5, 1.025, rf'$\xi_x=4\;\mathrm{{m}}^{{-1}},\ T={data["transverse"]["final_time"]:g}\;\mathrm{{s}},\ p=7,\ n_z=8$',
              ha='center', transform=axis.transAxes, fontsize=11)
    for index, value in enumerate(gains):
        axis.annotate(f'{value:.3f}', (index, value), xytext=(0, 7), textcoords='offset points', ha='center', fontsize=12)
    axis.grid(axis='y', which='major', alpha=.18, zorder=0)
    axis.text(.5, -.28, 'Reference-local subsystem only\n'+rf'IMEX2; $\Delta t={data["transverse"]["dt"]:g}$ s',
              transform=axis.transAxes, ha='center', va='top', fontsize=11, color=gray)

    axis = axes[1]
    constraints = data['necessary_gamma']['constraints']
    lower = [float(row['necessary_lower']) for row in constraints]
    upper = [float(row['necessary_upper']) for row in constraints]
    right = max(lower)*1.23
    axis.set_xlim(0, right)
    axis.set_ylim(-.45, 1.65)
    axis.set_yticks([1, 0], [r'$\beta=0.4$', r'$\beta=0.45$'])
    for y, lo, hi in zip((1, 0), lower, upper):
        axis.plot(hi, y+.15, marker='s', color=teal, markersize=7)
        axis.annotate('', (0, y+.15), (hi, y+.15), arrowprops={'arrowstyle': '->', 'color': teal, 'lw': 1.5})
        axis.text(hi+.08, y+.15, rf'$\Gamma\leq {hi:.3f}$', va='center', color=teal, fontsize=12)
        axis.plot(lo, y-.15, marker='o', color=red, markersize=7)
        axis.annotate('', (right, y-.15), (lo, y-.15), arrowprops={'arrowstyle': '->', 'color': red, 'lw': 1.5})
        axis.text(lo-.07, y-.15, rf'$\Gamma\geq {lo:.3f}$', ha='right', va='center', color=red, fontsize=12)
    axis.set_xlabel(r'Prescribed $\Gamma=\gamma/\nu_{t,0}$')
    axis.set_title('(b) Necessary bounds conflict', pad=28)
    axis.text(.5, 1.025, 'Same declared finite state box', ha='center', transform=axis.transAxes, fontsize=11)
    axis.legend(handles=[Line2D([], [], color=teal, marker='s', linestyle='', label='Covariance upper'),
                         Line2D([], [], color=red, marker='o', linestyle='', label='Momentum lower')],
                loc='upper center', frameon=False, fontsize=10, ncol=1, bbox_to_anchor=(.54, 1.00))
    axis.text(.5, -.28, 'Necessary inequalities have empty intersection\nNo sufficient region is plotted',
              transform=axis.transAxes, ha='center', va='top', fontsize=11, color=gray)

    axis = axes[2]
    records = data['initial_projection']['records']
    degree = [row['degree'] for row in records]
    error = [row['gradient_error'] for row in records]
    covariance = [row['inner_covariance'] for row in records]
    first, = axis.plot(degree, error, color=blue, marker='o', markersize=7, linewidth=1.7, label='Gradient coverage error')
    axis.set_yscale('log')
    axis.set_ylim(min(error)*.45, max(error)*2.5)
    axis.set_xticks(degree)
    axis.set_xlabel('Polynomial degree p')
    axis.set_ylabel(r'Max entry of $|A^*-A^*_0|$', color=blue)
    axis.tick_params(axis='y', colors=blue)
    axis.grid(axis='y', which='major', alpha=.18)
    other = axis.twinx()
    other.spines['right'].set_visible(True)
    second, = other.plot(degree, covariance, color=orange, marker='s', markersize=7, linewidth=1.7,
                          label='Sampled inner covariance')
    margin = (max(covariance)-min(covariance))*.17
    other.set_ylim(min(covariance)-margin, max(covariance)+margin)
    other.axhline(0, color=orange, linestyle='--', linewidth=1.0, alpha=.65)
    other.set_ylabel(r'Min $\lambda(R_A)$ [m$^2$/s$^2$]', color=orange)
    other.tick_params(axis='y', colors=orange)
    for index, (x, value) in enumerate(zip(degree, error)):
        offset, align = ((4, -18), 'left') if index == 0 else ((-4, -18), 'right') if index == len(degree)-1 else ((0, -18), 'center')
        axis.annotate(f'{value:.3g}', (x, value), xytext=offset, textcoords='offset points',
                      ha=align, color=blue, fontsize=10)
    for index, (x, value) in enumerate(zip(degree, covariance)):
        offset, align = ((6, 9), 'left') if index == 0 else ((-8, -17), 'right') if index == len(degree)-1 else ((0, 9), 'center')
        other.annotate(f'{value:+.3g}', (x, value), xytext=offset, textcoords='offset points',
                       ha=align, color=orange, fontsize=10)
    axis.set_title('(c) Initial reference projection', pad=28)
    parameters = data['initial_projection']['robust_parameters']
    axis.text(.5, 1.025, rf'$t=0,\ \beta={parameters[0]:g},\ \Gamma={parameters[5]:g}$',
              ha='center', transform=axis.transAxes, fontsize=11)
    axis.legend(handles=[first, second], loc='upper center', bbox_to_anchor=(.47, .99), frameon=False, fontsize=10)
    axis.text(.5, -.28, 'Zero-time resolution check only\np7 finite coverage is not an evolved trajectory',
              transform=axis.transAxes, ha='center', va='top', fontsize=11, color=gray)
    figure.text(.5, .035, 'No terrain accuracy claim. Local certificates and finite evaluation sets do not establish all-state or coupled stability.',
                ha='center', fontsize=11, color=gray)
    png, pdf = output/'m8-qualification.png', output/'m8-qualification.pdf'
    figure.savefig(png, dpi=220, facecolor='white', metadata={'Software': 'M8 evidence plot'})
    figure.savefig(pdf, facecolor='white', metadata={'Creator': 'plot-m8.py', 'CreationDate': None, 'ModDate': None})
    plt.close(figure)
    return [png, pdf]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--root', type=Path, help='Actual Linux M8 run root')
    source.add_argument('--data', type=Path, help='Previously extracted exact figure-data JSON')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--curate', type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific figure analysis and rendering run on the authoritative Linux host')
    args.output.mkdir(parents=True, exist_ok=False)
    if args.curate and not args.root:
        raise ValueError('Curation requires the original Linux run root')
    data = collect(args.root.resolve()) if args.root else json.loads(args.data.read_text())
    data_path = args.output/'figure-data.json'
    save_json(data_path, data)
    files = render(data, args.output)
    if args.curate:
        curate(args.root.resolve(), args.curate.resolve())
    receipt = dict(schema=1, kind='M8 standalone scientific figure', status='rendered',
        generated_utc=datetime.now(timezone.utc).isoformat(), analysis_only=True, solver_executed=False,
        command=[sys.executable, *sys.argv], plot_sha256=digest(__file__), figure_data_sha256=digest(data_path),
        source_commits=data['source_commits'], inputs=data['source_files'],
        cpu_affinity=sorted(os.sched_getaffinity(0)), numerical_threads=THREADS,
        versions=dict(python=sys.version, numpy=np.__version__, matplotlib=matplotlib.__version__),
        outputs=[dict(path=str(path), sha256=digest(path), bytes=path.stat().st_size) for path in files],
        curation=str(args.curate.resolve()) if args.curate else None,
        scope='Visual summary of already qualified/rejected evidence; rendering is not a new solver or model qualification')
    save_json(args.output/'figure-receipt.json', receipt)
    print(json.dumps(dict(status=receipt['status'], output=str(args.output),
                          outputs=receipt['outputs'], solver_executed=False), indent=2))


if __name__ == '__main__':
    main()
