#!/usr/bin/env python3
"""Record a CV3 bounded-design attempt without declaring calibration complete.

Read-only with respect to native runs. Reserved flat and hill data are never
loaded. Keep large native fields/logs in the Linux run root.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
from cao_tamura_experiment import ROOT, DATA, records, sha, write_csv, write_json


def bound(path):
    return dict(path=str(path), sha256=sha(path))


def last_native(folder, stem='case'):
    files = sorted(folder.glob('stage-*/case.rans.csv')) if (folder/'development.json').exists() else [folder/(stem+'.rans.csv')]
    available = [f for f in files if f.exists() and len(records(f))]
    if not available:
        return {}
    path = available[-1]
    rows = records(path)
    return dict(last_time=float(rows[-1]['time']), native_rows=len(rows), diagnostics=bound(path),
        last_min_k=float(rows[-1]['min_k']), last_min_epsilon=float(rows[-1]['min_eps']),
        last_strong_divergence=float(rows[-1]['divergence_l2']),
        last_weak_divergence=float(rows[-1]['weak_divergence_l2']),
        maximum_nonpositive_linear_images=max(float(r.get('nonpositive_linear_images',0))
            for f in available for r in records(f)))


def figures(conditional, run_root, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    points = records(DATA/'profiles.csv')
    fig, axes = plt.subplots(3, 2, figsize=(8, 8.7), constrained_layout=True)
    labels = (r'$U/U_\infty$', r'$-R_{xz}/U_\infty^2$', r'$\sigma_x/u_{\tau,\mathrm{source}}$')
    for column, surface in enumerate(('S','R')):
        item = conditional.get(surface+'/standard_original')
        for row, panel in enumerate(('a','b','c')):
            ax = axes[row, column]
            family = 'f3'+panel+'_'+surface
            target = [p for p in points if p['series_id']==family and float(p['weight'])>0]
            ax.errorbar([float(p['value']) for p in target], [.04*float(p['y_over_H']) for p in target],
                xerr=[float(p['value_bound']) for p in target], fmt='o', markersize=3,
                color='black', elinewidth=.6, label='CV1 digitization (value bounds)')
            if item:
                values = sorted(item['objective']['families'][family]['values'], key=lambda p:p['z'])
                values = [p for p in values if p['z']>=.022]
                ax.plot([p['value'] for p in values], [p['z'] for p in values], '-s',
                    markersize=3, color='#1769aa', label='Standard, original AWF')
                ax.set_title(surface+': '+item['candidate'] if row==0 else '')
            ax.axhspan(0., .022, color='.92', zorder=-1)
            ax.axhline(.022, color='.5', linestyle=':', linewidth=.8)
            ax.set_ylim(0., max(.3, max(.04*float(p['y_over_H']) for p in target)*1.03))
            ax.set_xlabel(labels[row]);ax.set_ylabel('Physical height (m)');ax.grid(alpha=.2)
            if panel=='b':ax.ticklabel_format(axis='x',style='sci',scilimits=(-2,2),useMathText=True)
        axes[0, column].legend(fontsize=7)
    fig.suptitle('Conditional standard-only design minima; no shared calibration freeze\n'
                 'Shaded heights excluded from fitting; numerical uncertainty pending CV4', fontsize=10)
    fig.savefig(output/'conditional-profiles.png', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True)
    qrows = records(run_root/'calibration/cv3/tangents-q.csv')
    axes[0].plot([float(r['normalized_shear']) for r in qrows],
        [1e6*float(r['min_transverse_tangent_h1e6']) for r in qrows], label='CLS native, 361 directions')
    axes[0].axhline(1e6*(3.8e-6+.09*.02**2/.1), label='Standard analytic', color='.4')
    axes[0].axhline(0., color='black', linewidth=.8);axes[0].set_xlim(0., 5.)
    axes[0].set_xlabel(r'$q=kG/\epsilon$');axes[0].set_ylabel(r'Minimum tangent ($10^{-6}$ viscosity units)')
    axes[0].set_title(r'Fixed $k=0.02$, $\epsilon=0.1$, $\nu=3.8\times10^{-6}$')
    for surface in ('S','R'):
        values = records(run_root/'calibration/cv3'/('tangents-'+surface+'.csv'))
        values = [r for r in values if r['id'].startswith(surface+'-measured-')]
        if not values:raise ValueError('Missing measured-rms inlet tangent family: '+surface)
        axes[1].plot([1e6*float(r['min_transverse_tangent_h1e6']) for r in values],
            [float(r['eta']) for r in values], label=surface+' inferred inlet')
    axes[1].axvline(0., color='black', linewidth=.8);axes[1].set_ylim(.011, .3)
    axes[1].set_xlabel(r'Minimum CLS tangent ($10^{-6}$ viscosity units)')
    axes[1].set_ylabel('Physical height (m)');axes[1].set_title('MeasuredRms / ShearLength family')
    for ax in axes:
        ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.savefig(output/'stress-tangent.png', dpi=180)
    plt.close(fig)


def record(campaign, run_root, build, output):
    output.mkdir(parents=True, exist_ok=False)
    selected = json.loads((campaign/'selection.json').read_text())
    recipe = json.loads((campaign/'recipe.json').read_text())
    cases, flat_rows, conditional = [], [], {}
    for candidate in selected['candidates']:
        surface = candidate['surface']
        for model, entry in candidate['models'].items():
            folder = Path(entry['directory'])
            state = json.loads((folder/'development.json').read_text())
            item = dict(surface=surface, candidate=candidate['id'], model=model,
                parameters=candidate['parameters'], status=state['status'], complete=state['complete'],
                development=bound(folder/'development.json'), native=last_native(folder),
                completed_development_segments=len(state['stages']), error=state.get('error'),
                final_drifts=[s['drift'] for s in state['stages'][-2:] if 'drift' in s])
            if state['complete']:
                objective = entry['objective']
                item.update(objective=objective, statistics=bound(folder/'statistics.csv'),
                    final_flow_times=state['stages'][-1]['elapsed_flow_times'],
                    maximum_free_residual=max(s['maximum_free_residual'] for s in state['stages']),
                    maximum_strip_identity_error=max(s['strip_momentum'].get('maximum_inventory_identity_error',0.) for s in state['stages']),
                    maximum_independent_strip_mass_defect=max(s['strip_continuity']['errors']['actual_flux_plus_side_flux'] for s in state['stages']))
                if surface=='R':
                    item['observational_displacement_scenarios']={}
                    for d in (.0025,.0035,.0045):
                        path=folder/f'objective-d{d:.4f}.json'
                        scenario=json.loads(path.read_text())
                        item['observational_displacement_scenarios'][str(d)]=dict(
                            source=bound(path), objective=scenario['objective'],
                            rms_score=scenario['families']['f3c_R']['score'],
                            rms_conditional_error=scenario['families']['f3c_R']['conditional_rmse'])
                if objective['objective'] is not None:
                    key = surface+'/'+model
                    if key not in conditional or objective['objective'] < conditional[key]['objective']['objective']:
                        conditional[key] = item
            cases.append(item)
            row = {k:item[k] for k in ('surface','candidate','model','status','complete','completed_development_segments')}
            row.update(candidate['parameters'])
            row.update({k:item['native'].get(k,'') for k in ('last_time','last_min_k','last_min_epsilon','last_strong_divergence')})
            if state['complete']:
                f = item['objective']['families']
                row.update(J=item['objective']['objective'], final_flow_times=item['final_flow_times'],
                    friction_ratio=f['friction']['value']/f['friction']['target'], delta99=f['thickness']['delta99'],
                    delta95=f['thickness']['delta95'],
                    U_profile_RMSE_over_Uinf=f['f3a_'+surface]['conditional_rmse'],
                    actual_ReH=item['objective']['incident']['actual_ReH'])
            flat_rows.append(row)
    write_csv(output/'design.csv', flat_rows, list(dict.fromkeys(k for r in flat_rows for k in r)))
    profiles = []
    for key, item in conditional.items():
        for family, values in item['objective']['families'].items():
            for v in values.get('values',[]):
                profiles.append(dict(surface=item['surface'], model=item['model'], candidate=item['candidate'],
                    family=family, **v))
    write_csv(output/'conditional-profiles.csv', profiles, list(profiles[0]) if profiles else ['surface'])
    controls = []
    for folder in sorted((run_root/'calibration/cv3/preflight').iterdir()):
        manifest, receipt = folder/'manifest.json', folder/'run.json'
        if not manifest.exists() or not receipt.exists():
            continue
        m, r = json.loads(manifest.read_text()), json.loads(receipt.read_text())
        xml = [f for f in folder.glob('*.xml') if sha(f)==r.get('case_sha256',m['case_sha256'])]
        if len(xml)!=1:
            controls.append(dict(name=folder.name, surface=m['config']['surface'], model=m['config']['model'],
                calibration_eligible=False, metadata_status='producing XML not uniquely present in this folder',
                producing_run=bound(receipt)))
            continue
        case=xml[0];root = ET.parse(case).getroot()
        order = int(root.find('CONDITIONS/TIMEINTEGRATIONSCHEME/ORDER').text)
        params = {n.text.split('=',1)[0].strip():n.text.split('=',1)[1].strip()
                  for n in root.findall('CONDITIONS/PARAMETERS/P')}
        floor = root.find("CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='0']")
        corrected = all(n.attrib.get('VALUE') == 'EpsilonWall' for n in floor
                        if n.attrib.get('VAR') in ('eps','eps_tilde'))
        log = folder/'run.log'
        errors = re.findall(r'RANS: [^\n]+', log.read_text(errors='replace')) if log.exists() else []
        controls.append(dict(name=folder.name, surface=m['config']['surface'], model=m['config']['model'],
            time_order=order, dt=m['config']['dt'], degree=m['config']['degree'],
            quadrature_points=root.find('EXPANSIONS/E').attrib['NUMPOINTS'],
            grad_div=params.get('RANSGradDiv','0'), startup_duration=params.get('RANSStartupDuration','0'),
            exit_code=r.get('exit_code'), timeout=r.get('timeout',False),
            complete_native_execution=r.get('exit_code') == 0 and not r.get('timeout',False),
            corrected_epsilon_boundary=corrected, supported_time_method=order in (1,2),
            calibration_eligible=False, error=errors[-1] if errors else '',
            **last_native(folder,case.stem), producing_run=bound(receipt), case=bound(case),
            solver_sha256=r.get('solver_sha256')))
    control_rows = [{k:v for k,v in r.items() if not isinstance(v,dict)} for r in controls]
    write_csv(output/'preflight.csv', control_rows, list(dict.fromkeys(k for r in control_rows for k in r)))
    regression = {}
    for name, file in [('implicit_order_and_history','regression-domain-implicit/results.json'),
                       ('native_inlet_startup_history','regression-cv3-delivery/evidence.json'),
                       ('grad_div_neutral_controls','regression-graddiv-neutral/evidence.json'),
                       ('stress_tangent_oracle','stress-tangent-oracle/results.json')]:
        path = run_root/'calibration/cv3'/file
        if path.exists():
            j = json.loads(path.read_text())
            summary = {k:j[k] for k in ('status','accepted','scope','solver_sha256','source_revision','binary_sha256','source_parent') if k in j}
            summary['result_names'] = list(j.get('results',{}))
            if 'comparisons' in j:summary['comparisons']=j['comparisons']
            summary['tests'] = {}
            for test, value in j.get('tests',{}).items():
                if isinstance(value,dict):
                    summary['tests'][test] = {k:value[k] for k in ('rates','maximum_scaled_difference') if k in value}
            if name=='stress_tangent_oracle':
                summary['maximum_scaled_error'] = {m:r['maximum_scaled_error'] for m,r in j['results'].items()}
                summary['rows'] = j['results']['cls']['rows']
            regression[name] = dict(**bound(path), summary=summary)
    tangent = {}
    for surface in ('S','R'):
        path = run_root/'calibration/cv3'/('tangents-'+surface+'.csv')
        rows = records(path)
        for family in ('old','measured'):
            values = [r for r in rows if r['id'].startswith(surface+'-'+family+'-')]
            if not values:
                continue
            minimum = min(values, key=lambda r:float(r['min_transverse_tangent_h1e6']))
            tangent[surface+'-'+family] = dict(source=bound(path), samples=len(values), minimum=minimum,
                maximum_absolute_normalized_shear=max(abs(float(r['normalized_shear'])) for r in values))
    post_domain=[]
    for path in sorted((run_root/'calibration/cv3/domain-image-targets').glob('*/*/cls_original/development.json')):
        j=json.loads(path.read_text())
        post_domain.append(dict(surface=j['config']['surface'],candidate=path.parents[1].name,
            status=j['status'],complete=j['complete'],error=j.get('error'),
            completed_development_segments=len(j['stages']),native=last_native(path.parent),
            development=bound(path),recipe=bound(path.parents[1]/'recipe.json')))
    source_files = ['solvers/IncNavierStokesSolver/EquationSystems/'+f for f in (
        'RANSVelocityCorrectionScheme.cpp','RANSVelocityCorrectionScheme.h','RANSMatrixFree.cpp','RANSMatrixFree.h')]
    source_files += ['solvers/IncNavierStokesSolver/Tests/RANS/'+f for f in (
        'cao_tamura_experiment.py','cao_tamura_calibration.py','cao_tamura_campaign.py',
        'cao_tamura_cv3_budgets.py','terrain_continuity.py','validate_cao_tamura_cv3.py',
        'validate_cao_tamura_calibration.py','validate_rans_implicit_advection.py',
        'validate_rans_stress_tangent.py','record_cao_tamura_cv3.py','CMakeLists.txt')]
    source_files += ['library/UnitTests/SolverUtils/RANSStressTangentAudit.cpp']
    figures(conditional,run_root,output)
    outcome = dict(milestone='CV3', status='INCOMPLETE', milestone_complete=False,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        reason=('No developed CLS candidate in the bounded design; no common S/R calibration can be frozen.'
                if not any(c['model']=='cls_original' and c['complete'] for c in cases) else
                'This bounded-design record does not seal wall transfers or a calibration; use the gated freeze workflow.'),
        scope='bounded flat-data design and numerical preflight; no reserved-flat or hill scores',
        source_parent=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source_sha256={n:sha(ROOT/n) for n in source_files},
        cv1_snapshot_sha256=sha(DATA/'snapshot.json'), selection=bound(campaign/'selection.json'),
        recipe=recipe, selected_common=selected['selected'],
        counts=dict(Counter(c['model']+':'+c['status'] for c in cases)),
        cases=cases, conditional_design_minima=conditional,post_domain_targets=post_domain,
        conditional_minima_are_frozen_calibrations=False, preflight=controls,
        regressions=regression, constitutive_diagnostic=tangent,
        build_receipt=json.loads((build/'build-receipt.json').read_text()),
        binaries={n:sha(build/n) for n in ('IncNavierStokesSolver','RANSWallFieldAudit','RANSObservationAudit','RANSFieldAudit')},
        stress_tangent_binary=bound(run_root/'build/RANSStressTangentAudit'),
        ctest=bound(run_root/'calibration/cv3/ctest-cv3-delivery.log'),
        artifacts={f.name:sha(f) for f in sorted(output.iterdir()) if f.is_file()},
        unexecuted=['shared S/R parameter freeze','developed ten-case wall-policy transfer',
            'reserved Fig.3d crosschecks','final common-calibration uncertainty assessment'],
        limitations=['No global canonical-model infeasibility claim.',
            'Local negative tangent alone does not establish the cause of the full solver rejection.',
            'Conditional standard minima are at the tested inlet-thickness bound; no interior/global optimum is identified.',
            'CV4 h/p/dt/domain/3D and observable numerical uncertainty remain unestablished.',
            'All preflight pilots, including exit-zero pilots, remain ineligible for calibration.'])
    write_json(output/'attempt.json',outcome)
    return outcome


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--run-root',type=Path,required=True)
    p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=record(args.campaign.resolve(),args.run_root.resolve(),args.build.resolve(),args.output.resolve())
    print(json.dumps(dict(status=result['status'],counts=result['counts'],selected_common=result['selected_common']),indent=2))


if __name__=='__main__':
    main()
