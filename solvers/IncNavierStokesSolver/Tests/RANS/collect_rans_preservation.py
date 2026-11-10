#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact evidence with explicit coverage and independent gate statuses.

Exit 0: required gates pass; 1: failed/missing required gate; 2: unresolved
unperturbed temporal observable. This collector never launches a solver.
"""
import argparse
import hashlib
import math
from pathlib import Path
import subprocess
import sys
from rans_qualification import write_json,coverage,audited_extrema,required_native_checks,native_check_coverage
from qualify_rans_preservation import load,MODELS,POLICIES,WALLS

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root;a.output.mkdir(parents=True,exist_ok=False)
    gates={};errors=[]
    def read(label,path,required='pass'):
        try:
            d=load(r/path)
            if required and d.get('status')!=required:errors.append(dict(input=label,status=d.get('status'),required=required))
            return d
        except Exception as error:errors.append(dict(input=label,failure=repr(error)));return {}
    manifest=read('runtime','runtime/source-build.json',None)
    kernel=read('independent-kernels','kernel-results.json');gates['independent_kernels']=kernel.get('status','fail')
    audit_names={'inactive_chi_zero','outside_spatial_support','active_invalid_cutoff','accepted_invalid_cutoff','inactive_invalid_scalar','inactive_invalid_constitutive','near_boundary_A95_forward','covariance_scale_ordinary','covariance_scale_small'}
    before=read('before-native','coherent-before-audits/audits.json');after=read('after-native','coherent-after-audits/audits.json');diagnostic=read('diagnostic','diagnostic-extended-checks.json')
    for label,d in (('before',before),('after',after)):
        d['coverage']=coverage(d.get('checks',[]),audit_names,key=lambda x:x['name'])
        if d['coverage']['status']!='pass' or any(x['status']!='pass' for x in d.get('checks',[])):errors.append(dict(input=label,reason='nine native audit controls failed/missing'))
    write_json(a.output/'audit-fixes.json',dict(before=before,after=after,eigen_and_A95=diagnostic,before_eigen=read('old-eigen','before/eigen.json',None),coherent_scale_controls=read('native-scale','covariance-native-scale-controls.json',None)))
    matrix=read('matrix','matrix/matrix.json');wanted={(m,p,w) for m in MODELS for p in POLICIES for w in WALLS}
    matrix_coverage=coverage(matrix.get('cases',[]),wanted)
    if matrix_coverage['status']!='pass':errors.append(dict(input='matrix-coverage',result=matrix_coverage))
    compact=[]
    for row in matrix.get('cases',[]):
        d=load(Path(row['path'])/'case.evolution.json');ext=audited_extrema(d)
        domains=[(s['time'],phase,x) for s in d['pimple_attempts'] for phase,x in s.items() if phase.endswith('_domain') and isinstance(x,dict)]
        cutoffs=[v for _,_,x in domains for v in x.get('cutoff_range',[]) if v>0]
        compact.append(dict(model=row['model'],kappa_policy=row['kappa_policy'],wall=row['wall'],execution=row['status'],time=d['time'],steps=d['steps'],
            all_state_stage_extrema=ext['over_all_audited_states_and_stages'],
            final_state={k:v for k,v in (ext['final_state'] or {}).items() if 'covariance' in k or 'minimum' in k or 'radicand' in k or 'cutoff' in k or k in ('coupled_norms','interface','common_velocity_k_dependent_epsilon_trace_defects')},
            physical_margin_inventory='actual native states plus declared dependent-trace targets; no tangent stencils; covariance provenance block identifies trace-target observations',
            initial_spatial_alpha_zero=d['raw_reference_norms'],final_time_discrete_stage=d['stage_audits'][-1],
            correction_budgets=dict(bound_calls=sum(s.get('bound_calls',0) for s in d['pimple_attempts']),signed_integral=0.,absolute_integral=0.),
            cutoff_positive_inventory_range=[min(cutoffs),max(cutoffs)] if cutoffs else [],mesh_sha256=row['xml_sha256'],raw_report=str(Path(row['path'])/'case.evolution.json')))
        if row['status']!='pass' or any(s.get('bound_calls',0) for s in d['pimple_attempts']):errors.append(dict(input=row['path'],reason='execution/correction gate'))
    beforeafter=read('bitwise-trajectories','before-after-trajectories.json')
    write_json(a.output/'matrix-results.json',dict(coverage=matrix_coverage,cases=compact,before_after=beforeafter))
    spatial=read('spatial','preservation-spatial/results.json');sc=[];weak=[]
    for folder in sorted((r/'preservation-spatial').iterdir()):
        path=folder/'metrics.json'
        if not path.exists():continue
        d=load(path);alpha=d['native_alpha_zero']
        d['actual_quadrature']['coefficient_initializer_wall_normal']=32 if folder.name.endswith(('-q1','-q2')) else d['configuration']['quadrature'][0]
        d['actual_quadrature']['accepted_referenceK_nested_integral']='same normal order as coefficient initializer, before the fixed Gauss32 radial projection' if '-Accepted-' in folder.name else 'not applicable (constant inner k)'
        sc.append(dict(name=folder.name,configuration=d['configuration'],actual_quadrature=d['actual_quadrature'],state_sha256=d['state_sha256'],
            representation=d['representation'],coupled_norms=alpha['coupled_norms'],k_inner_outer_shared=alpha['native_k_inner_outer_shared_interface_dual'],
            interface=alpha['interface'],dependent_trace=alpha['common_velocity_k_dependent_epsilon_trace_defects'],pressure_and_gauge=alpha['coupled_norms'][4:]))
        if 'analytical_weak_balance' in d:weak.append(dict(name=folder.name,**d['analytical_weak_balance']))
    expected_spatial={m+'-'+p+'-'+w+'-reference' for m in MODELS for p in POLICIES for w in WALLS}
    for m in MODELS:
        for policy in POLICIES:
            for w in WALLS:
                tag='-'.join((m,policy,w))
                expected_spatial.update(tag+'-'+s for s in (('q1','q2','h1','h2','p4','p5') if m=='E94E95' and policy=='Natural' else ('h1',)))
    spatial_coverage=coverage(spatial.get('checks',[]),expected_spatial,key=lambda x:x['name'])
    if spatial_coverage['status']!='pass' or any(x['status']!='pass' for x in spatial.get('checks',[])):errors.append(dict(input='spatial-coverage',result=spatial_coverage))
    tables={tag:{k:v for k,v in table.items() if k!='levels'} for tag,table in spatial.get('spatial_tables',{}).items()}
    write_json(a.output/'spatial-results.json',dict(status=spatial.get('status'),coverage=spatial_coverage,tables=tables,levels=sc))
    write_json(a.output/'analytical-weak-results.json',dict(status='pass' if len(weak)==14 and all(w['status']=='pass' for w in weak) else 'fail',results=weak))
    if len(weak)!=14 or any(w['status']!='pass' for w in weak):errors.append(dict(input='analytical-native-tests',reason='14 weak evaluations required'))
    temporal=[];perturb=[];temporal_runs=[];perturb_runs=[];timings=[]
    for model in MODELS:
        td=read(model+'-temporal','preservation-temporal-Q0-'+model+'/reanalysis.json',None)
        if td.get('status') not in ('pass','unresolved'):errors.append(dict(input=model+'-temporal',status=td.get('status')))
        temporal.extend(td.get('temporal_tables',[]))
        pd=read(model+'-perturbed','preservation-temporal-perturbed-'+model+'/results.json',None)
        if pd.get('status') not in ('pass','unresolved'):errors.append(dict(input=model+'-perturbed',status=pd.get('status')))
        if model=='A95':
            target=read('targeted-A95-EP-signal','preservation-temporal-perturbed-A95-EP-1pct/results.json')
            pd=dict(pd,temporal_tables=[t for t in pd.get('temporal_tables',[]) if t['wall']!='EP']+target.get('temporal_tables',[]))
        perturb.extend(pd.get('temporal_tables',[]))
        for prefix,container in (('Q0-',temporal_runs),('perturbed-',perturb_runs)):
            for f in sorted((r/('preservation-temporal-'+prefix+model)).glob('*/temporal-metrics.json')):
                d=load(f)
                timing=f.parent/'timing.jsonl'
                if timing.exists():
                    groups=dict(assembly=0.,solve=0.,diagnostic=0.,output=0.,other=0.)
                    for line in timing.read_text().splitlines():
                        import json
                        event=json.loads(line);name=event['name']
                        group='output' if name.startswith(('io.','output.')) else 'diagnostic' if 'audit' in name or 'domain' in name else 'assembly' if any(x in name for x in ('assembl','residual','evaluation','profile','source','transport')) else 'solve' if 'solve' in name or 'factor' in name else 'other'
                        groups[group]+=event['exclusive_wall_seconds']
                    timings.append(dict(case=str(f.parent),exclusive_seconds=groups,definition='sum exclusive instrumented scopes by label; startup, Python and unscoped work are separate'))
                container.append(dict(name=f.parent.name,**{k:v for k,v in d.items() if k not in ('extrema',)},extrema=d['extrema']['over_all_audited_states_and_stages']))
    for label,rows,expected in (('temporal',temporal,wanted),('perturbed',perturb,{(m,'Natural',w) for m in MODELS for w in WALLS})):
        check=coverage(rows,expected)
        if check['status']!='pass':errors.append(dict(input=label,coverage=check))
        for row in rows:
            if set(row['fields'])!={'u','v','k','epsilon_outer','Pi'} or any(f['status']=='fail' for f in row['fields'].values()):errors.append(dict(input=label,reason='field coverage or measured order failed',model=row['model'],wall=row['wall']))
    unresolved=[dict(model=t['model'],policy=t['kappa_policy'],wall=t['wall'],field=name,**value) for t in temporal for name,value in t['fields'].items() if value['status']=='unresolved_floor']
    write_json(a.output/'temporal-results.json',dict(execution='completed' if len(temporal_runs)==36 else 'missing',convergence='partially_unresolved' if unresolved else 'demonstrated',tables=temporal,runs=temporal_runs,unresolved=unresolved))
    targeted_metrics=[load(p) for p in (r/'preservation-temporal-perturbed-A95-EP-1pct').glob('*/temporal-metrics.json')]
    if len(targeted_metrics)!=3:errors.append(dict(input='targeted-perturbation',reason='three timesteps required'))
    write_json(a.output/'perturbation-results.json',dict(targeted_runs=targeted_metrics,status='pass' if len(perturb)==6 and all(t['status']=='pass' for t in perturb) else 'fail',tables=perturb,runs=perturb_runs))
    if len(temporal_runs)!=36 or len(perturb_runs)!=18:errors.append(dict(input='temporal-run-coverage',neutral=len(temporal_runs),perturbed=len(perturb_runs)))
    write_json(a.output/'timing-results.json',dict(cases=timings))
    uncertainty=read('uncertainty','uncertainty-final.json');write_json(a.output/'uncertainty-results.json',uncertainty)
    control=read('controls','continuation-controls/qualification.json');cv=native_check_coverage(control.get('checks',[]),control.get('expected_check_names'))
    if cv['status']!='pass':errors.append(dict(input='continuation-coverage',coverage=cv))
    extras=read('extra-regressions','extra-regressions.json');regression=read('standard-CLS','regression/baseline.json')
    write_json(a.output/'regressions.json',dict(continuation=control,coverage=cv,extras=extras,legacy=regression,
        command_line_failures=read('CLI-failures','command-line-failure-contract.json'),required_numerical_failure=read('synthetic-numerical-failure','deliberate-numerical-check/result.json')))
    production_changes=[]
    for path,digest in manifest.get('source_files',{}).items():
        if Path(path).suffix not in ('.cpp','.h') or path.startswith('library/UnitTests'):continue
        current=r/'source'/path
        if hashlib.sha256(current.read_bytes()).hexdigest()!=digest:production_changes.append(path)
    binary_changed=hashlib.sha256((r/'runtime/IncNavierStokesSolver').read_bytes()).hexdigest()!=manifest.get('solver',{}).get('sha256')
    if production_changes or binary_changed:errors.append(dict(input='production-provenance',changed_sources=production_changes,binary_changed=binary_changed))
    write_json(a.output/'source-build.json',manifest);write_json(a.output/'environment.json',read('environment','environment.json',None));write_json(a.output/'actual-loaded-libraries.json',read('loaded-libraries','runtime/actual-loaded-libraries.json',None))
    status='fail' if errors else 'unresolved' if unresolved else 'pass'
    gates.update(implementation_corrections='pass' if after.get('status')=='pass' and diagnostic.get('status')=='pass' else 'fail',analytical_native_tests='pass' if len(weak)==14 and all(w['status']=='pass' for w in weak) else 'fail',
        spatial=spatial.get('status'),quadrature_and_evaluation=uncertainty.get('status'),temporal_unperturbed='partially_unresolved' if unresolved else 'demonstrated',
        temporal_perturbed='pass' if len(perturb)==6 and all(t['status']=='pass' for t in perturb) else 'fail',restart_cache_rollback=control.get('status'),primary_trajectories_bitwise=beforeafter.get('status'))
    write_json(a.output/'overall-gates.json',dict(schema=1,status=status,exit_code=1 if errors else 2 if unresolved else 0,gates=gates,errors=errors,unresolved=unresolved,
        source_revision=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip(),runtime_source=manifest.get('revision'),compiled_CXX_revision=manifest.get('compiled_CXX_revision'),
        scope='bounded serial flat2D one-sweep, T=.0016; no long-time/physical/curved/3D/MPI/nonlinear qualification',root=str(r)))
    print(status);return 1 if errors else 2 if unresolved else 0
if __name__=='__main__':
    try:code=main()
    except Exception as error:
        if '--output' in sys.argv:
            directory=Path(sys.argv[sys.argv.index('--output')+1]);directory.mkdir(parents=True,exist_ok=True)
            write_json(directory/'collection-failure.json',dict(status='fail',failure=repr(error)))
        print('collection failed: '+repr(error),file=sys.stderr);code=1
    sys.exit(code)
