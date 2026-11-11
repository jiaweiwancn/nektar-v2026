#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT-Accepted-1 receipts; never relabel a failed physical pilot."""
import argparse,collections,csv,hashlib,json,shutil,subprocess
from pathlib import Path
from rans_qualification import write_json
from rans_request_manifest import sha


def load(path):return json.loads(Path(path).read_text())
def require(ok,reason):
    if not ok:raise ValueError(reason)
def canonical(x):return json.dumps(x,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root;s=r/'source';out=a.output;out.mkdir(parents=True,exist_ok=False)
    files={};mapped={}
    def bound(path):
        path=Path(path);files[str(path)]=sha(path);return load(path)
    def compact(value):
        if isinstance(value,list):return [compact(x) for x in value]
        if not isinstance(value,dict):return value
        result={}
        for k,v in value.items():
            if k=='actual_mapped_files':
                digest=hashlib.sha256(canonical(v)).hexdigest()
                mapped[digest]={p:h for p,h in v.items() if '.so' in Path(p).name or '/runtime-' in p}
                result['actual_mapped_files_digest']=digest
            else:result[k]=compact(v)
        return result
    feas=bound(r/'feasibility-v2/results.json');require(feas['status']=='pass','independent_feasibility_verification_failed')
    for name in ('request-manifest.json','completed-artifacts.json'):
        f=r/'feasibility-v2'/name;files[str(f)]=sha(f);shutil.copy2(f,out/('feasibility-'+name))
    flat=feas['flat'];curved=feas['curved'];require(len(flat)==180 and len({(x['model'],x['surface'],x['scenario'],x['H_m'],x['K_factor']) for x in flat})==180,'flat_coverage')
    def table(name,rows,fields):
        with (out/name).open('x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    common=['model','surface','H_m','H','hr','K','nu','Lambda','cutoff_existence_margin','cutoff','cutoff_over_H','zv','minimum_stretch','classification']
    def base(x):
        row={k:x.get(k) for k in common};row.update(hr_m=.04*x['hr'],q=x['coefficients']['q'],kappa=x['coefficients']['kappa'],
            sigma_epsilon=x['coefficients']['sigma_epsilon'],alpha=x['coefficients']['alpha'],c_mu_wall=x['coefficients']['c_mu_wall'],c_l_wall=x['coefficients']['c_l_wall'],
            minimum_k_over_K=x.get('reference_k',{}).get('minimum_over_K'),minimum_k_location=x.get('reference_k',{}).get('location'),
            Q_end=x.get('Q_end'),production_budget=x['finite_nu_production'],budget_defect=x.get('budget_defect'))
        return row
    rows=[]
    import numpy as np
    from validate_a93_a95_e94e95 import oracle
    for x in flat:
        row=base(x);row.update(scenario=x['scenario'],K_factor=x['K_factor'],J0=x['J0'],J2=x['J2'],
            ustar=x['counterfactual']['u'],matching_velocity=x['matching_velocity'],traction=x['total_reference_traction'],
            unshifted_z0_m=.04*x['unshifted_compatible_z0'],dissipation_plateau=x.get('dissipation_plateau'),
            flat_covariance_margin_over_K=x.get('flat_covariance_margin_over_K'),mechanical_work_defect=x.get('mechanical_work_defect'),
            quadrature_J0_uncertainty=x['quadrature_uncertainty']['J0'],minimum_k_quadrature_uncertainty=x['quadrature_uncertainty'].get('minimum_k_over_K'))
        G=np.zeros((3,3));G[0,2]=x['counterfactual']['u']/(x['coefficients']['kappa']*x['H'])
        law=oracle(x['model'],x['K'],x['counterfactual']['epsilon_amplitude']/x['H'],G)
        row['outer_reference_covariance_margin_over_K']=float(np.linalg.eigvalsh(law['R'])[0]/x['K'])
        rows.append(row)
    table('flat-feasibility.csv',rows,list(rows[0]))
    rows=[]
    for x in curved:
        row=base(x);row.update(x_over_Hh=x['x_over_Hh'],curvature=x['curvature'],geometric_classification=x['geometric_classification'],
            thin_layer_eligible=x['thin_layer_eligible'],within_sufficient_global_reach=x['within_sufficient_global_reach'],
            covariance=x['covariance'],frozen_column_eligible=x['classification']=='admissible_reference' and x['geometric_classification']=='admissible_geometry')
        rows.append(row)
    table('curved-columns.csv',rows,list(rows[0]))
    write_json(out/'feasibility-summary.json',dict(status=feas['status'],flat_count=len(flat),curved_count=len(curved),
        flat_classifications=dict(collections.Counter(x['classification'] for x in flat)),
        curved_classifications=dict(collections.Counter(x['classification'] for x in curved)),
        geometric_classifications=dict(collections.Counter(x['geometric_classification'] for x in curved)),
        incident=feas['incident'],fits=feas['fits'],smooth_8mm_counterexample=feas['smooth_8mm_counterexample'],
        negative_controls=feas['negative_controls'],global_reach_justification=feas['global_reach_justification'],
        native_flat_checks=[dict(model=x['model'],surface=x['surface'],scenario=x['scenario'],H_m=x['H_m'],K_factor=x['K_factor'],
            status=x['native']['status'],expected=x['native']['expected'],defects=x['native'].get('defects'),
            rejection=x['native'].get('native',{}).get('error')) for x in flat],
        decision=feas['decision']))
    receipts={}
    for name,path in dict(boundary_verification=r/'boundary-verification-v4/results.json',
        history_controls=r/'boundary-controls-v1/results.json',physical_budget_controls=r/'boundary-controls-supplement-v1/results.json',
        pilot=r/'pilot-v3/results.json',postmortem=r/'pilot-diagnostics-v1/results.json').items():
        receipts[name]=bound(path);write_json(out/(name+'.json'),compact(receipts[name]))
        for child in ('request-manifest.json','completed-artifacts.json'):
            f=path.parent/child
            if f.exists():files[str(f)]=sha(f);shutil.copy2(f,out/(name+'-'+child))
    witnesses=[]
    for name in ('A95-dt-0.002','A95-dt-0.001'):
        path=r/'pilot-v3'/name/'case.evolution.json';native=bound(path);initial=bound(path.parent.parent/(name+'-zero')/'case.evolution.json')
        require('negative square-root radicand' in native['failure'],'physical_failure_reason_changed')
        event=native['pimple_attempts'][-1]
        witnesses.append(dict(case=name,native_failure=native['failure'],attempted_time=event['time'],last_accepted_time=native['time'],
            accepted_steps=native['steps'],last_accepted_state=native['state'],rejected_state=native['rejected_state'],
            native_history=native['native_history'],native_history_times=native['native_history_times'],checkpoint=native['checkpoint'],
            initial_state=initial['state'],physics_signature_sha256=hashlib.sha256(native['signature'].encode()).hexdigest(),
            initial_boundary_receipt=native['initialization_receipt'],actual_failed_stage=event,
            independent_witness=receipts['postmortem']['witnesses'][name],report_sha256=sha(path),configuration_sha256=sha(path.parent/'case.xml')))
    write_json(out/'a95-failure-witnesses.json',witnesses)
    source=r/'pilot-diagnostics-v1/minimal-timed-reproducer/case.xml'
    # Portable reproduction changes only its output destination. The complete
    # tested original configuration digest remains bound in the postmortem.
    import xml.etree.ElementTree as ET
    root=ET.parse(source).getroot()
    for node in root.findall('./CONDITIONS/SOLVERINFO/I'):
        if node.get('PROPERTY')=='RANSWallCheckpoint':node.set('VALUE','checkpoint')
    ET.indent(root);ET.ElementTree(root).write(out/'minimal-a95-radicand.xml',encoding='utf-8',xml_declaration=True)
    hashes={}
    for name in ('RANSModels','AnalyticalWallFunction','EquilibriumAWFProfiles'):
        for ext in ('.h','.cpp'):
            path='library/SolverUtils/Turbulence/'+name+ext;old=subprocess.check_output(['git','show','b127d307:'+path],cwd=s)
            current=sha(s/path);require(hashlib.sha256(old).hexdigest()==current,'reused_kernel_source_changed:'+path);hashes[path]=current
    runtime=bound(r/'runtime-v2/complete-identity.json');write_json(out/'runtime-identity.json',runtime)
    write_json(out/'environment.json',dict(compiler=subprocess.check_output(['c++','--version'],text=True),
        system=subprocess.check_output(['uname','-a'],text=True),python=subprocess.check_output(['python3','--version'],text=True),
        numerical_libraries=dict(numpy=np.__version__,scipy=__import__('scipy').__version__),
        threads={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS')},
        runtime_ranks=1,compiled_MPI='ON',project_workers=16,compiler_workers=4,aggregate_limit=72,startup_reserve=16))
    for name in ('dataset-audit.json','feasibility-execution-v1.json','feasibility-execution-v2.json',
        'build-execution-v1.json','build-execution-v2.json','build-execution-v3.json','build-execution-v4.json','build-execution-v5.json',
        'native-pipeline-v2-execution.json','native-pipeline-v3-execution.json','pilot-execution-v2.json','postmortem-execution-v1.json','ctest-execution-v1.json','worker-guard-status.json'):
        f=r/name
        if f.exists():files[str(f)]=sha(f);shutil.copy2(f,out/name)
    old_source_scope=s/'docs/plans/a93-a95-e94e95-awf-directional-resolution-qualification-evidence/source-scope.json'
    write_json(out/'source-scope.json',dict(starting_revision='b127d307edb423f93a93aa9d861ce4b6a486573c',
        source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=s,text=True).strip(),
        source_status=subprocess.check_output(['git','status','--porcelain'],cwd=s,text=True),
        authoritative_environment='Linux serial, MPI-enabled build; all numerical libraries one thread',
        actual_runtime=r.as_posix()+'/runtime-v2',runtime_identity_sha256=sha(r/'runtime-v2/complete-identity.json'),
        unchanged_reused_kernel_sources=hashes,historical_directional_source_scope_sha256=sha(old_source_scope),
        historical_decisions='unchanged, no old evolved trajectory reused as open-floor evidence',
        reused='frozen dataset, unchanged kernels, sealed boundary/history controls; no old response or neutral-duration claim transferred',
        new='independent columns, open equations and verification, fresh zero-step/open histories, failed A95 pilot, A93/E94E95 smokes, read-only postmortem'))
    boundary_pass=all(receipts[n]['status']=='pass' for n in ('boundary_verification','history_controls','physical_budget_controls'))
    pilot=receipts['pilot'];require(pilot['status']=='fail','historical_pilot_failure_relabelled')
    decision=dict(schema='CT-Accepted-1-required-scope-v1',boundary_equations='pass' if boundary_pass else 'fail',
        independent_reference_feasibility='pass as a study, with explicit inadmissible physical candidates',
        A95_intended_rough_floor_pilot='failed: actual A95 radicand at post-pressure T=.004, planned T=.2 unreached',
        A93_smoke=pilot['cases']['A93-dt-0.001']['status'],E94E95_smoke=pilot['cases']['E94E95-dt-0.001']['status'],
        A95_common_prefix='diagnostic only, does not qualify planned endpoint or full advection transit',
        diagnostic_postmortem=receipts['postmortem']['status'],experimental_prediction_validation='not_attempted',
        native_curved_or_hill_evolution='not_attempted and guarded',readiness='small conservative flat prototype ready within verified scope; physical A95 precursor/benchmark not admitted',
        next_gate='diagnose finite-depth/delta constitutive-domain loss with signed pressure/strain and inlet/normal representation evidence before longer unforced development; use this small case to specify scalable conservative backend and curved coupling',
        historical_failures_and_eight_unresolved_epsilon_orders='unchanged')
    write_json(out/'decision.json',decision);write_json(out/'actual-loaded-library-maps.json',mapped)
    write_json(out/'raw-artifact-hashes.json',files)
    write_json(out/'hashes.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()})
    print(json.dumps(decision,indent=2));return 0 if boundary_pass and receipts['postmortem']['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
