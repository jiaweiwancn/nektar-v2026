#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable source-only comparisons; never lift/repair selected coefficients."""
import argparse
import math
from pathlib import Path
import shutil
from run_cls_epawf_picard import load,run,require,write_json,state_hash,sha
from run_cls_epawf_pimple import identity
from analyse_cls_epawf_constitutive_v2 import read_report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--upstream-root',type=Path,required=True)
    p.add_argument('--uniform-initial',type=Path);p.add_argument('--uniform-campaign',type=Path)
    p.add_argument('--uniform-steps',type=int,nargs='*',default=[8,16,32,160])
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    summary=dict(schema=1,probe=identity(a.probe.resolve()),states={},counterfactual_scope='same selected coefficient vector, geometry, time, sampled K, reference inputs; source only, not changed boundary/matching data')
    def replay(name,state,time,xml):
        folder=a.output/name;folder.mkdir();shutil.copy2(xml,folder/'case.xml')
        request=dict(operation='source_balance',state=state,history=state,time=time,alpha=0.,
            inner_stress_model='SaturatedR2',require_history_h0=False,
            stations=[-3.5,-2.625,-1.,1.60],heights=[0.,1e-6,1e-5,5e-5,1e-4,2e-4,5e-4,.001,.0025,.005,.01,.02,.05,.1,.25,.5,1.,2.])
        write_json(folder/'input.json',request)
        receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],folder,'probe.log',7200)
        write_json(folder/'receipt.json',receipt);require(receipt['exit_code']==0,'source readback failed: '+name)
        r=load(folder/'result.json');require(r['retained_state_unchanged'] and r['state_sha256']==state_hash(state),'readback changed state')
        valid=[v for v in r['points'] if v.get('valid')]
        require(all(abs(v['D_pre_identity_defect'])<1e-10*max(1.,abs(v['D_pre'])) for v in valid),'pre-destruction identity')
        require(all(abs(v['budget_identity_defect'])<1e-10*max(1.,abs(v['sum_w_D'])) for v in r['ray_budgets']),'ray normalization identity')
        points=[v for v in r['points'] if not v.get('valid') or (v.get('z') in (0.,.0001,.0025,.01,1.,2.))]
        compact=dict(time=time,state_sha256=state_hash(state),wall_policy=r['wall_policy'],input_sha256=sha(folder/'input.json'),
            result_sha256=sha(folder/'result.json'),xml_sha256=sha(folder/'case.xml'),points=points,ray_budgets=r['ray_budgets'],
            native_mass_solve_defect=r['native_mass_solve_defect'],component_load_sum_defect=r.get('component_load_sum_defect'),
            component_mass_rate_sum_defect=r.get('component_mass_rate_sum_defect'),diagnostic_mass_inverse_applications=r['diagnostic_mass_inverse_applications'])
        compact['fixed_native_quadrature_source_comparisons']=[{k:v for k,v in w.items() if k not in ('weak_source_residual','mass_inverted_source_rate')} for w in r['fixed_native_quadrature_source_comparisons']]
        summary['states'][name]=compact;write_json(a.output/'summary.partial.json',summary)
        return r
    r=a.upstream_root.resolve()
    for name,path in [('upstream-h0',r/'analysis16/initial'),('upstream-T.02',r/'monitor/step-000032/step32'),('upstream-T.045',r/'monitor/step-000072/step72')]:
        q=load(path/'input.json');replay(name,q['state'],q['time'],path/'case.xml')
    h=summary['states']['upstream-h0'];ground=[v for v in h['points'] if v.get('valid') and v['x']==1.6 and v['z']==0]
    for v in ground:
        expected_ze=2*v['c_length']*v['molecular_nu']/math.sqrt(v['sampled_K'])
        expected_ea=v['sampled_K']**1.5/(v['c_length']*expected_ze)
        require(abs(v['z_epsilon']/expected_ze-1)<1e-14 and abs(v['epsilon_A']/expected_ea-1)<1e-14,'original AWF cutoff/source mismatch')
        require(abs(v['local_k']/v['sampled_K']-1)<1e-14 and abs(v['limiter']-1)<1e-14 and abs(v['lambda']-1)<1e-14,'h0 comparison is not constant-k unit normalization')
    summary['initial_downstream_ground']=ground
    if a.uniform_initial:
        report=load(a.uniform_initial);replay('uniform-h0',report['accepted_states'][0],0.,a.uniform_initial.parent/'case.xml')
    if a.uniform_campaign:
        root=a.uniform_campaign.resolve()/'trajectory'
        rows={int(v['global_step']):v for v in map(__import__('json').loads,(root/'steps.jsonl').read_text().splitlines())}
        for n in a.uniform_steps:
            if n not in rows:continue
            row=rows[n];report=read_report(Path(row['segment']));start=int(report['steps'])-len(report['accepted_states'])+1
            replay('uniform-step'+str(n),report['accepted_states'][n-start],row['time'],Path(row['segment'])/'case.xml')
    summary['status']='read_only_source_diagnosis_completed';write_json(a.output/'summary.json',summary)
    print(a.output/'summary.json')


if __name__=='__main__':main()
