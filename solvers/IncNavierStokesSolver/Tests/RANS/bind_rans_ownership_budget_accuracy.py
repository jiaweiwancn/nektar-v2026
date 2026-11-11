#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind independent physical-budget accuracy to the exact corrected response."""
import argparse
from pathlib import Path
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest, sha
from rans_duration_contract import require, identity, MEMBERS, prefix
from rans_qualification import write_json
from rans_ownership_contract import budget_gate, OPERATOR


def bind(root,admission,cohort,dt,time,analysis_path,accuracy_path,output,baseline_cohort=None):
    r=root.resolve();ad=load(admission);cfg=ad['configuration'];require(cfg['model']=='A95' and cfg['wall']=='Accepted','budget_scope_not_full_A95_accepted_response')
    analysis=load(analysis_path);accuracy=load(accuracy_path);members={};trajectories={}
    for p in (analysis_path,accuracy_path):
        req=load(p.parent/'request-manifest.json')['request'];require(RequestManifest(p.parent,req,True).completed(p.parent),'unsealed_budget_source')
    for member in MEMBERS:
        members[member]={}
        for kind in ('baseline','perturbed'):
            label=((baseline_cohort or prefix(cfg['model'],cfg['policy'],cfg['wall'])+'-duration-'+str(dt)) if kind=='baseline' else cohort)+'-'+member
            folder=r/'campaigns'/label;rm=load(folder/'request-manifest.json');q=rm['request'];identity(q,cfg)
            require(q['numerical_operator']==OPERATOR,'old_operator_budget_history')
            if kind=='perturbed':require(q['velocity_amplitude']==2e-5,'foreign_budget_response_amplitude')
            else:require(q['velocity_amplitude']==0.,'perturbed_neutral_baseline')
            endpoint=folder/('T-'+format(time,'.4f'));require(RequestManifest(folder,q,True).completed(endpoint),'unsealed_budget_trajectory')
            report=load(endpoint/'case.evolution.json');require(report['numerical_operator']==OPERATOR,'old_operator_native_report')
            members[member][kind]=dict(trajectory_request_sha256=rm['sha256'],source_report_sha256=sha(endpoint/'case.evolution.json'))
            trajectories[kind,member]=dict(request=q,report=report,**members[member][kind])
    controls=set();inputs=load(analysis_path.parent/'request-manifest.json')['request']['inputs']
    expected={(kind,m,q,f) for kind in ('baseline','perturbed') for m in MEMBERS
              for q,f in (((1,16),(2,16),(2,32)) if m in MEMBERS[:2] else ((2,32),))}
    observed=set()
    for label,item in inputs.items():
        path=Path(item['source']);require(sha(path/'results.json')==item['result_sha256'] and sha(path/'request-manifest.json')==item['request_sha256'],'physical_budget_input_changed')
        req=load(path/'request-manifest.json')['request'];require(RequestManifest(path,req,True).completed(path),'unsealed_budget_leaf')
        identity(req['configuration'],cfg);require(abs(req['time']-time)<1e-13,'budget_leaf_horizon')
        found=[(kind,m) for (kind,m),v in trajectories.items() if all(req.get(k)==v[k] for k in ('trajectory_request_sha256','source_report_sha256'))]
        require(len(found)==1,'foreign_or_stale_budget_native_report')
        order=req['orders'][0];face=req['face_order'];require(order in (48,64) and face in (16,32),'foreign_budget_integration_control')
        key=(*found[0],1 if order==48 else 2,face);require(key not in observed,'duplicate_budget_leaf');observed.add(key);controls.add((order,face))
    require(observed==expected,'budget_leaf_companion_coverage')
    binding=dict(schema='NativeOwnershipBudgetBindingV1',numerical_operator=OPERATOR,configuration=cfg,
        cohort=cohort,amplitude=2e-5,time=time,members=members,analysis_sha256=sha(analysis_path),accuracy_sha256=sha(accuracy_path),
        controls=dict(volume_orders=[48,64],physical_face_orders=[16,32]))
    if baseline_cohort:
        binding.update(package_scope='DirectionalResolutionExperimentV1',baseline_cohort=baseline_cohort,
            member_spaces={m:{kind:value['request'].get('spatial_control') for (kind,mm),value in trajectories.items() if mm==m} for m in MEMBERS})
    gate=budget_gate(analysis,accuracy,binding,dict(configuration=cfg,cohort=cohort,time=time,members=members,
        analysis_sha256=sha(analysis_path),accuracy_sha256=sha(accuracy_path)))
    seal=RequestManifest(output,binding,False);value=dict(binding=binding,**gate)
    write_json(output/'results.json',value);seal.seal(output,['results.json']);return value


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('root','admission','analysis','accuracy','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--cohort',required=True);p.add_argument('--baseline-cohort');p.add_argument('--dt',type=float,required=True);p.add_argument('--time',type=float,default=.016)
    a=p.parse_args();value=bind(a.root,a.admission,a.cohort,a.dt,a.time,a.analysis,a.accuracy,a.output,a.baseline_cohort)
    print(value['status']);return 0 if value['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
