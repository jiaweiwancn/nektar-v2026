#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Physical flux quadrature and finite defects, separate from row reactions."""
import argparse
from pathlib import Path
from qualify_rans_preservation import load
from rans_duration_contract import require, MEMBERS
from rans_request_manifest import RequestManifest, sha
from rans_qualification import write_json
from rans_preservation_oracle import Reference


def whole_patch(receipt, ref):
    """Exclude the internal face from the physical external-boundary budget."""
    b=receipt['physical_budget'];parts=[b['equations'][name] for name in ('inner_k','outer_k')]
    interface=[q for q in b['physical_faces'] if q['axis']==2 and abs(q['coordinate'][2]-ref.H)<1e-13]
    outward={kind:sum(q['weight']*q['k_'+kind+'_flux'] for q in interface) for kind in ('frozen','current')}
    values={key:sum(p[key] for p in parts) for key in ('volume','storage_backward_difference','stage_source','instantaneous_source')}
    values['stage_physical_boundary_flux']=sum(p['stage_physical_boundary_flux'] for p in parts)-outward['frozen']
    values['instantaneous_physical_boundary_flux']=sum(p['instantaneous_physical_boundary_flux'] for p in parts)-outward['current']
    values['frozen_stage_physical_defect']=values['storage_backward_difference']-values['stage_source']-values['stage_physical_boundary_flux']
    values['instantaneous_physical_defect']=values['storage_backward_difference']-values['instantaneous_source']-values['instantaneous_physical_boundary_flux']
    values['frozen_operator_splitting']=values['instantaneous_physical_defect']-values['frozen_stage_physical_defect']
    values['interface_frozen_outward_flux_sum']=outward['frozen'];values['interface_current_outward_flux_sum']=outward['current']
    scale=ref.K*ref.us/ref.H*values['volume']
    return {key:value/scale for key,value in values.items() if key!='volume'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--cohort',required=True)
    p.add_argument('--time',type=float,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=a.root.resolve();receipts={};inputs={}
    for kind in ('perturbed','baseline'):
        for member in MEMBERS:
            versions=((1,16),(2,16),(2,32)) if member in MEMBERS[:2] else ((2,32),)
            for q,face in versions:
                label=a.cohort+'-'+kind+'-'+member+'-'+str(a.time)+'-Q'+str(q)+('-F32' if face==32 else '')
                path=r/'physical-budgets'/label;request=load(path/'request-manifest.json')['request']
                require(RequestManifest(path,request,True).completed(path),'physical_budget_unsealed:'+label)
                data=load(path/'results.json');require(data['status']=='pass','physical_budget_execution_failed:'+label)
                require(abs(data['time']-a.time)<1e-13 and data['physical_budget']['face_order']==face,
                        'physical_budget_requested_time_or_face_order_mismatch')
                receipts[kind,member,q,face]=data;inputs[label]=dict(source=str(path),result_sha256=sha(path/'results.json'),request_sha256=sha(path/'request-manifest.json'))
    cfg=next(iter(receipts.values()))['configuration']
    require(all(v['configuration']==cfg for v in receipts.values()),'physical_budget_configuration_mismatch')
    require(all(v['runtime']==next(iter(receipts.values()))['runtime'] for v in receipts.values()),
            'physical_budget_runtime_mismatch')
    ref=Reference(cfg['model'],cfg['policy'],cfg['wall'])
    for value in receipts.values():value['scaled']=dict(value['scaled'],whole_patch_k=whole_patch(value,ref))
    intent=dict(schema='PhysicalBudgetAnalysisV1',configuration=cfg,time=a.time,inputs=inputs,
                quadrature_absolute_floor=1e-9,quadrature_spatial_fraction=.01,
                temporal_absolute_floor=1e-9,temporal_spatial_fraction=.1,
                spatial_maximum_ratio=.85,residual_floor=5e-9)
    seal=RequestManifest(a.output,intent,False);quadrature=[];defects=[];responses={}
    for kind in ('perturbed','baseline'):
        for member in MEMBERS[:2]:
            versions=[receipts[kind,member,q,f] for q,f in ((1,16),(2,16),(2,32))]
            require(len({d['state_sha256'] for d in versions})==1 and
                    len({d['previous_state_sha256'] for d in versions})==1,
                    'physical_budget_frozen_coefficients_differ')
            other='B11' if member=='B5' else 'B5';comparison=receipts[kind,other,2,32]
            for equation,values in versions[2]['scaled'].items():
                for term,value in values.items():
                    gap=abs(value-comparison['scaled'][equation][term]);limit=max(1e-9,.01*gap)
                    changes=dict(volume=abs(versions[0]['scaled'][equation][term]-versions[1]['scaled'][equation][term]),
                                 face=abs(versions[1]['scaled'][equation][term]-value))
                    quadrature.append(dict(kind=kind,member=member,equation=equation,term=term,
                        changes=changes,limit=limit,status='pass' if max(changes.values())<=limit else 'unresolved'))
        for equation in ('inner_k','outer_k','outer_epsilon','whole_patch_k'):
            key='frozen_stage_physical_defect';values=[receipts[kind,m,2,32]['scaled'][equation][key] for m in MEMBERS]
            c,f=map(abs,values[:2]);gap=abs(values[0]-values[1]);delta=max(abs(values[0]-values[2]),abs(values[1]-values[3]))
            uncertainty=max((max(q['changes'].values()) for q in quadrature if q['kind']==kind and q['equation']==equation and q['term']==key),default=0.)
            floor=max(5e-9,5*uncertainty)
            spatial='unresolved_floor' if max(c,f)<=floor else 'pass' if f<=.85*c and c-f>5*uncertainty else 'failed_or_unresolved'
            tlimit=max(1e-9,.1*gap)
            defects.append(dict(kind=kind,equation=equation,signed_values=dict(zip(MEMBERS,values,strict=True)),
                fine_over_coarse=f/c if c else None,spatial_status=spatial,residual_floor=floor,
                quadrature_uncertainty=uncertainty,timestep_difference=delta,temporal_limit=tlimit,
                temporal_status='pass' if delta<=tlimit else 'unresolved',
                frozen_operator_splitting={m:receipts[kind,m,2,32]['scaled'][equation]['frozen_operator_splitting'] for m in MEMBERS}))
    for equation in ('inner_k','outer_k','outer_epsilon','whole_patch_k'):
        responses[equation]={m:{term:value-receipts['baseline',m,2,32]['scaled'][equation][term]
                               for term,value in receipts['perturbed',m,2,32]['scaled'][equation].items()} for m in MEMBERS}
    result=dict(status='pass' if all(q['status']=='pass' for q in quadrature) else 'unresolved',
        execution_status='completed',reconstruction_status='pass',configuration=cfg,time=a.time,
        quadrature_status='pass' if all(q['status']=='pass' for q in quadrature) else 'unresolved',
        quadrature=quadrature,finite_defects=defects,baseline_subtracted_budgets=responses,
        conservation_status='finite physical defects, with spatial and temporal classifications; no exact conservation theorem',
        scope='actual accepted storage and frozen-stage source, independently integrated physical face fluxes; no prescribed-row reaction used, dependent epsilon boundary storage retained; pairs bound uncertainty, not temporal order')
    write_json(a.output/'results.json',result);seal.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
