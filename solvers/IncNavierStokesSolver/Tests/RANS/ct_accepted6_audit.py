#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent physical scalar budgets and common physical CT6 observations.

Frozen-coefficient controls change integration only. No history is resumed,
signature edited, reaction substituted for flux, or state reset here.
"""
import argparse,json,math,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from ct_accepted6_cases import X,RADIAL,QUADRATURE
from ct_accepted6_operator import geometry
from ct_accepted_feasibility import parameters
from ct_accepted_open import probe,require
from ct_accepted_pilot import replay,budget_classification
from ct_accepted_budget_gate import physical_budgets
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha,canonical
import hashlib

CONTROLS=[('native',QUADRATURE,16),('q48f16',[48,12,16,16,16],16),
          ('q64f16',[64,16,20,20,20],16),('q64f32',[64,16,20,20,20],32)]


def common_grid(order=6):
    g,w=leggauss(order);points=[];weights=[]
    for a,b in zip(X[:-1],X[1:],strict=True):
        for c,d in zip(RADIAL[:-1],RADIAL[1:],strict=True):
            for i in range(order):
                for j in range(order):
                    x=(a+b)/2+(b-a)*g[i]/2;r=(c+d)/2+(d-c)*g[j]/2
                    metric,J,*_=geometry(x,r)
                    points.append([x,0.,r]);weights.append((b-a)*(d-c)*w[i]*w[j]*metric*J/4)
    n=len(points);sides=[0]*n
    for x in (-1.2,-1.,-.8,-.6,-.4,-.2,0.,.2,.4,.6,.8,1.,1.2):
        for side in (-1,1):points.append([x,0.,.5]);sides.append(side)
    return points,sides,weights,n


def budget(packet,model):
    p=parameters(model);o=packet['observation'];volume=math.fsum(q['weight'] for q in o['volume'])
    c=dict(model=model,K0=p['q']*.16,H=.5,x=X,knots=RADIAL,
           reference_column=dict(counterfactual=dict(epsilon_amplitude=.4**3/p['kappa'])))
    result=budget_classification(packet,c)
    for name,row in result['equations'].items():
        scale=(.4**3/(p['kappa']*.5) if name=='outer_epsilon' else p['q']*.16)*volume
        row.update(physical_scale=scale,scaled_frozen_defect=abs(row['frozen_stage_physical_defect'])/scale,
                   status='pass' if abs(row['frozen_stage_physical_defect'])/scale<=.02 else 'fail')
    result['physical_volume']=volume
    result['accuracy_status']='pass' if all(q['status']=='pass' for q in result['equations'].values()) else 'fail'
    result['scope']+='; curved physical volume is independently integrated a*J, not chart area'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('runtime','member','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();receipt=json.loads((a.member/'results.json').read_text())
    model=receipt['model'];initial=json.loads((a.member/'zero/initial.json').read_text())
    folder=a.member/'trajectory';report=json.loads((folder/'case.evolution.json').read_text())
    require(report['steps']>0,'no accepted step available for physical storage audit')
    points,sides,weights,n=common_grid()
    m=RequestManifest(a.output,dict(schema='CT6ActualStatePhysicalAuditV1',member_sha256=sha(a.member/'results.json'),
        native_report_sha256=sha(folder/'case.evolution.json'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),
        state_sha256=exact_vector_hash(report['state']),source_sha256=sha(Path(__file__)),controls=CONTROLS,
        physical_grid_sha256=hashlib.sha256(canonical([points,sides,weights])).hexdigest(),pressure='absolute Pi'))
    binding=dict(model=model,time=report['time'],dt=report['dt'],state_sha256=exact_vector_hash(report['state']),
                 native_report_sha256=sha(folder/'case.evolution.json'),configuration_sha256=sha(folder/'case.xml'))
    result=dict(schema='CT6ActualStatePhysicalAuditV1',budget_binding=binding,budgets={},
                physical_completion=receipt['physical_completion'],time=report['time'])
    i=next(i for i,t in enumerate(report['audit_times']) if abs(t-report['time'])<1e-13)
    require(i>0,'missing consecutive actual accepted states')
    previous=report['accepted_states'][i-1]
    attempt=next(q for q in reversed(report['pimple_attempts']) if q['scheduled_complete'])
    require(sha(attempt['k_snapshot'])==attempt['k_snapshot_sha256'],'changed actual pre-k frozen stage snapshot')
    snapshot=json.loads(Path(attempt['k_snapshot']).read_text())
    for label,orders,face in CONTROLS:
        sub=a.output/label;sub.mkdir();root=ET.parse(folder/'case.xml').getroot()
        for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,key,value)
        write(root,sub/'case.xml')
        req=replay(initial,report['state'],report['time'],[],export_basis=True,
                   export_physical_scalar_budget=True,previous_state=previous,stage_history=snapshot['history'],
                   export_frozen_scalar_reconciliation=True,physical_face_order=face)
        packet,e=probe(a.runtime,sub,req,'budget')
        require(e['exit_code']==0,'physical budget probe rejected '+label+':'+str(packet.get('failure')))
        o=packet['observation'];require(o['state']==report['state'] and o['identity']==initial['observation']['identity'],'frozen control changed effective state or verification identity')
        result['budgets'][label]=dict(budget(packet,model),binding=binding,execution=e,packet_sha256=sha(sub/'budget.json'),
            initialization_receipt=packet['initialization_receipt'] if 'initialization_receipt' in packet else None)
        write_json(a.output/'progress.json',result);print('CT6 physical budget',model,label,result['budgets'][label]['accuracy_status'],flush=True)
    result['budget_gate']=physical_budgets(result)
    sub=a.output/'common';sub.mkdir();write(ET.parse(folder/'case.xml').getroot(),sub/'case.xml')
    req=replay(initial,report['state'],report['time'],points,grid_sides=sides,
               export_test_basis=True,export_constitutive=True,export_manufactured=True)
    packet,e=probe(a.runtime,sub,req,'endpoint');require(e['exit_code']==0,'common physical observation probe failed')
    o=packet['observation'];g=np.asarray(o['grid_gradients'][:n]);w=np.asarray(weights)
    div=g[:,0,0]+g[:,1,2];mass=0.;by_face={}
    from ct_accepted_open import row_value
    for f in o['faces']:
        if f['interface']:continue
        U=np.array([row_value(f['velocity_k'],o['state'],o['offsets'][0])[0],0.,row_value(f['velocity_k'],o['state'],o['offsets'][1])[0]])
        load=f['weight']*float(U@f['area']);mass+=load
        key=str((f['axis'],f['side'])) if 'axis' in f else str(f['coordinate'])
        by_face[key]=by_face.get(key,0.)+load
    # Independently compare the scalar storage/source/face budget and continuity.
    mass_scale=.2*2.4;result['mass']=dict(outward_net=mass,scale=mass_scale,scaled_defect=abs(mass)/mass_scale,
        status='pass' if abs(mass)/mass_scale<=1e-8 else 'fail',by_face=by_face)
    result['common_observations']=dict(path=str(sub/'endpoint.json'),sha256=sha(sub/'endpoint.json'),execution=e,
        count=n,weights=weights,grid_sha256=hashlib.sha256(canonical([points,sides,weights])).hexdigest(),
        physical_volume=float(w.sum()),divergence_RMS=float(np.sqrt(w@(div*div)/w.sum())),divergence_maximum=float(max(abs(div))))
    result['status']='pass' if result['budget_gate']['status']==result['mass']['status']=='pass' else 'fail'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
