#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Short d0 RN implementation controls, not experimental preservation."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted7_top import make,request
from ct_accepted_open import run,probe,require,row_value
from ct_accepted_pilot import replay,budget_classification
from ct_accepted_budget_gate import physical_budgets
from ct_accepted6_audit import captured_k_matrix
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import write,set_parameter
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();manifest=RequestManifest(a.output,dict(schema='CT7RNHistoryBudgetControlsV1',runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source=sha(Path(__file__)),dt=1e-6,steps=2,scope='synthetic d0 implementation control only'))
    result=dict(schema='CT7RNHistoryBudgetControlsV1',checks={},models={},physical_experimental_admission=False)
    def evolve(name,model='A95',steps=2,**kwargs):
        folder=a.output/name;c=make(folder,a.runtime,model,steps=steps,**kwargs)
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log')
        path=folder/'case.evolution.json'
        r=json.loads(path.read_text()) if path.exists() else dict(failure=(folder/'solver.log').read_text()[-3000:],time=None)
        return folder,c,r,e
    try:
        full,c,r,e=evolve('full');require(e['exit_code']==0 and not r['failure'] and abs(r['time']-2e-6)<1e-16,'RN short evolution failed')
        first,_,one,e=evolve('first',steps=1);require(e['exit_code']==0,'RN one step failed')
        continued,_,end,e=evolve('continued',steps=1,restart=Path(one['checkpoint']));require(e['exit_code']==0 and same_history(r,end),'RN unchanged-method restart differs')
        result['checks']['restart']=dict(status='pass',state_sha256=exact_vector_hash(end['state']),report_sha256=sha(continued/'case.evolution.json'))
        folder=a.output/'rejected';make(folder,a.runtime,steps=2);root=ET.parse(folder/'case.xml').getroot();set_parameter(root,'RANSWallPimpleHardAuditRejectStep',2);write(root,folder/'case.xml')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');rejected=json.loads((folder/'case.evolution.json').read_text())
        require(e['exit_code']!=0 and 'deliberate' in rejected['failure'] and same_history(one,rejected),'RN rejected step did not roll back complete history')
        result['checks']['rollback']=dict(status='pass',reason=rejected['failure'],report_sha256=sha(folder/'case.evolution.json'))
        folder,_,wrong,e=evolve('old_boundary_restart',steps=0,restart=Path(one['checkpoint']),old_top=True)
        reason=wrong['failure']+(folder/'solver.log').read_text()[-2000:]
        require(e['exit_code']!=0 and 'signature' in reason.lower(),'old top history wrong rejection')
        result['checks']['old_boundary_history']=dict(status='pass',reason=reason)
        shifted,_,s,e=evolve('shifted',level=.125);require(e['exit_code']==0,'RN shifted pressure evolution failed')
        initial,pe=probe(a.runtime,full,request(c),'initial');require(pe['exit_code']==0,'RN initial observation failed')
        off=initial['observation']['offsets'];expected=np.zeros(len(r['state']))
        for q in initial['observation']['pressure_mass']:
            row=q['pressure'];expected[np.asarray(row['ids'],int)+off[4+q['region']]]=.125*np.asarray(row['constant'])
        defect=float(np.max(abs(np.asarray(s['state'])-np.asarray(r['state'])-expected)))
        require(defect<=2e-10,'RN traction-pressure shift changed evolved physical fields')
        result['checks']['pressure_level_evolution']=dict(status='pass',maximum_defect=defect)
        for model in ('A95','A93','E94E95'):
            if model=='A95':folder,construction,report=full,c,r
            else:
                folder,construction,report,e=evolve(model)
                if e['exit_code']!=0 or report['failure']:
                    result['models'][model]=dict(completion='fail',failure=report['failure'],report_sha256=sha(folder/'case.evolution.json'));continue
            ini,pe=probe(a.runtime,folder,request(construction),'budget-initial');require(pe['exit_code']==0,'budget initial reference rejected')
            attempt=report['pimple_attempts'][-1];path=Path(attempt['k_snapshot']);require(sha(path)==attempt['k_snapshot_sha256'],'changed actual frozen snapshot')
            k=json.loads(path.read_text());body=dict(budget_binding=dict(model=model,time=report['time'],dt=report['dt'],native_report_sha256=sha(folder/'case.evolution.json'),state_sha256=exact_vector_hash(report['state']),configuration_sha256=sha(folder/'case.xml')),budgets={})
            for label,orders,face in [('native',(32,8,12,12,12),16),('q48f16',(48,12,16,16,16),16),('q64f16',(64,16,20,20,20),16),('q64f32',(64,16,20,20,20),32)]:
                dest=folder/label;dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
                for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,key,value)
                write(root,dest/'case.xml');req=replay(ini,report['state'],report['time'],[],export_basis=True,export_physical_scalar_budget=True,
                    export_frozen_scalar_reconciliation=True,previous_state=report['accepted_states'][-2],stage_history=k['history'],actual_flow_stage=k['post_flow'],actual_pre_k_stage=k['pre_k'],physical_face_order=face)
                packet,pe=probe(a.runtime,dest,req,'budget');require(pe['exit_code']==0,'RN physical budget probe failed:'+str(packet.get('failure')))
                b=budget_classification(packet,construction);b.update(binding=body['budget_binding'],packet_sha256=sha(dest/'budget.json'),execution=pe)
                top=[v for v in packet['physical_scalar_budget']['physical_faces'] if v['region']==1 and v['axis']==2 and v['outward_side']==1]
                b['actual_top_diffusive_flux']=sum(v['weight']*v['epsilon_frozen_diffusive_outward'] for v in top)
                b['prescribed_top_diffusive_flux']=sum(v['weight']*float(np.linalg.norm(v['area']))*v['epsilon_diffusive_datum'] for v in packet['open_boundary']['top_faces'])
                b['top_flux_scope']='actual physical gradients and prescribed physical data reported separately; no substitution in budgets'
                body['budgets'][label]=b
                if label=='native':body['captured_k_matrix']=captured_k_matrix(packet,k)
            net=0.;o=ini['observation'];state=np.asarray(report['state'])
            for face in o['faces']:
                if face['interface']:continue
                u=np.array([row_value(face['velocity_k'],state,o['offsets'][b])[0] for b in (0,1)])
                net+=face['weight']*float(u@np.asarray(face['area'])[[0,2]])
            scale=construction['ustar']*construction['Z'];body['mass']=dict(outward_net=net,scale=scale,scaled_defect=abs(net)/scale,status='pass' if abs(net)/scale<=1e-8 else 'fail')
            body['budget_gate']=physical_budgets(body);write_json(folder/'physical-budgets.json',body)
            result['models'][model]=dict(completion='pass',time=report['time'],state_sha256=exact_vector_hash(report['state']),report_sha256=sha(folder/'case.evolution.json'),
                budget_gate=body['budget_gate'],mass=body['mass'],captured_k_matrix=body['captured_k_matrix'],physical_budget_receipt_sha256=sha(folder/'physical-budgets.json'))
            write_json(a.output/'progress.json',result)
        result['execution']='complete';result['status']='pass' if all(v['completion']=='pass' and v['budget_gate']['status']=='pass' and v['captured_k_matrix']['status']=='pass' and v['mass']['status']=='pass' for v in result['models'].values()) and len(result['models'])==3 else 'failed_or_incomplete_budget_accuracy'
    except Exception as error:result.update(status='fail',failure=repr(error))
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1

if __name__=='__main__':raise SystemExit(main())
