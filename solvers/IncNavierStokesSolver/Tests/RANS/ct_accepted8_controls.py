#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT8 displaced history/pressure/live-boundary controls on a small fixture."""
import argparse,copy,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted8_native import make,request,read,binding,col
from ct_accepted8_reference import POLICY,PROVIDER,FIT
from ct_accepted_open import run,probe,row_value,require
from ct_accepted_pilot import replay
from ct_accepted4_analysis import weak
from ct_accepted8_reference import reference
from ct_accepted_feasibility import write_json
from validate_a93_a95_e94e95 import oracle
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','reference','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');c=read(a.reference)['case'];c.update(Z=.5,knots=sorted(set([0.,c['H']/4,c['H']/2,c['roughness_input'],col(c).cutoff,c['H'],.3,.5])),degree=3,pressure_degree=1)
    m=RequestManifest(a.output,dict(schema='CT8DisplacedCompatibilityV1',case=c,reference=sha(a.reference),source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver')));result=dict(checks={},status='fail')
    def evolve(label,steps=2,restart=None,change=None):
        folder=a.output/label;make(folder,a.runtime,c,1e-6,steps,restart)
        if change:
            r=ET.parse(folder/'case.xml').getroot();change(r);write(r,folder/'case.xml')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');path=folder/'case.evolution.json'
        r=read(path) if path.exists() else dict(failure=(folder/'solver.log').read_text()[-3000:]);return folder,r,e
    try:
        folder,full,e=evolve('full');require(e['exit_code']==0 and not full['failure'],'displaced smoke rejected')
        first,one,e=evolve('one',1);require(e['exit_code']==0,'first step failed')
        end,continued,e=evolve('continued',1,Path(one['checkpoint']));require(e['exit_code']==0 and same_history(full,continued),'restart differs')
        result['checks']['restart']=dict(status='pass',state_sha256=exact_vector_hash(continued['state']),report_sha256=sha(end/'case.evolution.json'))
        dest,rejected,e=evolve('rollback',change=lambda r:set_parameter(r,'RANSWallPimpleHardAuditRejectStep',2))
        require(e['exit_code']!=0 and 'deliberate' in rejected['failure'] and same_history(one,rejected),'rejected state/history rollback differs')
        result['checks']['rollback']=dict(status='pass',failure=rejected['failure'],report_sha256=sha(dest/'case.evolution.json'))
        dest,wrong,e=evolve('stale_fit',0,Path(one['checkpoint']),lambda r:set_info(r,'RANSOpenFitPolicy','ForeignFitV1'))
        reason=wrong['failure']+(dest/'solver.log').read_text()[-1000:];require(e['exit_code']!=0 and 'signature' in reason.lower(),'stale fit history wrong rejection')
        result['checks']['stale_boundary_history']=dict(status='pass',reason=reason)
        dest,shifted,e=evolve('pressure_level',change=lambda r:set_parameter(r,'RANSOpenPressureLevel',.125));require(e['exit_code']==0,'pressure offset evolution rejected')
        ini,pe=probe(a.runtime,folder,request(c,True),'operator');require(pe['exit_code']==0,'operator fixture unavailable');binding(ini,c)
        off=ini['observation']['offsets'];expected=np.zeros(len(full['state']))
        for q in ini['observation']['pressure_mass']:
            row=q['pressure'];expected[np.array(row['ids'],int)+off[4+q['region']]]=.125*np.asarray(row['constant'])
        defect=float(max(abs(np.asarray(shifted['state'])-np.asarray(full['state'])-expected)));require(defect<2e-10,'absolute traction-pressure invariance')
        result['checks']['pressure_level']=dict(status='pass',maximum_defect=defect)
        for label,modify,word in [('wrong_model',lambda r:set_info(r,'TurbulenceModel','A93'),'Natural A95'),
            ('old_provider',lambda r:set_info(r,'RANSOpenReference','NeutralCompositeV1'),'displacement'),
            ('cache',lambda r:set_info(r,'RANSWallEvaluationReuse','ExactStepV1'),'caching')]:
            dest=a.output/label;make(dest,a.runtime,c,1e-6,0);r=ET.parse(dest/'case.xml').getroot();modify(r);write(r,dest/'case.xml')
            q,e=probe(a.runtime,dest,dict(operation='functions'),'negative');require(e['exit_code']!=0 and word.lower() in q.get('failure','').lower(),'wrong negative:'+label+str(q))
            result['checks'][label]=dict(status='pass',reason=q['failure'])
        fixed={int(i) for i,_ in ini['resolution']['constraints']};state=np.asarray(ini['observation']['state']).copy()
        ids=sorted({off[2]+i for f in ini['open_boundary']['top_faces'] for i,v in zip(f['velocity_k_row']['ids'],f['velocity_k_row']['phi'],strict=True) if abs(v)>1e-12}-fixed)
        state[ids[-1]]+=1e-3*c['K0'];req=replay(ini,state.tolist(),0.,[],export_basis=True)
        q,e=probe(a.runtime,folder,req,'top-k-changed');require(e['exit_code']==0,'admissible top-k change failed')
        old=[f['epsilon_diffusive_datum'] for f in ini['open_boundary']['top_faces']];new=[f['epsilon_diffusive_datum'] for f in q['open_boundary']['top_faces']];require(old==new,'prescribed physical flux changed with diffusivity')
        wrong=[];s=np.asarray(q['observation']['state'])
        for f in q['open_boundary']['top_faces']:
            G=np.zeros((3,3))
            for b,axis in enumerate((0,2)):G[axis]=row_value(f['velocity_k_row'],s,off[b])[1]
            k=row_value(f['velocity_k_row'],s,off[2])[0];eps=row_value(f['epsilon_row'],s,off[3])[0];law=oracle('A95',k,eps,G)
            wrong.append(-col(c).A/(c['Z']-c['displacement'])**2*law['nut']/1.3)
        signal=max(abs(x-y) for x,y in zip(wrong,new,strict=True));require(signal>1e-10,'nondiscriminating wrong evolving-diffusivity flux')
        result['checks']['physical_flux']=dict(status='pass',datum=new[0],wrong_evolving_diffusivity_signal=signal)
        result['checks']['model_history_binding']=dict(status='pass',signature=full['signature'],state_sha256=exact_vector_hash(full['state']),boundary=full['initialization_receipt'],scope='two tiny implementation steps, not strip preservation')
        result['status']='pass';result['execution']='complete'
    except Exception as error:result['failure']=repr(error)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
