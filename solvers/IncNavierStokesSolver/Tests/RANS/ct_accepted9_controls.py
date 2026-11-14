#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded dense/CSR all-row, consistent-mass, history and identity controls."""
import argparse,copy,json,os,resource,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.linalg import block_diag
import ct_accepted8_native as n
from ct_accepted7_top import make as make_zero,request as zero_request,classify
from ct_accepted_open import probe,run,require
from ct_accepted_pilot import replay
from ct_accepted3_pressure import pressure_action
from ct_accepted_feasibility import write_json
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity

BACKEND='AcceptedSerialSparseV1'

def enable(folder,zero=False):
    root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallLinearBackend',BACKEND);set_parameter(root,'FlexibleGMRES',1)
    if zero:set_info(root,'RANSWallSparseVerification','ZeroDisplacementV1')
    write(root,folder/'case.xml')

def maximum(x):return float(np.max(np.abs(x),initial=0.))

def compare_packets(dense,sparse,c):
    """All rows/actions, not a selection of modes or one random direction."""
    a=dense['open_frozen_systems'];b=sparse['open_frozen_systems'];checks={}
    for key in ('raw_residual','constrained_residual','external_constraints'):
        checks[key]=maximum(np.asarray(a[key])-np.asarray(b[key]))
    require(dense['observation']['offsets']==sparse['observation']['offsets'],'field numbering changed')
    require(dense['observation']['state']==sparse['observation']['state'],'initial coefficients changed')
    for name in ('volume','faces'):
        require(dense['observation'][name]==sparse['observation'][name],'native signed row/owner/side inventory changed:'+name)
    for i,(x,y) in enumerate(zip(a['masses'],b['masses'],strict=True)):checks['mass_'+str(i)]=maximum(np.asarray(x)-y)
    for x,y in zip(a['blocks'],b['blocks'],strict=True):
        require(x['block']==y['block'] and x['rows']==y['rows'],'free mask changed')
        for key in ('matrix','residual'):checks['block_'+str(x['block'])+'_'+key]=maximum(np.asarray(x[key])-y[key])
    for key in ('traction_load','diffusive_load','advection_load'):
        checks[key]=maximum(np.asarray(dense['open_boundary'][key])-sparse['open_boundary'][key])
    require(max(checks.values())<=2e-10,'dense/CSR complete row mismatch:'+str(checks))
    # Independent direct linear algebra checks the Nektar matrix-free GMRES
    # operation and both consistent-mass applications at alpha=1.
    packet=sparse['sparse_backend'];flow=a['blocks'][0];ids=np.asarray(flow['rows'],int);nv=packet['velocity_rows'];F=np.asarray(flow['matrix']);M=block_diag(*a['masses'],np.ones((1,1)))[np.ix_(ids[:nv],ids[:nv])]
    rm=np.asarray(flow['residual'][:nv]);rc=np.asarray(flow['residual'][nv:]);A=F[:nv,:nv];G=F[:nv,nv:];D=F[nv:,:nv];C=F[nv:,nv:]
    predictor=np.linalg.solve(A,-rm);bg=np.linalg.solve(M,G);S=D@bg-C;pressure=np.linalg.solve(S,D@predictor+rc);correction=bg@pressure
    checks.update(predictor=maximum(predictor-packet['predictor']),pressure_increment=maximum(pressure-packet['pressure_increment']),consistent_mass_correction=maximum(correction-packet['consistent_mass_correction']))
    require(max(checks.values())<=2e-10,'independent consistent-mass split correction mismatch')
    # Deliberately wrong implementations must be discriminated, not blessed
    # merely because the continuous equilibrium itself is close to balanced.
    lumped=D@(G/np.diag(M)[:,None])-C
    o=dense['observation'];off=o['offsets'];samples=o['pressure_mass'];V=np.zeros((len(samples),off[6]-off[4]));px=[]
    for i,q in enumerate(samples):
        pr=q['pressure'];V[i,np.asarray(pr['ids'],int)+off[4+q['region']]-off[4]]=pr['phi'];px.append(q['coordinate'][0])
    coefficient,_,rank,_=np.linalg.lstsq(V,px,rcond=None);require(rank==V.shape[1] and maximum(V@coefficient-px)<2e-10,'manufactured pressure x projection')
    direction=np.zeros(len(o['state']));direction[off[4]:off[6]]=coefficient
    _,rough,_=pressure_action(o,direction,c['roughness_input'],c['pressure_degree']);rough_signal=maximum(rough)
    # Pressure/roughness coverage is already independently reconstructed by
    # classify() for ALL pressure columns; compare omission to its native G.
    wrong_sign=maximum(S-(-D@bg-C));wrong_mass=maximum(S-lumped)
    stale_mask=maximum(F[0]-F[-1])
    require(min(wrong_sign,wrong_mass,stale_mask,rough_signal)>1e-9,'nondiscriminating pressure/sign/mass/free-mask controls')
    for x,y in zip(a['blocks'][1:],packet['scalars'],strict=True):
        require(x['block']==y['block'],'scalar order changed');increment=np.linalg.solve(np.asarray(x['matrix']),-np.asarray(x['residual']))
        checks['scalar_increment_'+str(x['block'])]=maximum(increment-y['diagnostic_increment'])
    require(max(checks.values())<=2e-10,'scalar inverse mismatch')
    return dict(status='pass',defects=checks,negative_roughness_pressure_omission_signal=rough_signal,negative_sign_signal=wrong_sign,negative_mass_replacement_signal=wrong_mass,negative_stale_mask_signal=stale_mask,
        flow_rows=len(ids),free_velocity=nv,free_pressure=len(ids)-nv,pressure_C=maximum(C),pressure_true_residual=packet['true_pressure_residual'],pressure_iterations=packet['pressure_iterations'],
        scope='complete small matrices/rows and independent consistent-mass linear correction; no governing extra solve')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','ct8','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    pre=n.read(a.ct8/'preflight-v6/results.json');c=copy.deepcopy(pre['case']);c.update(Z=.5,knots=sorted(set([0.,c['H']/4,c['H']/2,n.col(c).cutoff,c['roughness_input'],c['H'],.3,.5])))
    m=RequestManifest(a.output,dict(schema='CT9BackendControlsV1',source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),case=c,operator_tolerance=2e-10,trajectory_tolerance=1e-8))
    result=dict(schema='CT9BackendControlsV1',checks={},status='incomplete')
    def check(name,fn):
        try:result['checks'][name]=fn()
        except Exception as error:result['checks'][name]=dict(status='fail',failure=repr(error))
        write_json(a.output/'progress.json',result);print('CT9 control',name,result['checks'][name]['status'],flush=True)
    def operator(zero,degree):
        cc=copy.deepcopy(c);cc.update(degree=degree,pressure_degree=degree-2);packets=[];paths=[]
        for sparse in (False,True):
            dest=a.output/(('d0' if zero else 'displaced')+'-p'+str(degree)+('-sparse' if sparse else '-dense'))
            if zero:cc=make_zero(dest,a.runtime,degree=degree,interface_order=64);req=zero_request(cc)
            else:n.make(dest,a.runtime,cc,.002,0);req=n.request(cc,True)
            if sparse:enable(dest,zero);req.update(export_sparse_backend=True,sparse_alpha=1.)
            q,e=probe(a.runtime,dest,req,'operator');require(e['exit_code']==0,'native operator failed:'+str(q.get('failure')))
            independent=classify(q,cc,**({} if zero else dict(column=n.col(cc),reference_function=n.reference)))
            packets.append(q);paths.append(dict(path=str(dest/'operator.json'),sha256=sha(dest/'operator.json'),independent=independent,execution=e))
        return dict(compare_packets(*packets,cc),artifacts=paths)
    for zero,degree in ((True,3),(True,4),(False,3),(False,4)):check(('d0' if zero else 'displaced')+'-p'+str(degree),lambda z=zero,d=degree:operator(z,d))
    c.update(degree=3,pressure_degree=1)
    def evolve(label,steps=2,restart=None,modify=None,sparse=True):
        dest=a.output/label;n.make(dest,a.runtime,c,1e-6,steps,restart)
        if sparse:enable(dest)
        if modify:
            root=ET.parse(dest/'case.xml').getroot();modify(root);write(root,dest/'case.xml')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log');file=dest/'case.evolution.json'
        q=n.read(file) if file.exists() else dict(failure=(dest/'solver.log').read_text()[-3000:]);return dest,q,e
    def histories():
        dest,full,e=evolve('full');require(e['exit_code']==0 and not full['failure'],'sparse two-step failure:'+full.get('failure',''))
        first,one,e=evolve('one',1);require(e['exit_code']==0,'sparse first step failure')
        end,continued,e=evolve('continued',1,Path(one['checkpoint']));require(e['exit_code']==0 and same_history(full,continued),'same-backend history mismatch')
        failed,rejected,e=evolve('rollback',modify=lambda root:set_parameter(root,'RANSWallPimpleHardAuditRejectStep',2));require(e['exit_code']!=0 and 'deliberate' in rejected['failure'] and same_history(one,rejected),'rollback mismatch')
        dense,reference,e=evolve('dense-full',sparse=False);require(e['exit_code']==0,'dense rebuild smoke failed');defect=maximum(np.asarray(full['state'])-reference['state']);require(defect<1e-10,'backend tiny-trajectory mismatch')
        legacy=n.read(a.ct8/'controls-v2/full/case.evolution.json');require(same_history(reference,legacy),'rebuilt default dense route changed sealed CT8 state/history')
        negative=[]
        for name,modify,word,restart in [('backend_history',None,'signature',Path(one['checkpoint'])),('stale_fit',lambda root:set_info(root,'RANSOpenFitPolicy','ForeignFitV1'),'signature',Path(one['checkpoint'])),('wrong_model',lambda root:set_info(root,'TurbulenceModel','A93'),'A95',None),('unsupported_cache',lambda root:set_info(root,'RANSWallEvaluationReuse','ExactStepV1'),'caching',None)]:
            dest,q,e=evolve(name,0,restart,modify,sparse=name!='backend_history');message=q.get('failure','')+(dest/'solver.log').read_text()[-4000:];require(e['exit_code']!=0 and word.lower() in message.lower(),'wrong negative outcome:'+name+message);negative.append(dict(name=name,status='pass',reason=message[-2000:],execution=e))
        return dict(status='pass',uninterrupted_state=exact_vector_hash(full['state']),continued_state=exact_vector_hash(continued['state']),dense_sparse_state_maximum=defect,rebuilt_dense_CT8_history_bitwise_equal=True,legacy_native_report_sha256=sha(a.ct8/'controls-v2/full/case.evolution.json'),negative=negative,rollback_failure=rejected['failure'],native_reports={str(f):sha(f/'case.evolution.json') for f in (first,end,failed,dense)})
    check('restart_rollback_backend_identity',histories)
    result['status']='pass' if len(result['checks'])==5 and all(q['status']=='pass' for q in result['checks'].values()) else 'fail';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
