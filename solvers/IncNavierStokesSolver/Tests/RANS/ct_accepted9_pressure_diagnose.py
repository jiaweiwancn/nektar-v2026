#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent bounded solves of sealed CT8 frozen blocks, never evolution."""
import argparse,json,os,copy,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from ct_accepted9_pressure_run import sealed
from ct_accepted8_native import read
from ct_accepted_feasibility import write_json,require
from ct_accepted_open import probe
from ct_accepted_pilot import replay
from ct_accepted5_traction import profiles
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity

def load_matrix(info):
    f=Path(info['file']);require(sha(f)==info['sha256'],'changed frozen matrix');a=np.fromfile(f,dtype=np.float64);require(a.size==info['rows']*info['columns'],'frozen matrix dimensions');return a.reshape(info['rows'],info['columns'])

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct8','runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');pre=sealed(a.ct8/'preflight-v6');old=sealed(a.ct8/'strip-T008-v1');c=pre['case'];n=pre['global_count'];weights=np.asarray(pre['weights']);scale=c['ustar']**2
    ini=read(Path(pre['initial']['0.001']['directory'])/'initial.json');o=ini['observation'];off=o['offsets'];regions=np.asarray(o['grid_regions'][:n]);initial=np.asarray(o['state'])
    require(len(initial)==2957,'not the declared dense strip')
    m=RequestManifest(a.output,dict(schema='CT9FrozenPressureDiagnosisV1',preflight=sha(a.ct8/'preflight-v6/results.json'),pair=sha(a.ct8/'strip-T008-v1/results.json'),source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),scope='all solves diagnostic only; no governing field/history modification; same-state quadrature contrast'))
    results=dict(schema='CT9FrozenPressureDiagnosisV1',members={},frozen_quadrature={},status='pass');reports={};frames={}
    for name,row in old['members'].items():
        file=Path(row['directory'])/'case.evolution.json';require(sha(file)==row['native_report_sha256'],'changed sealed native report');native=read(file);reports[name]=native
        attempt=next(v for v in native['pimple_attempts'] if v['scheduled_complete']);file=Path(attempt['flow_snapshot']);require(sha(file)==attempt['flow_snapshot_sha256'],'changed first flow snapshot');s=read(file)
        F=load_matrix(s['frozen_flow_matrix']);M=load_matrix(s['free_velocity_mass']);nv=s['velocity_size'];rows=np.asarray(s['rows'],int);npres=len(rows)-nv;alpha=s['alpha'];A=F[:nv,:nv];G=F[:nv,nv:];D=F[nv:,:nv];C=F[nv:,nv:];rm=np.asarray(s['frozen_flow_residual'][:nv]);rc=np.asarray(s['frozen_flow_residual'][nv:]);mass=lu_factor(M)
        B=np.zeros((n,npres));inverse={v:i for i,v in enumerate(rows[nv:])}
        for i,(q,region) in enumerate(zip(ini['test_basis'][:n],regions,strict=True)):
            for j,value in zip(q['pressure']['ids'],q['pressure']['phi'],strict=True):B[i,inverse[off[4+int(region)]+j]]=value
        def norms(v):
            vals=B@v;out={}
            for label,mask in [('whole',np.ones(n,bool)),('inner',regions==0),('outer',regions==1)]:
                w=weights[mask];x=vals[mask];mean=float(w@x/w.sum());out[label]=dict(weighted_RMS=float(np.sqrt(w@(x*x)/w.sum())),signed_mean=mean,spatial_variation_RMS=float(np.sqrt(w@((x-mean)**2)/w.sum())),scale=scale)
            out['maximum_scaled_value']=float(max(abs(vals),default=0.))/scale;return out
        invG=lu_solve(mass,G);W=D@invG;S=W/alpha-C;predictor=lu_solve(lu_factor(A),-rm);fac=lu_factor(S)
        native_s=np.asarray(s['pressure_operator']);action_error=float(np.max(abs(S-native_s)))/max(1.,float(np.max(abs(S))))
        require(action_error<2e-10,'consistent-mass pressure action reconstruction')
        momentum=lu_solve(fac,D@predictor);continuity=lu_solve(fac,rc);increment=momentum+continuity;native_increment=np.asarray(s['pressure_increment'])
        reconstruction=norms(increment-native_increment);require(reconstruction['maximum_scaled_value']<1e-8,'pressure reconstruction exceeds declared uncertainty')
        exact=lu_solve(lu_factor(F),-np.asarray(s['frozen_flow_residual']));residual=F@exact+np.asarray(s['frozen_flow_residual']);require(float(max(abs(residual)))<2e-10,'exact frozen diagnostic residual')
        # Only velocity storage enters this diagnostic asymptotic hypothesis.
        # Dependent scalar lifting may differ by roundoff and is recorded,
        # without treating it as a changed momentum history.
        compatible=bool(np.array_equal(np.asarray(s['input'])[:off[2]],initial[:off[2]]) and np.array_equal(np.asarray(s['history'])[:off[2]],initial[:off[2]]));eligible=not np.any(C) and compatible
        limits={};sweep=[]
        if eligible:
            wf=lu_factor(W);lm=-lu_solve(wf,D@lu_solve(mass,rm));lc=alpha*lu_solve(wf,rc)
            limits=dict(momentum=norms(lm),continuity=norms(lc),sum=norms(lm+lc),finite_alpha_remainder=norms(increment-lm-lc))
            L=A-alpha*M
            for factor in (1,2,4,8):
                al=alpha*factor;v=lu_solve(lu_factor(L+al*M),-rm);pi=lu_solve(lu_factor(W/al),D@v+rc)
                sweep.append(dict(alpha=al,dt=1/al,pressure=norms(pi),continuity_limit=norms(al*lu_solve(wf,rc)),remainder=norms(pi-lm-al*lu_solve(wf,rc))))
        tr={}
        for label,delta in [('native_split',np.concatenate([np.asarray(s['velocity_increment']),native_increment])),('exact_frozen',exact)]:
            candidate=np.asarray(s['input']).copy();candidate[rows]+=delta
            try:
                v=profiles(ini,candidate,c);tr[label]={k:v[k] for k in ('mean_wall_shear','normalized_wall_traction_error','normalized_interface_traction_RMS','undefined_interface')}
            except Exception as error:tr[label]=dict(status='diagnostic_stress_unavailable',reason=str(error))
        results['members'][name]=dict(snapshot_sha256=attempt['flow_snapshot_sha256'],actual_alpha=alpha,free_velocity_rows=nv,free_pressure_rows=npres,all_free_pressure_ids=rows[nv:].tolist(),unused_gauge_in_system=bool(len(initial)-1 in rows),maximum_C=float(np.max(abs(C))),weak_continuity_maximum=float(max(abs(rc))),weak_continuity_vector=rc.tolist(),history_lift_velocity_compatible=compatible,full_input_initial_maximum=float(np.max(abs(np.asarray(s['input'])-initial))),full_history_initial_maximum=float(np.max(abs(np.asarray(s['history'])-initial))),pressure_action_scaled_defect=action_error,reconstruction=reconstruction,native_pressure=norms(native_increment),momentum_contribution=norms(momentum),continuity_contribution=norms(continuity),exact_coupled_pressure=norms(exact[nv:]),split_minus_exact=norms(increment-exact[nv:]),exact_frozen_residual=float(max(abs(residual))),limit_applicable=eligible,limit=limits,alpha_sweep=sweep,tractions=tr,scope='extra diagnostic factorizations, never governing solves or an admission shortcut')
        frames[name]=(F,M,A,G,D,C,rm,rc,rows,alpha,B,norms,W)
    # Changing quadrature does not change the represented state, its physical
    # coordinates, constraints, pressure level, or history. Alpha0 loads are
    # distinguished from the first frozen stage's storage/splitting load.
    F,M,A,G,D,C,rm,rc,rows,alpha,B,measure,W=frames['0.001'];mass=lu_factor(M);wf=lu_factor(W);nv=M.shape[0]
    for label,state in [('initial',initial),('evolved_dt001',np.asarray(reports['0.001']['state']))]:
        loads=[];receipt=[]
        for order in [(32,8,12,12,64),(48,12,12,16,64),(64,16,12,20,64)]:
            dest=a.output/(label+'-q'+str(order[0]));dest.mkdir();root=ET.parse(Path(old['members']['0.001']['directory'])/'case.xml').getroot()
            for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),order,strict=True):set_parameter(root,key,value)
            write(root,dest/'case.xml');t=0. if label=='initial' else reports['0.001']['time'];req=replay(ini,state.tolist(),t,[],export_resolution=True)
            packet,e=probe(a.runtime,dest,req,'resolution');require(e['exit_code']==0,'same-state nonlinear quadrature replay failed')
            require(packet['observation']['state']==state.tolist(),'quadrature replay changed coefficients');v=np.asarray(packet['resolution']['native']['spatial_unmasked']);loads.append(v);receipt.append(dict(orders=order,packet_sha256=sha(dest/'resolution.json'),state_sha256=exact_vector_hash(state),execution=e))
        changes=[]
        for v in loads[1:]:
            dr=v[rows]-loads[0][rows];pm=-lu_solve(wf,D@lu_solve(mass,dr[:nv]));pc=alpha*lu_solve(wf,dr[nv:]);changes.append(dict(free_momentum_load_maximum=float(max(abs(dr[:nv]))),free_continuity_load_maximum=float(max(abs(dr[nv:]))),diagnostic_pressure_momentum=measure(pm),diagnostic_pressure_continuity=measure(pc),diagnostic_pressure_sum=measure(pm+pc)))
        results['frozen_quadrature'][label]=dict(controls=receipt,changes_from_native=changes,scope='nonlinear governing quadrature at fixed coefficients; pressure sensitivity uses the stored first-stage W, not a new trajectory or observation-rule change')
    write_json(a.output/'results.json',results);m.seal(a.output,['results.json']);print('CT9 diagnostic reconstruction complete',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
