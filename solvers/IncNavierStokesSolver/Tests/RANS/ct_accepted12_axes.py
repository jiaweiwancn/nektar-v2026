#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Cartesian-axis pressure images and retained-stage physical tractions."""
import argparse,gc,json,os,resource
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted12_pressure import Frozen,maximum
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import require,row_value
from ct_accepted_feasibility import write_json,parameters
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct10','ct11','pressure','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    primary=sealed(a.pressure);target=sealed(a.ct11/'target-v1');ini=n.read(a.ct10/'target-preflight-v1/target/readiness.json');c=target['case'];weights=n.inventory(c)[2];off=ini['observation']['offsets'];nm=off[2]
    manifest=RequestManifest(a.output,dict(schema='CT12AxisAndStageV1',source=sha(Path(__file__)),pressure=sha(a.pressure/'results.json'),target=sha(a.ct11/'target-v1/results.json'),scope='read-only retained target operators/rows; no native solve or trajectory'))
    result=dict(schema='CT12AxisAndStageV1',status='incomplete',members={},stage_tractions={});vectors={};common={};feedback_common={};cache={}
    coefficients=parameters('A95');alpha=coefficients['alpha'];norm=None
    ends={round(f['coordinate'][0],14):f['velocity_k'] for f in ini['observation']['faces'] if f['interface'] and f['region']==1}
    for point,side,row in zip(ini['observation']['grid'],ini['observation']['grid_sides'],ini['test_basis'],strict=True):
        if point[2]==c['H'] and side==1:ends[round(point[0],14)]=row['velocity_k']
    def diffusion(K):
        require(K>0,'actual endpoint K invalid');rough=c['roughness_input']*np.sqrt(K)/(70*c['nu'])
        B=10.7*(1-rough**max(.5-.4*rough**.7,1-.79*rough**(-.28)))
        return c['nu']+alpha*np.sqrt(K)*max(0.,-c['nu']*B/np.sqrt(K))
    reference_gradient=np.zeros((3,3));reference_gradient[0,2]=c['ustar']/(coefficients['kappa']*(c['H']-c['displacement']))
    refQ=np.asarray(oracle('A95',c['K0'],(c['K0']/coefficients['q'])**1.5/(coefficients['kappa']*(c['H']-c['displacement'])),reference_gradient)['Q'])
    def external(state):
        s=np.asarray(state);lanes={}
        for f in ini['observation']['faces']:
            if f['interface']:continue
            area=np.asarray(f['area']);kind='outlet' if area[0]>0 and f['coordinate'][0]==c['x'][-1] else ('top' if area[2]>0 and f['coordinate'][2]==c['Z'] else None)
            if kind is None:continue
            region=f['region'];q=f['velocity_k'];G=np.zeros((3,3));u=np.zeros(3)
            for block,axis in ((0,0),(1,2)):u[axis],G[axis]=row_value(q,s,off[block])
            k=row_value(q,s,off[2])[0];pi=row_value(f['pressure'],s,off[4+region])[0]
            if region:
                ep=row_value(f['epsilon'],s,off[3])[0];law=oracle('A95',k,ep,G)
                T=law['nut']*(G+G.T-2*np.eye(3)*np.trace(G)/3)-np.asarray(law['Q'])-pi*np.eye(3)
            else:
                K=row_value(ends[round(f['coordinate'][0],14)],s,off[2])[0];T=diffusion(K)*(G+G.T)-pi*np.eye(3)
            normal=area/np.linalg.norm(area);tr=T@normal
            datum=np.array([refQ[2,2]-refQ[0,0] if region and kind=='outlet' else 0.,0.,c['ustar']**2]) if kind=='outlet' else np.array([c['ustar']**2,0.,0.])
            delta=tr-datum;mask=np.array([True,False,kind=='outlet']);measure=f['weight']*np.linalg.norm(area)
            lane=lanes.setdefault(kind+'/'+str(region),dict(measure=0.,sum=np.zeros(3),traction_sum=np.zeros(3),square=0.,pressure=0.,witness=None,count=0))
            lane['measure']+=measure;lane['sum']+=measure*delta;lane['traction_sum']+=measure*tr;lane['square']+=measure*float(delta[mask]@delta[mask]);lane['pressure']+=measure*pi;lane['count']+=1
            square=float(delta[mask]@delta[mask])
            if lane['witness'] is None or square>lane['witness']['square']:lane['witness']=dict(coordinate=f['coordinate'],geometry=q['geometry'],side='outward',square=square,represented_traction=tr.tolist(),prescribed_datum=datum.tolist() if kind=='outlet' else {'x':float(datum[0]),'z':'not prescribed; w is essential'},compared_signed_difference=delta[[0,2]].tolist() if kind=='outlet' else [float(delta[0])])
        return {key:dict(compared_components=['x','z'] if key.startswith('outlet') else ['x'],signed_mean_difference=(v['sum']/v['measure'])[[0,2]].tolist() if key.startswith('outlet') else [float(v['sum'][0]/v['measure'])],represented_mean_traction=(v['traction_sum']/v['measure']).tolist(),normalized_RMS=float(np.sqrt(v['square']/v['measure']))/c['ustar']**2,
            absolute_Pi_mean=v['pressure']/v['measure'],count=v['count'],witness=v['witness'],scope='physical constitutive traction minus prescribed datum; top compares tangential component only, no normal-traction condition; separate from the weak natural load') for key,v in lanes.items()}
    def ground(state):
        s=np.asarray(state);rows=[]
        for f in ini['observation']['faces']:
            if f['interface'] or f['coordinate'][2]!=0:continue
            x=f['coordinate'][0];e=ends[round(x,14)];K=row_value(e,s,off[2])[0];require(K>0,'actual sampled K invalid')
            D=diffusion(K)
            ux=row_value(f['velocity_k'],s,off[0])[1];w=row_value(f['velocity_k'],s,off[1])[1];S=ux[2]+w[0]
            pi=row_value(f['pressure'],s,off[4])[0]
            rows.append(dict(x=x,weight=f['weight'],K=K,D=D,S=S,tau=D*S,normal_traction=2*D*w[2]-pi,geometry=f['velocity_k']['geometry']))
        w=np.asarray([r['weight'] for r in rows]);t=np.asarray([r['tau'] for r in rows]);mean=float(w@t/w.sum())
        normal=np.asarray([r['normal_traction'] for r in rows])
        return dict(mean=mean,error=abs(mean-c['ustar']**2)/c['ustar']**2,standard_deviation=float(np.sqrt(w@((t-mean)**2)/w.sum())),
            normal_traction_mean=float(w@normal/w.sum()),normal_traction_RMS=float(np.sqrt(w@(normal*normal)/w.sum())),
            minimum=min(rows,key=lambda r:r['tau']),maximum=max(rows,key=lambda r:r['tau']),count=len(rows)),rows
    def interface(mom):
        pairs={}
        for f in mom['interface_fluxes']:
            if f['analytic']:continue
            pairs.setdefault(round(f['row']['coordinate'][0],14),[]).append(f)
        total=np.zeros(3);normal=0.;measure=0.;worst=None
        for x,fs in pairs.items():
            require(len(fs)==2 and {f['region'] for f in fs}=={0,1},'one-sided interface coverage')
            f,g=fs;d=np.asarray(f['traction'])/f['area'][2]-np.asarray(g['traction'])/g['area'][2];w=f['weight'];measure+=w;total+=w*d*d;normal+=w*d[2]
            if worst is None or float(d@d)>worst['square']:worst=dict(coordinate=f['row']['coordinate'],square=float(d@d),jump=d.tolist(),owners=[f['row']['geometry'],g['row']['geometry']])
        return dict(full_jump_RMS=float(np.sqrt(total.sum()/measure))/c['ustar']**2,normal_jump_RMS=float(np.sqrt(total[2]/measure))/c['ustar']**2,signed_normal_jump=normal/measure/c['ustar']**2,witness=worst,count=len(pairs),scope='same +z traction on both sides; outward faces were opposite; no interface flux added to volume assembly')
    for dt,member in primary['members'].items():
        accum={k:np.zeros(off[6]-off[4]) for k in ('streamwise','normal','continuity')};states={};checks=[];feedback=np.zeros(off[6]-off[4]);feedback_states={}
        for row in member['steps']:
            s=n.read(row['snapshot']['path']);require(sha(Path(row['snapshot']['path']))==row['snapshot']['sha256'],'changed flow packet')
            fr=Frozen(s,ini,c,weights,len(weights))
            if norm is None:norm=fr.images
            load=np.zeros(nm);load[fr.rows[:fr.nv]]=fr.rm
            parts={}
            for name,axis in [('streamwise',0),('normal',1)]:
                v=load.copy();v[:off[1] if axis else 0]=0.;v[off[2] if axis else off[1]:]=0.
                pv=fr.image(v,'axis/'+name);parts[name]=pv;accum[name]+=pv
            pc=fr.solve(fr.rc,'axis/continuity');parts['continuity']=pc;accum['continuity']+=pc
            defect=fr.measure(sum(parts.values(),np.zeros_like(pc))-np.asarray(s['pressure_increment']))
            require(defect['maximum_scaled_value']<1e-8,'Cartesian-axis signed reconstruction')
            native_pi=np.asarray(s['pressure_increment']);native_b=np.asarray(s['subtracted_mass_velocity_correction'])
            mass_defect=maximum(fr.alpha*fr.mul('M',native_b)-fr.mul('G',native_pi))
            action_defect=maximum(fr.mul('G',native_pi)-np.asarray(s['full_pressure_action'])[fr.rows[:fr.nv]])
            require(max(mass_defect,action_defect)<2e-10,'native complete pressure/mass action')
            old_pi=np.asarray(s['input'])[fr.rows[fr.nv:]]-np.asarray(ini['observation']['state'])[fr.rows[fr.nv:]]
            pi_rhs=-fr.mul('D',fr.direct('A',fr.mul('G',old_pi)))
            pi_feedback=fr.solve(pi_rhs,'stored-pressure-feedback');feedback+=pi_feedback
            repeat=None
            if round(row['time'],12)==.002:
                repeat=fr.measure(fr.solve(pi_rhs,'stored-pressure-feedback-tighter',2e-14)-pi_feedback)
                require(repeat['maximum_scaled_value']<1e-8,'stored-pressure feedback repeat')
            checks.append(dict(time=row['time'],stored_pressure_feedback=fr.measure(pi_feedback),feedback_tighter_repeat=repeat,axes={k:fr.measure(v) for k,v in parts.items()},signed_sum_error=defect,
                independent_native_mass_action_defect=mass_defect,complete_native_G_action_defect=action_defect,solves=fr.receipts))
            if round(row['time'],12) in (.001,.002):
                states[str(round(row['time'],12))]={k:v.copy() for k,v in accum.items()};feedback_states[str(round(row['time'],12))]=feedback.copy()
            result['members'][dt]=checks
            write_json(a.output/'progress.json',result)
            # S/W LinearOperator closures form Python reference cycles. Release
            # each complete factor set explicitly before constructing another.
            # Keep only the fixed physical observation object for comparisons.
            del fr;gc.collect()
        result['members'][dt]=checks;common[dt]=states;feedback_common[dt]=feedback_states
    for t in ('0.001','0.002'):
        delta={k:common['0.001'][t][k]-common['0.0005'][t][k] for k in accum};total=sum(delta.values(),np.zeros_like(pc));den=norm.inner(total,total)
        result.setdefault('common',{})[t]={k:dict(signed_alignment=norm.inner(v,total)/den,image=norm.measure(v)) for k,v in delta.items()};vectors[t]={k:v.tolist() for k,v in delta.items()}
        fp=feedback_common['0.001'][t]-feedback_common['0.0005'][t]
        result.setdefault('stored_pressure_feedback',{})[t]=dict(signed_alignment=norm.inner(fp,total)/den,image=norm.measure(fp),scope='G times actual stored pressure departure from unchanged initial Pi; contained in axis/component loads, not an extra independent contribution; no regauging')
    # Last fine step: actual input,corrected flow,pre-k and accepted fields.
    r=n.read(Path(target['members']['0.0005']['directory'])/'case.evolution.json');attempt=r['pimple_attempts'][-1];fs=n.read(attempt['flow_snapshot']);ks=n.read(attempt['k_snapshot'])
    result['initial_tractions']=dict(ground=ground(ini['observation']['state'])[0],external=external(ini['observation']['state']),state_sha256=exact_vector_hash(ini['observation']['state']))
    states={'input':fs['input'],'post_flow':ks['post_flow'],'pre_k':ks['pre_k'],'accepted':r['state']};previous=None;decomposition={}
    for name,state in states.items():
        stats,rows=ground(state);result['stage_tractions'][name]=dict(ground=stats,external=external(state),state_sha256=exact_vector_hash(state))
        if previous is not None:
            _,oldrows=previous;w=np.asarray([r['weight'] for r in rows]);terms=np.zeros(3);error=0.
            for old,new in zip(oldrows,rows,strict=True):
                ds=new['S']-old['S'];dd=new['D']-old['D'];v=np.array([old['D']*ds,old['S']*dd,dd*ds]);terms+=new['weight']*v;error=max(error,abs(new['tau']-old['tau']-v.sum()))
            require(error/c['ustar']**2<2e-10,'stage traction product decomposition')
            decomposition[previous[0]+'->'+name]=dict(gradient=terms[0]/w.sum(),endpoint_viscosity=terms[1]/w.sum(),cross=terms[2]/w.sum(),maximum_reconstruction_scaled=error/c['ustar']**2)
        previous=(name,rows)
    for label,record in primary['stage_loads'].items():
        path=Path(record['metadata']['path']);require(sha(path)==record['metadata']['sha256'],'changed native interface diagnostic');packet=n.read(path)
        result['stage_tractions'].setdefault(label,{})['interface']=interface(packet['ct10_momentum'])
    expected=target['members']['0.0005']['accepted_metrics'][-1]['traction']['normalized_wall_traction_error']
    require(abs(result['stage_tractions']['accepted']['ground']['error']-expected)<2e-10,'raw evolved ground traction disagreement')
    vm={k:int(v.split()[0]) for k,v in (line.split(':',1) for line in Path('/proc/self/status').read_text().splitlines() if line.startswith(('VmPeak:','VmHWM:')))}
    result.update(status='pass',ground_product_decomposition=decomposition,continuation_admitted=False,diagnostic_memory=dict(maximum_self_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**vm),factor_lifecycle='one Frozen factor set at a time; explicit del/gc after each step; fixed observation object retained')
    write_json(a.output/'axis-vectors.json',vectors);write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json','axis-vectors.json']);return 0
if __name__=='__main__':raise SystemExit(main())
