#!/usr/bin/env python3
"""Source-pinned Linux verification of the prescribed upstream transition.

Small manufactured runs verify the new coupling, not long-time CLS stability.
Raw meshes, states, logs and full-history checkpoints stay in the run directory.
"""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
THREADS={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from validate_cls_epawf_m4 import Runner as NativeRunner, state_difference
from validate_cls_epawf_m5 import spaces
from cls_epawf_m6_cases import transition_session,arc_offset
from cls_epawf_m3_cases import write,set_info,set_parameter

checks=0
errors={}
def require(value,message):
    global checks
    checks+=1
    if not value:raise RuntimeError(message)
def near(a,b,name,scale=1.,tolerance=1e-10):
    a,b=np.asarray(a),np.asarray(b)
    require(a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)),'invalid '+name)
    error=float(np.max(np.abs(a-b),initial=0.))/max(scale,float(np.max(np.abs(b),initial=0.)))
    errors[name]=max(errors.get(name,0.),error)
    require(error<=tolerance,f'{name}: {error} > {tolerance}')
    return error
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class Runner(NativeRunner):
    def probe(self,name,**kwargs):
        directory=self.output/name;directory.mkdir()
        xml=transition_session(self.template,directory/'case.xml',probe=True,**kwargs)
        self.command(name,'RANSCurvedWallProbe',[xml,directory/'probe.json'],kwargs.pop('ranks',1) if 'ranks' in kwargs else 1)
        return json.loads((directory/'probe.json').read_text())

def manufactured_error(raw,stationary=False):
    d=raw['dimension'];f=spaces(raw);off=raw['offsets'];time=raw['time']
    q,kappa=raw['parameters'][6:8];K0=q*.16;K=K0 if stationary else K0*(1+.02*math.sin(time))
    A=K**1.5/2.55;delta=(K/q)**1.5/kappa-A
    offset=arc_offset() if int(raw['map_parameters'][0])==4 else 0.
    square=np.zeros(d+4);volume=np.zeros(d+4)
    for b,field in enumerate(f):
        for e in field.elements:
            points,w=field.quadrature(e,24);ids,p,g=field.rows(e,points)
            actual=p@np.asarray(raw['state'])[off[b]+ids]
            if b<d or b>=d+2:exact=np.zeros_like(actual)
            elif b==d:exact=np.full_like(actual,K)
            else:
                s=np.clip((points[:,0]-offset+3)/.75,0,1)
                chi=1-10*s**3+15*s**4-6*s**5
                exact=(A+delta*chi)/points[:,2]
            square[b]+=np.dot(w,(actual-exact)**2);volume[b]+=np.sum(w)
    rms=np.sqrt(square/volume)
    values={'velocity':float(np.linalg.norm(rms[:d])),'k':float(rms[d]),'epsilon':float(rms[d+1]),
            'pressure':float(np.linalg.norm(rms[d+2:]))}
    require(raw['eta']==0,'transition/terrain anchor is nonzero')
    require(raw['upstream_transition'],'missing actual production transition')
    require(all(np.isfinite(list(values.values()))),'nonfinite manufactured error')
    trace_value=trace_gradient=measure=omission=0.
    z,w=np.polynomial.legendre.leggauss(32)
    for e in f[d+1].elements:
        if abs(e['lower'][2]-2)>1e-11:continue
        half=(e['upper'][0]-e['lower'][0])/2
        for xx,ww in zip(e['lower'][0]+half*(z+1),half*w):
            point=[xx,.37 if d==3 else 0.,2.]
            value,grad=f[d+1].evaluate(raw['state'],off[d+1],point,1)
            s=float(np.clip((xx-offset+3)/.75,0,1));chi=1-10*s**3+15*s**4-6*s**5
            derivative=-30*s*s*(1-s)**2/.75
            trace_value+=ww*(value-(A+delta*chi)/2)**2
            # The nonzero analytic tangential derivative occurs only on the
            # exactly flat fetch, where the physical and chart x axes agree.
            if xx-offset<=-2.:
                trace_gradient+=ww*(grad[0]-delta*derivative/2)**2
                omission+=ww*(delta*derivative/2)**2
            measure+=ww
    values.update(trace_value_rms=math.sqrt(trace_value/measure),trace_gradient_rms=math.sqrt(trace_gradient/measure),
                  omitted_grad_chi_rms=math.sqrt(omission/measure))
    require(values['omitted_grad_chi_rms']>1e-5,'manufactured trace has no detectable grad(chi) contribution')
    for sample in raw['transition_samples']:
        ell=sample['input'][0];gates=sample['gates']
        require(0<=gates[0]<=1,'analytic chi left [0,1]')
        if ell>=-2.25:near(gates[0],0,'no_downstream_reactivation')
        require(gates[2]==1 and gates[3]==1,'MMS left its declared state/admissibility unit plateau')
    return values

def stage_audit(raw):
    require(raw['eta']==0 and raw['upstream_transition'],'invalid terrain/transition stage policy')
    require(bool(raw['stage_audits']),'no native stage was solved')
    maximum=max(max(row[6:]) for row in raw['stage_audits'])
    require(maximum<1e-9,'unconverged transition coupled stage')
    near(np.sum(raw['regional_stage_loads'],axis=0),raw['final_stage_loads'],'regional_transition_stage_sum',tolerance=1e-10)
    require(abs(raw['state'][-1])<1e-9,'gauge multiplier would create a mass source')
    return {'maximum_scaled_residual':maximum,'unknowns':len(raw['state']),
            'boundary_surface_rms_correction':raw['boundary_flux_correction'],
            'energy_identity_error':abs(raw['energy_budget'][7])}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--mode',choices=['all','probes','smoke'],default='all')
    args=p.parse_args();require(platform.system()=='Linux','native numerical qualification belongs on Linux')
    run=Runner(args.build,args.output);sections={};source=Path(__file__).resolve().parents[4]
    summary={'schema':1,'milestone':'M6','started_utc':datetime.now(timezone.utc).isoformat(),
             'source_root':str(source),'thread_environment':THREADS,'status':'running','sections':sections}
    def save():
        summary.update(checks=checks,errors=errors,commands=run.commands)
        (run.output/'validation.json').write_text(json.dumps(summary,indent=2)+'\n')
    try:
        if args.mode in ('all','probes'):
            from cls_epawf_m6_operator_oracle import verify
            sections['operators']=[]
            for family,dimension in [('flat',2),('feet',2),('flat',3)]:
                raw=run.probe(f'operator-{family}-d{dimension}',family=family,dimension=dimension,quadrature=8)
                sections['operators'].append(verify(raw,near,require,directional=dimension==2));save()
        if args.mode in ('all','smoke'):
            sections['manufactured']=[];cases={}
            for family,degrees in [('flat',(3,4,5)),('feet',(3,4))]:
                records=[]
                for degree in degrees:
                    name=f'mms-{family}-p{degree}'
                    raw=run.evolution(name,generator=transition_session,family=family,degree=degree,steps=2,dt=.0025,quadrature=12)
                    value={'family':family,'degree':degree,'refine':1,'error':manufactured_error(raw),'stage':stage_audit(raw)}
                    records.append(value);sections['manufactured'].append(value);cases[(family,degree)]=raw;save()
                require(records[-1]['error']['epsilon']<records[0]['error']['epsilon']*.8,'outer epsilon p refinement failed: '+family)
                require(records[-1]['error']['trace_gradient_rms']<records[0]['error']['trace_gradient_rms']*.8,'trace gradient p refinement failed: '+family)
            refined=run.evolution('mms-flat-h2',generator=transition_session,family='flat',degree=3,refine=2,steps=2,dt=.0025)
            er=manufactured_error(refined);stage_audit(refined)
            require(er['epsilon']<sections['manufactured'][0]['error']['epsilon']*.8,'fixed-physical-position h refinement failed')
            sections['h_refinement']={'refine':2,'error':er};save()
            # Identical prescribed MMS and mesh: compare actual full/Schur
            # stages, collective ownership, and every native history slot.
            base=cases[('feet',3)]
            condensed=run.evolution('condensed-feet',generator=transition_session,family='feet',degree=3,steps=2,dt=.0025,condensed=True)
            stage_audit(condensed);state_difference(base,condensed,'M6_full_condensed')
            sections['parallel_restart']=[]
            for ranks in (2,4):
                other=run.evolution('mpi-feet-'+str(ranks),generator=transition_session,family='feet',degree=3,steps=2,dt=.0025,ranks=ranks,ownership_shift=1)
                stage_audit(other);state_difference(base,other,'M6_MPI_'+str(ranks))
                near(other['native_history'],base['native_history'],'M6_MPI_history_'+str(ranks),tolerance=1e-8)
                sections['parallel_restart'].append({'ranks':ranks,'stage':stage_audit(other)});save()
            checkpoint=run.output/'restart-first'/'state'
            first=run.evolution('restart-first',generator=transition_session,family='feet',degree=3,steps=1,dt=.0025,checkpoint=str(checkpoint))
            resumed=run.evolution('restart-last',generator=transition_session,family='feet',degree=3,steps=1,dt=.0025,
                                  restart=first['checkpoint'],ranks=4,ownership_shift=2)
            stage_audit(resumed);state_difference(base,resumed,'M6_full_history_restart')
            near(resumed['native_history'],base['native_history'],'M6_complete_restart_history',tolerance=1e-8)
            near(resumed['native_history_times'],base['native_history_times'],'M6_complete_restart_history_times',tolerance=1e-12)
            sections['parallel_restart'].append({'restart_ranks':4,'all_history_slots':len(resumed['native_history'])});save()
            # Actual native startup IMEX1 on the same newly coupled operator.
            one=run.evolution('imex1-transition',generator=transition_session,family='flat',degree=3,steps=2,order=1,dt=.0025)
            sections['imex1']={'error':manufactured_error(one),'stage':stage_audit(one)};save()
        summary.update(status='passed',finished_utc=datetime.now(timezone.utc).isoformat());save()
        print(json.dumps({'status':summary['status'],'checks':checks,'commands':len(run.commands)}),flush=True)
    except Exception as error:
        summary.update(status='failed',error=str(error),finished_utc=datetime.now(timezone.utc).isoformat());save();raise

if __name__=='__main__':main()
