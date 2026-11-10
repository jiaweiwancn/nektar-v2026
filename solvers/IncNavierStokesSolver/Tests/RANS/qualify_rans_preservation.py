#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Explicit bounded preservation jobs. Reports survive any failed gate.

Run spatial first; temporal never silently inherits a convergence conclusion.
Runtime binaries/libraries are frozen separately and identified by source-build.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from rans_qualification import THREADS,dependencies,write_json,coverage,compare_samples,exact_vector_hash,audited_extrema
dependencies('numpy','scipy')
import numpy as np
from rans_preservation_oracle import Reference,physical_grid,weak_balance,representation
from run_a93_a95_e94e95_matrix import case_xml,run_case
from run_cls_epawf_picard import sha
from rans_request_manifest import RequestManifest,runtime_identity
from run_cls_epawf_pimple import identity
import tempfile
from contextlib import nullcontext

Q=((32,8,8,8,8),(48,12,12,12,12),(64,16,16,16,16))
MODELS=('E94E95','A95','A93');POLICIES=('Natural','Prescribed');WALLS=('EP','Accepted')

def load(path):return json.loads(Path(path).read_text())
def nested(knots,level):
    for _ in range(level):knots=sorted(set(knots+[(a+b)/2 for a,b in zip(knots[:-1],knots[1:])]))
    return knots

def compact_packet(root):
    return {k:v for k,v in root.items() if k!='observation'} | {'observation':{k:v for k,v in root['observation'].items() if k not in ('volume','faces','pressure_mass')}}

def probe(executable,folder,request):
    write_json(folder/'request.json',request)
    with (folder/'probe.log').open('w') as f:
        child=subprocess.run([str(executable),'case.xml','request.json','packet.json'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=600,env=dict(os.environ,**THREADS))
    path=folder/'packet.json'
    if not path.exists():raise RuntimeError('probe produced no packet; exit '+str(child.returncode))
    full=load(path)
    if child.returncode or full.get('status')!='pass':raise RuntimeError('native read-only audit failed: '+json.dumps({k:v for k,v in full.items() if k!='observation'}))
    out=compact_packet(full);write_json(folder/'state-packet.json',out)
    return full

def archive_basis(folder):
    path=folder/'packet.json'
    with path.open('rb') as src,gzip.open(folder/'packet.json.gz','wb',compresslevel=3) as dst:
        import shutil;shutil.copyfileobj(src,dst)
    path.unlink()

def replay_request(base,grid,export_basis=False,state=None,time=0.):
    o=base['observation']
    return dict(operation='replay',identity=o['identity'],descriptor=o['resolved_descriptor'],reference=base['reference'],
                fixed_boundary=base['fixed_boundary'],state=o['state'] if state is None else state,time=time,grid=grid,export_basis=export_basis)

def residual_metrics(audit):
    # All are native unanchored alpha=0 receipts. Preserve original grouping.
    return {k:v for k,v in audit.items() if any(s in k for s in ('norm','dual','interface','trace','gauge','covariance','radicand','cutoff','divergence','traction'))}

def series_gate(values,uncertainty=1e-10):
    if len(values)!=3:return dict(status='fail',reason='three levels required',values=values)
    ratios=[values[i+1]/values[i] if values[i]>1e-30 else None for i in (0,1)]
    if max(values)<=5*uncertainty:return dict(status='unresolved_floor',values=values,uncertainty=uncertainty,ratios=ratios)
    ok=all(r is not None and r<=.85 for r in ratios) and all(abs(values[i]-values[i+1])>5*uncertainty for i in (0,1))
    return dict(status='pass' if ok else 'fail',values=values,uncertainty=uncertainty,ratios=ratios)

class IncompleteCampaign(Exception):pass

class Campaign:
    def __init__(self,args):
        self.r=args.root.resolve();self.phase=args.phase;self.models=tuple(args.models);self.policies=tuple(args.policies);self.walls=tuple(args.walls);self.epsilon_perturbation=args.epsilon_perturbation;self.temporal_quadrature=Q[args.quadrature_index];self.output=self.r/('preservation-'+args.phase+('-'+args.label if args.label else ''));self.horizon=args.horizon;self.boundary_policy=args.boundary_policy
        runtime=args.runtime_directory.resolve() if args.runtime_directory else self.r/'runtime'
        self.solver=runtime/'IncNavierStokesSolver';self.descriptor=runtime/'RANSNewModelsProbe';self.sampler=runtime/'RANSPreservationProbe'
        self.manifest=load(runtime/'source-build.json');self.result=dict(schema=1,phase=self.phase,checks=[],status='running',
            runtime_manifest_sha256=sha(runtime/'source-build.json'),runtime_source=self.manifest['revision'],ranks=1,threads=THREADS,
            source=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip())
        if sha(self.solver)!=self.manifest['solver']['sha256']:raise RuntimeError('frozen solver identity changed')
        os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(self.r/'tmp'))
        request=dict(verification_schema='PreservationCampaignV2',phase=self.phase,models=self.models,policies=self.policies,walls=self.walls,
            quadrature=list(self.temporal_quadrature),horizon=self.horizon,timesteps=[1e-4,5e-5,2.5e-5],integration='native IMEX1 PimpleSingleSweep consistent histories',
            epsilon_perturbation=self.epsilon_perturbation,boundary_policy=args.boundary_policy,boundary_identity=args.boundary_identity,
            runtime={name:runtime_identity(path) for name,path in [('solver',self.solver),('sampler',self.sampler),('descriptor',self.descriptor)]},initial_data=[])
        if self.output.exists() and not args.resume:raise FileExistsError('destination exists: explicit identical --resume required')
        if self.phase=='temporal':
            evaluation=load(self.r/'preservation-spatial/evaluation-grid.json')
            request['evaluation']=dict(grid=evaluation['grid'],weights=evaluation['weights'],sides=[0]*len(evaluation['grid']),export_basis=False)
        else:
            refs=[Reference(m,p,w) for m in self.models for p in self.policies for w in self.walls]
            panels=sorted(set(v for ref in refs for v in nested(ref.knots(),2)));g,w=physical_grid(panels,10,8)
            request['evaluation']=dict(grid=g,weights=w,sides=[0]*len(g),export_basis=True)
            request['spatial_series']=dict(reference_degree=3,h_levels=[0,1,2],p_degrees=[3,4,5],quadratures=Q,representation='native constrained neutral RadialProject')
        with nullcontext(tempfile.mkdtemp(prefix='request-preflight-',dir=self.r/'tmp')) as tmp:
            for m in self.models:
                for policy in self.policies:
                    for wall in self.walls:
                        ref=Reference(m,policy,wall);folder=Path(tmp)/('-'.join((m,policy,wall)))
                        self.config(folder,ref,ref.knots(),q=self.temporal_quadrature)
                        write_json(folder/'immutable-initialization-intent.json',dict(schema='ReadOnlyInitializationIntentV1',model=m,policy=policy,wall=wall,mesh=ref.knots(),degree=3,quadrature=list(self.temporal_quadrature),boundary_policy=self.boundary_policy,runtime=request['runtime'],operation='read-only neutral binding',time=0,steps=0))
                        full=probe(self.sampler,folder,dict(operation='neutral',export_basis=False,grid=[]))
                        request['initial_data'].append(dict(model=m,policy=policy,wall=wall,mesh=ref.knots(),degree=3,pressure_degree=1,
                            resolved=full['observation']['resolved_descriptor'],raw_sha256=exact_vector_hash(full['observation']['state']),
                            reference_sha256=exact_vector_hash(full['reference']),boundary=full['fixed_boundary'],representation='native constrained neutral projection'))
        self.admission=RequestManifest(self.output,request,args.resume)
        if args.resume and (self.output/'results.json').exists():self.result=load(self.output/'results.json')
        self.resume=args.resume;self.new_check_limit=args.new_check_limit;self.new_checks=0
        self.save()
    def save(self):write_json(self.output/'results.json',self.result)
    def check(self,name,fn):
        if self.resume and any(x['name']==name and x['status']=='pass' for x in self.result['checks']):
            folder=self.output/name.replace('-reference','-h0')
            if self.admission.completed(folder):return
        if self.new_check_limit is not None and self.new_checks>=self.new_check_limit:raise IncompleteCampaign('explicit new-check budget reached; identical resume required')
        self.new_checks+=1
        self.result['checks']=[x for x in self.result['checks'] if x['name']!=name]
        start=time.monotonic()
        try:
            value=fn();status=value.get('status','pass')
            self.result['checks'].append(dict(name=name,status=status,seconds=time.monotonic()-start,result=value))
        except Exception as error:
            self.result['checks'].append(dict(name=name,status='fail',seconds=time.monotonic()-start,failure=repr(error)))
        self.save();print(json.dumps({k:self.result['checks'][-1][k] for k in ('name','status','seconds')}),flush=True)
    def config(self,folder,ref,knots,p=3,q=Q[0],steps=0,dt=.0001,perturb=0.):
        return case_xml(folder,ref.model,ref.policy,ref.wall,steps=steps,dt=dt,descriptor_probe=self.descriptor,radial_knots=knots,degree=p,
                        quadrature=q,boundary_policy=self.boundary_policy,checkpoint_interval=max(steps,1),cache='ExactStepV1',perturb=perturb)
    def frozen(self,label,ref,knots,p=3,q=Q[0],base=None,basis=False,grid=None):
        folder=self.output/label
        config=dict(model=ref.model,policy=ref.policy,wall=ref.wall,knots=knots,degree=p,quadrature=list(q),runtime_manifest_sha256=self.result['runtime_manifest_sha256'])
        if self.resume and self.admission.completed(folder):
            metrics=load(folder/'metrics.json')
            if metrics.get('configuration')!=config:raise RuntimeError('existing frozen result identity mismatch')
            return load(folder/'state-packet.json'),metrics
        self.config(folder,ref,knots,p,q)
        request=dict(operation='neutral',export_basis=basis,grid=grid or []) if base is None else replay_request(base,grid or [],basis)
        full=probe(self.sampler,folder,request);packet=full['observation']
        metrics=dict(status='pass',configuration=config,execution_status='completed',physical_admissibility=full['physical_admissibility'],
            state_sha256=exact_vector_hash(packet['state']),reference_sha256=exact_vector_hash(full['reference']),
            actual_quadrature=packet['actual_quadrature'],residuals=residual_metrics(packet['audit_alpha_zero']),
            native_alpha_zero=packet['audit_alpha_zero'],source_XML_sha256=sha(folder/'case.xml'))
        if basis:metrics['analytical_weak_balance']=weak_balance(ref,packet)
        if grid:metrics['representation']=representation(ref,packet,self.grid_weights)
        write_json(folder/'metrics.json',metrics)
        self.admission.seal(folder,['case.xml','metrics.json','state-packet.json'])
        if basis:archive_basis(folder)
        else:(folder/'packet.json').unlink()
        return full,metrics
    def spatial(self):
        # A common evaluation grid spanning every refinement and reference breakpoint.
        refs=[Reference(m,p,w) for m in self.models for p in self.policies for w in self.walls]
        allknots=sorted(set(v for r in refs for v in nested(r.knots(),2)))
        grid,self.grid_weights=physical_grid(allknots,10,8)
        write_json(self.output/'evaluation-grid.json',dict(grid=grid,weights=self.grid_weights,order=10,xorder=8,panels=allknots))
        expected=[]
        for ref in refs:
            tag='-'.join((ref.model,ref.policy,ref.wall));expected.append(tag+'-reference')
            def initial(ref=ref,tag=tag):
                base,metrics=self.frozen(tag+'-h0',ref,ref.knots(),basis=True,grid=grid)
                write_json(self.output/(tag+'-base.json'),compact_packet(base))
                status=metrics['analytical_weak_balance']['status']
                return dict(status=status,model=ref.model,kappa_policy=ref.policy,wall=ref.wall,**{k:v for k,v in metrics.items() if k!='status'})
            self.check(tag+'-reference',initial)
            basepath=self.output/(tag+'-base.json')
            if not basepath.exists():continue
            base=load(basepath)
            if ref.model=='E94E95' and ref.policy=='Natural':
                for j in (1,2):
                    name=tag+'-q'+str(j);expected.append(name)
                    def quadrature(j=j,ref=ref,tag=tag,base=base):
                        _,metrics=self.frozen(tag+'-q'+str(j),ref,ref.knots(),q=Q[j],base=base,basis=j==2,grid=grid)
                        if metrics['state_sha256']!=exact_vector_hash(base['observation']['state']):raise AssertionError('quadrature changed frozen coefficients')
                        return metrics
                    self.check(name,quadrature)
                for level in (1,2):
                    name=tag+'-h'+str(level);expected.append(name)
                    self.check(name,lambda level=level,ref=ref,tag=tag:self.frozen(tag+'-h'+str(level),ref,nested(ref.knots(),level),q=Q[2],grid=grid)[1])
                for degree in (4,5):
                    name=tag+'-p'+str(degree);expected.append(name)
                    self.check(name,lambda degree=degree,ref=ref,tag=tag:self.frozen(tag+'-p'+str(degree),ref,ref.knots(),p=degree,q=Q[2],grid=grid)[1])
            else:
                name=tag+'-h1';expected.append(name)
                self.check(name,lambda ref=ref,tag=tag:self.frozen(tag+'-h1',ref,nested(ref.knots(),1),q=Q[2],grid=grid)[1])
        self.result['expected_check_names']=expected
        self.result['coverage']=coverage(self.result['checks'],set(expected),key=lambda x:x['name'])
        self.spatial_tables()
    def spatial_tables(self):
        tables={}
        for wall in self.walls:
            if 'E94E95' not in self.models or 'Natural' not in self.policies:continue
            tag='E94E95-Natural-'+wall
            paths=[self.output/(tag+s)/'metrics.json' for s in ('-h0','-h1','-h2','-p4','-p5','-q1','-q2')]
            if not all(p.exists() for p in paths):tables[tag]=dict(status='fail',reason='missing study levels');continue
            m=[load(p) for p in paths];byname={s:x for s,x in zip(('h0','h1','h2','p4','p5','q1','q2'),m,strict=True)}
            data={}
            for field in ('u','v','k','epsilon_outer','Pi'):
                for typ in ('values','gradient'):
                    def value(x):return x['representation']['fields'][field]['weighted_RMS'] if typ=='values' else x['representation']['gradient_weighted_RMS'][field]
                    qerr=abs(value(byname['q2'])-value(byname['q1']))
                    data[field+'_'+typ]=dict(h=series_gate([value(byname[k]) for k in ('h0','h1','h2')],max(1e-10,qerr)),
                        p=series_gate([value(byname[k]) for k in ('h0','p4','p5')],max(1e-10,qerr)),frozen_quadrature_change=qerr)
            # The Q2 replay supplies h0 at exactly the same quadrature as h/p.
            byname['h0_resolved']=byname['q2']
            spatial={}
            for key,indices in [('coupled_norms',list(range(8))),('native_k_inner_outer_shared_interface_dual',list(range(3)))]:
                labels=['momentum_x','momentum_z','global_k','outer_epsilon','continuity_inner','continuity_outer','gauge','dependent_trace'] if key=='coupled_norms' else ['inner_k','outer_k','shared_interface_k']
                for i,label in enumerate(labels):
                    def v(k):return byname[k]['native_alpha_zero'][key][i]
                    uncertainty=max(1e-10,abs(v('q2')-v('q1')))
                    h=series_gate([v(k) for k in ('h0_resolved','h1','h2')],uncertainty)
                    pp=series_gate([v(k) for k in ('h0_resolved','p4','p5')],uncertainty)
                    differences=[abs(v('h0_resolved')-v('h1')),abs(v('h1')-v('h2')),abs(v('h0_resolved')-v('p4')),abs(v('p4')-v('p5'))]
                    quadrature_ok=all(abs(v('q2')-v('q1'))<=max(1e-10,.01*d) for d in differences)
                    spatial[label]=dict(h=h,p=pp,quadrature_status='pass' if quadrature_ok else 'unresolved',quadrature_uncertainty=abs(v('q2')-v('q1')))
            # Preserve the full native arrays rather than guessing their labels.
            qa=byname['q1']['native_alpha_zero'];qb=byname['q2']['native_alpha_zero']
            qchanges={k:max((abs(a-b) for a,b in zip(v,qb[k],strict=True)),default=0.) for k,v in qa.items() if k in qb and len(v)==len(qb[k]) and ('norm' in k or 'dual' in k)}
            tables[tag]=dict(representation=data,native_spatial=spatial,quadrature_native_last_changes=qchanges,levels=byname,
                status='fail' if any(g['status']=='fail' for d in data.values() for g in (d['h'],d['p'])) or any(g['status']=='fail' for d in spatial.values() for g in (d['h'],d['p'])) or any(d['quadrature_status']!='pass' for d in spatial.values()) or any(x.get('analytical_weak_balance',{}).get('status','pass')!='pass' for x in m) else 'pass')
        self.result['spatial_tables']=tables;self.save()
    def temporal(self):
        # The first explicit temporal study uses the same fixed aligned p3 spaces.
        # Select a measured adequate rule, never only a wall-normal increase.
        spatial=load(self.r/'preservation-spatial/evaluation-grid.json');grid=spatial['grid'];weights=spatial['weights']
        self.grid_weights=weights;expected=[];tables=[]
        for model in self.models:
            for policy in self.policies:
                for wall in self.walls:
                    ref=Reference(model,policy,wall);tag='-'.join((model,policy,wall));runs=[]
                    for j,dt in enumerate((1e-4,5e-5,2.5e-5)):
                        name=tag+'-dt'+str(j);expected.append(name)
                        def trajectory(j=j,dt=dt,ref=ref,tag=tag):
                            folder=self.output/(tag+'-dt'+str(j));steps=round(self.horizon/dt)
                            env_value=str(folder/'timing.jsonl');os.environ['NEKTAR_RANS_TIMING_FILE']=env_value
                            if self.resume and self.admission.completed(folder):return load(folder/'temporal-metrics.json')
                            options=dict(steps=steps,dt=dt,descriptor_probe=self.descriptor,
                                radial_knots=ref.knots(),degree=3,quadrature=self.temporal_quadrature,cache='ExactStepV1',checkpoint_interval=steps,boundary_policy=self.boundary_policy,epsilon_perturbation=self.epsilon_perturbation)
                            row=load(folder/'result.json') if self.resume and self.admission.completed(folder,'trajectory-artifacts.json') else run_case(self.solver,folder,ref.model,ref.policy,ref.wall,**options)
                            if row['status']=='pass':self.admission.seal(folder,['case.xml','result.json','case.evolution.json'],'trajectory-artifacts.json')
                            if row['status']!='pass':return dict(status='fail',execution=row)
                            report=load(folder/'case.evolution.json')
                            # Reconstruct exactly the initializer on these same spaces.
                            initial=probe(self.sampler,folder,dict(operation='neutral',export_basis=False,grid=[]))
                            if exact_vector_hash(initial['observation']['state'])!=exact_vector_hash(report['fixed_reference']):raise AssertionError('probe/native neutral reference coefficients differ')
                            if not self.epsilon_perturbation and exact_vector_hash(initial['observation']['state'])!=exact_vector_hash(report['accepted_states'][0]):raise AssertionError('probe/native initial coefficients differ')
                            request=replay_request(initial,grid,state=report['state'],time=report['time'])
                            full=probe(self.sampler,folder,request)
                            packet=full['observation'];packet['grid_weights']=weights
                            result=dict(status='pass',execution_status='completed',time=report['time'],steps=report['steps'],dt=dt,
                                quadrature=list(self.temporal_quadrature),initial_sha256=exact_vector_hash(report['accepted_states'][0]),reference_sha256=exact_vector_hash(initial['reference']),
                                final_spatial_alpha_zero=residual_metrics(packet['audit_alpha_zero']),final_stage_residual=report['stage_audits'][-1],
                                refreshed_trace_defect_max=max(s['refreshed_trace_defect'] for s in report['pimple_attempts']),
                                source_splitting_defects={k:max(abs(s.get(k,0)) for s in report['pimple_attempts']) for k in ('momentum_splitting_identity_defect','k_same_state_source_storage_defect','epsilon_same_state_source_storage_defect')},
                                extrema=audited_extrema(report),drift=representation(ref,packet,weights))
                            result['initial_outer_epsilon_perturbation_over_epsilon_H']=self.epsilon_perturbation
                            if self.epsilon_perturbation:
                                unperturbed=load(self.r/('preservation-temporal-Q0-'+ref.model)/(tag+'-dt'+str(j))/'state-packet.json')['observation']
                                result['perturbation_response_vs_unperturbed']=compare_samples(packet,dict(unperturbed,grid_weights=weights),weights)
                            write_json(folder/'endpoint-samples.json',packet);write_json(folder/'temporal-metrics.json',result)
                            self.admission.seal(folder,['case.xml','result.json','case.evolution.json','temporal-metrics.json','endpoint-samples.json'])
                            (folder/'packet.json').unlink();return result
                        self.check(name,trajectory)
                        folder=self.output/name
                        if (folder/'endpoint-samples.json').exists():runs.append((load(folder/'endpoint-samples.json'),load(folder/'temporal-metrics.json')))
                    table=dict(model=model,kappa_policy=policy,wall=wall,status='fail',fields={})
                    if len(runs)==3:
                        if len({m['initial_sha256'] for _,m in runs})!=1:raise AssertionError('temporal triplet initial coefficient identities differ')
                        ab=compare_samples(runs[0][0],runs[1][0],weights);bc=compare_samples(runs[1][0],runs[2][0],weights)
                        for field in ('u','v','k','epsilon_outer','Pi'):
                            d1=ab['fields'][field]['weighted_RMS'];d2=bc['fields'][field]['weighted_RMS']
                            uncertainty=1e-10
                            measured=min(d1,d2)>5*uncertainty
                            order=math.log2(d1/d2) if measured else None
                            table['fields'][field]=dict(differences_weighted_RMS=[d1,d2],differences_L2=[ab['fields'][field]['L2'],bc['fields'][field]['L2']],
                                order=order,uncertainty_floor=uncertainty,status='unresolved_floor' if not measured else 'pass' if .8<=order<=1.2 else 'fail')
                        table['status']='fail' if any(x['status']=='fail' for x in table['fields'].values()) else 'pass'
                        table['convergence_status']='failed' if table['status']=='fail' else 'partially_unresolved' if any(x['status']=='unresolved_floor' for x in table['fields'].values()) else 'demonstrated'
                        table['initial_outer_epsilon_perturbation_over_epsilon_H']=self.epsilon_perturbation
                        table['initial_state_sha256']=runs[0][1]['initial_sha256'];table['terminal_time']=self.horizon
                    tables.append(table);self.result['temporal_tables']=tables;self.save()
        self.result['expected_check_names']=expected;self.result['coverage']=coverage(self.result['checks'],set(expected),key=lambda x:x['name'])
    def finish(self):
        statuses=[x['status'] for x in self.result['checks']]
        statuses.extend(t['status'] for t in self.result.get('temporal_tables',[]))
        statuses.extend(t['status'] for t in self.result.get('spatial_tables',{}).values())
        self.result['status']='pass' if statuses and all(s=='pass' for s in statuses) and self.result.get('coverage',{}).get('status')=='pass' else 'fail'
        if self.result['status']=='pass' and any(t.get('convergence_status')=='partially_unresolved' for t in self.result.get('temporal_tables',[])):self.result['status']='unresolved'
        self.result['execution_status']='completed' if self.result.get('coverage',{}).get('status')=='pass' and all(x['status']=='pass' for x in self.result['checks']) else 'failed_or_incomplete'
        if self.result['execution_status']=='completed':self.result.pop('reason',None)
        self.save();return 0 if self.result['status']=='pass' else 2 if self.result['status']=='unresolved' else 1

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--resume',action='store_true');p.add_argument('--new-check-limit',type=int);p.add_argument('--runtime-directory',type=Path);p.add_argument('--horizon',type=float,default=.0016);p.add_argument('--boundary-policy',default='Legacy',choices=('Legacy','FrozenReferenceV1'));p.add_argument('--boundary-identity',default='native-neutral-reference');p.add_argument('--phase',choices=('spatial','temporal'),required=True);p.add_argument('--models',nargs='+',choices=MODELS,default=list(MODELS));p.add_argument('--label',default='');p.add_argument('--policies',nargs='+',choices=POLICIES,default=list(POLICIES));p.add_argument('--walls',nargs='+',choices=WALLS,default=list(WALLS));p.add_argument('--epsilon-perturbation',type=float,default=0.);p.add_argument('--quadrature-index',type=int,choices=(0,1,2),default=2);a=p.parse_args()
    c=Campaign(a)
    try:getattr(c,a.phase)()
    except IncompleteCampaign as error:
        c.result['status']='incomplete';c.result['execution_status']='incomplete';c.result['reason']=str(error);c.save();return 2
    except Exception as error:c.result['checks'].append(dict(name='orchestration',status='fail',failure=repr(error)));c.save()
    return c.finish()
if __name__=='__main__':sys.exit(main())
