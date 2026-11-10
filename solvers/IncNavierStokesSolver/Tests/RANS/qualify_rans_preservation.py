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

class Campaign:
    def __init__(self,args):
        self.r=args.root.resolve();self.phase=args.phase;self.output=self.r/('preservation-'+args.phase);self.output.mkdir(exist_ok=True)
        self.solver=self.r/'runtime/IncNavierStokesSolver';self.descriptor=self.r/'runtime/RANSNewModelsProbe';self.sampler=self.r/'runtime/RANSPreservationProbe'
        self.manifest=load(self.r/'runtime/source-build.json');self.result=dict(schema=1,phase=self.phase,checks=[],status='running',
            runtime_manifest_sha256=sha(self.r/'runtime/source-build.json'),runtime_source=self.manifest['revision'],ranks=1,threads=THREADS,
            source=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip())
        if sha(self.solver)!=self.manifest['solver']['sha256']:raise RuntimeError('frozen solver identity changed')
        os.environ.update(THREADS,LD_LIBRARY_PATH=str(self.r/'runtime/lib'),TMPDIR=str(self.r/'tmp'))
        self.save()
    def save(self):write_json(self.output/'results.json',self.result)
    def check(self,name,fn):
        start=time.monotonic()
        try:
            value=fn();status=value.get('status','pass')
            self.result['checks'].append(dict(name=name,status=status,seconds=time.monotonic()-start,result=value))
        except Exception as error:
            self.result['checks'].append(dict(name=name,status='fail',seconds=time.monotonic()-start,failure=repr(error)))
        self.save();print(json.dumps({k:self.result['checks'][-1][k] for k in ('name','status','seconds')}),flush=True)
    def config(self,folder,ref,knots,p=3,q=Q[0],steps=0,dt=.0001,perturb=0.):
        return case_xml(folder,ref.model,ref.policy,ref.wall,steps=steps,dt=dt,descriptor_probe=self.descriptor,radial_knots=knots,degree=p,
                        quadrature=q,checkpoint_interval=max(steps,1),cache='ExactStepV1',perturb=perturb)
    def frozen(self,label,ref,knots,p=3,q=Q[0],base=None,basis=False,grid=None):
        folder=self.output/label;self.config(folder,ref,knots,p,q)
        request=dict(operation='neutral',export_basis=basis,grid=grid or []) if base is None else replay_request(base,grid or [],basis)
        full=probe(self.sampler,folder,request);packet=full['observation']
        metrics=dict(status='pass',execution_status='completed',physical_admissibility=full['physical_admissibility'],
            state_sha256=exact_vector_hash(packet['state']),reference_sha256=exact_vector_hash(full['reference']),
            actual_quadrature=packet['actual_quadrature'],residuals=residual_metrics(packet['audit_alpha_zero']),
            native_alpha_zero=packet['audit_alpha_zero'],source_XML_sha256=sha(folder/'case.xml'))
        if basis:metrics['analytical_weak_balance']=weak_balance(ref,packet)
        if grid:metrics['representation']=representation(ref,packet,self.grid_weights)
        write_json(folder/'metrics.json',metrics)
        if basis:archive_basis(folder)
        else:(folder/'packet.json').unlink()
        return full,metrics
    def spatial(self):
        # A common evaluation grid spanning every refinement and reference breakpoint.
        refs=[Reference(m,p,w) for m in MODELS for p in POLICIES for w in WALLS]
        allknots=sorted(set(v for r in refs for v in nested(r.knots(),2)))
        grid,self.grid_weights=physical_grid(allknots,10,8)
        write_json(self.output/'evaluation-grid.json',dict(grid=grid,weights=self.grid_weights,order=10,xorder=8,panels=allknots))
        expected=[]
        for ref in refs:
            tag='-'.join((ref.model,ref.policy,ref.wall));expected.append(tag+'-reference')
            def initial(ref=ref,tag=tag):
                base,metrics=self.frozen(tag+'-h0',ref,ref.knots(),basis=True,grid=grid)
                write_json(self.output/(tag+'-base.json'),compact_packet(base))
                status='pass' if metrics['analytical_weak_balance']['status']=='pass' else 'fail'
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
        for wall in WALLS:
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
            # Preserve the full native arrays rather than guessing their labels.
            qa=byname['q1']['native_alpha_zero'];qb=byname['q2']['native_alpha_zero']
            qchanges={k:max((abs(a-b) for a,b in zip(v,qb[k],strict=True)),default=0.) for k,v in qa.items() if k in qb and len(v)==len(qb[k]) and ('norm' in k or 'dual' in k)}
            tables[tag]=dict(representation=data,quadrature_native_last_changes=qchanges,levels=byname,
                status='pass' if all(x.get('analytical_weak_balance',{}).get('status','pass')=='pass' for x in m) else 'fail')
        self.result['spatial_tables']=tables;self.save()
    def temporal(self):
        # The first explicit temporal study uses the same fixed aligned p3 spaces.
        # Quadrature is Q2, never inferred from changing only normal order.
        spatial=load(self.r/'preservation-spatial/evaluation-grid.json');grid=spatial['grid'];weights=spatial['weights']
        self.grid_weights=weights;expected=[];tables=[]
        for model in MODELS:
            for policy in POLICIES:
                for wall in WALLS:
                    ref=Reference(model,policy,wall);tag='-'.join((model,policy,wall));runs=[]
                    for j,dt in enumerate((1e-4,5e-5,2.5e-5)):
                        name=tag+'-dt'+str(j);expected.append(name)
                        def trajectory(j=j,dt=dt,ref=ref,tag=tag):
                            folder=self.output/(tag+'-dt'+str(j));steps=round(.0016/dt)
                            env_value=str(folder/'timing.jsonl');os.environ['NEKTAR_RANS_TIMING_FILE']=env_value
                            row=run_case(self.solver,folder,ref.model,ref.policy,ref.wall,steps=steps,dt=dt,descriptor_probe=self.descriptor,
                                radial_knots=ref.knots(),degree=3,quadrature=Q[2],cache='ExactStepV1',checkpoint_interval=steps)
                            if row['status']!='pass':return dict(status='fail',execution=row)
                            report=load(folder/'case.evolution.json')
                            # Reconstruct exactly the initializer on these same spaces.
                            initial=probe(self.sampler,folder,dict(operation='neutral',export_basis=False,grid=[]))
                            if exact_vector_hash(initial['observation']['state'])!=exact_vector_hash(report['accepted_states'][0]):raise AssertionError('probe/native initial coefficients differ')
                            request=replay_request(initial,grid,state=report['state'],time=report['time'])
                            full=probe(self.sampler,folder,request)
                            packet=full['observation'];packet['grid_weights']=weights
                            result=dict(status='pass',execution_status='completed',time=report['time'],steps=report['steps'],dt=dt,
                                initial_sha256=exact_vector_hash(report['accepted_states'][0]),reference_sha256=exact_vector_hash(initial['reference']),
                                final_spatial_alpha_zero=residual_metrics(packet['audit_alpha_zero']),final_stage_residual=report['stage_audits'][-1],
                                refreshed_trace_defect_max=max(s['refreshed_trace_defect'] for s in report['pimple_attempts']),
                                source_splitting_defects={k:max(abs(s.get(k,0)) for s in report['pimple_attempts']) for k in ('momentum_splitting_identity_defect','k_same_state_source_storage_defect','epsilon_same_state_source_storage_defect')},
                                extrema=audited_extrema(report),drift=representation(ref,packet,weights))
                            write_json(folder/'endpoint-samples.json',packet);write_json(folder/'temporal-metrics.json',result)
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
                        table['initial_state_sha256']=runs[0][1]['initial_sha256'];table['terminal_time']=.0016
                    tables.append(table);self.result['temporal_tables']=tables;self.save()
        self.result['expected_check_names']=expected;self.result['coverage']=coverage(self.result['checks'],set(expected),key=lambda x:x['name'])
    def finish(self):
        statuses=[x['status'] for x in self.result['checks']]
        statuses.extend(t['status'] for t in self.result.get('temporal_tables',[]))
        self.result['status']='pass' if statuses and all(s=='pass' for s in statuses) and self.result.get('coverage',{}).get('status')=='pass' else 'fail'
        self.save();return 0 if self.result['status']=='pass' else 1

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--phase',choices=('spatial','temporal'),required=True);a=p.parse_args()
    c=Campaign(a)
    try:getattr(c,a.phase)()
    except Exception as error:c.result['checks'].append(dict(name='orchestration',status='fail',failure=repr(error)));c.save()
    return c.finish()
if __name__=='__main__':sys.exit(main())
