#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned high-order observations using a frozen producing-study engine.

No native attempt is repeated to repair quadrature. Existing case argv, XML,
code, runtime and output receipts remain bound to the frozen source checkout.
New governing work, if authorized by the screen gate, uses that same engine.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--frozen-source',type=Path,required=True);p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--retained',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True)
    p.add_argument('--mode',choices=('screen','select','full','transactions','report'),required=True);p.add_argument('--mesh',choices=('M1','M2'),default='M1');a=p.parse_args()
    frozen=a.frozen_source.resolve()/'solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(frozen))
    from run_wj5_native import Study as FrozenStudy
    from wj_immutable import artifacts,write_once,read,file_hash,digest,verify_artifacts,file_identity
    from wj_qualification import aggregate,component
    from run_cls_epawf_picard import require
    own=Path(__file__).resolve().parent;spec=importlib.util.spec_from_file_location('wj5_qualified_estimates',own/'wj5_measure.py');new_measure=importlib.util.module_from_spec(spec);spec.loader.exec_module(new_measure)
    a.solver=a.solver.resolve();a.probe=a.probe.resolve();summary=a.output.resolve()/a.policy/'summary.json'
    archive=None
    if summary.exists():
        data=summary.read_bytes();target=summary.parent/'summary-lineage'/(file_hash(summary)+'.json');target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():require(target.read_bytes()==data,'archived summary differs')
        else:
            with target.open('xb') as stream:stream.write(data)
        archive=file_identity(target)
    class Study(FrozenStudy):
        def measure(self,*args,order=8,xorder=4,**kwargs):
            return super().measure(*args,order=order,xorder=xorder,**kwargs)
        def observe_arithmetic(self,name,d,m,c):
            c['arithmetic']=self.arithmetic(name,d,m)
            c['boundary_uncertainty']=new_measure.boundary_uncertainty(m,c['arithmetic'],c['rule_settings']);return c
    s=Study(a)
    extra=artifacts([Path(__file__),own/'wj5_measure.py']);s.measure_code=dict(s.measure_code,measurement_orchestrator=extra)
    s.measure_version=digest(s.measure_code)[:12]
    design_path=s.root/('measurement-design-'+s.measure_version+'.json')
    design=dict(schema='WJ5ObservationRefinementV1',parent_design_sha256=file_hash(s.root/'design.json'),frozen_source=str(a.frozen_source.resolve()),
        frozen_revision=subprocess.check_output(['git','-C',str(a.frozen_source),'rev-parse','HEAD'],text=True).strip(),
        orchestrator_code=extra,baseline_order_z=8,baseline_order_x=4,independent_order_z=10,independent_order_x=6,
        absolute_floor=1e-10,reason='6->8 reference-error integration was unresolved; polynomial drift integration remained at roundoff; no governing rerun',
        native_mesh_dt_target_policy='unchanged original frozen design and exact cached native intents')
    write_once(design_path,design);s.s.setdefault('measurement_lineages',{})[s.measure_version]=dict(design=file_identity(design_path),parent_snapshot=archive,code=s.measure_code)
    s.s['components']['measurement_requalification']=component('NOT_RUN','final all-observation constraint-aware qualification required')
    try:
        if a.mode=='screen':
            mesh=a.mesh
            if mesh=='M2':
                siblings=[read(a.output.resolve()/q/'summary.json') for q in ('EP','Budget')]
                require(any(t.get('screens',{}).get('M1',{}).get('status')!='PASS' for t in siblings),'M2 escalation is not authorized')
            name='screen-'+mesh
            z,r,f=s.native(name+'-cp0',mesh=mesh);s.valid(z,r)
            require(len(z['state'])==dict(M1=2447,M2=2975)[mesh],'unexpected native custom mesh size')
            m,c0=s.measure(name+'-cp0-q8',z,f,mesh=mesh,assembly=True)
            d,r,f=s.native(name+'-32',32,z['checkpoint'],mesh=mesh,arithmetic_time=.001024);s.valid(d,r)
            m,c=s.measure(name+'-32-q8',d,f,mesh=mesh,assembly=True);s.observe_arithmetic(name+'-32-q8',d,m,c)
            qm,qc=s.measure(name+'-32-q10x6',d,f,mesh=mesh,order=10,xorder=6,diagnostics=False)
            changes={};leaves={}
            for kind in ('continuum_error','initial_drift'):
                for scope in ('full','central','constrained_boundary'):
                    for field,row in c['errors'][kind][scope]['fields'].items():
                        key=kind+'/'+scope+'/'+field;error=abs(row['L2']-qc['errors'][kind][scope]['fields'][field]['L2']);changes[key]=error
                        leaves[key]=component('PASS' if error<=1e-10 else 'UNRESOLVED',value=error,limit=1e-10)
            quad=aggregate(leaves);gate=aggregate(dict(startup=c0['errors']['accuracy_gate'],terminal=c['errors']['accuracy_gate'],quadrature=quad))
            entry=dict(initial=c0,terminal=c,unknowns=len(d['state']),run_seconds=r['seconds'],seconds_per_step=r['seconds']/32,
                quadrature=dict(order=[8,10],xorder=[4,6],absolute_floor=1e-10,changes=changes,gate=quad),status=gate['status'],
                measurement_lineage=file_identity(design_path),no_governing_attempt_repeated=True)
            s.s.setdefault('screens',{})[mesh]=entry;s.s['components']['local_screen']=component(entry['status'],mesh=mesh)
            failed=s.s['components'].get('screen',{})
            if 'observation quadrature sensitivity exceeds declared absolute floor' in failed.get('reason',''):
                s.s.setdefault('resolved_harness_failures',[]).append(dict(original_component=s.s['components'].pop('screen'),parent_snapshot=archive,
                    resolution='complete screen packet retained; independent higher observation order; original accuracy obstruction unchanged'))
        elif a.mode=='select':
            selection=s.select()
            if not selection['selected']:
                s.s['components']['selected_horizontal']=component('OBSTRUCTED','neither predeclared supported local mesh qualified at Tx; frozen T3 remains unreached')
                s.s['components']['selected_timestep']=component('NOT_RUN','no qualified selected mesh')
                s.s['components']['transactions']=component('NOT_RUN','no qualified selected trajectory')
        elif a.mode!='report':getattr(s,a.mode)()
        s.finish()
    except Exception as e:
        s.s['components']['requalified_'+a.mode]=component('RESOURCE_LIMIT' if 'RESOURCE_LIMIT' in str(e) else 'FAIL',str(e));s.s.setdefault('failures',[]).append(dict(mode='requalified_'+a.mode,reason=str(e),code=extra));s.finish();raise
    finally:s.ledger.finish_session(s.session);s.save()

if __name__=='__main__':main()
