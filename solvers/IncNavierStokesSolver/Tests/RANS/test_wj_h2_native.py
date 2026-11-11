#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Production-session open-boundary controls, strict passive checks and transactions."""
import argparse
import hashlib
import json
from pathlib import Path
import uuid
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from wj_h2_profiles import wall,NU,D0,TOP
from wj_h2_cases import case_xml
from cls_epawf_m3_cases import set_info,set_parameter,write
from validate_wj_h1 import run,require,save,sha
from compare_cls_epawf_pimple import same_history
from wj_h2_weak import independent_weak


def independent_boundary(report):
    d=report['open_boundary'];offsets=report['offsets'];n=len(report['state'])
    expected=np.zeros(n);measures={};mass=0;H=report['parameters'][1];h=report['parameters'][2]
    for r in d['boundary_rows']:
        normal=np.array(r['normal']);U=np.array(r['U']);G=np.array(r['G']);R=np.array(r['R']);k=r['k'];w=r['weight']
        role=r['role'];inner=r['inner'];un=U@normal;mass+=w*un
        measures[(role,inner)]=measures.get((role,inner),0)+w
        nt=wall(r['K'],H,h)['a']*max(r['chart'][2]-wall(r['K'],H,h)['zv'],0) if inner else 0
        T=-R+np.eye(3)*(2*k/3-r['pi'])+(NU*(G+G.T) if inner else 0)
        assert_allclose(r['actual_traction'],T@normal,rtol=2e-12,atol=3e-13)
        applied=np.array(r['applied_traction']);dk=NU+nt if inner else (572/6561)*k*k/r['epsilon']
        de=0 if inner else dk/1.3
        jk=-dk*np.dot(r['grad_k'],normal) if role=='inlet' else 0.
        je=-de*np.dot(r['grad_epsilon'],normal) if role=='inlet' and not inner else 0.
        if role=='top' and 'ReferenceTractionV1' in report['signature']:
            # The analytic control applies its nonzero outer epsilon flux.
            K=report['parameters'][5];je=(572/6561)*K*K/1.3/(TOP-D0)
        flux=np.array([U[0]*un-applied[0],U[2]*un-applied[2],k*un+jk,r['epsilon']*un+je])
        assert_allclose(r['flux'],flux,rtol=3e-12,atol=3e-13)
        for b in range(4):
            if b==3 and inner:continue
            ids=np.array(r['epsilon_ids'] if b==3 else r['ids'],dtype=int)+int(offsets[b])
            phi=r['epsilon_phi'] if b==3 else r['phi']
            expected[ids]+=w*np.array(phi)*flux[b]
    assert_allclose(d['boundary_load'],expected,rtol=3e-12,atol=3e-12)
    for role in ('inlet','outlet'):
        assert_allclose(measures[(role,True)],H,rtol=0,atol=2e-13)
        assert_allclose(measures[(role,False)],TOP-D0-H,rtol=0,atol=2e-12)
    assert_allclose(measures[('ground',True)],8,rtol=0,atol=2e-13)
    assert_allclose(measures[('top',False)],8,rtol=0,atol=2e-13)
    require(d['gauge_unknowns']==0 and d['boundary_flux_correction']==0,'open outlet retained gauge/ceiling repair')
    assert_allclose(mass,d['mass_flux'],rtol=0,atol=2e-13)
    return dict(max_load_error=float(np.max(abs(expected-d['boundary_load']))),mass_flux=mass,
        measures={str(k):v for k,v in measures.items()},rows=len(d['boundary_rows']))


def main():
    p=argparse.ArgumentParser();p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.output.resolve()/('run-'+uuid.uuid4().hex[:12]);root.mkdir(parents=True)
    summary=dict(status='STARTED',root=str(root),solver_sha256=sha(a.solver),probe_sha256=sha(a.probe),runs={})
    def native(label,steps=2,restart=None,reuse='Off',mutate=None,reason=None,reject=0):
        folder=root/label;xml=case_xml(folder,'analytic',steps,1e-4,reuse=reuse,restart=restart,reject=reject)
        if mutate:
            tree=ET.parse(xml);mutate(tree.getroot());write(tree.getroot(),xml)
        receipt=run([a.solver.resolve(),xml.resolve()],folder,timeout=600,expected_failure=reason)
        report=json.loads(xml.with_suffix('.evolution.json').read_text());summary['runs'][label]=dict(report=str(xml.with_suffix('.evolution.json')),receipt=receipt)
        return report,xml
    baseline,bxml=native('two')
    require(baseline['steps']==2 and baseline['failure']=='','analytic native control did not accept both steps')
    for r in baseline['pimple_attempts']:
        for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):
            require(r[k]==1,k)
        require(r['bound_activations']==0,'material correction in control')
    first,fx=native('one',1);second,sx=native('restart',1,first['checkpoint'])
    same_history(baseline,second)
    on,ox=native('reuse-on',2,reuse='ExactStepV1');same_history(baseline,on)
    for key in ('physical_samples','state_audits','wj_audit_packets'):
        require(baseline[key]==on[key],key+' cache changed physical/audit packets')
    immutable={str(f):sha(f) for f in (Path(first['checkpoint'])/'State.bin',fx.with_suffix('.evolution.json'))}
    frozen,zx=native('frozen',0,first['checkpoint']);same_history(first,frozen)
    require(frozen['initialization_lift_calls_iterations']==[0,0] and not frozen['pimple_attempts'],'read-only restart solved/lifted')
    rejected,rx=native('rejected',1,first['checkpoint'],reason='deliberate_PIMPLE_WJ_hard_audit_rejection',reject=2)
    same_history(first,rejected)
    require(rejected['rejected_state'] and rejected['rejected_state']!=first['state'],'rejected candidate not retained')
    require(all(sha(f)==h for f,h in immutable.items()),'transaction changed producer checkpoint/report')
    for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):
        require(rejected['pimple_attempts'][-1][k]==1,'rejected attempt accounting '+k)
    # Changes to actual physical boundary data must be signature-incompatible.
    def traction(r):r.find('./CONDITIONS/FUNCTION[@NAME="OutletTractionOuter"]/E[@VAR="u"]').set('VALUE','0.001')
    native('signature',0,first['checkpoint'],mutate=traction,reason='signature mismatch')
    def corner(r):r.find('./CONDITIONS/FUNCTION[@NAME="BoundaryValuesOuter"]/E[@VAR="epsilon"]').set('VALUE','0.1')
    native('corner',0,mutate=corner,reason='open_inlet_epsilon_corner_conflict')
    def backflow(r):r.find('./CONDITIONS/FUNCTION[@NAME="InitialConditionsOuter"]/E[@VAR="u"]').set('VALUE','-0.1')
    native('backflow',0,mutate=backflow,reason='open_boundary_backflow')
    checks=[]
    for label,report,xml,reuse in [('startup',None,bxml,False),('final',on,ox,True)]:
        if report is None:
            report,xml=native('startup',0,reuse='ExactStepV1')
        packet=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(xml.with_suffix('.evolution.json')),
            native_report_sha256=sha(xml.with_suffix('.evolution.json')),native_checkpoint=report['checkpoint'],
            checkpoint_sha256=sha(Path(report['checkpoint'])/'State.bin'),assembly=True,reuse=reuse,weak=True)
        inp=root/(label+'-checks-input.json');out=root/(label+'-checks.json');save(inp,packet)
        observer=root/(label+'-observer');observer.mkdir()
        run([a.probe.resolve(),xml.resolve(),inp,out],observer,timeout=900)
        c=json.loads(out.read_text());require(c['initialization_work']==[0,0] and c['governing_solves']==0,'passive check performed governing work')
        require(c['assembly_checks']['pressure_modes'][2]['action_max']>1e-4,'pressure constant still a nullspace')
        checks.append(dict(path=str(out),sha256=sha(out),boundary=independent_boundary(report),assembly=c['assembly_checks'],weak=independent_weak(c['weak_checks'],report)))
    summary.update(status='PASS',native_attempts=7,checks=checks,
        exact_off_on_and_restart=True,frozen_and_rejected_rollback=True,unchanged_producer=immutable)
    save(root/'summary.json',summary);print(json.dumps(dict(status='PASS',root=str(root),native_attempts=7)))

if __name__=='__main__':main()
