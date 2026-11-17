#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify full continuous-MMS steps with distinct moving cutoff surfaces."""
import argparse,json,math,os,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from scipy.optimize import brentq
import coupled_oracle as core
import coupled_mms_check as mms
import history_recovery_check as hr
import physical_cache_history_check as history
import implicit_algorithm_check as fields
require,read = core.require,history.read
TABLE = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'


def run(manifest,folder,table,deadline,out):
    require(manifest['schema'] == 'WJH27MovingCutoffMMSV1','wrong crossing manifest')
    root = hr.resolve(folder,manifest['root']);require(hr.identity(table)['sha256'] == TABLE,'reference changed')
    geometry_path,geometry_id = history.input_path(manifest['geometry'],root);geometry = read(geometry_path)
    face = geometry['cutoff_face'];require(0 < face < .4*.15,'crossing face is not below the collar')
    labels = {'poly8','poly12','hist8','hist12'};require(set(manifest['cases']) == labels,'incomplete variant/quadrature matrix')
    cases,results = {},{};ref = history.po.Reference(read(table))
    for label,spec in manifest['cases'].items():
        require(time.monotonic() < deadline,'crossing check exhausted admitted cap')
        path,ident = history.input_path(spec['native'],root);packet = read(path);data = packet['state']
        variant = 'Historical' if label.startswith('hist') else 'PolynomialScalarsV1'
        require(data['source_identity'] == manifest['source_identity'] and data['reference_identity'] == TABLE and
            data['scalar_quadrature_policy'] == 'CubicCutoffFacePanelsV1','crossing source/policy binding differs')
        require(packet['normal_main'] is True and packet['steady_mode'] is False and data['manufactured'] is True and
            data['mms_variant'] == variant and data['mms_body_amplitude'] == 1 and data['gauge'] is True and
            data['step'] == 2 and data['time'] == .04 and data['dt'] == .02,'crossing MMS case differs')
        require(data['implicit_solver'] == 'CoupledJFNKV1' and data['coupled_newton']['preconditioner_policy'] == 'CurrentDataDriftV1',
            'crossing solver/cache recipe differs')
        q = (8,12) if label.endswith('8') else (12,18)
        require((data['quadratureX'],data['quadratureR']) == q and len(data['elements']) == 8,'crossing quadrature/mesh differs')
        target = mms.IndependentMMS(read(table),1,variant);s = core.state(data['coefficients'],data['sizes'])
        oracle = core.Oracle(data,target,20,24,deadline);owners = []
        for e in data['elements']:
            xl,xr,rl,rr = e['bounds'];x,r,w = oracle.rule(e,s,data['time']);area = (xr-xl)*(rr-rl)
            require(np.all(w > 0) and abs(float(np.sum(w))-area) <= 1e-12*area,'independent crossing area/positive weights fail')
            if min(abs(rl-face),abs(rr-face)) > 1e-12: continue
            current = oracle.cutoff_x_splits(e,s,None,(xl,xr,face,face))
            exact = oracle.cutoff_x_splits(e,None,data['time'],(xl,xr,face,face))
            current = [x for x in current if abs(float(oracle.cutoff(oracle.trace(e,'k',np.asarray([x]),s))[0])-face) <= 1e-12*target.H]
            exact = [x for x in exact if abs(float(oracle.cutoff(target.functions['K'](x,target.H,data['time'])[0]))-face) <= 1e-12*target.H]
            if current or exact:
                require(len(current) == len(exact) == 1 and xl < current[0] < xr and xl < exact[0] < xr,'current/exact crossing owner inventory differs')
                owners.append(dict(owner=e['geometry'],bounds=e['bounds'],current_root=current[0],exact_root=exact[0],
                    root_difference=abs(current[0]-exact[0]),points=len(w),area=area))
        require(len(owners) == 2 and max(row['root_difference'] for row in owners) > 1e-10,
                'both adjacent owners with distinct current/exact surfaces must be exercised')
        report_path,report_id = history.input_path(spec['independent'],root);report = read(report_path)
        require(report['native_sha256'] == ident['sha256'] and report['table_sha256'] == TABLE and
            report['native_source_identity'] == manifest['source_identity'] and report['current_coupled_step_gate_pass'] is True and
            report['independent_quadrature_pass'] is True and report['nonzero_prescribed_body_coverage_pass'] is True,
            'full independent crossing equations/history/body gate failed')
        if label.endswith('8'):
            require(report['direct_reference_gate_pass'] is True,'crossing base requires an independently refreshed dense nonlinear reference')
        cases[label] = (data,s)
        results[label] = dict(native=ident,independent=report_id,owners=owners,current_equations=report['status'],
            direct_reference=report['direct_nonlinear']['status'],current_residual=report['independent_assemblies'][-1]['current_residual'])
    # Direct c(K_exact(x,t))-face roots over the complete declared interval.
    oracle = core.Oracle(cases['poly8'][0],mms.IndependentMMS(read(table),1,'PolynomialScalarsV1'),20,24,deadline)
    movement = []
    for when in (0.,.04,.2,.4):
        def residual(x):return float(oracle.cutoff(oracle.m.functions['K'](x,oracle.m.H,when)[0]))-face
        require(residual(0)*residual(1) < 0,'frozen face loses its interior exact crossing')
        point = brentq(residual,0,1,xtol=2e-14,rtol=4*np.finfo(float).eps)
        movement.append(dict(time=when,x=point,residual=residual(point)))
    require(abs(movement[0]['x']-movement[-1]['x']) >= .2,'exact moving crossing lacks meaningful support')
    differences = {}
    for name in ('poly','hist'):
        a,sa = cases[name+'8'];b,sb = cases[name+'12']
        for key in fields.EXACT_KEYS:
            if key not in ('quadratureX','quadratureR'):require(a[key] == b[key],'crossing q control changed '+key)
        delta = fields.coefficient_delta(sa,sb,a['sizes'],ref)
        require(all(row['scaled_MAX'] <= 1e-8 for row in delta.values()),'crossing q field difference exceeds fixed1e-8')
        differences[name] = delta
    out.update(status='PASS_MOVING_CUTOFF_COUPLED_MMS_ONLY',foundation_gate_pass=True,geometry=geometry_id,cases=results,
        exact_crossing_movement=movement,quadrature_coefficient_comparisons=differences,
        scope='Two continuous MMS families on actual crossing owners, with separate current/exact surfaces. Not a temporal-order, physical-reference, long-pilot or Task3 certificate.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','table','out'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True);a=p.parse_args();began=time.monotonic();require(not a.out.exists(),'crossing output is immutable')
    out=dict(schema='WJH27MovingCutoffMMSCheckV1',status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,
        manifest=hr.identity(a.manifest),checker=hr.identity(__file__),oracle=hr.identity(core.__file__),continuous_mms=hr.identity(mms.__file__),workers=1)
    try:
        require(sys.platform.startswith('linux') and math.isfinite(a.seconds) and 0 < a.seconds <= 60,'bounded Linux crossing assessment required')
        soft,hard=resource.getrlimit(resource.RLIMIT_AS);cap=min([8*1024**3]+[v for v in(soft,hard)if v != resource.RLIM_INFINITY]);resource.setrlimit(resource.RLIMIT_AS,(cap,hard))
        run(read(a.manifest),a.manifest.resolve().parent,a.table.resolve(),began+a.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds']=time.monotonic()-began
    with a.out.open('x')as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({key:out[key]for key in('status','reason','elapsed_seconds')if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__=='__main__':main()
