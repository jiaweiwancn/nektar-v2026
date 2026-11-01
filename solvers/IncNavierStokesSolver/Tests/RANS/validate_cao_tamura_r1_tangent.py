#!/usr/bin/env python3
"""Read-only 2D material-tangent diagnostic at immutable native checkpoints."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode=True
os.environ.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
import numpy as np
from numpy.polynomial import legendre as leg
from cao_tamura_experiment import ROOT,allocation,environment,execute,records,sha,write_csv,write_json
from record_cao_tamura_r1 import field_bound

GRAD=[f'g{i}{j}' for i in range(3) for j in range(3)]
INPUT=['id','eta','k','epsilon',*GRAD]
COARSE='min_transverse_tangent_h1e5'
FINE='min_transverse_tangent_h1e6'


def derivative(nodes):
    scale=(max(nodes)-min(nodes))/2
    x=(np.asarray(nodes)-np.mean(nodes))/scale
    weights=np.array([1./np.prod(x[i]-np.delete(x,i)) for i in range(len(x))])
    d=np.zeros((len(x),len(x)))
    for i in range(len(x)):
        for j in range(len(x)):
            if i!=j:d[i,j]=weights[j]/(weights[i]*(x[i]-x[j]))
        d[i,i]=-sum(d[i])
    return d/scale,x,scale


def samples(path,c):
    """Exact interpolation derivatives for the declared affine tensor grid."""
    groups={}
    for r in records(path):groups.setdefault(int(r['element']),[]).append(r)
    output=[];polynomial=0.;reconstruction=0.;gradient=0.
    for element,rows in sorted(groups.items()):
        rows.sort(key=lambda r:int(r['point']));n=math.isqrt(len(rows))
        if n*n!=len(rows) or [int(r['point']) for r in rows]!=list(range(n*n)):
            raise ValueError('Incomplete tensor-product element')
        x=np.array([float(r['x']) for r in rows]).reshape(n,n)
        z=np.array([float(r['y']) for r in rows]).reshape(n,n)
        if (np.max(abs(x-x[0,:]))>5e-13 or np.max(abs(z-z[:,0,None]))>5e-13):
            raise ValueError('This diagnostic requires the declared axis-aligned Plane mesh')
        dx,xi,sx=derivative(x[0,:]);dz,eta,sz=derivative(z[:,0])
        xx,zz=np.meshgrid(xi,eta)
        f=xx**2*zz**3+.2*xx*zz
        polynomial=max(polynomial,float(np.max(abs(sx*(f@dx.T)-(2*xx*zz**3+.2*zz)))),
                       float(np.max(abs(sz*(dz@f)-(3*xx**2*zz**2+.2*xx)))))
        velocities=[]
        vand=leg.legvander2d(xx,zz,[4,4]).reshape(n*n,25)
        for name in ('u','v'):
            u=np.array([float(r[name]) for r in rows]).reshape(n,n)
            ux=u@dx.T;uz=dz@u
            coeff=np.linalg.lstsq(vand,u.ravel(),rcond=None)[0].reshape(5,5)
            reconstruction=max(reconstruction,float(np.max(abs(leg.legval2d(xx,zz,coeff)-u)))/max(1.,float(np.max(abs(u)))))
            a=leg.legval2d(xx,zz,leg.legder(coeff,axis=0))/sx
            b=leg.legval2d(xx,zz,leg.legder(coeff,axis=1))/sz
            gradient=max(gradient,float(np.max(abs(a-ux)))/max(1.,float(np.max(abs(ux)))),
                         float(np.max(abs(b-uz)))/max(1.,float(np.max(abs(uz)))))
            velocities.append((ux,uz))
        selected={j*n+i for i in (n//2-1,n//2) for j in (n//2-1,n//2)}
        selected.update(q for q,r in enumerate(rows) if abs(float(r['x'])-c['xmin'])<1e-10
                        or any(abs(float(r['y'])-h)<1e-10 for h in (c['hm'],c['point_n'])))
        for q in sorted(selected):
            j,i=divmod(q,n);r=rows[q];g=np.zeros((3,3))
            g[0,0],g[0,2]=velocities[0][0][j,i],velocities[0][1][j,i]
            g[2,0],g[2,2]=velocities[1][0][j,i],velocities[1][1][j,i]
            for variant in ('actual','tracefree'):
                value=g.copy()
                if variant=='tracefree':
                    divergence=value[0,0]+value[2,2]
                    value[0,0]-=divergence/2;value[2,2]-=divergence/2
                row=dict(id=f'{element}:{q}:{variant}',eta=float(r['y']),k=float(r['k']),
                    epsilon=float(r.get('eps_tilde',r.get('eps'))),element=element,point=q,
                    x=float(r['x']),z=float(r['y']),variant=variant)
                row.update(zip(GRAD,value.ravel()));output.append(row)
    if polynomial>1e-10 or reconstruction>2e-10 or gradient>1e-8:
        raise ValueError(f'Gradient reconstruction failed: {polynomial}, {reconstruction}, {gradient}')
    return output,dict(polynomial_error=polynomial,p4_reconstruction_error=reconstruction,
        independent_gradient_error=gradient,elements=len(groups),samples=len(output))


def controls(exe,legacy,out):
    rows=[];old=[]
    reference=subprocess.check_output([str(exe),'--reference'],text=True).strip()
    q,kappa=map(float,reference.split(','))
    for i,shear in enumerate((0.,1.,10/3,8.)):
        r=dict(id=f'parallel{i}',eta=.02,k=.01,epsilon=.003)
        old.append(dict(r,shear=shear));r.update({g:0. for g in GRAD});r['g02']=shear;rows.append(r)
    r=dict(id='compatible-reference',eta=.1,k=.01,epsilon=.01/q)
    old.append(dict(r,shear=1.));r.update({g:0. for g in GRAD});r['g02']=1.;rows.append(r)
    write_csv(out/'control-general.csv',rows,INPUT)
    write_csv(out/'control-parallel.csv',old,['id','eta','k','epsilon','shear'])
    execute([exe,'cls',0,out/'control-general.csv',out/'control-general-result.csv'],out,'control-general')
    execute([legacy,'cls',0,out/'control-parallel.csv',out/'control-parallel-result.csv'],out,'control-parallel')
    a,b=[records(out/name) for name in ('control-general-result.csv','control-parallel-result.csv')]
    error=max(abs(float(x[k])-float(y[k]))/max(1e-12,abs(float(y[k]))) for x,y in zip(a,b) for k in (COARSE,FINE))
    if error>1e-8:raise ValueError('General-gradient probe disagrees with the legacy parallel oracle')
    return dict(status='PASS',parallel_relative_error=error,compatible_reference_q=q,
        compatible_reference_kappa=kappa,compatible_reference_highre_probe=next(r for r in a if r['id']=='compatible-reference'),
        reference_probe_scope='parallel canonical Q at arbitrary positive k/shear and zero resolved nu; not the complete reference trajectory')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('cases','kernel','legacy','build','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--cpu',default='79');a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    allocation(1,a.cpu)
    if len(a.cpu.split(','))!=1:raise ValueError('Use exactly one diagnostic worker')
    os.sched_setaffinity(0,{int(a.cpu)});build=environment(a.build)
    out=a.output.resolve()
    if out.is_relative_to(ROOT):raise ValueError('Store diagnostic data outside the source checkout')
    out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(__file__,out/'producer.py')
    cpp=Path(__file__).with_name('r1_general_stress_tangent.cpp');shutil.copy2(cpp,out/'producer.cpp')
    result=dict(status='RUNNING',scope='Local frozen constitutive tangents; not a developed case or global instability attribution',
        producer_sha256=sha(out/'producer.py'),cases={},extra_sources={str(cpp):sha(cpp),
            str(ROOT/'library/SolverUtils/Turbulence/RANSModels.cpp'):sha(ROOT/'library/SolverUtils/Turbulence/RANSModels.cpp'),
            str(ROOT/'library/SolverUtils/Turbulence/RANSModels.h'):sha(ROOT/'library/SolverUtils/Turbulence/RANSModels.h')},
        kernel=dict(path=str(a.kernel),sha256=sha(a.kernel)),legacy=dict(path=str(a.legacy),sha256=sha(a.legacy)),
        binaries={n:sha(build/n) for n in ('IncNavierStokesSolver','RANSWallFieldAudit')},
        sampling='four interior GLL points per element plus H/N/inlet traces, including one-sided copies; counts are not volume fractions')
    write_json(out/'evidence.json',result)
    try:
        exe=out/'general-tangent'
        compile_args=['c++','-std=c++17','-O3','-DNDEBUG','-Wall','-Wextra','-Wpedantic','-Werror',
            '-I'+str(ROOT/'library'),str(cpp),str(a.kernel),'-o',str(exe)]
        result['compile']=execute(compile_args,out,'compile',300);result['diagnostic_sha256']=sha(exe)
        result['controls']=controls(exe,a.legacy,out)
        summary=[]
        for label,item in json.loads(a.cases.read_text()).items():
            folder=Path(item['folder']);field=Path(item['field']);m=json.loads((folder/'manifest.json').read_text())
            start=json.loads((folder/'run-start.json').read_text())
            if (start['solver_sha256']!=result['binaries']['IncNavierStokesSolver']
                    or start['case_sha256']!=sha(folder/'case.xml') or m['case_sha256']!=start['case_sha256']
                    or start['geometry_sha256']!=sha(build/'geometry-lib/libSpatialDomains.so')):
                raise ValueError('Unbound native case/binary in '+label)
            c=m['config']
            if c['geometry']!='floor' or c['dimension']!=2:raise ValueError('Only declared affine 2D empty floors are supported')
            target=out/label;target.mkdir();binding=field_bound(field)
            initial=dict(stage_csv_sha256=sha(folder/'case.rans.csv'),manifest_sha256=sha(folder/'manifest.json'))
            execute([build/'RANSWallFieldAudit',folder/'case.xml',field,target/'full.csv'],target,'field')
            times=sorted({float(r['time']) for r in records(target/'full.csv')})
            if len(times)!=1 or not math.isfinite(times[0]) or ('time' in item and abs(times[0]-item['time'])>1e-10):
                raise ValueError('Wrong native diagnostic field time')
            data,checks=samples(target/'full.csv',c);write_csv(target/'samples.csv',data,list(data[0]))
            write_csv(target/'input.csv',[{k:r[k] for k in INPUT} for r in data],INPUT)
            nu=m['reference']['nu'] if c['diffusion']=='Molecular' else 0.
            metadata={r['id']:r for r in data};model_results={}
            for model in ('standard','cls'):
                file=target/(model+'.csv')
                execute([exe,model,nu,target/'input.csv',file],target,model,900)
                rows=records(file)
                if len(rows)!=len(data):raise ValueError('Incomplete material-tangent output')
                if model=='standard':
                    errors=[abs(float(r[FINE])-(nu+.09*metadata[r['id']]['k']**2/metadata[r['id']]['epsilon']))/
                        (nu+.09*metadata[r['id']]['k']**2/metadata[r['id']]['epsilon']) for r in rows]
                    if max(errors)>1e-6:raise ValueError('Standard transverse tangent does not equal nu+nu_t')
                    model_results['standard_control_error']=max(errors)
                for variant in ('actual','tracefree'):
                    selected=[r for r in rows if metadata[r['id']]['variant']==variant]
                    negative=[r for r in selected if max(float(r[COARSE]),float(r[FINE])) < -1e-12
                        and abs(float(r[COARSE])-float(r[FINE]))<=.02*max(abs(float(r[COARSE])),abs(float(r[FINE])))]
                    worst=min(selected,key=lambda r:float(r[FINE]))
                    value=dict(model=model,variant=variant,samples=len(selected),robust_negative=len(negative),
                        minimum=float(worst[FINE]),worst=dict(metadata[worst['id']],result=worst))
                    model_results[model+'/'+variant]=value
                    summary.append(dict(case=label,model=model,variant=variant,samples=len(selected),
                        robust_negative=len(negative),minimum=float(worst[FINE]),molecular_nu=nu))
            if (sha(folder/'case.rans.csv')!=initial['stage_csv_sha256']
                    or sha(folder/'manifest.json')!=initial['manifest_sha256'] or field_bound(field)!=binding):
                raise ValueError('Native diagnostic source changed during sampling')
            result['cases'][label]=dict(folder=str(folder),field=binding,**initial,
                time=times,resolved_nu=nu,physical_molecular_nu=m['reference']['nu'],run_start_sha256=sha(folder/'run-start.json'),
                gradient_checks=checks,results=model_results,
                source_sha256={str(f):sha(f) for f in (folder/'case.xml',target/'full.csv',target/'samples.csv',
                    target/'input.csv',target/'standard.csv',target/'cls.csv')})
        write_csv(out/'summary.csv',summary,list(summary[0]));result['compact_artifacts']={'summary.csv':sha(out/'summary.csv')}
        result['status']='COMPLETE'
    except Exception as error:
        result.update(status='FAIL',error=repr(error));write_json(out/'evidence.json',result);raise
    write_json(out/'evidence.json',result);print('Constitutive tangent diagnostic COMPLETE',flush=True)


if __name__=='__main__':main()
