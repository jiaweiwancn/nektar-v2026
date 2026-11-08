#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native frozen-stress scan and archived-stage attribution; zero PDE solves."""
import argparse
import csv
from datetime import datetime,timezone
import math
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from diagnose_cls_accepted_covariance import cls_tensor,STAGES
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,write_json

KAPPA=.4405539051866587


def tensor(k,e,g,beta=1.,gamma=0.,reference_shear=0.):
    R=cls_tensor(k,e,g);S=(g+g.T)*k/e;W=(g-g.T)*k/e
    q=max(np.linalg.norm(S),np.linalg.norm(W))/2**.5
    C=.3*(1. if q>8 else -math.expm1(-.36*math.exp(.75*q)))/(1+.35*q**1.5)
    ss=np.sum(S*S);oo=np.sum(W*W)
    cubic=k*C*(-10*C*C*(S@S@W-W@S@S)+(-5*C*C*ss+5*C*C*oo)*S)
    ref=np.zeros((3,3));ref[0,2]=ref[2,0]=reference_shear
    return R+(beta-1)*cubic-gamma*(g+g.T-ref)


def plot(out,operators,destination):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyse_cls_accepted_awf import dense_column
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    ax=axes[0,0]
    for name in ('refined','coarse','one'):
        s=out['cases'][name]['stages']
        ax.plot(range(5),[s[n]['minimum_outer'] for n in STAGES],marker='o',label=name+(' (mixed-side audit)' if name=='one' else ''))
    ax.axhline(0,color='k',lw=.6);ax.set_xticks(range(5),STAGES,rotation=15);ax.set_ylabel('Minimum covariance');ax.legend(fontsize=8)
    ax=axes[0,1]
    for law in ('Original','Regularized'):
        rows=[r for r in out['tangent_results'] if r['id']=='refined-post_flow-audit-'+law]
        ax.semilogx([r['h'] for r in rows],[r['min_work'] for r in rows],marker='o',label=law)
    ax.axhline(0,color='k',lw=.6);ax.set_xlabel('Physical gradient perturbation h');ax.set_ylabel('Minimum sampled tangent work');ax.legend()
    ax=axes[1,0];summary=load(operators/'summary.json')
    labels=sorted((k for k in summary['cases'] if k.startswith('Regularized-')),key=lambda k:summary['cases'][k]['time'])
    times=[summary['cases'][k]['time'] for k in labels]
    for b,name in enumerate(('u','w','k','epsilon')):
        ax.semilogy(times,[summary['cases'][k]['blocks'][b]['common_normalized'][0] for k in labels],marker='o',label=name)
    ax.set_xlabel('Physical time');ax.set_ylabel('Selected spatial residual / common scale');ax.legend()
    ax=axes[1,1]
    for label in labels:
        q=load(Path(summary['cases'][label]['full_free_loads_retained']))['diagnostics']['physical'];x,v=dense_column(q,-3.5)
        ax.plot(v['k'],v['z'],label='T='+format(q['time'],'.2g'))
    ax.set_ylim(0,.1);ax.set_xlabel('Transported k');ax.set_ylabel('z');ax.set_title('Cutoff-refined Regularized/Linear; x='+format(x,'.3g'));ax.legend()
    fig.savefig(destination,dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stages',type=Path,required=True)
    p.add_argument('--tangent',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--operators',type=Path,help='optional immutable strict operator readbacks for a compact scientific figure');a=p.parse_args()
    os.environ.update(THREADS);a.output.mkdir(exist_ok=False)
    queries={};rows=[];out=dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),governing_solves=0,
        tangent_executable_sha256=sha(a.tangent),scope='finite 721-angle x-z incompressible rank-one physical one-sided work scan at frozen k/epsilon; not a coupled stability theorem or a global positive-direction certificate',cases={})
    for case in ('refined','coarse','one'):
        packet=load(a.stages/case/'archived.json');data=load(a.stages/case/'stages/result.json')['diagnostics'];offset=[int(x) for x in data['committed']['offsets']]
        histories={stage:{str(b):state_hash(packet[stage][offset[b]:offset[b+1]]) for b in range(len(offset)-1)} for stage in STAGES}
        require(all(histories['committed'][str(b)]==histories['post_flow'][str(b)] for b in (2,3)), 'flow stage changed turbulence coefficients')
        entry=out['cases'][case]=dict(field_hashes=histories,flow_preserves_k_and_epsilon_bits=True,stages={},source_result_sha256=sha(a.stages/case/'stages/result.json'))
        for stage,d in data.items():
            w=d['minimum_outer'];entry['stages'][stage]=dict(minimum_outer=w['selected']['eigenvalues'][0],witness=w,
                one_sided_neighbors=d['one_sided_witnesses'],scalar_budget=d['scalar_budget'])
            if stage not in ('committed','post_flow','corrected'):continue
            choices=[('audit',w)]+[('physical-'+str(int(x['geometry'])),x) for x in d['one_sided_witnesses'] if not x['inner']]
            for role,point in choices:
                g=np.array(point['velocity_gradient']);k,e=point['k'],point['epsilon'];z=point['coordinate'][2]
                for law,beta,gamma in [('Original',1.,0.),('Regularized',.4,KAPPA*.4*z)]:
                    name=case+'-'+stage+'-'+role+'-'+law
                    expected=np.array(point['same_state_'+law]['tensor'])
                    require(np.max(abs(tensor(k,e,g,beta,gamma,.4/(KAPPA*z))-expected))<2e-12,'independent same-state variant tensor differs')
                    queries[name]=dict(case=case,stage=stage,role=role,law=law,k=k,epsilon=e,G=g.tolist(),distance=z,beta=beta,gamma=gamma,reference_shear=.4/(KAPPA*z),expected_tensor=expected.tolist())
        if case in ('refined','coarse'):
            left=data['committed']['minimum_outer'];right=data['post_flow']['minimum_outer']
            require(left['coordinate']==right['coordinate'] and left['geometry']==right['geometry'],
                'entrywise local attribution needs the same native one-sided witness')
            before=np.array(left['velocity_gradient']);after=np.array(right['velocity_gradient']);delta=after-before
            require(left['k']==right['k'] and left['epsilon']==right['epsilon'],'local gradient attribution requires frozen scalars')
            entry['gradient_attribution']=dict(scope='pointwise hybrid gradients, not represented coupled solutions or additive eigenvalue contributions',cases={})
            masks={'normal_diagonal':[(0,0),(2,2)],'both_shears':[(0,2),(2,0)],
                'G00':[(0,0)],'G02':[(0,2)],'G20':[(2,0)],'G22':[(2,2)]}
            for role,indices in masks.items():
                g=before.copy()
                for i,j in indices:g[i,j]+=delta[i,j]
                k,e=left['k'],left['epsilon'];z=left['coordinate'][2];R=tensor(k,e,g)
                name=case+'-hybrid-'+role+'-Original'
                queries[name]=dict(case=case,stage='hybrid',role=role,law='Original',k=k,epsilon=e,G=g.tolist(),
                    distance=z,beta=1.,gamma=0.,reference_shear=.4/(KAPPA*z),expected_tensor=R.tolist())
                entry['gradient_attribution']['cases'][role]=dict(covariance_minimum=float(np.linalg.eigvalsh(R)[0]),
                    divergence=float(np.trace(g)),tensor_trace=float(np.trace(R)),trace_minus_2k=float(np.trace(R)-2*k))
            z=data['corrected']['minimum_outer']['coordinate'][2];k=.5105080130598861;e=.4**3/(KAPPA*z);g=np.zeros((3,3));g[0,2]=.4/(KAPPA*z)
            for law,beta,gamma in [('Original',1.,0.),('Regularized',.4,KAPPA*.4*z)]:
                name=case+'-analytical-neutral-'+law;queries[name]=dict(case=case,stage='analytical-neutral',role='reference',law=law,k=k,epsilon=e,G=g.tolist(),distance=z,beta=beta,gamma=gamma,reference_shear=g[0,2],expected_tensor=tensor(k,e,g,beta,gamma,g[0,2]).tolist())
    source=a.output/'input.csv'
    with source.open('w',newline='') as stream:
        writer=csv.writer(stream,lineterminator='\n');writer.writerow(['id','distance','k','epsilon','beta','gamma','reference_shear',*['g'+str(i)+str(j) for i in range(3) for j in range(3)]])
        for name,q in queries.items():writer.writerow([name,q['distance'],q['k'],q['epsilon'],q['beta'],q['gamma'],q['reference_shear'],*np.array(q['G']).ravel()])
    command=[str(a.tangent.resolve()),'--general',str(source.resolve()),str((a.output/'native.csv').resolve())];start=time.monotonic()
    q=subprocess.run(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300);(a.output/'native.log').write_bytes(q.stdout);require(q.returncode==0,'native tangent scan failed')
    for row in csv.DictReader((a.output/'native.csv').open()):
        name=row['id'];p=queries[name];h=float(row['h']);angle=float(row['angle']);side=int(row['side']);normal=np.array([math.cos(angle),0.,math.sin(angle)]);polarization=np.array([-math.sin(angle),0.,math.cos(angle)])
        g=np.array(p['G']);delta=np.outer(polarization,normal);base=tensor(p['k'],p['epsilon'],g,p['beta'],p['gamma'],p['reference_shear']);trial=tensor(p['k'],p['epsilon'],g+side*h*delta,p['beta'],p['gamma'],p['reference_shear'])
        work=-float(np.sum(delta*(trial-base)))/(side*h);native=float(row['min_work']);matrix=np.array([[float(row['R'+str(i)+str(j)]) for j in range(3)] for i in range(3)])
        require(np.max(abs(matrix-base))<2e-12 and abs(work-native)<2e-7*max(1,float(row['nut'])),'independent physical tangent/work differs')
        row={k:float(v) if k!='id' else v for k,v in row.items()};row.update(independent_work=work,work_difference=abs(work-native),base_covariance_minimum=float(np.linalg.eigvalsh(matrix)[0]))
        rows.append(row)
    out.update(status='passed',queries=queries,tangent_results=rows,receipt=dict(command=command,seconds=time.monotonic()-start,input_sha256=sha(source),output_sha256=sha(a.output/'native.csv')),interpretation='A negative converged work in a sampled incompressible direction is a concrete frozen principal-part violation. Positive values on this finite scan do not prove all-direction or coupled-system stability. Covariance and tangent signs are reported separately.')
    if a.operators:
        figure=a.output/'covariance-and-evolution.png';plot(out,a.operators,figure)
        out['figure']=dict(path=str(figure),sha256=sha(figure))
    write_json(a.output/'native.json',out);print(a.output/'native.json')


if __name__=='__main__':main()
