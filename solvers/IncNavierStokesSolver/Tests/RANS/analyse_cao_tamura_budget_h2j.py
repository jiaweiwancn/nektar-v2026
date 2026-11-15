#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent enlarged-space replay and compact H2j numerical evidence."""
import argparse,time,os
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_h2i_convex import Context
from cao_tamura_h2h_scalar import Budget
from analyse_cao_tamura_budget_h2i import replay,tensor_check
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha,THREADS

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();os.environ.update(THREADS);start=time.monotonic();r=a.root;f=load(r/'prepare/frozen.json');ctx=Context(Path(f['context']));search=load(r/'joint-search/result.json');a.output.mkdir(exist_ok=False);count=0;proofs=[];native={}
    old=load(a.previous/'joint-context/context.json')
    for key in ('geometry_identity','inventory_identity','reference_sha256','defect_sha256','patch_signature'):require(old[key]==ctx.doc[key],'old physical identity mismatch '+key)
    for anchor in ('A','B'):
        P=ctx.problem(anchor,search['retention']['theta'])
        for ph in (v for v in search['phases'] if v['anchor']==anchor):
            initial=P.scan(np.zeros(P.dim),ph['result']['delta'],Budget(seconds=60));count+=1;proof=replay(P,ph,initial);actual=[]
            for m in proof['limiting_constraints']:
                m=dict(m)
                if m['kind']=='epsilon_bound':
                    outer=np.flatnonzero(~P.inner);m['point']=int(outer[m['hard_row']%len(outer)]);m['bound']='upper' if m['hard_row']<len(outer) else 'lower'
                if 'point' in m:
                    i=m['point'];m.update(chart=ctx.block('chart')[i].tolist(),physical=ctx.block('physical')[i].tolist(),geometry=int(P.geometry[i]),memberships=ctx.arrays['members'][ctx.arrays['members'][:,1]==i].tolist())
                    if m['kind']=='PSD':row,_=P.rayleigh_cut(i,m['vector'],ph['result']['delta']);m['new_velocity_cut_norm']=float(np.linalg.norm(row[41:P.nv]))
                actual.append(m)
            proof['limiting_constraints']=actual;proofs.append(proof)
        receipt=load(r/('native-'+anchor)/'receipt.json');d=load(Path(receipt['output']));candidate=load(r/'joint-search'/(anchor+'-actual')/'best-candidate.json');s=P.scan(np.array(candidate['joint_coefficients']),1e-4*P.K0,Budget(seconds=60));count+=1;seen=d['seen_sampling']
        for order,v in seen.items():
            for name in ('inner_covariance_witness','outer_covariance_witness'):
                w=v[name];o=w['observation'];ids=np.flatnonzero(np.all(ctx.block('chart')==o['coordinate'],axis=1)&(P.geometry==o['geometry']));require(len(ids)==1,'native witness not unique');check=tensor_check(w,P.nt[ids[0]] if o['inner'] else None);require(check['tensor_error']<1e-12,'native tensor replay failed');v[name+'_independent']=check
        au=d['audit'];native[anchor]=dict(state_sha256=candidate['state_sha256'],anchor_k_sha256=candidate['anchor_k_sha256'],passed=receipt['passed'],gates=receipt['gates'],native_covariance=au['accepted_awf_covariance'],native_inner_point=au['accepted_awf_inner_minimum_point'],native_outer_point=au['accepted_awf_outer_minimum_point'],native_inner_provenance=au['native_inner_minimum_geometry_element_side_block'],native_outer_provenance=au['native_outer_minimum_geometry_element_side_block'],native_GL16=d['native_divergence_GL16'],strong_continuity=au.get('strong_continuity'),scalar_budget=au['accepted_awf_scalar_domain'],coupled_norms=au['coupled_norms'],seen=seen,native_budget=d['native_budget'],reference_sha256=d['reference_sha256'],defect_sha256=d['defect_sha256'],output=dict(path=receipt['output'],sha256=sha(Path(receipt['output']))),new_velocity_coefficients=candidate['joint_coefficients'][41:P.nv])
    atomic(a.output/'obstruction.json',dict(classification=search['classification'],replays=proofs,independent_complete_scans=count,independent_seconds=time.monotonic()-start,scope='new 44-velocity/eight-epsilon space, exact fixed-k anchors and declared sampled physical bounds; not global infeasibility'))
    require(search['budget']['complete_scans']+count<=120 and search['budget']['seconds']+time.monotonic()-start<600,'aggregate replay budget exceeded')
    atomic(a.output/'native.json',dict(diagnostic_states=native,repeatable_checks=load(r/'candidate-checks.json'),fresh_hill_GL28=0,physical_hill_attempts=0,physical_hill_accepted_steps=0))
    flat={}
    prior=load(a.previous/'flat-H1/results.json')
    for model,item in load(r/'flat/results.json').items():
        d=load(Path(item['output']))['context'];flat[model]={k:d[k] for k in ('base_sha256','state_sha256','setup','scalar_domain','fields_changed','baseline_audit','scalar_H1_audit','joint_audit')}
        for key in ('baseline','scalar_H1','projected'):
            v=d[key];flat[model][key]={k:v[k] for k in ('selected_spatial_norms','analytical_free_load_norms','radial_errors','one_sided_derivative_jumps','radial_global_free_sum_error','global_components')};flat[model][key]['errors']=v['physical_reference']['flat_reference_errors'];flat[model][key]['signed_components']={}
            for name,q in v['radial_components'].items():
                term=name.split(':')[-1];flat[model][key]['signed_components'].setdefault(term,np.zeros(4));flat[model][key]['signed_components'][term]+=q['signed_dual_projection_on_total']
            flat[model][key]['signed_components']={k:v.tolist() for k,v in flat[model][key]['signed_components'].items()}
        require(d['scalar_H1']['state']==prior[model]['data']['projected']['state'],'old scalar H1 field changed')
        base=np.array(d['baseline']['state'],dtype='<f8');v=np.array(d['projected']['state'],dtype='<f8');oldv=np.array(d['scalar_H1']['state'],dtype='<f8');require(base[2280:].tobytes()==v[2280:].tobytes() and oldv[1520:2280].tobytes()==v[1520:2280].tobytes(),'flat epsilon/pressure/k control changed');flat[model]['preserved_k_epsilon_pressure_bits']=True
    atomic(a.output/'flat.json',dict(models=flat,adopted=False,physical_steps=0,weak_continuity_enforced_in_KKT=True))
    atomic(a.output/'construction.json',dict(frozen=f,primary_patch=load(r/'prepare/level1/construction.json'),complete_response_bits_preserved=True,physical_identities_preserved=True))
    times={}
    for p in r.glob('*usage.txt'):
        try:times[p.name]=dict(zip(('wall_seconds','user_seconds','system_seconds','peak_RSS_KiB'),map(float,p.read_text().split())))
        except ValueError:pass
    atomic(a.output/'cost.json',dict(measurements=times,construction_numerical_seconds=f['numerical_seconds'],context_setup_seconds=f['total_setup_seconds'],optimizer=search['budget'],independent_replay_scans=count,total_search_and_replay_scans=search['budget']['complete_scans']+count,independent_replay_seconds=time.monotonic()-start,physical_steps=0))

if __name__=='__main__':main()
