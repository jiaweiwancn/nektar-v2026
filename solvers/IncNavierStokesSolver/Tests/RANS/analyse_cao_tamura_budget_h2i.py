#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent saved-cut dual replay and compact H2i diagnostic reduction.

No LP is solved here. A nonpositive or inconsistent replay never becomes an
infeasibility claim. The finite-box correction permits imperfect stationarity
and omitted tiny multipliers while retaining a rigorous necessary relaxation.
"""
import argparse, collections, os, time
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_h2i_convex import Context
from cao_tamura_h2h_scalar import Budget
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, THREADS, state_hash


def replay(P, phase, initial):
    result=phase['result']; proof=result['proof']; delta=result['delta']
    A=np.array(proof['rows']); b=np.array(proof['rhs']); y=np.array(proof['multipliers'])
    require(np.all(y<=0), 'dual multiplier sign')
    outer=np.flatnonzero(~P.inner); nout=len(outer)
    error=0.; witness=[]
    for row,rhs,m in zip(A,b,proof['cut_metadata']):
        kind=m['kind']
        if kind=='epsilon_bound':
            j=m['hard_row']; i=int(outer[j%nout]); sign=1 if j<nout else -1
            r=np.r_[np.zeros(P.nv),sign*P.E[i]/P.ec[i],0.]
            h=(result['epsilon_band'][1]-P.e[i]/P.ec[i]) if sign==1 else (P.e[i]/P.ec[i]-result['epsilon_band'][0])
            scale=max(1.,float(np.max(abs(r))),abs(h))
        elif kind=='PSD':
            r,h=P.rayleigh_cut(m['point'],m['vector'],delta); scale=m['scale']
        elif kind=='response_cap':
            r,h=P.cap_cut(m['point'],m['component'],m['sign']); scale=m['scale']
        else:
            Q=P.mass if kind=='mass' else P.div[{16:0,20:1,24:2}[int(kind.split('GL')[1])]]
            limit=.12**2 if kind=='mass' else P.div_limits[{16:0,20:1,24:2}[int(kind.split('GL')[1])]]
            r,h=P.quadratic_cut(Q,limit,np.r_[m['tangent'],np.zeros(P.ne)]); scale=m['scale']
        error=max(error,float(np.max(abs(row-r/scale))),abs(rhs-h/scale))
        witness.append(dict(m,dual_multiplier=float(y[len(witness)])))
    require(error<2e-13, 'saved necessary cuts differ from immutable data')
    c=np.r_[np.zeros(P.dim),1.]
    bounds=np.array(P.bounds+[(0.,initial['slack']+1.)],np.longdouble)
    # inf over a containing finite box is <= inf over all physical constraints.
    # No numerical stationarity equality or missing bound multiplier is assumed.
    Al=A.astype(np.longdouble); bl=b.astype(np.longdouble); yl=y.astype(np.longdouble)
    d=c.astype(np.longdouble)-Al.T@yl
    lower=bl@yl+np.sum(d*np.where(d>=0,bounds[:,0],bounds[:,1]))
    box=np.max(abs(bounds),axis=1)
    arithmetic=32*np.finfo(float).eps*(1+np.sum(abs(bl*yl))+np.sum(abs(Al)*abs(yl[:,None])*box))
    allowance=result['proof']['receipt']['certificate']['native_data_rounding_allowance']
    checked=float(lower-arithmetic-allowance)
    primal=result['proof']['receipt']['certificate']['primal']
    require(checked>1e-7 and abs(primal-float(lower))<1e-7, 'independent limited dual replay failed')
    return dict(anchor=phase['anchor'],buffered=phase['buffered'],cut_reconstruction_error=error,
                lower_bound_before_allowances=float(lower),arithmetic_allowance=float(arithmetic),
                native_data_allowance=allowance,checked_lower_bound=checked,
                primal=primal,gap_before_allowances=primal-float(lower),
                stationarity_inf=float(np.max(abs(d))),dual_signs=True,
                active_kinds=dict(collections.Counter(m['kind'] for m in proof['cut_metadata'])),
                limiting_constraints=witness,verified=True,
                scope='finite sampled transformed-tensor slack; fixed-k anchor, 41 velocity/8 epsilon directions and original physical caps; not global physical infeasibility')


def tensor_check(w, nt=None):
    o=w['observation']; k=o['k']; e=o['epsilon']; G=np.array(o['gradient']); R=np.array(w['tensor'])
    if nt is None:
        require(not o['inner'], 'inner replay requires independently saved auxiliary viscosity')
        nt=.09*k*k/e
    iso=2*k/3*np.eye(3); linear=-nt*(G+G.T); independent=iso+linear
    ev,V=np.linalg.eigh(independent); v=np.array(w['critical_eigenvector'])
    return dict(tensor_error=float(np.max(abs(independent-R))),eigenvalue_error=abs(float(ev[0])-w['eigenvalues'][0]),
                isotropic=iso.tolist(),linear=linear.tolist(),nonlinear_zero=True,regularization_zero=True,
                eigenvalues=ev.tolist(),Rayleigh=float(v@independent@v),eigenvector_residual=float(np.linalg.norm(independent@v-ev[0]*v)),
                production=float(np.sum((nt*(G+G.T))*G)),eddy_viscosity=nt)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
    start=time.monotonic();os.environ.update(THREADS);r=a.root;previous=load(r/'independent-reduction.json') if (r/'independent-reduction.json').exists() else {};ctx=Context(r/'joint-context/context.json');search=load(r/'joint-search/result.json')
    theta=search['retention']['theta'];out=dict(version='H2iIndependentReductionV1',replays=[],native={},flat={},full_scans=0)
    matrices={}
    for anchor in ('A','B'):
        P=ctx.problem(anchor,theta); initial={}
        for ph in (x for x in search['phases'] if x['anchor']==anchor):
            initial=P.scan(np.zeros(P.dim),ph['result']['delta'],Budget(seconds=60));out['full_scans']+=1
            out['replays'].append(replay(P,ph,initial))
        candidate=load(r/'joint-search'/(anchor+'-actual')/'best-candidate.json')
        receipt=load(r/('native-'+anchor)/'receipt.json');d=load(Path(receipt['output']))
        s=P.scan(np.array(candidate['joint_coefficients']),1e-4*P.K0,Budget(seconds=60));out['full_scans']+=1
        require(receipt['state_sha256']==candidate['state_sha256'], 'native candidate binding')
        seen={}
        for order,v in d['seen_sampling'].items():
            seen[order]={k:x for k,x in v.items() if k not in ('state','samples')}
            for key in ('inner_covariance_witness','outer_covariance_witness'):
                w=v[key];o=w['observation'];points=ctx.arrays['points'];mask=np.all(ctx.block('chart')==o['coordinate'],axis=1)&(P.geometry==o['geometry'])
                ids=np.flatnonzero(mask);require(len(ids)==1, 'independent witness location is ambiguous')
                nt=P.nt[ids[0]] if o['inner'] else None
                seen[order][key+'_independent']=tensor_check(w,nt)
                require(seen[order][key+'_independent']['tensor_error']<1e-12,'independent native tensor replay')
        au=d['audit'];out['native'][anchor]=dict(candidate_sha256=candidate['state_sha256'],anchor_k_sha256=candidate['anchor_k_sha256'],
            passed=receipt['passed'],gates=receipt['gates'],native_covariance=au['accepted_awf_covariance'],
            scalar_budget=au['accepted_awf_scalar_domain'],residual_and_constraints=au['coupled_norms'],
            native_inner_witness=au['accepted_awf_inner_minimum_point'],native_outer_witness=au['accepted_awf_outer_minimum_point'],
            native_inner_provenance=au['native_inner_minimum_geometry_element_side_block'],native_outer_provenance=au['native_outer_minimum_geometry_element_side_block'],
            seen=seen,seconds=receipt['seconds'],optimizer_true_scan={k:s[k].tolist() if isinstance(s[k],np.ndarray) else s[k] for k in ('inner_min','outer_min','mass_squared','strong_divergence_squared','velocity_max','epsilon_ratio','hard_violation')})
        # A bounded native accessibility packet: extrema plus important sided
        # dual rows, all recomputed at this actual diagnostic state.
        selected=[]
        for mask in (P.cov&P.inner,P.cov&~P.inner):
            ids=np.flatnonzero(mask);selected.append(int(ids[np.argmin(s['lambda_covariance'][ids])]))
        sided=[]
        chart=ctx.block('chart'); foot=np.all(chart==[.23549580558569316,0,.0065],axis=1)
        for geom in (35,41):
            ids=np.flatnonzero(foot&(P.geometry==geom)&P.cov)
            for i in ids:
                G=np.array([[s['gradient'][i,0],0,s['gradient'][i,1]],[0,0,0],[s['gradient'][i,2],0,s['gradient'][i,3]]]);nu=P.nt[i] if P.inner[i] else .09*P.k[i]**2/s['epsilon'][i];R=2*P.k[i]/3*np.eye(3)-nu*(G+G.T);ev,V=np.linalg.eigh(R)
                sided.append(dict(point=int(i),chart=chart[i].tolist(),physical=ctx.block('physical')[i].tolist(),geometry=geom,inner=bool(P.inner[i]),dependent_epsilon=bool(ctx.block('dependent')[i]),gradient=G.tolist(),k=float(P.k[i]),epsilon=float(s['epsilon'][i]),nu=float(nu),tensor=R.tolist(),eigenvalues=ev.tolist(),eigenvector=V[:,0].tolist(),memberships=ctx.arrays['members'][ctx.arrays['members'][:,1]==i].tolist()))
            actual=[int(i) for i in ids if not ctx.block('dependent')[i]]
            if actual:selected.append(actual[0])
        out['native'][anchor]['old_protected_matching_sides']=sided
        proof=next(x['result']['proof'] for x in search['phases'] if x['anchor']==anchor and not x['buffered'])
        ranked=sorted(zip(proof['cut_metadata'],proof['multipliers']),key=lambda v:abs(v[1]),reverse=True)
        for m,y in ranked:
            if m['kind']=='PSD' and m['point'] not in selected:selected.append(m['point'])
            if len(selected)==6:break
        witnesses=[]
        for i in selected:
            G=np.array([[s['gradient'][i,0],0,s['gradient'][i,1]],[0,0,0],[s['gradient'][i,2],0,s['gradient'][i,3]]]);nu=P.nt[i] if P.inner[i] else .09*P.k[i]**2/s['epsilon'][i];ev,V=np.linalg.eigh(2*P.k[i]/3*np.eye(3)-nu*(G+G.T))
            witnesses.append(dict(chart=ctx.block('chart')[i].tolist(),geometry=int(P.geometry[i]),eigenvector=V[:,0].tolist(),point=int(i),covariance=float(ev[0]),memberships=ctx.arrays['members'][ctx.arrays['members'][:,1]==i].tolist()))
        packet=dict(operation='h2i_accessibility',expected_model='StandardKEpsilon',inner_stress_model='Linear',state=candidate['state'],state_sha256=candidate['state_sha256'],base_sha256=candidate['base_sha256'],basis_file=ctx.doc['basis_file'],basis_sha256=ctx.doc['basis_sha256'],witnesses=witnesses)
        atomic(r/('accessibility-'+anchor+'-input.json'),packet)
    flat=load(r/'flat-H1/results.json')
    for model,item in flat.items():
        d=item['data'];out['flat'][model]={k:d[k] for k in ('base_sha256','state_sha256','H1_length','setup','scalar_domain','fields_changed')}
        for key in ('baseline','projected'):out['flat'][model][key]={k:v for k,v in d[key].items() if k not in ('state','reference','defect')}
    total_scans=previous.get('cumulative_replay_scans',previous.get('full_scans',0))+out['full_scans']
    replay_seconds=previous.get('cumulative_replay_seconds',previous.get('seconds',0))+time.monotonic()-start
    require(search['budget']['complete_scans']+total_scans<=120 and search['budget']['seconds']+replay_seconds<=600,'aggregate diagnostic/optimizer budget exhausted')
    out.update(seconds=time.monotonic()-start,cumulative_replay_scans=total_scans,cumulative_replay_seconds=replay_seconds,optimizer_aggregate=dict(search['budget'],diagnostic_replay_scans=total_scans),
               context_sha256=sha(ctx.path),search_sha256=sha(r/'joint-search/result.json'),physical_attempts=0,fresh_GL28=0)
    atomic(r/'independent-reduction.json',out)

if __name__=='__main__':main()
