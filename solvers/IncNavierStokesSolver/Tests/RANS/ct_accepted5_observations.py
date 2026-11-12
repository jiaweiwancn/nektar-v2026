#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned CT4 observation correction without refitting or native evolution."""
import argparse,ast,csv,hashlib,json,math,os,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from scipy.optimize import brentq,minimize_scalar
from ct_accepted5_geometry import closest_normal,hill,slope,A
from ct_accepted4_reference import write_json
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha

def legacy():
    g=['git','-c','index.threads=1','-c','core.preloadIndex=false','show','b9882937a5096aa9c73219dd3b6288cd8a346983:solvers/IncNavierStokesSolver/Tests/RANS/ct_accepted4_reference.py']
    source=subprocess.check_output(g,cwd=Path(__file__).resolve().parents[4],text=True)
    tree=ast.parse(source);node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='normal_distance')
    ns=dict(np=np,math=math,brentq=brentq);exec(compile(ast.Module(body=[node],type_ignores=[]),'reviewed_CT4_function','exec'),ns)
    return ns['normal_distance'],hashlib.sha256(ast.get_source_segment(source,node).encode()).hexdigest()

def controls():
    old,digest=legacy();w=old(0.,1.25);v=closest_normal(0.,1.25)
    require(abs(w['distance']-2.7950849718747373)<1e-12 and w['footpoint_x']==-2.5,'legacy_crest_witness_not_reproduced')
    require(v['footpoint_x']==0 and v['distance']==.25,'corrected_crest_witness_failed')
    cases=[]
    def check(x,z,r=None,t=None):
        q=closest_normal(x,z);require(q['distance'] is not None,'unexpected_inside_control')
        require(q['distance']<=z-hill(x)+2e-12,'vertical_candidate_bound')
        mirror=closest_normal(-x,z);require(abs(q['distance']-mirror['distance'])<2e-12 and abs(q['footpoint_x']+mirror['footpoint_x'])<2e-9,'left_right_symmetry')
        if r is not None:require(abs(q['distance']-r)<2e-11,'independent_normal_distance')
        if t is not None:require(abs(q['footpoint_x']-t)<2e-9,'independent_footpoint')
        cases.append(dict(x=x,z=z,expected_distance=r,expected_foot=t,result=q));return q
    for z in (1.,1+1e-9,1.01,1.25,2.,8.):check(0.,z,z-1,0.)
    for x in (-5.,-3.,3.,5.):
        for z in (0.,.05,.175,.4):check(x,z,z,x)
    for t in (-2.5,-2.499999999,-2.4,-1.25,-.3,.3,1.25,2.4,2.499999999,2.5):
        for r in (0.,.01,.15):
            hp=slope(t);n=math.hypot(1,hp);x=t-r*hp/n;z=hill(t)+r/n
            # Zero-offset values computed independently can differ by one ulp
            # from hill(x); retain the exact graph for on-surface controls.
            if r==0:z=hill(x)
            check(x,z,r,t)
    continuity=[check(x,1.25) for x in (0.,1e-10,1e-8,1e-6,1e-4)]
    require(abs(continuity[1]['distance']-.25)<1e-12,'crest_continuity')
    # Analytic F=F'=0, a double stationary root that need not change sign.
    t=2.;hpp=-2*A*A*math.cos(2*A*t);z=hill(t)+(1+slope(t)**2)/hpp;x=t-slope(t)*(z-hill(t));q=check(x,z)
    require(any(abs(t-u)<2e-9 for u in q['normal_roots']),'tangent_root_omitted')
    # Independent distance minimization, not the stationarity polynomial.
    oracle=[]
    for x,z in [(.001,1.25),(.73,1.5),(2.3,1.7),(2.6,2.),(x,z),(-2.5,2.),(0.,8.)]:
        knots=np.linspace(-2.5,2.5,65);candidates=[(z*z,x)] if abs(x)>=2.5 else []
        fun=lambda t:(t-x)**2+(hill(t)-z)**2
        for lo,hi in zip(knots[:-1],knots[1:]):
            m=minimize_scalar(fun,bounds=(lo,hi),method='bounded',options={'xatol':1e-14});candidates.extend([(m.fun,m.x),(fun(lo),lo),(fun(hi),hi)])
        d2,foot=min(candidates);native=closest_normal(x,z);require(abs(math.sqrt(d2)-native['distance'])<2e-10,'independent_global_distance_minimum')
        oracle.append(dict(x=x,z=z,independent_distance=math.sqrt(d2),root_method_distance=native['distance'],independent_foot=foot))
    return dict(status='pass',reviewed_function_sha256=digest,old_witness=w,corrected_witness=v,analytic_controls=cases,independent_minimization=oracle,
        completeness='all stationary roots isolated by analytic critical points; tangent and endpoint roots included; binary64 near-coalescent multiplicities explicitly delimited in method receipts')

def classify(q,mm):return 'unavailable_or_inside' if q['distance'] is None else 'inner' if q['distance']<mm/40 else 'outer'

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('data','ct4-audit','references','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();m=RequestManifest(a.output,dict(schema='CT5ObservationAuditV1',sources={f:sha(Path(__file__).with_name(f)) for f in ('ct_accepted5_observations.py','ct_accepted5_geometry.py','ct_accepted4_reference.py')},
        original_audit_sha256=sha(a.ct4_audit),frozen_reference_sha256=sha(a.references),dataset_sha256=sha(a.data/'snapshot.json')))
    result=dict(schema='CT5ObservationAuditV1',status='fail',controls=controls(),height_audit=[],corrections=[],Fig7=[],new_native_evolution=False)
    old=json.loads(a.ct4_audit.read_text());rows=list(csv.DictReader((a.data/'profiles.csv').open()));series={q['id']:q for q in json.loads((a.data/'series.json').read_text())}
    for q in old['height_audit']:
        x=q['x_over_Hh'];z=q['global_z_over_Hh'];bound=q['vertical_extraction_bound_over_Hh'];new=dict(q,nominal=closest_normal(x,z))
        for key,zz in [('lower_vertical_envelope',z-bound if bound is not None else None),('upper_vertical_envelope',z+bound if bound is not None else None)]:new[key]=closest_normal(x,zz) if zz is not None else None
        new['interface_nominal']={str(mm):classify(new['nominal'],mm) for mm in (5.5,6.,6.5)};result['height_audit'].append(new)
        nd,od=new['nominal']['distance'],q['nominal']['distance']
        if (nd is None)!=(od is None) or (nd is not None and abs(nd-od)>2e-10):
            result['corrections'].append(dict(series=q['series'],point_id=q['point_id'],x=x,z=z,old=q['nominal'],new=new['nominal'],old_classification=q['interface_nominal'],new_classification=new['interface_nominal']))
    for r in rows:
        if not r['series_id'].startswith('f7') or r['quality']!='visible':continue
        sid=r['series_id'];z=1+6.25*float(r['axis_y']);bound=6.25*float(r['axis_y_bound'])
        q=dict(series=sid,point_id=r['point_id'],x_over_Hh=0.,global_z_over_Hh=z,vertical_extraction_bound_over_Hh=bound,
            nominal=closest_normal(0,z),lower_vertical_envelope=closest_normal(0,z-bound),upper_vertical_envelope=closest_normal(0,z+bound),
            role=r['role'],quantity=series[sid]['quantity'],value=float(r['value']),value_bound=float(r['value_bound']),weight=float(r['weight']),
            fixed_physical_probe_height=False,coordinate='global z/Hh=1+6.25*(yprime/delta); extraction bounds scaled by6.25',position_uncertainty='not reported; extraction bound is not experimental/probe uncertainty')
        q['interface_nominal']={str(mm):classify(q['nominal'],mm) for mm in (5.5,6.,6.5)};result['Fig7'].append(q);result['height_audit'].append(q)
    changed=result['corrections'];require(len(changed)==40 and all(q['x']==0 for q in changed),'unexpected_CT4_nominal_correction_coverage')
    flips={str(mm):[dict(series=q['series'],point_id=q['point_id']) for q in changed if '_SS' in q['series'] and q['old_classification'][str(mm)]=='outer' and q['new_classification'][str(mm)]=='inner'] for mm in (5.5,6.,6.5)}
    require(all(len(v)==4 for v in flips.values()) and len(result['Fig7'])==36,'expected_smooth_flips_or_Fig7_coverage_missing')
    lower={sid:min(q['lower_vertical_envelope']['distance']*40 for q in result['Fig7'] if q['series']==sid) for sid in {q['series'] for q in result['Fig7']}}
    require(9.33<lower['f7a_RR']<9.38 and 9.33<lower['f7b_RR']<9.38,'Fig7_rough_bound_conversion')
    rough=[q for q in result['height_audit'] if '_RR' in q['series'] and q['series'].startswith(('f5','f7','f10','f12'))]
    limiter=min(rough,key=lambda q:(q['lower_vertical_envelope'] or q['nominal'])['distance'] if (q['lower_vertical_envelope'] or q['nominal'])['distance'] is not None else math.inf)
    dist=(limiter['lower_vertical_envelope'] or limiter['nominal'])['distance']*40
    require(limiter['series'].startswith('f12_RR') and abs(dist-7)<1e-9,'rough_limiting_height_changed')
    refs=json.loads(a.references.read_text());ratios={key:dict(inner_auxiliary_nu_over_nut_H=c['nu_over_nut_H'],outer_nu_over_kappa_utau_H=c['nu']/(c['coefficients']['kappa']*c['utau']*c['H']),physical_outer_molecular_diffusion='zero by stipulated equation regime; this ratio is an omitted-scale indicator') for key,c in refs['cases'].items()}
    result.update(status='pass',corrected_CT4_nominal_distances=len(changed),smooth_outer_to_inner=flips,old_coverage=len(old['height_audit']),new_coverage=len(result['height_audit']),Fig7_rough_lower_heights_mm=lower,
        limiting_rough_turbulence_height_mm=dist,limiting_rough_point=dict(series=limiter['series'],point_id=limiter['point_id']),high_Re_indicators=ratios,
        Fig11='physical heights unavailable where delta7.5 unspecified',uncertainty='frozen coordinate/velocity extraction envelopes retained; nominal delta and assumed displacement separate; no invented experimental/probe uncertainty',
        physical_inputs='not refitted; CT4 fixed reference artifact unchanged',original_audit_scope='CT4 sealed record retained; closest-distance/classification claims corrected by this version; mathematical wall and native CT4 trajectory decisions unchanged')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT5 observation audit pass',len(changed),len(result['Fig7']));return 0
if __name__=='__main__':raise SystemExit(main())
