#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT2 two-axis deterministic extraction envelopes and fixed-fit wall heights.
No statistical confidence interval, experimental error model or native terrain.
Invoke through the owned CT launcher.
"""
import argparse,csv,itertools,json,math,os,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import Column,parameters,fit_roughness,native_compare,write_json
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha,runtime_identity

def interpolate(z,v,target,logarithmic):
    x=np.log(z-3.5) if logarithmic else z;t=math.log(target-3.5) if logarithmic else target
    j=max(0,min(len(z)-2,int(np.searchsorted(z,target))-1));a=(t-x[j])/(x[j+1]-x[j])
    return float((1-a)*v[j]+a*v[j+1]),[j,j+1]

def envelope(z,v,dz,dv,target,logarithmic,height=True,velocity=True):
    # All adjacent bracketing pairs and polytope vertices, including z_i=target.
    # The data's disjoint height intervals guarantee preserved source ordering.
    rows=[]
    for i in range(len(z)-1):
        js=(i,i+1);choices=[]
        for j in js:
            zz=[z[j]-dz[j],z[j]+dz[j]] if height else [z[j]]
            if height and z[j]-dz[j]<=target<=z[j]+dz[j]:zz.append(target)
            vv=[v[j]-dv[j],v[j]+dv[j]] if velocity else [v[j]]
            choices.append(list(itertools.product(zz,vv)))
        for a,b in itertools.product(*choices):
            if not a[0]<=target<=b[0] or a[0]>=b[0]:continue
            zz=z.copy();vv=v.copy();zz[i],vv[i]=a;zz[i+1],vv[i+1]=b
            val,ids=interpolate(zz,vv,target,logarithmic)
            rows.append(dict(value=val,bracket=ids,coordinates_mm=zz.tolist(),values=vv.tolist()))
    require(rows,'no_possible_extraction_bracket')
    return dict(lower=min(rows,key=lambda r:r['value']),upper=max(rows,key=lambda r:r['value']),vertices=len(rows),possible_brackets=sorted(set(tuple(r['bracket']) for r in rows)))

def scenarios(data):
    rows=[r for r in csv.DictReader((data/'profiles.csv').open()) if r['series_id']=='f3a_R' and float(r['weight'])>0]
    rows.sort(key=lambda r:float(r['axis_y']));z=np.array([float(r['axis_y']) for r in rows]);v=np.array([float(r['value']) for r in rows]);dz=np.array([float(r['axis_y_bound']) for r in rows]);dv=np.array([float(r['value_bound']) for r in rows])
    require(np.all(z[1:]-dz[1:]>z[:-1]+dz[:-1]),'coordinate_intervals_overlap_requires_order_polytope_extension')
    cases=[dict(name='central',coordinates_mm=z.tolist(),values=v.tolist())];bounds={}
    for label,h,b in [('height_only',True,False),('velocity_only',False,True),('both_axes',True,True)]:
        a=envelope(z,v,dz,dv,16.,True,h,b);c=envelope(z,v,dz,dv,40.,False,h,b);bounds[label]=dict(calibration16=a,normalization40=c)
        # The two brackets use disjoint point sets here; each joint corner is
        # attainable, not an independence assumption for random errors.
        used16=set(sum([list(x) for x in a['possible_brackets']],[]));used40=set(sum([list(x) for x in c['possible_brackets']],[]))
        require(not used16&used40,'shared_anchor_points_require_joint_polytope')
        for s,t in itertools.product(('lower','upper'),repeat=2):
            zz=z.copy();vv=v.copy()
            for j in used16:zz[j]=a[s]['coordinates_mm'][j];vv[j]=a[s]['values'][j]
            for j in used40:zz[j]=c[t]['coordinates_mm'][j];vv[j]=c[t]['values'][j]
            cases.append(dict(name=label+'-16'+s+'-40'+t,coordinates_mm=zz.tolist(),values=vv.tolist()))
    for sz,sv in itertools.product((-1,1),repeat=2):cases.append(dict(name='coherent_axes-'+str(sz)+'-'+str(sv),coordinates_mm=(z+sz*dz).tolist(),values=(v+sv*dv).tolist()))
    for c in cases:
        zz=np.array(c['coordinates_mm']);vv=np.array(c['values']);a,ia=interpolate(zz,vv,16.,True);b,ib=interpolate(zz,vv,40.,False)
        c.update(U16_over_Uinf=a,U40_over_Uinf=b,U16_over_UH=a/b,Uinf_over_UH=1/b,utau_over_UH=.051/b,bracket16=[rows[i]['point_id'] for i in ia],bracket40=[rows[i]['point_id'] for i in ib])
    return cases,bounds

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--data',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    prior_contract=json.loads((a.previous/'source/docs/plans/cao-tamura-accepted-1-evidence/contract.json').read_text())
    for name,digest in prior_contract['frozen_dataset']['hashes'].items():require(sha(a.data/name)==digest,'frozen_dataset_changed:'+name)
    req=dict(schema='CT2WallHeightExtractionV1',data_snapshot_sha256=sha(a.data/'snapshot.json'),uncertainty_sha256=sha(a.data/'uncertainty.json'),previous_sha256=sha(a.previous/'feasibility-v2/results.json'),code_sha256=sha(Path(__file__)),native_probe=runtime_identity(a.runtime/'RANSNewModelsProbe'),heights_mm=[5,8,11],scenario_count=17,quadrature=[48,96],thin_thresholds=[.1,.25],scope='independent frozen reference columns, no curved velocity/covariance certificate')
    manifest=RequestManifest(a.output,req);cases,bounds=scenarios(a.data)
    result=dict(schema='CT2WallOperatingHeightV1',diagnostic_execution='pass',status='pass',extraction_bounds=bounds,scenarios=cases,flat=[],curved=[],fits={},historical_smooth='CT1 8mm negative-k results unchanged; SS preparation unresolved',uncertainty='deterministic extraction envelope only. Assumed d=3.5mm fixed; delta/Re/utau physical and measurement uncertainties unquantified. Coherent shifts account for a possible shared axis placement; conservative per-point bounds do not assert statistical independence.')
    process=subprocess.Popen([str(a.runtime/'RANSNewModelsProbe')],stdin=subprocess.PIPE,stdout=subprocess.PIPE,text=True)
    stations=sorted(set(np.linspace(-2.5,2.5,33).tolist()+[-3.,-2.5-1e-8,-2.5+1e-8,-1.25,0.,1.25,2.5-1e-8,2.5+1e-8,3.]))
    previous=json.loads((a.previous/'feasibility-v2/results.json').read_text());nu=1/10500
    # CT1 used the stored y/H column (rounded separately from axis_y). Keep
    # its central normalization/calibration EXACT; retain the tiny column
    # representation difference explicitly rather than retuning central hr.
    central=cases[0];ref=previous['incident']['R']
    central['axis_column_rounding_comparison']={k:central[k] for k in ('U16_over_Uinf','U40_over_Uinf','U16_over_UH','utau_over_UH')}
    central.update(U16_over_Uinf=ref['anchors']['0.016']['value'],U40_over_Uinf=ref['U_H_over_Uinf']['value'],
        U16_over_UH=ref['anchors']['0.016']['value']*ref['Uinf_over_U_H'],Uinf_over_UH=ref['Uinf_over_U_H'],utau_over_UH=ref['ustar'])
    central['central_input_policy']='sealed CT1 exact values; other deterministic scenarios use frozen axis columns and both transformed extraction bounds'
    try:
        for case in cases:
            for model in ('A93','A95','E94E95'):
                key=case['name']+':'+model;us=case['utau_over_UH'];K=parameters(model)['q']*us*us
                fit=fit_roughness(model,nu,us,.4,case['U16_over_UH'],'velocity');result['fits'][key]=fit
                if fit['status']!='fit':continue
                hr=fit['hr']
                if case['name']=='central':require(abs(hr-previous['fits'][model+':model-velocity-fit']['hr'])<2e-12,'central_roughness_changed')
                for mm in (5,8,11):
                    H=mm/40;c=Column(model,K,H,hr,nu);row=c.result();tight=Column(model,K,H,hr,nu,order=96).result()
                    row.update(scenario=case['name'],height_mm=mm,hr_mm=hr*40,hr_lt_H=hr<H,interface_over_cube_height=mm/5,quadrature_uncertainty=dict(cutoff_over_H=abs(c.cutoff-tight['cutoff'])/H if c.cutoff else 0.,minimum_k_over_K=abs(row['reference_k']['minimum_over_K']-tight['reference_k']['minimum_over_K']) if c.cutoff else 0.))
                    require(max(row['quadrature_uncertainty'].values())<2e-8,'flat_quad_unresolved')
                    if case['name']=='central':row['native']=native_compare(process,c,row)
                    result['flat'].append(row)
                    for x in stations:
                        slope=-math.pi/5*math.sin(math.pi*x/2.5) if abs(x)<2.5 else 0.;second=-math.pi**2/12.5*math.cos(math.pi*x/2.5) if abs(x)<2.5 else 0.;curv=second/(1+slope*slope)**1.5
                        c=Column(model,K,H,hr,nu,curv);r=c.result();ratio=abs(curv)*H
                        r.update(scenario=case['name'],height_mm=mm,x_over_Hh=x,hr_mm=hr*40,hr_lt_H=hr<H,thin_ratio=ratio,thin_0p1=ratio<=.1,thin_0p25=ratio<=.25,global_reach=12.5/math.pi**2,global_reach_valid=H<12.5/math.pi**2,curved_covariance='not certified without actual Cartesian velocity/strain')
                        if x==0. or x in (-1.25,1.25) or abs(abs(x)-2.5)<2e-8:
                            t=Column(model,K,H,hr,nu,curv,96).result()
                            r['quadrature_uncertainty']=dict(cutoff_over_H=abs(c.cutoff-t['cutoff'])/H if c.cutoff else 0.,minimum_k_over_K=abs(r['reference_k']['minimum_over_K']-t['reference_k']['minimum_over_K']) if c.cutoff else 0.)
                            require(max(r['quadrature_uncertainty'].values())<2e-8,'metric_quad_unresolved')
                        result['curved'].append(r)
                write_json(a.output/'progress.json',dict(fits=len(result['fits']),flat=len(result['flat']),curved=len(result['curved'])))
        summary={}
        for model in ('A93','A95','E94E95'):
            for mm in (5,8,11):
                f=[r for r in result['flat'] if r['model']==model and r['height_mm']==mm];c=[r for r in result['curved'] if r['model']==model and r['height_mm']==mm]
                summary[model+':'+str(mm)]=dict(hr_mm_range=[min(r['hr_mm'] for r in f),max(r['hr_mm'] for r in f)],flat_admissible_scenarios=sum(r['classification']=='admissible_reference' for r in f),scenarios=len(f),minimum_column_k_over_K=min(r['reference_k']['minimum_over_K'] for r in c if r['cutoff'] is not None),maximum_thin_ratio=max(r['thin_ratio'] for r in c),thin_0p1=all(r['thin_0p1'] for r in c),thin_0p25=all(r['thin_0p25'] for r in c),all_global_reach=all(r['global_reach_valid'] for r in c),all_stretches_valid=all(r['minimum_stretch']>=.2 for r in c),all_reference_columns_admissible=all(r['classification']=='admissible_reference' for r in c),interfaces='5mm at nominal cube tops; 8/11mm above cubes but roughness-sublayer homogenization not demonstrated; wall origin unchanged')
        result['summary']=summary
    finally:process.stdin.close();process.wait()
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']);print('pass',len(result['flat']),len(result['curved']));return 0
if __name__=='__main__':raise SystemExit(main())
