#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Unused flat stress checks and physical normal-height extraction audit."""
import argparse,csv,json,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted4_reference import normal_distance,write_json
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('references','data','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();m=RequestManifest(a.output,dict(schema='CT4UnusedObservationsV1',reference_sha256=sha(a.references),data_sha256=sha(a.data/'snapshot.json'),source_sha256=sha(Path(__file__))))
    audit=subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('cao_tamura_data.py')),'--data',str(a.data),'--report',str(a.output/'data-audit.json')],capture_output=True,text=True)
    require(audit.returncode==0,'frozen_dataset_auditor_failed')
    refs=json.loads(a.references.read_text());rows=list(csv.DictReader((a.data/'profiles.csv').open()));checks={}
    for model in refs['fits']:
        c=refs['cases'][model+'-6.0'];rms=c['outer_equilibrium']['RMS_over_utau'];cases={}
        for series,predicted in [('f3b_R',.051**2),('f3c_R',rms[0]),('f3d_R',rms[2])]:
            data=[r for r in rows if r['series_id']==series and r['quality']=='visible'];weights=np.array([float(r['weight']) for r in data]);d=np.array([predicted-float(r['value']) for r in data])
            cases[series]=dict(predicted=predicted,weighted_RMS_discrepancy=float(np.sqrt(weights@(d*d)/weights.sum())),
                observed_range=[min(float(r['value']) for r in data),max(float(r['value']) for r in data)],
                points=[dict(id=r['point_id'],observed=float(r['value']),predicted=predicted,residual=float(v),extraction_bound=float(r['value_bound'])) for r,v in zip(data,d,strict=True)],
                role='unused CT4 check; original frozen calibration/cross-check role unchanged; no stress fitted')
        checks[model]=dict(outer_RMS_over_utau=dict(x=rms[0],spanwise=rms[1],vertical=rms[2]),
            inner_reference='all normal variances=2*k_w(z)/3; k_w is nonconstant, so inner RMS depends on normal distance; interface covariance need not be continuous',checks=cases)
    heights=[]
    for r in rows:
        if r['quality']!='visible' or not r['series_id'].startswith(('f4','f5','f8','f10','f12')) or not r['x_over_H']:continue
        x=float(r['x_over_H']);fixed=r['series_id'].startswith('f12')
        if fixed:z=.002/.04 if '_SS' in r['series_id'] else .007/.04;bound=None
        elif r['y_over_H']:z=float(r['y_over_H']);bound=float(r['axis_y_bound'])
        else:continue
        central=normal_distance(x,z)
        low=normal_distance(x,z-bound) if bound is not None else None;high=normal_distance(x,z+bound) if bound is not None else None
        heights.append(dict(series=r['series_id'],point_id=r['point_id'],x_over_Hh=x,global_z_over_Hh=z,
            nominal=central,vertical_extraction_bound_over_Hh=bound,lower_vertical_envelope=low,upper_vertical_envelope=high,
            fixed_physical_probe_height=fixed,position_uncertainty='not reported; plotted velocity bound is not a probe-height bound',
            interface_nominal={str(mm):('unavailable_or_inside' if central['distance'] is None else 'inner' if central['distance']<mm/40 else 'outer') for mm in (5.5,6.,6.5)}))
    summary={}
    for mm in (5.5,6.,6.5):
        samples=[q for q in heights if q['series'].startswith(('f5','f10','f12'))]
        summary[str(mm)]=dict(total=len(samples),inner=sum(q['interface_nominal'][str(mm)]=='inner' for q in samples),
            unavailable_or_inside=sum(q['interface_nominal'][str(mm)]=='unavailable_or_inside' for q in samples),
            scope='nominal closest-normal observation mapping only; C1 foot limits and possible multiple normal roots retained; no native curved capability or roughness-sublayer certificate')
    duplicate={}
    for model,f in refs['fits'].items():
        points=[];num=den=0.
        for r in rows:
            if r['series_id']!='f3e_R' or r['quality']!='visible':continue
            # Published axis is (y-d) in mm; the frozen physical-height column
            # already adds the authors' assumed d. Never displace the log law.
            z=float(r['y_over_H']);require(z>6.5/40,'duplicate_point_below_all_ct4_outer_interfaces')
            from ct_accepted4_reference import parameters
            pred=f['B']/f['ustar']+np.log(z)/parameters(model)['kappa'];delta=pred-float(r['value']);w=float(r['weight']);num+=w*delta*delta;den+=w
            points.append(dict(point_id=r['point_id'],physical_z_over_Hh=z,observed_U_over_utau=float(r['value']),predicted_U_over_utau=pred,
                residual=delta,extraction_bound=float(r['value_bound']),weight=w))
        require(len(points)==15,'missing_duplicate_flat_mean_coverage')
        duplicate[model]=dict(points=points,weighted_RMS_discrepancy=float(np.sqrt(num/den)),
            role='unchanged duplicate_crosscheck; not part of the CT4 four-point fit; original Fig3a/Fig3e normalization discrepancy remains visible')
    result=dict(schema='CT4UnusedObservationsV1',diagnostic_execution='pass',data_audit=json.loads((a.output/'data-audit.json').read_text()),
        flat_turbulence=checks,height_audit=heights,turbulence_height_summary=summary,
        duplicate_flat_mean=duplicate,
        flat_height_envelopes=refs['observation_heights']['flat'],Fig11='no physical height invented from unspecified delta75',
        normalization='reported utau/Uinf held fixed; no fitting of K or stresses; dimensional Fig12 RMS not scored without independent velocity scale',
        uncertainty='coordinate and velocity extraction are deterministic envelopes; shared-source dependence in CT4 fit/40mm normalization retained; probe, measurement, displacement and physical parameter uncertainty remain separately unquantified')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json','data-audit.json']);print('CT4 unused observation audit pass');return 0
if __name__=='__main__':raise SystemExit(main())
