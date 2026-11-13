#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finish declared H11 local component comparisons from verified scalar records."""
import argparse
import csv
from pathlib import Path
import math
import numpy as np
from run_wj_h11 import load,save
from publish_wj_h11 import table
from study_wj_h9 import identity
from study_wj_h6 import configuration
from wj_h11_analysis import ratio


def run(root):
    out=root/'local-comparisons-v1';out.mkdir();pub=root/'publication-v1';c=configuration(6.)
    scales=dict(U=1.,W=1.,G=1.,k=c.K,epsilon=c.A/c.H,R=c.K,pi=c.u_tau**2,mechanical_pressure=c.u_tau**2,margin=1.,transverse=1.)
    source=root/'science/common-v1/fixed-owner-history.json';fixed=load(source);groups={}
    for v in fixed:
        d=v['record'];p=d['physical'];o=d['owner'];step=round(v['time']/.01)
        if abs(v['time']-step*.01)>128*np.finfo(float).eps:continue
        key=(step,tuple(float(x).hex()for x in p),o['geometry'],o['normal_side'],o['x_side']);groups.setdefault(key,{})[v['branch']]=d
    rows=[]
    for key,values in groups.items():
        if set(values)!=set(('coarse','fine','quarter')):continue
        data=[values[b]for b in ('coarse','fine','quarter')];p=data[0]['physical'];o=data[0]['owner']
        for field in ('U','G','k','epsilon','R','pi','mechanical_pressure','margin','transverse'):
            arrays=[np.asarray(d[field])for d in data]
            for component in np.ndindex(arrays[0].shape)if arrays[0].shape else [()]:
                a,b,d=[float(v[component])for v in arrays];label=field+''.join(str(x)for x in component)
                scale=scales.get(field,1.);rows.append(dict(time=key[0]*.01,x=p[0],physical_z=p[2],physical_z_mm=40*p[2],geometry=o['geometry'],normal_side=o['normal_side'],x_side=o['x_side'],quantity=label,values=[a,b,d],signed_D1=a-b,signed_D2=b-d,scale=scale,**ratio(abs(a-b),abs(b-d),4e-11*scale)))
        projectors=[]
        for d in data:
            xi=np.array([math.cos(d['angle']),0,math.sin(d['angle'])]);projectors.append(np.outer(xi,xi))
        rows.append(dict(time=key[0]*.01,x=p[0],physical_z=p[2],physical_z_mm=40*p[2],geometry=o['geometry'],normal_side=o['normal_side'],quantity='minimizing_direction_projector',scale=1.,**ratio(float(np.linalg.norm(projectors[0]-projectors[1])),float(np.linalg.norm(projectors[1]-projectors[2])),4e-11)))
    table(out/'local-component-pairs.csv',rows)
    pairs=[]
    with(pub/'field-pairs.csv').open(newline='')as f:
        for row in csv.DictReader(f):
            denominator=float(row['scale'])*(math.sqrt(float(row['measure']))if row['norm']=='L2'else 1.)
            pairs.append(dict(row,normalized_D1=float(row['D1'])/denominator,normalized_D2=float(row['D2'])/denominator,normalized_floor=float(row['floor'])/denominator))
    table(out/'field-pairs-with-scales.csv',pairs)
    normalized=[]
    measures=load(pub/'Q1-continuity.json')[0]['regions'];volume={'inner':measures['4']['measure'],'outer':measures['5']['measure']};volume['full']=volume['inner']+volume['outer']
    for record in load(pub/'errors-drift.json'):
        result={k:record[k]for k in ('branch','step','time')}
        for kind in ('reference_error','compatible_CP0_drift'):
            result[kind]={}
            for scope,fields in record[kind].items():
                result[kind][scope]={}
                for field,values in fields.items():
                    measure=volume['outer']if field=='epsilon'else volume[scope];scale=scales[field]
                    result[kind][scope][field]=dict(RMS_over_amplitude=values['L2']/(scale*math.sqrt(measure)),sampled_max_over_amplitude=values['sampled_max']/scale,amplitude=scale,measure=measure)
        normalized.append(result)
    save(out/'normalized-errors-drift.json',normalized)
    save(out/'receipt.json',dict(status='PASS',source_fixed=identity(source),source_pairs=identity(pub/'field-pairs.csv'),source_errors=identity(pub/'errors-drift.json'),
        source=identity(__file__),scales=scales,physical_measures=volume,local_rows=len(rows),field_rows=len(pairs),scope='only declared same-owner/field comparisons and normalization of existing scalar bytes; no state evaluation/native work; direction uses sign-invariant xi tensor xi'))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()]))
    print('H11 local/scaled scalar tables PASS',len(rows),len(pairs),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();run(a.root)
