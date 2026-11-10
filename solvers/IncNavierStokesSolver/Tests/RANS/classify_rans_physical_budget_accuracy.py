#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Classify finite physical-budget accuracy separately from reconstruction."""
import argparse,json,math
from pathlib import Path
from rans_duration_contract import require,MEMBERS
from rans_request_manifest import RequestManifest,sha
from rans_qualification import write_json


def classify_row(row):
    values=row['signed_values'];require(set(values)==set(MEMBERS),'budget_companion_missing')
    require(all(math.isfinite(v) for v in values.values()),'nonfinite_physical_budget')
    c,f=map(abs,[values['B5'],values['B11']]);u=row['quadrature_uncertainty']
    require(math.isfinite(u) and u>=0,'invalid_budget_uncertainty')
    floor=max(5e-9,5*u);gap=abs(values['B5']-values['B11'])
    delta=max(abs(values[m]-values[m+'-dt2']) for m in MEMBERS[:2]);limit=max(1e-9,.1*gap)
    spatial='unresolved_floor' if max(c,f)<=floor else 'pass' if f<=.85*c and c-f>5*u else 'failed_or_unresolved'
    return dict(kind=row['kind'],equation=row['equation'],signed_values=values,
        spatial_status=spatial,residual_floor=floor,maximum_ratio=.85,
        fine_over_coarse=f/c if c else None,timestep_difference=delta,temporal_limit=limit,
        temporal_status='pass' if delta<=limit else 'unresolved')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();data=json.loads(a.source.read_text())
    manifest=json.loads((a.source.parent/'request-manifest.json').read_text())['request']
    require(RequestManifest(a.source.parent,manifest,True).completed(a.source.parent),'physical_budget_analysis_unsealed')
    require(data['execution_status']=='completed' and data['reconstruction_status']=='pass','physical_budget_reconstruction_incomplete')
    expected={(kind,equation) for kind in ('perturbed','baseline') for equation in ('inner_k','outer_k','outer_epsilon','whole_patch_k')}
    rows=[classify_row(row) for row in data['finite_defects']]
    require(len(rows)==len(expected) and {(r['kind'],r['equation']) for r in rows}==expected,'physical_budget_accuracy_coverage')
    quadrature=bool(data['quadrature']) and all(max(q['changes'].values())<=q['limit'] for q in data['quadrature'])
    status='pass' if quadrature and all(r['spatial_status'] in ('pass','unresolved_floor') and r['temporal_status']=='pass' for r in rows) else 'failed_or_unresolved'
    request=dict(schema='PhysicalBudgetAccuracyClassificationV1',source_sha256=sha(a.source),source_request_sha256=sha(a.source.parent/'request-manifest.json'),
        criteria='unchanged previously declared .85 spatial, max(5e-9,5*quadrature) floor, max(1e-9,.1*spatial difference) timestep bound')
    out=RequestManifest(a.output,request,False)
    result=dict(status=status,execution_status='completed',reconstruction_status='pass',configuration=data['configuration'],time=data['time'],
        quadrature_status='pass' if quadrature else 'unresolved',checks=rows,
        source=str(a.source),source_sha256=sha(a.source),
        scope='accuracy/convergence of finite physical flux defects; separate from successful reconstruction and prescribed-row reaction identities; original receipt remains unchanged')
    write_json(a.output/'results.json',result);out.seal(a.output,['results.json']);print(status)
    return 0 if status=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
