#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate local spatial and temporal statuses from sealed CT3 primitives."""
import argparse,copy,json,math,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from rans_request_manifest import RequestManifest,sha
from ct_accepted_open import require


def classify(result):
    out=copy.deepcopy(result)
    for contrast in out['contrasts'].values():
        evolved=[]
        for dt,rows in contrast['comparisons'].items():
            for q in rows:
                v=q['local'];v['original_combined_status']=v['status']
                if 'maximum_absolute_radicand_difference' not in v:
                    v['spatial_status']=v['status'];v['temporal_status']='unavailable_actual_domain'
                    continue
                d,m=v['maximum_absolute_radicand_difference'],v['minimum_margin_both_spaces']
                require(math.isfinite(d) and math.isfinite(m),'nonfinite_required_local_comparison')
                v['spatial_status']='failed_actual_domain' if m<0 else 'unresolved_floor' if d<=1e-10 else 'pass' if m>0 and d<=.25*m else 'unresolved_large_local_spatial_difference'
                if dt=='0.001':
                    td=v['maximum_dt_radicand_difference'];limit=max(1e-10,.1*d)
                    v['temporal_status']='pass' if math.isfinite(td) and td<=limit else 'unresolved_temporal_uncertainty'
                else:v['temporal_status']='not_primary_fine_dt_comparison'
                v['status']='pass' if v['spatial_status']==v['temporal_status']=='pass' else 'unresolved'
                if dt=='0.001' and q['time']>0:
                    for field,u in q['field_spatial_interpretation_uncertainty'].items():
                        limit=max(1e-7,.1*q['fields']['fields'][field]['scaled_weighted_RMS'])
                        u['limit']=limit;u['status']='pass' if math.isfinite(u['maximum_dt_difference']) and u['maximum_dt_difference']<=limit else 'unresolved_temporal_uncertainty'
                    evolved.append(q)
        contrast['status']='bounded_on_declared_grid' if evolved and all(q['fields']['status']=='pass' and q['local']['status']=='pass' and
            all(v['status']=='pass' for v in q['field_spatial_interpretation_uncertainty'].values()) for q in evolved) else 'unresolved_spatial_sensitivity'
    out['evolved_spatial_status']='bounded_contrasts_only' if out['contrasts'] and all(q['status']=='bounded_on_declared_grid' for q in out['contrasts'].values()) else 'unresolved'
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    m=RequestManifest(a.output,dict(schema='CT3IndependentLocalClassificationV1',input_sha256=sha(a.input),source_sha256=sha(Path(__file__)),criteria='unchanged predeclared quarter-margin and tenth-spatial-difference/floors'))
    result=classify(json.loads(a.input.read_text()));result['classification_scope']='read-only reclassification of the same sealed primitives; both failed criteria remain visible; no field/threshold change'
    # A positive global field status cannot hide a failed local spatial test,
    # even when the local timestep-uncertainty test also fails.
    q=result['contrasts']['degree_same_partition']['comparisons']['0.001'][-1]['local']
    require(q['spatial_status']!='pass','expected_recorded_local_spatial_limitation_missing')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print(result['evolved_spatial_status']);return 0
if __name__=='__main__':raise SystemExit(main())
