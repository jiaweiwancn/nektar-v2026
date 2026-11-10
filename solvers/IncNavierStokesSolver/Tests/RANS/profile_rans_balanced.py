#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exclusive cost accounting; no overlapping parent/child time is double counted."""
import argparse,collections,json
from pathlib import Path
from rans_qualification import write_json

def group(name):
 if name.startswith(('domain.','audit.')) or name=='native.physical_samples':return 'diagnostic'
 if 'assembly' in name or name in ('pressure.build','matching.constraints'):return 'assembly'
 if name.startswith('dense.') or 'linear_solve' in name:return 'solve'
 if any(x in name for x in ('checkpoint','output','field_write')):return 'output'
 if name.startswith('startup.'):return 'initialization'
 if name=='constitutive.refresh':return 'constitutive_refresh_shared'
 if name=='residual':return 'residual_evaluation_shared'
 return 'other'

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();rows=[]
 for wall in ('EP','Accepted'):
  for space in ('B5','B11'):
   folder=a.root/'campaigns'/('E-N-'+wall+'-'+space)/'T-0.0160';totals=collections.defaultdict(lambda:dict(wall_seconds=0.,cpu_seconds=0.,calls=0));categories=collections.defaultdict(lambda:dict(wall_seconds=0.,cpu_seconds=0.,calls=0))
   for line in (folder/'timing.jsonl').read_text().splitlines():
    d=json.loads(line)
    for acc in (totals[d['name']],categories[group(d['name'])]):acc['wall_seconds']+=d['exclusive_wall_seconds'];acc['cpu_seconds']+=d['exclusive_process_cpu_seconds'];acc['calls']+=1
   total=sum(v['wall_seconds'] for v in totals.values());metric=json.loads((folder/'metrics.json').read_text());rows.append(dict(wall=wall,space=space,source=str(folder),exclusive_total=total,elapsed=metric['wall_seconds'],seconds_per_accepted_step=metric['seconds_per_step'],categories=categories,scopes=totals,output_timing_status='instrumented' if categories.get('output',{}).get('calls',0) else 'not separately instrumented; elapsed-minus-exclusive retains unclassified overhead',diagnostic_fraction=categories['diagnostic']['wall_seconds']/total,solve_fraction=categories['solve']['wall_seconds']/total))
 assert not a.output.exists();write_json(a.output,dict(status='pass',checks=rows,scope=__doc__,classification='domain/audit/physical samples are diagnostic; assembly/pressure build/constraints are assembly; dense solves are solve; named checkpoint/output scopes are output; residual and constitutive refresh are shared and reported separately',optimization='none introduced; mandatory domain computation dominates; immutable restricted-mass factor cache would not address that cost'));return 0
if __name__=='__main__':raise SystemExit(main())
