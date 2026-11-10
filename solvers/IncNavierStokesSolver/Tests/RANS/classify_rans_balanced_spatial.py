#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Frozen initialization gate; full and fixed-common-test defects both count."""
import argparse
from pathlib import Path
from qualify_rans_preservation import load
from rans_preservation_oracle import Reference
from run_rans_neutral_pilot import scaled,FIELDS
from rans_qualification import write_json

def classify(root,wall,model='E94E95',policy='Natural',screen=None):
    rows={x['name']:x for x in load(root/'diagnosis'/(screen or ('balanced-'+wall))/'results.json')['checks']};base=rows.get('BaseP5') or next(x for x in load(root/'diagnosis'/('initial-'+wall)/'results.json')['checks'] if x['name']=='BaseP5');c,f=(rows[n] for n in ('Balanced5V1','Balanced11V1'));checks={}
    for eq in (('u','k','w') if model in ('A93','A95') else ('u','k')):
        native=lambda x:x['spatial']['native_k_partition_mass_dual']['inner'] if eq=='k' else x['spatial']['native_full_mass_dual'][eq]
        for kind,fn in [('full_native',native),('fixed_common',lambda x:x['common_test']['fields'][eq]['fine_restricted_common_mass_dual'])]:
            v=[fn(x) for x in (base,c,f)];ratio=v[2]/v[1];checks[eq+'-'+kind]=dict(base=v[0],candidate=v[1],fine=v[2],fine_over_candidate=ratio,limit=.85,scaled_floor=1e-8,status='pass' if (v[1]<1e-8 or ratio<=.85) and v[2]<v[0] and abs(v[1]-v[2])>5e-10 else 'fail')
    ref=Reference(model,policy,wall);errors=[scaled(x['representation'],ref) for x in (c,f)]
    for index,name in enumerate(('candidate','fine')):
        checks[name+'-field_accuracy']=dict(status='pass' if all(errors[index][q]['weighted_RMS']<=(1e-4 if q=='v' else 1e-3) for q in FIELDS) else 'fail',errors=errors[index])
        checks[name+'-evolution_limit']=dict(status='pass' if (c,f)[index]['space']['dof']<=3000 else 'fail',unknowns=(c,f)[index]['space']['dof'],limit=3000)
    return dict(status='pass' if all(x['status']=='pass' for x in checks.values()) else 'fail',checks=checks,spaces=[c['space'],f['space']],other_equations={name:x['spatial']['native_full_mass_dual'] for name,x in [('candidate',c),('fine',f)]},interfaces={name:x['spatial']['native_k_partition_mass_dual'] for name,x in [('candidate',c),('fine',f)]},scope='initialization gate only; no evolved spatial or temporal conclusion')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--wall',choices=('EP','Accepted'),required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--model',default='E94E95');p.add_argument('--policy',default='Natural');p.add_argument('--screen');a=p.parse_args();assert not a.output.exists();r=classify(a.root,a.wall,a.model,a.policy,a.screen);write_json(a.output,r);print(r['status']);return 0 if r['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
