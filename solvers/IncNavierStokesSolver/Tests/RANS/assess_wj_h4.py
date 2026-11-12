#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact existing H4 physical boundary rows; no reconstruction or solve."""
import argparse
import csv
import gc
import gzip
from pathlib import Path
from run_wj_h4 import load,save,sha
from wj_h3_analysis import identity


def extract(root):
    out=root/'boundary-supplement-v1';out.mkdir();view=load(root/'native-observation-view.json')
    assert sha(view['path'])==view['sha256'];inputs=[];count=0;interfaces=[]
    with gzip.open(out/'physical-boundaries.csv.gz','wt',newline='')as f:
        writer=None
        for name,v in load(view['path']).items():
            native=Path(v['native']['path']);assert sha(native)==v['native']['sha256'];p=load(native)
            assert p['governing_solves']==0 and p['initialization_work']==[0,0]
            inputs.append(v['native'])
            for row in p['open_boundary']['boundary_rows']:
                r=dict(case=name,time=p['open_boundary']['time'],region='inner'if row['inner']else 'outer',role=row['role'],
                    geometry=row['geometry'],element=row['element'],x=row['physical'][0],z=row['physical'][2],chart_rho=row['chart'][2],
                    x_side='lower'if row['role']=='inlet'else 'upper'if row['role']=='outlet'else 'interior',
                    normal_side='lower'if row['role']=='ground'else 'upper'if row['role']=='top'else 'interior',
                    k=row['k'],epsilon_field=None if row['inner']else row['epsilon'],epsilon_scope='absent:inner'if row['inner']else 'outer_physical',
                    pi=row['pi'],sampled_K=row['K'],cutoff=row['cutoff'],physical_weight=row['weight'])
                for a,axis in enumerate('xyz'):
                    r['U'+axis]=row['U'][a];r['normal_'+axis]=row['normal'][a]
                    r['actual_traction_'+axis]=row['actual_traction'][a];r['applied_traction_'+axis]=row['applied_traction'][a]
                    for b,other in enumerate('xyz'):r['G'+axis+other]=row['G'][a][b];r['R'+axis+other]=row['R'][a][b]
                for j,block in enumerate(('u','w','k','epsilon')):
                    r['prescribed_weak_flux_'+block]=None if row['inner']and block=='epsilon'else row['flux'][j]
                if writer is None:writer=csv.DictWriter(f,fieldnames=list(r));writer.writeheader()
                writer.writerow(r);count+=1
            interfaces.extend(dict(case=name,time=p['open_boundary']['time'],**row)for row in p['open_boundary']['interface'])
            del p;gc.collect()
    save(out/'native-interface-rows.json',interfaces)
    save(out/'manifest.json',dict(source=identity(__file__),view=view,inputs=inputs,rows=count,
        outputs=[identity(out/'physical-boundaries.csv.gz'),identity(out/'native-interface-rows.json')],
        role='exact existing native physical-row data; outlet actual/applied tractions distinct; no inner epsilon field; no governing work'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);extract(p.parse_args().root)
