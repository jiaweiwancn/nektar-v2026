#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Supplement sealed history controls with physical MMS budgets and corners."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import manufactured_fixture,fixture,functions,probe,run,require
from ct_accepted_pilot import replay
from rans_request_manifest import RequestManifest,runtime_identity,sha
from run_cls_epawf_picard import write


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--base',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');base=json.loads(a.base.read_text());require(base['status']=='pass','base_controls_failed')
    old=json.loads((a.base.parent/'request-manifest.json').read_text())['request'];identity=runtime_identity(a.runtime/'IncNavierStokesSolver');require(identity==old['runtime'],'foreign_base_controls_runtime')
    seal=json.loads((a.base.parent/'completed-artifacts.json').read_text());require(all(sha(a.base.parent/n)==h for n,h in seal['files'].items()),'base_controls_changed')
    manifest=RequestManifest(a.output,dict(schema='CTAcceptedOpenNativeControlsSupplementV1',runtime=identity,
        base_control_path=str(a.base),base_control_sha256=sha(a.base),dt=.001,steps=2,
        physical_budget_limit=2e-9,source='manufactured verification only, never the physical pilot'))
    result=dict(status='fail',checks=base['checks'].copy(),reused=dict(path=str(a.base),sha256=sha(a.base)),executions=[])
    try:
        folder=a.output/'manufactured';manufactured_fixture(folder,a.runtime,steps=2,dt=.001)
        initial,e=probe(a.runtime,folder,dict(operation='functions',grid=[],export_basis=True),'initial');result['executions'].append(e);require(e['exit_code']==0,'MMS_initialization_failed')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');result['executions'].append(e)
        r=json.loads((folder/'case.evolution.json').read_text());require(e['exit_code']==0 and not r['failure'] and abs(r['time']-.002)<1e-13,'MMS_native_advancement_failed')
        snapshot=json.loads(Path(r['pimple_attempts'][-1]['k_snapshot']).read_text())
        request=replay(initial,r['state'],r['time'],[],export_basis=True,export_physical_scalar_budget=True,
            previous_state=r['accepted_states'][-2],stage_history=snapshot['history'],export_frozen_scalar_reconciliation=True,physical_face_order=32)
        packet,e=probe(a.runtime,folder,request,'budget');result['executions'].append(e);require(e['exit_code']==0,'MMS_budget_replay_failed')
        b=packet['physical_scalar_budget'];defects={name:abs(row['frozen_stage_physical_defect']) for name,row in b['equations'].items()}
        rec={name:max(abs(x) for x in row['reconstruction_defect']) for name,row in b['reconciliation']['assembled'].items()}
        require(max(defects.values())<=2e-9 and max(rec.values())<=2e-10,'independent_MMS_physical_budget_or_weak_reconstruction')
        o=packet['observation'];off=np.array(o['offsets'],int);state=np.array(o['state']);mass=0.
        for f in o['faces']:
            if f['interface']:continue
            q=f['velocity_k'];ids=np.array(q['ids'],int);u=[np.array(q['phi'])@state[ids+off[j]] for j in (0,1)]
            mass+=f['weight']*(u[0]*f['area'][0]+u[1]*f['area'][2])
        require(abs(mass)<=1e-8,'MMS_independent_mass_flux')
        result['checks']['independent_physical_storage_source_flux']=dict(status='pass',physical_defects=defects,
            weak_load_reconstruction_defects=rec,mass_outward=mass,physical_budget=b,
            scope='actual consecutive accepted states and frozen one-sweep stages, independently evaluated gradient/advective face fluxes; no row reaction substituted; manufactured verification only')
        folder=a.output/'bad_inlet_corner';path=fixture(folder,a.runtime);root=ET.parse(path).getroot()
        f=root.find("./CONDITIONS/FUNCTION[@NAME='BoundaryValuesOuter']/E[@VAR='eps_tilde']");f.set('VALUE','1');write(root,path)
        packet,e=probe(a.runtime,folder,dict(operation='functions',grid=[]),'corner');result['executions'].append(e)
        require(e['exit_code']!=0 and 'incompatible' in packet.get('failure','') and 'epsilon' in packet.get('failure',''),'inlet_dependent_corner_wrong_rejection')
        result['checks']['inlet_dependent_corner']=dict(status='pass',reason=packet['failure'])
        result['status']='pass'
    except Exception as error:result['failure']=repr(error)
    write_json(a.output/'results.json',result)
    if result['status']=='pass':manifest.seal(a.output,['results.json']+[str(f.relative_to(a.output)) for f in a.output.rglob('case.xml')])
    print(result['status'],result.get('failure',''));return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
