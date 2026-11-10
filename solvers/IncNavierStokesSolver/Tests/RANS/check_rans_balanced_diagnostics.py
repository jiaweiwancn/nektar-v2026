#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native canonical-observer eligibility and unchanged positive-x controls."""
import argparse,os,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import THREADS,write_json,exact_vector_hash
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load,Q
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info,write

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--previous-root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'diagnostic-regression';out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'));checks=[]
 for wall in ('EP','Accepted'):
  folder=out/wall;ref=Reference('E94E95','Natural',wall);case_xml(folder,'E94E95','Natural',wall,steps=4,dt=.0004,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=ref.knots(),degree=5,quadrature=Q[1],boundary_policy='FrozenReferenceV1')
  root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallPilotAccuracy','NeutralV1');write(root,folder/'case.xml')
  with (folder/'solver.log').open('w') as f:c=subprocess.run([str(r/'runtime/IncNavierStokesSolver'),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=600)
  new=load(folder/'case.evolution.json');old=load(a.previous_root/'p5/campaigns'/('E-N-'+wall+'-h0p5')/'T-0.0016/case.evolution.json')
  assert c.returncode==0 and not new['failure'] and new['time']==old['time']==.0016
  equal={k:(exact_vector_hash(new[k])==exact_vector_hash(old[k])) for k in ('state','initial','fixed_reference','native_history_times')}
  equal['accepted_states']=len(new['accepted_states'])==len(old['accepted_states']) and all(exact_vector_hash(x)==exact_vector_hash(y) for x,y in zip(new['accepted_states'],old['accepted_states'],strict=True))
  equal['native_history']=len(new['native_history'])==len(old['native_history']) and all(exact_vector_hash(x)==exact_vector_hash(y) for x,y in zip(new['native_history'],old['native_history'],strict=True))
  equal['state_audits']=new['state_audits']==old['state_audits'];assert all(equal.values()),equal
  assert 'not drift from a perturbed initial state' in new['pilot_accuracy_semantics']
  checks.append(dict(name='positive-x-'+wall,status='pass',bitwise_equal=equal,old_source=str(a.previous_root),new_source=str(folder),state_sha256=exact_vector_hash(new['state']),semantics=new['pilot_accuracy_semantics']))
  negative=out/(wall+'-negative-direction');negative.mkdir();root=ET.parse(folder/'case.xml').getroot()
  parameters=root.find('./CONDITIONS/PARAMETERS')
  for q in list(parameters):
   if q.text.split('=')[0].strip()=='RANSWallDirectionX':parameters.remove(q)
  ET.SubElement(parameters,'P').text='RANSWallDirectionX = -1'
  # A new output destination; no existing checkpoint is changed.
  set_info(root,'RANSWallCheckpoint',str(negative/'checkpoint'));write(root,negative/'case.xml')
  with (negative/'solver.log').open('w') as f:c=subprocess.run([str(r/'runtime/IncNavierStokesSolver'),'case.xml'],cwd=negative,stdout=f,stderr=subprocess.STDOUT,timeout=60)
  log=(negative/'solver.log').read_text();reason='NeutralV1 accuracy observer requires canonical positive-x reference direction';assert c.returncode!=0 and reason in log
  checks.append(dict(name='negative-direction-'+wall,status='pass',exit_code=c.returncode,reason=reason))
 write_json(out/'results.json',dict(status='pass',checks=checks));print('pass');return 0
if __name__=='__main__':sys.exit(main())
