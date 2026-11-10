#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Deterministic native CTest wrapper; retains every input on success/failure."""
import argparse,os,shutil,subprocess,sys,tempfile
from pathlib import Path
from rans_qualification import THREADS,write_json
from continue_cls_epawf_pimple import freeze

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--sampler',type=Path,required=True);p.add_argument('--descriptor-probe',type=Path,required=True);a=p.parse_args();os.environ.update(THREADS)
 r=Path(tempfile.mkdtemp(prefix='rans-frozen-native-',dir=os.environ.get('TMPDIR')));(r/'tmp').mkdir();freeze(a.solver.resolve(),r/'runtime')
 shutil.copy2(a.sampler,r/'runtime/RANSPreservationProbe');shutil.copy2(a.descriptor_probe,r/'runtime/RANSNewModelsProbe');rows=[]
 for script in ('check_rans_frozen_initial_history.py','check_rans_interface_sides.py'):
  command=[sys.executable,'-B',str(Path(__file__).with_name(script)),'--root',str(r)]
  with (r/(script+'.log')).open('w') as f:code=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT).returncode
  rows.append(dict(script=script,exit_code=code))
 status='pass' if all(x['exit_code']==0 for x in rows) else 'fail';write_json(r/'ctest-result.json',dict(status=status,checks=rows));print(status,str(r));return 0 if status=='pass' else 1
if __name__=='__main__':sys.exit(main())
