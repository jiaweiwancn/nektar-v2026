#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Zero-step native witness separating fixed-neutral change from actual drift."""
import argparse,os,subprocess
from pathlib import Path
from types import SimpleNamespace
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from run_rans_neutral_pilot import configure
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load
from rans_qualification import THREADS,write_json,exact_vector_hash

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();out=r/'observer-perturbed-baseline';out.mkdir(exist_ok=False);os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
 args=SimpleNamespace(root=r,model='A95',policy='Natural',wall='EP',dt=.0004,degree=5,level=0,refinement='Base',columns=1,quadrature_index=1,cache='ExactStepV1',velocity_amplitude=2e-5,mesh_file=None)
 configure(out,args,Reference('A95','Natural','EP'),0)
 with (out/'solver.log').open('w') as f:c=subprocess.run([str(r/'runtime/IncNavierStokesSolver'),'case.xml'],cwd=out,stdout=f,stderr=subprocess.STDOUT,timeout=300)
 d=load(out/'case.evolution.json');values=d['state_audits'][0]['pilot_accuracy_E_D_wall_gauge'];assert c.returncode==0 and not d['failure'] and d['time']==0 and d['steps']==0
 assert exact_vector_hash(d['state'])==exact_vector_hash(d['initial']) and max(values[5:7])>1e-8
 assert 'not drift from a perturbed initial state' in d['pilot_accuracy_semantics']
 write_json(out/'results.json',dict(status='pass',performed_steps=0,actual_initial_state_sha256=exact_vector_hash(d['initial']),state_sha256=exact_vector_hash(d['state']),actual_drift_from_initial=0.,native_change_from_fixed_neutral=values[5:10],analytical_reference_error=values[:5],native_semantics=d['pilot_accuracy_semantics']));return 0
if __name__=='__main__':raise SystemExit(main())
