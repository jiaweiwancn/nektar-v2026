#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify the immutable evidence supporting the selected duration admission."""
import argparse,hashlib,json
from pathlib import Path
from rans_request_manifest import canonical
from qualify_rans_preservation import load
from rans_qualification import write_json

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--verify-existing',action='store_true');a=p.parse_args();r=a.root.resolve();names=['initial-gate-EP.json','initial-gate-Accepted.json','comparisons/E94E95-Natural-EP-0.0016.json','comparisons/E94E95-Natural-EP-0.016.json','comparisons/E94E95-Natural-Accepted-0.0016.json','comparisons/E94E95-Natural-Accepted-0.016.json','balanced-quadrature/results.json','balanced-fine-quadrature/results.json','evolved-quadrature-0.016/results.json','evaluation-uncertainty/results.json','diagnostic-regression/results.json','observer-perturbed-baseline-v2/results.json','common-test-selector-regression/results.json','evidence-contract/results.json']
 names += ['evolved-audits/balanced-E-N-'+w+'-'+v+'/results.json' for w in ('EP','Accepted') for v in ('B5','B11')];proof=[]
 for name in names:
  path=r/name;d=load(path);assert d['status']=='pass',(name,d['status']);proof.append(dict(path=name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),status=d['status']))
 expected=dict(schema='BalancedDurationAdmissionV1',status='pass',configurations=[['E94E95','Natural','Accepted']],reason='Nonconstant accepted wall is the limiting qualified E94E95 scalar case: larger full/common inner-k defects and reference error than EP; both E walls pass paired spatial and both-member dt uncertainty gates. A95/A93 qualification continues independently; no duration claim for those configurations.',proof=proof,cohort=dict(dt=.004,half_dt=.002,milestones=[.016,.1,.5],target_theta=.1,independent_initial_histories=True,unchanged_method_restart=True,optional_targeted_refinement_dt=.002),worker_allocation='initially two duration controllers beside at most four short-qualification controllers; numeric libraries one thread; aggregate <=72')
 path=r/'duration-admission.json'
 if a.verify_existing:assert canonical(load(path))==canonical(expected),'admission evidence changed'
 else:assert not path.exists();write_json(path,expected)
 output=r/'duration-admission-verification.json';assert not output.exists();write_json(output,dict(status='pass',admission_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),required_proofs=len(proof),scope='selected E94E95 Natural Accepted cohort only; independent model qualification retains its own gates'));return 0
if __name__=='__main__':raise SystemExit(main())
