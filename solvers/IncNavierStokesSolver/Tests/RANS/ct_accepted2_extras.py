#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent high-order initial normal-load control and computed-flow replay."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from ct_accepted_feasibility import write_json
from ct_accepted2_diagnosis import analytic,make_case
from ct_accepted2_storage import artifact
from ct_accepted_open import require,run
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,write
from rans_request_manifest import RequestManifest,sha,runtime_identity

def high_order(xml,projection,order):
    root=ET.parse(xml).getroot();z=np.array(projection['nodes']);p=np.array(projection['independent_nodal_projection']);x,w=leggauss(order);total=0.
    for j,(a,b) in enumerate(zip(z[:-1],z[1:])):
        q=(a+b)/2+(b-a)*x/2;projected=p[j]*(b-q)/(b-a)+p[j+1]*(q-a)/(b-a)
        derivative=-2*(q-6.25)/1.5**2*np.exp(-((q-6.25)/1.5)**2)
        # Analytical Qzz=-Pi on the exact 1D supplied initial field.
        total+=(b-a)/2*np.dot(w,(analytic(root,'Outer','p',q)-projected)*derivative)
    return float(total*.25**3/6)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();rt=r/'runtime';os.environ['LD_LIBRARY_PATH']=str(rt/'lib')
    manifest=RequestManifest(a.output,dict(schema='CT2IndependentProjectionAndReplayV1',runtime=runtime_identity(rt/'IncNavierStokesSolver'),source_sha256=sha(Path(__file__)),normal_quadrature=[48,96],checks='unchanged supplied Pi=-Qzz; no source correction; diagnostic snapshots Off and timing only'))
    result=dict(status='pass',projection_controls={})
    for layout in ('B','X','Zplus','XZplus','Zminus'):
        prep=r/('companion-initial-v1' if layout=='Zminus' else 'initial-v2')
        q=json.loads((prep/'results.json').read_text())['layouts'][layout]['0.002']['pressure_projection'];xml=prep/'work'/(layout+'-0.002')/'case.xml'
        lo,hi=[high_order(xml,q,n) for n in (48,96)];native=q['signed_common_test']['sum']
        result['projection_controls'][layout]=dict(native_common_normal_load=native,independent48=lo,independent96=hi,quadrature_difference=abs(lo-hi),native_vs_tight_defect=abs(native-hi),scope='initial reference field only, fixed phi=x*(.25-x)*exp(-((z-6.25)/1.5)^2); independent analytical Q and represented P1 pressure; not a perturbed/evolved equilibrium assertion')
    # Read the actual frozen binary arrays, preserving their generating row IDs.
    report=json.loads((r/'matrix-v1/work/B-0.002/case.evolution.json').read_text());last=report['pimple_attempts'][-1]
    f=json.loads(artifact(last['flow_snapshot'],last['flow_snapshot_sha256']).read_text())
    def matrix(name):
        q=f[name];return np.fromfile(artifact(q['file'],q['sha256']),dtype=np.float64).reshape(q['rows'],q['columns'])
    A=matrix('frozen_flow_matrix');M=matrix('free_velocity_mass');nv=f['velocity_size'];G=A[:nv,nv:];D=A[nv:,:nv];C=A[nv:,nv:]
    p=np.array(f['predictor_increment']);dz=np.array(f['pressure_increment']);bg=np.array(f['subtracted_mass_velocity_correction']);res=np.array(f['frozen_flow_residual']);du=np.array(f['velocity_increment'])
    defects=dict(predictor=float(np.max(abs(A[:nv,:nv]@p+res[:nv]))),mass_correction=float(np.max(abs(f['alpha']*M@bg-G@dz))),pressure=float(np.max(abs(np.array(f['pressure_operator'])@dz-D@p-res[nv:]))),continuity=float(np.max(abs(D@du+C@dz+res[nv:]))),splitting=float(np.max(abs(A@np.r_[du,dz]+res-np.r_[G@dz-A[:nv,:nv]@bg,np.zeros(len(dz))]))))
    require(max(defects.values())<2e-10,'binary_computed_flow_replay_failed')
    result['binary_flow_reconstruction']=dict(status='pass',maximum_absolute_defects=defects,snapshot_sha256=last['flow_snapshot_sha256'],additional_governing_solves=0,method='matrix-vector checks of saved computed increments, no new linear system solved')
    source=r/'matrix-v1/work/B-0.002';c=json.loads((source/'construction.json').read_text());folder=a.output/'snapshots-off-timed';make_case(folder,source/'case.xml',c,'B',.002,2)
    root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallPimpleSnapshots','Off');write(root,folder/'case.xml')
    os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([rt/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
    new=json.loads((folder/'case.evolution.json').read_text())
    require(e['exit_code']!=0 and same_history(report,new) and report['rejected_state']==new['rejected_state'] and report['failure']==new['failure'],'snapshot_or_timing_instrumentation_changed_state_history_failure')
    timing={}
    for line in (folder/'timing.jsonl').read_text().splitlines():
        t=json.loads(line);row=timing.setdefault(t['name'],dict(exclusive_process_cpu_seconds=0.,exclusive_wall_seconds=0.,calls=0))
        row['exclusive_process_cpu_seconds']+=t['exclusive_process_cpu_seconds'];row['exclusive_wall_seconds']+=t['exclusive_wall_seconds'];row['calls']+=1
    result['instrumentation_compatibility']=dict(status='pass',accepted_and_rejected_states_and_native_history='bitwise equal',failure='unchanged negative A95 radicand at attempt .004',snapshot_off=True,execution=e,timing=timing,report_sha256=sha(folder/'case.evolution.json'),expected_process_exit='nonzero physical rejection, successful regression')
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']);print('extras pass');return 0
if __name__=='__main__':raise SystemExit(main())
