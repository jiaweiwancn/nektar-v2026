#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Execution-only adapters for the IMMUTABLE H13 six-millimetre target.

No calibration, roughness inversion, profile, height-family or mesh generation.
"""
import copy
import hashlib
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from cls_epawf_m3_cases import set_info, set_parameter, write
from wj_h13_domain import topology

EXPECTED=dict(Nx=14,Ni=16,No=34,state=24688)
QUALIFICATION=('mass','operators','sweep','pressure','bound','initializer','restart','rollback','cache','runtime','dense_default')
INFO={'RANSWallLinearBackend','RANSWallMappedResourceDesign','RANSWallReport','RANSWallCheckpoint','RANSWallRestart',
      'RANSWallPimpleSnapshotDirectory','RANSWallPimpleArithmetic','RANSWallPimpleKBoundTraceFile','RANSWallPimpleKBoundDumpFile'}
PARAMETERS={'NumSteps','TimeStep','IO_CheckSteps','IO_InfoSteps'}

def file_id(path):
    p=Path(path);return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def native_command(admission,phase,argv):
    if phase not in ('target_cp0','pilot')or admission.get('status')!='PASS':raise ValueError('closed target command')
    if list(map(str,argv))!=admission['commands'][phase]:raise ValueError('unbound target command')
    for v in admission['command_inputs'][phase]:
        if file_id(v['path'])!=v:raise ValueError('changed immutable target command input')

class OwnerIndex:
    """Candidate index only; the established owner validator defines every row."""
    def __init__(self,report):
        import numpy as np
        self.report=report;self.mesh=report['fields'][0]['mesh']
        self.axes=[sorted({e[k][a]for e in self.mesh for k in('lower','upper')})for a in(0,2)]
        for z in self.axes:
            if min(np.diff(z))<=256*np.finfo(float).eps*max(1.,max(abs(v)for v in z)):
                raise ValueError('owner index needs separated frozen tensor intervals')
        ix={v:i for i,v in enumerate(self.axes[0])};iz={v:i for i,v in enumerate(self.axes[1])}
        self.cells={}
        for j,e in enumerate(self.mesh):
            key=(ix[e['lower'][0]],iz[e['lower'][2]])
            if key in self.cells:raise ValueError('duplicate tensor owner cell')
            self.cells[key]=(j,e)
    def __call__(self,p):
        from bisect import bisect_left
        from wj_h4_points import chart,owner_point
        q=chart(self.report,p);indices=[]
        for z,a in zip(self.axes,(0,2)):
            j=bisect_left(z,q[a]);indices.append(range(max(0,j-2),min(len(z)-1,j+2)))
        rows=[]
        for i in indices[0]:
            for k in indices[1]:
                item=self.cells.get((i,k))
                if item is None:continue
                order,e=item
                try:o=owner_point(self.report,e,p)
                except ValueError:continue
                rows.append((order,o))
        if not rows:raise ValueError('physical point has no native owner')
        return [o for _,o in sorted(rows,key=lambda v:v[0])]

def protected(root):
    r=copy.deepcopy(root)
    for parent,tag in (('./CONDITIONS/SOLVERINFO','I'),('./CONDITIONS/PARAMETERS','P')):
        node=r.find(parent)
        for v in list(node):
            key=v.get('PROPERTY')if tag=='I'else v.text.split('=')[0].strip()
            if key in (INFO if tag=='I'else PARAMETERS):node.remove(v)
    # Whitespace-only XML serialization and comments have no physical role.
    for node in r.iter():
        if node.text is not None and not node.text.strip():node.text=None
        if node.tail is not None and not node.tail.strip():node.tail=None
    return ET.tostring(r)

def validate_pair(original,execution):
    old=ET.parse(original).getroot()if not isinstance(original,ET.Element)else original
    new=ET.parse(execution).getroot()if not isinstance(execution,ET.Element)else execution
    if protected(old)!=protected(new):raise ValueError('execution diff changes protected physical/mesh/Function/Time data')
    p={v.text.split('=')[0].strip():v.text.split('=',1)[1].strip()for v in new.findall('./CONDITIONS/PARAMETERS/P')}
    if float(p['Time'])!=0. or float(p['RANSWallHeight']).hex()!=float(.15).hex():raise ValueError('frozen six-millimetre geometry/time differs')
    if not 0<=int(p['NumSteps'])<=4 or int(p['IO_CheckSteps'])!=1:raise ValueError('H14 zero/four-step/checkpoint ceiling differs')
    return dict(status='PASS',protected_sha256=hashlib.sha256(protected(old)).hexdigest(),
        original=file_id(original)if not isinstance(original,ET.Element)else None,
        execution=file_id(execution)if not isinstance(execution,ET.Element)else None,
        changed_execution_parameters={k:p[k]for k in PARAMETERS if k in p},no_height_or_profile_generation=True)

def clone(original,output,*,dt,steps=0,restart='None',design='None',backend='MappedBlockV1'):
    if backend not in ('Dense','MappedBlockV1'):raise ValueError('unknown backend')
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    root=ET.parse(original).getroot()
    for k,v in dict(RANSWallLinearBackend=backend,RANSWallMappedResourceDesign=design,RANSWallReport=str(output/'native.json'),
                    RANSWallCheckpoint=str(output/'CP'),RANSWallRestart=str(restart),RANSWallPimpleSnapshotDirectory=str(output/'snapshots'),
                    RANSWallPimpleArithmetic='Off',RANSWallPimpleKBoundDumpFile='None').items():set_info(root,k,v)
    for k,v in dict(TimeStep=dt,NumSteps=steps,IO_CheckSteps=1).items():set_parameter(root,k,v)
    path=output/'case.xml';write(root,path);diff=validate_pair(original,path)
    (output/'execution-diff.json').write_text(json.dumps(diff,indent=2)+'\n')
    return path

def target_preflight(original,execution,manifest,qualification,build,allocation):
    diff=validate_pair(original,execution)
    if manifest['primary_reference']['H_w'].hex()!=float(.15).hex():raise ValueError('different H13 primary')
    if {k:manifest['counts'][k]for k in EXPECTED}!=EXPECTED:raise ValueError('different H13 target dimensions')
    if topology(execution)['sha256']!=manifest['shared_topology']['sha256']:raise ValueError('different H13 target topology')
    missing=[k for k in QUALIFICATION if qualification.get(k)!='PASS']
    reasons=list(missing)
    if build.get('status')!='PASS':reasons.append('build')
    if allocation.get('status')!='PASS' or allocation.get('projected_total_seconds',float('inf'))>7200:reasons.append('resources')
    return dict(status='PASS'if not reasons else'OBSTRUCTED',native_launch_allowed=not reasons,
        height_mm=6,state_size=24688,missing_gates=reasons,execution_diff=diff,
        parent_H13_capacity_status='OBSTRUCTED; immutable historical dense-route result',
        backend='MappedBlockV1',fresh_state_required=True,old_strip_checkpoint_substitution=False,
        scientific_reference='unchanged H13 full-height neutral; physical interface9.5mm',
        field_spaces=manifest['counts'],actual_bandwidths='native permutation/factor audits must remain within predeclared allocation')
