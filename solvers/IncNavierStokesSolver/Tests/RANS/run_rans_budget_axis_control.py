#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Two predeclared, under-limit spatial diagnostics after failed budget accuracy.

This wraps the unchanged pilot for a new destination. It never changes an
existing request, trajectory, boundary artifact, timestep or qualification gate.
The original B5/B11 accuracy decision is not replaced by these controls.
"""
import argparse
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest, sha
from rans_duration_contract import require, identity
from rans_ownership_contract import OPERATOR
import run_rans_neutral_pilot as pilot


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','hypothesis','binding','mesh-file'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--control',choices=('NormalEndpointsV1','StreamwiseBaseV1','StreamwiseP6V1'),required=True)
    p.add_argument('--label',required=True)
    p.add_argument('--dt',type=float,choices=(.002,.001),required=True)
    p.add_argument('--amplitude',type=float,choices=(0.,2e-5),required=True)
    p.add_argument('--initialize-only',action='store_true')
    p.add_argument('--resume',action='store_true')
    a=p.parse_args();r=a.root.resolve();cfg=dict(model='A95',policy='Natural',wall='Accepted')
    failed=load(a.binding);binding=failed['binding'];identity(binding['configuration'],cfg)
    require(binding['numerical_operator']==OPERATOR,'axis_control_old_operator')
    require(binding['amplitude']==2e-5 and abs(binding['time']-.016)<=1e-13,'axis_control_foreign_response')
    require(failed['reconstruction_status']=='pass' and failed['accuracy_status']=='failed_or_unresolved',
            'axis_control_requires_completed_failed_accuracy')
    br=load(a.binding.parent/'request-manifest.json')['request']
    require(RequestManifest(a.binding.parent,br,True).completed(a.binding.parent),'axis_control_unsealed_failure')
    h=load(a.hypothesis);identity(h['configuration'],cfg)
    require(h['schema']=='PhysicalBudgetAxisHypothesisV1' and h['control']==a.control,'axis_control_hypothesis_identity')
    require(h['mesh_sha256']==sha(a.mesh_file),'axis_control_predeclared_mesh_changed')
    require(h['source_primary_budget_sha256']==sha(Path(h['source_primary_budget'])),
            'axis_control_diagnosis_source_changed')
    base_path=r/'meshes/A95-N-Accepted/Balanced5V1.json';base=load(base_path)
    require(h['base_mesh_sha256']==sha(base_path),'axis_control_base_mesh_changed')
    spec=load(a.mesh_file)
    field_degrees=[6,5] if a.control=='StreamwiseP6V1' else [5,5]
    pressure_degrees=[v-2 for v in field_degrees]
    knots=sorted(set(base['knots']+[.000625,3.75])) if a.control=='NormalEndpointsV1' else base['knots']
    require(spec['knots']==knots and spec['degree']==5 and spec['field_degrees']==field_degrees and
            spec['pressure_degrees']==pressure_degrees and spec['control']==a.control,'axis_control_space_mismatch')
    expected={'StreamwiseP6V1':2769,'StreamwiseBaseV1':2355,'NormalEndpointsV1':2589}[a.control]
    require(h['expected_unknowns']==expected and expected<=3000,'axis_control_backend_limit')
    sample=pilot.probe
    def checked_probe(*args,**kwargs):
        result=sample(*args,**kwargs)
        require(result['observation']['numerical_operator']==OPERATOR,'axis_control_runtime_operator_mismatch')
        require(sum(result['observation']['sizes'])==expected,'axis_control_actual_unknown_count_mismatch')
        for b,elements in enumerate(result['observation']['field_basis_configuration']):
            expected_degrees=pressure_degrees if b>=4 else field_degrees
            for element in elements:
                require([axis[1]-1 for axis in element['axes_type_modes_pointsType_points']]==expected_degrees,
                        'axis_control_consumed_directional_basis_mismatch')
        return result
    pilot.probe=checked_probe
    initialize=pilot.configure
    def configure(folder,args,ref,steps,restart=None):
        initialize(folder,args,ref,steps,restart)
        if a.control=='StreamwiseP6V1':
            xml=folder/'case.xml';root=ET.parse(xml).getroot()
            for e in root.findall('./EXPANSIONS/E'):
                fields=e.get('FIELDS','').split(',')
                degrees=pressure_degrees if fields==['p'] else field_degrees
                require(e.get('BASISTYPE')=='Modified_A,Modified_A','axis_control_basis_changed')
                require(e.get('NUMPOINTS')=='8,8','axis_control_native_points_changed')
                e.set('NUMMODES',','.join(str(v+1) for v in degrees))
            pilot.write(root,xml)
    pilot.configure=configure
    factory=pilot.RequestManifest
    class BoundRequest(factory):
        def __init__(self,directory,request,resume=False):
            request=dict(request,spatial_control=dict(schema='PhysicalBudgetAxisControlV1',control=a.control,
                field_degrees=field_degrees,pressure_degrees=pressure_degrees,expected_unknowns=expected,
                hypothesis_sha256=sha(a.hypothesis),failed_accuracy_binding_sha256=sha(a.binding),
                driver_sha256=sha(Path(__file__)),scope='diagnosis only; original required budget-accuracy decision retained'))
            super().__init__(directory,request,resume)
    pilot.RequestManifest=BoundRequest
    sys.argv=[str(Path(pilot.__file__)), '--root',str(r),'--label',a.label,'--model','A95','--policy','Natural',
        '--wall','Accepted','--mesh-file',str(a.mesh_file),'--degree','5','--dt',str(a.dt),
        '--velocity-amplitude',str(a.amplitude),'--milestones','.016','--max-time','.016','--timeout','14400']
    if a.initialize_only:sys.argv+=['--initialize-only']
    if a.resume:sys.argv+=['--resume']
    code=pilot.main()
    metrics=load(r/'campaigns'/a.label/'initial-metrics.json')
    require(sum(metrics['dof'])==expected,'axis_control_actual_unknown_count_mismatch')
    return code


if __name__=='__main__':raise SystemExit(main())
