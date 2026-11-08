#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused case-policy and acceptance-contract checks; no numerical solver."""
from copy import deepcopy
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

from cls_epawf_m8_r1_cases import parameters,solver_info
from run_cls_accepted_awf import CASES,FORMULATION,INVENTORY,WORK,case_xml,check_report
from run_cls_epawf_picard import load,require


def rejects(fn,word):
    try:fn()
    except ValueError as e:
        require(word in str(e),str(e));return
    raise AssertionError('expected rejection: '+word)


def main():
    definitions={}
    with tempfile.TemporaryDirectory(prefix='accepted-awf-cases-') as tmp:
        root=Path(tmp)
        for case,(H,partition) in CASES.items():
            xml=case_xml(case,root/case,8);tree=ET.parse(xml).getroot();p=parameters(tree);s=solver_info(tree)
            require(s['RANSWALLFORMULATION']==FORMULATION and s['RANSWALLKSAMPLING']=='InterfaceTrace','wrong formulation/functional')
            require(float(p['RANSWALLSAMPLEHEIGHT'])==float(p['RANSWALLTAILSTART'])==float(p['RANSWALLHEIGHT'])==H,'endpoint/tail heights differ')
            require(s['RANSWALLTRANSITION']=='Fixed' and float(p['RANSWALLCHI'])==0 and
                s['RANSWALLDESTRUCTION']=='Original','legacy transition/cap leaked into new source')
            require(s['RANSWALLVARIANT']=='Original' and float(p['RANSWALLCUBICBETA'])==1 and
                float(p['RANSWALLREGULARIZATIONSCALE'])==0 and s['RANSWALLINNERSTRESS']=='Linear','faithful stress/outer model changed')
            require(s['RANSWALLPIMPLEKINVENTORY']==INVENTORY and s['RANSWALLPIMPLEKBOUNDSOLVER']=='RetainedPool','inventory/finite solver mismatch')
            definition=load(xml.parent/'definition.json');definitions[case]=definition
            require(definition['work']==WORK and definition['z_K']==definition['z_m']==H,'work/geometry policy differs')
            nodes=[list(map(float,x.text.split())) for x in tree.findall('GEOMETRY/VERTEX/V')]
            radial=sorted({x[1] for x in nodes});x=sorted({x[0] for x in nodes})
            require(radial==definition['radial_knots'] and x==definition['tangential_knots'],'actual XML mesh differs from declared mesh')
            if partition=='one':require([z for z in radial if z<=H]==[0.,H],'one-element layer still contains internal knots')
            require(p['SLOPE']=='Kappa*sqrt(K0/Qstar)' and p['YAPSCALE']=='K0^1.5/(2.5*E0)','reference helper changed calibrated slope or physical Yap scale')
            require({f.get('NAME') for f in tree.findall('CONDITIONS/FUNCTION')}=={'InitialPerturbationInner','InitialPerturbationOuter'},'unexpected forcing/boundary functions')
        for a,b in [('Accepted-AWF-2G','Accepted-AWF-2O'),('Accepted-AWF-05G','Accepted-AWF-05O')]:
            da,db=definitions[a],definitions[b];H=da['H']
            require([z for z in da['radial_knots'] if z>=H]==[z for z in db['radial_knots'] if z>=H] and
                da['tangential_knots']==db['tangential_knots'],'paired outer/tangential meshes differ')
        require(all([z for z in d['radial_knots'] if z>=2]==[z for z in definitions['Accepted-AWF-2G']['radial_knots'] if z>=2] for d in definitions.values()),'common mesh above2 changed')
        for case in ('Accepted-AWF-2G','Accepted-AWF-05G'):
            xml=case_xml(case,root/(case+'-cutoff-v1'),0,inner_refinement='cutoff-v1');d=load(xml.parent/'definition.json');base=definitions[case]
            require(d['radial_knots']==sorted(set(base['radial_knots']+[.0025,.005,.02,.04])),'focused refinement lost original physical knots')
            require([z for z in d['radial_knots'] if z>=d['H']]==[z for z in base['radial_knots'] if z>=d['H']],'inner follow-up changed outer mesh')
            require(d['model_choices']==base['model_choices'] and d['z_K']==d['z_m']==base['H'] and
                d['tensor_unknown_estimate']<=3000 and d['inner_refinement']=='cutoff-v1','focused representation changed model or pilot scope')
        rejects(lambda:case_xml('Accepted-AWF-05O',root/'hidden-refinement',0,inner_refinement='cutoff-v1'),'graded follow-up')
        xml=case_xml('Accepted-AWF-2G',root/'outer-pressure-p2',0,inner_refinement='cutoff-v1',outer_refinement='pressure-p2-v1',checkpoint_every=1)
        d=load(xml.parent/'definition.json');tree=ET.parse(xml).getroot()
        require(d['tensor_unknown_estimate']==2866 and d['pressure_degree']==dict(inner=1,outer=2),
            'outer pressure control changed pilot estimate or pressure degrees')
        require(d['radial_knots']==sorted(set(definitions['Accepted-AWF-2G']['radial_knots']+[.0025,.005,.02,.04])) and
            d['model_choices']==definitions['Accepted-AWF-2G']['model_choices'],'outer pressure control changed mesh or constitutive model')
        rows=tree.findall('./EXPANSIONS/E[@FIELDS="p"]')
        require(tree.find('./GEOMETRY/DOMAIN').text=='C[1,2]' and
            tree.find('./EXPANSIONS/E').get('COMPOSITE')=='C[1,2]',
            'regional expansions must use composites declared in DOMAIN')
        require([(row.get('COMPOSITE'),row.get('NUMMODES')) for row in rows]==[('C[1]','2,2'),('C[2]','3,3')],
            'regional pressure expansion did not preserve inner p1 and select outer p2')
        require(parameters(tree)['IO_CHECKSTEPS']=='1','focused control lacks exact last-accepted checkpoints')
        rejects(lambda:case_xml('Accepted-AWF-05G',root/'hidden-outer-control',0,outer_refinement='pressure-p2-v1'),'cutoff-resolved H2')
        rejects(lambda:case_xml('Accepted-AWF-2G',root/'hidden-fallback',0,companion='regularized'),'preserved faithful-model blocker')
    minimal=dict(temporal_coupling='PimpleSingleSweep',order=1,wall_policy=FORMULATION,eta=0.,
        k_bound_mode='MassProjection',inner_stress_model='Linear',implicit_callback_count=1,stage_audits=[],
        pimple_attempts=[dict(time=.000625,scheduled_complete=True,additional_governing_equation_solves=0,
            outer_sweeps=1,momentum_predictor_solves=1,pressure_system_solves=1,epsilon_transport_solves=1,
            k_transport_solves=1,bound_calls=1,k_bound=dict(success=False,failure='genuine projection rejection'))])
    # Four scheduled solves can complete without an accepted numerical step.
    require(not check_report(minimal)[0]['numerically_accepted'],'a rejected candidate was relabelled accepted')
    accepted=deepcopy(minimal);accepted['stage_audits']=[[.000625]];row=accepted['pimple_attempts'][0]
    row.update(mass_inverse_applications=2,pressure_gauge_size=2,additional_correction_or_restoration_solves=1,
        bound_kkt_solves=1,corrected_continuity_gauge=[0.]*3,boundary_defect=0.,momentum_splitting_identity_defect=0.,
        imposed_trace_defect=0.,trace_mass_column_max=1.,retained_trace_column_max=1.,epsilon_frozen_linear_residual=0.,
        k_frozen_linear_residual=0.,refreshed_original_block_norms=[5.]*8,
        k_bound=dict(success=True,activations=1,finite_pool_complete=True,primal_floor_defect=0.,ground_floor_defect=0.,
            work_controls=dict(finite_qp_solver='RetainedPool')))
    require(check_report(accepted)[0]['numerically_accepted'],'one-sweep acceptance incorrectly required original nonlinear convergence')
    broken=deepcopy(accepted);broken['pimple_attempts'][0]['k_bound']['primal_floor_defect']=1e-12
    rejects(lambda:check_report(broken),'complete finite/physical feasibility')
    broken=deepcopy(accepted);broken['pimple_attempts'][0]['k_transport_solves']=2
    rejects(lambda:check_report(broken),'extra scheduled pass')
    print('accepted-AWF case definitions, exact endpoint semantics, paired meshes, and rejection/acceptance distinction passed')


if __name__=='__main__':main()
