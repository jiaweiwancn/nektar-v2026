#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One fresh native initialization with the first inner cell split at .0025.

This is a zero-time representation/source probe, not a finer-mesh trajectory
or a coefficient transfer/restart. It cannot establish time/spatial convergence.
"""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,run,write_json,require,state_hash,sha,write
from run_cls_epawf_pimple import prepare,check_report
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r2_cases import RADIAL,X


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--baseline-sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output=a.output.resolve();xml=prepare(a.output,.000625,0,True,checkpoint=True,k_bound='MassProjection',
        floor_relative=1e-8,inventory='ConstitutiveV2',wall_policy='UniformEPFlat')
    root=ET.parse(xml).getroot();tensor_mesh(root,[0.,.0025]+RADIAL['outer6'][1:],3,2,6.,1.,points=[8,8],xknots=X);write(root,xml)
    receipt=run([str(a.solver.resolve()),'case.xml'],a.output,'solver.log',7200);write_json(a.output/'receipt.json',receipt)
    require(receipt['exit_code']==0,'refined native initialization failed')
    r=load(a.output/'case.evolution.json');check_report(r)
    require(r['steps']==r['implicit_callback_count']==0 and not r['initial_k_bound_check']['activations'],'zero-time probe changed k or advanced')
    request=dict(operation='source_balance',state=r['state'],history=r['state'],time=0.,alpha=0.,
        inner_stress_model='SaturatedR2',require_history_h0=False,stations=[-3.5,-2.625,-1.,1.6],
        heights=[0.,1e-6,1e-5,1e-4,.0005,.001,.0025,.005,.0075,.01,.02,.1,1.,2.])
    write_json(a.output/'input.json',request)
    receipt=run([str(a.probe.resolve()),'case.xml','input.json','source.json'],a.output,'source.log',7200)
    require(receipt['exit_code']==0,'refined source readback failed');s=load(a.output/'source.json')
    baseline=load(a.baseline_sources)
    def station(source):return [v for v in source['points'] if v.get('valid') and v['variant']==-1 and v['x']==1.6]
    keys=('z','local_k','sampled_K','production','reference_production','source','represented_gradient','analytical_gradient','selected_native_spatial_rate')
    rows={name:[{k:v[k] for k in keys} for v in station(source)] for name,source in [('original_h0',baseline),('refined_h0',s)]}
    require(s['retained_state_unchanged'] and all(v['chi']==1 for v in s['points'] if v.get('valid') and v['variant']==-1),'refined policy changed')
    summary=dict(status='native_initialization_source_probe_completed',time=0.,production_steps=0,governing_transport_solves=0,
        mesh_change='split original [0,.01] first cell into [0,.0025] and [.0025,.01]; original p3 and all physical inputs retained',
        initialization='fresh native analytical neutral family; no checkpoint migration or same-polynomial claim',
        state_sha256=state_hash(r['state']),xml_sha256=sha(xml),initialization_lift_calls_iterations=r['initialization_lift_calls_iterations'],
        initial_raw_norms=r['raw_reference_norms'],points=rows,component_load_sum_defect=s['component_load_sum_defect'],
        component_mass_rate_sum_defect=s['component_mass_rate_sum_defect'],native_mass_solve_defect=s['native_mass_solve_defect'],
        limitation='zero-time representation/source comparison only; no mesh-converged evolution claim')
    write_json(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
