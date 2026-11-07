#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Hold native one-element h0 fixed while changing source integration order."""
import argparse
import os
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from run_cls_epawf_picard import load,write_json,run,require,state_hash,THREADS
from run_cls_epawf_height_comparison import case_xml


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case',choices=('U05-one',))
    p.add_argument('--initial-report',type=Path)
    a=p.parse_args();require(bool(a.case)==bool(a.initial_report),'explicit case and retained report are paired')
    os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False);out={}
    def probe(case,order,operation,state,**extra):
        folder=a.output/(case+'-'+str(order)+'-'+operation)
        xml=case_xml(case,folder,0,normal=order,inventory='ConstitutiveV3' if a.case else 'ConstitutiveV2')
        write_json(folder/'input.json',dict(operation=operation,state=state,history=state,alpha=0.,time=0.,require_history_h0=False,inner_stress_model='SaturatedR2',**extra))
        rc=run([str(a.probe),str(xml),'input.json','result.json'],folder,'probe.log',1800);write_json(folder/'receipt.json',rc)
        require(rc['exit_code']==0,'fixed-coefficient quadrature probe failed');return load(folder/'result.json')
    inputs=[(a.case,a.initial_report)] if a.case else [('D',a.root/'screening/D/initial/case.evolution.json'),('H',a.root/'focused/initial-H/case.evolution.json')]
    for case,initial in inputs:
        d=load(initial);s=d['state'];low=probe(case,16,'height_comparison',s);high=probe(case,32,'height_comparison',s)
        require(low['state_sha256']==high['state_sha256']==state_hash(s),'quadrature replay changed selected state')
        require(low['representation_L2_velocity_gradient_production_referenceProduction'][:2]==high['representation_L2_velocity_gradient_production_referenceProduction'][:2],'integration changed represented velocity')
        mass=probe(case,16,'contact',s,h=[],export_mass=True);M=np.array(mass['free_k_mass']);factor=cho_factor(M,lower=True)
        rows=np.array(mass['free_k_rows'],dtype=int)-int(mass['offsets'][2]);sources={}
        for name,x in low['weak_k_component_loads'].items():
            x=np.array(x)[rows];y=np.array(high['weak_k_component_loads'][name])[rows];e=y-x
            norm=lambda v:float(np.sqrt(max(0.,v@cho_solve(factor,v))))
            sources[name]=dict(all_free_modes_mass_dual_difference=norm(e),baseline_mass_dual=norm(x),maximum_load_difference=float(max(abs(e),default=0.)))
        out[case]=dict(actual_h0_sha256=state_hash(s),normal_orders=[16,32],state_unchanged=True,
            wall_policy=d.get('wall_policy'),z_m=d.get('height_comparison_initial',{}).get('z_m'),
            physical_mass_solve_residual=mass['mass_solve_residual'],weak_loads=sources,
            dimensional_residual_norms=[low['dimensional_mass_dual_spatial_first4'],high['dimensional_mass_dual_spatial_first4']],
            fixed_bands=[low['fixed_bands'],high['fixed_bands']],production_governing_solves=0,diagnostic_k_transport_solves=mass['diagnostic_k_solves'],
            policy='fixed coefficients; same cutoff/roughness/geometric integration panels; mean remains fixed native GL8 dz/H')
        write_json(a.output/'summary.json',out)


if __name__=='__main__':main()
