#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only component norms at each height's actual initialized state."""
import argparse
import math
import os
from pathlib import Path
from run_cls_epawf_picard import load,run,require,write_json,THREADS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screening',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+',default=['C','E','G']);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);os.environ.update(THREADS)
    out={}
    for c in a.cases:
        root=a.screening/c/'initial';initial=load(root/'case.evolution.json');q=initial['height_comparison_final'];folder=a.output/c;folder.mkdir()
        write_json(folder/'input.json',dict(operation='resolution_components',state=initial['state'],history=initial['state'],time=0.,alpha=0.,compute_scalar_rate=False,
            require_history_h0=False,inner_stress_model='SaturatedR2'))
        rc=run([str(a.probe),str(root/'case.xml'),'input.json','components.json'],folder,'probe.log',1800);write_json(folder/'receipt.json',rc)
        require(rc['exit_code']==0,'outer component replay failed '+c);d=load(folder/'components.json');H=q['H'];kappa=d['parameters'][7];E0=.4**3/(kappa*H)
        conversion=E0*.4/H*math.sqrt(6*(4-H))
        require(max(map(abs,d['native']['external_lift_delta']),default=0.)==0.,'component readback changed essential data')
        out[c]=dict(H=H,source_state_sha256=q['state_sha256'],epsilon_component_mass_dual={k:v[3]*conversion for k,v in d['norms'].items() if k.startswith('term_') and v[3]},
            epsilon_original_mass_dual=q['dimensional_mass_dual_spatial_first4'][3],all_source_vectors=str(folder/'components.json'),
            coefficients_unchanged=True,epsilon_trace_retained=True,production_governing_solves=0)
        write_json(a.output/'summary.json',out)


if __name__=='__main__':main()
