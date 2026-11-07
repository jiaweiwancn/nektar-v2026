#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact frozen-step increment balance in shared physical bands, using native M_ff.

These projected discrete loads are distinct from continuum source/flux estimates.
There are no transport solves in this readback.
"""
import argparse
import os
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from run_cls_epawf_picard import load,run,require,write_json,sha,THREADS


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True);os.environ.update(THREADS);out={}
    for case in 'ABCDEFGH':
        root=a.root/'screening'/case
        if not (root/'result.json').exists():continue
        initial=load(root/'initial/case.evolution.json');folder=a.output/case;folder.mkdir(exist_ok=True)
        if not (folder/'mass-and-bands.json').exists():
            write_json(folder/'input.json',dict(operation='contact',state=initial['state'],history=initial['state'],time=0.,alpha=0.,h=[],export_mass=True,
                require_history_h0=False,inner_stress_model='SaturatedR2'))
            rc=run([str(a.probe),str(root/'initial/case.xml'),'input.json','mass-and-bands.json'],folder,'probe.log',1800);write_json(folder/'receipt.json',rc)
            require(rc['exit_code']==0,'physical mass/band-row replay failed '+case)
        mass=load(folder/'mass-and-bands.json');require(mass['state_sha256']==initial['height_comparison_final']['state_sha256'],'cached mass replay belongs to a different case')
        M=np.array(mass['free_k_mass']);rows=np.array(mass['free_k_rows'],dtype=int)
        offset=int(mass['offsets'][2]);size=int(mass['offsets'][3])-offset;index=rows-offset;factor=cho_factor(M,lower=True)
        bands=mass['fixed_physical_band_rows']['regions'];ell=np.array([b['coefficient_integral_row'] for b in bands]);free=ell[:,index]
        dual=cho_solve(factor,free.T).T
        loads=initial['height_comparison_final']['weak_k_component_loads'];initial_components={};component_sum=np.zeros(len(rows))
        for name in ('conservative_advection','k_diffusion','inner_k_production','inner_k_destruction','outer_k_production','outer_k_destruction','prescribed_scalar_source'):
            v=np.array(loads[name])[index];component_sum+=v
            initial_components[name]=dict(all_free_modes_mass_dual_norm=float(np.sqrt(max(0.,v@cho_solve(factor,v)))),
                projected_band_rates=(-dual@v).tolist(),maximum_free_load=float(max(abs(v),default=0.)))
        require(max(abs(component_sum-np.array(loads['spatial_sum'])[index]),default=0.)<1e-12,'initial physical components do not sum; possible double counting')
        totals={k:np.zeros(len(bands)) for k in ('input_load','implicit_increment','bound','actual','closure')};steps=[]
        report=load(root/'trajectory/case.evolution.json')
        for attempt in report['pimple_attempts']:
            if attempt['time']>report['time']+1e-13:continue
            packet_path=Path(attempt['k_snapshot']);require(sha(packet_path)==attempt['k_snapshot_sha256'],'native packet differs')
            q=load(packet_path);require(np.array_equal(rows,np.array(q['transport_rows'],dtype=int)),'native free-row definition differs')
            require(sha(Path(q['matrix_file']))==q['matrix_sha256'],'frozen transport matrix differs')
            A=np.fromfile(q['matrix_file'],dtype='<f8').reshape((len(rows),len(rows)));alpha=q['alpha'];h=1/alpha
            pre=np.array(q['pre_k']);raw=np.array(q['raw']);corrected=np.array(q['corrected']);delta=np.array(q['raw_increment']);residual=np.array(q['transport_residual'])
            repair=(corrected-raw)[rows];actual=(corrected-pre)[rows]
            require(max(abs((raw-pre)[rows]-delta),default=0.)<1e-14,'raw increment mismatch')
            prescribed=sorted(set(range(offset,offset+size))-set(rows))
            require(np.array_equal(pre[prescribed],raw[prescribed]) and np.array_equal(raw[prescribed],corrected[prescribed]),'prescribed k changed')
            require(max(abs(pre[offset:offset+size]-np.array(q['history'])[offset:offset+size]),default=0.)<1e-13,'k history differs; storage must not be called spatial input')
            linear=residual+A@delta;implicit=A@delta-alpha*(M@delta)
            values=dict(input_load=-h*(dual@residual),implicit_increment=-h*(dual@implicit),bound=free@repair,actual=free@actual)
            values['closure']=values['actual']-values['input_load']-values['implicit_increment']-values['bound']
            added=dual@(A@repair);storage=alpha*(free@repair)
            corrected_identity=dual@(residual+A@(delta+repair))-added
            require(max(abs(linear),default=0.)<1e-12 and max(abs(values['closure']),default=0.)<1e-12 and max(abs(corrected_identity),default=0.)<1e-10,'frozen physical-band identity failed')
            for k,v in values.items():totals[k]+=v
            steps.append(dict(time=q['time'],increments={k:v.tolist() for k,v in values.items()},added_frozen_load_rates=added.tolist(),storage_reaction_rates=storage.tolist(),
                raw_linear_residual=float(max(abs(linear),default=0.)),corrected_identity_defect=float(max(abs(corrected_identity),default=0.))))
        initial_k=np.array(initial['state'])[offset:offset+size];final_k=np.array(report['state'])[offset:offset+size]
        endpoint=ell@(final_k-initial_k)
        require(max(abs(endpoint-totals['actual']),default=0.)<1e-12,'step budget did not telescope to physical band storage')
        out[case]=dict(bands=[dict(name=b['name'],lower=b['lower'],upper=b['upper']) for b in bands],steps=steps,totals={k:v.tolist() for k,v in totals.items()},
            initialized_weak_components=initial_components,initialized_projected_band_rates=(-dual@component_sum).tolist(),
            endpoint_integral_change=endpoint.tolist(),production_governing_solves=0,diagnostic_k_transport_solves=mass['diagnostic_k_solves'],
            mass_source='native PicardKMass on identical free rows; saved matrix is A_k, not M_ff',
            meaning='all totals have units of integrated k; input=-h ell^T Mff^-1 r_spatial, implicit=-h ell^T Mff^-1 L Delta_raw; correction=ell^T delta_bound')
        write_json(a.output/'summary.json',out)


if __name__=='__main__':main()
