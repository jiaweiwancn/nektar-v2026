#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Adjacent-dt physical comparisons for retained late-state native trajectories."""
import argparse,time,resource
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic
from wp2_all_state import Rows,US
from wp2_common_quadrature import resample
from wp2_qualification import common_differences,quarter_budget_differences
from wp2_independent_reference import evaluate as analytical_reference


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--progress',type=Path)
    a=parser.parse_args();r=a.root;started=time.monotonic();progress_path=a.progress or r/'late-standard-dt005/progress.json';progress=load(progress_path)
    require(progress['import_validation']['passed'],'unqualified local import')
    old_path=r/'preservation-recovery-standard/step368/case.evolution.json';context=load(old_path)
    coarse={round(t,12):(t,s) for t,s in zip(context['audit_times'],context['accepted_states'])}
    fine={};bindings={str(old_path):sha(old_path),str(progress_path):sha(progress_path)}
    for label,entry in progress['results'].items():
        if label!='fine-initial' and not label.startswith('fine-step'):continue
        path=Path(entry['report']);require(sha(path)==entry['report_sha256'],'fine report changed');d=load(path);bindings[str(path)]=sha(path)
        for t,s in zip(d['audit_times'],d['accepted_states']):
            key=round(t,12)
            if key in fine:require(state_hash(fine[key][1])==state_hash(s),'restart boundary state differs')
            fine[key]=(t,s)
    rows=Rows(r/'geometry-rows-selected/result.json');seed=fine[min(fine)][1];K0=context['parameters'][5];H=context['parameters'][1]
    def evaluate(state):
        v=rows.tuple(np.asarray(state),context)
        dense=resample(dict(context,wp2_common_samples=dict(context['wp2_common_samples'],values=v['values']),
                            physical_samples=dict(wall_K=v['endpoint_K'])),28)
        sample=dict(columns=context['wp2_common_samples']['columns'],values=dense['values'],weights=dense['weights'].tolist(),regions=dense['regions'].tolist())
        # Keep the original analytical normalization. Neither the represented
        # cold reference nor the late seed replaces that denominator.
        sample['reference_values']=reference_values
        return sample,dense
    rv=rows.tuple(np.asarray(context['fixed_reference']),context)
    reference_dense=resample(dict(context,wp2_common_samples=dict(context['wp2_common_samples'],values=rv['values']),physical_samples=dict(wall_K=rv['endpoint_K'])),28)
    reference_values=analytical_reference(context,reference_dense['physical'])['values']
    seed_sample,_=evaluate(seed);history=[]
    for key in sorted(set(coarse)&set(fine)):
        tc,sc=coarse[key];tf,sf=fine[key];require(abs(tc-tf)<2e-13,'common times differ; no interpolation allowed')
        left,dl=evaluate(sc);right,dr=evaluate(sf)
        difference=common_differences(left,right,US,K0,H)
        cr=common_differences(left,seed_sample,US,K0,H);fr=common_differences(right,seed_sample,US,K0,H)
        response={}
        for name in ['U','W','k','epsilon','pi','mechanical_pressure']:
            den=fr[name]['dimensional_L2'];response[name]=None if den<=1e-12 else difference[name]['dimensional_L2']/den
        tensors={}
        for side in [0,1]:
            mask=dl['regions']==side;w=dl['weights'][mask]
            for name in ['R','T']:tensors[str(side)+'_'+name]=float(np.sqrt(np.sum(w[:,None,None]*(dl[name][mask]-dr[name][mask])**2)/sum(w))/US**2)
        minima=[]
        for data in [dl,dr]:
            R=data['R'];v=np.minimum(R[:,1,1],.5*(R[:,0,0]+R[:,2,2]-np.hypot(R[:,0,0]-R[:,2,2],2*R[:,0,2])))
            minima.append([float(min(v[data['regions']==side])) for side in [0,1]])
        history.append(dict(coarse_actual_time=tc,fine_actual_time=tf,coarse_state_sha256=state_hash(sc),fine_state_sha256=state_hash(sf),
            differences=difference,coarse_response_from_common_seed=cr,fine_response_from_common_seed=fr,difference_fraction_of_fine_response=response,
            quarter_budget=quarter_budget_differences(left,right,US,K0,H),full_tensor_RMS_difference_over_ustar2=tensors,
            full_tensor_quarter_budget_passed=max(tensors.values())<=.0025,GL28_covariance_inner_outer_coarse_fine=minima))
    atomic(a.output,dict(version='WP2LateCommonStateComparisonV2',fine_native_status=progress['status'],history=history,
        common_seed_sha256=state_hash(seed),fixed_reference_sha256=state_hash(context['fixed_reference']),bindings=bindings,
        response_denominator_floor_dimensional_L2=1e-12,pressure='absolute, no gauge realignment',quadrature='seen GL28 physical union partition, exact affine coefficient-derived fields; dependent inner epsilon excluded from transported field norms',
        source_sha256=sha(Path(__file__)),seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Local common-state adjacent timestep differences only; no temporal order, global accumulated-error or recovery qualification claim. No PDE solves or changes to retained states.'))


if __name__=='__main__':main()
