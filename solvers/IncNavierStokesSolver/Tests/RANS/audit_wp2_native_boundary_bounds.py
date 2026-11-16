#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native all-time wall/interface/continuity comparison bounds for WP2."""
import argparse,math
from pathlib import Path
from analyze_wp2_approach import paths,K0
from cao_tamura_h2d_state import load
from wp2_all_state import US
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic


def collect(root,model,phase):
    progress=root/f'preservation-{phase}-{model}/progress.json';p=load(progress);out={};sources=[];identity_error=0.
    require(p['status']=='completed_bounded_trajectory','boundary trajectory incomplete')
    for file in paths(root,model,phase,p):
        d=load(file);sources.append(dict(path=str(file),sha256=sha(file)));volume=math.fsum(d['wp2_common_samples']['weights'])
        require(len(d['accepted_states'])==len(d['wp2_per_state_monitoring'])==len(d['state_audits']),'native boundary inventory lengths differ')
        for state,monitor,audit in zip(d['accepted_states'],d['wp2_per_state_monitoring'],d['state_audits']):
            metrics=monitor['physical_checks'];g=metrics['unanchored_free_block_norms'];interface=audit['interface'];key=round(monitor['time'],12)
            independent=math.hypot(*audit['strong_continuity'][:2])/math.sqrt(volume)*.0065/US
            identity_error=max(identity_error,abs(independent-metrics['strong_divergence_RMS_over_ustar_H']))
            row=dict(time=monitor['time'],state_sha256=state_hash(state),wall_RMS_error=metrics['wall_shear_RMS_relative_error'],
                     wall_mean=metrics['wall_mean_shear'],interface_normal_RMS=abs(interface[2])/US**2,
                     interface_full_RMS=abs(interface[3])/US**2,interface_k_flux_RMS=abs(interface[5])/(US*K0[model]),
                     strong=independent,weak=max(g[4:7]),
                     net_flux=abs(metrics['boundary_status']['net_volume_flux'])/abs(next(f[3] for f in metrics['boundary_status']['faces'] if f[0]==0 and f[1]==-1)))
            if key in out:require(out[key]['state_sha256']==row['state_sha256'],'boundary restart changed coefficients')
            else:out[key]=row
    require(identity_error<=2e-12,'native regional L2 and physical RMS strong-divergence identities differ')
    return out,dict(trajectory=str(progress),sha256=sha(progress),sources=sources,strong_divergence_identity_max_absolute_error=identity_error)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();result={}
    for model in ['CLS','standard']:
        data={};bindings={}
        for phase in ['coarse','fine','confirmation']:data[phase],bindings[phase]=collect(a.root,model,phase)
        comparisons={}
        for name,right in [('timestep','fine'),('mesh','confirmation')]:
            history=[]
            for key in sorted(set(data['coarse'])&set(data[right])):
                left,other=data['coarse'][key],data[right][key];require(abs(left['time']-other['time'])<2e-14,'boundary common time differs')
                bounds={k:left[k]+other[k] for k in ['wall_RMS_error','interface_normal_RMS','interface_full_RMS','interface_k_flux_RMS','strong','weak','net_flux']}
                limits={k:.0025 for k in bounds};limits.update(strong=.00025,weak=2.5e-10,net_flux=2.5e-10)
                gates={k:dict(value=v,limit=limits[k],passed=v<=limits[k]) for k,v in bounds.items()}
                history.append(dict(left_actual_time=left['time'],right_actual_time=other['time'],gates=gates,
                    mean_wall_shear_difference=abs(left['wall_mean']-other['wall_mean'])/US**2,passed=all(v['passed'] for v in gates.values())))
            comparisons[name]=dict(common_native_times=len(history),history=history,passed=all(x['passed'] for x in history))
        result[model]=dict(bindings=bindings,comparisons=comparisons)
    atomic(a.output,dict(version='WP2NativeBoundaryTriangleBoundsV1',models=result,
        interpretation='same physical faces and zero interface reference jump: triangle inequality bounds each paired RMS change; conservative sufficient bounds, not additive error percentages',
        source_sha256=sha(Path(__file__)),governing_solves=0))


if __name__=='__main__':main()
