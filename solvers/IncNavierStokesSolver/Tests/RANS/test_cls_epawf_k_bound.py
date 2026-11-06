#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed-horizon bounded experiments; no horizon shortening on failure.

Every case is independent from actual h0 unless an explicit checkpoint is
supplied. Phase snapshots and nonlinear-qualification exits stay separate.
"""
import argparse
import json
import math
from pathlib import Path
import platform
import subprocess
import struct
import xml.etree.ElementTree as ET

from run_cls_epawf_picard import load,run,sha,state_hash,require,write_json,set_info,set_parameter,write
from run_cls_epawf_pimple import ROOT,identity,prepare
from compare_cls_epawf_pimple import differences,same_history,describe


def frozen_accounting(packet, attempt):
    """Independent readback of the saved transport matrix, without a solve."""
    n = len(packet['transport_rows'])
    payload = Path(packet['matrix_file']).read_bytes()
    require(len(payload) == 8*n*n, 'frozen k matrix shape changed')
    matrix = struct.unpack('<'+str(n*n)+'d', payload)
    raw_delta = packet['raw_increment']
    correction = [packet['corrected'][int(i)]-packet['raw'][int(i)] for i in packet['transport_rows']]
    raw = []; added = []
    for i in range(n):
        raw.append(packet['transport_residual'][i]+math.fsum(matrix[i*n+j]*raw_delta[j] for j in range(n)))
        added.append(math.fsum(matrix[i*n+j]*correction[j] for j in range(n)))
    maximum = lambda values: max(map(abs, values), default=0.)
    result = dict(raw=maximum(raw), added=maximum(added),
                  corrected=maximum(x+y for x,y in zip(raw,added)),
                  raw_vector_readback_defect=maximum(x-y for x,y in zip(raw,packet['raw_linear_residual'])))
    require(result['raw_vector_readback_defect'] < 1e-12, 'raw frozen residual packet differs')
    if 'k_corrected_frozen_linear_residual' in attempt:
        require(abs(result['corrected']-attempt['k_corrected_frozen_linear_residual']) < 1e-12 and
                abs(result['added']-attempt['k_added_frozen_residual']) < 1e-12, 'corrected frozen residual accounting differs')
    return result


def replay(probe,xml,packet,output,*,floor=None,fixture=None):
    output.mkdir(parents=True,exist_ok=False)
    tree=ET.parse(xml).getroot();set_info(tree,'RANSWallPimpleKBound','MassProjection')
    if floor is not None:set_parameter(tree,'RANSWallPimpleKFloorRelative',floor)
    write(tree,output/'case.xml')
    if fixture:
        initial=packet
        request=dict(operation='k_bound',fixture=fixture,state=initial['state'],pre_k=initial['state'],
            history=initial['state'],time=0.,alpha=0.,require_history_h0=False,inner_stress_model='SaturatedR2')
    else:
        request=dict(operation='k_bound',state=packet['raw'],pre_k=packet['pre_k'],history=packet['history'],
            time=packet['time'],alpha=packet['alpha'],require_history_h0=False,inner_stress_model='SaturatedR2')
    write_json(output/'input.json',request)
    receipt=run([str(probe),str(output/'case.xml'),'input.json','result.json'],output,'probe.log',900)
    require(receipt['exit_code']==0,'native k-bound replay failed')
    result=load(output/'result.json');bound=result['k_bound']
    require(result['non_k_change']==result['essential_k_change']==0,'bound changed non-k or essential data')
    first,last=map(int,(result['offsets'][2],result['offsets'][3]))
    require(state_hash(result['raw'][:first]+result['raw'][last:])==state_hash(result['corrected'][:first]+result['corrected'][last:]),'non-k bitwise invariance failed')
    require(result['shared_cg_jump']<1e-12,'bound broke native CG continuity')
    require(abs(result['quadrature_correction_mass_norm']-bound['mass_norm_correction'])<1e-12 and
            abs(result['quadrature_correction_integral']-bound['integral_total'])<1e-12,'mass/integral accounting differs')
    write_json(output/'receipt.json',receipt)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--final-time',type=float,default=.0025)
    p.add_argument('--dts',type=float,nargs='+',default=[.000625,.0003125,.00015625])
    p.add_argument('--floor',type=float,default=1e-8)
    p.add_argument('--skip-baseline',action='store_true')
    a=p.parse_args();require(platform.system()=='Linux','Linux required')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    a.solver=a.solver.resolve();a.probe=a.probe.resolve()
    runner=Path(__file__).with_name('run_cls_epawf_pimple.py')
    summary=dict(schema=1,source=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        solver=identity(a.solver),probe=identity(a.probe),final_time=a.final_time,cases={},status='running')
    def save():write_json(a.output/'summary.json',summary)
    def case(name,dt,mode='MassProjection',floor=a.floor,final=a.final_time,*extra):
        output=a.output/name
        cmd=['python3','-B',str(runner),'--solver',str(a.solver),'--output',str(output),'--dt',str(dt),
             '--final-time',str(final),'--diagnostic','--k-bound',mode,'--floor-relative',str(floor),
             '--snapshots','--budget-time',str(a.final_time),*map(str,extra)]
        with (a.output/(name+'.driver.log')).open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
        r=load(output/'case.evolution.json');m=load(output/'manifest.json')
        summary['cases'][name]=dict(command=cmd,status=m['status'],accepted_steps=r['steps'],time=r['time'],failure=r['failure'],
            seconds=m['execution']['seconds'],xml_sha256=sha(output/'case.xml'),bound_integral=sum(x.get('k_bound',{}).get('integral_total',0.)
                for x in m['attempts'] if x['numerically_accepted']),activations=sum(x.get('bound_activations',0) for x in m['attempts'] if x['numerically_accepted']))
        save();return r
    save()
    baseline=[];bounded=[]
    for i,dt in enumerate(a.dts):
        baseline.append(case(f'off-{i}',dt,'Off') if not a.skip_baseline else None)
        bounded.append(case(f'bound-{i}',dt))
    summary['fixed_time_differences']=[]
    for i in range(len(bounded)-1):
        x,y=bounded[i:i+2]
        if abs(x['time']-a.final_time)<1e-13 and abs(y['time']-a.final_time)<1e-13:
            summary['fixed_time_differences'].append(dict(pair=[i,i+1],physical_L2=differences(x,y)))
    summary['status']='fixed_horizon_cases_recorded';save()
    print(json.dumps(summary['cases'],indent=2))


if __name__=='__main__':main()
