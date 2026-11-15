#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small primitive tables from sealed CT11 evidence; no new native run."""
import argparse,json,os
from pathlib import Path
from collections import defaultdict
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted9_pressure_run import sealed
from ct_accepted_feasibility import write_json,require
from rans_request_manifest import RequestManifest,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('qualified','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();q=a.qualified;sealed(q,'decision.json');read=lambda x:json.loads((q/x).read_text())
    m=RequestManifest(a.output,dict(schema='CT11PrimitiveTablesV1',qualified_seal=sha(q/'completed-artifacts.json'),source=sha(Path(__file__)),scope='derived tables of sealed results and hash-checked exclusive timing traces only; no new native trajectory'))
    eq=read('equivalence.json');cap=read('capacity.json');controls=read('controls.json');target=read('target.json');cost=read('cost.json');out=dict(schema='CT11PrimitiveTablesV1',observer=dict(cases=len(eq['checks']),controls=sum(len(v) for v in eq['checks'].values()),maximum_scaled_primitive_difference=max(r['all_primitive_difference']['maximum_scaled_difference'] for v in eq['checks'].values() for r in v.values()),maximum_scaled_reconstruction=max(x['scaled_defect'] for v in eq['checks'].values() for r in v.values() for x in r['classification']['reconstruction'].values()),maximum_scaled_advective_IBP=max(x['scaled_defect'] for v in controls['postchecks']['cases'].values() for r in v.values() for x in r['advective_IBP'].values()),prescribed_epsilon_storage_signal=controls['postchecks']['actual_prescribed_epsilon_storage_signal']),capacity={},target={},temporal=target['temporal'],traction_temporal=target['traction_temporal'],regional_temporal_pressure=target['regional_temporal_pressure'],exclusive_accepted_cost={},resources=read('decision.json')['resources'])
    for k,r in cap['controls'].items():out['capacity'][k]=dict(execution=r['execution'],output_bytes=r['output_bytes'],row_counts={n:r['volume_observer'][n] for n in ('inner_points','outer_points')},maximum_scaled_reconstruction=max(x['scaled_defect'] for x in r['frozen_load_reconstruction'].values()))
    for k,r in target['members'].items():
        last=r['accepted_metrics'][-1];out['target'][k]=dict(time=r['time'],native_CPU=r['execution']['process_cpu_seconds'],native_wall=r['execution']['wall_seconds'],peak_RSS_KiB=r['execution']['peak_RSS_KiB'],initialization_CPU=r['zero_step_execution']['process_cpu_seconds'],drift={n:v['scaled'] for n,v in last['drift'].items()},absolute_reference_error={n:v['scaled'] for n,v in last['absolute_reference_error'].items()},ground=last['traction']['normalized_wall_traction_error'],interface=last['traction']['normalized_interface_traction_RMS'],mass=last['mass']['scaled'],budgets={n:{e:v['scaled_frozen_defect'] for e,v in b['equations'].items()} for n,b in r['budgets']['budgets'].items()},all_audited_minima=r['domain']['native_audited_extrema']['over_all_audited_states_and_stages'],independent_inventory=r['domain']['fixed'],stage_linear_checks=r['actual_stage_checks'],regional_pressure=r['regional_pressure'],corrections=r['corrections'],extra_governing_solves=r['extra_solves'])
        f=Path(cost[k]['timing_path']);require(sha(f)==cost[k]['timing_sha256'],'changed accepted timing trace');names=defaultdict(lambda:dict(calls=0,CPU=0.,wall=0.));steps={}
        for line in f.read_text().splitlines():
            t=json.loads(line)
            if t['step']<=0:continue
            v=names[t['name']];v['calls']+=t.get('calls',1);v['CPU']+=t['exclusive_process_cpu_seconds'];v['wall']+=t['exclusive_wall_seconds'];steps[t['step']]=steps.get(t['step'],0.)+t['exclusive_process_cpu_seconds']
        require(len(steps)==len(r['accepted_metrics'])-1,'accepted timing count');out['exclusive_accepted_cost'][k]=dict(steps=steps,categories=dict(names),native_nonstep_CPU=r['execution']['process_cpu_seconds']-sum(steps.values()),scope='exclusive accepted-step context only; initialization/output outside context separate; no nested inclusive totals')
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print(json.dumps(dict(observer=out['observer'],pressure=out['temporal']['Pi'],capacity={k:{n:v['execution'][n] for n in ('process_cpu_seconds','wall_seconds','peak_RSS_KiB','peak_virtual_KiB')} for k,v in out['capacity'].items()},target={k:{n:v[n] for n in ('drift','ground','interface','mass','native_CPU','peak_RSS_KiB')} for k,v in out['target'].items()}),indent=2),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
