#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind the fresh CT15 trajectories to complete, unchanged CT14 raw witnesses."""
import argparse,gc,resource
from pathlib import Path
from rans_request_manifest import RequestManifest,sha,canonical
from rans_qualification import exact_vector_hash
from ct_accepted_feasibility import write_json
from ct_accepted_open import require
import ct_accepted8_native as n

def numerical(value):
    """Exclude source paths/schema prose, preserve every numerical tree field."""
    if isinstance(value,dict):
        return {k:numerical(v) for k,v in value.items()
                if not isinstance(v,str) or k.endswith('sha256')}
    if isinstance(value,list):return [numerical(v) for v in value]
    return value

def stages(report):
    result=[]
    for attempt in report['pimple_attempts']:
        if not attempt.get('scheduled_complete'):continue
        row={}
        for kind in ('flow','k'):
            path=Path(attempt[kind+'_snapshot'])
            require(sha(path)==attempt[kind+'_snapshot_sha256'],'changed actual '+kind+' snapshot')
            packet=n.read(path)
            if kind=='k':require(sha(Path(packet['matrix_file']))==packet['matrix_sha256'],'changed actual k matrix')
            else:
                for v in packet.values():
                    if isinstance(v,dict) and 'file' in v and 'sha256' in v:
                        require(sha(Path(v['file']))==v['sha256'],'changed flow stage matrix')
            import hashlib
            row[kind]=hashlib.sha256(canonical(numerical(packet))).hexdigest()
        result.append(row)
    return result

def main():
    p=argparse.ArgumentParser()
    for key in ('ct14','static','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--pair',type=Path);p.add_argument('--static-proof',type=Path)
    p.add_argument('--phase',choices=['static','full'],default='full')
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    require(a.phase=='static' or a.pair is not None,'full comparison needs the fresh pair')
    m=RequestManifest(a.output,dict(schema='CT15CompleteOracleComparisonV1',phase=a.phase,source=sha(Path(__file__)),oracle=sha(a.ct14/'pair-v2/results.json'),candidate=sha(a.pair/'results.json') if a.pair else None,static=sha(a.static/'results.json'),static_proof=sha(a.static_proof/'results.json') if a.static_proof else None))
    result=dict(status='incomplete',static={},members={},criterion='binary64 numerical arrays identical, stronger than the prospective scaled tolerances; source paths/component/schema prose excluded')
    if a.static_proof:
        from ct_accepted9_pressure_run import sealed
        prior=sealed(a.static_proof);request=n.read(a.static_proof/'request-manifest.json')['request']
        require(prior['status']=='pass' and request['static']==sha(a.static/'results.json') and request['oracle']==sha(a.ct14/'pair-v2/results.json'),'foreign/incomplete static comparison')
        result['static']=prior['static']
    for label in (() if a.static_proof else ('target','alpha2000','q48','q64')):
        old=n.read(a.ct14/'static-v1'/label/'frozen.json');new=n.read(a.static/label/'frozen.json')
        require(numerical(old['sparse_backend'])==numerical(new['sparse_backend']),'selected sparse action changed:'+label)
        if label=='target':require(old['test_basis']==new['test_basis'],'selected immutable FE rows changed')
        result['static'][label]=dict(status='pass',all_sparse_backend_numerical_data_equal=True,immutable_FE_rows_equal=label=='target',oracle_sha256=sha(a.ct14/'static-v1'/label/'frozen.json'),candidate_sha256=sha(a.static/label/'frozen.json'))
        del old,new;gc.collect()
    for dt in (() if a.phase=='static' else ('0.00025','0.000125')):
        old=n.read(a.ct14/'pair-v2'/dt/'case.evolution.json');new=n.read(a.pair/dt/'case.evolution.json')
        require(old['audit_times']==new['audit_times'],'accepted coverage changed')
        hashes=[exact_vector_hash(s) for s in old['accepted_states']]
        require(hashes==[exact_vector_hash(s) for s in new['accepted_states']],'actual initial/accepted binary64 coefficients changed')
        before=stages(old);after=stages(new);require(before==after,'actual flow/epsilon/k numerical stage or canonical matrix changed')
        result['members'][dt]=dict(status='pass',times=old['audit_times'],accepted_binary64_hashes=hashes,actual_stage_hashes=after,all_initial_stage_and_accepted_numerical_data_equal=True,maximum_scaled_action_difference=0.,maximum_scaled_state_absolute_pressure_and_traction_difference=0.,oracle_report_sha256=sha(a.ct14/'pair-v2'/dt/'case.evolution.json'),candidate_report_sha256=sha(a.pair/dt/'case.evolution.json'))
        del old,new;gc.collect()
    result['status']='pass';write_json(a.output/'results.json',result);m.seal(a.output,['results.json'])
if __name__=='__main__':main()
