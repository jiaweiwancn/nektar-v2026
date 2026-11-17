#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check both cold imports after an actual owned scientific/recovery rotation."""
import argparse
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import history_recovery_check as hr
import physical_cache_history_check as history
import physical_oracle as po
import scientific_accumulator as accumulator
import scientific_retention_check as retention
import terminal_output_check as terminal

require,read = po.require,history.read


def run(manifest,folder,deadline,out):
    require(manifest['schema'] == 'WJH27PostRetentionColdImportsV1','wrong cold-import manifest')
    root = hr.resolve(folder,manifest['root']);assessment_path,assessment_id = history.input_path(manifest['assessment'],root)
    assessment = read(assessment_path)
    require(assessment['status'] == 'PASS_SCIENTIFIC_REDUCTION_AND_ROTATION_ONLY' and assessment['raw_payload_rotation_performed'] is True,
        'scientific rotation is unverified')
    receipt_path,receipt_id = history.input_path(assessment['retention']['receipt'],root);receipt = read(receipt_path)
    require(receipt['status'] == 'PASS_ACTUAL_OWNED_ROTATION','actual owned rotation failed')
    for item in receipt['removed']:require(not Path(item['path']).exists(),'removed owned payload reappeared')
    for item in receipt['retained']:require(hr.identity(Path(item['path'])) == item,'retained owned payload changed')
    for item in receipt['immutable_indices']:require(hr.identity(Path(item['path'])) == item,'original native index changed')
    require(hr.identity(Path(receipt['archive']['path'])) == receipt['archive'],'protected pre-retention archive changed')
    cp,_ = history.input_path(manifest['control'],root);control = read(cp)['state']
    evidence = {}
    require(set(manifest['reloads']) == {'9','12'},'both retained recoveries must be tested')
    for label,spec in manifest['reloads'].items():
        retention.clock(deadline);step = int(label)
        result_path,result_id = history.input_path(spec['result'],root);packet = read(result_path);data = packet['state']
        _,session_id,parameters,solver = terminal.session(spec['session'],root)
        history.configuration(data,control)
        require(packet['restart'] is True and data['step'] == step and data['time'] == step*.005 and
            parameters['NumSteps'] == '0' and parameters['IO_CheckSteps'] == '0' and parameters['WJScientificSteps'] == '0' and
            solver['WJTerminalOutput'] == 'False','zero-step cold-import recipe differs')
        incoming,input_id = history.input_path(spec['input'],root)
        require(input_id in receipt['retained'] and incoming.name == 'epoch3.recovery-'+label+'.bin' and
            hr.resolve(root,solver['WJRestart']) == incoming,'cold import is not the retained recovery')
        states,binary = history.binary_history(incoming,data,.005,deadline,input_id['sha256'])
        require(binary['step'] == step,'retained recovery clock differs')
        exported = history.exported(data,states,reload=True)
        prefix = history.replace_suffix(result_path,'.result.json','');final = Path(str(prefix)+'.final.bin')
        require(final.read_bytes() == incoming.read_bytes(),'cold zero-step serialization changed the full binary history')
        log = Path(str(prefix)+'.steps.jsonl');require(log.read_bytes() == b'','cold zero-step process invented an accepted step')
        require(packet['complete_step_seconds'] == [] and packet['complete_advancing_seconds'] == 0 and
            data['coupled_newton']['refreshes'] == data['coupled_newton']['reuses'] == 0 and data['rebuilds'] == 0,
            'cold zero-step process claims advancement/preconditioner work')
        index_path = Path(str(prefix)+'.science-index.json');index = read(index_path)
        require(index['run_complete'] is True and index['initial_step'] == index['final_step'] == step and index['frames'] == [] and
            index['terminal_output'] is False,'cold zero-step scientific index differs')
        require(all(index[key] == data[key] for key in ('source_identity','reference_identity','signature','dt')),'cold import index identity differs')
        require(packet['scientific_frame_seconds'] == 0 and 0 < packet['scientific_index_seconds'] <= packet['final_export_seconds'],
            'cold zero-step export timing differs')
        evidence[label] = dict(result=result_id,session=session_id,input=binary,output=hr.identity(final),log=hr.identity(log),
            scientific_index=hr.identity(index_path),complete_binary_history_exact=True,accepted_steps=0,
            advancing_residual_evidence_claimed=False,history_roles=list(exported),
            process_seconds=packet['process_seconds'] if 'process_seconds' in packet else None,
            construction_seconds=packet['construction_seconds'],initialisation_seconds=packet['initialisation_seconds'],
            final_export_seconds=packet['final_export_seconds'])
    final_generation,_ = history.input_path(manifest['final_generation'],root);generation,_ = accumulator.read_generation(final_generation)
    require(generation['window_complete'] is True and generation['weight'] == .0364,'retained scientific statistics no longer reload')
    by_kind = {}
    for item in receipt['retained']+receipt['removed']:
        name = Path(item['path']).name
        kind = 'scientific_frame' if '.science-' in name else 'recovery' if '.recovery-' in name else 'final_checkpoint' if name.endswith('.final.bin') else 'statistics_npz'
        by_kind.setdefault(kind,[]).append(item['bytes'])
    files = [p for p in root.rglob('*') if p.is_file()]
    protected = [p for p in files if (root/'protected-retention-v2') in p.parents]
    working = [p for p in files if p.name.startswith(('epoch1.','epoch2.','epoch3.')) or p.parent == final_generation.parent]
    line_bytes = max(len(line)+1 for i in [1,2,3] for line in (root/('epoch'+str(i)+'.steps.jsonl')).read_bytes().splitlines())
    index_bytes = max(item['bytes'] for item in receipt['immutable_indices'])
    batch_steps = 500
    projected_payloads = (batch_steps+2)*max(by_kind['scientific_frame'])+(math.ceil(batch_steps/3)+2)*max(by_kind['recovery'])+3*max(by_kind['final_checkpoint'])+3*max(by_kind['statistics_npz'])
    out['storage'] = dict(actual_run_directory_bytes=sum(p.stat().st_size for p in files),
        protected_retention_evidence_bytes=sum(p.stat().st_size for p in protected),
        retained_selected_payload_bytes=receipt['retained_payload_bytes'],removed_selected_payload_bytes=receipt['removed_payload_bytes'],
        epoch_native_and_V2_generation_working_bytes=sum(p.stat().st_size for p in working),
        maximum_measured_bytes={key:max(values) for key,values in by_kind.items()},
        maximum_measured_step_log_bytes=line_bytes,maximum_four_frame_index_bytes=index_bytes,
        illustrative_batch_steps=batch_steps,illustrative_raw_payload_peak_bytes=projected_payloads,
        illustrative_current_batch_log_index_bytes=batch_steps*line_bytes+math.ceil(batch_steps/4)*index_bytes,
        formula='(N+2)*SCImax+(ceil(N/3)+2)*RECmax+3*FINALmax+3*STATmax; logs/indices separately',
        forecast_scope='Measured 700-cell every-step fixture only; input/XML/FLD/metadata and qualification evidence remain separately charged; no selected production resolution or rigorous future-byte bound',
        V2_generation_metadata='Per-batch provenance plus previous-generation hash; cumulative JSON records retained, older array payloads rotated',
        qualification_archive_is_production_saving=False,production_storage_qualified=False)
    out.update(status='PASS_BOUNDED_SCIENTIFIC_RETENTION_FOUNDATION_ONLY',foundation_gate_pass=True,assessment=assessment_id,
        retention_receipt=receipt_id,reloads=evidence,retained_statistics=hr.identity(final_generation),
        actual_rotation=True,both_retained_recoveries_cold_imported=True,monitoring_escalation_qualified=False,
        selected_recipe_qualified=False,full_H27_qualification=False,H27_D_pass=False,TASK3_READY=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','out'):parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True);args = parser.parse_args();began = time.monotonic()
    require(not args.out.exists(),'post-retention assessment is immutable')
    out = dict(schema='WJH27PostRetentionColdImportCheckV1',status='FAIL',foundation_gate_pass=False,full_H27_qualification=False,
        H27_D_pass=False,TASK3_READY=False,manifest=hr.identity(args.manifest),workers=1,
        sources={key:hr.identity(path) for key,path in dict(checker=__file__,retention=retention.__file__,accumulator=accumulator.__file__,history=history.__file__,JSON_reader=hr.__file__).items()})
    try:
        require(sys.platform.startswith('linux') and sys.byteorder == 'little' and math.isfinite(args.seconds) and 0 < args.seconds <= 60,
            'bounded authoritative Linux check required')
        soft,hard = resource.getrlimit(resource.RLIMIT_AS);cap = min([8*1024**3]+[v for v in (soft,hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS,(cap,hard));run(read(args.manifest),args.manifest.resolve().parent,began+args.seconds,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds'] = time.monotonic()-began;out['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    retention.durable_json(args.out,out);print(json.dumps({key:out[key] for key in ('status','reason','elapsed_seconds') if key in out}))
    if not out['foundation_gate_pass']:raise SystemExit(2)


if __name__ == '__main__':main()
