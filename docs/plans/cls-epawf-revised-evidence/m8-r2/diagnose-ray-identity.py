#!/usr/bin/env python3
"""Diagnose one retained terminal moving-ray identity failure; never repair it.

Exactly one Analyzer.domain call is allowed. Only its exact geometry-identity
ValueError is inspected, through the pinned domain traceback frame. Pairing
unmatched coordinates is diagnostic metadata: the original gate stays false.
No weak residual, Jacobian, endpoint, retry or native executable is evaluated.
"""
from collections import Counter, defaultdict
from decimal import Decimal, localcontext
import importlib.util
import json
from pathlib import Path
import platform
import sys
import argparse
import os
import time
import traceback

FAILURE = 'Native/independent moving-ray inventory geometry differs'
FAILED_REPORT_SHA = '9208cf90b8d6bd7f11f0f46d117c5e9b9d79cac640f8518fecf49855465d0aae'
ORIGINAL_COORDINATE_TOLERANCE = 2e-11
ORIGINAL_WEIGHT_TOLERANCE = 2e-11


def metadata_module():
    path = Path(__file__).with_name('analyze-terminal-record.py')
    spec = importlib.util.spec_from_file_location('r2_terminal_metadata',path)
    module = importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module,path


def descriptor(row):
    return (row['inventory'],int(row['geometry']),int(row['block']),int(row['side']),
            bool(row['strict']),bool(row['dependent_trace']))


def take_difference(rows, excess, key):
    remaining=excess.copy();result=[]
    for index,row in enumerate(rows):
        identity=key(row)
        if remaining[identity]>0:
            result.append((index,row));remaining[identity]-=1
    if any(remaining.values()):
        raise ValueError('Counter excess could not be extracted without altering original rows')
    return result


def scaled_errors(left,right,np):
    left,right=np.asarray(left,float),np.asarray(right,float)
    if left.size==0:
        return dict(count=0,maximum_absolute_difference=0.,maximum_componentwise_scaled_difference=0.)
    if left.shape!=right.shape or not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise ValueError('Invalid retained diagnostic comparison arrays')
    return dict(count=int(left.shape[0]),maximum_absolute_difference=float(np.max(np.abs(left-right))),
        maximum_componentwise_scaled_difference=float(np.max(np.abs(left-right)/np.maximum(1.,np.abs(right)))))


def bin_witness(a,b,axis):
    import math
    x,y=float(a['coordinate'][axis]),float(b['coordinate'][axis])
    rx,ry=round(x,12),round(y,12)
    if rx==ry:
        return None
    dx,dy=Decimal.from_float(x),Decimal.from_float(y)
    with localcontext() as context:
        context.prec=200
        bx,by=Decimal(str(rx)),Decimal(str(ry));mid=(bx+by)/2
        native_distance=str(abs(dx-mid));independent_distance=str(abs(dy-mid))
    ulp=max(math.ulp(x),math.ulp(y))
    return dict(axis=axis,native=x,independent=y,native_hex=x.hex(),independent_hex=y.hex(),
        native_round12=rx,independent_round12=ry,absolute_difference=abs(x-y),
        difference_in_maximum_endpoint_ulps=abs(x-y)/ulp,
        adjacent_decimal_bins=abs(bx-by)==Decimal('1e-12'),
        exact_native_decimal=str(dx),exact_independent_decimal=str(dy),decimal_bin_midpoint=str(mid),
        exact_values_straddle_midpoint=min(dx,dy)<=mid<=max(dx,dy),
        native_distance_to_midpoint=native_distance,independent_distance_to_midpoint=independent_distance)


def pair_unmatched(native,independent,point_key,np):
    """Minimum total L-infinity coordinate distance within exact descriptors."""
    from scipy.optimize import linear_sum_assignment
    left,right=defaultdict(list),defaultdict(list)
    for item in native:left[descriptor(item[1])].append(item)
    for item in independent:right[descriptor(item[1])].append(item)
    pairs,remaining,groups=[],[],[]
    for identity in sorted(set(left)|set(right)):
        a,b=left[identity],right[identity]
        if not a or not b:
            remaining.extend(dict(side='native' if a else 'independent',sorted_index=index,row=row)
                for index,row in (a or b))
            groups.append(dict(descriptor=list(identity),native_count=len(a),independent_count=len(b),paired=0,
                complete_pairing=False,unique_nearest_assignment_proven=False))
            continue
        if len(a)*len(b)>4_000_000:
            remaining.extend(dict(side=side,sorted_index=index,row=row)
                for side,values in (('native',a),('independent',b)) for index,row in values)
            groups.append(dict(descriptor=list(identity),native_count=len(a),independent_count=len(b),paired=0,
                complete_pairing=False,unique_nearest_assignment_proven=False,
                pairing_skipped_for_metadata_memory_cap=True,cost_matrix_entry_limit=4_000_000))
            continue
        ca=np.asarray([row['coordinate'] for _,row in a],float)
        cb=np.asarray([row['coordinate'] for _,row in b],float)
        distance=np.max(np.abs(ca[:,None,:]-cb[None,:,:]),axis=2)
        ai,bi=linear_sum_assignment(distance)
        used_a,used_b=set(map(int,ai)),set(map(int,bi))
        unique=len(a)==len(b)
        group_pairs=[]
        for i,j in zip(ai,bi):
            i,j=int(i),int(j);cost=float(distance[i,j]);row_min=float(np.min(distance[i]));column_min=float(np.min(distance[:,j]))
            row_ties=int(np.count_nonzero(distance[i]==row_min));column_ties=int(np.count_nonzero(distance[:,j]==column_min))
            nearest=cost==row_min and cost==column_min
            unique=unique and nearest and row_ties==column_ties==1
            coordinate=scaled_errors([a[i][1]['coordinate']],[b[j][1]['coordinate']],np)
            weight=scaled_errors([a[i][1]['weight']],[b[j][1]['weight']],np)
            alternatives=np.delete(distance[i],j)
            observations=dict(descriptor=list(identity),native_sorted_index=a[i][0],independent_sorted_index=b[j][0],
                native_point_key=list(point_key(a[i][1])),independent_point_key=list(point_key(b[j][1])),
                native_row=a[i][1],independent_row=b[j][1],assigned_maximum_coordinate_distance=cost,
                coordinate_difference=coordinate,weight_difference=weight,
                assigned_is_both_row_and_column_nearest=nearest,row_minimum_distance=row_min,column_minimum_distance=column_min,
                equal_row_minimum_candidate_count=row_ties,equal_column_minimum_candidate_count=column_ties,
                minimum_nonassigned_row_distance=float(np.min(alternatives)) if len(alternatives) else None,
                local_equal_minimum_ambiguity=row_ties>1 or column_ties>1,
                original_coordinate_tolerance=ORIGINAL_COORDINATE_TOLERANCE,original_weight_tolerance=ORIGINAL_WEIGHT_TOLERANCE,
                bin_witnesses=[item for axis in range(3) if (item:=bin_witness(a[i][1],b[j][1],axis)) is not None])
            pairs.append(observations);group_pairs.append(len(pairs)-1)
        for side,values,used in (('native',a,used_a),('independent',b,used_b)):
            remaining.extend(dict(side=side,sorted_index=index,row=row) for k,(index,row) in enumerate(values) if k not in used)
        groups.append(dict(descriptor=list(identity),native_count=len(a),independent_count=len(b),paired=len(ai),
            pair_indices=group_pairs,complete_pairing=len(a)==len(b),unique_nearest_assignment_proven=bool(unique),
            uniqueness_scope='Sufficient strict row/column-nearest criterion; unresolved ties are retained'))
    coordinates=scaled_errors([p['native_row']['coordinate'] for p in pairs],[p['independent_row']['coordinate'] for p in pairs],np)
    weights=scaled_errors([p['native_row']['weight'] for p in pairs],[p['independent_row']['weight'] for p in pairs],np)
    return dict(method='Minimum total L-infinity coordinate distance, one-to-one within the exact six-field structural descriptor; no weights/values used to choose pairs',
        paired_count=len(pairs),remaining_unmatched=remaining,remaining_unmatched_count=len(remaining),groups=groups,pairs=pairs,
        paired_coordinate_errors=coordinates,paired_weight_errors=weights,
        coordinate_difference_exceeds_original_tolerance=coordinates['maximum_componentwise_scaled_difference']>ORIGINAL_COORDINATE_TOLERANCE,
        weight_difference_exceeds_original_tolerance=weights['maximum_componentwise_scaled_difference']>ORIGINAL_WEIGHT_TOLERANCE,
        ambiguity_present=any(p['local_equal_minimum_ambiguity'] for p in pairs),
        all_groups_have_unique_nearest_assignment_proven=all(g['unique_nearest_assignment_proven'] for g in groups),
        uniqueness_scope='Only the extracted Counter excess. Exact-key duplicate call identities and original call ordering are not proved by this diagnostic assignment',
        acceptance_gate_changed=False,original_geometry_identity_gate=False,diagnostic_pairing_is_not_qualification=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','preflight','preflight-session','stage-session','trace','runtime-pin',
                 'failure-report','failure-receipt','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--expected-commit',required=True)
    parser.add_argument('--record-index',type=int,default=0)
    parser.add_argument('--cpu',type=int,default=49)
    args=parser.parse_args();T,terminal_path=metadata_module()
    T.require(platform.system()=='Linux','Retained ray-identity diagnosis belongs to authoritative Linux')
    args.source=args.source.resolve();args.output.mkdir(parents=True,exist_ok=False);started=time.monotonic()
    helper=args.source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_stage_analysis.py'
    receipt=dict(schema=1,status='preparing',success=False,native_calls=0,weak_calls=0,domain_calls=0,
        Jacobian_or_endpoint_calls=0,started_utc=T.stamp(),command_argv=sys.argv,wrapper_sha256=T.sha(__file__),
        source_commit=T.git(args.source,'rev-parse','HEAD'),source_clean_start=not bool(T.git(args.source,'status','--porcelain')),
        cpu_affinity=sorted(os.sched_getaffinity(0)),selected_cpu=args.cpu,
        numerical_threads={key:os.environ.get(key) for key in T.THREADS},declared_external_capacity_seconds=900,
        inputs={},failures=[],original_geometry_identity_gate=False,original_failed_analysis_promoted=False,
        scope='Terminal record only: diagnose exact rounded-key identity failure using one original domain call; no gate change, weak/J/native call, path qualification or physical conclusion')
    def bind(path,expected=None):
        item=T.file_identity(path)
        T.require(expected is None or item['sha256']==expected,'Changed bound input: '+item['path'])
        old=receipt['inputs'].setdefault(item['path'],item)
        T.require(old==item,'Input changed during diagnosis: '+item['path'])
    def checkpoint(name,result):
        T.save(args.output/(name+'.json'),result)
        receipt['latest_completed_checkpoint']=T.file_identity(args.output/(name+'.json'))
        T.save(args.output/'receipt.json',receipt)
    T.save(args.output/'receipt.json',receipt)
    try:
        T.require(receipt['source_commit']==args.expected_commit and receipt['source_clean_start'],'Unpinned driver source')
        T.require(args.cpu>=0 and receipt['cpu_affinity']==[args.cpu] and all(v=='1' for v in receipt['numerical_threads'].values()),'Wrong CPU/thread allocation')
        for path in (args.preflight,args.preflight_session,args.stage_session,args.runtime_pin,args.failure_report,args.failure_receipt,helper,terminal_path,Path(__file__)):bind(path)
        T.require(T.sha(args.failure_report)==FAILED_REPORT_SHA,'Require the exact retained failed full auxiliary analysis')
        failed=json.loads(args.failure_report.read_text());old=json.loads(args.failure_receipt.read_text())
        T.require(failed['analysis_success'] is False and old['success'] is False and old['producing_source_commit']==T.PRODUCER,
            'Historical failed gate/provenance was not retained')
        record,prefix=T.snapshot_record(args.trace,args.record_index,args.output);bind(prefix['snapshot'],prefix['selected_sha256'])
        receipt['active_trace_prefix']=prefix
        retained=old['active_trace_prefix']
        T.require(prefix['selected_sha256']==retained['selected_sha256'] and prefix['prefix_sha256']==retained['prefix_sha256'],
            'Selected completed record is not byte-identical to the retained failed analysis')
        pin=json.loads(args.runtime_pin.read_text());producing=Path(pin['source_root'])
        T.require(pin['success'] and pin['clean'] and pin['source_commit']==T.PRODUCER and
            T.git(producing,'rev-parse','HEAD')==T.PRODUCER and not T.git(producing,'status','--porcelain'),'Producing source/runtime differs from original9f20')
        key='build/IncNavierStokesSolver';bind(pin['runtime_resolved_paths'][key],pin['runtime_artifacts'][key])
        receipt.update(producing_source_commit=T.PRODUCER,producing_solver=T.file_identity(pin['runtime_resolved_paths'][key]),
            producing_runtime_sections={key:T.canonical(pin[key]) for key in T.RUNTIME_KEYS})
        T.require(receipt['producing_runtime_sections']==old['producing_runtime_sections'],'Producing fingerprint sections differ from the retained failure')
        sys.dont_write_bytecode=True;sys.path.insert(0,str(helper.parent))
        import numpy as np
        import cls_epawf_m8_r2_stage_analysis as analysis
        import cls_epawf_m8_r2_reuse as reuse
        historical_source=reuse.option(old,'source');historical_helper=historical_source/reuse.HELPER
        T.require(T.git(historical_source,'rev-parse','HEAD')==old['source_commit'] and
            not T.git(historical_source,'status','--porcelain'),'Historical failed-analysis source pin changed')
        expected=old['imported_source_files_start'][str(historical_helper)]
        bind(historical_helper,expected)
        T.require(reuse.scientific_core(historical_helper)==reuse.scientific_core(helper),'Original/current scientific AST differs')
        for path,module in [(Path(m.__file__).resolve(),m) for m in list(sys.modules.values()) if getattr(m,'__file__',None)]:
            if path.is_relative_to(args.source):
                bind(path)
                relative=str(path.relative_to(args.source))
                if relative.startswith('solvers/') and relative not in (reuse.HELPER,reuse.METADATA):
                    prior=historical_source/relative
                    T.require(str(prior) in old['imported_source_files_start'] and T.sha(path)==old['imported_source_files_start'][str(prior)],
                        'Scientific import differs from retained failure: '+relative)
        raw=json.loads(args.preflight.read_text());binding=T.validate_record(raw,record,args.preflight_session,args.stage_session,analysis)
        analyzer=analysis.Analyzer(dict(raw,stage_test_target_time=record['time']),args.stage_session,args.source)
        receipt.update(status='running',domain_call_started_utc=T.stamp(),domain_calls=1)
        checkpoint('input-binding',dict(binding=binding,static_coverage=analyzer.static_coverage,record_snapshot=prefix,
            original_failed_report_sha256=FAILED_REPORT_SHA,scientific_AST_sha256=reuse.scientific_core(helper)))
        try:
            analyzer.domain(np.asarray(record['state'],float),record['domain'],retain=True)
        except ValueError as error:
            if str(error)!=FAILURE:
                raise
            frames=[];tb=error.__traceback__
            while tb is not None:
                if tb.tb_frame.f_code.co_name=='domain' and Path(tb.tb_frame.f_code.co_filename).resolve()==helper.resolve():frames.append(tb)
                tb=tb.tb_next
            T.require(len(frames)==1 and frames[0].tb_frame.f_locals.get('self') is analyzer,'Failure does not belong to the expected original domain frame')
            local=frames[0].tb_frame.f_locals
            a,b=list(local['a']),list(local['b']);partial=dict(local['result'])
            T.require(local['native'] is record['domain'] and partial['admissible'] is True,'Wrong failure point or pre-identity domain state')
            T.require(Counter(map(id,local['records']))==Counter(map(id,b)),
                'Captured sorted independent rows do not preserve the actual domain records')
            receipt.update(expected_failure_reproduced=True,domain_call_finished_utc=T.stamp(),
                caught_failure=dict(type=type(error).__name__,message=str(error),filename=frames[0].tb_frame.f_code.co_filename,
                    line=frames[0].tb_lineno,traceback=traceback.format_exc()))
            del frames,tb,local
        else:
            raise ValueError('Exact historical geometry-identity failure was not reproduced; no gate is promoted')
        T.require(len(a)==len(b),'Earlier passed common inventory length check no longer holds')
        ca,cb=Counter(map(analyzer.point_key,a)),Counter(map(analyzer.point_key,b))
        excess_a,excess_b=ca-cb,cb-ca
        T.require(bool(excess_a or excess_b),'Historical rounded-key failure reproduced with no multiset differences; sorted ordering requires separate review')
        unmatched_a=take_difference(a,excess_a,analyzer.point_key);unmatched_b=take_difference(b,excess_b,analyzer.point_key)
        all_coordinates=scaled_errors([row['coordinate'] for row in a],[row['coordinate'] for row in b],np)
        all_weights=scaled_errors([row['weight'] for row in a],[row['weight'] for row in b],np)
        minimum_keys=('minimum_dimensionless_margin','minimum_inner_k','minimum_outer_k','minimum_outer_epsilon',
            'minimum_sampled_K','minimum_projected_trace_epsilon','minimum_b','minimum_lambda','cutoff_range')
        minima={key:dict(independent_pre_identity=partial.get(key),native=record['domain'].get(key)) for key in minimum_keys}
        unmatched=dict(original_geometry_identity_gate=False,failure=receipt['caught_failure'],native_count=len(a),independent_count=len(b),
            exact_structural_fields=['inventory','geometry','block','side','strict','dependent_trace'],
            counter_difference_native=[dict(point_key=list(key),count=count) for key,count in sorted(excess_a.items())],
            counter_difference_independent=[dict(point_key=list(key),count=count) for key,count in sorted(excess_b.items())],
            unmatched_native_count=sum(excess_a.values()),unmatched_independent_count=sum(excess_b.values()),
            matched_rounded_key_count=sum((ca&cb).values()),original_sorted_order_coordinate_errors=all_coordinates,
            original_sorted_order_weight_errors=all_weights,original_coordinate_tolerance=ORIGINAL_COORDINATE_TOLERANCE,
            original_weight_tolerance=ORIGINAL_WEIGHT_TOLERANCE,minimum_scalars_already_computed=minima,
            prior_comparison_errors=dict(analyzer.audit.errors),
            prior_checks_completed=analyzer.audit.checks,
            coordinate_and_weight_near_calls_after_identity_gate_were_not_reached=True)
        checkpoint('unmatched-observations',unmatched)
        import scipy
        from scipy.optimize import _lsap
        bind(_lsap.__file__)
        receipt['assignment_runtime']=dict(scipy_version=scipy.__version__,
            numpy_version=np.__version__,assignment_image=T.file_identity(_lsap.__file__))
        pairing=pair_unmatched(unmatched_a,unmatched_b,analyzer.point_key,np)
        result=dict(schema=1,success=True,diagnosis_completed=True,original_geometry_identity_gate=False,
            original_failed_full_analysis_remains_failed=True,minimum_scalars_already_computed=minima,
            unmatched=unmatched,pairing=pairing,domain_calls=1,weak_calls=0,Jacobian_or_endpoint_calls=0,native_calls=0,
            interpretation='Counter and one-to-one coordinate comparisons diagnose the retained rounded-key rejection; no original threshold or identity gate is changed',
            no_physical_or_path_qualification=True)
        checkpoint('validation',result);receipt['success']=True
    except Exception as error:
        receipt['failures'].append(dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
    try:
        for item in receipt['inputs'].values():T.require(T.file_identity(item['path'])==item,'Changed input: '+item['path'])
        if 'active_trace_prefix' in receipt:
            prefix=receipt['active_trace_prefix'];prefix['prefix_sha256_end']=T.prefix_sha(args.trace,prefix['prefix_bytes'])
            T.require(prefix['prefix_sha256_end']==prefix['prefix_sha256'],'Analyzed complete native prefix changed')
        receipt.update(source_commit_end=T.git(args.source,'rev-parse','HEAD'),source_clean_end=not bool(T.git(args.source,'status','--porcelain')))
        T.require(receipt['source_commit_end']==receipt['source_commit'] and receipt['source_clean_end'],'Driver source changed')
    except Exception as error:
        receipt['success']=False;receipt['failures'].append(dict(type='ProvenanceError',message=str(error)))
    if receipt['failures']:receipt['success']=False
    receipt.update(status='completed' if receipt['success'] else 'failed',finished_utc=T.stamp(),seconds=time.monotonic()-started)
    T.save(args.output/'receipt.json',receipt)
    print(json.dumps({key:receipt.get(key) for key in ('success','status','domain_calls','seconds','failures')}),flush=True)
    return 0 if receipt['success'] else 1


if __name__=='__main__':
    raise SystemExit(main())
