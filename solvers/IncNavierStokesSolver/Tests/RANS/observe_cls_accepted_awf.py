#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only accepted-AWF observer with compact incremental correction budgets.

Read immutable completed native blocks. The observer never repairs fields,
advances a native history, retries a rejection, or changes production progress.
"""
import argparse
import datetime
from itertools import combinations
import json
import os
from pathlib import Path
import time

from analyse_cls_accepted_awf import compact,differences,reference_errors,normalized,report,band_key
from continue_cls_epawf_pimple import CampaignLock,atomic,load_progress,process_identity
from run_cls_epawf_picard import load,require,sha


class MetricLedger:
    """Idempotent observer records; distinguish reports from real checkpoints."""
    def __init__(self,path):
        self.path=path;self.records={}
        if not path.exists():return
        lines=path.read_text().splitlines();repair=False
        for i,line in enumerate(lines):
            try:r=json.loads(line)
            except json.JSONDecodeError:
                require(i==len(lines)-1,'corrupt interior observer record');repair=True;break
            n=int(r['global_step'])
            if n in self.records:
                require(self.records[n]['source_identity']==r['source_identity'],'conflicting observer identity')
                repair=True;continue
            self.records[n]=r
        if repair:
            temporary=path.with_suffix('.reconcile.tmp')
            with temporary.open('w') as f:
                for r in self.records.values():f.write(json.dumps(r,separators=(',',':'))+'\n')
                f.flush();os.fsync(f.fileno())
            os.replace(temporary,path)
    def append(self,record):
        n=record['global_step']
        if n in self.records:
            require(self.records[n]['source_identity']==record['source_identity'],'conflicting accepted observer identity')
            return
        with self.path.open('a') as f:
            f.write(json.dumps(record,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno())
        self.records[n]=record


def step_metric(step,attempt,source_identity):
    b=attempt.get('k_bound',{});integrals=b.get('integrals_inner_outer_near_positive_negative',[0.]*5)
    return dict(global_step=int(step),source_identity=source_identity,time=attempt['time'],
        raw_mass_norm=b.get('mass_norm_raw_step'),bound_mass_norm=b.get('mass_norm_correction'),
        accepted_mass_norm=b.get('mass_norm_accepted_step'),raw_bound_inner_product=b.get('mass_inner_product_raw_bound'),
        increment_norm_identity_defect=b.get('mass_norm_increment_identity_defect'),correction_raw_ratio=b.get('correction_to_step_ratio'),
        integral_signed=b.get('integral_total',0.),integral_positive=integrals[3],integral_negative=integrals[4],
        integral_absolute=integrals[3]-integrals[4],bound_activations=b.get('activations',0),
        boundary_floor_defect=b.get('primal_floor_defect'),ground_floor_defect=b.get('ground_floor_defect'),
        bound_seconds=b.get('seconds_total',0.),kkt_solves=b.get('kkt_solves',0),
        inventory_rebuilds=b.get('inventory_rebuilds',0),mass_factor_applications=b.get('mass_factor_applications',0),
        raw_linear_residual=attempt.get('k_frozen_linear_residual'),corrected_frozen_residual=attempt.get('k_corrected_frozen_linear_residual'),
        added_A_delta=attempt.get('k_added_frozen_residual'),storage_alpha_M_delta=attempt.get('k_added_storage_residual'),
        matching_defect=attempt.get('refreshed_trace_defect'),minimum=b.get('corrected_minimum'),ground=b.get('corrected_ground'),
        physical_bands=[{k:r[k] for k in ('name','x_region','x_bounds','lower','upper','correction_signed',
            'correction_positive','correction_negative','correction_absolute') if k in r}
            for r in (b.get('physical_bands') or {}).get('regions',[])],
        solve_counts={k:attempt.get(k,0) for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves',
            'epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')})


def budget(records,first,last,dt):
    selected=[r for n,r in records.items() if first<n<=last];ids={r['global_step'] for r in selected}
    duration=(last-first)*dt;bands={}
    out=dict(first_step=first,last_step=last,duration=duration,complete_coverage=ids==set(range(first+1,last+1)),
        covered_steps=len(selected),signed=sum(r['integral_signed'] for r in selected),
        positive=sum(r['integral_positive'] for r in selected),negative=sum(r['integral_negative'] for r in selected),
        absolute=sum(r['integral_absolute'] for r in selected),activations=sum(r['bound_activations'] for r in selected),
        mass_norm_sum=sum(r['bound_mass_norm'] or 0. for r in selected),
        maximum_correction_raw_ratio=max((r['correction_raw_ratio'] or 0. for r in selected),default=0.),
        maximum_increment_norm_identity_defect=max((abs(r['increment_norm_identity_defect'] or 0.) for r in selected),default=0.),
        bound_wall_seconds=sum(r['bound_seconds'] for r in selected),fixed_bands=bands)
    for r in selected:
        for b in r['physical_bands']:
            name=band_key(b);acc=bands.setdefault(name,dict(lower=b['lower'],upper=b['upper'],
                x_bounds=b.get('x_bounds'),signed=0.,positive=0.,negative=0.,absolute=0.))
            for k in ('signed','positive','negative','absolute'):acc[k]+=b['correction_'+k]
    if duration>0:out['mean_correction_rate']={k:out[k]/duration for k in ('signed','positive','negative','absolute')}
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cases',nargs='+',required=True)
    p.add_argument('--campaigns',nargs='+',type=Path,required=True);p.add_argument('--initial-reports',nargs='+',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--once',action='store_true')
    p.add_argument('--interval',type=float,default=15.);a=p.parse_args()
    require(len(a.cases)==len(a.campaigns)==len(a.initial_reports) and len(set(a.cases))==len(a.cases),'case/initial/campaign inventory differs')
    require(1<=a.interval<=60,'bounded observer period required')
    a.output.mkdir(parents=True,exist_ok=True);lock=CampaignLock(a.output);states={};initial={};metrics={};digests={};cursors={};prefix_done=set();report_paths={}
    observations={};windows={};paths=dict(zip(a.cases,a.campaigns));cursorfile=a.output/'cursors.json'
    if cursorfile.exists():cursors=load(cursorfile)
    previous=load(a.output/'summary.json') if (a.output/'summary.json').exists() else {}
    digests=previous.get('immutable_report_sha256',{})
    paired_history=previous.get('paired',{})
    for case,path in zip(a.cases,a.initial_reports):
        require(Path(case).name==case and case not in ('.','..'),'case identifier must be one filename component')
        d=report(path);require(d['wall_policy']=='CLSCalibratedBudgetAWF','observer requires explicit new wall formulation')
        if not d.get('accepted_awf_initial',{}).get('common_samples'):
            atomic(a.output/'summary.json',dict(status='initial_diagnostic_unavailable',case=case,
                query=d.get('accepted_awf_initial'),production_governing_solves=0,
                interpretation='optional observer stopped; native campaign outcome remains in its own progress record'))
            return
        initial[case]=d['accepted_awf_initial'];states[case]={0:initial[case]};report_paths[case]={}
        prior=previous.get('cases',{}).get(case,{})
        observations[case]={int(n):v for n,v in prior.get('observations',{}).items()};windows[case]=prior.get('fixed_physical_windows',{})
        metrics[case]=MetricLedger(a.output/(case+'-steps.jsonl'));digests[str(path)]=sha(path)
    while True:
        result=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='observing',
            observer=process_identity(os.getpid()),production_governing_solves=0,cases={},paired=paired_history,
            scope='immutable native readbacks; no global polynomial positivity or steady-convergence claim',
            budgets='exact discrete repair integrals and storage diagnostics; physical source/flux estimates remain separate')
        terminal=True
        for case,root in paths.items():
            q=load_progress(root/'progress.json');front=int(q['accepted_steps']);dt=q['dt'];prefix=int(q.get('prefix_steps',0));ledger=metrics[case]
            require(q.get('case_family')=='AcceptedAWFV1','campaign family differs from accepted model')
            terminal=terminal and q['status'] in ('numerically_rejected','targets_completed')
            # Import an immutable short prefix once. If its report covers only
            # a late block, earlier correction budgets remain explicitly unknown.
            if q.get('prefix_report') and prefix and case not in prefix_done:
                prefix_file=Path(q['prefix_report']);prefix_hash=sha(prefix_file)
                require(prefix_hash==q['prefix_report_sha256'],'immutable prefix report changed')
                prefix_report=report(prefix_file);accepted={x[0] for x in prefix_report['stage_audits']}
                for t in prefix_report['pimple_attempts']:
                    if t['time'] not in accepted:continue
                    n=round(t['time']/dt)
                    ledger.append(step_metric(n,t,dict(kind='immutable_native_prefix_report',sha256=prefix_hash,native_step=n)))
                if prefix_report.get('accepted_awf_final',{}).get('common_samples'):
                    n=int(prefix_report['steps']);s=prefix_report['accepted_awf_final'];states[case][n]=s;report_paths[case][n]=prefix_file
                    digests[str(prefix_file)]=prefix_hash
                    if n not in observations[case]:
                        change=differences(s,initial[case]);observations[case][n]=dict(time=s['time'],diagnostic=compact(s),
                            drift_from_actual_h0=change,common_normalized_drift=normalized(change),
                            error_to_own_reference=reference_errors(s),cumulative_correction=budget(ledger.records,0,n,dt))
                prefix_done.add(case)
            path=root/'steps.jsonl';cursor=int(cursors.get(case,0))
            if path.exists():
                with path.open('rb') as f:
                    require(cursor<=path.stat().st_size,'source ledger shrank; inspect before reconciling observer')
                    f.seek(cursor)
                    while True:
                        position=f.tell();line=f.readline()
                        if not line or not line.endswith(b'\n'):break
                        rec=json.loads(line);n=int(rec['global_step'])
                        if n>front:break
                        t=rec['attempt'];require(t['scheduled_complete'] and not t['failure'],'unaccepted source record in native ledger')
                        counts=(t['outer_sweeps'],t['momentum_predictor_solves'],t['pressure_system_solves'],
                            t['epsilon_transport_solves'],t['k_transport_solves'],t['bound_calls'])
                        require(counts==(1,1,1,1,1,1) and t['additional_governing_equation_solves']==0,'changed production solve schedule')
                        ledger.append(step_metric(n,t,dict(kind='native_checkpoint',sha256=rec['checkpoint_sha256'])));cursor=f.tell()
            cursors[case]=cursor
            keep_from=max(0,front-round(.12/dt))
            for seg in q['segments']:
                n=int(seg['end_step']);f=Path(seg['checkpoint']).parent/'case.evolution.json.gz'
                if n>160 and n%32 and n%80 and n not in q['requested_targets']:continue
                report_paths[case][n]=f
                if n in states[case] or (n in observations[case] and n<keep_from) or not f.exists():continue
                try:d=report(f)
                except (OSError,EOFError,json.JSONDecodeError):continue
                if not d.get('accepted_awf_final',{}).get('common_samples'):
                    observations[case][n]=dict(time=d['time'],invalid_optional_query=d.get('accepted_awf_final'))
                    continue
                require(d['steps']==n and abs(d['accepted_awf_final']['time']-d['time'])<1e-13,'completed native block identity differs')
                s=d['accepted_awf_final'];states[case][n]=s;digests[str(f)]=sha(f)
                if n not in observations[case]:
                    change=differences(s,initial[case]);observations[case][n]=dict(time=s['time'],diagnostic=compact(s),
                        drift_from_actual_h0=change,common_normalized_drift=normalized(change),
                        error_to_own_reference=reference_errors(s),
                        cumulative_correction=budget(ledger.records,0,n,dt))
                for duration in (.02,.05,.1):
                    steps=round(duration/dt);old=n-steps
                    if old not in states[case]:continue
                    before=states[case][old];require(abs(s['time']-before['time']-duration)<1e-11,'fixed physical window differs')
                    windows[case][str(old)+'-'+str(n)]=dict(start=before['time'],end=s['time'],duration=duration,
                        dimensional_field_change=differences(s,before),
                        wall_stress_change=s['mean_wall_stress']-before['mean_wall_stress'],
                        exact_correction_budget=budget(ledger.records,old,n,dt))
                # Bound peak memory even when first observing an already long
                # campaign. Completed compact observations/windows remain.
                for old in list(states[case]):
                    if old and old<n-round(.12/dt):del states[case][old]
            oldest=max(0,max(states[case])-round(.12/dt))
            for n in list(states[case]):
                if n and n<oldest:del states[case][n]
            source=initial[case]['common_samples'];bulk=sum(w*u for w,u in zip(source['weights'],source['u']))/sum(source['weights']);flow=6/bulk
            driver=process_identity(q['pid']) if q.get('pid') else None
            expected_native=q.get('native_identity');native=process_identity(expected_native['pid']) if expected_native else None
            supervisor_file=root.parent/'status.json';supervisor=None;supervisor_matches=None
            if supervisor_file.exists():
                sv=load_progress(supervisor_file);expected_supervisor=sv.get('supervisor_identity')
                supervisor=process_identity(expected_supervisor['pid']) if expected_supervisor else None
                supervisor_matches=bool(supervisor and supervisor==expected_supervisor)
            recent=[seg['seconds']/(seg['end_step']-seg['start_step']) for seg in q['segments'][-4:]
                if seg['end_step']>seg['start_step'] and not seg['failure']]
            seconds_per_step=sum(recent)/len(recent) if recent else None
            result['cases'][case]=dict(status=q['status'],accepted_steps=front,accepted_time=q['accepted_time'],
                checkpoint=q['checkpoint'],checkpoint_sha256=q['checkpoint_sha256'],completed_targets=q['completed_targets'],next_target=q.get('next_target'),
                failure=q.get('failure'),case_definition=q.get('case_definition'),
                source_revision=q['source_revision'],frozen_executable_sha256=q['executable_sha256'],
                driver_identity=driver,driver_matches=bool(driver and str(driver['start'])==str(q.get('process_start'))),
                native_identity=native,native_matches=bool(native and native==expected_native),
                supervisor_identity=supervisor,supervisor_matches=supervisor_matches,
                progress_age_seconds=time.time()-q['updated_unix'],recent_wall_seconds_per_accepted_step=seconds_per_step,
                estimated_hours_to_next_target=(q['next_target']-front)*seconds_per_step/3600
                    if seconds_per_step is not None and q.get('next_target') else None,
                cumulative_signed_correction_from_production=q['cumulative_bound_integral'],
                elapsed_wall_times=q['accepted_time']/(initial[case]['H']/.4),elapsed_initial_flow_times=q['accepted_time']/flow,
                wall_time_scale=initial[case]['H']/.4,initial_flow_through_time=flow,
                latest_step=ledger.records[max(ledger.records)] if ledger.records else None,
                observations=observations[case],fixed_physical_windows=windows[case])
        for left,right in combinations(a.cases,2):
            pair=result['paired'].setdefault(left+' vs '+right,{})
            valid=lambda c:{n for n,o in observations[c].items() if not o.get('invalid_optional_query')}
            common=({0}|valid(left))&({0}|valid(right))
            def retained(case,n):
                if n in states[case]:return states[case][n]
                f=report_paths[case].get(n)
                if f is None or not f.exists():return None
                return report(f)['accepted_awf_final']
            for n in sorted(common):
                if str(n) in pair:continue
                astate,bstate=retained(left,n),retained(right,n)
                if astate is None or bstate is None:continue
                require(abs(astate['time']-bstate['time'])<1e-12,'paired observations have different physical times')
                pair[str(n)]=dict(time=astate['time'],dimensional_field_difference=differences(astate,bstate))
        result['immutable_report_sha256']=digests;result['status']='terminal_observed' if terminal else 'observing'
        atomic(cursorfile,cursors);atomic(a.output/'summary.json',result)
        if terminal or a.once:return
        time.sleep(a.interval)


if __name__=='__main__':main()
