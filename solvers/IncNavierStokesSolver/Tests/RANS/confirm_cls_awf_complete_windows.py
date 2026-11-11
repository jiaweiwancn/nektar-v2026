#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Confirm every steady-window field on native element-resolved physical panels."""
import argparse
from datetime import datetime, timezone
import math, os
from pathlib import Path
from analyse_cls_accepted_awf import normalized
from analyse_cls_epawf_height_comparison import differences
from cls_awf_native_metrics import flat_partition, gl4_differences
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint
from qualify_cls_awf_neutral import observe
from run_cls_epawf_picard import THREADS, load, require, sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('common-panel','probe','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--case',required=True);a=p.parse_args();os.environ.update(THREADS)
    a.output.mkdir(parents=True,exist_ok=False);src=load(a.common_panel/'summary.json');case=src['cases'][a.case];policy=src['policy']
    out=dict(status='running',captured_UTC=datetime.now(timezone.utc).isoformat(),case=a.case,
        source_sha256=sha(a.common_panel/'summary.json'),probe_sha256=sha(a.probe),script_sha256=sha(Path(__file__)),
        policy=policy,observations={},windows=[],governing_solves=0,physical_commits=0,native_history_advanced=False)
    data={};x=None
    for step in case['selected_steps']:
        old=case['observations'][str(int(step))];cp=Path(old['receipt']['checkpoint']);saved=checkpoint(cp)
        xml=Path(old['receipt']['source_report']).parent/'case.xml';xp,zp=flat_partition(xml)
        require(x is None or x==xp,'native tangential partition changed');x=xp
        selection=dict(checkpoint=str(cp),checkpoint_sha256=old['receipt']['checkpoint_sha256'],xml=str(xml),
            native_step=saved['steps'],native_time=saved['time'],source_report=old['receipt']['source_report'],H=.5,balance=False)
        row,q=observe(a.probe,selection,a.output/f'step-{int(step):06d}');data[int(step)]=q
        require(all(z in q['common_reporting_knots'] for z in zp),'native vertical interface omitted')
        require(abs(math.fsum(q['common_samples']['epsilon_weights'])-21)<1e-11,'full epsilon area differs')
        out['observations'][str(int(step))]=row;atomic(a.output/'summary.json',out)
    for first,last,old in zip(case['selected_steps'],case['selected_steps'][1:],case['windows']):
        left,right=data[int(first)],data[int(last)];dim=differences(right,left);norm=normalized(dim)
        independent=gl4_differences(right,left,x)
        errors={k:abs(independent['dimensional'][k]-dim[k]) for k in norm}
        require(max(errors.values())<2e-12,'GL4/GL16 physical norm disagreement')
        field_pass=all(norm[k]<=v for k,v in policy['normalized_window_L2_change_maximum'].items())
        other=all(old[k] for k in ('residuals_pass','matching_pass','traction_pass','bound_supply_pass'))
        out['windows'].append(dict(start=left['time'],end=right['time'],old_normalized_field_change=old['normalized_field_change'],
            confirmed_normalized_field_change=norm,dimensional_field_change=dim,independent_gl4=independent,
            gl4_gl16_absolute_error=errors,fields_pass=field_pass,nonfield_gates_pass=other,
            complete_screen_pass=field_pass and other))
    count=sum(w['complete_screen_pass'] for w in out['windows'])
    out.update(status='passed',passing_windows=count,classification='candidate_steady_selected_operator_at_declared_screen' if count==3 else 'not_steady_by_declared_screen',completed_UTC=datetime.now(timezone.utc).isoformat())
    atomic(a.output/'summary.json',out);print(a.output/'summary.json')


if __name__=='__main__':main()
