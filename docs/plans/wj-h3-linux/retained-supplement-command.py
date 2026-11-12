#!/usr/bin/env python3
"""Additional passive base-T=.08 integration and unapplied tensor diagnostics."""
from pathlib import Path
import sys
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h3 import load,save,sha
from wj_h3_analysis import actual,identity,quadrature,state_context,principal
from validate_wj1 import Probe

root=Path(sys.argv[1]);out=root/'retained-supplement-v1';out.mkdir(exist_ok=False)
h=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111')
progress=load(h/'study-v1/progress.json');record=progress['runs']['base-pilot-on']
rp=actual(record['report']);assert sha(rp)==record['report_sha256'];r=load(rp)
s,hist,a,t,p=state_context(r);q=quadrature(r,s,hist,a)
idx=load(h/'observations-v3/index.json');entry=idx['base08-nativeweak'];f=actual(entry['native'])
assert sha(f)==entry['native_sha256'];native=load(f)['weak_checks']
assert_allclose(q['views'][8]['components']['outer_momentum_stress'],native['native_terms']['outer_momentum_stress'],rtol=2e-11,atol=2e-13)
assert q['smallest_adequate']==12 and q['reference_status']=='PASS'
original_errors={}
for key,views in [('base08',q['views']),('refined-fine08',load(root/'retained-v1/half08.json')['quadrature']['views'])]:
    views={int(k):v for k,v in views.items()}
    err=float(np.max(abs(np.array(views[8]['components']['outer_momentum_stress'])-views[24]['components']['outer_momentum_stress'])))
    old=idx[key+'-q8']['independent_weak']
    expected=old['component_max_errors']['outer_momentum_stress']
    assert abs(err-expected)<2e-13
    original_errors[key]=dict(original_q8_vs_q24_error=expected,independently_reconstructed_error=err,
        original_regional_gate=old['conditions']['regional_sum'],original_status=old['status'])
rr=load(actual(progress['runs']['refined-pilot']['report']));w=rr['physical_range_events'][-1]['wj_packet']['principal_witnesses']['in_plane']
G=np.array(w['raw_G_row_major']).reshape(3,3);shear=np.zeros((3,3));shear[0,2]=G[0,2]
kernel=actual('/media/jay/新加卷/nektar-v2026-runs/wj-h2-20261111/build-v1/WJKernelProbe');probe=Probe(str(kernel));diagnostics={}
for name,g in [('pure_shear_diagnostic_only',shear),('planar_tracefree_diagnostic_only',G-np.trace(G)*np.diag([1.,0,1.])/2)]:
    v=principal(w['k'],w['epsilon'],g);direction=np.outer(v['polarization'],v['xi'])
    c=probe.closure(w['k'],w['epsilon'],g);d=probe.closure(w['k'],w['epsilon'],g,(0,0,direction))
    form=float(-np.sum(direction*d['R'])/c['scalar'][1]);assert_allclose(form,v['in_plane'],rtol=1e-9,atol=2e-11)
    diagnostics[name]=dict(G=g.tolist(),principal=v,compiled_form=form,applied_to_state=False)
probe.close();assert diagnostics['pure_shear_diagnostic_only']['compiled_form']<0
save(out/'base08-full-views.json',q)
result=dict(status='PASS',role='post-control passive confirmation; fixed q_star12 unchanged',producer=identity(rp),
    callback_history_sha256=p['stage_history_sha256'],same_node_native_packet=identity(f),
    q_star=q['smallest_adequate'],reference_status=q['reference_status'],comparisons=q['comparisons'],
    original_H2_errors=original_errors,unapplied_counterfactuals=diagnostics,governing_solves=0,
    full_views=identity(out/'base08-full-views.json'),command_sha256=sha(__file__))
save(out/'summary.json',result);print(result['original_H2_errors'],diagnostics,flush=True)
