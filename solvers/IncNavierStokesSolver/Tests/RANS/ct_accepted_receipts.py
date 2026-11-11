#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Audit consumed inputs and preserve complete small native stage receipts."""
import argparse,hashlib,json,subprocess
from pathlib import Path
from rans_qualification import write_json,exact_vector_hash
from rans_request_manifest import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root;rows={}
    for name in ('A95-dt-0.002','A95-dt-0.001','A93-dt-0.001','E94E95-dt-0.001'):
        folder=r/'pilot-v3'/name;zero=r/'pilot-v3'/(name+'-zero')
        n=json.loads((folder/'case.evolution.json').read_text());z=json.loads((zero/'case.evolution.json').read_text());q=json.loads((zero/'initial.json').read_text())
        b=q['initialization_receipt'];c=n['initialization_receipt']
        assert z['state']==n['accepted_states'][0]==q['observation']['state']
        assert z['signature']==n['signature']
        for x in (b,c):
            assert x['raw_initial_sha256']==x['effective_initial_sha256']==x['fixed_reference_sha256']
            assert x['explicit_boundary_sha256']==x['implicit_boundary_sha256']
            assert x['raw_effective_max_difference']==x['source_target_boundary_max_difference']==0
            assert x['open_explicit_implicit_agreement']
        for k in ('raw_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256',
                  'open_traction_load_sha256','open_diffusive_load_sha256','open_advective_load_sha256','live_fixed_boundary'):
            assert b[k]==c[k]
        assert exact_vector_hash(z['state'])==b['raw_initial_sha256']
        assert all(x.get('bound_calls',0)==0 and x.get('additional_correction_or_restoration_solves',0)==0 for x in n['pimple_attempts'])
        rows[name]=dict(consumed_input_status='pass',descriptor=q['observation']['resolved_descriptor'],
            initial_state_sha256=b['raw_initial_sha256'],explicit_implicit_boundary_sha256=b['explicit_boundary_sha256'],
            fixed_reference_sha256=b['fixed_reference_sha256'],outlet_loads_sha256={k:b[k] for k in ('open_traction_load_sha256','open_diffusive_load_sha256','open_advective_load_sha256')},
            raw_effective_difference=0,source_target_boundary_difference=0,zero_step_vs_consumed_coefficients='bitwise_equal',
            actual_quadrature=q['observation']['actual_quadrature'],basis=q['observation']['field_basis_configuration'],
            report_sha256=sha(folder/'case.evolution.json'),initial_report_sha256=sha(zero/'case.evolution.json'),initial_packet_sha256=sha(zero/'initial.json'),
            signature_sha256=hashlib.sha256(n['signature'].encode()).hexdigest(),
            audit_times=n['audit_times'],state_audits=n['state_audits'],stage_audits=n['stage_audits'],pimple_attempts=n['pimple_attempts'])
    assert rows['A95-dt-0.002']['initial_state_sha256']==rows['A95-dt-0.001']['initial_state_sha256']
    old=json.loads((r/'boundary-controls-v1/closed_old/case.evolution.json').read_text())
    new=json.loads((r/'boundary-controls-v1/closed_new/case.evolution.json').read_text())
    assert old['signature']==new['signature']
    sources={}
    for path in ('solvers/IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.cpp',
                 'solvers/IncNavierStokesSolver/EquationSystems/RANSHistory.cpp'):
        before=subprocess.check_output(['git','show','b127d307:'+path],cwd=r/'source')
        assert hashlib.sha256(before).hexdigest()==sha(r/'source'/path);sources[path]=sha(r/'source'/path)
    result=dict(schema='CTAcceptedConsumedInputStageReceiptsV1',status='pass',cases=rows,
        closed_signature_status='bitwise_equal',closed_signature_sha256=hashlib.sha256(new['signature'].encode()).hexdigest(),
        unchanged_native_history_reader_sources=sources,pilot_qualification='A95 failed; this receipt does not change its status')
    write_json(a.output,result);print('pass');return 0
if __name__=='__main__':raise SystemExit(main())
