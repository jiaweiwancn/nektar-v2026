#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify the selected immutable H9 audit bytes and close its existing ledger."""
import argparse
import json
from pathlib import Path
import re
import time
from run_wj_h9 import load,save,size
from study_wj_h9 import identity,check


def verify(root,publication):
    out=root/'closure';out.mkdir();imports=load(root/'imports.json');verified=[]
    for v in load(publication/'seal.json')['files']:verified.append(check(v))
    for item in imports['reports'].values():
        for key in ('report','xml','original_report'):verified.append(check(item[key]))
    for item in imports['attempts']:
        for key in ('packet','matrix','original_packet','original_matrix','checkpoint'):verified.append(check(item[key]))
    for item in imports['anchors'].values():verified.append(check(item['native']))
    for key in ('inventory','records','local','H8_accounting'):verified.append(check(imports[key]))
    for key in ('ELF','wrapper'):verified.append(check(imports['runtime'][key]))
    for v in imports['runtime']['private_libraries']:verified.append(check(v))
    bindings=load(publication/'phase-bindings.json');expected={('coarse',1),('fine',1),('fine',2)}|{('quarter',i)for i in range(1,5)}
    assert {(v['branch'],v['step'])for v in bindings}==expected and len(bindings)==7
    assert all(v['P_C_bitwise']and v['post_flow_pre_k_bitwise']and v['raw_corrected_bitwise']and v['flow_A_U_pi_bitwise']for v in bindings)
    moments=load(publication/'all-Q1-moments.json')
    assert len(moments)==7
    for v in moments:
        assert sorted(len(x['C'])for x in v['regions'].values())==[34,70]
        for x in v['regions'].values():
            assert max(x[k]for k in ('linear_identity_max','repair_max','post_flow_max'))<2e-15
    telescopes=load(root/'science-v2/telescopes.json')
    for v in telescopes:assert v.get('identity_max',v.get('maximum_identity_defect',float('inf')))<2e-13
    compiled=load(publication/'compiled-checks.json');assert len(compiled)==12
    assert max(v['error']/v['nu_eff']for v in compiled)<2e-11
    assert max(abs(v['implicit_dN_error'])for v in compiled)<2e-11
    gates=load(publication/'component-gates.json');assert gates['status']=='PASS'
    assert all(gates[k]==0 for k in ('new_governing','new_CP0','new_strict','new_builds'))
    logs={};tests=set()
    for name in ('focused-tests-v4','focused-h9-final'):
        log=root/'operations'/name/'command.log';text=log.read_text();assert '\nOK\n'in text
        tests.update(re.findall(r'^(test_\w+) \(([^)]+)\) \.\.\. ok$',text,re.M))
        logs[name]=dict(log=identity(log),request=identity(root/'operations'/name/'request.json'),receipt=identity(root/'operations'/name/'receipt.json'))
    assert len(tests)==51
    matrices=load(publication/'matrix-checks.json');assert len(matrices)==7 and all(v['factorizations']==v['solves']==0 for v in matrices)
    book=load(root/'accounting.json');assert all(book[k]==0 for k in ('governing_attempts','cp0_reservations','strict_readbacks','builds'))
    forbidden={'IncNavierStokesSolver','RANSWallResolutionProbe','cmake','ninja','make'};operations=[]
    for name in book['entries']:
        folder=root/'operations'/name;request=load(folder/'request.json');assert not(forbidden&{Path(x).name for x in request['argv']})
        if book['entries'][name]['status']=='RESERVED':
            assert name=='closure-verification','only this verifier may be running'
            continue
        maps=load(folder/'process-maps.json')
        assert not any(forbidden&{Path(x).name for x in v['command'].split()}for v in maps.values())
        operations.append(dict(label=name,request=identity(folder/'request.json'),receipt=identity(folder/'receipt.json'),log=identity(folder/'command.log'),maps=identity(folder/'process-maps.json')))
    save(out/'verification.json',dict(status='PASS',selected_input_and_publication_hashes=verified,operations=operations,executed_tests=logs,unique_methods=51,methods=[list(v)for v in sorted(tests)],
        governing_attempts=0,new_CP0=0,new_strict=0,builds=0,attempts_audited=7,full_Q1_tests_per_attempt=[34,70],
        completed_batch=dict(k_matvecs=21,k_factorizations=0,k_solves=0,diagnostic_mass_factorizations=52,diagnostic_mass_solves=28,kernel_queries=24,weak_evaluations=24),
        failed_prefix=dict(k_matvecs=3,diagnostic_mass_factorizations=52,diagnostic_mass_solves=4,projected_test_fits=134,kernel_queries=0,
            counting_scope='inferred from preserved exception location after coarse matrix/moment evaluation; no completed attempt published'),
        intermediate_status='diagnostic only; no inherited accepted-state or strict-reader eligibility',
        native_anchors='existing H8 read only; no new native view',loaded_runtime=identity(publication/'runtime.json'),
        no_governing_schur_or_momentum_assembly=True,no_saved_predictor_correction_split=True,source=identity(__file__)))
    print('PASS immutable hashes, 51 unique methods, seven attempts, zero new native work',flush=True)


def finalize(root,publication):
    out=root/'closure';assert load(out/'verification.json')['status']=='PASS';book=load(root/'accounting.json');assert not book.get('closed')
    assert all('seconds'in v and v['status']!='RESERVED'for v in book['entries'].values()),'unknown reservations remain charged; do not close as known'
    actual=sum(v['seconds']for v in book['entries'].values());science=sum(v['seconds']for k,v in book['entries'].items()if k.startswith('science-batch'))
    # Charge full setup and closure reserves in addition to timed child work.
    # These cover parent orchestration, inspections, transfer, source/figure QA,
    # versioned packaging, failure triage and the final scalar receipt work.
    total=book['setup_seconds']+actual+book['closure_reserve'];assert total<=600 and science<=360
    peak=max(v['peak_tree_rss_bytes']for v in book['entries'].values());output=size(root)
    assert peak<=4*1024**3 and output<=2*1024**3
    book.update(closed=True,status='PASS',measured_operation_seconds=actual,charged_seconds=total,science_batch_seconds=science,
        closure_charge_seconds=book['closure_reserve'],closure_charge_scope='full conservative parent/closure reserve charged, not reclaimed',peak_process_tree_RSS=peak,new_root_bytes_before_final_receipts=output)
    save(root/'accounting.json',book);save(out/'accounting-final.json',book)
    save(out/'bridge.json',dict(publication_seal=identity(publication/'seal.json'),verification=identity(out/'verification.json'),accounting=identity(out/'accounting-final.json'),
        original_H8_accounting=check(load(root/'imports.json')['H8_accounting']),runtime=identity(publication/'runtime.json'),
        verifier_receipt=identity(root/'operations/closure-verification/receipt.json'),
        verifier_maps=identity(root/'operations/closure-verification/process-maps.json'),
        source_versions=[identity(p)for directory in [root/'source',*sorted(root.glob('observer-v*'))]for p in sorted(directory.glob('*.py'))],
        superseded=[identity(root/'operations'/name/'receipt.json')for name,v in book['entries'].items()if v['status']!='PASS'],
        scalar_reclassification=dict(original=identity(root/'science-v2/summary.json'),corrected=identity(root/'assessment-v1/summary.json'),reason='physical additive terms distinguished from diagnostic buckets; same load bytes and tolerances'),
        limits=dict(seconds=600,science_seconds=360,RSS=4*1024**3,output=2*1024**3,workers=72),scope='remote-heavy hashes verified on Linux; transport will independently verify compact bytes on Mac'))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()],status='PASS'))
    print(json.dumps(dict(status=book['status'],charged=total,science=science,peak_RSS=peak,new_output=output)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--publication',type=Path,required=True);p.add_argument('--finalize',action='store_true');a=p.parse_args()
    (finalize if a.finalize else verify)(a.root,a.publication)
