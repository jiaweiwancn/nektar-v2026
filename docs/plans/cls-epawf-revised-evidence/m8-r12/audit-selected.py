#!/usr/bin/env python3
"""Audit exactly the selector-chosen R11 trial, reusing an identical complete audit."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time

sys.dont_write_bytecode=True
THREADS=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
         'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--policy',choices=('L1','Lf','N1','Nf'),required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--comparison',type=Path)
    a=p.parse_args()
    assert platform.system()=='Linux' and len(os.sched_getaffinity(0))==1
    assert all(os.environ.get(k)=='1' for k in THREADS)
    e=a.epoch.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r4_audit')
    driver=module(evidence/'m8-r11/experiment.py','r11_selected')
    comparison=a.comparison or e/'comparison.json'
    C=h.load(comparison)
    records=[x for x in C['events'] if x['kind']=='r8_trial']
    policy=next(x for x in C['policies'] if x['name']==a.policy)
    chosen=driver.selected_trial(policy,records)
    started=time.monotonic()
    if chosen is None:
        h.save(a.output,dict(policy=a.policy,selected=None,status='no_selected_native_candidate',new_native_trials=0))
        return
    # Reuse an independently completed R12 audit only for an identical tuple.
    # This also covers Nf selecting its unchanged, newly audited Lf seed.
    for prior_path in sorted((e/'audits').glob('*.json')):
        if prior_path.resolve()==a.output.resolve():continue
        prior=h.load(prior_path)
        if prior.get('status')!='independently_audited':continue
        if (all(prior.get(k,{}).get('sha256')==chosen[k]['sha256'] for k in ('state','residual')) and
                prior.get('history',{}).get('sha256')==C['initial']['history']['sha256']):
            assert C['native_history_verified']
            h.save(a.output,dict(policy=a.policy,status='reused_exact_R12_audit',state=chosen['state'],
                residual=chosen['residual'],history=C['initial']['history'],fraction=chosen['fraction'],
                audit=h.pin(prior_path),seconds=time.monotonic()-started,new_native_trials=0,script=h.pin(__file__)))
            return
    terminal=e/'terminal-audit.json'
    if terminal.exists():
        audit=h.load(terminal)
        retained=C['selected']
        if audit.get('success') and all(chosen[k]['sha256']==retained[k]['sha256'] for k in ('state','residual')):
            assert C['native_history_verified']
            h.save(a.output,dict(policy=a.policy,state=chosen['state'],residual=chosen['residual'],
                history=C['initial']['history'],fraction=chosen['fraction'],audit=h.pin(terminal),
                status='reused_retained_complete_audit',prior_phases=audit['phases'],
                seconds=time.monotonic()-started,new_native_trials=0,script=h.pin(__file__)))
            return
    r10=h.load(evidence/'m8-r10/decision-ledger.json')
    old=h.load(h.bound(r10['final']['snapshot']))
    if chosen['state']['sha256']==old['state']['sha256'] and chosen['residual']['sha256']==old['residual']['sha256'] and C['initial']['history']['sha256']==old['history']['sha256']:
        prior=h.load(h.bound(r10['final']['audit']))
        assert prior['success'] and C['native_history_verified']
        h.save(a.output,dict(policy=a.policy,state=chosen['state'],residual=chosen['residual'],
            history=C['initial']['history'],fraction=chosen['fraction'],audit=r10['final']['audit'],
            status='reused_exact_complete_audit',prior_phases=prior['phases'],seconds=time.monotonic()-started,
            new_native_trials=0,script=h.pin(__file__)))
        return
    sys.path.insert(0,str(e/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
    import numpy as np
    from decimal import Decimal
    import cls_epawf_m8_r2_stage_analysis as analysis
    from cls_epawf_m8_r3_resolution import ResolutionAudit
    receipt=h.load(h.bound(C['receipt']))
    manifest=h.load(h.bound(receipt['manifest']))
    raw=h.load(h.bound(manifest['current_preflight']))

    class Capture(analysis.Analyzer):
        range_inventory=None
        range_domain=None

        def domain(self,state,native=None,retain=False):
            if native is self.range_inventory and native is not None:
                return self.range_domain
            return super().domain(state,native,retain)

        def domain_result(self,admitted,phase_name,records,counts,margins,minimum,native):
            self.retained_rows=records
            return super().domain_result(admitted,phase_name,records,counts,margins,minimum,native)

    analyzer=Capture(dict(raw,stage_test_target_time=.0025),h.bound(receipt['session']),e/'source')
    metric=ResolutionAudit(raw,e/'source')
    state=np.array(h.vector(chosen['state']))
    history=np.array(h.vector(C['initial']['history']))
    native=np.array(h.vector(chosen['residual']))
    domain=analyzer.domain(state,chosen['scalar'],retain=True)
    assert domain['admissible']
    analyzer.range_inventory={'observations':analyzer.retained_rows}
    analyzer.range_domain=domain
    constraints=analyzer.constraints(state)
    effective=analyzer.external_lift(state)
    loads,quadrature=analysis.weak_loads(raw,analyzer.fields,effective)
    unmasked=loads+analyzer.time_mass(effective,history,400.)
    independent=unmasked.copy()
    for row,target in constraints.items():independent[row]=state[row]-target
    error=float(np.max(np.abs(native-independent)/np.maximum(1.,np.abs(independent))))
    assert error<=2e-10
    physical=analyzer.physical(state,history,400.,unmasked=unmasked,constraints=constraints,
                              native_domain=analyzer.range_inventory)
    common=physical['common_inventory_physical_qualification']
    assert common['qualified_on_common_inventory']
    _,E,D,redistribution,*_=physical['source_energy']['robust_source']
    source_error=E-D-redistribution
    energy_error=physical['energy']['direct_remainder']
    assert max(abs(source_error),abs(redistribution),abs(energy_error))<3e-10
    wall=[]
    for precision in (80,110):
        rows=h.wall_rows(raw,state,precision)
        row=min(rows,key=lambda r:Decimal(r['minimum']['value']))
        wall.append(dict(row['minimum'],geometry=row['geometry']))
    assert Decimal(wall[-1]['value'])>0 and abs(Decimal(wall[0]['value'])-Decimal(wall[1]['value']))<Decimal('1e-65')
    h.save(a.output,dict(policy=a.policy,status='independently_audited',state=chosen['state'],
        residual=chosen['residual'],history=C['initial']['history'],fraction=chosen['fraction'],
        native_blocks=metric.norms(native,constraints).tolist(),independent_blocks=metric.norms(independent,constraints).tolist(),
        difference_blocks=metric.norms(native-independent,constraints).tolist(),component_value_error=error,
        domain={k:v for k,v in domain.items() if k not in ('independent_trace','independent_branch_context','profile_branch_pattern')},
        trace=domain['independent_trace'],ground=wall[-1],
        common=common,accounting=dict(source_identity_error=source_error,redistribution=redistribution,energy_error=energy_error),
        energy=physical['energy'],traction=physical['full_total_traction_and_k_flux'],
        seconds=time.monotonic()-started,new_native_trials=0,producing_revision=receipt['source_commit'],
        script=h.pin(__file__),comparison=h.pin(comparison)))


if __name__=='__main__':
    main()
