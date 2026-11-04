#!/usr/bin/env python3
"""Prepare bound local R6 replays and audit the three cut representations."""
import argparse
from decimal import Decimal, localcontext
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import xml.etree.ElementTree as ET

THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
sys.dont_write_bytecode = True
for name in THREADS:
    os.environ[name] = '1'


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def prepare(epoch, cases, h, policy):
    source = epoch / 'source'
    old = h.load(source / 'docs/plans/cls-epawf-revised-evidence/m8-r5/decision-ledger.json')
    pilot, local = old['pilot'], old['local']
    rr = h.load(h.bound(pilot['receipt']))
    lr = h.load(h.bound(local['receipt']))
    old_manifest = h.load(h.bound(rr['manifest']))
    runtime_file = epoch / 'runtime-fingerprint.json'
    runtime = h.load(runtime_file)
    exe = h.pin(epoch / 'build/IncNavierStokesSolver')
    h.require(runtime['success'] and runtime['runtime_artifacts']['build/IncNavierStokesSolver'] == exe['sha256'],
              'Wrong R6 executable/runtime')
    c5 = h.load(h.bound(pilot['collection']))
    local5 = h.load(h.bound(lr['durable_collection']['file']))
    for case, collection, index in (('terminal', c5, 1), ('control', local5, 0), ('pilot', c5, 1)):
        if case == 'pilot' and policy != 'RefinedR6':
            continue
        base = collection['completed_iterates'][index]
        current, = (row for row in collection['jacobians'] if row['iteration'] == base['iteration'])
        matrix = current['event']['matrix']
        h.require(current['status'] == 'complete' and current['iteration'] == base['iteration'] == index and
            matrix == (pilot['matrices'][1]['matrix'] if index else old['replay']['current_terminal_matrix']),
            'Matrix does not belong to the retained state')
        h.bound(matrix)
        for key in ('snapshot','state','residual','constraints','history'):
            h.bound(base[key])
        h.require(base['history']['sha256'] == h.predictor.H0 and base['alpha'] == 400 and
            base['time'] == .0025 and not base['history_advanced'], 'Wrong retained history or stage')
        if index == 1:
            h.require(base['snapshot'] == pilot['last_snapshot']['snapshot'], 'Wrong terminal snapshot')
        out = cases / case
        out.mkdir(parents=True)
        tree = ET.parse(h.bound(rr['session']))
        info = tree.find('.//SOLVERINFO')
        selectors = {x.get('PROPERTY'):x for x in info}
        if 'RANSWallR5Updates' in selectors:
            info.remove(selectors['RANSWallR5Updates'])
        selectors['RANSWallStageTest'].set('VALUE','ExactTargetR3' if case == 'pilot' else 'DirectionExperimentR5')
        selectors['RANSWallR3PredictorFile'].set('VALUE',str(out/'predictor.json'))
        selectors['RANSWallR3ProgressFile'].set('VALUE',str(out/'case.r3-progress.jsonl'))
        ET.SubElement(info,'I',PROPERTY='RANSWallDirectionAccuracy',VALUE=policy)
        xml=out/'case.xml'
        tree.write(xml,encoding='unicode')
        h.require(h.predictor.protected_xml(xml) == old_manifest['protected_physical_xml_sha256'], 'Changed frozen target')
        binding = dict(matrix=matrix, state_sha256=base['state']['sha256'],
            history_sha256=base['history']['sha256'],residual_sha256=base['residual']['sha256'],
            assembly_json=json.dumps(current['event']['details']['assembly']))
        manifest = dict(old_manifest, state=base['state'], history=base['history'],
            producing_source_commit=runtime['source_commit'], current_executable=exe,
            current_target_xml=h.pin(xml), current_runtime=h.pin(runtime_file),
            origin=dict(alpha=400,time=.0025,converged=False,inner_stress_model='SaturatedR2',
                iteration=index,snapshot=base['snapshot'],residual=base['residual'],constraints=base['constraints']))
        if case != 'pilot':
            manifest['frozen_action']=binding
        h.save(out/'predictor.json',manifest)
        if case != 'pilot':
            probe_xml=Path(rr['session']['path']).with_name('case-m4-native-rank0.xml')
            probe=ET.parse(probe_xml)
            selector=probe.find(".//SOLVERINFO/I[@PROPERTY='RANSWallGlobalization']")
            selector.set('VALUE','JointRestorationR5')
            probe.write(out/'probe.xml',encoding='unicode')
            h.save(out/'packet.json',dict(state=h.vector(base['state']),history=h.vector(base['history']),
                alpha=400,time=.0025,only_external_direction=True))


def high_precision_active(normals, rhs, z, number_active, np):
    """Independent 60-digit Gram solve only for the small retained active set."""
    with localcontext() as ctx:
        ctx.prec=60
        n = [[Decimal(float(x)) for x in row] for row in normals]
        b = [Decimal(float(x)) for x in rhs]
        count = number_active
        # Unit row scaling keeps the Gram calculation in the QR coordinates.
        length=[sum(x*x for x in row).sqrt() for row in n]
        unit=[[x/l for x in row] for l,row in zip(length,n)]
        gram=[[sum(a*b for a,b in zip(row,col)) for col in unit[:count]]
              +[target/l] for row,target,l in zip(unit[:count],b,length)]
        for j in range(count):
            p=max(range(j,count),key=lambda i:abs(gram[i][j]))
            gram[j],gram[p]=gram[p],gram[j]
            g=gram[j][j]
            if not g:
                return dict(status='rank_deficient_at_60_digits')
            gram[j]=[x/g for x in gram[j]]
            for i in range(count):
                if i != j:
                    g=gram[i][j]
                    gram[i]=[v-g*w for v,w in zip(gram[i],gram[j])]
        dual=[row[-1] for row in gram]
        polished=[sum(mu*n_j for mu,n_j in zip(dual,col)) for col in zip(*unit[:count])]
        gaps=[sum(a*b for a,b in zip(row,polished))-target for row,target in zip(n,b)]
        return dict(status='solved',minimum_scaled_multiplier=float(min(dual)),
            max_abs_active_gap=float(max(map(abs,gaps[:count]))),
            entering_gap=float(gaps[-1]) if len(gaps)>count else None,
            maximum_z_difference=float(max(abs(x-Decimal(float(y))) for x,y in zip(polished,z))))


def audit(epoch, cases, h, kind):
    import numpy as np
    case=cases/kind
    receipt=h.load(cases/(kind+'-run')/'receipt.json')
    col=h.load(h.bound(receipt['durable_collection']['file']))
    man=h.load(h.bound(receipt['manifest']))
    base,=col['completed_iterates']
    assert all(base[key]['sha256'] == (man[key] if key in ('state','history') else
        man['origin'][key])['sha256'] for key in ('state','history','constraints','residual'))
    current,=col['jacobians']
    assert current['event']['matrix']['sha256']==man['frozen_action']['matrix']['sha256']
    assert current['event']['details']['loaded_bound_matrix']
    size=2165
    ld=np.longdouble
    vec=lambda pin: np.fromfile(h.bound(pin),dtype='<f8').astype(ld)
    matrix=np.memmap(h.bound(current['event']['matrix']),dtype='<f8',mode='r',shape=(size,size)).astype(ld)
    residual=vec(base['residual'])
    constraints=vec(base['constraints']).reshape(-1,2)[:,0].astype(int)
    events=[json.loads(line) for line in Path(col['progress']['path']).read_text().splitlines()]
    findings=[]
    for event in (row for row in events if row['kind']=='r6_cut_audit'):
        z,d,newton,Bz=(vec(event[key]) for key in ('z','direction','newton','Bz'))
        product=matrix@d
        error=product+residual-Bz
        correction=d-newton
        inverse_error=matrix@correction-Bz
        rows=[]
        n=[]
        for item in event['rows']:
            a,normal,x=(vec(item[key]) for key in ('a','n','transpose_solution'))
            n.append(normal)
            dual_error=a-matrix.T@x
            coefficient=np.dot(a,d)-ld(item['rhs'])
            transformed=np.dot(normal,z)-ld(item['transformed_rhs'])
            rows.append(dict(item,independent_coefficient_gap=float(coefficient),
                independent_transformed_gap=float(transformed),
                coefficient_minus_transformed=float(coefficient-transformed),
                transpose_residual_inf=float(np.max(np.abs(dual_error))),
                transpose_dot_error=float(np.dot(dual_error,correction)),
                forward_dot_error=float(np.dot(x,inverse_error)),
                metric_dot_error=float(np.dot(x,Bz)-np.dot(normal,z)),
                shift_roundoff=float(np.dot(a,newton)-ld(item['rhs'])+ld(item['transformed_rhs']))))
        lengths=np.array([row['normal_length'] for row in rows])
        active=[i for i,row in enumerate(rows) if not row['incoming']]
        unit=np.array(n,dtype=float)/lengths[:,None]
        singular=np.linalg.svd(unit[active],compute_uv=False)
        polish=high_precision_active(n,[row['transformed_rhs'] for row in rows],z,len(active),np)
        dual=np.array([rows[i]['multiplier_scaled'] for i in active])
        gaps=np.array([rows[i]['independent_transformed_gap']/lengths[i] for i in active])
        history=event['exchanges']
        findings.append(dict(name=event['name'],status=event['status'],exchange=event['exchange'],
            active_count=len(active),active_condition=float(singular[0]/singular[-1]),
            minimum_qr_diagonal=min(event['qr_diagonal']),minimum_scaled_multiplier=float(min(dual)),
            stationarity_l2=float(np.linalg.norm(np.asarray(z,float)-dual@unit[active])),
            complementarity=float(max(abs(gaps*dual))),
            forward_Jd_residual_inf=float(max(abs(error))),forward_Jcorrection_residual_inf=float(max(abs(inverse_error))),
            Jnewton_residual_inf=float(max(abs(matrix@newton+residual))),
            algebraic_residual_inf=float(max(abs(error[constraints]))),
            exact_recurrent_entries=sum(bool(x.get('exact_previous_alias_ids')) for x in history),
            minimum_projected_length=min((x['projected_length'] for x in history if 'projected_length' in x),default=None),
            original_tolerance=event['original_tolerance'],worst_native_violation=event['worst_native_violation'],
            worst_extended_violation=event['worst_extended_violation'],maximum_scan_difference=event['maximum_scan_difference'],
            high_precision_active=polish,rows=rows))
    h.require(findings,'No detailed R6 cut replay')
    h.save(case/'cut-audit.json',dict(matrix=current['event']['matrix'],base=base['snapshot'],
        native_progress=col['progress'],findings=findings,terminal=col['terminal']))
    for result in findings:
        print({key:result[key] for key in result if key!='rows'})
        for row in result['rows']:
            if row['incoming']:
                print({k:v for k,v in row.items() if not isinstance(v,dict)})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=('prepare','audit','qualify','budget'))
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--cases',type=Path,required=True)
    p.add_argument('--policy',default='ReplayR6',choices=('ReplayR6','RefinedR6'))
    p.add_argument('--kind',default='terminal',choices=('terminal','control'))
    a=p.parse_args()
    if platform.system()!='Linux':
        raise ValueError('Use the designated Linux host')
    epoch,cases=a.epoch.resolve(),a.cases.resolve()
    evidence=epoch/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r4_prepare')
    if a.mode=='prepare':
        prepare(epoch,cases,h,a.policy)
    elif a.mode=='audit':
        audit(epoch,cases,h,a.kind)
    else:
        r5=module(evidence/'m8-r5/compare.py','r5_compare')
        if a.mode=='qualify':
            r5.qualify(epoch,cases,h,49,case=a.kind)
        else:
            replay=h.load(cases/'terminal/cut-audit.json')
            h.require(any(x['status']=='passed_original_cut_scan' for x in replay['findings']),
                      'No verified original-row direction')
            old=h.load(evidence/'m8-r5/decision-ledger.json')
            previous=h.load(h.bound(old['pilot']['collection']))['terminal']['profiling']['complete_action_assembly']
            r5.stage_budget(cases,h,case='terminal',
                bound_action_seconds=1.2*previous['seconds']/previous['calls'])


if __name__=='__main__':
    main()
