#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Two exact saved H20 owners: sparse replay and signed endpoint algebra only."""
import argparse
import csv
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import time
from run_wj_h4 import THREADS, load
os.environ.update(THREADS)  # before numerical imports
import numpy as np
from wj_h4_points import owner_point, field, vector_hash
from wj_h12_phase import local, decompositions, identity, cancellation
from wj_h15_observer import atomic, identity as pin, verify, bits
from wj_h17_observer import equal, native_load
from wj_h20_observer import resolve, read_anchor
from wj_h20_analysis import DT, TIMES, absolute_step, selected, ratio
from wj_h19_analysis import fixed_action

OLD = Path('/home/jay/build/wj-h20-three-level-resume-20261115')
POINT = [15., 0., .25607517888672315]
OWNERS = ((251, 'upper'), (265, 'lower'))
PAIRS = (('C', 'M'), ('M', 'F'), ('C', 'F'))
ROLES = tuple((t, z, absolute_step(z, t), geo, side)
              for t in TIMES for z in DT for geo, side in OWNERS)
NAMES = ('U', 'G', 'devG', 'k', 'epsilon', 'pi', 'R', 'tau', 'B', 'd', 'nu_eff',
         'a', 'g', 'mhat', 'Praw', 'N', 'F_N', 'transverse')
K = .008809015145515192
ANCHOR_HASHES = dict(state='631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077',
                     fixed_reference='9e2af3379f05a2ffbadc2cee491aefcab33900674c8b7ca6af380c922399d849',
                     fixed_defect='cf3dba7cf72ffc915d727bae7aad2b12c60eff2f1bb077159e9ee95a71268cbd')


def scale(name):
    if name in ('k', 'R'):
        return K
    if name == 'epsilon':
        return .051**3 / (.381917863546733 * .15)
    if name == 'pi':
        return .051**2
    if name in ('d', 'nu_eff'):
        return .001
    return 1.


def plain(v):
    if isinstance(v, np.ndarray):
        return plain(v.tolist())
    if isinstance(v, np.generic):
        return plain(v.item())
    if isinstance(v, float) and not np.isfinite(v):
        return None  # only descriptive unresolved ratios; no field replacement
    if isinstance(v, dict):
        return {k: plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [plain(x) for x in v]
    return v


def write(path, value):
    atomic(path, plain(value))


def checked_owner(structure, record, index):
    geo, side = OWNERS[index]
    equal(record['physical'], POINT, 'exact physical coordinate bits')
    o = record['owner']
    if (o['geometry'], o['normal_side'], o['x_side'], o['region']) != (geo, side, 'upper', 'outer'):
        raise ValueError('wrong exact outlet owner')
    e = next(v for v in structure['fields'][0]['mesh'] if v['geometry'] == geo)
    equal(owner_point(structure, e, POINT), o, 'owner bounds, sides and roundoff snap')
    return o


def linear(structure, state, record, index):
    o = checked_owner(structure, record, index)
    U = np.zeros(3)
    G = np.zeros((3, 3))
    for block, axis in ((0, 0), (1, 2)):
        U[axis], G[axis] = field(structure, state, block, o)
    k, gk = field(structure, state, 2, o)
    e, ge = field(structure, state, 3, o)
    pi, gp = field(structure, state, 5, o)
    return dict(U=U, G=G, devG=G-np.trace(G)*np.eye(3)/3, trace_G=np.trace(G),
                k=k, epsilon=e, pi=pi, mechanical_pressure=pi-2*k/3,
                grad_k=gk, grad_epsilon=ge, grad_pi=gp)


def endpoint_window(c0, c1, f0, f1, name):
    """True endpoint increments; never labelled flow/k/epsilon native phases."""
    if any(name not in v for v in (c0, c1, f0, f1)):
        return dict(status='NOT_EVALUABLE')
    incoming = np.asarray(c0[name]) - f0[name]
    coarse = np.asarray(c1[name]) - c0[name]
    finer = np.asarray(f1[name]) - f0[name]
    end = np.asarray(c1[name]) - f1[name]
    net = coarse - finer
    return dict(status='PASS',incoming=incoming,coarser_endpoint_increment=coarse,
                minus_finer_endpoint_increment=-finer,net_update=net,endpoint=end,
                identity=identity(end-incoming-coarse+finer,[end,incoming,coarse,finer],scale(name)),
                branch_cancellation=cancellation(incoming,[coarse,-finer],end,scale(name)),
                net_cancellation=cancellation(incoming,[net],end,scale(name)))


def comparison(actual, expected, s):
    err = float(np.max(abs(np.asarray(actual)-expected)))
    assert np.isfinite(err) and err <= 2e-11*s, 'published numerical anchor differs'
    return dict(absolute_error=err,scaled_error=err/s,tolerance=2e-11*s)


def complete_structure(bounded, original):
    """Bounded field metadata and original model metadata have separate pins."""
    for name in ('fields','offsets','geometry_signature'):
        equal(bounded[name],original[name],'unchanged original '+name)
    result=dict(bounded)
    for name in ('parameters','map_parameters','quadrature_orders'):
        if name in bounded:
            equal(bounded[name],original[name],'original '+name)
        result[name]=original[name]
    return result


def geometry_compatible(a,b):
    # Structure bytes contain dt-bearing signature/intent as well as FE rows.
    # Exclude only these execution identity fields from spatial equality.
    for name in (set(a)|set(b))-{'intent','physical_signature_sha256'}:
        equal(a[name],b[name],'same spatial structure '+name)


def prepare(root):
    index = load(root/'controls/h20-index.json')
    pub = load(verify(index['publication']))
    verify(index['closure']); closed = load(verify(index['accounting']))
    assert closed['closed'], 'H20 allocation must remain closed'
    pins = {Path(v['path']).relative_to(OLD/'publication-v2').as_posix():v for v in pub['files']}
    exported = load(verify(pins['fixed-probe-triplets.json']))
    inputs = load(verify(pins['controls/inputs--selected.json']))
    heavy = load(verify(pins['heavy-manifest.json']))
    heavies = {v['path']:v for v in heavy['files']}
    selected_anchors = []
    owners = None
    for triplet in exported:
        assert triplet['time'] in TIMES
        records = [r for r in triplet['records'] if r['name']=='outer_principal']
        assert len(records)==2
        records = [next(r for r in records if r['owner']['geometry']==geo) for geo,_ in OWNERS]
        if owners is None:
            owners = [{k:r[k] for k in ('physical','owner')} for r in records]
        else:
            equal([{k:r[k] for k in ('physical','owner')} for r in records],owners,'fixed owners across times')
        for z in DT:
            n = absolute_step(z,triplet['time'])
            a = resolve(OLD,inputs,z,triplet['time'])
            # Newly computed anchors are sealed both by H20's publication and
            # by its selected heavy manifest, not just a mutable live path.
            if str(n) not in inputs['states'].get(z,{}):
                sealed = pins[f'{z}/step-{n}.qualification-commit.json']
                verify(sealed)
                assert a['qualified_commit']['sha256']==sealed['sha256']
                for key in ('report','checkpoint','callback'):
                    assert a[key]==heavies[a[key]['path']]
            for value in a.values():
                if isinstance(value,dict) and set(('path','bytes','sha256')) <= value.keys():
                    verify(value)
            selected_anchors.append(dict(level=z,time=triplet['time'],step=n,anchor=a))
    assert len(selected_anchors)==12 and len(ROLES)==24
    q=load(selected_anchors[-1]['anchor']['qualified_commit']['path'])
    structure_pin=pin(OLD/'science/F/monitor/structure.json')
    assert structure_pin['sha256']==q['binding']['structure_sha256']
    bounded=load(verify(structure_pin))
    original=native_load(verify(inputs['original']['producer']))
    assert vector_hash(original['state'])==ANCHOR_HASHES['state']
    structure=complete_structure(bounded,original)
    for j,o in enumerate(owners): checked_owner(structure,o,j)
    structure_bindings=[]
    for item in selected_anchors:
        q=load(item['anchor']['qualified_commit']['path'])
        expected=q['binding'].get('structure_sha256')
        if expected is None:
            # M2's H18 B startup verification is a separate audit of the
            # original A2 completed producer, and has the older compact schema.
            assert (item['level'],item['time'],q['schema'])==('M',.005,'WJStartupVerificationCommitV1')
            structure_bindings.append(dict(level='M',time=.005,mode='H20 shared spatial support plus original H14 field-map identity',
                 sources=[structure_pin,inputs['original']['producer']],original_producer_and_own_callback_retained=True))
        else:
            own=pin(Path(item['anchor']['report']['path']).parent/'structure.json')
            assert own['sha256']==expected
            geometry_compatible(bounded,load(verify(own)))
            structure_bindings.append(dict(level=item['level'],time=item['time'],mode='own sealed structure, exact common geometry',source=own))
    expected=load(root/'controls/helper-pins.json')
    for name,digest in expected.items():
        module=importlib.import_module(name)
        assert pin(module.__file__)['sha256']==digest,'changed numerical helper '+name
    imports=[pin(m.__file__) for _,m in sorted(sys.modules.items())
             if getattr(m,'__file__',None) and '/Tests/RANS/' in str(m.__file__)]
    write(root/'controls/mapped-structure.json',structure)
    write(root/'controls/selected.json',dict(schema='WJH21InputsV1',H20_index=index,
          H20_ledger_closed=True,structure=pin(root/'controls/mapped-structure.json'),
          structure_inputs=[structure_pin,inputs['original']['producer']],structure_bindings=structure_bindings,owners=owners,states=selected_anchors,
          exported_records=pins['fixed-probe-triplets.json'],imported_helpers=imports,
          source_commit=load(root/'controls/source.json'),active_analysis=pin(__file__),roles=ROLES,
          permitted_aliases='H20 publication copies have identical hashes to live selected producer/commit bytes'))
    write(root/'controls/python-runtime.json',dict(python=sys.version,executable=pin(sys.executable),
          numpy=np.__version__,numpy_module=pin(np.__file__),threads={k:os.environ[k] for k in THREADS},
          loaded_libraries=[pin(p) for p in sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
              if len(line.split())>=6 and line.split()[-1].startswith('/') and Path(line.split()[-1]).is_file()})]))
    print('prepared 12 immutable accepted anchors, 24 owner roles',flush=True)


def science(root):
    started=time.monotonic();selection=load(root/'controls/selected.json')
    for v in selection['imported_helpers']:verify(v)
    structure=load(verify(selection['structure']))
    published=load(verify(selection['exported_records']))
    out=root/'science-v1';out.mkdir(exist_ok=False)
    records=[];lineage=[];replays=[];shared=[];costs=[]
    for i,item in enumerate(selection['states']):
        began=time.monotonic();z=item['level'];t=item['time'];anchor=item['anchor']
        r,audit=read_anchor(anchor,z,t)
        c=selected(r,z,t)
        assert r['original_qh0_state_sha256']==ANCHOR_HASHES['state']
        for name in ('fixed_reference','fixed_defect'):
            assert vector_hash(r[name])==ANCHOR_HASHES[name]
        lineage.append(dict(level=z,time=t,absolute_step=c['step'],dt=r['dt'],time_bits=bits(c['time']),alpha=c['alpha'],
             state_sha256=vector_hash(r['state']),callback_history_sha256=vector_hash(c['history']),
             native_history_sha256=[vector_hash(v)for v in r['native_history']],native_history_times=r['native_history_times'],
             fixed_reference_sha256=vector_hash(r['fixed_reference']),fixed_defect_sha256=vector_hash(r['fixed_defect']),
             physical_signature_sha256=hashlib.sha256(r['signature'].encode()).hexdigest(),anchor=anchor,
             role='accepted endpoint with genuine own producer; never an intermediate phase'))
        pub=next(v for v in published if v['time']==t)
        for j,owner in enumerate(selection['owners']):
            v=local(linear(structure,r['state'],owner,j))
            assert v['constitutive_status']=='DEFINED','saved endpoint closure undefined'
            old=next(vv for vv in pub['records'] if vv['name']=='outer_principal' and vv['owner']['geometry']==OWNERS[j][0])
            checks={name:comparison(v[name],old['states'][z][name],scale(name)) for name in ('U','G','k','epsilon','R','pi','nu_eff')}
            for name,alias in (('mhat','margin'),('Praw','production'),('transverse','transverse'),('covariance','covariance'),('N','root')):
                checks[name]=comparison(v[name],old['states'][z][alias],scale(name))
            replays.append(dict(level=z,time=t,geometry=OWNERS[j][0],checks=checks))
            records.append(dict(level=z,time=t,step=c['step'],geometry=OWNERS[j][0],owner=owner,values=v))
            for source in DT:
                angle=old['states'][source]['angle']
                values=fixed_action(np.array([v['k']]),np.array([v['epsilon']]),v['G'][None],np.array([angle]))
                expected=old['principal']['theta_'+source]['values'][z]
                comparison(values['normalized'][0],expected['normalized'],1.)
                comparison(values['dimensional'][0],expected['dimensional'],.001)
                shared.append(dict(time=t,geometry=OWNERS[j][0],angle_from=source,level=z,
                                  values={k:x[0]for k,x in values.items()},published_error=float(values['normalized'][0]-expected['normalized'])))
        elapsed=time.monotonic()-began;costs.append(elapsed)
        if i==0:
            # This retained first state is the pilot; no duplicate evaluation.
            remaining=11*elapsed;forecast=1.25*(remaining+30)
            assert time.monotonic()-started+forecast<180,'RESOURCE_LIMIT science remainder'
            write(out/'pilot-admission.json',dict(measured_first_state_seconds=elapsed,remaining_state_estimate_seconds=remaining,
                  analysis_figure_allowance_seconds=30,headroom=1.25,forecast_seconds=forecast,closure_reserved_seconds=90,status='PASS'))
        print(z,c['step'],t,'sparse replay PASS',elapsed,flush=True)
    lookup={(r['time'],r['level'],r['geometry']):r['values']for r in records}
    pairs=[];windows=[];increments=[];summary=[]
    for geo,_ in OWNERS:
        for t in TIMES:
            a,b,c=(lookup[t,z,geo]for z in DT)
            summary.append(dict(time=t,geometry=geo,minimized=ratio(abs(a['mhat']-b['mhat']),abs(b['mhat']-c['mhat'])),
                fixed=ratio(abs(a['a']-b['a']),abs(b['a']-c['a'])),signed_minimized_CM=a['mhat']-b['mhat'],
                signed_minimized_MF=b['mhat']-c['mhat'],signed_fixed_CM=a['a']-b['a'],signed_fixed_MF=b['a']-c['a']))
            for coarse,finer in PAIRS:
                v=decompositions(lookup[t,finer,geo],lookup[t,coarse,geo])
                pairs.append(dict(time=t,geometry=geo,pair=coarse+'-'+finer,**v))
        for lo,hi in zip(TIMES[:-1],TIMES[1:]):
            for z in DT:
                increments.append(dict(geometry=geo,level=z,start=lo,end=hi,
                    decomposition=decompositions(lookup[lo,z,geo],lookup[hi,z,geo])))
            for coarse,finer in PAIRS:
                windows.append(dict(geometry=geo,pair=coarse+'-'+finer,start=lo,end=hi,
                    endpoint_quantities={name:endpoint_window(lookup[lo,coarse,geo],lookup[hi,coarse,geo],
                         lookup[lo,finer,geo],lookup[hi,finer,geo],name)for name in NAMES}))
    write(out/'lineage.json',lineage);write(out/'endpoint-values.json',records)
    write(out/'H20-replay.json',replays);write(out/'shared-directions.json',shared)
    write(out/'pair-decompositions.json',pairs);write(out/'endpoint-increments.json',increments)
    write(out/'endpoint-windows.json',windows);write(out/'principal-summary.json',summary)
    write(out/'counts.json',dict(accepted_states=12,endpoint_owner_roles=24,shared_direction_actions=72,
          pair_decompositions=24,branch_endpoint_interval_decompositions=18,paired_endpoint_windows=18,
          full_volume_rows=0,native_processes=0,kernel_requests=0,governing_attempts=0,cp0=0,strict_readbacks=0,
          builds=0,saved_governing_matrix_operations=0,mass_solves=0,phase_or_predictor_packets='UNAVAILABLE_NOT_INFERRED',
          implicit_root_tangents='Python complex-step and independent angular calculation at actual saved states',
          state_replay_seconds=costs,total_science_seconds=time.monotonic()-started))
    print('all 24 endpoint-owner roles and signed identities PASS',flush=True)


def publish(root):
    out=root/'publication-v1';out.mkdir(exist_ok=False)
    science=root/'science-v1'
    import shutil
    for p in sorted(science.glob('*.json')):shutil.copyfile(p,out/p.name)
    for p in sorted((root/'controls').glob('*.json')):shutil.copyfile(p,out/('control--'+p.name))
    summary=load(science/'principal-summary.json');pairs=load(science/'pair-decompositions.json')
    windows=load(science/'endpoint-windows.json');values=load(science/'endpoint-values.json')
    def table(name,rows):
        with (out/name).open('x')as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    table('principal-differences.csv',[dict(time=r['time'],geometry=r['geometry'],E1_min=r['signed_minimized_CM'],E2_min=r['signed_minimized_MF'],
        rho_min=r['minimized']['rho'],E1_fixed=r['signed_fixed_CM'],E2_fixed=r['signed_fixed_MF'],rho_fixed=r['fixed']['rho'])for r in summary])
    table('action-decomposition.csv',[{k:r[k]for k in ('time','geometry','pair','delta_d','delta_inverse_nu','dimensional_term',
         'normalization_term','delta_a','delta_g','delta_mhat')}for r in pairs])
    rows=[]
    for r in windows:
        for name,v in r['endpoint_quantities'].items():
            shape=np.shape(v['endpoint'])
            for index in np.ndindex(shape) if shape else [()]:
                get=lambda x:float(np.asarray(x)[index])
                rows.append(dict(geometry=r['geometry'],pair=r['pair'],start=r['start'],end=r['end'],quantity=name,component=str(index),
                    incoming=get(v['incoming']),coarser_increment=get(v['coarser_endpoint_increment']),minus_finer_increment=get(v['minus_finer_endpoint_increment']),
                    net_update=get(v['net_update']),endpoint=get(v['endpoint']),identity_residual=get(v['identity']['residual']),
                    net_cancellation=None if v['net_cancellation']['cancellation'] is None else np.asarray(v['net_cancellation']['cancellation'],dtype=object)[index]))
    table('signed-endpoint-windows.csv',rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9})
    fig,axes=plt.subplots(1,2,figsize=(10,3.7),sharey=True)
    for ax,(geo,side)in zip(axes,OWNERS):
        rr=[r for r in summary if r['geometry']==geo]
        for key,label,style in (('signed_minimized_CM','C-M minimized','o-'),('signed_minimized_MF','M-F minimized','s-'),
                                ('signed_fixed_CM','C-M fixed H*','o--'),('signed_fixed_MF','M-F fixed H*','s--')):
            ax.plot([r['time']for r in rr],[r[key]*1e6 for r in rr],style,label=label,linewidth=1.2,markersize=4)
        ax.axhline(0,color='0.4',linewidth=.6);ax.set_title(f'Geometry {geo}, {side} normal owner')
        ax.set_xlabel('Physical time (H_h / U_infinity)');ax.set_xticks(TIMES);ax.grid(alpha=.2)
    axes[0].set_ylabel('Signed normalized-action difference (10^-6)');axes[1].legend(fontsize=8)
    fig.suptitle('Outlet x=15, physical height 10.2430071555 mm; separate one-sided traces')
    fig.tight_layout();fig.savefig(out/'signed-owner-history.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for ax,(geo,side)in zip(axes,OWNERS):
        r=next(v for v in windows if v['geometry']==geo and v['pair']=='M-F' and v['start']==.015)
        for j,name in enumerate(('a','g','mhat')):
            v=r['endpoint_quantities'][name]
            ax.bar(np.arange(3)+j*.23-.23,[v['incoming']*1e6,v['net_update']*1e6,v['endpoint']*1e6],width=.23,label=name)
        ax.axhline(0,color='0.3',linewidth=.6);ax.set_xticks(range(3),['Incoming .015','Net update .015–.020','Endpoint .020'],rotation=12)
        ax.set_title(f'Geometry {geo}, {side} normal owner');ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Signed M-F difference (10^-6)');axes[1].legend(title='a: fixed; g: gap; mhat: min',fontsize=8)
    fig.suptitle('Terminal cancellation: inherited and newly accumulated endpoint differences')
    fig.tight_layout();fig.savefig(out/'terminal-cancellation.png',dpi=160);plt.close(fig)
    write(out/'component-gates.json',dict(input_and_endpoint_replay='PASS',tensor_and_normalization_identities='PASS',
        endpoint_telescopes='PASS',new_native_qualification='NOT_RUN',H20_uniform_temporal_contraction='NOT_DEMONSTRATED_UNCHANGED',
        underlying_PDE_stage_causality='UNRESOLVED_NO_SAVED_PHASE_PACKETS',scope='24 accepted endpoint-owner records only'))
    write(out/'accounting-publication-snapshot.json',load(root/'accounting.json'))
    write(out/'manifest.json',dict(files=[pin(p)for p in sorted(out.iterdir())if p.is_file()],
         scope='compact H21 publication; OPEN accounting snapshot is not final closure'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--mode',choices=('prepare','science','publish'),required=True);a=p.parse_args()
    assert a.root.name=='wj-h21-outlet-saved-state-20261115'
    dict(prepare=prepare,science=science,publish=publish)[a.mode](a.root)
