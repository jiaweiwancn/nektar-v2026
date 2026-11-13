#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify H10's immutable native matrix and close its existing allocation."""
import argparse
from pathlib import Path
import re
import shlex
import subprocess
from run_wj_h10 import load,save,size
from study_wj_h9 import identity,check
from study_wj_h10 import selected_analysis
from wj_h10_analysis import MATRIX,QUALITY
from wj_h9_local import same
from wj_h4_points import vector_hash


def verify(root):
    out=root/'closure';out.mkdir();pub=root/'publication-v1';verified=[]
    for v in load(pub/'seal.json')['files']:verified.append(check(v))
    for v in load(pub/'manifest.json')['heavy']:check(v)
    for v in load(root/'science/input-intent.json')['records']:
        verified.extend(check(v[k])for k in ('original','copy'))
    gates=load(pub/'component-gates.json');assert gates['status']=='PASS'
    assert (gates['new_CP0'],gates['new_governing'],gates['new_strict'])==(3,3,6)
    runtime=load(root/'science/runtime-build.json')
    for v in runtime['executables']:
        for k in ('ELF','wrapper'):verified.append(check(v[k]))
    for v in runtime['private_libraries']:verified.append(check(v))
    verified.extend(check(runtime[k])for k in ('kernel','geometry','configured_commands'))
    bindings=[];pre=None;unique=set();methods=set();test_logs={}
    for branch,dt in MATRIX:
        original=load(root/f'science/imports/{branch}.report.json')
        cp0=load(root/f'science/cases/{branch}/cp0/case.evolution.json')
        assert same(cp0['initial_velocity_diagnostic']['pre_state'],original['accepted_states'][0])
        if pre is not None:
            assert all(same(cp0[k],pre[k])for k in ('state','fixed_reference','fixed_defect'))
        pre=cp0
        for phase in ('cp0','run'):
            v=selected_analysis(root,branch,phase);assert v['status']=='PASS'
            verified.extend(check(v[k])for k in ('producer','checkpoint'));r=load(v['producer']['path'])
            assert len(r['state'])==2328 and not r['failure']
            assert r['dt']==dt and r['steps']==(phase=='run')
            assert r['initial_velocity_calls']==(phase=='cp0')
            assert phase=='cp0'or r['initialization_lift_calls_iterations']==[0,0]
            unique.add(vector_hash(r['state']))
            view=load(root/f'science/views/{branch}-{phase}/result.json');assert view['status']=='PASS'
            assert view['anchor']['status']==view['quadrature']['status']=='PASS'
            n=load(view['native']['path']);ctx=view['context']
            assert n['initial_velocity_mode']=='FlatInletTraceV1'and n['initial_velocity_calls']==0
            assert n['initialization_work']==[0,0]and n['governing_solves']==0
            assert n['state_sha256']==ctx['state_sha256']==vector_hash(r['state'])
            assert n['report_sha256']==ctx['report']['sha256']==v['producer']['sha256']
            assert n['checkpoint_sha256']==ctx['checkpoint']['sha256']==v['checkpoint']['sha256']
            assert n['callback_history_sha256']==ctx['history_sha256']
            assert n['weak_checks']['alpha_mass']==ctx['alpha']and n['open_boundary']['time']==ctx['time']
            for side,values in view['native_Q1_moments'].items():
                assert len(values)=={'4':34,'5':70}[side]and max(map(abs,values))<=2e-15
            sample_label=branch+('-step'if phase=='run'else '-cp0')
            for label in (sample_label,f'{branch}-{phase}-strict'):
                samples=load(root/'operations'/label/'native-exec-samples.json');assert samples
                assert all(s['status']=='PASS'and s['libc']==[str(root/'runtime/glibc/libc.so.6')]for s in samples)
            if phase=='run':
                ar=load(v['arithmetic']['path']);assert ar['status']=='PASS'
                assert all(z['status']=='PASS'for z in ar['checks'].values())
                assert ar['binding']['P_C_bitwise']and ar['binding']['raw_corrected_bitwise']
            bindings.append(dict(branch=branch,phase=phase,producer=v['producer'],checkpoint=v['checkpoint'],strict=check(view['native']),
                state=n['state_sha256'],history=n['callback_history_sha256'],alpha=ctx['alpha'],time=ctx['time'],initializer_calls=r['initial_velocity_calls']))
    common=load(root/'science/final-common/result.json');assert common['status']=='PASS'and len(unique)==4
    assert all(v['min_in_plane']>=QUALITY for v in common['states'].values())
    for label in ('focused-tests-v2','focused-h10-audit-recovery'):
        folder=root/'operations'/label;text=(folder/'command.log').read_text();assert '\nOK\n'in text
        found=re.findall(r'^(test_\w+) \(([^)]+)\)',text,re.M)
        assert len(found)==int(re.search(r'Ran (\d+) tests',text).group(1));methods.update(found)
        test_logs[label]={k:identity(folder/k)for k in ('command.log','request.json','receipt.json')}
    assert len(methods)==63
    # Preserve the compiler's actual dependency files, not just configured .cpp names.
    dependencies=set();depfiles=[]
    for p in sorted((root/'build').rglob('*.o.d')):
        depfiles.append(identity(p));text=p.read_text().replace('\\\n',' ')
        dependencies.update(Path(v)for v in shlex.split(text.partition(':')[2])if Path(v).is_file())
    compiler=Path('/usr/bin/c++').resolve()
    save(out/'compiled-inputs.json',dict(compiler=identity(compiler),compiler_version=subprocess.check_output([str(compiler),'--version'],text=True),
        compile_commands=identity(root/'build/compile_commands.json'),cache=identity(root/'build/CMakeCache.txt'),
        native_source_repair=check(identity(root/'native-source-repair.json')),dependencies=[identity(p)for p in sorted(dependencies)],dependency_files=depfiles,
        scope='governor ca82fdb; reader ca82fdb plus committed9447f29 synthetic fixture header; runtime and Python revisions remain separate'))
    save(out/'verification.json',dict(status='PASS',verified_compact_and_imports=verified,native_bindings=bindings,unique_physical_states=4,
        methods=[list(v)for v in sorted(methods)],unique_passing_methods=63,test_logs=test_logs,
        native_fixture=check(identity(root/'science/fixture-v2/native.json')),compiled_inputs=identity(out/'compiled-inputs.json'),
        new_CP0s=3,new_governing_attempts=3,new_strict_readbacks=6,new_accepted_steps=3,bound_calls=3,bound_activations=0,
        passive_pressure_mass_algebra=dict(completed_scans_and_final_common_factorizations=364,completed_solves=32,
            failed_coarse_scan_factorizations=52,failed_coarse_scan_solves=2,
            scope='Rows/MomentAudit algebra only; native mass-dual diagnostics are separately labelled, never governing solves'),
        kernel_queries=14,scope='all heavy science hashes verified on Linux; Mac independently verifies compact transport only',source=identity(__file__)))
    print('H10 verification PASS: 63 methods, 3 CP0s, 3 steps, 6 strict views',flush=True)


def finalize(root):
    out=root/'closure';assert load(out/'verification.json')['status']=='PASS';book=load(root/'accounting.json');assert not book.get('closed')
    assert not book.get('h10_halt')and not book.get('coverage_halt')
    assert all('seconds'in v and v['status']!='RESERVED'for v in book['entries'].values())
    actual=sum(v['seconds']for v in book['entries'].values());build=60+sum(v['seconds']for v in book['entries'].values()if v['phase']=='build')
    science=sum(v['seconds']for v in book['entries'].values()if v['phase']=='science');total=60+actual+300
    assert build<=2400 and science<=900 and total<=3600
    peak=max(v['peak_tree_rss_bytes']for v in book['entries'].values());build_bytes=sum(size(root/p)for p in ('source','build','runtime'));science_bytes=size(root/'science')
    assert peak<=4*1024**3 and build_bytes<=8*1024**3 and science_bytes<=2*1024**3
    assert(book['cp0_reservations'],book['governing_attempts'],book['strict_readbacks'])==(3,3,6)
    operations=[]
    for label in book['entries']:
        p=root/'operations'/label;operations.append(dict(label=label,**{k:identity(p/(k+'.json'))for k in ('request','receipt')},log=identity(p/'command.log'),maps=identity(p/'process-maps.json')))
    book.update(closed=True,status='PASS',measured_operation_seconds=actual,charged_seconds=total,build_phase_including_setup=build,science_seconds=science,
        fully_charged_closure_seconds=300,peak_process_tree_RSS=peak,build_source_runtime_bytes=build_bytes,science_bytes=science_bytes,
        closure_scope='full300s includes parent inspections, source/figure QA, verification orchestration, transport and final receipt work; not reclaimed')
    save(root/'accounting.json',book);save(out/'accounting-final.json',book)
    pub=root/'publication-v1';save(out/'bridge.json',dict(status='PASS',publication_seal=identity(pub/'seal.json'),verification=identity(out/'verification.json'),
        final_accounting=identity(out/'accounting-final.json'),operations=operations,setup=identity(root/'setup.json'),
        source_versions=[identity(p)for directory in sorted(root.glob('observer-v*'))if directory.is_dir()for p in sorted(directory.glob('*.py'))]+[identity(p)for p in sorted(root.glob('observer-v*.py'))],
        superseded=[identity(root/'operations'/k/'receipt.json')for k,v in book['entries'].items()if v['status']!='PASS'],
        original_H9_publication_seal=check(identity(root/'source/docs/plans/wj-h9-linux/publication-v3/seal.json')),
        resource_transport='remote heavy verification complete; compact transport hash is independently checked on Mac',
        historical_outcomes='H8 early common-time contraction UNRESOLVED; old ledgers and all previous outcomes unchanged',source=identity(__file__)))
    save(out/'seal.json',dict(status='PASS',files=[identity(p)for p in sorted(out.iterdir())if p.is_file()]))
    print(dict(status='PASS',charged_seconds=total,build_including_setup=build,science=science,peak_RSS=peak,build_bytes=build_bytes,science_bytes=science_bytes),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--finalize',action='store_true');a=p.parse_args();(finalize if a.finalize else verify)(a.root)
