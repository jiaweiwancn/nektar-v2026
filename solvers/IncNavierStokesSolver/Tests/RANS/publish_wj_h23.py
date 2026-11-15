#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Truthful compact H23 stopped-prefix publication from immutable saved bytes."""
import argparse,csv,gzip,json,math,os,shutil
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash,bits
from wj_h17_observer import native_load
from wj_h23_continuation import DT,START,END,STRIDE,TIMES,ORDER,absolute_step,qualified
from publish_wj_h20_resumed import table,saved_arrays
from run_wj_h23 import usage
from run_wj_h15 import size

def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');selected=load(root/'inputs/selected.json');copies=[];heavy=[]
    def copy(path,name):
        path=Path(path)
        if not path.exists():return
        q=out/name;q.parent.mkdir(parents=True,exist_ok=True)
        if path.suffix=='.json'and path.stat().st_size>100000:
            q=Path(str(q)+'.gz')
            with gzip.GzipFile(filename=str(q),mode='wb',mtime=0)as f:f.write(path.read_bytes())
            assert gzip.decompress(q.read_bytes())==path.read_bytes();method='lossless gzip'
        else:shutil.copyfile(path,q);method='byte copy'
        copies.append(dict(original=identity(path),stored=identity(q),transformation=method))
    outcome=load(root/'science/M/outcome.json')
    assert outcome['native_exit']is None and outcome['producer']is None and outcome['error']
    assert not(root/'science/M/native-launch.json').exists()and(root/'science/M/native.log').stat().st_size==0
    assert not list((root/'science').glob('*/monitor/step-*.completion.json'))
    assert book['native_reservations']==['M','strict_M']and book['native_halt']
    strict=load(root/'science/M/final-strict/result.json');assert strict['status']=='PASS'and strict['selected_absolute_step']==8 and strict['genuine_alpha']==400
    preflight=load(root/'science/preflight/result.json');assert preflight['status']==preflight['comparison']['status']=='PASS'
    counts=dict(governor_coordinator_launches=1,actual_governor_processes=0,actual_governing_attempts=0,native_completions=0,
        new_startup_verification_commits=0,new_qualified_endpoints=0,governor_scientific_audits=0,observer_scientific_events=0,
        observer_participating_state_assessments=0,passive_preflight_state_assessments=3,protected_strict_reads=1,
        strict_roles=['original qualified H20 M8/.02, genuine alpha400'],remaining_strict_reads='NOT_RUN: C/F governors not admitted',
        maximum_attempt_capacity_reserved=8,driver_unknown_capacity_charge_retained=8,actual_zero_established_by_outcome_and_absent_native_launch=True,
        new_CP0=0,bootstrap=0,spatial_initializations=0,history_productions=0,native_builds=0,normal_phase_vector_matrix_captures=0,
        diagnostic_mass_solves=0,saved_matrix_operations=0,new_scientific_time=0)
    atomic(out/'execution-counts.json',counts)
    planned=[]
    for z in ORDER:
        for n in range(START[z]+1,END[z]+1):
            planned.append(dict(level=z,dt=DT[z],absolute_step=n,time=n*DT[z],alpha=1/DT[z],status='NOT_RUN',
                planned_scientific=n%STRIDE[z]==0,actual_online='NOT_RUN',actual_native_full='NOT_RUN',actual_independent='NOT_RUN',qualified_commit=False,
                reason='M prelaunch active-import compatibility failure; monotone stop, no governing retry'))
    table(out/'all-step-timeline.csv',planned)
    audits=[]
    for z in ORDER:
        for n in [START[z],*range(START[z]+STRIDE[z],END[z]+1,STRIDE[z])]:
            audits.append(dict(level=z,step=n,time=n*DT[z],role='startup'if n==START[z]else'completed',status='NOT_RUN',
                reason='observer cache import failed before native startup'if z=='M'and n==START[z]else'governor not launched'))
    table(out/'scheduled-audits.csv',audits)
    table(out/'new-time-pairs.csv',[dict(time=t,C_step=absolute_step('C',t),M_step=absolute_step('M',t),F_step=absolute_step('F',t),status='NOT_RUN',
        physical='NOT_RUN',absolute_sensitivity='NOT_RUN',contraction='NOT_RUN')for t in TIMES])
    bindings=[]
    for z in ORDER:
        r,cp=qualified(selected['anchors'][z],z,START[z]);bindings.append(dict(level=z,step=START[z],time=r['time'],time_bits=bits(r['time']),dt=DT[z],alpha=r['last_stage_alpha'],
            coupled_state_sha256=vector_hash(r['state']),checkpoint=selected['anchors'][z]['checkpoint'],producer=selected['anchors'][z]['report'],
            genuine_callback=selected['anchors'][z]['callback'],native_history_hashes=[vector_hash(v)for v in r['native_history']],native_history_times=r['native_history_times'],
            reference_sha256=vector_hash(r['fixed_reference']),defect_sha256=vector_hash(r['fixed_defect']),original_qh0=r['original_qh0_state_sha256'],
            original_qualification=selected['anchors'][z]['qualified_commit'],new_H23_startup_verification='NOT_RUN',eta=cp['eta']))
        for pin in selected['anchors'][z].values():heavy.append(pin)
    atomic(out/'start-state-bindings.json',bindings)
    preflight_rows=[];q1=[]
    for z,p in preflight['physical'].items():
        preflight_rows.append(dict(level=z,role='passive retained H20 .02, not H23 evolution',summary=p['summary']))
        for region,m in p['all_regional_Q1_moments'].items():
            for i,value in enumerate(m['moments']):q1.append(dict(level=z,time=.02,region=region,test=i,signed_moment=value))
    atomic(out/'retained-02-physical.json',preflight_rows);table(out/'retained-02-signed-Q1.csv',q1)
    fields=[];ratios=[];errors=[]
    for pair,values in preflight['comparison']['pairs'].items():
        for region,data in values['temporal'].items():
            for name,m in data.items():
                for measure in ('rms','maximum'):fields.append(dict(time=.02,pair=pair,region=region,field=name,metric=measure,**m[measure]))
    diagnostic=preflight['comparison']['contraction']
    for kind,group in (('field',diagnostic['fields']),('principal',diagnostic['principal'])):
        for name,regions in group.items():
            for region,data in regions.items():
                for metric in ('RMS','maximum'):ratios.append(dict(time=.02,group=kind,quantity=name,region=region,metric=metric,**data[metric]))
    for z,kinds in preflight['errors'].items():
        for kind,regions in kinds.items():
            for region,values in regions.items():
                for name,m in values.items():errors.append(dict(level=z,time=.02,kind=kind,region=region,field=name,RMS=m['rms']['value'],maximum=m['maximum']['value']))
    table(out/'retained-02-all-pair-screens.csv',fields);table(out/'retained-02-contraction.csv',ratios);table(out/'retained-02-reference-drift-change.csv',errors)
    atomic(out/'retained-02-fixed-owners.json',diagnostic['fixed_probes']);atomic(out/'retained-02-interfaces.json',preflight['interfaces'])
    # Load each immutable NPZ array once. No evaluation, reconstruction or solve.
    arrays=saved_arrays(verify(preflight['signed_arrays']));inventory=load(verify(preflight['inventory']))
    atomic(out/'retained-02-inventory-summary.json',dict(points=len(inventory['points']),owners=len(inventory['owned']),physical_rows=inventory['physical_rows'],
        quadrature_rows=inventory['quadrature_rows'],inventory=preflight['inventory'],signed_arrays=preflight['signed_arrays'],scope='saved retained triplet support only'))
    heavy.extend((preflight['inventory'],preflight['signed_arrays']))
    # Explicitly retain nonzero incoming signed component differences.
    incoming=[];w=arrays['weights']
    for a,b in (('C','M'),('M','F'),('C','F')):
        for name in ('U','G','k','R','pi','production','K','rc','destruction'):
            delta=arrays[a+'_'+name]-arrays[b+'_'+name]
            incoming.append(dict(pair=a+'-'+b,field=name,mean=(np.tensordot(w,delta,axes=1)/w.sum()).tolist(),
                minimum=np.min(delta,axis=0).tolist(),maximum=np.max(delta,axis=0).tolist(),scope='incoming .02 difference, not a new .02->.04 result'))
    atomic(out/'incoming-02-signed-components.json',incoming);del arrays
    for name in ('inputs/selected.json','inputs/historical-indices.json','inputs/historical-immutability.json','source-manifest.json','source-manifest-v1.json','source-manifest-v2.json','execution-config.json','preflight-host.json','admission-M.json'):
        copy(root/name,'controls/'+name.replace('/','--'))
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','target-cache-compatibility.json','setup.json','preflight/result.json','passive-binding-repair-v1/result.json','passive-binding-repair-v1/compatibility.json'):
        copy(root/'science'/name,'controls/'+name.replace('/','--'))
    for z in ORDER:
        for name in ('case.xml','execution-diff.json','intent.json','outcome.json','observer.log','native.log'):
            copy(root/'science'/z/name,z+'/'+name)
    for name in ('request.json','result.json'):copy(root/'science/M/final-strict'/name,'M/final-strict/'+name)
    heavy.append(identity(root/'science/M/final-strict/native.json'))
    for p in sorted((root/'operations').glob('*/receipt.json')):
        copy(p,'operations/'+p.parent.name+'/receipt.json')
        for n in ('request.json','native-exec-samples.json','command.log'):
            if p.parent.name!='stopped-publication-v1':copy(p.parent/n,'operations/'+p.parent.name+'/'+n)
    for p in (root/'passive-repair-v1').glob('*'):copy(p,'source/passive-repair-v1/'+p.name)
    for n in ('bootstrap_wj_h23.py','run_wj_h23.py','study_wj_h23.py','wj_h23_continuation.py','wj_h23_observer.py','test_wj_h23_continuation.py'):
        copy(root/'source-v3/solvers/IncNavierStokesSolver/Tests/RANS'/n,'source/executed-v3/'+n)
    copy(Path(__file__),'source/'+Path(__file__).name)
    # Verify actual executed reader maps. No governor native-exec evidence exists.
    build=load(root/'science/runtime-build-v1.json');expected={v['ELF']['path']:v['ELF']for v in build['executables']}
    deps={v['actual']['path']:v['actual']for v in load(root/'science/runtime-dependencies-v1.json')['libraries']};maps=[];libs={}
    for path in (root/'operations').glob('*/native-exec-samples.json'):
        for v in load(path):
            assert v['status']=='PASS'
            for n in v['native']:verify(expected[n]);assert Path(n).name=='RANSWallResolutionProbe.elf'
            for n in v['files']:
                if '.so'in Path(n).name or Path(n).name.startswith('ld-linux'):libs[n]=identity(verify(deps[n]))
            maps.append(dict(operation=path.parent.name,pid=v['pid'],native=v['native'],libc=v['libc']))
    assert maps
    atomic(out/'runtime-execution.json',dict(status='PASS_READER_ONLY',compiled_source_capsule=build['native_source_revision'],
        static_governor_dependencies='VERIFIED in setup; governor execution NOT_RUN',reader_samples=maps,libraries=list(libs.values())))
    costs=dict(status='NO_NEW_ORDINARY_ADVANCE_OR_AUDIT_MEASUREMENT',governing_samples=0,conditional_H23_forecast='NOT_EVALUABLE',
        actual_operation_usage=usage(book),protected_read_seconds=book['entries']['protected-strict-M8-v1']['seconds'],
        passive_triplet_preflight_seconds=book['entries']['retained-triplet-preflight-v1']['seconds'],
        historical_capture_Off_priors=dict(H20_ordinary_M_seconds=41.35,H20_ordinary_F_seconds=41.86,H20_conditional_days=dict(C=49.9,M=98.2,F=197.3),
            H20_measured_encoding_GB=dict(C=194,M=274,F=429),scope='historical conditional, not H23 measurements or duration admission'),
        horizon=505.02466521278404,terminal_window=[447.8483832973273,505.02466521278404],from_time=.04,
        hypothetical_remaining_from_04={z:dict(remaining_physical_time=505.02466521278404-.04,steps=math.ceil(505.02466521278404/DT[z])-END[z],
            ceiling_absolute_time=math.ceil(505.02466521278404/DT[z])*DT[z],recovery_stride=round(.025/DT[z]),regular_stride=round(.25/DT[z]),terminal_stride=round(.125/DT[z]))for z in DT},
        limitations='No H23 ordinary/warm/monitoring/retention sample; .04 was not reached; cadence/rotation/namespace/compression/production budgets unqualified')
    atomic(out/'conditional-cost-storage.json',costs)
    storage=dict(actual_new_science_bytes=size(root/'science'),source_copies_bytes=sum(size(p)for p in root.iterdir()if p.is_dir()and p.name.startswith(('source','passive-repair'))),
        actual_recoveries=0,planned_recoveries=28,all_recoveries_retained=True,required_recovery_output='NOT_RUN',
        selected_historical_checkpoint_packages={z:{k:v['bytes']for k,v in a.items()}for z,a in selected['anchors'].items()},
        fixed_runtime='reused in place, 0 copied bytes',governor_logs_bytes=(root/'science/M/native.log').stat().st_size,
        compressed_compact_copies=[dict(input_bytes=v['original']['bytes'],output_bytes=v['stored']['bytes'])for v in copies if v['transformation']=='lossless gzip'])
    atomic(out/'storage.json',storage)
    gates=dict(status='OBSTRUCTED_EXECUTION_COMPATIBILITY_BEFORE_GOVERNOR',input_lineage='PASS',retained_02_preflight='PASS_EXACT',
        focused_controls='PASS; earlier selector failure retained',actual_native_trajectory='NOT_RUN',new_startup_parity='NOT_RUN',
        new_physical_qualification='NOT_RUN',new_absolute_sensitivity='NOT_RUN',new_contraction='NOT_RUN',
        protected_M8_read='PASS_ORIGINAL_H20_ANCHOR',passive_binding_repair='PASS_COMPONENT_ONLY',repaired_governor='NOT_RUN',
        all_closed_history_unchanged=True,counts=counts,affordability='NOT_ADMITTED')
    atomic(out/'component-gates.json',gates)
    atomic(out/'heavy-manifest.json',dict(files=heavy));atomic(out/'copy-bindings.json',copies)
    # A signed figure of the incoming state only, with explicit nonexecution labels.
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    probes=[p for p in diagnostic['fixed_probes']if p['name']=='outer_principal']
    fig,axes=plt.subplots(1,len(probes),figsize=(5*len(probes),3.5),squeeze=False,layout='constrained')
    for ax,p in zip(axes[0],probes):
        d=p['principal']['minimized'];ax.bar(['C-M','M-F'],[d['signed_E1'],d['signed_E2']],color=['#276a9f','#bd5837']);ax.axhline(0,color='black',lw=.8)
        ax.set(title=f"Geometry {p['owner']['geometry']} {p['owner']['normal_side']}",ylabel='Signed normalized principal difference at .02')
    fig.suptitle('Retained H20 input: x=15, z=10.2430071555 mm\nH23 continuation NOT_RUN; no .025-.04 result')
    fig.savefig(out/'retained-02-signed-owners.png',dpi=160);plt.close(fig)
    atomic(out/'accounting-publication-snapshot.json',book)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='stopped-prefix publication, completed logs only; closure separately sealed'))
    print(json.dumps(gates),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
