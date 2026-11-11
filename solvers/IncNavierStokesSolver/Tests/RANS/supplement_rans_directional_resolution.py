#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finish the unchanged admitted experiment after a preflight runtime-alias failure."""
import argparse,os,resource,subprocess,sys,time
from pathlib import Path
from qualify_rans_preservation import load
from rans_request_manifest import RequestManifest,runtime_identity,sha
from rans_duration_contract import require,MEMBERS
from rans_qualification import THREADS,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--old-root',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();old=a.old_root.resolve();s=Path(__file__).parent;original=r/'directional-resolution-v1';runtime=r/'runtime'
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'),PYTHONDONTWRITEBYTECODE='1',NEKTAR_RANS_PROBE_TIMEOUT='14400')
    original_request=load(original/'request-manifest.json')['request'];require(RequestManifest(original,original_request,True).completed(original),'directional_source_job_unsealed')
    failure=load(original/'results.json');require(failure.get('failure')=="ValueError('directional_execution_or_required_gate:neutral-audit')",'directional_supplement_wrong_failure')
    require('admitted_response_diagnostic_runtime_changed' in (original/'neutral-audit.log').read_text(),'directional_supplement_wrong_preflight_reason')
    canonical=runtime_identity(runtime/'RANSPreservationProbe');alias=runtime_identity(r/'runtime-resolution-v3/RANSPreservationProbe')
    require(canonical['sha256']==alias['sha256'] and canonical['libraries']==alias['libraries'] and
        (r/'runtime-resolution-v3').resolve()==runtime,'runtime_alias_not_same_bytes_and_target')
    admission=original/'admission.json';ad=load(admission);require(canonical==ad['duration_diagnostic_runtime'],'directional_admitted_runtime_bytes_changed')
    tag='A95-N-Accepted-directional-dt-0.002';base=[tag+'-neutral-'+m for m in MEMBERS];pert=[tag+'-response-'+m for m in MEMBERS]
    native_requests={label:sha(r/'campaigns'/label/'request-manifest.json') for label in (*base,*pert)}
    out=r/'directional-resolution-supplement-v1'
    intent=dict(schema='DirectionalResolutionExecutionSupplementV1',original_failed_job_sha256=sha(original/'results.json'),
        original_request_sha256=sha(original/'request-manifest.json'),admission_sha256=sha(admission),
        native_requests=native_requests,time=.016,amplitude=2e-5,timesteps=[.002,.001],runtime=canonical,
        alias_identity=alias,fix='canonicalize the same executable directory before strict identity checking; no runtime, boundary, state, method or gate changes',
        scope='reuse four sealed neutral endpoints and eight initialized bindings; continue only the still-zero-step perturbed cases')
    seal=RequestManifest(out,intent);events=[];failures=[]
    def run(name,cmd,required=True):
        before=resource.getrusage(resource.RUSAGE_CHILDREN);start=time.monotonic()
        with (out/(name+'.log')).open('x') as stream:code=subprocess.run(list(map(str,cmd)),stdout=stream,stderr=subprocess.STDOUT).returncode
        after=resource.getrusage(resource.RUSAGE_CHILDREN);events.append(dict(name=name,command=list(map(str,cmd)),exit_code=code,
            wall_seconds=time.monotonic()-start,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
        if code:failures.append(name)
        write_json(out/'progress.json',dict(events=events,failures=failures))
        if required:require(code==0,'directional_supplement_required_gate:'+name)
        return code
    try:
        # This wrapper records the failed CLI preflight; the original result and
        # log remain unchanged. The existing audit supplement validates these
        # four source identities before gathering any independent native audit.
        previous=r/'response-audits'/(tag+'-neutral-0.016');old_request=load(previous/'request-manifest.json')['request']
        wrapper=out/'failed-preflight';wrapper_request=dict(old_request,schema='AuditRuntimeAliasExecutionFailureV1',
            original_failed_manifest_sha256=sha(previous/'request-manifest.json'),source_log_sha256=sha(original/'neutral-audit.log'))
        wrapped=RequestManifest(wrapper,wrapper_request)
        write_json(wrapper/'results.json',dict(status='fail',configuration=ad['configuration'],time=.016,
            failure='runtime alias path mismatch before native audit',performed_native_audits=0))
        wrapped.seal(wrapper,['results.json'])
        run('neutral-audit',[sys.executable,'-B',s/'audit_rans_response.py','--root',r,'--admission',admission,'--cohort-prefix',tag+'-neutral','--time','.016',
            '--receipt-tag','canonical-v1','--supplement-from',wrapper/'results.json'])
        binding=old/'physical-budget-bindings/A95-N-Accepted-response-2e-05-dt-0.002-0.016/results.json'
        for i,label in enumerate(pert):
            folder=r/'campaigns'/label;q=load(folder/'request-manifest.json')['request']
            require(RequestManifest(folder,q,True).completed(folder/'actual-initial'),'supplement_actual_initial_not_sealed')
            control='StreamwiseBaseV1' if i in (0,2) else 'StreamwiseP6V1'
            run('evolve-'+label,[sys.executable,'-B',s/'run_rans_budget_axis_control.py','--root',r,'--binding',binding,
                '--hypothesis',original/(control+'-hypothesis.json'),'--mesh-file',original/(control+'-mesh.json'),
                '--control',control,'--label',label,'--dt','.001' if i>=2 else '.002','--amplitude','2e-5','--resume'])
            require(sha(folder/'request-manifest.json')==native_requests[label],'supplement_mutated_native_request')
        run('response-comparison',[sys.executable,'-B',s/'compare_rans_nonlinear_response.py','--root',r,'--perturbed',*pert,'--baseline',*base,
            '--time','.016','--output-directory',out/'response-comparison'],False)
        run('response-audit',[sys.executable,'-B',s/'audit_rans_response.py','--root',r,'--admission',admission,'--cohort-prefix',tag+'-response','--time','.016'],False)
        for kind,labels in (('baseline',base),('perturbed',pert)):
            for i,(member,label) in enumerate(zip(MEMBERS,labels,strict=True)):
                for q,face in (((1,16),(2,16),(2,32)) if i<2 else ((2,32),)):
                    output=tag+'-'+kind+'-'+member+'-0.016-Q'+str(q)+('-F32' if face==32 else '')
                    run('budget-'+kind+'-'+member+'-'+str(q)+'-'+str(face),[sys.executable,'-B',s/'audit_rans_physical_scalar_budget.py',
                        '--root',r,'--label',label,'--time','.016','--output-label',output,'--quadrature-index',q,'--face-order',face,'--runtime-name','runtime'],False)
        analysis=out/'budget-analysis'
        run('budget-analysis',[sys.executable,'-B',s/'analyze_rans_physical_budgets.py','--root',r,'--cohort',tag,'--time','.016','--output',analysis],False)
        run('budget-accuracy',[sys.executable,'-B',s/'classify_rans_physical_budget_accuracy.py','--source',analysis/'results.json','--output',out/'budget-accuracy'],False)
        run('budget-binding',[sys.executable,'-B',s/'bind_rans_ownership_budget_accuracy.py','--root',r,'--admission',admission,
            '--cohort',tag+'-response','--baseline-cohort',tag+'-neutral','--dt','.002','--time','.016','--analysis',analysis/'results.json',
            '--accuracy',out/'budget-accuracy/results.json','--output',out/'budget-binding'],False)
        result=dict(status='pass' if not failures else 'failed_or_unresolved',events=events,failures=failures,
            configuration=ad['configuration'],time=.016,original_failed_receipt_unchanged=True,
            response_comparison=load(out/'response-comparison/results.json'),physical_budget_analysis=load(analysis/'results.json'),
            physical_budget_accuracy=load(out/'budget-accuracy/results.json'),physical_budget_binding=load(out/'budget-binding/results.json'),
            response_extension_eligible=False,scope=intent['scope'])
    except Exception as error:result=dict(status='failed_or_unresolved',failure=repr(error),events=events,failures=failures,
        original_failed_receipt_unchanged=True,response_extension_eligible=False)
    write_json(out/'results.json',result);seal.seal(out,[p for p in out.rglob('*') if p.is_file() and p.name not in ('request-manifest.json','completed-artifacts.json')])
    print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
