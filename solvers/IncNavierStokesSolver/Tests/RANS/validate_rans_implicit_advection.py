#!/usr/bin/env python3
"""Linux implicit-advection transport/history verification, separate from fitting."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
import subprocess
from pathlib import Path
import sys
from threading import Lock
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
from cao_tamura_experiment import ROOT, environment, execute, allocation, records, sha, write_json
from generate_m3 import generate
from generate_m4 import generate_cls, info, parameter


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpus',default='9,10,11,12');p.add_argument('--resume',action='store_true');args=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Authoritative validation uses Linux')
    build=environment(args.build);out=args.output.resolve();out.mkdir(parents=True,exist_ok=args.resume)
    cpus=args.cpus.split(',');allocation(4,args.cpus)
    previous=json.loads((out/'results.json').read_text()) if args.resume and (out/'results.json').exists() else {}
    result=dict(source_revision=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        solver_sha256=sha(build/'IncNavierStokesSolver'),threads={k:os.environ[k] for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS')},
        tests={},maximum_workers=4,source_hashes={str(Path(__file__).relative_to(ROOT)):sha(__file__)})
    def persist():write_json(out/'results.json',result)
    # SymPy's process-global assumption caches are not thread-safe. Solver
    # executions may run concurrently; symbolic case construction must not.
    generator_lock=Lock()
    def solve(name,model,dt,steps,cpu,modify=None,ranks=1,expected_time=None):
        folder=out/name;folder.mkdir(exist_ok=args.resume);case=folder/'case.xml'
        generator=generate if model=='standard' else generate_cls
        with generator_lock:
            generator(case,kind='mms',order=2,dt=dt,steps=steps,modes=9 if model=='standard' else 7,points=14)
        tree=ET.parse(case);root=tree.getroot()
        info(root,'RANSLinearSolver','MatrixFree')
        info(root,'RANSAdvectionTreatment','Implicit')
        info(root,'GlobalSysSoln','DirectStaticCond' if ranks==1 else 'IterativeStaticCond')
        for key,value in dict(RANSNonlinearTolerance=1e-12,RANSPressureTolerance=1e-12,
            RANSMaxIterations=10000,IterativeSolverTolerance=1e-13,IO_CheckSteps=steps).items():parameter(root,key,value)
        if modify:modify(root)
        tree.write(case,encoding='utf-8',xml_declaration=True)
        saved=folder/'verification.json'
        if args.resume and saved.exists():
            cached=json.loads(saved.read_text())
            if cached['case_sha256']==sha(case) and previous.get('solver_sha256')==result['solver_sha256']:
                assert abs(cached['audit']['time']-(expected_time if expected_time is not None else dt*steps))<1e-11
                return dict(folder=str(folder),audit=cached['audit'])
        receipt=execute(allocation(ranks,cpu)+[build/'IncNavierStokesSolver',case.name,'--use-scotch'],folder,'solver',900)
        stages=records(folder/'case.rans.csv')
        assert stages and abs(float(stages[-1]['time'])-(expected_time if expected_time is not None else dt*steps))<1e-11
        assert all(float(r['min_k'])>0 and float(r['min_eps'])>0 and float(r['guard_activations'])==0 and float(r['weak_divergence_l2'])<1e-7 for r in stages)
        execute(allocation(1,cpu.split(',')[0])+[build/'RANSFieldAudit',case,folder/'case.fld',folder/'case.fld'],folder,'audit')
        text=(folder/'audit.log').read_text();audit=json.loads(text[text.index('{"time"'):])
        write_json(folder/'verification.json',dict(receipt=receipt,case_sha256=sha(case),audit=audit,stages=stages))
        return dict(folder=str(folder),audit=audit)
    for model in ('standard','cls'):
        dts=(.05,.025,.0125,.00625)
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures=[pool.submit(solve,f'{model}-dt{dt:g}',model,dt,round(.2/dt),cpu) for dt,cpu in zip(dts,cpus)]
            runs=[f.result() for f in futures]
        fields=('u','v','k','eps' if model=='standard' else 'eps_tilde','p')
        rates={f:[math.log(a['audit']['fields'][f]['l2_error']/b['audit']['fields'][f]['l2_error'],2) for a,b in zip(runs,runs[1:])] for f in fields}
        result['tests'][model+'-temporal']=dict(runs=runs,rates=rates)
        persist();print(model,'temporal',rates,flush=True)
        # Report all rates; a velocity/scalar failure rejects second-order use.
        if not all(min(rates[f][-2:])>1.75 for f in fields if f!='p'):
            raise RuntimeError('Implicit-advection transport order verification failed')
        full=Path(runs[1]['folder'])
        for prefix in (1,2,5):
            name=f'{model}-prefix{prefix}'
            part=solve(name,model,.025,prefix,cpus[0])
            def restart(root):
                group=root.find("CONDITIONS/FUNCTION[@NAME='InitialConditions']")
                group.clear();group.set('NAME','InitialConditions')
                for v in root.findall('CONDITIONS/VARIABLES/V'):
                    ET.SubElement(group,'F',VAR=v.text,FILE=str(Path(part['folder'])/'case.fld'))
                info(root,'RANSRestartMode','History')
            resumed=solve(f'{model}-restart{prefix}',model,.025,8-prefix,args.cpus,restart,4,expected_time=.2)
            folder=Path(resumed['folder'])
            execute(allocation(1,cpus[0])+[build/'RANSFieldAudit',full/'case.xml',full/'case.fld',folder/'case.fld'],folder,'compare')
            text=(folder/'compare.log').read_text();audit=json.loads(text[text.index('{"time"'):])
            error=max(v['difference_linf']/max(1,abs(v['min']),abs(v['max'])) for v in audit['fields'].values())
            result['tests'][f'{model}-restart{prefix}']=dict(maximum_scaled_difference=error,audit=audit,ranks_prefix=1,ranks_continuation=4)
            persist();assert error<1e-6
    result['accepted']=True;persist();print('Implicit-advection transport/history checks PASS',flush=True)


if __name__=='__main__':main()
