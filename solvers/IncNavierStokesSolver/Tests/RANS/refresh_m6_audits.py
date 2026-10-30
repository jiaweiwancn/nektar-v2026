#!/usr/bin/env python3
"""Re-audit saved M6 fields after an audit-only change, preserving solver evidence.

The solver hash must be unchanged, and every complete-field audit must reproduce
its previous time and field JSON exactly. History publication paths can differ
in the deliberate backup test. Original reports are retained beside refreshed ones.
This does not permit reusing results after a solver or fixture change.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition,message):
    if not condition:raise RuntimeError(message)


def numerical(value):return {key:value[key] for key in ('time','fields')}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--build',type=Path,required=True)
    args=parser.parse_args();require(sys.platform.startswith('linux'),'Use the authoritative Linux environment')
    run=args.run.resolve();build=args.build.resolve();audit=build/'RANSFieldAudit'
    binaries={p.name:sha(p) for p in (build/'IncNavierStokesSolver',audit)}
    roots=[run/name for name in ('final-study','final-verification','final-columns','final-robust','final-standard',
                                'final-outflow-equivalence','final-regressions')]
    reports=[p for root in roots for pattern in ('results*.json','validation.json') for p in root.rglob(pattern)]
    used_cases={p.parent/name/'case.json' for p in reports
                for name in json.loads(p.read_text()).get('cases',{})}
    used_cases.add(run/'final-columns/column_settle/case.json')
    commands=[];checked={};previous={}
    env=os.environ.copy();env.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})

    def evaluate(command,cwd,expected):
        command=[str(audit)]+list(map(str,command[1:]));key=(str(cwd),tuple(command))
        if key in checked:return checked[key]
        result=subprocess.run(command,cwd=cwd,env=env,capture_output=True,text=True,timeout=300)
        require(result.returncode==0,'Re-audit failed: '+str(cwd))
        raw=result.stdout+result.stderr;value=json.loads(raw[raw.index('{"time"'):])
        require(numerical(value)==numerical(expected),'Complete-field audit changed: '+str(cwd))
        log=cwd/('audit-refresh-'+str(len(commands))+'.log');log.write_text(raw)
        commands.append(dict(command=command,cwd=str(cwd),log_sha256=sha(log),complete_field_json_unchanged=True))
        checked[key]=value;return value

    def preserve(path):
        backup=path.with_name(path.name+'.before-audit-refresh')
        if not backup.exists():backup.write_bytes(path.read_bytes())
        return dict(path=str(backup),sha256=sha(backup))

    # Refresh every case cache, including the completed settling column whose
    # suite stopped before publishing a report. Negative solver tests need no audit.
    for root in roots:
        print('Refreshing saved cases in '+str(root),flush=True)
        for path in sorted(root.rglob('case.json')):
            if path not in used_cases:continue
            case=json.loads(path.read_text());fingerprint=case['fingerprint']
            binarykey='binaries' if 'binaries' in fingerprint else 'binary_sha256';old=fingerprint[binarykey]
            require(old['IncNavierStokesSolver']==binaries['IncNavierStokesSolver'],'Solver changed: '+str(path))
            if old==binaries:continue
            directory=path.parent;name=directory.name
            fieldkey='audit' if 'audit' in case else 'field_audit' if 'field_audit' in case else None
            if fieldkey:
                command=[audit,name+'.xml',name+'.fld',name+'.fld']
                if (directory/'profiles.csv').exists():command+=['profiles.csv']
                case[fieldkey]=evaluate(command,directory,case[fieldkey])
            backup=preserve(path);case['audit_refresh']=dict(previous=backup,previous_binaries=old,binaries=binaries)
            fingerprint[binarykey]=binaries
            for key in ('artifacts','artifacts_sha256'):
                for name in list(case.get(key,{})):
                    if (directory/name).is_file():case[key][name]=sha(directory/name)
            path.write_text(json.dumps(case,indent=2)+'\n')

    for path in sorted(reports):
        print('Refreshing report '+str(path),flush=True)
        report=json.loads(path.read_text());key='binaries' if 'binaries' in report else 'binary_sha256';old=report[key]
        require(old['IncNavierStokesSolver']==binaries['IncNavierStokesSolver'],'Solver changed: '+str(path))
        if old==binaries:continue
        # Recheck every comparison command too. The old JSON in its log must
        # agree exactly; solver output and reference-monitor CSVs are untouched.
        for item in report['commands']:
            command=item['command']
            if Path(command[0]).name!='RANSFieldAudit':continue
            cwd=Path(item['cwd']);candidates=[]
            for log in cwd.glob('*.log'):
                if log.name.startswith('audit-refresh-'):continue
                raw=log.read_text(errors='replace')
                if '{"time"' in raw:
                    try:candidates.append(json.loads(raw[raw.index('{"time"'):]))
                    except json.JSONDecodeError:pass
            require(candidates,'Missing previous audit JSON: '+str(cwd))
            result=subprocess.run([str(audit)]+command[1:],cwd=cwd,env=env,capture_output=True,text=True,timeout=300)
            require(result.returncode==0,'Comparison re-audit failed: '+str(cwd))
            value=json.loads(result.stdout[result.stdout.index('{"time"'):])
            require(any(numerical(value)==numerical(v) for v in candidates),'Comparison audit changed: '+str(cwd))
            evaluate(command,cwd,value)
        previous[str(path)]=preserve(path)
        report['audit_refresh']=dict(previous=previous[str(path)],previous_binaries=old,
            source_sha256=sha(Path(__file__).resolve().parents[4]/'library/UnitTests/SolverUtils/RANSFieldAudit.cpp'))
        report[key]=binaries
        for name,case in report.get('cases',{}).items():
            cache=path.parent/name/'case.json'
            if cache.exists():report['cases'][name]=json.loads(cache.read_text())
        path.write_text(json.dumps(report,indent=2)+'\n')
    receipt=dict(status='passed',utc=datetime.now(timezone.utc).isoformat(),binaries=binaries,
        commands=commands,previous_reports=previous,solver_outputs_reused=True)
    (run/'audit-refresh.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(status='passed',complete_field_audits=len(commands),reports=len(previous)),indent=2))


if __name__=='__main__':main()
