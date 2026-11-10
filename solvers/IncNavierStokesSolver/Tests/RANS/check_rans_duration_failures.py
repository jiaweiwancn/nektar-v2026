#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Preserve actual CLI negative fixtures and their intended rejection reasons."""
import argparse,copy,json,subprocess,sys
from pathlib import Path
from test_rans_duration_contract import DurationContract


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=False);fixture=DurationContract();fixture.setUp();good=fixture.good;variants=[('positive',good,None)]
    def changed(name,reason,fn):
        d=copy.deepcopy(good);fn(d);variants.append((name,d,reason))
    changed('missing-time','missing_required_time',lambda d:d['spatial_audits']['B11']['checks'].pop(1))
    changed('wrong-model','configuration_mismatch',lambda d:d.update(configuration=dict(d['configuration'],model='E94E95')))
    changed('wrong-wall','configuration_mismatch',lambda d:d.update(configuration=dict(d['configuration'],wall='EP')))
    changed('missing-space','missing_companion',lambda d:d['companions'].pop('B11'))
    changed('missing-timestep','missing_companion',lambda d:d['companions'].pop('B11-dt2'))
    changed('short-observer','scalar_horizon_incomplete',lambda d:d['companions']['B11-dt2']['scalar_coverage'].update(last_time=.1))
    changed('bad-evolved-reduction','evolved_reduction_failed',lambda d:d['spatial_audits']['B11']['checks'][0]['spatial']['native_full_mass_dual'].update(u=1.1e-3))
    (out/'request.json').write_text(json.dumps(fixture.config));rows=[]
    for name,receipt,reason in variants:
        path=out/(name+'-receipt.json');path.write_text(json.dumps(receipt));dest=out/(name+'-result.json');cmd=[sys.executable,'-B',str(Path(__file__).with_name('classify_rans_duration_bundle.py')),'--request',str(out/'request.json'),'--receipt',str(path),'--output',str(dest)];code=subprocess.run(cmd).returncode;result=json.loads(dest.read_text());matched=(code==0 and result['status']=='pass') if reason is None else code!=0 and any(reason in e for e in result.get('failures',[]));rows.append(dict(name=name,status='pass' if matched else 'fail',actual_exit_code=code,expected_reason=reason,result=result,command=cmd))
    result=dict(status='pass' if all(x['status']=='pass' for x in rows) else 'fail',checks=rows);(out/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
