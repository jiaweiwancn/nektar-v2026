from pathlib import Path
import json,hashlib
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/reference-inverse-v19';src=r/'source-v19';prefix='solvers/IncNavierStokesSolver/EquationSystems/';lib='library/SolverUtils/Turbulence/'
common=[lib+n+s for n in('WJEARSM','RANSModels','AnalyticalWallFunction')for s in('.h','.cpp')]
groups={'mixed_BC_gauge_Q3_force':(8,common+[prefix+n+s for n in('WJNativeFlow','WJNativeMixedFlow')for s in('.h','.cpp')]+[prefix+'WJNativeFailures.h']),
 'constitutive_and_planar_admissibility':(7,common+[prefix+n+s for n in('WJNativeFlow','WJNativeWallV2','WJNativeConstitutive','WJNativeAdmissibility')for s in('.h','.cpp')])}
matched={}
for name,(version,paths)in groups.items():
 before=r/('source-v'+str(version));a=json.loads((r/('source-v'+str(version)+'.json')).read_text());b=json.loads((r/'source-v19.json').read_text());am={v['path']:v for v in a['files']};bm={v['path']:v for v in b['files']};rows=[]
 for path in paths:
  old,new=identity(before/path),identity(src/path);assert old['sha256']==new['sha256']==am[str(before/path)]['sha256']==bm[str(src/path)]['sha256'],path
  rows.append(dict(relative_path=path,old=old,current=new))
 matched[name]=dict(source_version=version,old_manifest=identity(r/('source-v'+str(version)+'.json')),current_manifest=identity(r/'source-v19.json'),relevant_numerical_sources_match=True,files=rows)
reports={}
for case in('closed','slip','physicaloutlet'):
 p=r/('science/mixed-'+case+'-oracle-v8.json');d=json.loads(p.read_text());raw=Path(d['input']);assert d['input_sha256']==identity(raw)['sha256']and d['status']=='PASS_TINY_MIXED_BC_Q3_FOUNDATION_ONLY'
 reports[case]=dict(report=identity(p),native=identity(raw),maximum_operator_error=max(v['relative_error']for v in d['checks'].values()),gauge=d['gauge'],unaugmented_null_modes=d['unaugmented_null_modes'],native_v_map=d['native_v_map'],direct_solution_relative_error=d['direct_solution_relative_error'],known_scalar=d['known_scalar'],mutations=d['mutations'])
for name,rawname in [('admissibility-independent-v2','admissibility-native-v2'),('constitutive-independent-v1','constitutive-native-v1')]:
 p=r/'science'/(name+'.json');raw=r/'science'/(rawname+'.json');d=json.loads(p.read_text());n=json.loads(raw.read_text());assert d['status']=='PASS'and d['native_sha256']==identity(raw)['sha256'];reports[name]=dict(report=identity(p),native=identity(raw),states=len(d['states']))
 if name.startswith('admissibility'):reports[name].update(negative_covariance_count=n['negative_covariance_count'],negative_principal_count=n['negative_principal_count'],unsupported_spanwise_rejected=n['unsupported_spanwise_rejected'])
 else:reports[name].update({k:d[k]for k in('maximum_value_parity_error','maximum_direction_parity_error','maximum_joint_FD_error','value_fallbacks','direction_fallbacks')})
new=json.loads((o/'independent.json').read_text());assert new['inverse_guard_gate_pass']is True
atomic(o/'reuse-audit.json',dict(schema='WJH27VerificationReuseAuditV1',status='PASS_REUSED_NUMERICAL_FOUNDATIONS_ONLY',source_identity=json.loads((r/'source-v19.json').read_text())['expected_source_identity'],groups=matched,legacy_reports=reports,new_source19_inverse_guard_report=identity(o/'independent.json'),
 full_H27_qualification=False,TASK3_READY=False,scope='Hash-bound reuse of unchanged numerical implementations and their existing authoritative Linux tests. New source19 inverse guards checked separately. No transfer of kernel timings, physical mesh/dt, long pilots, production retry/statistics/storage or complete cost.'))
print('SOURCE19_REUSE_AUDIT_PASSED')
