#!/usr/bin/env python3
"""Bind the full M6 evidence to its localized physical-orientation guard fix.

This is a narrow, mechanically checked source transition, not a general
exception to the strict final documentation-only runtime fingerprint.
"""
import argparse,hashlib,json,os,shutil,subprocess
from datetime import datetime,timezone
from pathlib import Path
os.environ['PYTHONDONTWRITEBYTECODE']='1'
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
p.add_argument('--original-audit',required=True);p.add_argument('--guard',required=True)
p.add_argument('--guard-fingerprint',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=a.root;source=root/'source';original=root/a.original_audit;guard=root/a.guard
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
old=read(original/'integration.json');ori=read(guard/'qualification/validation.json');gr=read(guard/'receipt.json')
before=read(original/'qualification-fingerprint.json');after=read(root/a.guard_fingerprint)
assert old['status']=='verified_before_orientation_guard' and ori['success'] and gr['status']=='passed'
assert ori['source_clean_start'] and ori['source_clean_end']
base=old['source_commit'];revision=gr['source_commit']
assert ori['source_commit']==ori['source_commit_end']==revision==after['source_commit']==git('rev-parse','HEAD')
assert gr['source_end_commit']==revision and not gr['source_end_status'] and not git('status','--porcelain')
assert after['success'] and after['clean'] and before['source_commit']==base
assert gr['main_commit']==gr['main_end_commit']==old['preserved_linux_main_commit'] and not gr['main_end_status']
assert len(ori['constructors'])==6 and sum(c['expected_accept'] for c in ori['constructors'])==4
assert len(ori['operator_replays'])==3 and all(r['original_mathematical_fields_unchanged'] for r in ori['operator_replays'])
for replay in ori['operator_replays']:
 assert sha(Path(replay['original_report']))==replay['original_report_sha256']
 assert sha(Path(replay['replay_report']))==replay['replay_report_sha256']
changed=[name for name in set(before['source_inputs'])|set(after['source_inputs']) if before['source_inputs'].get(name)!=after['source_inputs'].get(name)]
production='solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionTransition.cpp'
probe='library/UnitTests/SolverUtils/RANSCurvedWallProbe.cpp'
cmake='solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt'
test='solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m6_orientation.py'
assert set(changed)=={production,probe,cmake,test},changed
def prior(name):return subprocess.check_output(['git','-C',str(source),'show',base+':'+name],text=True)
old_guard='    Require(dx>0,"the prescribed upstream coordinate is oriented along positive q1");\n'
new_guard=('    const auto incoming=m_mapping.At({lower[0],lower[1],0});\n'
 '    const double incomingX=dx*incoming.tangent1[0]+dy*incoming.tangent2[0];\n'
 '    Require(incomingX>0,"the physical incoming flow must increase the prescribed wall coordinate");\n')
assert prior(production).count(old_guard)==1
assert (source/production).read_text()==prior(production).replace(old_guard,new_guard)
anchor='        const auto mapping=LoadRANSWallNormalMap(serial,d,parameters.geometry.height);RANSWallEvolutionPatch patch(serial,graph,world->GetComm(),parameters,mapping);\n'
following='        const double K=parameters.model.Q()*.16,H=parameters.geometry.height,Z=patch.upper[2],time=.2,alpha=5.3;'
text=(source/probe).read_text();start=text.index(anchor)+len(anchor);end=text.index(following,start)
instrumentation=text[start:end]
assert instrumentation.startswith('        int constructorOnly=0;') and 'if(constructorOnly)' in instrumentation
assert text[:start]+text[end:]==prior(probe)
text=(source/cmake).read_text();start=text.index('                add_test(NAME RANSCLSEPAWFM6OrientationOracle')
end=text.index('                add_test(NAME RANSCLSEPAWFM5NoteOracle',start)
registration=text[start:end]
expected=('                add_test(NAME RANSCLSEPAWFM6OrientationOracle COMMAND ${Python3_EXECUTABLE}\n'
 '                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m6_orientation.py\n'
 '                    --build $<TARGET_FILE_DIR:RANSCurvedWallProbe>\n'
 '                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m6-orientation)\n'
 '                set_tests_properties(RANSCLSEPAWFM6OrientationOracle PROPERTIES PROCESSORS 1 TIMEOUT 1800\n'
 '                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")\n')
assert registration==expected and text[:start]+text[end:]==prior(cmake)
assert test not in before['source_inputs']
runtime_changed=[name for name in set(before['runtime_artifacts'])|set(after['runtime_artifacts']) if before['runtime_artifacts'].get(name)!=after['runtime_artifacts'].get(name)]
allowed={'build/IncNavierStokesSolver','build/RANSEquilibriumEvolutionProbe','build/RANSCurvedWallProbe'}
assert set(runtime_changed)<=allowed,runtime_changed
assert all(before[key]==after[key] for key in ('runtime_bindings','runtime_environment','python_packages','tool_images'))
names={r['name'] for r in json.loads(subprocess.check_output(['ctest','--test-dir',str(root/'build'),'--show-only=json-v1'],text=True))['tests']}
coverage=dict(old['registered_test_coverage'])
coverage['RANSCLSEPAWFM6OrientationOracle']={'status':'passed','execution':'direct registered orientation entry point, plus three complete operator replays',
 'source_commit':revision,'receipt':str(guard/'qualification/validation.json'),'checks':ori['checks']}
assert set(coverage)==names and len(names)==31
receipt=dict(old,status='closed',source_commit=revision,full_qualification_commit=base,
 collected_utc=datetime.now(timezone.utc).isoformat(),registered_test_coverage=coverage,total_registered_tests=31,
 orientation_checks=ori['checks'],orientation_source_transition={'from':base,'to':revision,
  'changed_numerical_inputs':sorted(changed),'changed_runtime_images':sorted(runtime_changed),
  'mechanically_verified':'only physical incoming-orientation validation, optional constructor-only probe, and its registered oracle added',
  'equations_time_integrator_and_fixed_mode_paths_unchanged':True,
  'operator_replays':ori['operator_replays'],
  'reuse_scope':'full prior qualification on its admitted flat/cosine cases; all 28 earlier regressions use Fixed mode and do not call the changed guard'},
 fingerprint=str(root/a.guard_fingerprint),fingerprint_sha256=sha(root/a.guard_fingerprint))
a.output.mkdir(parents=True,exist_ok=False)
for path in original.iterdir():
 if path.is_file():shutil.copyfile(path,a.output/path.name)
shutil.copyfile(original/'integration.json',a.output/'pre-guard-integration.json')
shutil.copyfile(root/a.guard_fingerprint,a.output/'guard-fingerprint.json')
shutil.copyfile(guard/'qualification/validation.json',a.output/'orientation-verification.json')
shutil.copyfile(guard/'receipt.json',a.output/'orientation-receipt.json')
index=read(original/'raw-artifact-index.json')
for path in sorted(guard.rglob('*')):
 if path.is_file():index.append({'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
(a.output/'raw-artifact-index.json').write_text(json.dumps(index,indent=2)+'\n')
(a.output/'integration.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'closed','source_commit':revision,'full_qualification_commit':base,'registered_tests':31,
                  'orientation_checks':ori['checks'],'replayed_operators':len(ori['operator_replays']),'raw_artifacts':len(index)}))
