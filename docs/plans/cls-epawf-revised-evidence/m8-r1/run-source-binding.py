"""Bind committed final helpers/current C++ to the frozen tested build.

No rebuild or trajectory. Historical and controlled source/build worktrees
stay fixed; a separate detached checkout holds the committed audit helpers.
Added evidence/scripts legitimately change an all-input manifest, so this
receipt compares compiled inputs explicitly rather than claiming otherwise.
"""
import hashlib, json, os, subprocess
from pathlib import Path

r1 = Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r1')
root = r1/'controlled'; frozen = root/'source'; build = root/'build'
out = r1/'source-binding-001'; out.mkdir(exist_ok=False)
main = Path('/media/jay/新加卷/nektar-v2026')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
git = lambda p,*args: subprocess.check_output(['git','-C',str(p),*args],text=True).strip()
assert git(main,'rev-parse','HEAD') == 'f4b70baec10016dd098005bc5a5507b87e4a57bf'
assert not git(main,'status','--porcelain') and not git(frozen,'status','--porcelain')
assert git(frozen,'rev-parse','HEAD') == '360cbb7ac728b575cdf9af21e6d0963a1d88acb9'
bundle = r1/'helpers.bundle'
subprocess.run(['git','-C',str(main),'fetch','--no-tags',str(bundle),'HEAD'],check=True)
revision = git(main,'rev-parse','FETCH_HEAD')
assert revision.startswith('f53bc61')
source = out/'source'
subprocess.run(['git','-C',str(main),'worktree','add','--detach',str(source),revision],check=True)
compiled = []
for name in git(frozen,'ls-files').splitlines():
    p = Path(name)
    if p.suffix.lower() in ('.c','.cc','.cpp','.cxx','.h','.hpp','.hxx','.cmake') or p.name == 'CMakeLists.txt':
        compiled.append(name)
assert compiled
differences = [name for name in compiled if not (source/name).is_file() or sha(source/name) != sha(frozen/name)]
assert not differences, differences
added_compiled = [name for name in git(source,'diff','--name-only',git(frozen,'rev-parse','HEAD'),revision).splitlines()
    if Path(name).suffix.lower() in ('.c','.cc','.cpp','.cxx','.h','.hpp','.hxx','.cmake') or Path(name).name == 'CMakeLists.txt']
assert not added_compiled, added_compiled
helpers = {str(p.relative_to(source)):sha(p) for p in sorted(
    (source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('cls_epawf_m8_r1_*.py'))}
assert len(helpers)==6
pin = root/'controls-001/runtime-before.json'; before = json.loads(pin.read_text())
env = dict(os.environ)
for key,value in before['runtime_environment'].items():
    if value is None: env.pop(key,None)
    else: env[key]=value
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS',
            'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'): assert env[key]=='1'
argv = ['python3','-B',str(frozen/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'),
    '--root',str(root),'--output',str(out/'runtime-final.json'),'--compare',str(pin)]
with (out/'runtime.log').open('w') as stream:
    code=subprocess.run(argv,cwd=out,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=240).returncode
assert code==0
receipt=dict(schema=1,status='passed',runner_sha256=sha(__file__),bundle_sha256=sha(bundle),
    helper_source_commit=revision,helper_source_clean=not bool(git(source,'status','--porcelain')),
    compiled_source_commit=git(frozen,'rev-parse','HEAD'),compiled_source_clean=True,
    compiled_input_count=len(compiled),compiled_inputs_byte_identical=True,
    compiled_manifest_sha256=hashlib.sha256(json.dumps({name:sha(source/name) for name in compiled},sort_keys=True,separators=(',',':')).encode()).hexdigest(),
    helpers=helpers,mathematical_residual_changed=False,
    qualification='No new C++ law or production solver policy; final helpers differ from historical epoch by explicit new files.',
    all_input_manifest_identity_claimed=False,
    actual_runtime_pin=dict(path=str(out/'runtime-final.json'),sha256=sha(out/'runtime-final.json')),
    compared_runtime_pin=dict(path=str(pin),sha256=sha(pin)),runtime_command=argv,
    runtime_log_sha256=sha(out/'runtime.log'),runtime_exit_code=code,
    main_checkout_preserved=git(main,'rev-parse','HEAD')=='f4b70baec10016dd098005bc5a5507b87e4a57bf' and not git(main,'status','--porcelain'))
assert receipt['helper_source_clean'] and receipt['main_checkout_preserved']
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(status=receipt['status'],compiled_inputs=len(compiled),helpers=len(helpers),source_commit=revision)))
