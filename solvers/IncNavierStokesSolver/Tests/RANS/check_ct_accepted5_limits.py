#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT5 package-wide half-step limit and unchanged numerical classifiers."""
import argparse,ast,json,os,subprocess,tempfile
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted5_run import require_unused_halfstep
from ct_accepted_open import require
from rans_request_manifest import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    path=Path(__file__).with_name('ct_accepted5_run.py');old=subprocess.check_output(['git','-c','index.threads=1','-c','core.preloadIndex=false','show','acafc765:solvers/IncNavierStokesSolver/Tests/RANS/ct_accepted5_run.py'],cwd=path.resolve().parents[4],text=True)
    def classifiers(source):return {x.name:ast.dump(x,include_attributes=False) for x in ast.parse(source).body if isinstance(x,ast.FunctionDef) and x.name in ('basic_gate','traction_gate','physical')}
    require(classifiers(old)==classifiers(path.read_text()),'changed_classification_primitives')
    with tempfile.TemporaryDirectory(prefix='CT5-halfstep-',dir=os.environ['TMPDIR']) as folder:
        root=Path(folder);require_unused_halfstep(root);(root/'initial-quarter').mkdir()
        try:require_unused_halfstep(root)
        except (ValueError,RuntimeError) as e:require('single_conditional_halfstep_already_consumed' in str(e),'wrong_halfstep_rejection')
        else:raise RuntimeError('second_halfstep_admitted')
    out=dict(status='pass',unused_control='admitted',second_control='expected_package_wide_rejection',classifiers='AST-identical to qualified acafc765',source_sha256=sha(path),test_sha256=sha(Path(__file__)),scope='orchestration-only; no native trajectory or physical input changed')
    with a.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print('CT5 package half-step limit pass');return 0


if __name__=='__main__':raise SystemExit(main())
