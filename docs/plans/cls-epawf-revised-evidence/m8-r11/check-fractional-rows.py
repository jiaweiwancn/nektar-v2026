#!/usr/bin/env python3
"""Read-only original-C checks on the retained unfinished Lf/Nf candidates."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode=True


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    assert len(os.sched_getaffinity(0))==1
    e=a.epoch.resolve()
    evidence=e/'source/docs/plans/cls-epawf-revised-evidence'
    h=module(evidence/'m8-r4/prepare.py','r11_fractional')
    C=h.load(e/'comparison-v2.json')
    receipt=h.load(h.bound(C['receipt']))
    manifest=h.load(h.bound(receipt['manifest']))
    sys.path.insert(0,str(e/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
    import numpy as np
    from cls_epawf_m8_r3_resolution import ResolutionAudit
    raw=h.load(h.bound(manifest['current_preflight']))
    audit=ResolutionAudit(raw,e/'source')
    R=np.array(h.vector(C['initial']['residual']))
    J=np.fromfile(h.bound(C['matrix']),dtype='<f8').reshape(len(R),-1)
    constraints={int(i):v for i,v in np.fromfile(h.bound(C['initial']['constraints']),dtype='<f8').reshape(-1,2)}
    results=[]
    for row in C['events']:
        if row['kind']!='r6_cut_audit' or row['status']=='passed_original_cut_scan':
            continue
        d=np.array(h.vector(row['direction']))
        Jd=J@d
        native=R+.125*Jd
        full=R+Jd
        wrong=R/.125+Jd
        full_norm=audit.norms(full,constraints)
        assert max(full_norm[raw['dimension']+2:])<=5e-12
        fractional_error=native-(1.-.125)*R
        fractional_norm=audit.norms(fractional_error,constraints)
        assert max(fractional_norm[raw['dimension']+2:])<=5e-12
        results.append(dict(name=row['name'],unqualified_candidate=True,direction=row['direction'],
            full_algebraic_blocks=full_norm[raw['dimension']+2:].tolist(),
            native_fractional_algebraic_blocks=audit.norms(native,constraints)[raw['dimension']+2:].tolist(),
            fractional_C_identity_error=fractional_norm[raw['dimension']+2:].tolist(),
            incorrectly_scaled_C_blocks=audit.norms(wrong,constraints)[raw['dimension']+2:].tolist()))
    h.save(a.output,dict(comparison=h.pin(e/'comparison-v2.json'),script=h.pin(__file__),
        new_native_trials=0,new_source_values=0,results=results))


if __name__=='__main__':
    main()
