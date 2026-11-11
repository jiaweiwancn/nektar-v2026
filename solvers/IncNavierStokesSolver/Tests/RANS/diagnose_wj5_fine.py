#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only context from retained WJ-4 L2/Nx1 at Tx; not an hx gate."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--frozen-source',type=Path,required=True);p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--retained',type=Path,required=True);p.add_argument('--policy',choices=('EP','Budget'),required=True);a=p.parse_args();a.mode='report'
    frozen=a.frozen_source.resolve()/'solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(frozen))
    from run_wj5_native import Study,read_execution_intent
    from wj2_cases import normal_knots
    from wj_immutable import artifacts,verify_artifacts,read,file_hash,write_once,digest
    from run_cls_epawf_picard import require
    own=Path(__file__).resolve().parent;spec=importlib.util.spec_from_file_location('wj5_fine_estimates',own/'wj5_measure.py');estimates=importlib.util.module_from_spec(spec);spec.loader.exec_module(estimates)
    a.solver=a.solver.resolve();a.probe=a.probe.resolve();s=Study(a)
    try:
        folder=a.retained.resolve()/'study'/a.policy/'longer-coarse-32';intent=read_execution_intent(folder/'execution-intent.json');receipt=read(folder/'receipt.json')
        require(digest(intent)==receipt['intent_sha256'],'retained fine intent changed');verify_artifacts(intent['inputs']);verify_artifacts(receipt['artifacts'])
        for code in intent['code'].values():verify_artifacts(code)
        for identity in intent['runtime'].values():
            require(file_hash(identity['path'])==identity['sha256'],'retained fine runtime changed')
            for p,h in identity['libraries'].items():require(file_hash(p)==h,'retained fine dependency changed')
        d=read(folder/'case.evolution.json');require(len(d['state'])==2493 and d['steps']==32 and abs(d['time']-.001024)<2e-15,'retained fine shape/time differs')
        code=artifacts([Path(__file__),own/'wj5_measure.py']);s.measure_code=dict(s.measure_code,fine_diagnostic=code);s.measure_version=digest(s.measure_code)[:12]
        knots=normal_knots(2)
        write_once(s.root/('fine-diagnostic-design-'+s.measure_version+'.json'),dict(schema='WJ5FineContextV1',role='read-only retained uniform L2/Nx1 context; joint x/z comparison, not pure hx qualification',
            actual_zknots=knots,actual_xknots=[0.,1.],report=artifacts([folder/'case.evolution.json',folder/'receipt.json',folder/'execution-intent.json']),code=code,governing_attempts=0))
        # This extra key is only an observation lookup. The immutable governing
        # mesh design is unchanged; no native case is generated from this key.
        s.design['meshes']['L2_diagnostic']=knots
        m,c=s.measure('retained-L2Nx1-32',d,folder,mesh='L2_diagnostic',nx=1,assembly=True,order=8,xorder=4)
        c['arithmetic']=s.arithmetic('retained-L2Nx1-32',d,m);c['boundary_uncertainty']=estimates.boundary_uncertainty(m,c['arithmetic'],c['rule_settings'])
        c['comparison_role']='joint normal/horizontal diagnostic context; not selected-mesh preservation or horizontal qualification'
        s.s.setdefault('retained_fine',{})['retained-L2Nx1-32']=c;s.finish()
    finally:s.ledger.finish_session(s.session);s.save()

if __name__=='__main__':main()
