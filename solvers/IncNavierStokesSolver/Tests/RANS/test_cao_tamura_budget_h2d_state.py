#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Signed-zero parser and native immutable-row rejection, no physical advance."""
import argparse,json,struct,os
from pathlib import Path
from cao_tamura_h2d_state import load,rebind_fixed,bits
from run_cao_tamura_budget_h2d import case
from diagnose_cls_accepted_covariance import run_probe
from continue_cls_epawf_pimple import atomic,freeze
from run_cls_epawf_picard import THREADS,sha,require

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','probe','preflight','h1'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);out=a.root/'signed-contract-v1';out.mkdir(exist_ok=False)
    f=out/'parser.json';f.write_text('{"state":[-0,0,-0.0,1]}')
    v=load(f)['state'];require(bits(v[0])==bits(v[2])!=bits(v[1]),'integer -0 was normalized')
    base=dict(feasibility=dict(base=[-0.,1.,2.],free_velocity_rows=[1]));c=dict(state=[0.,3.,2.]);restored,rows=rebind_fixed(base,c)
    require(rows==[0] and bits(restored[0])==bits(-0.) and restored[1]==3.,'signed-zero rebind changed allowed values')
    try:rebind_fixed(base,dict(state=[.1,3.,2.]));raise RuntimeError('nonzero frozen-row error accepted')
    except ValueError:pass
    probe,runtime=freeze(a.probe,out/'runtime');xml=case(out/'case',a.h1,load(a.preflight));old=load(a.root/'search-v1/best-candidate.json')
    text,receipt=run_probe(probe,xml,dict(operation='h2d_candidate',inner_stress_model='Linear',xml_sha256=sha(xml),state=old['state'],state_sha256=old['state_sha256'],base_sha256=old['base_sha256']),out/'rejected-native',success=False)
    require('changed scalar or pressure bits' in text,'native did not reject numerically equal changed bits')
    atomic(out/'progress.json',dict(status='passed',runtime=runtime,receipt=receipt,tests=dict(parser_preserves_integer_negative_zero=True,rebind_frozen_bits_only=True,fixed_numeric_change_rejected=True,native_fixed_signed_bits_rejection=True),governing_solves=0,physical_steps=0))
if __name__=='__main__':main()
