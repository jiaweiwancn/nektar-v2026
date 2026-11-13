#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded H2e setup and independent native certification (Linux only)."""
import argparse, json, os, platform, time, resource
from pathlib import Path
from cao_tamura_h2d_state import load
from run_cao_tamura_budget_h2d import case
from run_cls_epawf_picard import THREADS,sha,require
from continue_cls_epawf_pimple import atomic,freeze
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import state_hash

def context(a):
    import numpy as np
    a.output.mkdir(parents=True,exist_ok=False)
    exe,runtime=freeze(a.probe.resolve(),a.output/'runtime')
    pre=load(a.preflight);xml=case(a.output/'case',a.h1,pre)
    # Only decode the bounded immutable prefix; never materialize 789 MB JSON.
    with a.old_context.open() as f:prefix=f.read(4*1024*1024)
    decoder=json.JSONDecoder(parse_int=lambda t:-0.0 if t=='-0' else int(t))
    def field(name):return decoder.raw_decode(prefix,prefix.index('"'+name+'":')+len(name)+3)[0]
    base=field('base');basis=np.array(field('basis'),dtype='<f8');require(basis.shape==(32,len(base)),'saved basis dimensions')
    folder=a.output/'arrays';folder.mkdir();basis.tofile(a.output/'old-basis.f64')
    p=dict(operation='h2e_context',inner_stress_model='Linear',xml_sha256=sha(xml),base_sha256=state_hash(base),basis_file=str(a.output/'old-basis.f64'),basis_sha256=sha(a.output/'old-basis.f64'),context_directory=str(folder))
    result,receipt=run_probe(exe,xml,p,a.output/'native-context')
    result=load(Path(receipt['output']));result.update(runtime=runtime,receipt=receipt,xml=str(xml),xml_sha256=sha(xml),array_directory=str(folder),old_context=str(a.old_context),preflight=str(a.preflight),preflight_sha256=sha(a.preflight),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    atomic(a.output/'context.json',result)

def certify(a):
    a.output.mkdir(parents=True,exist_ok=False)
    d=load(a.context);c=load(a.candidate);xml=Path(d['xml']);exe=Path(d['runtime']['path'])
    if a.probe:exe,d['runtime']=freeze(a.probe.resolve(),a.output/'runtime')
    p=dict(operation='h2e_candidate',inner_stress_model='Linear',xml_sha256=sha(xml),state=c['state'],state_sha256=c['state_sha256'],base_sha256=c['base_sha256'],coupled=c.get('stage')=='B')
    result,receipt=run_probe(exe,xml,p,a.output/'native')
    result=load(Path(receipt['output']));extra=result['extra_sampling'];audit=result['audit'];K0=load(a.preflight)['heights'][-1]['K0']
    # Native complete acceptance and independent initializer constraints remain
    # separate. Audit field names are retained verbatim in the receipt.
    gates=dict(native_admissible=bool(audit['admissible'][0]) and bool(audit['accepted_awf_hard_admissible'][0]),native_budget=audit['accepted_awf_scalar_domain'][8]<=1e-10,native_algebraic_boundary=max(audit['coupled_norms'][4:])<1e-9,extra_covariance=extra['minimum_covariance']>=-1e-10*K0,strong_divergence=extra['strong_divergence_L2']<=c['divergence_limit']+1e-8,velocity_L2=extra['velocity_change_L2']<=.12*(1+1e-10),sampled_velocity_increment=extra['maximum_sampled_velocity_increment']<=.5+1e-10,optimizer_interior_margin=c['min_lambda_over_K0']>=1e-4)
    out=dict(version='H2eNativeSeedReceiptV1',gates=gates,passed=all(v for k,v in gates.items() if k!='optimizer_interior_margin'),state_sha256=c['state_sha256'],base_sha256=c['base_sha256'],candidate_sha256=sha(a.candidate),runtime=d['runtime'],probe_output_file=receipt['output'],probe_output_sha256=sha(Path(receipt['output'])),receipt=receipt,extra_sampling=extra,audit=audit,contract=result['outer_initializer_contract'],geometry_identity=result['geometry_identity'],inventory_identity=result['inventory_identity'],patch_signature=result['patch_signature'],reference_sha256=result['reference_sha256'],defect_sha256=result['defect_sha256'],physical_steps=0,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    atomic(a.output/'certification.json',out)

def certificate_v2(native,receipt,base,runtime,probe,provenance):
    """Issue only after complete checks; never promote a rejected hill seed."""
    require(receipt['passed'],'cannot certify rejected native/extra initializer')
    require(native['audit']['accepted_awf_hard_admissible']==[1] and native['audit']['admissible']==[1],'cannot certify failed native covariance/scalar inventory')
    require(native['audit']['accepted_awf_scalar_domain'][8]<=1e-10 and max(native['audit']['coupled_norms'][4:])<1e-9,'cannot certify failed budget/algebraic/boundary constraints')
    require(native['extra_sampling']['minimum_covariance']>=-1e-10*provenance['K0'],'cannot certify failed independent covariance')
    require(native['extra_sampling']['strong_divergence_L2']<=provenance['divergence_limit']+1e-8,'cannot certify excessive strong divergence')
    require(native['extra_sampling']['velocity_change_L2']<=.12*(1+1e-10) and native['extra_sampling']['maximum_sampled_velocity_increment']<=.5+1e-10,'cannot certify excessive velocity response')
    require(state_hash(base)==native['reference_sha256'],'certificate reference is not the immutable base')
    contract=native['outer_initializer_contract'];rows={k:contract[k] for k in ('velocity_rows','outer_k_rows','outer_epsilon_rows')}
    allowed=set(sum(rows.values(),[]));fixed=[v for i,v in enumerate(base) if i not in allowed]
    return dict(version='CertifiedNeutralStartupV2',physical_t0=0,complete_native_admissible=True,state=native['state'],state_sha256=native['state_sha256'],immutable_base=base,immutable_base_sha256=state_hash(base),fixed_reference=native['fixed_reference'],fixed_reference_sha256=native['reference_sha256'],fixed_defect=native['fixed_defect'],fixed_defect_sha256=native['defect_sha256'],patch_signature=native['patch_signature'],geometry_identity=native['geometry_identity'],inventory_identity=native['inventory_identity'],allowed_update_contract='OuterScalarInitializerV1',allowed_rows=rows,frozen_base_rows_sha256=state_hash(fixed),scalar_envelope_version='OuterRatioV1:0.5:2:zero-inlet-interface',basis_sha256=provenance['basis_sha256'],context_sha256=provenance['context_sha256'],prescribed_boundary_sha256=state_hash(base),runtime=dict(sha256=runtime['sha256'],library_path=str(Path(runtime['path']).parent/'lib'),libraries=runtime['libraries']),validation_probe=probe,native_receipt_file=provenance['native_receipt_file'],native_receipt_sha256=sha(Path(provenance['native_receipt_file'])),source_revision=provenance['source_revision'],source_hashes=provenance['source_hashes'])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=['context','certify'])
    for n in ('probe','h1','preflight','old-context','context','candidate'):p.add_argument('--'+n,type=Path)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();require(platform.system()=='Linux','authoritative Linux only');os.environ.update(THREADS)
    (context if a.phase=='context' else certify)(a)
if __name__=='__main__':main()
