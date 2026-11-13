#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H2f evidence, source/runtime identities and indexed Linux artifacts."""
import argparse,datetime,hashlib,json,pathlib,platform,subprocess,xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import sha,require

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--preflight',type=pathlib.Path,required=True);a=p.parse_args();r=a.root;s=r/'source';out=r/'evidence';out.mkdir(exist_ok=False)
    compact=load(r/'compact.json');assert compact['physical_hill_steps']==0
    startup=load(r/'changed-startup/verification.json');checks=load(r/'focused-checks.json');follow=load(r/'startup-followup/checks.json');require(startup['status']=='passed_nonzero_affine_fixture' and checks['passed'] and follow['passed'],'final focused qualification incomplete')
    records={'native.json':compact,'verification.json':dict(focused=checks,startup=follow,helper_replay=load(r/'helper-replay.json'),mass_metric_QP_passed='passed' in (r/'qp-checks.log').read_text(),branch_tests_passed='OK' in (r/'branch-checks.log').read_text()),'startup.json':startup,'stage-A.json':load(r/'stage-A-final/progress.json'),'stage-B.json':load(r/'stage-B/progress.json'),'hybrid.json':load(r/'certify-hybrid/certification.json'),'sparse-mass.json':load(r/'sparse-mass-final/sparse-mass.json'),'symbolic-mass.json':load(r/'symbolic-mass.json'),'zero-replay.json':load(r/'zero-replay.json')}
    for name,d in records.items():atomic(out/name,d)
    base=pathlib.Path('/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2d-20261113/replay-final-v1/PhysicalHeightInterfaceV2-probe/result.json');d=load(base);atomic(out/'baseline-correction.json',dict(state_sha256=d['state_sha256'],native_receipt=str(base),sha256=sha(base),covariance=d['audit']['accepted_awf_covariance'],PDE_solves=0,read_only=True))
    names=['library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionComposite.cpp']
    scripts=s/'solvers/IncNavierStokesSolver/Tests/RANS'
    names += [str(p.relative_to(s)) for p in sorted(scripts.glob('*h2f*.py')) if not p.name.startswith('._')]+['solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_initializer_context.py','solvers/IncNavierStokesSolver/Tests/RANS/prototype_cao_tamura_sparse_mass.py','solvers/IncNavierStokesSolver/Tests/RANS/symbolic_cao_tamura_mass.py']
    source_hashes={n:sha(s/n) for n in dict.fromkeys(names)}
    dataset=s/'docs/plans/cao-tamura-validation-data';data_hashes={p.name:sha(p) for p in dataset.iterdir() if p.is_file()}
    context=load(r/'enriched-context/context.json');xml=pathlib.Path(context['xml']);tree=ET.parse(xml).getroot();params=[e.text for e in tree.findall('./CONDITIONS/PARAMETERS/P')];info={e.attrib['PROPERTY']:e.attrib['VALUE'] for e in tree.findall('./CONDITIONS/SOLVERINFO/I')}
    atomic(out/'case-definition.json',dict(parent='H2e fixed primary PhysicalHeightInterfaceV2',xml=str(xml),xml_sha256=sha(xml),solver_info=info,parameters=params,coupled_slots=2804,physical_hill_steps=0,initializer_direction_counts=context['context']['direction_counts'],base_sha256=context['context']['base_sha256'],geometry_identity=context['geometry_identity'],inventory_identity=context['inventory_identity']))
    runtime={name:load(r/name/'identity.json')['frozen'] for name in ('pilot-solver','pilot-probe')}
    timestamp=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()
    atomic(out/'source-revisions.json',dict(version='H2fSourceRuntimeV1',capture_Asia_Shanghai=timestamp,parent='e93c0fe8f52fee76d23d1ee23c30e6359d47f7e7',source_hashes=source_hashes,policy_sha256=sha(s/'docs/plans/cls-cao-tamura-budget-h2f-evidence/policy.json'),preflight_sha256=sha(a.preflight),frozen_dataset_hashes=data_hashes,platform=platform.platform(),mount=subprocess.check_output(['findmnt','-n','-o','SOURCE,UUID,TARGET','--target',str(r)],text=True).strip(),compiler_flags='Release -fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS',build_workers=4,native_ranks=1,library_threads=1,aggregate_worker_limit=72,runtime=runtime,context_runtime=context['runtime'],changed_fixture_runtime=startup['runtime'],changed_fixture_probe=startup['probe'],stage_A_probe=load(r/'certify-A/certification.json')['runtime'],stage_B_probe=load(r/'certify-B/certification.json')['runtime'],note='Frozen runtime variants and source overlays remain indexed. Probe-only GL20 regional output fields were added after the first context export. Strict schema checks do not alter valid numerical records. The A stopping-label correction preserves its original receipt and exact coefficients.'))
    files=[]
    for p in sorted(r.rglob('*')):
        if not p.is_file() or p.is_symlink():continue
        rel=p.relative_to(r).as_posix()
        if rel.startswith(('source/','build/')):continue
        if p.name=='State.bin' or p.suffix in ('.json','.jsonl','.gz','.log','.bin') or p.name in ('usage.txt',):files.append(dict(path=rel,bytes=p.stat().st_size,sha256=sha(p)))
    for folder in ('native-data','native-extra','native-combined','zero-context','enriched-context'):
        d=load(r/folder/('native-data.json' if folder.startswith('native') else 'context.json'))
        for n,v in d['context']['files'].items():p=pathlib.Path(d['array_directory'])/v['file'];files.append(dict(path=p.relative_to(r).as_posix(),bytes=p.stat().st_size,sha256=v['sha256'],shape=v['shape'],dtype=v['dtype']))
    atomic(out/'artifact-index.json',dict(version='H2fArtifactsV1',root=str(r),capture_Asia_Shanghai=timestamp,files=files,complete_native_hill_certifications=3,physical_hill_attempts=0,physical_hill_steps=0,retention='candidate vectors, input packets, exact native reports, histories, checkpoints, all frozen runtimes and source overlays preserved; other projects untouched'))
if __name__=='__main__':main()
