#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable compact CT2 export: complete classifications, tables and raw hashes."""
import argparse,csv,hashlib,json,shutil,subprocess
from pathlib import Path
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();out=a.output.resolve();out.mkdir(exist_ok=False);copied={}
    sources={'rough-pressure-control.json':'rough-v1/results.json','primary-initial.json':'initial-v2/results.json','companion-initial.json':'companion-initial-v1/results.json','primary-trajectories.json':'recovered-primary-v1/results.json','companion-trajectories.json':'companion-v1/results.json','primary-audits.json':'primary-audits-v1/results.json','companion-audits.json':'companion-audits-v1/results.json','boundary-verification.json':'boundary-verification-v1/results.json','history-controls.json':'history-controls-v1/results.json','shared-model-smokes.json':'shared-smokes-v1/results.json','continuation02.json':'continue02-v1/results.json','continuation02-audits.json':'continue02-audits-v1/results.json','continuation2.json':'continue2-v1/results.json','continuation2-audits.json':'continue2-audits-v1/results.json','extras.json':'extras-v1/results.json','decision.json':'final-v1/results.json','runtime-build-identity.json':'runtime-identity.json','relocation-runtime-bridge.json':'relocation-runtime-bridge.json','build-execution.json':'build-execution.json','post-reboot-execution.json':'post-reboot-execution.json','followup-execution.json':'followup-execution.json','dataset-audit.json':'dataset-audit.json'}
    required={'rough-pressure-control.json','primary-initial.json','companion-initial.json','primary-trajectories.json','companion-trajectories.json','primary-audits.json','companion-audits.json','boundary-verification.json','history-controls.json','shared-model-smokes.json','extras.json','decision.json','runtime-build-identity.json','relocation-runtime-bridge.json','dataset-audit.json'}
    for name,source in sources.items():
        f=r/source
        if not f.exists():
            if name in required:raise ValueError('missing_required_compact_body:'+name)
            continue
        data=json.loads(f.read_text())
        # Coordinates/weights reside in heavy input packets; pin the one common
        # grid by hash and counts rather than duplicate it in both init receipts.
        if name.endswith('initial.json'):
            data['fixed_grid_sha256']=hashlib.sha256(json.dumps(data.pop('fixed_grid'),separators=(',',':')).encode()).hexdigest()
            data['global_weights_sha256']=hashlib.sha256(json.dumps(data.pop('global_weights'),separators=(',',':')).encode()).hexdigest()
        if name=='decision.json':
            # Full bodies are separately exported, with exact raw receipt hashes.
            data['selected_continuations']={k:dict(status=v.get('physical_completion'),receipt_sha256=sha(r/k/'results.json'),complete_body='continuation02.json' if k=='continue02-v1' else 'continuation2.json') for k,v in data['selected_continuations'].items()}
            data['resource_compliance'].pop('highwater_samples',None)
        dump(out/name,data);copied[name]=dict(source=str(f),raw_sha256=sha(f))
    # Preserve immutable requests/seals, including the original failed collector.
    manifests={}
    for f in sorted(r.glob('*/request-manifest.json')):
        manifests[str(f.relative_to(r))]=dict(sha256=sha(f),body=json.loads(f.read_text()))
    for f in sorted(r.glob('*/completed-artifacts.json')):
        manifests[str(f.relative_to(r))]=dict(sha256=sha(f),body=json.loads(f.read_text()))
    dump(out/'requests-and-seals.json',manifests)
    wall=json.loads((r/'wall-heights-v2/results.json').read_text())
    summary={k:v for k,v in wall.items() if k not in ('flat','curved')};summary['full_results_sha256']=sha(r/'wall-heights-v2/results.json');dump(out/'wall-height-summary.json',summary)
    def table(name,rows,keys):
        with (out/name).open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader()
            for row in rows:
                result={k:row.get(k) for k in keys}
                for k,v in result.items():
                    if isinstance(v,(list,dict)):result[k]=json.dumps(v,separators=(',',':'))
                writer.writerow(result)
    flat=[]
    for row in wall['flat']:
        q=dict(row);q['minimum_k_over_K']=row.get('reference_k',{}).get('minimum_over_K');q['minimum_k_location']=row.get('reference_k',{}).get('location');q['stationary_Q_roots']=row.get('reference_k',{}).get('stationary_Q_roots');flat.append(q)
    table('flat-height-feasibility.csv',flat,['scenario','model','height_mm','hr_mm','hr_lt_H','K','nu','Lambda','cutoff','cutoff_over_H','cutoff_existence_margin','J0','J2','finite_nu_production','matching_velocity','total_reference_traction','unshifted_compatible_z0','minimum_k_over_K','minimum_k_location','stationary_Q_roots','flat_covariance_margin_over_K','flat_covariance_minimum_location','budget_defect','mechanical_work_defect','classification','quadrature_uncertainty','coefficients'])
    curved=[]
    for row in wall['curved']:
        q=dict(row);q['minimum_k_over_K']=row.get('reference_k',{}).get('minimum_over_K');q['minimum_k_location']=row.get('reference_k',{}).get('location');q['stationary_Q_roots']=row.get('reference_k',{}).get('stationary_Q_roots');curved.append(q)
    table('curved-height-feasibility.csv',curved,['scenario','model','height_mm','x_over_Hh','hr_mm','hr_lt_H','curvature','minimum_stretch','global_reach','global_reach_valid','thin_ratio','thin_0p1','thin_0p25','Lambda','cutoff','cutoff_existence_margin','finite_nu_production','matching_velocity','total_reference_traction','minimum_k_over_K','minimum_k_location','stationary_Q_roots','Q_end','budget_defect','classification','quadrature_uncertainty','curved_covariance'])
    # Full high-water ownership attribution and separate outer command timings.
    highwater=[];execution={}
    for d in sorted(r.glob('launch-*')):
        if (d/'execution.json').exists():execution[d.name]=json.loads((d/'execution.json').read_text())
        if (d/'worker-guard-highwater.jsonl').exists():
            highwater.extend(dict(control=d.name,**json.loads(line)) for line in (d/'worker-guard-highwater.jsonl').read_text().splitlines())
    dump(out/'resource-highwater.json',highwater);dump(out/'owned-executions.json',execution)
    hashes={}
    campaigns=['initial-v1','initial-v2','matrix-v1','companion-initial-v1','companion-v1','continue02-v1','continue2-v1','boundary-verification-v1','history-controls-v1','shared-smokes-v1','primary-audits-v1','companion-audits-v1','continue02-audits-v1','continue2-audits-v1','extras-v1']
    for name in campaigns:
        for f in sorted((r/name).rglob('*')):
            if f.is_file() and (f.suffix in ('.xml','.json','.bin','.fld') or any('checkpoint' in part for part in f.parts)):hashes[str(f.relative_to(r))]=sha(f)
    dump(out/'raw-artifact-hashes.json',hashes)
    source=r/'source-after-reboot';cpp=source/'solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp'
    dump(out/'source-scope.json',dict(starting_revision='eedb3738017f46f9c71b857bf814213ad070d7f7',build_revision=json.loads((r/'runtime-identity.json').read_text())['source_revision'],CPP_diagnostic_commit=subprocess.check_output(['git','rev-parse','3c6b826'],cwd=source,text=True).strip(),export_source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=source,text=True).strip(),source_status=subprocess.check_output(['git','status','--porcelain'],cwd=source,text=True),Pimple_cpp_sha256=sha(cpp),CPP_changes='computed-stage diagnostics only; no constitutive, momentum/scalar/outlet/pressure operator equation change',diagnostic_budget_change='physical geometry replaces hardcoded .25 normalization; baseline geometry unchanged',historical_results='unchanged',reused='CT1 functions, original kernel/analytical evidence via byte identities, pinned runtime across explicit mount bridge; primary trajectories reused read-only after collector failure'))
    dump(out/'export-provenance.json',copied)
    dump(out/'hashes.json',{f.name:sha(f) for f in out.iterdir() if f.is_file()})
    print('exported',len(list(out.iterdir())),'files')
    return 0
if __name__=='__main__':raise SystemExit(main())
