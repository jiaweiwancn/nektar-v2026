#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Measured quadrature and evaluation uncertainty; no continuum changes."""
import argparse
from pathlib import Path
import shutil
import sys
from rans_qualification import write_json,compare_samples,exact_vector_hash
from qualify_rans_preservation import load,probe,replay_request,Q
from rans_preservation_oracle import Reference,physical_grid,representation
from run_a93_a95_e94e95_matrix import run_case

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root
    import os
    os.environ['LD_LIBRARY_PATH']=str(r/'runtime/lib');os.environ['TMPDIR']=str(r/'tmp')
    output=r/'uncertainty';output.mkdir(exist_ok=False);grid=load(r/'preservation-spatial/evaluation-grid.json');weights=grid['weights'];result=dict(checks=[],status='running')
    def save():write_json(output/'results.json',result)
    for model in ('E94E95','A95','A93'):
        for wall in ('EP','Accepted'):
            name=model+'-Natural-'+wall+'-dt0'
            try:
                ref=Reference(model,'Natural',wall)
                lowdir=r/('preservation-temporal-Q0-'+model)/name
                lo=load(lowdir/'state-packet.json');lor=load(lowdir/'case.evolution.json')
                if wall=='EP':
                    highdir=r/('preservation-temporal-'+model)/name;hi=load(highdir/'state-packet.json');hir=load(highdir/'case.evolution.json')
                else:
                    highdir=output/name
                    row=run_case(r/'runtime/IncNavierStokesSolver',highdir,model,'Natural',wall,steps=16,dt=.0001,descriptor_probe=r/'runtime/RANSNewModelsProbe',
                                 radial_knots=ref.knots(),degree=3,quadrature=Q[2],cache='ExactStepV1',checkpoint_interval=16)
                    if row['status']!='pass':raise RuntimeError('Q2 control failed')
                    hir=load(highdir/'case.evolution.json')
                    hi=probe(r/'runtime/RANSPreservationProbe',highdir,replay_request(lo,grid['grid'],state=hir['state'],time=hir['time']))
                assert exact_vector_hash(lor['accepted_states'][0])==exact_vector_hash(hir['accepted_states'][0]),'quadrature controls start at different coefficient states'
                assert exact_vector_hash(hi['observation']['state'])==exact_vector_hash(hir['state']),'high-order sample packet is not final native state'
                comparison=compare_samples(dict(hi['observation'],grid_weights=weights),dict(lo['observation'],grid_weights=weights),weights)
                maximum=max(v['weighted_RMS'] for v in comparison['fields'].values())
                result['checks'].append(dict(name='Q2-Q0-'+model+'-'+wall,status='pass' if maximum<1e-10 else 'fail',comparison=comparison,
                    declaration='all field RMS differences below predeclared 1e-10 uncertainty floor',high=str(highdir),low=str(lowdir)))
            except Exception as error:result['checks'].append(dict(name='Q2-Q0-'+model+'-'+wall,status='fail',failure=repr(error)))
            save()
    for wall in ('EP','Accepted'):
        try:
            ref=Reference('E94E95','Natural',wall);tag='E94E95-Natural-'+wall;directory=r/'preservation-spatial';base=load(directory/(tag+'-base.json'));base_metrics=load(directory/(tag+'-h0')/'metrics.json')
            dense,dense_weights=physical_grid(grid['panels'],14,10);folder=output/('evaluation-'+wall);folder.mkdir();shutil.copy2(directory/(tag+'-h0')/'case.xml',folder/'case.xml')
            hi=probe(r/'runtime/RANSPreservationProbe',folder,replay_request(base,dense));higher=representation(ref,hi['observation'],dense_weights)
            changes={name:abs(base_metrics['representation']['fields'][name]['weighted_RMS']-higher['fields'][name]['weighted_RMS']) for name in higher['fields']}
            gradients={name:abs(base_metrics['representation']['gradient_weighted_RMS'][name]-value) for name,value in higher['gradient_weighted_RMS'].items()}
            result['checks'].append(dict(name='evaluation-GL10-GL14-'+wall,status='pass' if max([*changes.values(),*gradients.values()])<1e-10 else 'fail',
                value_RMS_changes=changes,gradient_RMS_changes=gradients,lower_orders=[10,8],higher_orders=[14,10],same_state_sha256=exact_vector_hash(base['observation']['state'])))
        except Exception as error:result['checks'].append(dict(name='evaluation-GL10-GL14-'+wall,status='fail',failure=repr(error)))
        save()
    result['mass_quadrature']=dict(actual_global_flat_rule='Gauss max(8,p+3) per direction, from RANSWallGlobalField::Mass/Quadrature; patch override affects regional/diagnostic integration only',
        degrees=[3,4,5],maximum_mass_integrand_degree=10,minimum_exact_Gauss_degree=15,status='exact_polynomial_rule',
        temporal_comparison_rule='common panel grid contains every element boundary; GL10 integrates squared p<=5 polynomial differences exactly in exact arithmetic',
        projection_rule='fixed independent Gauss32 geometric subpanels; its coefficient approximation belongs to representation, not frozen residual quadrature')
    expected={'Q2-Q0-'+m+'-'+w for m in ('E94E95','A95','A93') for w in ('EP','Accepted')}|{'evaluation-GL10-GL14-'+w for w in ('EP','Accepted')}
    from rans_qualification import coverage
    result['coverage']=coverage(result['checks'],expected,key=lambda x:x['name']);result['status']='pass' if result['coverage']['status']=='pass' and all(x['status']=='pass' for x in result['checks']) else 'fail';save();print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
