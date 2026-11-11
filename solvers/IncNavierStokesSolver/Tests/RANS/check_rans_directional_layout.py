#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native tensor rows, faces, masses, continuous trace and finite mixed rank."""
import argparse,os,sys,tempfile,subprocess,xml.etree.ElementTree as ET
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.linalg import block_diag,cholesky,solve_triangular,svdvals
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import write
from qualify_rans_preservation import Q,load,replay_request
from diagnose_rans_balanced_spatial import run_probe
from rans_qualification import write_json
from rans_request_manifest import RequestManifest,runtime_identity
from rans_duration_contract import require


def classify(packet):
    o=packet['observation'];v=packet['tensor_validation'];mixed=packet['mixed_space'];checks=[]
    require(v['schema'] in ('NativeDirectionalTensorChecksV1','NativeDirectionalTensorChecksV2') and len(v['fields'])==6,'tensor_check_coverage')
    for row in v['fields']:
        b=int(row['block']);degrees=row['axis_degrees'];n=len(o['field_basis_configuration'][b])
        require(row['mode_checks']==25*n*(degrees[0]+1)*(degrees[1]+1),'individual_mode_coverage')
        limits=np.array([3e-12,3e-12,3e-12,2e-10,2e-10,3e-12,5e-12,5e-12])
        d=np.array(row['defects_value_dx_dz_dxz_dzz_mass_polynomial_faces'])
        require(d.shape==limits.shape and np.isfinite(d).all() and (d<=limits).all(),'native_tensor_row_or_face_mismatch:'+str(b))
        m=np.array(row['mass']);require(np.max(abs(m-m.T))<=3e-12,'native_mass_not_symmetric');cholesky(m,lower=True)
        checks.append(dict(block=b,defects=d.tolist(),limits=limits.tolist(),mass_positive=True,mode_checks=row['mode_checks']))
        if v['schema']=='NativeDirectionalTensorChecksV2':
            require(row['outer_native_restriction_maximum_defect']<=5e-12 and
                len(row['outer_native_restriction_exports'])==(3 if b<3 else 0),'outer_native_restriction_or_export_failed')
            checks[-1]['outer_native_restriction_maximum_defect']=row['outer_native_restriction_maximum_defect']
    for ground in v['continuous_ground_fixtures']:
        g=ground['diagnostic'];require(g['admissible']==ground['expected_admissible'],'native_continuous_ground_trace:'+ground['name'])
        require(g['admissible']==(g['minimum']>=g['arithmetic_scale']),'ground_arithmetic_classification')
        if ground['name']=='narrow_negative_sextic':require(abs(g['point'][0]-.123456789)<2e-6 and g['minimum']<0 and ground['ordinary_81_point_minimum']>0,'narrow_ground_minimum_missing')
        if ground['name']=='interior_arithmetic_margin':require(abs(g['minimum'])<=g['arithmetic_scale'] and g['arithmetic_scale']>0,'ground_arithmetic_margin_reason_missing')
    if v['fields'][0]['axis_degrees'][0]==6:require(len(v['continuous_ground_fixtures'])==10,'degree_six_ground_coverage')
    a=np.array(mixed['matrix']);rows=np.array(mixed['rows'],int);nv=int(mixed['velocity_size']);np_=len(rows)-nv-1
    offset=np.array(o['offsets'],int);sizes=np.array(o['sizes'],int)
    masses=[np.array(r['mass']) for r in v['fields']]
    full=block_diag(*masses[:2],*masses[4:6],np.ones((1,1)))
    ids=[]
    for row in rows:
        if row<offset[2]:ids.append(row)
        elif row<offset[5]:ids.append(sizes[0]+sizes[1]+row-offset[4])
        elif row<offset[6]:ids.append(sizes[0]+sizes[1]+sizes[4]+row-offset[5])
        else:ids.append(len(full)-1)
    gram=full[np.ix_(ids,ids)];factor=cholesky(gram,lower=True)
    white=solve_triangular(factor,a,lower=True);white=solve_triangular(factor,white.T,lower=True).T
    spectra={name:svdvals(m) for name,m in [('mixed_with_gauge',white),('divergence',white[nv:nv+np_,:nv]),('pressure_action',white[:nv,nv:nv+np_])]}
    ranks={name:int(np.count_nonzero(s>1e-11*s[0])) for name,s in spectra.items()}
    require(ranks['mixed_with_gauge']==len(rows),'additional_mixed_null_mode')
    require(ranks['divergence']==np_-1 and ranks['pressure_action']==np_-1,'additional_pressure_or_divergence_null_mode')
    constant=np.array(mixed['pressure_constant_with_zero_gauge']);raw=a@constant
    constant_defect=float(np.max(abs(raw[:nv]),initial=0.))
    require(constant_defect<=2e-12,'global_pressure_constant_has_velocity_action')
    # The implemented one-sweep pressure inverse uses the free consistent
    # velocity mass, not the momentum block. Inspect that actual Schur system.
    mv=gram[:nv,:nv];g=a[:nv,nv:];d=a[nv:,:nv];c=a[nv:,nv:]
    schur=d@np.linalg.solve(mv,g)-c
    pf=cholesky(gram[nv:,nv:],lower=True)
    sw=solve_triangular(pf,schur,lower=True)
    sw=solve_triangular(pf,sw.T,lower=True).T
    pressure_spectra={name:svdvals(m) for name,m in
        [('pimple_with_gauge',sw),('pimple_without_gauge',sw[:np_,:np_])]}
    pressure_ranks={name:int(np.count_nonzero(s>1e-11*s[0])) for name,s in pressure_spectra.items()}
    require(pressure_ranks['pimple_with_gauge']==np_+1 and
            pressure_ranks['pimple_without_gauge']==np_-1,'additional_pimple_pressure_null_mode')
    return dict(status='pass',tensor_checks=checks,continuous_ground=v['continuous_ground_fixtures'],
        finite_mixed_rank=dict(ranks=ranks,expected=dict(mixed_with_gauge=len(rows),divergence=np_-1,pressure_action=np_-1),
            singular_values={k:s.tolist() for k,s in spectra.items()},relative_rank_threshold=1e-11,
            pimple_pressure_ranks=pressure_ranks,
            pimple_pressure_singular_values={k:s.tolist() for k,s in pressure_spectra.items()},
            common_pressure_constant_action_defect=constant_defect,
            scope='small actual assembled two-region fixture and common gauge; no general inf-sup claim'),
        radial_projection=[r['radial_projection_samples'] for r in v['fields']],unknowns=sum(o['sizes']))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--output',type=Path);a=p.parse_args()
    if a.output is None:a.output=Path(tempfile.mkdtemp(prefix='rans-directional-native-',dir=os.environ.get('TMPDIR')))/'verification'
    if (a.runtime/'lib').exists():os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    layout={'P3':([3,3],[1,1]),'P4':([4,4],[2,2]),'P5':([5,5],[3,3]),'X6Z5':([6,5],[4,3])}
    manifest=RequestManifest(a.output,dict(schema='DirectionalNativeQualificationV2',layouts=layout,knots=[0.,2.,4.],
        fixture='read-only projected scalar/pressure seed, velocity multiplied by .01; its own fixed boundary and reference, no evolution or neutral-equilibrium claim',
        runtime={n:runtime_identity(a.runtime/n) for n in ('RANSPreservationProbe','RANSNewModelsProbe')},
        tensor_limits=[3e-12,3e-12,3e-12,2e-10,2e-10,3e-12,5e-12,5e-12],rank_relative_limit=1e-11))
    rows={};files=[]
    try:
        for name,(degree,pressure) in layout.items():
            folder=a.output/name;case_xml(folder,'A95','Natural','Accepted',steps=0,degree=degree[1],quadrature=Q[1],
                descriptor_probe=a.runtime/'RANSNewModelsProbe',radial_knots=[0.,2.,4.],boundary_policy='FrozenReferenceV1')
            root=ET.parse(folder/'case.xml').getroot()
            for e in root.findall('./EXPANSIONS/E'):e.set('NUMMODES',','.join(str(x+1) for x in (pressure if e.get('FIELDS')=='p' else degree)))
            write(root,folder/'case.xml')
            # A two-element logarithmic interpolant is not an admissible neutral
            # pilot. Use its coefficients only as a construction seed. Retain its
            # physical outcome; the actual read-only fixture must pass all gates.
            write_json(folder/'seed-request.json',dict(operation='neutral',grid=[],export_basis=False))
            with (folder/'seed.log').open('x') as stream:
                code=subprocess.run([str(a.runtime/'RANSPreservationProbe'),'case.xml','seed-request.json','seed-packet.json'],cwd=folder,stdout=stream,stderr=subprocess.STDOUT).returncode
            seed=load(folder/'seed-packet.json');require(seed.get('execution')=='completed','tensor_fixture_seed_construction_failed')
            write_json(folder/'seed-outcome.json',dict(exit_code=code,status=seed['status'],scope='unqualified coefficient construction only; no history or evolution'))
            state=np.array(seed['reference']);end=int(seed['observation']['offsets'][2]);state[:end]*=.01
            request=replay_request(seed,[],state=state.tolist());request['reference']=state.tolist()
            request['fixed_boundary']=[[row,value*(.01 if row<end else 1.)] for row,value in seed['fixed_boundary']]
            request.update(export_tensor_validation=True,export_mixed_space=True)
            packet,_=run_probe(a.runtime,folder,request)
            for b,row in enumerate(packet['tensor_validation']['fields']):
                require(row['axis_degrees']==[* (degree if b<4 else pressure),0],'requested_axis_orders_not_constructed')
            files += [folder/'case.xml',folder/'request.json',folder/'packet.json.gz',folder/'probe.log'];rows[name]=classify(packet)
        radial=max(abs(np.array(x)-np.array(y)).max() for x,y in zip(rows['P5']['radial_projection'],rows['X6Z5']['radial_projection'],strict=True))
        require(radial<=3e-12,'radial_projection_depends_on_streamwise_order')
        result=dict(status='pass',layouts=rows,radial_projection_P5_X6Z5_max_difference=float(radial))
    except Exception as e:result=dict(status='fail',failure=repr(e),layouts=rows)
    write_json(a.output/'results.json',result)
    files=[p for p in a.output.rglob('*') if p.is_file() and p.name not in ('completed-artifacts.json','request-manifest.json')]
    manifest.seal(a.output,files);print(result['status'],str(a.output.resolve()));return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
