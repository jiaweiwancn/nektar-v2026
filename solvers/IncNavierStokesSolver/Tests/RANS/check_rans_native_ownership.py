#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Production quadrature ownership, native Stokes and retained frozen witnesses.

All rows come from Refresh/ResidualResolution in the real native patch. Expected
volumes, constants and Stokes identities are independent geometry identities.
Counterfactual test-basis replacement is never used by this qualification.
"""
import argparse
import gzip
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
import numpy as np
from rans_request_manifest import RequestManifest, sha, runtime_identity
from rans_preservation_oracle import Reference, weak_balance
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info, set_parameter, write

OPERATOR = 'GeneratedVolumeOwnerV1'
LEGACY = 'LegacyCoordinateToleranceV1'


def load(path):
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as f:
        return json.load(f)


def primitives(packet, ref):
    o = packet['observation']
    assert o['numerical_operator'] == OPERATOR
    elements = {e['geometry']: e for e in o['elements']}
    n = o['sizes'][0]
    global_volume = np.zeros((3, n), dtype=np.longdouble)
    global_face = np.zeros_like(global_volume)
    local = {g: np.zeros_like(global_volume) for g in elements}
    volumes = {g: 0. for g in elements}
    wrong = []; unity = 0.; gradient = 0.; roundoff = 0; tiny = 0
    panels = set()
    for q in o['volume']:
        row = q['velocity_k']; g = row['geometry']; e = elements[g]
        x = q['coordinate']; w = q['weight']; ids = np.asarray(row['ids'], int)
        if q['generating_geometry'] != g or any(
                not e['lower'][a] <= x[a] <= e['upper'][a] for a in (0, 2)):
            wrong.append(dict(point=x, selected=g, generating=q['generating_geometry']))
        assert q['pressure']['geometry'] == g
        if q['region'] == 1: assert q['epsilon']['geometry'] == g
        a, b = q['generating_radial_panel']; assert b > a
        if q['region'] == 0:
            assert e['lower'][2] <= a <= b <= e['upper'][2]
            panels.add((a, b))
            tiny += int(b-a < 1e-11)
        c = np.asarray(row['constant'], np.longdouble)
        phi = np.asarray(row['phi'], np.longdouble)
        derivative = np.asarray(row['derivative'], np.longdouble)
        unity = max(unity, float(abs(phi@c-1)))
        for axis in (0, 2):
            gradient = max(gradient, float(abs(derivative[axis]@c)) *
                           (e['upper'][axis]-e['lower'][axis]))
        roundoff += int(row['endpoint_roundoff_axes'] != 0)
        v = np.longdouble(w)*derivative
        global_volume[:, ids] += v; local[g][:, ids] += v
        volumes[g] += w
    for f in o['faces']:
        ids = np.asarray(f['velocity_k']['ids'], int)
        action = (np.longdouble(f['weight']) * np.asarray(f['area'], np.longdouble)[:, None]
                  * np.asarray(f['velocity_k']['phi'], np.longdouble))
        global_face[:, ids] += action
        if f['area'][0] != 0: local[f['velocity_k']['geometry']][:, ids] -= action
    # All native radial faces, including internal faces, are exported one-sided.
    for f in packet['resolution']['one_sided_faces']:
        ids = np.asarray(f['test']['ids'], int)
        local[f['geometry']][2, ids] -= (np.longdouble(f['weight'])*f['outward_normal_z'] *
                                        np.asarray(f['test']['phi'], np.longdouble))
    global_stokes = float(np.max(np.abs(global_volume-global_face)))
    element_stokes = max(float(np.max(np.abs(a))) for a in local.values())
    volume_error = max(abs(volumes[g]-(e['upper'][0]-e['lower'][0]) *
                           (e['upper'][2]-e['lower'][2])) for g,e in elements.items())
    fixed = {int(i) for i,v in packet['resolution']['constraints']}
    free_flux = 0.
    for block,axis in ((0,0),(1,2)):
        ids = [i for i in range(n) if o['offsets'][block]+i not in fixed]
        free_flux = max(free_flux, float(np.max(np.abs(global_volume[axis,ids]), initial=0)))
    pressure_volume = 0.
    for q in o['pressure_mass']:
        pressure_volume += q['weight']*float(np.dot(q['pressure']['phi'], q['pressure']['constant']))
    exact_volume = sum((e['upper'][0]-e['lower'][0])*(e['upper'][2]-e['lower'][2])
                       for e in elements.values())
    weak = weak_balance(ref, o)
    limits = dict(unity=2e-12, scaled_constant_gradient=2e-12,
                  global_stokes=2e-12, element_stokes=2e-12,
                  element_volume=2e-12, free_constant_pressure_flux=2e-12, pressure_gauge_volume=2e-12)
    values = dict(unity=unity, scaled_constant_gradient=gradient,
                  global_stokes=global_stokes, element_stokes=element_stokes,
                  element_volume=volume_error, free_constant_pressure_flux=free_flux,
                  pressure_gauge_volume=abs(pressure_volume-exact_volume))
    ok = not wrong and all(values[k] <= limits[k] for k in limits) and weak['status']=='pass'
    return dict(status='pass' if ok else 'fail', ownership_errors=wrong,
                primitive_defects=values, primitive_limits=limits,
                observed_local_endpoint_roundoff_rows=roundoff,
                tiny_panel_points=tiny, minimum_panel_width=min(b-a for a,b in panels),
                analytical_native_weak_balance=weak,
                numerical_operator=OPERATOR, performed_time_steps=packet['performed_time_steps'],
                scope='actual production rows; analytical weak balance is distinct from represented-state residual')


def permute_geometry(path):
    root = ET.parse(path).getroot(); quads = root.findall('./GEOMETRY/ELEMENT/Q')
    ids = sorted(int(q.get('ID')) for q in quads); mapping = dict(zip(ids, reversed(ids)))
    for q in quads: q.set('ID', str(mapping[int(q.get('ID'))]))
    def replacement(match):
        values = []
        for word in match.group(1).split(','):
            pair = word.split('-'); values += list(range(int(pair[0]), int(pair[-1])+1))
        return 'Q['+','.join(str(mapping[i]) for i in sorted(values))+']'
    for c in root.findall('./GEOMETRY/COMPOSITE/C'):
        c.text = re.sub(r'Q\[([^]]+)\]', replacement, c.text)
    write(root, path)
    return mapping


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--descriptor-probe', type=Path, required=True)
    p.add_argument('--old-probe', type=Path)
    p.add_argument('--previous-root', type=Path)
    p.add_argument('--output', type=Path)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--directional', action='store_true', help='qualified uniform (6,5)/(4,3) route only')
    a = p.parse_args()
    out = a.output.resolve() if a.output else Path(tempfile.mkdtemp(prefix='rans-owned-volume-',dir=os.environ.get('TMPDIR')))/'evidence'
    a.probe=a.probe.resolve();a.descriptor_probe=a.descriptor_probe.resolve()
    if a.old_probe:a.old_probe=a.old_probe.resolve()
    if a.previous_root:a.previous_root=a.previous_root.resolve()
    if a.directional and (a.old_probe or a.previous_root):raise ValueError('directional layout has no legacy operator witness')
    fixture = Path(__file__).resolve().parents[4]/'docs/plans/a93-a95-e94e95-awf-nonlinear-duration-response-qualification-evidence/minimal-ownership/case.xml'
    frozen={}
    if a.previous_root:
        for wall in ('Accepted','EP'):
            source=a.previous_root/('evolved-audits/response-A95-N-'+wall+'-response-2e-05-dt-0.002-0.016-B5/T-0.0160/packet.json.gz')
            frozen[wall]=dict(packet_sha256=sha(source),configuration_sha256=sha(source.parent/'case.xml'))
    intent = dict(schema='NativeOwnershipGeometryContractV1', numerical_operator=OPERATOR,
                  probe=runtime_identity(a.probe), descriptor=runtime_identity(a.descriptor_probe),
                  old_probe=runtime_identity(a.old_probe) if a.old_probe else None,
                  previous_root=str(a.previous_root), source_sha256=sha(Path(__file__)),
                  fixture_sha256=sha(fixture),frozen_inputs=frozen,
                  primitive_absolute_tolerance=2e-12,
                  weak_tolerances='unchanged native-test oracle: scaled 1e-9 or absolute 1e-12',
                  method='production Refresh rows; native cutoff used as fixed mesh node; no row replacement')
    if a.directional:intent['directional_layout']='UniformDirectionalModifiedA2DV1:6,5:4,3'
    seal = RequestManifest(out, intent, a.resume)
    if seal.completed(out): return 0 if load(out/'results.json')['status']=='pass' else 1
    results = []; sealed_files = []; ref = Reference('A95','Natural','Accepted')

    def run(name, request, source_xml=None, us=.4, node=None, reverse=False, old=False, expected_error=None):
        folder = out/name; folder.mkdir(exist_ok=False)
        if source_xml:
            shutil.copy2(source_xml,folder/'case.xml')
        else:
            folder.rmdir()
            knots = sorted(set([0.,.0025,.005,.01,node or ref.zc,.02,.04,.1,.5,1.,2.,3.,4.]))
            case_xml(folder,'A95','Natural','Accepted',steps=0,degree=3,radial_knots=knots,
                     quadrature=(48,12,12,12,12),descriptor_probe=a.descriptor_probe,
                     boundary_policy='FrozenReferenceV1',cache='ExactStepV1')
        root=ET.parse(folder/'case.xml').getroot()
        set_info(root,'RANSWallRestart','None');set_info(root,'RANSWallCheckpoint',str(folder/'unused-checkpoint'))
        set_parameter(root,'RANSWallFrictionVelocity',us);set_parameter(root,'K0',ref.q*us*us)
        if a.directional:
            for e in root.findall('./EXPANSIONS/E'):
                e.set('NUMMODES','5,4' if e.get('FIELDS')=='p' else '7,6')
        write(root,folder/'case.xml')
        mapping=permute_geometry(folder/'case.xml') if reverse else None
        write_json(folder/'request.json',request)
        executable=a.old_probe if old else a.probe
        started=time.monotonic()
        with (folder/'probe.log').open('x') as f:
            code=subprocess.run([str(executable),'case.xml','request.json','packet.json'],cwd=folder,
                                stdout=f,stderr=subprocess.STDOUT,timeout=14400).returncode
        item=dict(name=name,exit_code=code,seconds=time.monotonic()-started,
                  configuration_sha256=sha(folder/'case.xml'),request_sha256=sha(folder/'request.json'),
                  geometry_permutation=mapping)
        if expected_error:
            item.update(status='pass' if code!=0 and expected_error in (folder/'probe.log').read_text() else 'fail',
                        expected_rejection=expected_error)
            results.append(item);write_json(out/'progress.json',dict(checks=results));return None
        if code or not (folder/'packet.json').exists():
            item.update(status='fail',failure=(folder/'probe.log').read_text()[-2500:]);results.append(item)
            write_json(out/'progress.json',dict(checks=results));return None
        body=load(folder/'packet.json');assert body['performed_time_steps']==0
        item['packet_sha256']=sha(folder/'packet.json')
        if old:
            obs=body['observation'];point=obs['grid'][0]
            owner=[e['geometry'] for e in obs['elements'] if all(e['lower'][d]<point[d]<e['upper'][d] for d in (0,2))]
            selected=obs['grid_geometry_by_field'][0][0]
            item.update(status='pass' if owner==[5] and selected==4 else 'fail',
                        native_qualification='fail',expected_failure='native_geometry_not_containing_point',
                        selected_geometry=selected,strict_containing_geometry=owner)
        elif request.get('export_basis'):
            item.update(primitives(body,Reference('A95','Natural','Accepted',us)))
        else:
            obs=body['observation'];point=obs['grid'][0]
            owner=[e['geometry'] for e in obs['elements'] if all(e['lower'][d]<point[d]<e['upper'][d] for d in (0,2))]
            item.update(status='pass' if obs['grid_geometry_by_field'][0][0] in owner else 'fail',
                        selected_geometry=obs['grid_geometry_by_field'][0][0],strict_containing_geometry=owner)
        if not old:
            assert body['numerical_operator']==OPERATOR
            if request.get('export_basis'):
                actual=body['observation']['profiles'][0]['cutoff']
                item['cutoff']=actual;item['fixed_mesh_node']=node or ref.zc
                item['cutoff_minus_node']=actual-item['fixed_mesh_node']
        results.append(item)
        # Preserve the original serialized packet, without retaining a second uncompressed copy.
        with gzip.open(folder/'packet.json.gz','wb',compresslevel=3) as f:f.write((folder/'packet.json').read_bytes())
        (folder/'packet.json').unlink();sealed_files.append(str((folder/'packet.json.gz').relative_to(out)))
        write_json(out/'progress.json',dict(checks=results))
        print(name,item['status'],flush=True)
        return body

    point=[-.7699026741943047,0.,.014912396183309551]
    request=dict(operation='neutral',grid=[point],grid_sides=[0],export_test_basis=True,grid_test_geometries=[5])
    if a.old_probe: run('old-minimal',request,source_xml=fixture,old=True)
    run('corrected-minimal',request,source_xml=fixture)
    run('wrong-explicit-owner',dict(request,grid_test_geometries=[4]),source_xml=fixture,
        expected_error='outside its recorded native element')
    # Reversed IDs: let the strict locator identify the correct physical owner.
    run('reversed-ids-minimal',dict(request,grid_test_geometries=[]),source_xml=fixture,reverse=True)
    volume_request=dict(operation='neutral',grid=[],export_basis=True,export_resolution=True)
    seed=run('cutoff-seed',volume_request)
    if seed:
        node=seed['observation']['profiles'][0]['cutoff']
        for label,us in [('coincident',.4),('lower-input',.4*(1-1e-10)),('upper-input',.4*(1+1e-10))]:
            run('panel-'+label,volume_request,us=us,node=node)
        run('reversed-ids-volume',volume_request,node=node,reverse=True)
        variants=[r for r in results if r['name'].startswith('panel-')]
        signs=[r.get('cutoff_minus_node') for r in variants]
        results.append(dict(name='cutoff-position-coverage',status='pass' if len(signs)==3 and signs[0]==0 and min(signs)<0<max(signs) and
                            all(r['tiny_panel_points']>0 for r in variants[1:]) else 'fail',
                            displacements=signs,scope='fixed mesh, physical input perturbations only; no snapped coordinates or removed panels'))
        trace=run('exact-and-near-interface',dict(volume_request,
            grid=[[-.37,0.,2.],[-.37,0.,2.],[-.37,0.,2.-1e-13],[-.37,0.,2.+1e-13]],
            grid_sides=[-1,1,0,0]),node=node)
        if trace:
            o=trace['observation'];geometries={e['geometry']:e for e in o['elements']}
            values=np.array(o['grid_values']);jump=float(np.max(abs(values[0,[0,1,2]]-values[1,[0,1,2]])))
            owners=o['grid_geometry_by_field'];regions=[0,1,0,1]
            owner_ok=all((geometries[g[0]]['lower'][2]>=2.)==bool(region) and
                g[0]==g[1]==g[2]==g[4] and (g[3]==g[0] if region else g[3]==-1)
                for g,region in zip(owners,regions,strict=True))
            results.append(dict(name='exact-side-and-nearby-owner-separation',status='pass' if owner_ok and jump<=2e-11 else 'fail',
                common_velocity_k_jump=jump,actual_field_owners=owners,declared_sides=o['grid_sides'],
                one_sided_velocity_and_k_gradients=o['grid_gradients'],scope='no gradient or eddy-viscosity continuity imposed'))
    if a.previous_root:
        for wall in ('Accepted','EP'):
            source=a.previous_root/('evolved-audits/response-A95-N-'+wall+'-response-2e-05-dt-0.002-0.016-B5/T-0.0160/packet.json.gz')
            old=load(source);o=old['observation']
            req=dict(operation='legacy_operator_replay',identity=o['identity'],descriptor=o['resolved_descriptor'],
                     state=o['state'],reference=old['reference'],fixed_boundary=old['fixed_boundary'],time=o['time'],
                     source_numerical_operator=LEGACY,target_numerical_operator=OPERATOR,
                     grid=[],export_basis=True,export_resolution=True)
            # EP has a distinct reference; apply its unchanged oracle explicitly.
            folder=out/('frozen-'+wall);folder.mkdir();shutil.copy2(source.parent/'case.xml',folder/'case.xml')
            write_json(folder/'request.json',req)
            with (folder/'probe.log').open('x') as f:
                code=subprocess.run([str(a.probe),'case.xml','request.json','packet.json'],cwd=folder,
                                    stdout=f,stderr=subprocess.STDOUT,timeout=14400).returncode
            item=dict(name='frozen-'+wall,exit_code=code,source_packet_sha256=sha(source))
            if code==0:
                body=load(folder/'packet.json');item.update(primitives(body,Reference('A95','Natural',wall)))
                old_points=np.asarray([q['coordinate'] for q in o['volume']],dtype='<f8')
                new_points=np.asarray([q['coordinate'] for q in body['observation']['volume']],dtype='<f8')
                old_weights=np.asarray([q['weight'] for q in o['volume']],dtype='<f8')
                new_weights=np.asarray([q['weight'] for q in body['observation']['volume']],dtype='<f8')
                same=old_points.tobytes()==new_points.tobytes() and old_weights.tobytes()==new_weights.tobytes()
                item.update(coordinates_and_weights_bitwise_unchanged=same,source_numerical_operator=body['source_numerical_operator'],
                            target_numerical_operator=body['numerical_operator'],trajectory_requalified=False)
                if not same:item['status']='fail'
                with gzip.open(folder/'packet.json.gz','wb',compresslevel=3) as f:f.write((folder/'packet.json').read_bytes())
                (folder/'packet.json').unlink();sealed_files.append(str((folder/'packet.json.gz').relative_to(out)))
            else:item.update(status='fail',failure=(folder/'probe.log').read_text()[-2500:])
            results.append(item);write_json(out/'progress.json',dict(checks=results));print(item['name'],item['status'],flush=True)
    result=dict(status='pass' if results and all(r['status']=='pass' for r in results) else 'fail',
                contract=intent,checks=results,
                scope='native production geometry and analytical weak identities only; no evolution qualification or reference residual subtraction')
    sealed_files += [str(p.relative_to(out)) for pattern in ('*/case.xml','*/request.json','*/probe.log') for p in out.glob(pattern)]
    write_json(out/'results.json',result);seal.seal(out,['results.json',*sealed_files])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
