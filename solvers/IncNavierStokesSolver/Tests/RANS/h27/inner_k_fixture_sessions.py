#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Create isolated physical sessions with the retained v2 boundary data."""

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

from coupled_fixtures import generate, mesh
from inner_k_fixture_oracle import FixtureReference, atomic_out
from inner_k_mechanism import read
from physical_oracle import TABLE_SHA, identity, require


SOURCE19 = '56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'


def axes(path):
    coordinates = [tuple(map(float, v.text.split())) for v in ET.parse(path).findall('GEOMETRY/VERTEX/V')]
    return tuple(sorted({p[i] for p in coordinates}) for i in (0,1))


def edit_session(template, path, options, parameters):
    doc = ET.parse(template)
    si = doc.find('CONDITIONS/SOLVERINFO')
    for key,value in options.items():
        item = si.find(f"I[@PROPERTY='{key}']")
        if item is None:
            item = ET.SubElement(si,'I',PROPERTY=key)
        item.set('VALUE',str(value))
    params = doc.find('CONDITIONS/PARAMETERS')
    for key,value in parameters.items():
        # Replace every existing occurrence: some retained inputs have
        # the same Krylov parameter twice, with identical values.
        matches = [p for p in params if p.text.split('=',1)[0].strip()==key]
        if not matches:
            matches = [ET.SubElement(params,'P')]
        for item in matches:
            item.text = f'{key} = {value}'
    with path.open('xb') as stream:
        doc.write(stream,encoding='utf-8')
    with Path(str(path)+'.opt').open('xb') as stream:
        stream.write(Path(str(template)+'.opt').read_bytes())
    return identity(path)


def f2(root, folder):
    science = root/'science'
    table = science/'reference-tables-v1.json'
    require(identity(table)['sha256']==TABLE_SHA,'boundary table changed')
    ref = FixtureReference(read(table))
    out = science/folder
    out.mkdir()
    base = science/'cp2-refinement-v19/base.xml'
    old = science/'cp2-refinement-v19/base.result.json'
    parity = edit_session(base,out/'parity.xml',
                          dict(WJInitialFields=old,WJInitialSourceIdentity=SOURCE19,WJOutputPrefix=out/'parity-native'),
                          dict(WJSteadySeconds=65))
    physical = science/'cutoff-panels-v19/physical.xml'
    x, z = axes(physical)
    d0 = z[0]
    require(abs(d0-.0875)<1e-13 and (len(x)-1)*(len(z)-1)==700 and
            all(s in x for s in (-2.5,0.,2.5)),'unexpected retained diagnostic geometry')
    c = float(ref.cutoff(ref.Kin))
    face = d0+c
    nearby = [s for s in z if abs(face-s)<=1e-8]
    require(len(nearby)<=1 and (not nearby or abs(face-nearby[0])<=1e-12),
            'ambiguous archived historical cutoff face')
    radii = sorted((set(z)-set(nearby)) | {face})
    template = out/'actions-template.xml'
    mesh(template,x,radii,quadrature=12)
    case = out/'actions.xml'
    generate(case,table,out/'actions-native',mode='Physical',solve='SteadyReference',
             dt=.05,steps=0,backend='FastFull3DV1',quadrature_x=12,quadrature_r=16,
             steady_seconds=60,template=template)
    return out,dict(parity=parity,parity_baseline=identity(base),parity_seed=identity(old),
                    action=identity(case),action_baseline=identity(physical),axis_size=[len(x),len(radii)],
                    cutoff_face_z=face,retained_matching_face=nearby,table=identity(table))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase',choices=['f2'])
    p.add_argument('root',type=Path)
    p.add_argument('--folder',default='inner-k-f2',choices=['inner-k-f2','inner-k-f2-r1'])
    args = p.parse_args()
    folder, report = f2(args.root,args.folder)
    atomic_out(folder/'sessions.json',report)


if __name__=='__main__':
    main()
