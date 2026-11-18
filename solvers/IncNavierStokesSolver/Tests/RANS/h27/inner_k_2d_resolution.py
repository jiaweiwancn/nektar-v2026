#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Run the measured directional, quadrature and algebraic F3 controls."""

import argparse
import json
from pathlib import Path

from inner_k_2d_runs import APP, BOUNDARY, INPUT, OUT, PY, admission, audit, operation, read
from resources import identity


def prepare(args):
    from coupled_fixtures import generate, mesh
    from inner_k_fixture_sessions import axes, edit_session
    from inner_k_fixture_oracle import atomic_out

    base = OUT/'base.xml'
    x, z = axes(base)
    spec = read(args.recipe)
    forecast = json.loads(args.forecast)
    case = spec['case']
    x2 = sorted(set(x)|set(spec['extra_x']))
    z2 = sorted(set(z)|set(spec.get('extra_z',[])))
    assert len(x2)>len(x) and x2[0]==x[0] and x2[-1]==x[-1]
    assert z2[0]==z[0] and z2[-1]==z[-1]
    widths = [b-a for a,b in zip(x2,x2[1:])]
    grading = max(max(a/b,b/a) for a,b in zip(widths,widths[1:]))
    assert grading<=2.01 and (len(x2)-1)*(len(z2)-1)<=5808
    caps = spec['caps']
    old = OUT/'base/steady.json'
    topology = OUT/f'{case}-template.xml'
    path = OUT/f'{case}.xml'
    mesh(topology,x2,z2,quadrature=12)
    source = read(args.gate)['diagnostic_source']
    generate(path,BOUNDARY,OUT/case,mode='Physical',solve='SteadyReference',
             dt=.05,steps=0,backend='FastFull3DV1',quadrature_x=8,quadrature_r=12,
             steady_seconds=caps['directional-steady']-45,template=topology,
             initial_fields=old,initial_source_identity=source)
    rows = {case:identity(path)}
    for label, overrides in (
            ('quadrature',dict(WJQuadratureX=12,WJQuadratureR=16)),
            ('algebraic',dict(WJLinearTolerance=3e-13,WJCoupledTolerance=3e-12))):
        tag = f'{case}-{label}'
        rows[tag] = edit_session(path,OUT/f'{tag}.xml',
            dict(WJOutputPrefix=OUT/tag,WJInitialFields=OUT/case/'steady.json'),
            {**overrides,'WJSteadySeconds':caps[label+'-steady']-25})
    atomic_out(OUT/f'{case}-inputs.json',dict(schema='WJInnerK2DResolutionInputsV1',
        case=case, recipe=identity(args.recipe), base_current=identity(OUT/'base-current.json'),
        base_seed=identity(old), same_physical_loads_as=identity(base), sessions=rows,
        x=x2, physical_z=z2, cells=[len(x2)-1,len(z2)-1], maximum_adjacent_x_ratio=grading,
        slab_rule='All coarse faces remain in the common-owner inventory.',
        admission=forecast))


def run(args):
    spec = read(args.recipe)
    base = read(OUT/'base-current.json')
    old = read(OUT/'base-inputs.json')
    source19 = read(Path(old['v2_audit']['path']))
    old_rms = source19['core']['physical_envelope']['domain/inner/k']['RMS']
    assert base['current_pass'] and base['physical_domain']['domain/inner/k']['RMS']+.001 < old_rms
    assert spec['base_current_sha256']==identity(OUT/'base-current.json')['sha256']
    case, caps = spec['case'], spec['caps']
    assert case.isidentifier()
    forecast = admission(caps,spec['memory_forecast_bytes'],args.gate)
    operation(case+'-inputs','analysis',caps['inputs'],*PY,__file__,'prepare',args.gate,args.recipe,json.dumps(forecast))
    base_native = OUT/'base/steady.json'
    native = OUT/case/'steady.json'
    operation(case+'-steady','science',caps['directional-steady'],APP,'steady-fixture',OUT/f'{case}.xml',OUT/case)
    audit(case+'-current',caps['directional-current'],native,OUT/f'{case}-current.json',args.gate,32,'--localize')
    audit(case+'-common-owners',caps['directional-compare'],base_native,OUT/f'base-to-{case}.json',args.gate,20,
          '--other',native,'--control','directional')
    for label in ('quadrature','algebraic'):
        tag = f'{case}-{label}'
        result = OUT/tag/'steady.json'
        operation(tag+'-steady','science',caps[label+'-steady'],APP,'steady-fixture',OUT/f'{tag}.xml',OUT/tag)
        audit(tag+'-current',caps[label+'-current'],result,OUT/f'{tag}-current.json',args.gate,20)
        audit(tag+'-common-owners',caps[label+'-compare'],native,OUT/f'{case}-to-{tag}.json',args.gate,20,
              '--other',result,'--control',label)
    print('INNER_K_F3_RESOLUTION_CONTROLS_FINISHED',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('run','prepare'))
    parser.add_argument('gate',type=Path)
    parser.add_argument('recipe',type=Path)
    parser.add_argument('forecast',nargs='?')
    args = parser.parse_args()
    (run if args.stage=='run' else prepare)(args)


if __name__=='__main__':
    main()
