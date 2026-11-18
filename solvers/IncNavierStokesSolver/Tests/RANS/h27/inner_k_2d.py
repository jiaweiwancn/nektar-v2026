#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent current-state and resolution evidence for the steady inner-k fixture."""

import argparse
from pathlib import Path
import time

import physical_oracle as po
from inner_k_fixture_oracle import FixtureReference, atomic_out, current_subset, state, verify_guards
from inner_k_mechanism import read
from physical_reference_localization import measure
import numpy as np


GUARDS = ('minimum_k','minimum_epsilon','minimum_covariance_over_k','minimum_principal_over_mu')


def inner_k_rows(norms):
    return {key:value for key,value in norms.items() if key.endswith('inner/k')}


def checked_state(path, gate, reference):
    data, info = state(path)
    po.require(gate['F2']=='PASS' and info['diagnostic_source']==gate['diagnostic_source'] and
               info['kernel_source']==gate['kernel_source'] and data['reference_identity']==reference.identity
               and data['inner_k_transport']==reference.inner_k_transport and data['mms_body_amplitude']==0
               and data['boundary_identity']==reference.boundary_identity,
               'native steady candidate is not bound to the tested F2 operators')
    return data, info


def audit(args, reference, gate, deadline):
    data, info = checked_state(args.native,gate,reference)
    oracle = po.PhysicalOracle(data,reference,args.order,deadline)
    core = oracle.assemble()
    boundary = oracle.boundary_checks()
    profiles = oracle.profiles()
    interfaces = oracle.interfaces()
    loc = measure(data,reference,args.order,deadline,sample_extrema=True) if args.localize else None
    minima = {key:data[key] for key in GUARDS}
    guard_error = None
    try:
        verify_guards(minima,core,boundary)
    except ValueError as error:
        guard_error = str(error)
    envelope = core['physical_envelope']
    peak = None if loc is None else max((p['sampled_peak'] for p in loc['columns']),
                                         key=lambda row: abs(row['signed_over_Kin']))
    report = dict(schema='WJInnerK2DStateAssessmentV1', native=po.identity(args.native), state_sha256=data['state_sha256'],
                  diagnostic_source=data['source_identity'], kernel_source=info['kernel_source'],
                  reference_identity=data['reference_identity'], boundary_identity=data['boundary_identity'],
                  native_settings={k:data[k] for k in ('dt','quadratureX','quadratureR','steady_velocity_preconditioner',
                                                       'linear_tolerance','nonlinear_tolerance','initial_field_transfer')},
                  native_statistics={k:info[k] for k in ('construction_seconds','initialization_seconds','steady_seconds',
                                                             'newton_steps','krylov_iterations','jacobian_actions','true_residual')},
                  native_minima=minima, guard_error=guard_error, current=current_subset(core), boundary=boundary,
                  budget=core['budgets'], owner_interfaces=interfaces, quadrature_order=args.order,
                  physical_domain=envelope, physical_stations=profiles['norms'], station_inventory=profiles['stations'],
                  ground_stations=profiles['ground_stations'], localization_peak=peak,
                  inner_k_pass=all(v['status']=='MET' for v in inner_k_rows({**envelope,**profiles['norms']}).values())
                              and (peak is None or abs(peak['signed_over_Kin'])<=po.ENVELOPES['k'][1]),
                  other_physical_envelopes_pass=all(value['status']=='MET' for key,value in
                      {**envelope,**profiles['norms']}.items() if not key.endswith('/k'))
                              and all(abs(row['scaled_deviation'])<=.10 for row in profiles['ground_stations']),
                  current_pass=guard_error is None and core['steady_residual_status']=='MET',
                  peak_station_samples=len(profiles['rows']))
    if args.localize:
        root = args.out.parent
        stem = args.out.stem
        for tag, packet in (('localization',loc),('profiles',profiles),('weak-and-physical-core',core)):
            path = root/f'{stem}-{tag}.json'
            atomic_out(path,packet)
            report[tag] = po.identity(path)
    return report


def compare(args, reference, gate, deadline):
    a = checked_state(args.native,gate,reference)[0]
    if args.control=='baseline':
        b, metadata = po.load_native(args.other)
        other_ref = po.Reference(read(args.table))
        # The independent topology check snaps trusted domain endpoints.
        # Compare the two original inventories before either is normalized.
        po.require(po.identity(args.other)==gate['v2'] and metadata['steady_solver']['converged']
                   and a['elements']==b['elements'] and a['geometry']==b['geometry']
                   and a['dt']==b['dt'] and a['quadratureX']==b['quadratureX'] and a['quadratureR']==b['quadratureR'],
                   'v2 comparison is not mesh- and settings-matched')
        po.require(a['essential_flow_lift']==b['essential_flow_lift'] and a['essential_k_lift']==b['essential_k_lift'],
                   'v2 and fixture inlet or ground lifts changed')
    else:
        b = checked_state(args.other,gate,reference)[0]
        other_ref = reference
    a_oracle = po.PhysicalOracle(a,reference,args.order,deadline)
    b_oracle = po.PhysicalOracle(b,other_ref,args.order,deadline)
    result = a_oracle.compare(b_oracle)
    if args.control=='baseline':
        result['status'] = 'MEASURED'
        result['scope'] = 'Matched mesh and boundary physical difference between two transport laws.'
        for row in result['norms'].values():
            del row['status'], row['limits']
    k = inner_k_rows(result['norms'])
    limits = {'RMS':.001,'MAX':.001} if args.control=='directional' else {'RMS':.0001,'MAX':.0001}
    k_pass = all(v[key]<=limit for v in k.values() for key,limit in limits.items())
    result.update(schema='WJInnerK2DCommonOwnerAssessmentV1', native=po.identity(args.native),
                  other=po.identity(args.other), distinct_reference_laws=args.control=='baseline',
                  physical_control_only=args.control=='baseline', kind=args.control, quadrature_order=args.order,
                  inner_k_numerical_pass=None if args.control=='baseline' else k_pass,
                  inner_k_numerical_limits=None if args.control=='baseline' else limits)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('table','gate','native','out'):
        parser.add_argument('--'+key,type=Path,required=True)
    parser.add_argument('--seconds',type=float,required=True)
    parser.add_argument('--order',type=int,choices=(20,32),required=True)
    parser.add_argument('--localize',action='store_true')
    parser.add_argument('--other',type=Path)
    parser.add_argument('--control',choices=('baseline','quadrature','algebraic','directional'))
    args = parser.parse_args()
    po.require(args.seconds>0 and args.out.parent.is_dir() and not args.out.exists() and
               bool(args.other)==bool(args.control) and not (args.other and args.localize) and
               po.identity(args.table)['sha256']==po.TABLE_SHA, 'invalid bounded analysis request')
    reference, gate = FixtureReference(read(args.table)), read(args.gate)
    deadline = time.monotonic()+args.seconds
    result = compare(args,reference,gate,deadline) if args.other else audit(args,reference,gate,deadline)
    result['oracle'] = po.identity(__file__)
    result['F2_gate'] = po.identity(args.gate)
    atomic_out(args.out,result)
    print(f'INNER_K_2D_ASSESSMENT {args.out.name}',flush=True)


if __name__=='__main__':
    main()
