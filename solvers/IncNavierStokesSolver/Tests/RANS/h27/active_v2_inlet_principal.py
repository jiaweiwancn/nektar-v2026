#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Locate the supported-principal obstruction at the fixed inlet/H trace."""

import argparse
from pathlib import Path
import sys
import time

import physical_oracle as po
from active_v2_inlet import Column, independent, make_profile
from active_v2_inlet_obstruction import fixed_shear, two_ended
from inner_k_columns import planar_response
from inner_k_fixture_oracle import atomic_out
from inner_k_mechanism import read
from reference_oracle import Q_PEAK, shear_law
import numpy as np
from scipy.optimize import brentq


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('table','branch','out'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--seconds',type=float,required=True)
    args = p.parse_args()
    table = po.identity(args.table)
    po.require(sys.platform.startswith('linux') and args.seconds>0 and
               table['sha256']==po.TABLE_SHA,'requires confirmed Linux calibration')
    prior = read(args.branch)
    po.require(prior['status']=='BLOCKED' and prior['table']==table,'branch audit binding differs')
    ref = po.Reference(read(args.table))
    column = Column(ref,time.monotonic()+args.seconds)
    args.out.mkdir()
    k, E, aux = ref.Kin, ref.A/ref.H, float(ref.aux(ref.Kin,ref.H))
    floor = 1e-8

    def witness(q):
        G = q*E/k
        values = [np.array([value]) for value in (k,E,G,aux,1.)]
        report, _ = planar_response(*values,ref.nu)
        shear = ref.nu*G+k*shear_law(q)[0]
        return dict(q=q,shear=shear,**report)

    probes = [witness(float(q)) for q in np.linspace(.001,Q_PEAK,193)]
    key = 'minimum_principal_over_mu'
    crossings = [(lo['q'],hi['q']) for lo,hi in zip(probes[:-1],probes[1:])
                 if (lo[key]-floor)*(hi[key]-floor)<0]
    po.require(len(crossings)==1,'H principal crossing is not uniquely bracketed in the finite scan')
    q_cross = brentq(lambda q:witness(q)[key]-floor,*crossings[0],xtol=2e-11)
    critical = witness(q_cross)
    controls = []
    guess = prior['last_column']['Q_H_over_A']
    for tag,q in (('below',q_cross-2e-5),('above',q_cross+2e-5)):
        trace = witness(q)
        t = trace['shear']/ref.tau
        cases = fixed_shear(column,t,guess,tight=True)
        po.require(len(cases['solutions'])==1,'principal crossing column has no unique sampled bracket')
        item = cases['solutions'][0]
        x = [t,item['Q_H_over_A']]
        state,sample = column.integrate(x,tight=True)
        anchor = item['U_H']
        assessment,fields = independent(column,x,state,sample,anchor,
                                       require_admissible=False,forward_control=False)
        report = dict(label=tag,H=trace,trials=cases,current=assessment,U_H=anchor,
                      gap=prior['target_U_H']-anchor)
        if tag=='below' and assessment['sampled_admissibility_pass']:
            report['independent_match'] = two_ended(column,t,item['k_ground'],x[1],anchor,sample)
            try:
                artifact, report['witness_representation'] = make_profile(column,sample,anchor,fields,18)
            except ValueError as error:
                report['witness_representation'] = dict(status='NOT_MET',reason=str(error))
            else:
                artifact.update(schema='WJH27UnmatchedColumnWitnessV1', U_H=anchor,
                                original_outer_target=prior['target_U_H'],T_in=trace['shear'],Q_H=assessment['Q_H'])
                path = args.out/'unmatched-witness-polynomials.json'
                atomic_out(path,artifact)
                report['diagnostic_polynomial']=po.identity(path)
        controls.append(report)
    po.require(controls[0]['H'][key]>floor and controls[1]['H'][key]<floor and
               all(v['gap']>0 for v in controls),'principal crossing did not corroborate the obstruction')
    atomic_out(args.out/'principal-assessment.json',dict(schema='WJH27MatchedInletPrincipalObstructionV1',
        prior=po.identity(args.branch),table=table,first_sampled_H_crossing=critical,H_bracket=crossings[0],
        controls=controls,scope='Fixed K and E at H. Finite q scan plus directional stationary points of the supported planar quartic.'))
    atomic_out(args.out/'H-principal-scan.json',probes)
    print('ACTIVE_V2_INLET_PRINCIPAL_OBSTRUCTION_CHECKED',flush=True)


if __name__=='__main__':
    main()
