#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""INACTIVE v3 counterfactual scalar-diffusion diagnostic; never runtime data."""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import fast_reference as fr
import physical_oracle as po
import numpy as np


def identity(path):
    path = Path(path)
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('table', 'proposal', 'out'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    args = parser.parse_args()
    start = time.monotonic()
    report = dict(schema='WJH27InactiveAuxDiffusionReferenceControlV1',
        status='FAIL', proposal_only=True, active_model_changed=False,
        runtime_table_produced=False, physical_reference_qualified=False,
        full_H27_qualification=False, TASK3_READY=False,
        scope='Finite counterfactual columns only; no two-dimensional reference, evolving PDE, covariance/principal or success inference',
        proposal=identity(args.proposal), existing_v2_table=identity(args.table),
        sources=dict(helper=identity(__file__), generator=identity(fr.__file__),
                     complete_IVP=identity(fr.oracle.__file__), inlet=identity(po.__file__)),
        factors=[.9, .95, 1., 1.05, 1.1], rows=[])
    error = None
    try:
        assert sys.platform.startswith('linux') and 0 < args.seconds <= 180
        assert report['existing_v2_table']['sha256'] == po.TABLE_SHA
        gen = fr.Generator(args.seconds)
        # The existing independent Column already names this alternate
        # mathematical diffusion law. No active PDE/table/option is changed.
        gen.p = replace(gen.p, inner_diffusion='auxiliary')
        assert gen.p.sigma_k == 1
        report['candidate_parameters'] = asdict(gen.p)
        report['equation'] = 'Dkw=nu+nuA(K,r); unchanged approved v2 parallel shear and raw production'
        for factor in report['factors']:
            fast = gen.fast(factor)
            row = gen.full_validation(factor, fast['logc'])
            row.update(roughness_intercept=fast['column'].zv,
                conditional_velocity_anchor=fast['anchor'],
                conditional_production_over_A=fast['production_over_A'],
                zero_flux_logcutoff_slope_error=abs(row['complete_IVP_logcutoff_slope']-1))
            report['rows'].append(row)
            assert 0 < row['c_over_H'] < .4 and row['roughness_intercept'] < 0
            assert max(abs(row['ground_flux_over_A']), abs(row['independent_budget_over_A']),
                       row['epsilon_corner_relative_error']) <= 1e-10
            assert row['zero_flux_logcutoff_slope_error'] <= 1e-8
        column, trial = gen.full_cache[float(1.).hex()]
        reference = po.Reference(json.loads(args.table.read_text()))
        c3, c2 = trial['c'], float(reference.cutoff(np.array([reference.Kin]))[0])
        radii = sorted(set(np.linspace(0, column.p.H, 129).tolist()+
                           [c2, c3, reference.old_c, .4*column.p.H, .75*column.p.H]))
        r = np.asarray(radii)
        old = reference.inlet(r)['k']
        candidate = np.asarray([column.p.K*trial['sample'](v/column.p.H)[0] for v in r])
        diffusion3 = np.asarray([column.p.nu+column.auxiliary(v/column.p.H) for v in r])
        destruction2 = column.A/np.maximum(r, c2)
        diffusion2 = column.p.nu+column.p.cd*old*old/(column.p.sigma_k*destruction2)
        report['inlet_diagnostic'] = dict(r=r.tolist(), physical_z_mm=((r+reference.d0)*40).tolist(),
            old_inlet_k_over_Kin=(old/reference.Kin).tolist(),
            candidate_counterfactual_k_over_Kin=(candidate/reference.Kin).tolist(),
            existing_v2_Dkw_at_inlet=diffusion2.tolist(), candidate_Dkw=diffusion3.tolist(),
            D_v2_over_candidate=(diffusion2/diffusion3).tolist(),
            old_inlet_cutoff_over_H=reference.old_c/column.p.H,
            active_v2_reference_cutoff_over_H=c2/column.p.H,
            candidate_reference_cutoff_over_H=c3/column.p.H,
            ground_D_ratio=float(diffusion2[0]/diffusion3[0]),
            ground_k_departure_over_Kin=float((candidate[0]-old[0])/reference.Kin),
            sampled_k_departure_MAX_over_Kin=float(np.max(abs(candidate-old))/reference.Kin),
            interpretation='Counterfactual versus unchanged inlet only; does not predict or qualify actual 2D development')
        report['conditional_integrations'] = gen.trials
        report['status'] = 'PASS_PROPOSAL_REFERENCE_CONTROLS_ONLY'
    except Exception as caught:
        error = caught
        report.update(exception=type(caught).__name__, reason=str(caught))
    report['elapsed_seconds'] = time.monotonic()-start
    with args.out.open('x') as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'], elapsed_seconds=report['elapsed_seconds'],
                         ground_D_ratio=report.get('inlet_diagnostic', {}).get('ground_D_ratio'),
                         candidate_ground_k_departure=report.get('inlet_diagnostic', {}).get('ground_k_departure_over_Kin'))))
    if error is not None:
        raise error


if __name__ == '__main__':
    main()
