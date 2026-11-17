#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native row-schema compatibility and pressure-only polynomial fixture."""
import argparse
from copy import deepcopy
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import require, sha
from wp2_all_state import Rows


def main():
    parser = argparse.ArgumentParser()
    for key in ('root', 'probe', 'output'):
        parser.add_argument('--'+key, type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    old = load(args.root/'geometry-rows-selected/result.json')
    source = load(args.root/'geometry-rows-selected/input.json')
    xml = Path(load(args.root/'geometry-rows-selected/receipt.json')['command'][1])
    exported = {}
    for name, version in (('default', 'P1V1'), ('p1_v2', 'RegionalP2V1')):
        packet = dict(source, pressure_rows=version)
        data, receipt = run_probe(args.probe, xml, packet, args.output/name)
        for key, value in old['arrays'].items():
            require(value['sha256'] == data['arrays'][key]['sha256'], 'archived row bytes changed: '+key)
        exported[name] = dict(path=receipt['output'], receipt=receipt)
    context_path = args.root/'preservation-recovery-standard/step120/case.evolution.json'
    context = load(context_path)
    baseline = Rows(args.root/'geometry-rows-selected/result.json')
    successor = Rows(Path(exported['p1_v2']['path']))
    left, right = baseline.tuple(context['state'], context), successor.tuple(context['state'], context)
    for key in ('values', 'R', 'T', 'G', 'gk', 'ge', 'gp', 'endpoint_K'):
        require(left[key].tobytes() == right[key].tobytes(), 'zero-enrichment pressure replay changed '+key)
    native_replay = successor.check_native(context)
    candidates = list((args.root/'pressure-case-generation-checks/standard/inner-p2').glob('*.xml'))
    require(len(candidates) == 1, 'prospective pressure XML is ambiguous')
    p2_xml = candidates[0]
    p2_packet = dict(operation='wp2_geometry_rows', inner_stress_model='Linear', expected_model='StandardKEpsilon',
        xml_sha256=sha(p2_xml), z_partition=old['physical_z_partition'], pressure_rows='RegionalP2V1',
        pressure_polynomial_fixture=True)
    data, receipt = run_probe(args.probe, p2_xml, p2_packet, args.output/'mixed_p2')
    rows = Rows(Path(receipt['output']))
    require(rows.size == 9093 and data['sizes'][-3:-1] == [611, 140], 'native pressure dimensions differ from prediction')
    require(np.array_equal(rows.a['p_counts'][:, 0], np.where(rows.inner, 9, 4)), 'regional row widths are wrong')
    fixture = data['pressure_polynomial_fixture']
    require(fixture['governing_solves'] == 0 and not fixture['candidate_state'], 'pressure fixture masquerades as physical initialization')
    state = np.zeros(rows.size)
    for region in fixture['regions']:
        offset = int(region['offset'])
        state[offset:offset+len(region['coefficients'])] = region['coefficients']
    x = (rows.physical[:, 0]+.2)/.8
    height = np.where(rows.inner, .0065, .8-.0065)
    z = (rows.physical[:, 2]-np.where(rows.inner, 0., .0065))/height
    quadratic = rows.inner.astype(float)
    exact = .003+.001*x+.0004*z+.0002*x*z+quadratic*(.0007*x*x+.0003*z*z)
    dx = (.001+.0002*z+quadratic*.0014*x)/.8
    dz = (.0004+.0002*x+quadratic*.0006*z)/height
    errors = [float(max(abs(rows.evaluate(state, 'p', kind)-target)))
              for kind, target in [('phi', exact), ('dx', dx), ('dz', dz)]]
    require(errors[0] < 2e-12 and max(errors[1:]) < 2e-9, 'independent pressure polynomial/gradient replay failed')
    bad = deepcopy(p2_packet)
    bad.update(pressure_rows='P1V1', pressure_polynomial_fixture=False)
    log, negative = run_probe(args.probe, p2_xml, bad, args.output/'reject_mislabelled_p2', success=False)
    require('unsupported native pressure degree' in log, 'wrong rejection reason for mislabelled P2')
    require(not list((args.output/'reject_mislabelled_p2').glob('*.bin')), 'invalid schema wrote partial binary rows')
    atomic(args.output/'result.json', dict(version='WP2RegionalPressureRowsChecksV1', passed=True,
        default_arrays_bitwise_archived=True, zero_enrichment_fields_and_tensors_bitwise=True,
        native_tuple_replay=native_replay, native_mixed_dimensions=data['sizes'], native_slots=rows.size,
        actual_pressure_widths=[9, 4], independent_polynomial_errors=errors,
        native_polynomial_errors=[{k:v for k,v in region.items() if k != 'coefficients'} for region in fixture['regions']],
        auxiliary_pressure_mass_projections=2, governing_solves=0, physical_candidates=0,
        exports=exported, mixed_receipt=receipt, mislabelled_schema_rejected=negative,
        source_sha256=sha(Path(__file__)), probe_sha256=sha(args.probe),
        scope='Pressure-only synthetic projection/row verification. Native cold-state audit, complete pressure action, '
              'rank/nullspace, startup and recovery are UNEXECUTED for the new pressure representation.'))


if __name__ == '__main__':
    main()
