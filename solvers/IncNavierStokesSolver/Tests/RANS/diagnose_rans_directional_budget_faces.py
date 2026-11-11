#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only physical-face decomposition; require unchanged sealed budget totals."""
import argparse, copy, gzip, os, shutil
from pathlib import Path
from qualify_rans_preservation import load
from rans_duration_contract import MEMBERS, require
from rans_request_manifest import RequestManifest, runtime_identity, sha
from rans_qualification import THREADS, exact_vector_hash, write_json
from rans_preservation_oracle import Reference
from diagnose_rans_balanced_spatial import run_probe
from analyze_rans_physical_budgets import whole_patch


COMPONENTS = ('k_advective_outward', 'k_frozen_diffusive_outward',
              'k_current_diffusive_outward', 'epsilon_advective_outward',
              'epsilon_frozen_diffusive_outward', 'epsilon_current_diffusive_outward')


def decompose(receipt, ref):
    budget = receipt['physical_budget']
    grouped = {}
    max_point_defect = 0.
    for face in budget['physical_faces']:
        require(face['component_schema'] == 'SeparatePhysicalFaceFluxV1', 'face_component_schema')
        region, axis, side = (face[k] for k in ('region', 'axis', 'outward_side'))
        interface = axis == 2 and abs(face['coordinate'][2]-ref.H) < 1e-13
        key = f'region-{region}:axis-{axis}:side-{side}:interface-{int(interface)}'
        row = grouped.setdefault(key, dict(region=region, axis=axis, side=side,
            internal_interface=interface, samples=0, weight=0., **dict.fromkeys(COMPONENTS, 0.)))
        row['samples'] += 1
        row['weight'] += face['weight']
        for field in ('k', 'epsilon') if region else ('k',):
            for stage in ('frozen', 'current'):
                independent = face[field+'_'+stage+'_diffusive_outward']-face[field+'_advective_outward']
                max_point_defect = max(max_point_defect, abs(independent-face[field+'_'+stage+'_flux']))
        for term in COMPONENTS:
            if face[term] is not None:
                row[term] += face['weight']*face[term]
    equations = {}
    max_integral_defect = 0.
    for equation in ('inner_k', 'outer_k', 'outer_epsilon', 'whole_patch_k'):
        whole = equation == 'whole_patch_k'
        region = 0 if equation == 'inner_k' else 1
        field = 'epsilon' if equation == 'outer_epsilon' else 'k'
        rows = [v for v in grouped.values() if (whole or v['region'] == region)
                and not (whole and v['internal_interface'])]
        values = whole_patch(receipt, ref) if whole else receipt['scaled'][equation]
        volume = sum(budget['equations'][q]['volume'] for q in ('inner_k', 'outer_k')) if whole else budget['equations'][equation]['volume']
        scale = (ref.A/ref.H if field == 'epsilon' else ref.K)*ref.us/ref.H*volume
        adv = sum(v[field+'_advective_outward'] for v in rows)/scale
        frozen = sum(v[field+'_frozen_diffusive_outward'] for v in rows)/scale
        current = sum(v[field+'_current_diffusive_outward'] for v in rows)/scale
        defect = max(abs(frozen-adv-values['stage_physical_boundary_flux']),
                     abs(current-adv-values['instantaneous_physical_boundary_flux']))
        max_integral_defect = max(max_integral_defect, defect)
        equations[equation] = dict(values, advective_outward=adv,
            frozen_diffusive_outward=frozen, current_diffusive_outward=current,
            physical_flux_component_reconstruction_defect=defect,
            scale=scale, physical_face_groups=[key for key, value in grouped.items() if value in rows])
    require(max_point_defect <= 5e-14, 'physical_face_components_do_not_reproduce_native_flux')
    require(max_integral_defect <= 5e-12, 'physical_face_groups_do_not_reproduce_native_budget')
    return dict(equations=equations, physical_face_groups=grouped,
        maximum_point_flux_defect=max_point_defect, maximum_scaled_integral_defect=max_integral_defect)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--cohort', required=True)
    p.add_argument('--runtime-name', required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); r = a.root.resolve(); runtime = r/a.runtime_name
    os.environ.update(THREADS, LD_LIBRARY_PATH=str(runtime/'lib'), TMPDIR=str(r/'tmp'),
                      PYTHONDONTWRITEBYTECODE='1')
    inputs = {}
    for kind in ('baseline', 'perturbed'):
        for member in MEMBERS:
            leaf = r/'physical-budgets'/(a.cohort+'-'+kind+'-'+member+'-0.016-Q2-F32')
            source = load(leaf/'request-manifest.json')['request']
            require(RequestManifest(leaf, source, True).completed(leaf), 'source_budget_not_sealed')
            inputs[kind+':'+member] = dict(source=str(leaf), result_sha256=sha(leaf/'results.json'),
                packet_sha256=sha(leaf/'packet.json.gz'), request_sha256=sha(leaf/'request.json'))
    intent = dict(schema='DirectionalPhysicalFaceDecompositionV1', time=.016,
        inputs=inputs, diagnostic_runtime=runtime_identity(runtime/'RANSPreservationProbe'),
        policy='same sealed states, coefficients, consumed boundaries and physical points; extra read-only component output only; no new evolution or changed budget classification')
    seal = RequestManifest(a.output, intent)
    results = {}
    for key, item in inputs.items():
        leaf = Path(item['source']); source = load(leaf/'results.json')
        ref = Reference(*(source['configuration'][name] for name in ('model', 'policy', 'wall')))
        out = a.output/key.replace(':', '-'); out.mkdir()
        shutil.copy2(leaf/'case.xml', out/'case.xml')
        request = load(leaf/'request.json')
        packet, seconds = run_probe(runtime, out, request)
        require(exact_vector_hash(packet['observation']['state']) == source['state_sha256'], 'face_audit_changed_state')
        extra = packet['physical_scalar_budget']
        stripped = copy.deepcopy(extra)
        for face in stripped['physical_faces']:
            for component in (*COMPONENTS, 'component_schema'):
                del face[component]
        require(stripped == source['physical_budget'], 'probe_only_face_export_changed_old_budget')
        detailed = dict(source, physical_budget=extra)
        results[key] = dict(decompose(detailed, ref), seconds=seconds,
            source=item, packet_sha256=sha(out/'packet.json.gz'), state_sha256=source['state_sha256'],
            previous_state_sha256=source['previous_state_sha256'], unchanged_old_budget_exact=True)
    responses = {}
    for member in MEMBERS:
        before, after = (results[kind+':'+member] for kind in ('baseline', 'perturbed'))
        responses[member] = {equation: {name: value-before['equations'][equation][name]
            for name, value in values.items() if isinstance(value, (float, int)) and name != 'scale'}
            for equation, values in after['equations'].items()}
    result = dict(status='pass', reconstruction_status='pass', accuracy_status='not_reclassified',
        results=results, baseline_subtracted_components=responses,
        scope='independently evaluated physical advective/diffusive faces, including nonzero natural-wall gradients; interface sides separate and excluded from whole-patch external totals; original total-budget accuracy is unchanged')
    write_json(a.output/'results.json', result)
    seal.seal(a.output, [q for q in a.output.rglob('*') if q.is_file() and q.name not in ('request-manifest.json', 'completed-artifacts.json')])
    print('physical-face decomposition reproduced every sealed budget exactly')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
