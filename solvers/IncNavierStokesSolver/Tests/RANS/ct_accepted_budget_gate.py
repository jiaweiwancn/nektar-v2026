# SPDX-License-Identifier: MIT
"""Complete physical-budget primitive admission; labels never supply coverage."""
import math

CONTROLS = {'native', 'q48f16', 'q64f16', 'q64f32'}
EQUATIONS = {'inner_k', 'outer_k', 'outer_epsilon', 'whole_k'}
RECONSTRUCTION = {'k', 'outer_epsilon'}


def finite_numbers(value):
    if isinstance(value, dict):return all(finite_numbers(v) for v in value.values())
    if isinstance(value, (list, tuple)):return all(finite_numbers(v) for v in value)
    return math.isfinite(value) if isinstance(value, (float, int)) else True


def physical_budgets(row):
    reasons, controls, values = [], {}, {}
    budgets = row.get('budgets', {})
    if set(budgets) != CONTROLS:
        reasons.append('missing_or_foreign_budget_control')
    for label, budget in budgets.items():
        errors = []
        if 'budget_binding' in row and budget.get('binding') != row['budget_binding']:
            errors.append('state_configuration_budget_identity_mismatch')
        if set(budget.get('equations', {})) != EQUATIONS:
            errors.append('missing_or_foreign_budget_equation')
        if set(budget.get('reconstruction', {})) != RECONSTRUCTION:
            errors.append('missing_or_foreign_reconstruction')
        for name, q in budget.get('equations', {}).items():
            if not finite_numbers(q):errors.append('nonfinite_equation_primitive:' + name)
            keys = ('physical_scale', 'storage_backward_difference', 'stage_source',
                    'stage_physical_boundary_flux', 'frozen_stage_physical_defect')
            if any(not isinstance(q.get(k), (float, int)) or not math.isfinite(q[k]) for k in keys):
                errors.append('missing_or_nonfinite_equation_primitive:' + name)
                continue
            scale = q['physical_scale']
            if scale <= 0:
                errors.append('nonpositive_physical_scale:' + name)
                continue
            defect = q['storage_backward_difference'] - q['stage_source'] - q['stage_physical_boundary_flux']
            if abs(defect-q['frozen_stage_physical_defect']) > 2e-10*max(1., scale):
                errors.append('inconsistent_physical_defect:' + name)
            if abs(defect)/scale > .02:
                errors.append('physical_budget_accuracy:' + name)
            values.setdefault(name, []).append(defect/scale)
        for name, q in budget.get('reconstruction', {}).items():
            if not finite_numbers(q):errors.append('nonfinite_reconstruction_primitive:' + name)
            v = q.get('scaled_defect')
            if 'normalization' in q:
                normal=q['normalization'];maximum=q.get('maximum_defect')
                if not isinstance(normal,(float,int)) or not math.isfinite(normal) or normal<=0 or not isinstance(maximum,(float,int)) or not math.isfinite(maximum) or maximum<0:
                    errors.append('invalid_reconstruction_scale:' + name)
                elif not isinstance(v,(float,int)) or abs(maximum/normal-v)>2e-14*max(1.,abs(v)):
                    errors.append('inconsistent_reconstruction_primitives:' + name)
                elif maximum/normal>2e-10:errors.append('reconstruction_accuracy:' + name)
            if not isinstance(v, (float, int)) or not math.isfinite(v) or v < 0 or v > 2e-10:
                errors.append('missing_or_failed_reconstruction_primitive:' + name)
        controls[label] = dict(status='pass' if not errors else 'fail', reasons=errors)
        reasons.extend(label+':'+e for e in errors)
    spreads = {n: max(v)-min(v) for n, v in values.items()}
    if any(v > .02 for v in spreads.values()):
        reasons.append('quadrature_or_face_uncertainty')
    return dict(status='pass' if not reasons else 'fail', reasons=reasons, controls=controls,
                scaled_signed_defect_spreads=spreads, uncertainty_limit=.02,
                required_controls=sorted(CONTROLS), required_equations=sorted(EQUATIONS),
                required_reconstruction=sorted(RECONSTRUCTION))


def paired_coverage(members, horizon, expected):
    """Bind required dts, time/state/config/runtime identities supplied by caller.

    Each member carries actual primitives and the separately derived expected
    identity. Initial vectors may differ between spaces, never within a dt pair.
    """
    reasons = []
    if set(members) != {'0.002', '0.001'}:
        reasons.append('missing_or_foreign_timestep_member')
    initials = []
    for dt, row in members.items():
        t = row.get('time')
        if not isinstance(t, (int, float)) or not math.isfinite(t) or abs(t-horizon) > 1e-13:
            reasons.append(dt+':incomplete_required_endpoint')
        if row.get('failure'):
            reasons.append(dt+':native_rejection')
        if row.get('identity') != expected.get(dt):
            reasons.append(dt+':state_configuration_identity_mismatch')
        initials.append(row.get('initial_state_sha256'))
        if physical_budgets(row)['status'] != 'pass':
            reasons.append(dt+':physical_budget_primitives')
    if not initials or not all(initials) or len(set(initials)) != 1:
        reasons.append('same_space_initial_coefficients_differ_or_missing')
    return dict(status='pass' if not reasons else 'fail', reasons=reasons,
                required_horizon=horizon, required_dt=[.002, .001])
