# SPDX-License-Identifier: MIT
"""Finite cumulative WP2 resource amendments; the original clock never resets."""
import hashlib,json
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha


def wall_limit(root):
    path=root/'resource-extension-13h.json'
    if not path.exists():return 43200
    amendment=load(path)
    require(amendment['version']=='WP2CumulativeResourceAmendmentV1' and
            amendment['original_wall_seconds']==43200 and amendment['cumulative_wall_seconds']==46800,
            'unsupported WP2 resource amendment')
    require(amendment['origin_sha256']==sha(root/'origin.json'),'WP2 clock origin changed')
    original=load(root/'source/docs/plans/cls-cao-tamura-budget-wp2-evidence/policy.json')
    expected=hashlib.sha256(json.dumps(original['acceptance'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
    require(amendment['unchanged_acceptance_sha256']==expected,'WP2 acceptance criteria changed')
    successor=root/'resource-extension-16h.json'
    if successor.exists():
        extension=load(successor)
        require(extension['version']=='WP2CumulativeResourceAmendmentV2' and
                extension['original_wall_seconds']==43200 and
                extension['previous_cumulative_wall_seconds']==46800 and
                extension['cumulative_wall_seconds']==57600 and
                extension['additional_physical_attempt_limit']==180,
                'unsupported WP2 global-control amendment')
        require(extension['predecessor_sha256']==sha(path),'WP2 predecessor amendment changed')
        require(extension['origin_sha256']==sha(root/'origin.json'),'WP2 clock origin changed')
        require(extension['unchanged_acceptance_sha256']==expected,'WP2 acceptance criteria changed')
        confirmation=root/'resource-extension-24h.json'
        if confirmation.exists():
            control=load(confirmation)
            require(control['version']=='WP2CumulativeResourceAmendmentV3' and
                    control['original_wall_seconds']==43200 and control['previous_cumulative_wall_seconds']==57600 and
                    control['cumulative_wall_seconds']==86400 and control['additional_physical_attempt_limit']==120 and
                    control['reserved_attempts_before']==1646,
                    'unsupported WP2 confirmation-recovery amendment')
            require(control['predecessor_sha256']==sha(successor),'WP2 predecessor amendment changed')
            require(control['origin_sha256']==sha(root/'origin.json'),'WP2 clock origin changed')
            require(control['unchanged_acceptance_sha256']==expected,'WP2 acceptance criteria changed')
            pressure=root/'resource-extension-pressure-28h.json'
            if pressure.exists():
                last=load(pressure)
                require(last['version']=='WP2CumulativeResourceAmendmentV4' and
                        last['original_wall_seconds']==43200 and last['previous_cumulative_wall_seconds']==86400 and
                        last['cumulative_wall_seconds']==100800 and last['reserved_attempts_before']==1766 and
                        last['additional_physical_attempt_limit']==1102,
                        'unsupported WP2 pressure resource amendment')
                require(last['predecessor_sha256']==sha(confirmation) and last['origin_sha256']==sha(root/'origin.json'),
                        'WP2 pressure resource predecessor or clock changed')
                require(last['unchanged_acceptance_sha256']==expected and
                        last['pressure_scope_sha256']==sha(root/'pressure-control-scope.json'),
                        'WP2 pressure resource acceptance or numerical scope changed')
                return last['cumulative_wall_seconds']
            return control['cumulative_wall_seconds']
        return extension['cumulative_wall_seconds']
    return amendment['cumulative_wall_seconds']


def pressure_scope(root):
    """Prospective numerical control under the unchanged cumulative deadline."""
    path=root/'pressure-control-scope.json'
    if not path.exists():return None
    d=load(path);previous=root/'resource-extension-24h.json'
    require(d['version']=='WP2PressureControlScopeV1' and
            d['cumulative_wall_seconds']==86400 and wall_limit(root) in (86400,100800) and
            d['reserved_attempts_before']==1766 and d['additional_physical_attempt_limit']==742 and
            d['spatial_case_limit']==40 and d['physical_attempt_limit_unchanged']==6000,
            'unsupported WP2 pressure-control scope')
    require(d['predecessor_sha256']==sha(previous) and d['origin_sha256']==sha(root/'origin.json'),
            'WP2 pressure-control origin or predecessor changed')
    require(d['unchanged_acceptance_sha256']==load(previous)['unchanged_acceptance_sha256'],
            'WP2 pressure-control acceptance changed')
    if (root/'resource-extension-pressure-28h.json').exists():
        last=load(root/'resource-extension-pressure-28h.json')
        d=dict(d,additional_physical_attempt_limit=last['additional_physical_attempt_limit'],
               effective_cumulative_wall_seconds=last['cumulative_wall_seconds'])
    return d
