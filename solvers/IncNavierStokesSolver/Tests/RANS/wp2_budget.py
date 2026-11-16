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
        return extension['cumulative_wall_seconds']
    return amendment['cumulative_wall_seconds']
