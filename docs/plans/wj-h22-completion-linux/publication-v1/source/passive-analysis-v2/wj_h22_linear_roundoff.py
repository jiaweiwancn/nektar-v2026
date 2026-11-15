# SPDX-License-Identifier: MIT
"""Propagate a certified FE-gradient roundoff bound through the H22 linear map.

128 eps and the original nonzero scales are unchanged. A composition cannot
discard the error already certified for its reconstructed input. Direct and
propagated certificates are both retained, including direct-bound exceedances.
"""
import numpy as np

def dev_bound(bound):
    b=np.asarray(bound)
    assert b.shape==(3,3)and np.all(np.isfinite(b))and np.all(b>=0)
    return b+np.eye(3)*np.trace(b)/3

def propagated_identity(residual,terms,inherited,amplitude=1.):
    r=np.asarray(residual);incoming=np.asarray(inherited)
    assert np.all(np.isfinite(incoming))and np.all(incoming>=0)
    direct=128*np.finfo(float).eps*(amplitude+sum(abs(np.asarray(x))for x in terms))
    bound=direct+incoming
    if np.any(~np.isfinite(r))or np.any(abs(r)>bound):
        raise ValueError('linear composed identity exceeds propagated primitive roundoff bound')
    return dict(status='PASS_PROPAGATED_PRIMITIVE_BOUNDS',residual=r.tolist(),bound=bound.tolist(),
        direct_bound=direct.tolist(),inherited_bound=incoming.tolist(),maximum_error=float(np.max(abs(r))),
        maximum_bound=float(np.max(bound)),direct_bound_exceeded_count=int(np.count_nonzero(abs(r)>direct)),
        interpretation='unchanged 128 eps primitive certificates plus their linear propagation; not a physical comparison floor')
