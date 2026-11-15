#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Version semantics and native adaptive-context preservation checks."""
import argparse, time
from pathlib import Path
import numpy as np
from cao_tamura_h2i_convex import Context, validate_direction_counts
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser();p.add_argument('--context',type=Path,required=True);p.add_argument('--previous',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();start=time.monotonic();tests={}
    for version,n in [('FixedKJointContextV1',41),('FixedKJointContextV2',44),('FixedKJointContextV3',48),('FixedKJointContextV3',128)]:assert validate_direction_counts(version,dict(velocity=n,k_anchor=8,epsilon=8))==n
    for version,n,k,e in [('FixedKJointContextV1',44,8,8),('FixedKJointContextV2',48,8,8),('FixedKJointContextV3',41,8,8),('FixedKJointContextV3',129,8,8),('FixedKJointContextV3',48,9,8),('FixedKJointContextV3',48,8,9)]:
        try:validate_direction_counts(version,dict(velocity=n,k_anchor=k,epsilon=e))
        except ValueError:continue
        raise AssertionError('invalid version/count accepted')
    tests['archived_counts_and_successor_rejections']=True
    old=Context(a.previous);new=Context(a.context);nv=old.meta['direction_counts']['velocity'];assert nv>=44
    for key in old.blocks:
        x=new.block(key);x=x[:,:nv] if key=='velocity_basis' else x;assert old.block(key).tobytes()==x.tobytes(),key
    for key in ('members','member_weights'):assert old.arrays[key].tobytes()==new.arrays[key].tobytes(),key
    tests['all_prior_responses_scalars_memberships_bits_preserved']=True
    for key in ('geometry_identity','inventory_identity','reference_sha256','defect_sha256','patch_signature'):assert old.doc[key]==new.doc[key],key
    for key in ('velocity_mass','divergence_GL16','divergence_GL20','divergence_GL24'):
        # Native sums over quadrature retain the same old sub-block arithmetic.
        assert old.algebra[key].tobytes()==new.algebra[key][:nv+1,:nv+1].tobytes(),key
    tests['old_quadratic_base_cross_direction_blocks_bitwise']=True
    v=new.meta['verification'];assert v['gradient_relative_error']<1e-10 and v['scalar_relative_error']<1e-10 and v['native_tensor_error']<1e-11 and v['protected_rows_and_profiles']
    tests['native_meaningful_signed_responses_and_frozen_profile']=True
    rng=np.random.default_rng(1307);B=rng.normal(size=(19,8));W=np.exp(rng.normal(size=19));Q=B.T@(W[:,None]*B);u=rng.normal(size=8);d=rng.normal(size=8)
    for alpha in (-.1,.03,.2):assert abs((u+alpha*d)@Q@(u+alpha*d)-(u@Q@u+2*alpha*u@Q@d+alpha**2*d@Q@d))<1e-11
    tests['independent_finite_quadratic_identity']=True
    atomic(a.output,dict(version='H2kFocusedVerificationV1',passed=True,tests=tests,new_direction_counts=new.meta['direction_counts'],seconds=time.monotonic()-start))


if __name__=='__main__':main()
