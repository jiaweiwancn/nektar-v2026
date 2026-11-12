# SPDX-License-Identifier: MIT
"""Explicit read-only bridge for the user-directed post-reboot disk relocation."""
import json
from pathlib import Path
from rans_request_manifest import sha
OLD='/media/jay/新加卷/'
NEW='/media/jay/新加卷1/'

def artifact(path,expected_sha256):
    p=Path(path)
    if not p.exists() and str(p).startswith(OLD):p=Path(NEW+str(p)[len(OLD):])
    if sha(p)!=expected_sha256:raise ValueError('relocated artifact bytes differ: '+str(p))
    return p

def canonical_location(value):
    # Only the approved disk prefix is relocated. Every digest, other path,
    # library member and physics/method identity remains exact.
    if isinstance(value,str):return OLD+value[len(NEW):] if value.startswith(NEW) else value
    if isinstance(value,dict):return {canonical_location(k):canonical_location(v) for k,v in value.items()}
    if isinstance(value,list):return [canonical_location(v) for v in value]
    return value
