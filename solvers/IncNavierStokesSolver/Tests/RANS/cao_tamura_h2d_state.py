# SPDX-License-Identifier: MIT
"""Binary64 numerical JSON contract: the native writer can emit integer -0."""
import json,struct
from pathlib import Path

def load(path):
    with Path(path).open() as f:return json.load(f,parse_int=lambda s:-0.0 if s=='-0' else int(s))

def bits(value):return struct.pack('<d',value)

def rebind_fixed(base_document,candidate):
    """Repair only serialization of frozen rows; preserve optimized free values.

    No optimization, projection, Lift or change in any numerical fixed value.
    The native reader independently checks every restored fixed coefficient.
    """
    f=base_document['feasibility'];base=f['base'];free=set(map(int,f['free_velocity_rows']));state=candidate['state'][:];changed=[]
    for i,b in enumerate(base):
        if i not in free:
            if state[i]!=b:raise ValueError('fixed numerical value changed '+str(i))
            if bits(state[i])!=bits(b):changed.append(i)
            state[i]=b
    return state,changed
