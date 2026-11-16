# SPDX-License-Identifier: MIT
"""Explicit scope for immutable H-plus trace addresses; no current values."""
import hashlib,json,math,sys
from ct_accepted_open import require

_active = None


def positive_h_trace_index(initial,x,height):
    catalog = _active
    if catalog is not None and catalog.initial is initial and catalog.height == height:
        catalog.queries += 1
        try:
            i = catalog.indices[x]
        except KeyError:
            raise StopIteration from None
        o = initial['observation']
        require(o['grid'][i] == [x,0.,height] and o['grid_sides'][i] == 1,
                'stale immutable positive-H trace index')
        return i
    o = initial['observation']
    return next(i for i,(q,side) in enumerate(zip(o['grid'],o['grid_sides'],strict=True))
                if q == [x,0.,height] and side == 1)


class TraceIndex:
    policy = 'CT16-ExactTraceIndexV1'
    maximum_grid_rows = 262144
    maximum_entries = 8192
    maximum_index_bytes = 64*1024**2

    def __init__(self,initial,height):
        self.initial = initial; self.height = height; self.queries = 0; self.previous = None
        o = initial['observation']; grid = o['grid']; sides = o['grid_sides']
        require(math.isfinite(height) and height > 0 and len(grid) == len(sides) == len(initial['test_basis'])
                and len(grid) <= self.maximum_grid_rows, 'invalid or oversized immutable trace packet')
        self.indices = {}
        for i,(q,side) in enumerate(zip(grid,sides,strict=True)):
            require(len(q) == 3 and all(math.isfinite(x) for x in q), 'nonfinite immutable trace coordinate')
            if q[1] == 0 and q[2] == height and side == 1:
                self.indices.setdefault(q[0],i)  # same first matching row as next()
        require(len(self.indices) <= self.maximum_entries, 'immutable trace entry bound')
        self.bytes = sys.getsizeof(self.indices)+sum(sys.getsizeof(k)+sys.getsizeof(v) for k,v in self.indices.items())
        require(self.bytes <= self.maximum_index_bytes, 'immutable trace byte bound')
        self.geometry_sha256 = self.geometry_digest()

    def geometry_digest(self):
        o = self.initial['observation']
        return hashlib.sha256(json.dumps([o['grid'],o['grid_sides'],len(self.initial['test_basis']),self.height],
                             separators=(',',':'),allow_nan=False).encode()).hexdigest()

    def __enter__(self):
        global _active
        self.previous = _active; _active = self
        return self

    def __exit__(self,*args):
        global _active
        try:
            require(_active is self and self.geometry_digest() == self.geometry_sha256,
                    'immutable trace geometry changed inside validation scope')
        finally:
            _active = self.previous
