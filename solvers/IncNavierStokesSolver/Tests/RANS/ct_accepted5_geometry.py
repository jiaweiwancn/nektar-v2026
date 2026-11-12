# SPDX-License-Identifier: MIT
"""Closest point to the exact compact cosine-squared graph, in Hh units.

F(t)=t-x-h'(t)(z-h(t)) is half the distance-squared derivative.
For a=pi/5 and c=cos(2at), F'=1+a²+a²(2z-1)c-2a²c².
Its critical points can therefore be isolated without a sampling grid.
Every open subinterval is monotone; endpoint/critical zeros include tangent
roots. This gives all stationary candidates in exact arithmetic. Binary64
coalescence below the stated root tolerance is explicitly unresolved as a
multiplicity question; it does not turn a vertical distance into a clamp.
"""
import math,sys
from scipy.optimize import brentq

L=2.5
A=math.pi/5
EPS=sys.float_info.epsilon

def hill(t):return math.cos(A*t)**2 if abs(t)<L else 0.
def slope(t):return -A*math.sin(2*A*t) if abs(t)<L else 0.

def unique(values,tol=64*EPS):
    out=[]
    for v in sorted(values):
        if not out or abs(v-out[-1])>tol*max(1.,abs(v),abs(out[-1])):out.append(v)
    return out

def closest_normal(x,z):
    if not math.isfinite(x) or not math.isfinite(z):raise ValueError('nonfinite graph coordinate')
    vertical=z-hill(x)
    if vertical<0:return dict(status='inside_terrain',distance=None,method='ClosestGraphNormalV2')
    b=2*z-1;C=1+1/A**2
    if not math.isfinite(b):raise ValueError('graph coordinate exceeds binary64 isolation range')
    d=math.hypot(b,math.sqrt(8*C));q=(b+math.copysign(d,b))/2
    critical=[-L,0.,L];cos_roundoff=0
    for c in (q/2,-C/q):
        if -1<=c<=1:pass
        elif abs(abs(c)-1)<=32*EPS:c=math.copysign(1.,c);cos_roundoff+=1
        else:continue
        t=math.acos(c)/(2*A);critical.extend((-t,t))
    cuts=unique(critical);roots=[]
    def f(t):return t-x-slope(t)*(z-hill(t))
    def tolerance(t):return 64*EPS*(1+abs(t)+abs(x)+abs(slope(t)*(z-hill(t))))
    for t in cuts:
        if abs(f(t))<=tolerance(t):roots.append(t)
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        fl,fh=f(lo),f(hi)
        if fl*fh<0:roots.append(brentq(f,lo,hi,xtol=5e-15,rtol=8*EPS))
    roots=unique(roots)
    candidates=unique([-L,L,*roots,*([x] if abs(x)>=L else [])])
    t=min(candidates,key=lambda t:(t-x)**2+(hill(t)-z)**2)
    distance=math.hypot(t-x,hill(t)-z)
    bound=128*EPS*max(1.,abs(x),abs(z),distance,vertical)
    if distance>vertical+bound:raise ArithmeticError('closest-graph distance exceeds vertical candidate; root isolation failed')
    return dict(status='foot_limit' if abs(t)==L else 'closest_normal',footpoint_x=t,distance=distance,vertical=vertical,
        normal_root_count=len(roots),normal_roots=roots,critical_points=cuts,candidate_count=len(candidates),
        stationarity_residual=f(t),cosine_endpoint_roundoff=cos_roundoff,method='ClosestGraphNormalV2',
        arithmetic_scope='all monotone stationarity intervals; roots coalescing within64eps*coordinate scale are deduplicated; no separate multiplicity certificate below binary64 resolution')
