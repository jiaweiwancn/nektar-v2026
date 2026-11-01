"""Laboratory-frame mean observables; methods do not use experimental targets."""
import math


def speedup(hill, floor, scale):
    if not all(math.isfinite(v) for v in (hill, floor, scale)) or scale <= 0:
        raise ValueError('Nonfinite speedup input/scale')
    if abs(floor) <= 1e-10*scale:
        return dict(available=False, reason='near_zero_flat_denominator')
    return dict(available=True, value=hill/floor-1)


def recovery_zeros(points):
    """Negative-to-positive recovery; never interpolate across missing probes."""
    points=sorted(points,key=lambda p:p[0])
    if len({x for x,_ in points})!=len(points):
        raise ValueError('Duplicate recovery coordinates')
    roots=[]
    for (x,a),(y,b) in zip(points,points[1:]):
        if a is None or b is None:
            continue
        if not all(math.isfinite(v) for v in (x,y,a,b)):
            raise ValueError('Nonfinite recovery trace')
        if a < 0 <= b:
            roots.append(x-a*(y-x)/(b-a))
    return roots


def intercept(points):
    """Unweighted declared tail fit x=a+b*y; intercept is at physical ground."""
    n=len(points)
    if n<2 or len({y for y,_ in points})!=n:
        raise ValueError('At least two distinct tail heights are required')
    mean_y=sum(y for y,_ in points)/n;mean_x=sum(x for _,x in points)/n
    slope=sum((y-mean_y)*(x-mean_x) for y,x in points)/sum((y-mean_y)**2 for y,_ in points)
    return mean_x-slope*mean_y


def reattachment(rows, hill_height=.04):
    """Ground extrapolation of U=0 dividing points, with two/three-tail sensitivity.

    This is distinct from a finite-probe recovery or a wall-shear zero. The
    lowest two available heights are the downstream tail; three are required
    for the declared method sensitivity. Multiple recoveries stay ambiguous.
    """
    groups={}
    for row in rows:
        if not row['id'].startswith('dividing:'):
            continue
        x,z=float(row['x']),float(row['z'])
        if x<=.1 or z<=0 or float(row['span'])!=0:
            raise ValueError('Dividing probes require positive vertical height on the flat wake x>0.1, span=0')
        groups.setdefault(z,[]).append((x,None if row['u']=='' else float(row['u'])))
    dividing=[];tail=[]
    for z,points in sorted(groups.items()):
        roots=recovery_zeros(points)
        dividing.append(dict(height_m=z,recovery_x_m=roots,missing_probes=sum(v is None for _,v in points)))
        if len(roots)==1:
            tail.append((z,roots[0]))
    result=dict(available=False,dividing_points=dividing,
                method='linear x(y) ground intercept using the two lowest unique-recovery heights; three-height sensitivity',
                interpretation='ground-extrapolated dividing line, not finite-probe zero, wall-shear zero or bubble length')
    # An unavailable/ambiguous low tail cannot be replaced silently by higher probes.
    if len(groups)<3 or len(tail)<3 or any(len(d['recovery_x_m'])!=1 or d['missing_probes'] for d in dividing[:3]):
        result['reason']='three complete, uniquely recovering low tail heights are required'
        return result
    two,three=intercept(tail[:2]),intercept(tail[:3])
    result.update(available=True,x_R_m=two,x_R_over_H=two/hill_height,
                  three_tail_x_R_over_H=three/hill_height,
                  tail_method_difference_over_H=abs(two-three)/hill_height,
                  tail_heights_m=[y for y,_ in tail[:3]])
    return result
