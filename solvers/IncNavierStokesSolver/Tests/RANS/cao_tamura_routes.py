"""Explicit experimental wall routes and independent native mask/coverage audit."""
import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET

ROUTES=('EPTransition','UniformEP','Original','CompatibleLog')


def resolve_route(model, route=None):
    ep=model.endswith(('_fixed','_switched'))
    expected='CompatibleLog' if model=='standard' else 'Original' if model.endswith('_original') else 'EPTransition'
    route=expected if route is None else route
    if route not in ROUTES or (ep and route not in ('EPTransition','UniformEP')) or (not ep and route!=expected):
        raise ValueError('Wall route is incompatible with the model/policy: '+str(route))
    return route


def envelope(x, begin, end):
    if not math.isfinite(begin+end) or end<=begin:raise ValueError('Invalid transition endpoints')
    r=max(0.,min(1.,(x-begin)/(end-begin)))
    return 1.-10*r**3+15*r**4-6*r**5


def route_metadata(case):
    root=ET.parse(case).getroot()
    info={n.attrib['PROPERTY']:n.attrib['VALUE'] for n in root.findall('CONDITIONS/SOLVERINFO/I')}
    params={n.text.split('=',1)[0].strip():float(n.text.split('=',1)[1])
            for n in root.findall('CONDITIONS/PARAMETERS/P')}
    wall=info.get('RANSWallModel','None');spatial=info.get('RANSAWFSpatialWeight','None')
    if wall=='AWFEquilibrium':
        if spatial not in ('None','Ramp'):raise ValueError('Unsupported experimental envelope')
        kind='EPTransition' if spatial=='Ramp' else 'UniformEP'
    else:kind='Original' if wall=='AWFOriginal' else 'CompatibleLog'
    if info.get('RANSExperimentalRoute',kind)!=kind:raise ValueError('Declared route differs from the actual XML envelope')
    out=dict(version='cao-tamura-wall-route-v1',kind=kind,wall_model=wall,
        weight_mode=info.get('RANSAWFWeight','Switched') if wall=='AWFEquilibrium' else None,
        fixed_weight=params.get('RANSAWFFixedWeight',1.) if wall=='AWFEquilibrium' else None,
        failure_policy=info.get('RANSAWFFailurePolicy','UseOriginal'))
    if kind=='EPTransition':
        a,b=params['RANSAWFTransitionStart'],params['RANSAWFTransitionEnd']
        envelope(a,a,b);out.update(begin=a,end=b)
    return out


def coverage(case):
    """Strict R1 fixed-EP audit; no silent fallback in the requested region.

    Counts retain duplicated one-sided native traces and are not area weights.
    Prescribed downstream zero correction is explicitly distinguished from a
    failed upstream correction. Check every accepted post-startup native time.
    """
    case=Path(case);route=route_metadata(case)
    if route['kind'] not in ('EPTransition','UniformEP') or route['weight_mode']!='Fixed':
        raise ValueError('R1 coverage requires an explicitly fixed EP route')
    groups={}
    for file in case.parent.glob(case.stem+'.terrain-points-r*.csv'):
        with file.open(newline='') as stream:
            for r in csv.DictReader(stream):
                if r['phase']=='1':groups.setdefault(float(r['time']),[]).append(r)
    if not groups:raise ValueError('No accepted native wall samples')
    expected=set();partial=0
    with case.with_suffix('.rans.csv').open(newline='') as stream:
        for r in csv.DictReader(stream):
            try:valid=None not in r and all(math.isfinite(float(v)) for v in r.values())
            except (TypeError,ValueError):valid=False
            if valid:expected.add(float(r['time']))
            else:partial+=1
    if not expected:raise ValueError('No complete finite accepted solver stages')
    if expected-set(groups):raise ValueError('Missing wall coverage at accepted solver times: '+str(sorted(expected-set(groups))))
    # An aborted process can finish its wall output immediately before its
    # stage summary is written. Such extra frames cannot certify acceptance.
    ignored=sorted(set(groups)-expected)
    groups={t:rows for t,rows in groups.items() if t in expected}
    maxima=dict(mask=0.,requested=0.,active=0.);fallback=0;counts={}
    for time,rows in groups.items():
        bins=dict(upstream=0,transition=0,downstream=0)
        for r in rows:
            x=float(r['x'])
            if not math.isfinite(x):raise ValueError('Nonfinite native wall position')
            w=envelope(x,route['begin'],route['end']) if route['kind']=='EPTransition' else 1.
            request=w*route['fixed_weight'];actual=float(r['chi'])
            values=[float(r[k]) for k in ('mask','requested_chi','chi')]
            if not all(math.isfinite(v) for v in values):raise ValueError('Nonfinite wall coverage')
            maxima['mask']=max(maxima['mask'],abs(values[0]-w))
            maxima['requested']=max(maxima['requested'],abs(values[1]-request))
            maxima['active']=max(maxima['active'],abs(actual-request))
            if request>1e-12 and int(r['reason'])>=2:fallback+=1
            bins['upstream' if w>=1.-1e-12 else 'downstream' if w<=1e-12 else 'transition']+=1
        if not bins['upstream'] or (route['kind']=='EPTransition' and (not bins['transition'] or not bins['downstream'])):
            raise ValueError('Native grid does not sample every wall-route region')
        counts[str(time)]=bins
    if abs(route['fixed_weight']-1.)>1e-12 or max(maxima.values())>5e-9 or fallback:
        raise ValueError('Fixed EP route lost its prescribed correction: '+str(dict(errors=maxima,fallback=fallback)))
    return dict(route=route,status='PASS',accepted_times=len(groups),time_begin=min(groups),time_end=max(groups),
        expected_solver_times=len(expected),ignored_nonstage_wall_times=ignored,incomplete_solver_rows=partial,
        maximum_errors=maxima,fallback_native_traces=fallback,region_counts=counts,
        rule='actual chi equals spatial envelope times fixed weight at every accepted trace; no admissibility fallback',
        weighting='native counts including one-sided copies, not an area fraction')
