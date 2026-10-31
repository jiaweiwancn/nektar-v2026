"""Independent collocation and side-flux audit for uniform terrain fixtures.

Uses saved wall scalars, analytical geometry and independently constructed GLL
operators; never calls the solver's surface differentiation or integration API.
"""
import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def gll(n):
    degree=n-1
    def legendre(x):
        a,b=1.,x
        for k in range(2,degree+1): a,b=b,((2*k-1)*x*b-(k-1)*a)/k
        return b,a
    nodes=[-math.cos(math.pi*i/degree) for i in range(n)]
    for i in range(1,n-1):
        x=nodes[i]
        for _ in range(30):
            p,previous=legendre(x)
            first=degree*(x*p-previous)/(x*x-1)
            second=(2*x*first-degree*(degree+1)*p)/(1-x*x)
            change=first/second;x-=change
            if abs(change)<2e-16: break
        nodes[i]=x
    barycentric=[1/math.prod(x-y for j,y in enumerate(nodes) if i!=j) for i,x in enumerate(nodes)]
    matrix=[]
    for i,x in enumerate(nodes):
        row=[barycentric[j]/(barycentric[i]*(x-y)) if j!=i else 0. for j,y in enumerate(nodes)]
        row[i]=-sum(row);matrix.append(row)
    weights=[2/(degree*(degree+1)*legendre(x)[0]**2) for x in nodes]
    return nodes,weights,matrix


def audit(directory):
    directory=Path(directory);root=ET.parse(directory/'neutral.xml').getroot()
    params={p.text.split('=',1)[0].strip():float(p.text.split('=',1)[1]) for p in root.findall('CONDITIONS/PARAMETERS/P')}
    info={p.attrib['PROPERTY']:p.attrib['VALUE'] for p in root.findall('CONDITIONS/SOLVERINFO/I')}
    dim=int(root.find('GEOMETRY').attrib['DIM'])
    count=int(root.find('EXPANSIONS/E').attrib['NUMPOINTS'].split(',')[0])
    _,weights,D=gll(count);face_size=count**(dim-1)
    all_points=[]
    for path in directory.glob('neutral.terrain-points-r*.csv'):
        all_points += [{k:float(v) for k,v in r.items()} for r in csv.DictReader(path.open()) if float(r['phase'])==1]
    if not all_points: raise RuntimeError('No accepted wall points for continuity audit')
    groups={}
    for row in all_points: groups.setdefault((row['time'],int(row['point'])//face_size),[]).append(row)
    log=info.get('RANSWallModel','None')=='None'
    errors=dict(target=0.,normalP=0.,target_flux_plus_side_flux=0.,actual_flux_plus_side_flux=0.,geometry=0.)
    budgets={};maximum_transpiration=0.
    centre,half=params.get('RANSTerrainCentre',0.),params.get('RANSTerrainHalfLength',.1)
    amplitude=params.get('RANSTerrainHeight',0.) if info.get('RANSWallGeometry')=='CosineHill' else 0.
    for (time,_),points in groups.items():
        points.sort(key=lambda r:r['point'])
        if len(points)!=face_size: raise RuntimeError('Incomplete native surface face')
        def at(i,j=0): return points[i+count*j]
        stream_axis=0
        if dim==3 and abs(at(count-1)['x']-at(0)['x'])<abs(at(count-1)['span']-at(0)['span']): stream_axis=1
        def point(s,y=0): return at(s,y) if stream_axis==0 else at(y,s)
        dx=(point(count-1)['x']-point(0)['x'])/2
        dy=(point(0,count-1)['span']-point(0,0)['span'])/2 if dim==3 else 1.
        if abs(dx)<1e-12 or (dim==3 and abs(dy)<1e-12): raise RuntimeError('Invalid surface coordinate orientation')
        middle=(point(0)['x']+point(count-1)['x'])/2
        local_target=local_actual=side=0.
        for j in range(count if dim==3 else 1):
            side += params['H']*(point(count-1,j)['uP']-point(0,j)['uP'])*math.copysign(1,dx)*(abs(dy)*weights[j] if dim==3 else 1)
            for i in range(count):
                r=point(i,j);x=r['x'];angle=math.pi*(x-centre)/half
                slope=params.get('RANSTerrainSlope',0.) if info.get('RANSWallGeometry')=='Plane' else 0.
                curvature=0.
                if amplitude and abs(middle-centre)<half:
                    slope=-amplitude*math.pi/(2*half)*math.sin(angle)
                    curvature=-amplitude*(math.pi/half)**2/2*math.cos(angle)/(1+slope*slope)**1.5
                arc=math.hypot(1,slope);H,P=r['H'],r['P'];JH,JP=1-curvature*H,1-curvature*P
                V=lambda z:z-curvature*z*z/2
                ds=sum(D[i][k]*point(k,j)['uP'] for k in range(count))/(dx*arc)
                dz=sum(D[j][k]*point(i,k)['spanP'] for k in range(count))/dy if dim==3 else 0.
                target=0. if log else -(H*ds+V(H)*dz)/JH
                normalP=0. if log else -(P*ds+V(P)*dz)/JP
                errors['target']=max(errors['target'],abs(target-r['normal_target']))
                if not log: errors['normalP']=max(errors['normalP'],abs(normalP-r['wP']))
                errors['geometry']=max(errors['geometry'],abs(JH-r['JH']),abs(V(H)-r['V']))
                measure=abs(dx)*arc*weights[i]*(abs(dy)*weights[j] if dim==3 else 1)
                local_target += r['normal_target']*JH*measure
                local_actual += (r['normal_target']+r['normal_error'])*JH*measure
                maximum_transpiration=max(maximum_transpiration,abs(r['normal_target']))
        if dim==3:
            for i in range(count):
                r=point(i,0);c=r['curvature'];arc=math.hypot(1,params.get('RANSTerrainSlope',0.))
                if amplitude and abs(middle-centre)<half:
                    arc=math.hypot(1,-amplitude*math.pi/(2*half)*math.sin(math.pi*(r['x']-centre)/half))
                side += (r['H']-c*r['H']**2/2)*(point(i,count-1)['spanP']-point(i,0)['spanP'])*math.copysign(1,dy)*abs(dx)*arc*weights[i]
        if log: side=0.
        budget=budgets.setdefault(time,dict(target_flux=0.,actual_flux=0.,side_flux=0.))
        budget['target_flux']+=local_target;budget['actual_flux']+=local_actual;budget['side_flux']+=side
    for budget in budgets.values():
        errors['target_flux_plus_side_flux']=max(errors['target_flux_plus_side_flux'],abs(budget['target_flux']+budget['side_flux']))
        errors['actual_flux_plus_side_flux']=max(errors['actual_flux_plus_side_flux'],abs(budget['actual_flux']+budget['side_flux']))
    # The collocation formula and integrated conservation are independent gates.
    # Geometry interpolation is finite-order, so retain an explicit geometry
    # tolerance while keeping the physical flux defect at the declared 1e-7.
    if max(errors['target'],errors['normalP'])>5e-8 or errors['geometry']>1e-8 or errors['target_flux_plus_side_flux']>5e-9 or errors['actual_flux_plus_side_flux']>1e-7:
        raise RuntimeError('Independent curved strip-continuity audit failed: '+str(errors))
    return dict(errors=errors,maximum_transpiration=maximum_transpiration,budgets=budgets)
