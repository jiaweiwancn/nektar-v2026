#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Full-size matched hill/flat mesh and guarded native capacity preflight."""
import argparse,copy,json,math,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import mpmath as mp
from generate_m3 import mesh,element
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_parameter,set_info,write
from ct_accepted_open import probe,require
from rans_request_manifest import RequestManifest,sha,runtime_identity

L=2.5;Z=20.;NC=.3;NR=1.;H=.15
X=[x/.04 for x in (-1.,-.8,-.68,-.56,-.52,-.48,-.44,-.40,-.36,-.32,-.28,-.24,-.20,-.16,-.12,-.1,-.075,-.05,-.025,0.,.025,.05,.075,.1,.14,.2,.3,.45,.6,.8,1.,1.2)]


def graph(q,amplitude=1.):
    if abs(q)>=L:return 0.,0.,0.
    t=math.pi*q/L
    return amplitude*(1+math.cos(t))/2,-amplitude*math.pi/(2*L)*math.sin(t),-amplitude*math.pi**2/(2*L**2)*math.cos(t)


def mapping(q,r,amplitude):
    h,p,_=graph(q,amplitude);a=math.hypot(1,p);normal=(q-r*p/a,h+r/a)
    vertical=(q,h+(Z-h)*r/Z)
    if r<=NC:return normal
    if r>=NR:return vertical
    t=(r-NC)/(NR-NC);b=t**3*(10-15*t+6*t*t)
    return tuple((1-b)*x+b*y for x,y in zip(normal,vertical,strict=True))


def certificate():
    iv=mp.iv;mp.mp.dps=40;leaves=[]
    def interval(lo,hi):return iv.mpf([lo,hi])
    def lower(x):return float(x.a)
    def upper(x):return float(x.b)
    def bounds(q0,q1,r0,r1):
        q=interval(q0,q1);r=interval(r0,r1);pi=iv.pi
        h=(1+iv.cos(pi*q/L))/2;p=-pi/(2*L)*iv.sin(pi*q/L);pp=-pi*pi/(2*L*L)*iv.cos(pi*q/L)
        a=iv.sqrt(1+p*p);J=1-r*pp/(a*a*a);v=1-h/Z
        lo=(r0-NC)/(NR-NC);hi=(r1-NC)/(NR-NC)
        smooth=lambda t:t**3*(10-15*t+6*t*t)
        rate=lambda t:30*t*t*(1-t)**2/(NR-NC)
        b=interval(math.nextafter(smooth(lo),-math.inf),math.nextafter(smooth(hi),math.inf))
        bp=interval(math.nextafter(min(rate(lo),rate(hi)),-math.inf),math.nextafter(1.875/(NR-NC) if lo<=.5<=hi else max(rate(lo),rate(hi)),math.inf))
        xq=(1-b)*J+b;zq=p*((1-b)*J+b*v)
        xr=p/a*(r*bp-(1-b));zr=(1-b)/a+b*v+r*bp*(v-1/a)
        return lower(xq*zr-zq*xr),lower(xq),lower(J)
    def box(q0,q1,r0,r1,depth=0):
        det,xq,j=bounds(q0,q1,r0,r1)
        if det>0 and xq>0:leaves.append((det,xq,j));return
        require(depth<18,'map interval certificate unresolved; no geometry admission')
        if (q1-q0)/.1>(r1-r0)/.05:
            mid=(q0+q1)/2;box(q0,mid,r0,r1,depth+1);box(mid,q1,r0,r1,depth+1)
        else:
            mid=(r0+r1)/2;box(q0,q1,r0,mid,depth+1);box(q0,q1,mid,r1,depth+1)
    for i in range(50):
        for j in range(14):box(-L+2*L*i/50,-L+2*L*(i+1)/50,NC+(NR-NC)*j/14,NC+(NR-NC)*(j+1)/14)
    return dict(policy='C1CollarRoofV1',normal_collar=NC,transition_end=NR,roof=Z,
        interval_transition_cells=len(leaves),minimum_interval_Jacobian=min(v[0] for v in leaves),minimum_interval_xq=min(v[1] for v in leaves),
        collar_sufficient_Jacobian_bound=1-NC*math.pi**2/(2*L**2),outer_vertical_Jacobian_bound=1-1/Z,
        pure_normal_at_roof_foot_limit=1-Z*math.pi**2/(2*L**2),
        injectivity='positive oriented one-sided analytic charts, x_q>0 and simple fixed inlet/roof/outlet/ground boundary; matching vertical C1-foot branches; degree-one global inverse on each conforming subdomain',
        scope='analytic map interval certificate, not a certificate for the native polynomial geometry interpolant or coupled curved evolution')


def curved_geometry(root,radial,amplitude):
    root.remove(root.find('GEOMETRY'));mesh(root,X,radial);g=root.find('GEOMETRY')
    seeds={int(v.get('ID')):tuple(map(float,v.text.split())) for v in g.find('VERTEX')}
    def point(seed):
        x,z=mapping(seed[0],seed[1],amplitude);return x,z,0.
    for v in g.find('VERTEX'):v.text=' '.join(f'{x:.17g}' for x in point(seeds[int(v.get('ID'))]))
    if amplitude:
        curves=element(g,'CURVED')
        for i,edge in enumerate(g.find('EDGE')):
            left,right=[seeds[int(k)] for k in edge.text.split()];values=[]
            for j in range(16):
                t=j/15;seed=tuple((1-t)*a+t*b for a,b in zip(left,right,strict=True));values+=point(seed)
            element(curves,'E',' '.join(f'{v:.17g}' for v in values),ID=i,EDGEID=int(edge.get('ID')),NUMPOINTS=16,TYPE='PolyEvenlySpaced')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','reference','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();manifest=RequestManifest(a.output,dict(schema='CT7TerrainFlatPreflightV1',runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__)),reference=sha(a.reference)))
    radial=[0.,H/64,H/32,H/16,H/8,H/4,H/2,H,NC,NR,Z];r=H
    while r*1.25<Z:r*=1.25;radial.append(r)
    radial=sorted(set(radial));cert=certificate();cases={}
    refs=json.loads(a.reference.read_text())
    refdef=dict(schema='CT7ProspectivePhysicalReferenceV1',fit_source_sha256=sha(a.reference),fits=refs['fits'],
        wall_origin=0.,displacement=.0875,reference_height=H,normalization='Hh=.04m, original U(Hh); no cap',
        outer='U=utau/kappa*log((z-d)/z0),k=q*utau^2,epsilon=utau^3/(kappa*(z-d)); roof20Hh',
        inner='requires admitted separately versioned displaced Accepted policy; unavailable for failed original roughness matches',
        boundary='AcceptedRichardsNorrisFlatV1 physical data; proposed displaced extension not yet native qualified',physical_admission=False)
    write_json(a.output/'reference-definition.json',refdef)
    for degree in (3,4):
        for kind,amp in [('flat',0.),('terrain',1.)]:
            folder=a.output/(kind+'-p'+str(degree));xml=case_xml(folder,'A95','Natural','Accepted',steps=0,descriptor_probe=a.runtime/'RANSNewModelsProbe',degree=degree,radial_knots=radial,xknots=X)
            root=ET.parse(xml).getroot();curved_geometry(root,radial,amp)
            set_parameter(root,'RANSWallHeight',H);set_parameter(root,'RANSWallSampleHeight',H);set_parameter(root,'RANSWallTailStart',H)
            set_parameter(root,'RANSOpenDisplacement',.0875);set_info(root,'RANSWallBoundaryPolicy','AcceptedRichardsNorrisFlatV1');set_info(root,'RANSOpenTop','RichardsNorris2011V1')
            for f in list(root.findall('CONDITIONS/FUNCTION')):root.find('CONDITIONS').remove(f)
            write(root,xml);row=dict(xml_sha256=sha(xml),degree=degree,pressure_degree=degree-2,geometry=kind,reference_definition_sha256=sha(a.output/'reference-definition.json'))
            if kind=='flat':
                packet,e=probe(a.runtime,folder,dict(operation='native_layout_inventory'),'count')
                require(e['exit_code']==0,'actual native flat count failed');row.update(native_count=packet,execution=e,count_packet_sha256=sha(folder/'count.json'))
            geometry={}
            for order in (16,32):
                dest=folder/('geometry-q'+str(order));dest.mkdir();g=copy.deepcopy(root)
                for exp in g.findall('EXPANSIONS/E'):exp.set('NUMPOINTS',str(order)+','+str(order))
                write(g,dest/'case.xml');packet,e=probe(a.runtime,dest,dict(operation='mapped_mesh_geometry_inventory'),'geometry')
                require(e['exit_code']==0,'native mesh Jacobian failed:'+kind)
                expected=1100. if not amp else 1097.5
                geometry[str(order)]=dict(packet,execution=e,expected_physical_volume=expected,volume_defect=abs(packet['physical_volume']-expected),packet_sha256=sha(dest/'geometry.json'))
            row['native_geometry']=geometry;cases[kind+'-p'+str(degree)]=row
            write_json(a.output/'progress.json',dict(certificate=cert,cases=cases));print('CT7 target preflight',kind,degree,'complete',flush=True)
    require(all(not v['native_count']['within_unchanged_3000_guard'] for k,v in cases.items() if k.startswith('flat')),'unexpected capacity: amend prelaunch contract before any evolution')
    result=dict(schema='CT7TerrainFlatPreflightV1',status='pass',definitions=dict(x=X,radial=radial,roof=Z,H=H,collar=NC,transition_end=NR,reference=sha(a.output/'reference-definition.json')),
        analytic_geometry_certificate=cert,cases=cases,full_domain_preservation='blocked_by_native_capacity_and_complete_physical_reference',native_curved_coupling='unsupported; analytic/native geometry preflight is not coupled-operator admission')
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json','reference-definition.json'])
    return 0

if __name__=='__main__':raise SystemExit(main())
