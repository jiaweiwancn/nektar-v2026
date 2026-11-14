#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only CT8 observation, retained temporal vectors and capacity controls."""
import argparse,copy,csv,json,math,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import write_json,require,parameters
from ct_accepted5_geometry import closest_normal
from ct_accepted7_domain import X,Z,NC,NR,certificate,curved_geometry
from ct_accepted_open import probe,row_value
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_parameter,write
from ct_accepted6_compare import observations
from rans_request_manifest import RequestManifest,sha,runtime_identity

def read(p):return json.loads(Path(p).read_text())

def observation_audit(data,case):
    H=case['H'];records=list(csv.DictReader((data/'profiles.csv').open()));out=[]
    # Reuse the sealed V2 complete closest-point audit, including envelopes
    # and all Fig7 ratios; frozen source roles and uncertainty remain unchanged.
    audit_path=Path(__file__).resolve().parents[4]/'docs/plans/cao-tamura-accepted-5-evidence/qualified-v1/observation-audit.json'
    # Repository evidence index records the actual export name; caller may
    # locate it without inventing or altering observation coordinates.
    candidates=list((Path(__file__).resolve().parents[4]/'docs/plans/cao-tamura-accepted-5-evidence').rglob('*observation*.json'))
    base=None
    for path in candidates:
        q=read(path)
        if isinstance(q,dict) and q.get('schema')=='CT5ObservationAuditV1':base=q;audit_path=path;break
    require(base is not None,'missing corrected CT5 observation audit')
    for q in base['height_audit']:
        r=copy.deepcopy(q);d=r['nominal'].get('distance');vals=[]
        def scan(x):
            if isinstance(x,dict):
                if x.get('distance') is not None:vals.append(float(x['distance']))
                for k,v in x.items():
                    if k!='nominal':scan(v)
            elif isinstance(x,list):
                for v in x:scan(v)
        scan(r)
        if d is not None:vals.append(d)
        r['CT8_nominal']='unavailable' if d is None else 'inner_AWF_coupled_comparison' if d<H else 'outer'
        r['CT8_extraction_range']=dict(minimum=min(vals),maximum=max(vals)) if vals else None
        r['CT8_interface_overlapping']=bool(vals and min(vals)<=H<=max(vals));out.append(r)
    rr=[r for r in out if r['series'].startswith('f12_RR')];over=[r for r in out if r['series'].startswith('f5b_RR') and r['CT8_interface_overlapping']]
    require(rr and all(r['CT8_nominal']=='inner_AWF_coupled_comparison' for r in rr),'7mm RR observations not classified inner')
    require(any(abs(abs(float(r['x_over_Hh']))-1.25)<1e-10 for r in over),'missing RR +/-1.25 interface overlap')
    p=parameters('A95');us=case['ustar'];nu=case['nu'];d=case['displacement']
    return dict(schema='CT8ObservationAuditV1',status='pass',source=sha(audit_path),height_A_mm=7.5,rows=out,
        interface_overlapping=[dict(series=q['series'],point_id=q['point_id'],x=q.get('x_over_Hh'),range=q['CT8_extraction_range']) for q in over],
        indicators=dict(H_plus=H*us/nu,zeta_plus=(H-d)*us/nu,outer_nu_over_nut=nu/(p['kappa']*us*(H-d)),
            cube_clearance_mm=2.5,curvature_height_product=math.pi**2/12.5*H,metric_minimum_stretch=1-math.pi**2/12.5*H,
            historical_point1_screen_pass=False),
        Fig12='7mm RR observations are inside the AWF and require actual transported-state stresses, not outer equilibrium stress',
        Fig11='physical heights unavailable where delta7.5 unspecified',scope='closest-normal graph and frozen extraction envelopes; no invented physical probe uncertainty; positive metric is not curved accuracy')

def temporal(old,ct7):
    members=[old/'evolution-v2/A95-p4-dt002',old/'evolution-v3/A95-p4-dt001',ct7/'ct6-finer-p4']
    audits=[old/'physical-audits-v3/A95-p4-dt002',old/'physical-audits-v3/A95-p4-dt001',ct7/'ct6-finer-p4-audit']
    records=[observations(x,y) for x,y in zip(members,audits,strict=True)]
    reports=[read(x/'trajectory/case.evolution.json') for x in members]
    require(len({len(x['state']) for x in reports})==1 and all(r['initial']==reports[0]['initial'] for r in reports),'temporal coefficient spaces/initials differ')
    # Native volume rows are the actual physical mass quadrature, independent
    # from the common observation RMS used by the sealed temporal classifier.
    ini=Path(read(members[-1]/'results.json')['initial_packet']['path']);require(sha(ini)==read(members[-1]/'results.json')['initial_packet']['sha256'],'changed native mass-row packet')
    packet=read(ini);o=packet['observation'];off=o['offsets'];a=np.asarray(reports[0]['state'])-reports[1]['state'];b=np.asarray(reports[1]['state'])-reports[2]['state'];out={}
    for block,name in enumerate(('U','w','k','outer_epsilon','Pi')):
        aa=bb=ab=rr=0.
        for q in o['volume']:
            if block==3 and not q['region']:continue
            row=q['epsilon'] if block==3 else q['pressure'] if block==4 else q['velocity_k'];ob=off[4+q['region']] if block==4 else off[block]
            av=row_value(row,a,ob)[0];bv=row_value(row,b,ob)[0];weight=q['weight']
            aa+=weight*av*av;bb+=weight*bv*bv;ab+=weight*av*bv;rr+=weight*(av-2*bv)**2
        out[name]=dict(mass_norm_a=math.sqrt(aa),mass_norm_b=math.sqrt(bb),alignment=ab/math.sqrt(aa*bb),
            a_minus_2b_mass_norm=math.sqrt(rr),relative_to_a=math.sqrt(rr/aa),ratio=math.sqrt(aa/bb),
            physical_mass_definition='sum actual native volume weights times represented difference squared; pressure includes both regional blocks, epsilon only outer')
    absolute={}
    for dt,r in zip((.002,.001,.0005),records,strict=True):
        delta=r[3]-r[4];w=r[5];regions=r[6];rows={}
        for j,name in enumerate(('U','w','k','outer_epsilon','Pi')):
            mask=regions==1 if j==3 else np.ones(len(w),bool);ww=w[mask];e=delta[mask,j];mean=float(ww@e/ww.sum())
            rows[name]=dict(absolute_RMS=float(np.sqrt(ww@(e*e)/ww.sum())),mean_error=mean,variation_RMS=float(np.sqrt(ww@((e-mean)**2)/ww.sum())))
        absolute[str(dt)]=rows
    return dict(schema='CT8RetainedTemporalVectorDiagnosisV1',execution='pass',differences=out,absolute_manufactured_errors=absolute,
        binding=[dict(member=str(m),receipt=sha(m/'results.json'),audit=sha(q/'results.json'),native_report=sha(m/'trajectory/case.evolution.json')) for m,q in zip(members,audits,strict=True)],
        historical_gate='CT6/CT7 failed upper-band U/w classifications unchanged; no fresh governing trajectory or prospective threshold selected from these results',
        interpretation='aligned decreasing differences and a-2b diagnose finite-step contamination; do not prove asymptotic order or remove the remaining joined-map risk')

def domain(rt,dest,case):
    dest.mkdir();H=case['H'];radial=[0.,H/64,H/32,H/16,H/8,H/4,H/2,H,NC,NR,Z];r=H
    while r*1.25<Z:r*=1.25;radial.append(r)
    radial=sorted(set(radial));out=dict(schema='CT8TargetDomainInventoryV1',radial_grid=radial,x=X,height=H,roof=Z,geometry_certificate=certificate(),cases={})
    for degree in (3,4):
        for kind,amp in [('flat',0.),('terrain',1.)]:
            folder=dest/(kind+'-p'+str(degree));xml=case_xml(folder,'A95','Natural','Accepted',steps=0,descriptor_probe=rt/'RANSNewModelsProbe',degree=degree,radial_knots=radial,xknots=X)
            root=ET.parse(xml).getroot();curved_geometry(root,radial,amp)
            pars=root.find('CONDITIONS/PARAMETERS')
            for n in list(pars):pars.remove(n)
            for k,v in dict(NumSteps=0,TimeStep=.001,Time=0.,RANSWallHeight=H,RANSWallSampleHeight=H,RANSWallTailStart=H).items():set_parameter(root,k,v)
            for f in list(root.findall('CONDITIONS/FUNCTION')):root.find('CONDITIONS').remove(f)
            write(root,xml);row=dict(mesh_sha256=sha(xml),degree=degree,physical_reference=case,geometry_only=True)
            if kind=='flat':
                q,e=probe(rt,folder,dict(operation='native_layout_inventory'),'count');require(e['exit_code']==0,'native capacity count failed');row['native_count']=q
            row['native_geometry']={}
            for order in (16,32):
                sub=folder/('geometry-'+str(order));sub.mkdir();r=copy.deepcopy(root)
                for ex in r.findall('EXPANSIONS/E'):ex.set('NUMPOINTS',str(order)+','+str(order))
                write(r,sub/'case.xml');q,e=probe(rt,sub,dict(operation='mapped_mesh_geometry_inventory'),'geometry');require(e['exit_code']==0,'geometry preflight failed');row['native_geometry'][str(order)]=dict(q,execution=e,packet_sha256=sha(sub/'geometry.json'))
            out['cases'][kind+'-p'+str(degree)]=row
    out.update(status='preflight_complete',full_domain_evolution='blocked by unchanged3000slot guard',joined_map_coupled_evolution='not qualified');return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','reference','old','ct7','data','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();m=RequestManifest(a.output,dict(schema='CT8IndependentAssessmentV1',source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),reference=sha(a.reference)))
    c=read(a.reference)['case'];results=dict(observations=observation_audit(a.data,c),temporal=temporal(a.old,a.ct7),domain=domain(a.runtime,a.output/'domain',c))
    for name,r in results.items():write_json(a.output/(name+'.json'),r)
    write_json(a.output/'results.json',dict(status='pass',scope='read-only diagnostics and geometry/count inventory only',files={name:sha(a.output/(name+'.json')) for name in results}));m.seal(a.output,['results.json',*(name+'.json' for name in results)]);return 0
if __name__=='__main__':raise SystemExit(main())
