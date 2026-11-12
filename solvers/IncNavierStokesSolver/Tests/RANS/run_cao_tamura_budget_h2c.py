#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded analytical-neutral H2c family. Historical finite cases are not restarts."""
import argparse, csv, hashlib, json, math, os, platform, subprocess, time
from pathlib import Path
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic, freeze, CampaignLock
from run_cao_tamura_budget_h2a import prepare
from cls_epawf_m4_cases import tensor_mesh
from diagnose_cls_accepted_covariance import run_probe, checkpoint, common_packet
from run_cls_epawf_picard import THREADS,load,require,sha,set_info,set_parameter,write,state_hash
from cao_tamura_budget_h1 import DATA, NU, USTAR, UINF

X=[-.2,-.1,-.05,0.,.05,.1,.6]

def case(folder,h1,preflight,H,geometry='hill',steps=0,dt=2e-6,restart=None,audit=True,fixture=False,seed="TangentialNeutralV1"):
    h=next(v for v in preflight['heights'] if abs(v['H']-H)<1e-10)
    xml=prepare(folder,h1,steps,dt,restart,not fixture)
    t=ET.parse(xml).getroot()
    # Negative viscosity intercept is permitted. Knots follow the positive
    # molecular/eddy and cutoff scales, not the previous finite-family knot.
    rapid=preflight['nu']/(preflight['cMu_wall']*preflight['cLength_wall']*math.sqrt(h['K0']))-h['viscosity_intercept']
    require(0 < rapid/8 < .00005,'new physical rapid-variation knot outside first panel')
    z=[0.,rapid/8,.00005,.0002,h['cutoff'],h['h_s'],H,.012,.025,.05,.1,.2,.4,.8]
    tensor_mesh(t,z,3,2,.8,1.,points=[16,16],xknots=X)
    for k,v in dict(RANSWallReference='AnalyticalNeutralInflowV1',RANSWallRegularizedReference='AnalyticalNeutralInflowV1',
        RANSWallCompositeSeed=seed,RANSWallNeutralFitSHA256=hashlib.sha256(json.dumps(preflight['fit'],sort_keys=True).encode()).hexdigest(),
        RANSWallGeometry='LegacyV1' if geometry=='affine' else 'CompositeHillV1',
        RANSWallMap='Flat' if geometry=='affine' else 'CosineHill',RANSWallPimpleKInventory='AcceptedBudgetV2',
        RANSWallAuditReuse='ExactObservationV1' if audit else 'Off').items():set_info(t,k,v)
    # Remove unused finite-content identity entirely from the new family.
    info=t.find('./CONDITIONS/SOLVERINFO')
    for e in list(info):
        if e.get('PROPERTY') in ('RANSWallReferenceFile','RANSWallReferenceSHA256'):info.remove(e)
    for k,v in dict(RANSWallHeight=H,RANSWallTailStart=H,RANSWallSampleHeight=H,
        RANSWallRoughness=h['h_s'],RANSWallNeutralZ0=preflight['fit'][0]['z0_eff'],
        RANSWallMapAmplitude=.04 if geometry=='hill' else 0.,RANSWallMapLengthX=.1,
        Kinvis=preflight['nu'],RANSWallFrictionVelocity=preflight['u_tau'],
        K0=h['K0'],RANSWallNormalOrder=16,RANSWallTangentialOrder=16).items():set_parameter(t,k,v)
    write(t,xml)
    ni=sum(v<=H for v in z[1:]);no=len(z)-1-ni;nx=len(X)-1
    sizes=[(3*nx+1)*(3*(ni+no)+1)]*3+[(3*nx+1)*(3*no+1),(nx+1)*(ni+1),(nx+1)*(no+1),1]
    require(sum(sizes)<=3000,'H2c dense size guard exceeded')
    definition=dict(version='CaoTamuraBudgetH2cV1',H_A=H,h_s=h['h_s'],z_m=H,z_K=H,
        inlet='AnalyticalNeutralInflowV1 global unshifted height; no constant tail',fit_sha256=hashlib.sha256(json.dumps(preflight['fit'],sort_keys=True).encode()).hexdigest(),
        geometry=geometry,boundary='CompatibilityV1' if fixture else 'RoughFlatOpenV1',
        mesh=dict(x=X,z=z,degree=3,regional_pressure_degree=1,points=[16,16],first_knot=rapid/8,rapid_viscosity_scale=rapid),slots=sum(sizes),sizes=sizes,
        model='Regularized CLS beta=.4 Gamma=1 gamma0=0; Linear inner; calibrated budget AWF',
        inventory='AcceptedBudgetV2',dt=dt,physical_steps=steps,fresh_history=restart is None,
        initialization=seed+'; native analytical/AWF radial then physical streamfunction and native continuity setup',
        measurement_qualification='exploratory: matching plane only 1.1-1.3 cube heights; no established roughness-sublayer clearance',
        checkpoint=str(restart) if restart else None,xml_sha256=sha(xml))
    atomic(folder/'h2c-definition.json',definition);atomic(folder/'definition.json',definition)
    return xml

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('solver','probe','preflight','h1','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--phase',choices=['screen','primary','flat','flat-finalize','native-seed','hill','cache'],default='screen')
    p.add_argument('--gate',type=Path)
    p.add_argument('--flat-root',type=Path)
    p.add_argument('--seed',choices=['TangentialNeutralV1','PhysicalHeightInterfaceV1'],default='TangentialNeutralV1')
    a=p.parse_args();require(platform.system()=='Linux','authoritative Linux required')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output)
    os.environ.update(THREADS);os.environ.update(VECLIB_MAXIMUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
    solver,image=freeze(a.solver.resolve(),a.output/'runtime');probe,pimage=freeze(a.probe.resolve(),a.output/'probe-runtime')
    preflight=load(a.preflight);out=dict(status='started',runtime=image,probe=pimage,preflight_sha256=sha(a.preflight),results={},rank=1,threads=THREADS)
    def save(**kw):out.update(kw);atomic(a.output/'progress.json',out)
    def native(name,H,geometry='affine',steps=0,dt=2e-6,cp=None,fixture=False,mutate=None):
        xml=case(a.output/name,a.h1,preflight,H,geometry,steps,dt,cp,fixture=fixture,seed=a.seed)
        if mutate:
            t=ET.parse(xml).getroot();mutate(t);write(t,xml)
        env=dict(os.environ,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(xml.parent),NEKTAR_RANS_TIMING_FILE=str(xml.parent/'timing.jsonl'))
        cmd=['/usr/bin/time','-f','%e %U %S %M','-o',str(xml.parent/'usage.txt'),str(solver),str(xml)]
        start=time.monotonic()
        with (xml.parent/'solver.log').open('w') as f:code=subprocess.run(cmd,cwd=xml.parent,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=1800).returncode
        report=xml.with_suffix('.evolution.json');q=load(report) if report.exists() else None
        receipt=dict(code=code,seconds=time.monotonic()-start,usage=(xml.parent/'usage.txt').read_text(),command=cmd,
            xml_sha256=sha(xml),report=str(report),report_sha256=sha(report) if report.exists() else None,
            failure=q.get('failure') if q else (xml.parent/'solver.log').read_text()[-3000:])
        if q:
            for row in q.get('pimple_attempts',[]):
                if row.get('accepted',not row.get('rejected',False)):
                    for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):require(row[k]==1,'wrong schedule '+k)
                    require(row['solved_epsilon_retained_bitwise'],'epsilon overwritten')
        out['results'][name]=dict(receipt=receipt,report=str(report));save();return q
    save()
    try:
        if a.phase in ('screen','primary','cache'):
            heights=[.0065] if a.phase in ('cache','primary') else [.0055,.006,.0065]
            for H in heights:
                for geom in ('affine','hill') if a.phase=='screen' else ('hill',):
                    name=str(round(H*10000))+'-'+geom
                    xml=case(a.output/name,a.h1,preflight,H,geom,audit=a.phase!='cache',seed=a.seed)
                    d,r=run_probe(probe,xml,dict(operation='h2c_seed',inner_stress_model='Linear',xml_sha256=sha(xml)),a.output/(name+'-probe'))
                    out['results'][name]=dict(receipt=r,state_sha256=d['state_sha256'],audit=d['audit'],setup=d['seed_setup'],checks=d['checks'],
                        moving_foot_fixture=d['moving_foot_fixture'],scalar_inventory=d['scalar_inventory'],seed_constitutive=d['seed_constitutive'],target=d['seed_target_witnesses'],budget=d['native_budget'])
                    save()
            save(status='screen_completed')
        elif a.phase in ('flat','flat-finalize'):
            from compare_cls_epawf_pimple import same_history
            source=a.output
            if a.phase=='flat-finalize':
                require(a.flat_root,'immutable completed flat root required')
                previous=load(a.flat_root/'progress.json')
                require(previous['status']=='gate_stop' and 'finite ground reaction requires' in previous['failure'] and previous['tests']==dict(split_restart_bitwise=True,post_assignment_rollback=True),
                    'only the saved readback-only common-lifting stop may be finalized')
                source=a.flat_root
                q=load(source/'flat-h0/case.evolution.json');r=load(source/'flat8/case.evolution.json')
                first=load(source/'flat4/case.evolution.json');second=load(source/'flat4-restart/case.evolution.json')
                for label in ('flat-h0','flat8','flat4','flat4-restart'):
                    old=previous['results'][label]['receipt'];require(old['code']==0 and sha(Path(old['report']))==old['report_sha256'],'retained accepted report identity changed')
                require(not q['failure'] and not r['failure'],'retained trajectory is not accepted')
                out['retained_accepted_trajectory']=str(source);out['new_physical_steps']=0
                out['results'].update({k:v for k,v in previous['results'].items() if k in ('flat-h0','flat8','flat4','flat4-restart','flat-assigned-reject')})
                cp=Path(q['checkpoint'])
            else:
                q=native('flat-h0',.0065,fixture=True);require(q and not q['failure'],'flat initialization failed')
                cp=Path(q['checkpoint']);r=native('flat8',.0065,steps=8,cp=cp,fixture=True)
                require(r and not r['failure'],'flat compatibility trajectory rejected')
                first=native('flat4',.0065,steps=4,cp=cp,fixture=True)
                second=native('flat4-restart',.0065,steps=4,cp=Path(first['checkpoint']),fixture=True)
            require(same_history(r,second),'analytical flat restart history/state mismatch')
            parent_sha=sha(cp/'State.bin')
            bad=load(source/'flat-assigned-reject/case.evolution.json') if a.phase=='flat-finalize' else native('flat-assigned-reject',.0065,steps=1,cp=cp,fixture=True,mutate=lambda t:set_parameter(t,'RANSWallPimplePostAssignRejectStep',1))
            require(same_history(q,bad) and sha(cp/'State.bin')==parent_sha,'analytical rollback changed parent/state/history')
            out['tests']=dict(split_restart_bitwise=True,post_assignment_rollback=True)
            for label,report in [('flat-h0',q),('flat8',r)]:
                parent=Path(report['checkpoint']);seed=checkpoint(parent)
                packet=dict(common_packet(parent),operation='h2c_checkpoint',inner_stress_model='Linear',
                    state=seed['coupled'],history=seed['coupled'],alpha=0.,time=seed['time'],
                    native_report_file=str(source/label/'case.evolution.json'),native_report_sha256=sha(source/label/'case.evolution.json'))
                d,receipt=run_probe(probe,source/label/'case.xml',packet,a.output/(label+'-readback'))
                out['results'][label]['readback']=dict(receipt=receipt,diagnostic=d)
                save()
            out['flat_reference_preservation'] = dict(initial=q['finite_inflow_reference'],final=r['finite_inflow_reference'],
                physical_samples_initial=q['physical_samples'],physical_samples_final=r['physical_samples'],
                result='absolute predeclared targets assessed in collector; no response-relative drift screen')
            for key,mutation in [('fit',lambda t:set_info(t,'RANSWallNeutralFitSHA256','0'*64)),
                                 ('origin',lambda t:set_parameter(t,'RANSWallReferenceOriginZ',-1e-5)),
                                 ('inventory',lambda t:set_info(t,'RANSWallPimpleKInventory','AcceptedBudgetV1'))]:
                bad=native('reject-'+key,.0065,cp=cp,fixture=True,mutate=mutation)
                require(not bad or bad['failure'],'incompatible case accepted')
            save(status='completed_flat8')
        elif a.phase=='native-seed':
            q=native('hill-h0',.0065,'hill')
            save(status='initialization_rejected' if not q or q['failure'] else 'initialized_zero_steps')
        else:
            require(a.gate,'screen gate required');g=load(a.gate/'progress.json')
            require(g['results']['65-hill']['audit']['accepted_awf_hard_admissible']==[1],'hill seed gate failed; no evolving retry')
            # Implemented only once the native gate has passed. Each dt uses
            # its own fresh initialization and integrator histories.
            zero=native('hill-h0',.0065,'hill');finezero=native('hill-half-h0',.0065,'hill',dt=1e-6)
            require(state_hash(zero['state'])==state_hash(finezero['state']),'dt changed mapped seed')
            c=native('hill8',.0065,'hill',8,cp=Path(zero['checkpoint']));f=native('hill16',.0065,'hill',16,dt=1e-6,cp=Path(finezero['checkpoint']))
            require(c and f and not c['failure'] and not f['failure'],'hill leg rejected')
            from compare_cls_epawf_pimple import same_history
            first=native('hill4',.0065,'hill',4,cp=Path(zero['checkpoint']));second=native('hill4restart',.0065,'hill',4,cp=Path(first['checkpoint']))
            require(same_history(c,second),'mapped split restart differs')
            bad=native('hill-rollback',.0065,'hill',1,cp=Path(zero['checkpoint']),mutate=lambda t:set_parameter(t,'RANSWallPimplePostAssignRejectStep',1))
            require(same_history(zero,bad),'mapped rollback differs')
            save(status='completed_hill_pair')
    except Exception as e:save(status='gate_stop',failure=str(e));raise
    finally:lock.stream.close()

if __name__=='__main__':main()
