#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Actual-target finite-step pressure reconstruction, never a governing update.

Uses captured CSR blocks, sparse consistent mass and native aggregate loads.
No dense Schur, M^-1G columns or per-element global vectors are constructed.
"""
import argparse,gc,inspect,json,os,resource,shutil,time,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu,LinearOperator,gmres
import ct_accepted8_native as n
from ct_accepted10_diagnose import matrix,PressureImages
from ct_accepted9_pressure_run import sealed
from ct_accepted_pilot import replay
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity
from run_cls_epawf_picard import set_parameter,write
PROGRESS=None


def maximum(v):
    return float(np.max(np.abs(v),initial=0.))


def smx(m):
    return maximum(m.data)


def factor(m):
    lu=splu(m.tocsc())
    require(lu.L.nnz+lu.U.nnz<=40000000,'diagnostic sparse factor exceeds40M')
    return lu


class Frozen:
    def __init__(self,s,initial,c,weights,count):
        self.s=s;self.nv=s['velocity_size'];self.rows=np.asarray(s['rows'],int)
        self.images=PressureImages(s,initial,c,weights,count)
        self.M=self.images.M;self.ml=self.images.lu
        self.G=self.images.G;self.D=self.images.D;self.C=self.images.C
        self.F=matrix(s['frozen_flow_matrix']);nv=self.nv
        self.A=self.F[:nv,:nv].tocsr();self.al=factor(self.A)
        # Match the native CSR extended accumulation and three-step direct
        # residual refinement. This matters for the tiny continuity RHS, whose
        # relative accuracy cannot be inferred from a large-RHS total solve.
        self.extended={key:m.astype(np.longdouble) for key,m in
                       (('M',self.M),('G',self.G),('D',self.D),('C',self.C),('A',self.A))}
        self.alpha=s['alpha'];self.rc=np.asarray(s['frozen_flow_residual'][nv:])
        self.rm=np.asarray(s['frozen_flow_residual'][:nv]);self.receipts=[]
        self.S=LinearOperator(self.C.shape,matvec=lambda p:self.mul('D',self.direct('M',self.mul('G',p)))/self.alpha-self.mul('C',p))
        pc=(self.D@diags(1/self.M.diagonal())@self.G/self.alpha-self.C).tocsc()
        self.pl=factor(pc);self.P=LinearOperator(self.S.shape,matvec=self.pl.solve)
        self.blocks=dict(G=smx(self.F[:nv,nv:]-self.G),D=smx(self.F[nv:,:nv]-self.D),C=smx(self.F[nv:,nv:]-self.C))
        require(max(self.blocks.values())<=2e-10 and not np.any(self.C.data),'captured full blocks/C disagree')
        require(len(initial['observation']['state'])-1 not in self.rows,'unused gauge entered solve')

    def mul(self,key,v):
        return np.asarray(self.extended[key]@np.asarray(v,dtype=np.longdouble),dtype=np.float64)

    def direct(self,key,rhs):
        lu=self.ml if key=='M' else self.al;rhs=np.asarray(rhs);x=lu.solve(rhs)
        bound=2e-12*maximum(rhs)+1e-25
        for attempt in range(4):
            residual=self.mul(key,x)-rhs
            if maximum(residual)<=bound:return x
            require(attempt<3,'diagnostic direct solve true residual:'+key)
            x-=lu.solve(residual)
        raise ValueError('unreachable diagnostic refinement')

    def solve(self,rhs,label,rtol=2e-12,operator=None):
        rhs=np.asarray(rhs);mag=float(np.linalg.norm(rhs))
        if mag==0:return np.zeros(self.S.shape[0])
        op=self.S if operator is None else operator;ticks=[];started=time.process_time()
        absolute=1e-27 if rtol<2e-12 else 1e-25
        if PROGRESS:write_json(PROGRESS,dict(label=label,rhs_norm=mag,rhs_maximum=maximum(rhs),relative_tolerance=rtol,absolute_tolerance=absolute,phase='begin'))
        print('CT12 diagnostic solve',label,'rhs',mag,flush=True)
        kw={'rtol' if 'rtol' in inspect.signature(gmres).parameters else 'tol':rtol}
        # scipy counts restart cycles; 12x80 stays below the native component's
        # 1000-iteration work ceiling. Tiny continuity noise also needs the
        # explicit UNSCALED absolute tolerance, not relative-only normalization.
        p,info=gmres(op,rhs/mag,M=self.P,atol=absolute/mag,restart=80,maxiter=12,
                     callback=lambda x:ticks.append(float(x)),callback_type='pr_norm',**kw)
        p*=mag;res=op@p-rhs;relative=float(np.linalg.norm(res))/mag
        receipt=dict(label=label,iterations=len(ticks),relative_residual=relative,maximum_true_residual=maximum(res),
            absolute_tolerance=absolute,rhs_maximum=maximum(rhs),CPU=time.process_time()-started,info=int(info),
            requested_scipy_stopping_test_met=bool(info==0),physical_image=self.measure(p))
        self.receipts.append(receipt)
        if PROGRESS:write_json(PROGRESS,dict(receipt,phase='end'))
        # A restart-limit label is not an error bound. Retain it explicitly,
        # require the independently recomputed unscaled residual, and then
        # qualify the physical image by the mandatory tighter/signed-sum tests.
        require(info>=0 and maximum(res)<=2e-10*maximum(rhs)+absolute,'true pressure residual:'+label+':'+str(receipt))
        return p

    def image(self,load,label,rtol=2e-12):
        v=np.asarray(load)[self.rows[:self.nv]]
        return self.solve(-self.mul('D',self.direct('A',v)),label,rtol)

    def mass_image(self,load,label,rtol=2e-12):
        v=np.asarray(load)[self.rows[:self.nv]]
        return self.solve(-self.mul('D',self.direct('M',v))/self.alpha,label,rtol)

    def measure(self,p):return self.images.measure(p)
    def inner(self,a,b):return self.images.inner(a,b)


def main():
    global PROGRESS
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct10','ct11','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--reuse-packets',type=Path,help='Exact immutable native packets from a failed diagnostic; never reuse derived classifications')
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    target=sealed(a.ct11/'target-v1');pre=sealed(a.ct10/'target-preflight-v1')
    require(target['qualification']['reasons']==['temporal:Pi'],'historical pair changed')
    initial_path=a.ct10/'target-preflight-v1/target/readiness.json'
    require(sha(initial_path)==pre['packet_sha256'],'changed target observation/test packet')
    initial=n.read(initial_path);c=target['case'];require(c==pre['mesh'],'foreign target mesh')
    weights=n.inventory(c)[2];count=len(weights);off=initial['observation']['offsets']
    init=np.asarray(initial['observation']['state']);require(len(init)==60907,'not actual target')
    require(sha(a.runtime/'RANSPreservationProbe')=='20236b14eb4a90dea0f2a9a38b92e41135a83e9241c978fb265e185098f586e8','foreign read-only runtime')
    reused={}
    if a.reuse_packets:
        old=n.read(a.reuse_packets/'request-manifest.json')['request']
        require(old['target']==sha(a.ct11/'target-v1/results.json') and old['preflight']==sha(a.ct10/'target-preflight-v1/results.json') and
                old['initial']==sha(initial_path) and old['runtime']==runtime_identity(a.runtime/'RANSPreservationProbe'),'foreign native packet request')
        for folder in sorted(a.reuse_packets.iterdir()):
            if folder.is_dir() and (folder/'loads.json').exists():
                reused[folder.name]={name:sha(folder/name) for name in ('case.xml','loads-request.json','loads-execution.json','loads.json')}
    manifest=RequestManifest(a.output,dict(schema='CT12TargetPressureV1',contract=sha(a.contract),source=sha(Path(__file__)),
        target=sha(a.ct11/'target-v1/results.json'),preflight=sha(a.ct10/'target-preflight-v1/results.json'),
        reused_native_packets=dict(path=str(a.reuse_packets),files=reused),initial=sha(initial_path),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),scope='read-only captured target operators/states; extra diagnostic solves only'))
    PROGRESS=a.output/'last-solve.json'
    initial_mom=initial['ct10_momentum'];actual0={k:np.asarray(v) for k,v in initial_mom['actual'].items()}
    analytic0={k:np.asarray(v) for k,v in initial_mom['analytic'].items()};nm=initial_mom['momentum_size']
    cache={exact_vector_hash(init):(initial_mom,dict(path=str(initial_path),sha256=sha(initial_path),reused=True))}
    export_count=0;reused_count=0

    def loads(state,t,template,label):
        nonlocal export_count,reused_count
        h=exact_vector_hash(state)
        if h in cache:return cache[h]
        req=replay(initial,list(state),t,[],export_faces_only=True,export_ct10_momentum=True,ct10_export_points=False)
        if label in reused:
            dest=a.reuse_packets/label
            require(sha(dest/'case.xml')==sha(template) and n.read(dest/'loads-request.json')==req,'foreign reused native state/configuration')
            q=n.read(dest/'loads.json');e=n.read(dest/'loads-execution.json');reused_count+=1
        else:
            export_count+=1;require(export_count<=10,'prospective aggregate export count exceeded')
            dest=a.output/label;dest.mkdir();shutil.copy2(template,dest/'case.xml')
            print('CT12 native load',label,t,flush=True)
            q,e=probe(a.runtime,dest,req,'loads')
        require(e['exit_code']==0,'native target load export:'+str(q.get('failure')))
        require(q['observation']['state']==list(state),'read-only observer changed input')
        m=q['ct10_momentum'];require(m['maximum_reconstruction_error']<=2e-10,'native all-row load reconstruction')
        require(m['maximum_input_lift_difference']==0.,'nonzero effective lift requires separate attribution')
        meta=dict(path=str(dest/'loads.json'),sha256=sha(dest/'loads.json'),state_sha256=h,execution=e,
            all_row_reconstruction=m['maximum_reconstruction_error'],effective_lift=m['maximum_input_lift_difference'],
            prescribed_reaction_maximum=maximum(m['prescribed_reactions']),shared_trace_ids=len(m['shared_velocity_trace_ids']))
        cache[h]=(m,meta);return cache[h]

    out=dict(schema='CT12TargetPressureV1',status='incomplete',members={},common={},stage_loads={},scope='diagnostic reconstruction is separate from unchanged CT11 pressure failure')
    vectors={};totals={};accum_by_time={};snapshots={};bands=None;negatives={}
    for dt,member in target['members'].items():
        folder=Path(member['directory']);file=folder/'case.evolution.json'
        require(sha(file)==member['native_report_sha256'],'changed sealed native report')
        report=n.read(file);accum={};tel=np.zeros(off[6]-off[4]);rows_out=[];times={}
        for i,attempt in enumerate(report['pimple_attempts']):
            f=Path(attempt['flow_snapshot']);require(sha(f)==attempt['flow_snapshot_sha256'],'changed retained flow snapshot')
            s=n.read(f);state=np.asarray(s['input']);fr=Frozen(s,initial,c,weights,count)
            print('CT12 frozen reconstruction',dt,s['time'],flush=True)
            nv=fr.nv;ids=fr.rows[nv:];require(np.array_equal(ids,np.arange(off[4],off[6])),'pressure numbering/absolute free mask changed')
            mom,meta=loads(state,s['time']-float(dt),folder/'case.xml','input-'+dt+'-'+str(i))
            actual={k:np.asarray(v) for k,v in mom['actual'].items()}
            require(set(actual)==set(actual0),'component coverage changed')
            spatial=sum(actual.values(),np.zeros(nm));storage=fr.rm-spatial[fr.rows[:nv]]
            storage_global=np.zeros(nm);storage_global[fr.rows[:nv]]=storage
            free_difference=state[fr.rows[:nv]]-np.asarray(s['history'])[fr.rows[:nv]]
            storage_expected=s['alpha']*(fr.M@free_difference)
            pred=-fr.direct('A',fr.rm);rhs=fr.mul('D',pred)+fr.rc
            dpi=fr.solve(rhs,'total');tight=fr.solve(rhs,'total-tighter',2e-14)
            native=np.asarray(s['pressure_increment']);error=fr.measure(dpi-native);unc=fr.measure(tight-dpi)
            require(error['maximum_scaled_value']<1e-8 and unc['maximum_scaled_value']<1e-8,'pressure reconstruction/repeat uncertainty')
            continuity=fr.solve(fr.rc,'continuity');st=fr.image(storage_global,'storage-lifting')
            continuity_unc=fr.measure(fr.solve(fr.rc,'continuity-tighter',2e-14)-continuity)
            storage_unc=fr.measure(fr.image(storage_global,'storage-lifting-tighter',2e-14)-st)
            require(max(continuity_unc['maximum_scaled_value'],storage_unc['maximum_scaled_value'])<1e-8,'small-RHS image uncertainty')
            signed={};per_component={}
            for key,v in actual.items():
                pa=fr.image(analytic0[key],'continuous/'+key)
                pr=fr.image(actual0[key]-analytic0[key],'representation/'+key)
                pc=fr.image(v,'current/'+key)
                signed['continuous/'+key]=pa;signed['representation/'+key]=pr;signed['evolution/'+key]=pc-pa-pr
                per_component[key]=dict(current=fr.measure(pc),initial_representation=fr.measure(pr),evolution=fr.measure(pc-pa-pr))
            signed.update(continuity=continuity,storage_lifting=st)
            summed=sum(signed.values(),np.zeros_like(dpi));sum_error=fr.measure(summed-native)
            require(sum_error['maximum_scaled_value']<1e-8,'signed finite-step pressure sum')
            mass_limit=fr.mass_image(spatial,'current-small-step-image')
            tel+=native
            accepted=np.asarray(report['accepted_states'][i+1]);absolute=accepted[ids]-init[ids]
            telescope=fr.measure(tel-absolute)
            require(telescope['maximum_scaled_value']<1e-8,'absolute pressure increment telescope')
            kf=Path(attempt['k_snapshot']);require(sha(kf)==attempt['k_snapshot_sha256'],'changed actual pre-k packet');ks=n.read(kf)
            stage_pressure=dict(input_to_previous_accepted=maximum(state[ids]-np.asarray(report['accepted_states'][i])[ids]),
                corrected_to_input_plus_increment=maximum(np.asarray(s['corrected_flow'])[ids]-state[ids]-native),
                post_flow_to_corrected=maximum(np.asarray(ks['post_flow'])[ids]-np.asarray(s['corrected_flow'])[ids]),
                pre_k_to_post_flow=maximum(np.asarray(ks['pre_k'])[ids]-np.asarray(ks['post_flow'])[ids]),
                accepted_to_pre_k=maximum(accepted[ids]-np.asarray(ks['pre_k'])[ids]))
            require(max(stage_pressure.values())<=2e-10,'unattributed pressure stage change')
            row=dict(time=s['time'],snapshot=dict(path=str(f),sha256=sha(f)),load=meta,
                alpha=s['alpha'],C_maximum=smx(fr.C),unused_gauge=False,free_velocity=nv,free_pressure=len(ids),
                block_reconstruction=fr.blocks,predictor_reconstruction=maximum(pred-np.asarray(s['predictor_increment'])),
                initial_weak_continuity_maximum=maximum(fr.rc),native_increment=fr.measure(native),
                reconstructed_increment=fr.measure(dpi),reconstruction=error,tighter_repeat=unc,signed_sum_error=sum_error,
                continuity=fr.measure(continuity),storage_and_lifting=fr.measure(st),
                continuity_tighter=continuity_unc,storage_lifting_tighter=storage_unc,
                backward_velocity_storage_maximum=maximum(storage_expected),storage_lift_reconstruction=maximum(storage-storage_expected),
                small_step_momentum_image=fr.measure(mass_limit),finite_alpha_remainder=fr.measure(dpi-continuity-mass_limit),
                telescope=telescope,absolute_pressure_change=fr.measure(absolute),stage_pressure=stage_pressure,
                trace={k:attempt[k] for k in ('imposed_trace_defect','refreshed_trace_defect','trace_physical_storage_contribution_max')},
                component_images=per_component,predictor_factor_nonzeros=fr.al.L.nnz+fr.al.U.nnz)
            if (dt=='0.001' and i==0) or (dt=='0.0005' and i==len(report['pimple_attempts'])-1):
                repeats={}
                for key in ('outer/volume/nonlinear_normal','outer/volume/advection','outer/outlet/natural_advection'):
                    pv=signed['continuous/'+key]+signed['representation/'+key]+signed['evolution/'+key]
                    q=fr.image(actual[key],'component-tighter/'+key,2e-14);err=fr.measure(q-pv)
                    require(err['maximum_scaled_value']<1e-8,'component pressure-image tighter repeat:'+key)
                    repeats[key]=err
                row['component_tighter_repeats']=repeats
            if not out['members']:
                # Discriminating actions, no alternate governing update.
                direction=np.sin(np.arange(len(ids))*.371)*c['ustar']**2
                right=fr.S@direction;wrong=fr.D@(fr.G@direction/fr.M.diagonal())/fr.alpha-fr.C@direction
                neg=maximum(wrong-right);require(neg>2e-10,'diagonal-mass control has no action signal')
                omitted=fr.G+fr.D.T;rough_signal=maximum(omitted@direction)
                require(rough_signal>2e-10,'rough-pressure action control has no signal')
                negatives=dict(diagonal_governing_mass=dict(maximum_action_defect=neg,limit=2e-10),omitted_rough_pressure=dict(maximum_action_defect=rough_signal,limit=2e-10))
            if dt=='0.001' and i==0:
                # Exact frozen coupled correction through an A-inverse action,
                # without a global coupled factor or a dense Schur matrix.
                exactop=LinearOperator(fr.S.shape,matvec=lambda p:fr.mul('D',fr.direct('A',fr.mul('G',p)))-fr.mul('C',p))
                exactpi=fr.solve(rhs,'exact-frozen-coupled',operator=exactop)
                exact_repeat=fr.measure(fr.solve(rhs,'exact-frozen-coupled-tighter',rtol=2e-14,operator=exactop)-exactpi)
                require(exact_repeat['maximum_scaled_value']<1e-8,'exact coupled pressure repeat uncertainty')
                exactu=pred-fr.direct('A',fr.mul('G',exactpi));res=fr.F@np.r_[exactu,exactpi]+np.r_[fr.rm,fr.rc]
                require(maximum(res)<2e-10,'exact coupled diagnostic all-row residual')
                row['exact_coupled']=dict(pressure=fr.measure(exactpi),split_difference=fr.measure(dpi-exactpi),maximum_residual=maximum(res),tighter_repeat=exact_repeat)
                # Stream one radial band vector at a time from local sparse supports.
                bounds={str(e['geometry']):(e['lower'][2],e['upper'][2]) for e in mom['elements']}
                static=fr.mass_image(np.asarray(mom['total_actual'])-np.asarray(mom['total_analytic']),'initial-representation-total')
                den=fr.inner(static,static);require(den>0,'unresolved initial representation signal')
                bandrows=[];bandvectors=[]
                for band in sorted(set(bounds.values())):
                    load=np.zeros(nm)
                    for loc in mom['local_sparse']:
                        if bounds[loc['geometry']]!=band:continue
                        for term in loc['terms'].values():np.add.at(load,np.asarray(term['ids'],int),np.asarray(term['values'])*(-1 if loc['analytic'] else 1))
                    pv=fr.mass_image(load,'band/'+str(band));bandvectors.append(pv)
                    bandrows.append(dict(lower=band[0],upper=band[1],signed_alignment=fr.inner(pv,static)/den,image=fr.measure(pv),maximum_load=maximum(load)))
                berr=fr.measure(sum(bandvectors,np.zeros_like(static))-static);require(berr['maximum_scaled_value']<1e-8,'radial-band signed sum')
                bands=dict(rows=bandrows,signed_sum_error=berr,total_image=fr.measure(static))
            for key,v in signed.items():accum[key]=accum.get(key,np.zeros_like(v))+v
            if round(s['time'],12) in (.001,.002):times[str(round(s['time'],12))]={k:v.copy() for k,v in accum.items()}
            row['solve_checks']=fr.receipts;rows_out.append(row)
            write_json(a.output/'progress.json',dict(member=dt,completed=len(rows_out),last=row,scope='partial diagnostic,not qualified'))
            # Native component sum/faces/owners remain in external bound packets.
            if i==len(report['pimple_attempts'])-1:
                endmom,endmeta=loads(accepted,s['time'],folder/'case.xml','accepted-'+dt)
                out['stage_loads'][dt+'/accepted']=dict(metadata=endmeta,component_load_changes={k:maximum(np.asarray(v)-np.asarray(mom['actual'][k])) for k,v in endmom['actual'].items()})
                if dt=='0.0005':
                    last=mom
                    for name in ('post_flow','pre_k'):
                        m,meta_stage=loads(np.asarray(ks[name]),s['time'],folder/'case.xml',name+'-'+dt)
                        out['stage_loads'][dt+'/'+name]=dict(metadata=meta_stage,component_load_changes={k:maximum(np.asarray(v)-np.asarray(last['actual'][k])) for k,v in m['actual'].items()});last=m
            out['members'][dt]=dict(steps=rows_out,report_sha256=sha(file),actual_initial_sha256=exact_vector_hash(init),final_absolute_pressure=fr.measure(absolute))
            # Persist derived work before any later diagnostic can fail. Native
            # packets stay immutable and are rebound by request/state hashes.
            write_json(a.output/'partial-results.json',out)
            write_json(a.output/'partial-bands.json',bands)
            np.savez(a.output/('pressure-vectors-'+dt+'.npz'),**{t+'|'+k:v for t,terms in times.items() for k,v in terms.items()})
            snapshots[dt]=s;del fr;gc.collect()
        accum_by_time[dt]=times
        vectors[dt]={t:{k:v.tolist() for k,v in terms.items()} for t,terms in times.items()}
    # Fixed physical pressure norm and actual free pressure space for signed pairing.
    fr=Frozen(snapshots['0.0005'],initial,c,weights,count)
    # Nonlinear governing-quadrature control at identical stored coefficients;
    # observation rows and the pressure-image operator remain fixed.
    fine=target['members']['0.0005'];file=Path(fine['directory'])/'case.evolution.json';report=n.read(file)
    endstate=report['state'];native_mom,_=cache[exact_vector_hash(endstate)];quad=[]
    for orders in ((48,12,12,16,64),(64,16,12,20,64)):
        dest=a.output/('held-q'+str(orders[0]));dest.mkdir();root=ET.parse(Path(fine['directory'])/'case.xml').getroot()
        for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,key,value)
        write(root,dest/'case.xml');req=replay(initial,endstate,.002,[],export_faces_only=True,export_ct10_momentum=True,ct10_export_points=False)
        label=dest.name
        if label in reused:
            olddest=a.reuse_packets/label
            require(sha(olddest/'case.xml')==sha(dest/'case.xml') and n.read(olddest/'loads-request.json')==req,'foreign held-quadrature request')
            q=n.read(olddest/'loads.json');e=n.read(olddest/'loads-execution.json');packet=olddest/'loads.json'
        else:
            q,e=probe(a.runtime,dest,req,'loads');packet=dest/'loads.json'
        if e['exit_code']!=0:
            # Preserve the failed high-order capacity control. It cannot become
            # passing quadrature evidence, but must not erase the six completed
            # pressure reconstructions or trigger another native export.
            require(q.get('failure')=='std::bad_alloc','unexpected held-state native failure')
            quad.append(dict(status='blocked',orders=orders,execution=e,packet=str(packet),packet_sha256=sha(packet),failure=q['failure'],address_space_GiB=16,scope='unresolved governing-quadrature capacity; no physical continuation admission'))
            continue
        require(q['observation']['state']==endstate,'held target quadrature input changed')
        m=q['ct10_momentum'];delta=np.asarray(m['total_actual'])-np.asarray(native_mom['total_actual'])
        pv=fr.image(delta,'held-q'+str(orders[0]));err=fr.measure(pv)
        require(m['maximum_reconstruction_error']<=2e-10 and err['maximum_scaled_value']<1e-8,'held-state quadrature uncertainty')
        quad.append(dict(status='pass',orders=orders,execution=e,packet=str(packet),packet_sha256=sha(packet),maximum_load_change=maximum(delta),pressure_image=err,reconstruction=m['maximum_reconstruction_error']))
    for t in sorted(set(accum_by_time['0.001'])&set(accum_by_time['0.0005'])):
        aa=accum_by_time['0.001'][t];bb=accum_by_time['0.0005'][t]
        delta={k:aa[k]-bb[k] for k in aa};total=sum(delta.values(),np.zeros(fr.S.shape[0]));den=fr.inner(total,total)
        require(den>0,'unresolved common-time pressure signal')
        table={k:dict(signed_alignment=fr.inner(v,total)/den,pressure=fr.measure(v)) for k,v in delta.items()}
        direct=next(x for x in target['common_times'] if abs(x['time']-float(t))<1e-13)['fields']['Pi']
        rms=fr.measure(total)['whole']['weighted_RMS'];require(abs(rms-direct['weighted_RMS'])/c['ustar']**2<1e-8,'paired pressure sum differs from sealed physical norm')
        groups={}
        for key,v in delta.items():
            group=key.split('/')[0];groups[group]=groups.get(group,np.zeros_like(v))+v
        out['common'][t]=dict(pressure=fr.measure(total),sealed_CT11=direct,components=table,
            groups={k:dict(signed_alignment=fr.inner(v,total)/den,pressure=fr.measure(v)) for k,v in groups.items()})
    out.update(status='pass' if all(q['status']=='pass' for q in quad) else 'reconstruction_pass_quadrature_blocked',reconstruction_status='pass',
        quadrature_status='pass' if all(q['status']=='pass' for q in quad) else 'blocked',bands=bands,negative_actions=negatives,native_export_count=export_count,reused_native_export_count=reused_count,held_quadrature=quad,
        pressure_convention='absolute total-traction level; common means and variations retained without regauging',
        historical_qualification=target['qualification'],continuation_admitted=False)
    write_json(a.output/'pressure-vectors.json',vectors);write_json(a.output/'results.json',out)
    manifest.seal(a.output,['results.json','pressure-vectors.json'])
    print('CT12 target pressure reconstruction complete',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
