#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only CT10 signed native loads, independent tensors and pressure images.

Run with ct_accepted3_command.py. No state, boundary, checkpoint or governing
operator is changed. Local loads have sparse native support; Schur actions use
consistent mass solves without materializing M^-1 G or a dense Schur matrix.
"""
import argparse,copy,math,os,resource,shutil
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.sparse import csr_matrix,diags
from scipy.sparse.linalg import splu,LinearOperator,gmres
import ct_accepted8_native as n
from ct_accepted8_reference import column,reference
from validate_a93_a95_e94e95 import oracle
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import probe,require,row_value
from ct_accepted_pilot import replay
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity


def matrix(info):
    p=Path(info['file']);require(sha(p)==info['sha256'],'changed captured matrix')
    with p.open('rb') as f:
        require(f.read(8)==b'RANSCSR1','not a native CSR packet')
        nr,nc,nz=map(int,np.fromfile(f,dtype='<u8',count=3))
        ptr=np.fromfile(f,dtype='<u8',count=nr+1);ids=np.fromfile(f,dtype='<u4',count=nz);v=np.fromfile(f,dtype='<f8',count=nz)
        require((nr,nc)==(info['rows'],info['columns']) and len(v)==nz and len(ptr)==nr+1 and not f.read(1),'CSR shape/length mismatch')
        require(np.isfinite(v).all() and ptr[0]==0 and ptr[-1]==nz and np.all(ptr[1:]>=ptr[:-1]) and np.all(ids<nc),'invalid CSR primitives')
    return csr_matrix((v,ids,ptr),shape=(nr,nc))


class PressureImages:
    def __init__(self,snapshot,initial,c,weights,count):
        self.s=snapshot;self.c=c;self.scale=c['ustar']**2;self.rows=np.asarray(snapshot['rows'],int);self.nv=snapshot['velocity_size']
        self.M=matrix(snapshot['free_velocity_mass']);self.lu=splu(self.M.tocsc())
        op=snapshot['pressure_operator'];self.G=matrix(op['G']);self.D=matrix(op['D']);self.C=matrix(op['C'])
        require(not np.any(self.C.data),'CT10 W image requires actual C=0')
        require(len(initial['observation']['state'])-1 not in self.rows,'unused gauge entered physical pressure system')
        inv={int(v):i for i,v in enumerate(self.rows[self.nv:])};rr=[];cc=[];vv=[];off=initial['observation']['offsets']
        self.regions=np.asarray(initial['observation']['grid_regions'][:count]);self.weights=np.asarray(weights)
        require(count==len(weights),'pressure observation weight count')
        for i,(q,r) in enumerate(zip(initial['test_basis'][:count],self.regions,strict=True)):
            for j,v in zip(q['pressure']['ids'],q['pressure']['phi'],strict=True):rr.append(i);cc.append(inv[off[4+int(r)]+j]);vv.append(v)
        self.B=csr_matrix((vv,(rr,cc)),shape=(count,len(inv)))
        ir=[];ic=[];iv=[];self.ifaces=[f for f in initial['observation']['faces'] if f['interface']]
        for i,f in enumerate(self.ifaces):
            for j,v in zip(f['pressure']['ids'],f['pressure']['phi'],strict=True):ir.append(i);ic.append(inv[off[4+f['region']]+j]);iv.append(v)
        self.IB=csr_matrix((iv,(ir,ic)),shape=(len(self.ifaces),len(inv)))
        self.W=LinearOperator((len(inv),len(inv)),matvec=lambda p:self.D@self.lu.solve(self.G@p))
        pc=(self.D@diags(1/self.M.diagonal())@self.G).tocsc();self.pc=splu(pc)
        self.P=LinearOperator(self.W.shape,matvec=self.pc.solve);self.receipts=[]

    def measure(self,p):
        v=self.B@p;out={}
        for name,mask in [('whole',np.ones(len(v),bool)),('inner',self.regions==0),('outer',self.regions==1)]:
            w=self.weights[mask];x=v[mask];mean=float(w@x/w.sum())
            out[name]=dict(weighted_RMS=float(np.sqrt(w@(x*x)/w.sum())),signed_mean=mean,spatial_variation_RMS=float(np.sqrt(w@((x-mean)**2)/w.sum())),scale=self.scale)
        out['maximum_scaled_value']=float(max(abs(v),default=0))/self.scale
        iv=self.IB@p;pairs={}
        for f,x in zip(self.ifaces,iv,strict=True):pairs.setdefault(round(f['coordinate'][0],14),[]).append((f['region'],float(x),f['weight']))
        jumps=[]
        for qs in pairs.values():
            require(len(qs)==2 and {q[0] for q in qs}=={0,1},'pressure image interface side coverage')
            inn=next(q for q in qs if not q[0]);outer=next(q for q in qs if q[0]);jumps.append((inn[1]-outer[1],inn[2]))
        out['interface_inner_minus_outer']=dict(signed_mean=sum(x*w for x,w in jumps)/sum(w for _,w in jumps),weighted_RMS=math.sqrt(sum(x*x*w for x,w in jumps)/sum(w for _,w in jumps)),scale=self.scale)
        return out

    def inner(self,a,b):return float((self.B@a)@(self.weights*(self.B@b)))

    def image(self,global_load,label,rtol=2e-12):
        load=np.asarray(global_load)[self.rows[:self.nv]];rhs=-(self.D@self.lu.solve(load));mag=float(np.linalg.norm(rhs))
        if mag==0:return np.zeros(self.W.shape[0])
        ticks=[]
        x,info=gmres(self.W,rhs/mag,M=self.P,rtol=rtol,atol=0.,restart=80,maxiter=1000,callback=lambda v:ticks.append(v),callback_type='pr_norm')
        x*=mag;defect=self.W@x-rhs
        rel=float(np.linalg.norm(defect))/mag
        require(info==0 and rel<2e-10,'diagnostic true Schur residual:'+label+':'+str(rel))
        mass_defect=self.M@self.lu.solve(load)-load
        self.receipts.append(dict(label=label,iterations=len(ticks),rhs_norm=mag,true_relative_residual=rel,maximum_mass_defect=float(max(abs(mass_defect),default=0))))
        return x


def tensors_and_loads(packet,initial,c):
    """Independent constitutive/profile formulas applied to exported native rows."""
    p=packet['ct10_momentum'];state=np.asarray(initial['observation']['state']);off=p['offsets'];nm=p['momentum_size']
    truth=column(c['K0'],c['H'],c['roughness_input'],c['nu'],c['displacement']);ref=lru_cache(None)(lambda z,r:reference(truth,z,bool(r)))
    cols=lru_cache(None)(lambda K:column(K,c['H'],c['roughness_input'],c['nu'],c['displacement']))
    loads=[{key:np.zeros(nm) for key in p['actual']} for _ in range(2)];defects={};partition=0.;gradconst=0.;volume=[0.,0.]
    def put(version,key,q,values,gradient,weight):
        ids=np.asarray(q['ids'],int);phi=np.asarray(q['phi']);dphi=np.asarray(q['derivative'])
        for a,axis in enumerate((0,2)):
            v=weight*(phi*values[a]+dphi.T@gradient[a]);np.add.at(loads[version][key],off[a]+ids,v)
    def check(name,a,b):defects[name]=max(defects.get(name,0.),float(np.max(abs(np.asarray(a)-np.asarray(b)))))
    samplers=[]
    for q in p['samplers']:
        v=q['outer_velocity_k'];pr=q['outer_pressure'];K,gk=row_value(v,state,off[2]);_,gp=row_value(pr,state,off[5]);U,gu=row_value(v,state,off[0])
        fp=np.array([gp[0]-2*gk[0]/3,0.,0.]);fi=np.array([2*U*gu[0],0.,0.]);check('rough_pressure',fp,q['pressure_force']);check('rough_inertia',fi,q['inertial_force']);check('endpoint_K',K,q['K']);samplers.append((K,fp,fi))
    for point in p['volume_points']:
        r=point['region'];q=point['velocity_k'];z=q['coordinate'][2];weight=point['weight'];reg='outer/' if r else 'inner/'
        partition=max(partition,abs(sum(q['phi'])-1));gradconst=max(gradconst,float(max(abs(np.sum(q['derivative'],axis=1)))));volume[r]+=weight
        for version,record in enumerate(point['actual_analytic']):
            if version:
                value=ref(z,r);u=np.array([value[0],0.,0.]);G=np.zeros((3,3));G[0,2]=value[4];k,ep,pi=value[1:4];nt=0. if r else float(truth.D(z)-c['nu']);fp=fi=np.zeros(3)
            else:
                u=np.zeros(3);G=np.zeros((3,3))
                for block,axis in enumerate((0,2)):u[axis],G[axis]=row_value(q,state,off[block])
                k,_=row_value(q,state,off[2]);pi,_=row_value(point['pressure'],state,off[4+r])
                if r:ep,_=row_value(point['epsilon_row'],state,off[3]);nt=0.;fp=fi=np.zeros(3)
                else:
                    K,fp,fi=samplers[point['sampler']];cl=cols(K);ep=cl.A/(max(z,cl.cutoff)-c['displacement']);nt=float(cl.D(z)-c['nu'])
            mol=np.zeros((3,3));lin=np.zeros((3,3));nl=np.zeros((3,3))
            if r:
                law=oracle('A95',k,ep,G);nl=-law['Q'];lin=2*law['nut']*((G+G.T)/2-np.eye(3)*np.trace(G)/3)
            else:mol=c['nu']*(G+G.T);lin=nt*(G+G.T)
            for key,val in dict(u=u,G=G,k=k,epsilon=ep,Pi=pi,nt_aux=nt,molecular=mol,linear=lin,nonlinear=nl).items():check(key,val,record[key])
            vv=[np.zeros(3),np.zeros(3)];adv=[-u[axis]*u for axis in (0,2)]
            for a,axis in enumerate((0,2)):vv[a][axis]=-pi
            put(version,reg+'volume/augmented_pressure',q,[0,0],vv,weight);put(version,reg+'volume/advection',q,[0,0],adv,weight)
            for label,T in [('molecular',mol),('linear',lin),('nonlinear',nl)]:
                normal=np.zeros((2,3));shear=np.zeros((2,3))
                for a,axis in enumerate((0,2)):
                    normal[a,axis]=T[axis,axis];shear[a]=T[axis];shear[a,axis]=0.
                put(version,reg+'volume/'+label+'_normal',q,[0,0],normal,weight);put(version,reg+'volume/'+label+'_shear',q,[0,0],shear,weight)
            for label,f in [('pressure',fp),('inertia',fi)]:put(version,reg+'volume/rough_'+label,q,[-f[0],-f[2]] if not r and z<=c['roughness_input'] else [0,0],np.zeros((2,3)),weight)
    Q=oracle('A95',truth.K,truth.A/(truth.H-truth.disp),np.array([[0,0,truth.um/(truth.p['kappa']*(truth.H-truth.disp))],[0,0,0],[0,0,0]]))['Q']
    for f in p['natural_faces']:
        q=f['row'];r=f['region'];reg='outer/' if r else 'inner/';area=np.asarray(f['area']);z=q['coordinate'][2];kind=f['kind'];weight=f['weight']
        tr=np.array([Q[2,2]-Q[0,0] if r else 0.,truth.t0]) if kind=='outlet' else np.array([truth.t0,0.])
        for version in (0,1):
            u=np.array([ref(z,r)[0],0.,0.]) if version else np.array([row_value(q,state,off[0])[0],0.,row_value(q,state,off[1])[0]])
            put(version,reg+kind+'/natural_traction',q,-np.linalg.norm(area)*tr,np.zeros((2,3)),weight)
            put(version,reg+kind+'/natural_advection',q,u[[0,2]]*(u@area),np.zeros((2,3)),weight)
    comparisons={}
    for version,name in enumerate(('actual','analytic')):
        err={key:float(max(abs(v-np.asarray(p[name][key])))) for key,v in loads[version].items()};comparisons[name]=dict(maximum_defect=max(err.values()),components=err)
    require(max(defects.values())<2e-9 and max(x['maximum_defect'] for x in comparisons.values())<2e-10,'independent tensor/load disagreement')
    return dict(status='pass',primitive_defects=defects,load_reconstruction=comparisons,partition_defect=partition,constant_gradient_defect=gradconst,integrated_volume=volume,expected_volume=[(c['x'][-1]-c['x'][0])*c['H'],(c['x'][-1]-c['x'][0])*(c['Z']-c['H'])],scope='independent profile quadrature and spectral-tensor formulas; native test rows, owners, coordinates and weights unchanged')


def attribution(p,images,c):
    actual={k:np.asarray(v) for k,v in p['actual'].items()};analytic={k:np.asarray(v) for k,v in p['analytic'].items()};delta={k:actual[k]-analytic[k] for k in actual}
    total={name:sum(vals.values(),np.zeros(p['momentum_size'])) for name,vals in [('actual',actual),('analytic',analytic),('delta',delta)]}
    error=float(max(abs(total['actual']-np.asarray(p['native_raw']))));require(error<=2e-10,'native signed sum failed')
    frozen_error=float(max(abs(total['actual'][images.rows[:images.nv]]-np.asarray(images.s['frozen_flow_residual'][:images.nv]))))
    require(frozen_error<=2e-10,'initial original momentum differs from actual frozen momentum load')
    full={name:images.image(v,name) for name,v in total.items()};den=images.inner(full['delta'],full['delta']);require(den>0,'no resolved representation pressure image')
    vectors={};table={};off=p['offsets']
    for key in actual:
        pa=images.image(actual[key],'actual/'+key);pb=images.image(analytic[key],'analytic/'+key);dv=pa-pb;vectors[key]=dv
        table[key]=dict(actual=images.measure(pa),analytic=images.measure(pb),difference=images.measure(dv),signed_alignment=images.inner(dv,full['delta'])/den,maximum_load_difference=float(max(abs(delta[key]))),free_signed_load_sum=float(np.sum(delta[key][images.rows[:images.nv]])))
    sumerr=images.measure(sum(vectors.values(),np.zeros_like(full['delta']))-full['delta'])
    require(sumerr['maximum_scaled_value']<1e-8,'signed pressure image sum disagreement')
    # Per-axis image decomposition makes normal versus streamwise sources explicit.
    axes={}
    for name,lo,hi in [('streamwise',0,off[1]),('normal',off[1],p['momentum_size'])]:
        v=total['delta'].copy();v[:lo]=0.;v[hi:]=0.;pv=images.image(v,name);axes[name]=dict(image=images.measure(pv),signed_alignment=images.inner(pv,full['delta'])/den)
    bounds={str(e['geometry']):(e['lower'][2],e['upper'][2]) for e in p['elements']};bands={};local_sum=[np.zeros(p['momentum_size']),np.zeros(p['momentum_size'])]
    for local in p['local_sparse']:
        version=int(local['analytic']);band=bounds[local['geometry']];v=bands.setdefault(band,[np.zeros(p['momentum_size']),np.zeros(p['momentum_size'])])[version]
        for term in local['terms'].values():
            ids=np.asarray(term['ids'],int);values=np.asarray(term['values']);np.add.at(v,ids,values);np.add.at(local_sum[version],ids,values)
    for version,name in enumerate(('actual','analytic')):require(float(max(abs(local_sum[version]-total[name])))<=2e-10,'local sparse supports do not reconstruct global '+name)
    rows=[];band_vectors=[]
    for band,(b,a) in sorted(bands.items()):
        pv=images.image(b-a,'band/'+str(band));band_vectors.append(pv);rows.append(dict(lower=band[0],upper=band[1],signed_alignment=images.inner(pv,full['delta'])/den,pressure=images.measure(pv),load_maximum=float(max(abs(b-a)))))
    banderr=images.measure(sum(band_vectors,np.zeros_like(full['delta']))-full['delta']);require(banderr['maximum_scaled_value']<1e-8,'band image sum disagreement')
    tight=images.image(total['delta'],'delta-tighter',rtol=2e-14);unc=images.measure(tight-full['delta']);require(unc['maximum_scaled_value']<1e-8,'pressure solve uncertainty too large')
    gram=[[images.inner(a,b)/den for b in band_vectors] for a in band_vectors]
    signal=np.array([abs(q['signed_alignment']) for q in rows]);threshold=.8*signal.sum();choices=[]
    for i in range(len(rows)):
        for j in range(i,len(rows)):
            if signal[i:j+1].sum()>=threshold:choices.append((j-i+1,i,j));break
    width,i,j=min(choices);uniform=width>math.ceil(len(rows)/3);selected=list(range(len(rows))) if uniform else list(range(i,j+1))
    knots=sorted(set(c['knots']+[(rows[k]['lower']+rows[k]['upper'])/2 for k in selected]))
    selection=dict(schema='CT10SingleRadialContrastV1',status='selected_before_evolution',rule='contract80percent_absolute_signed_alignment,smallest_contiguous_band;uniform_if_more_than_one_third',uniform=uniform,diagnosed_band=[rows[i]['lower'],rows[j]['upper']],selected_bands=[[rows[k]['lower'],rows[k]['upper']] for k in selected],knots=knots,old_knots=c['knots'],degree=4,pressure_degree=2,continuous_physical_inputs_unchanged=True,signal_captured=float(signal[i:j+1].sum()/signal.sum()))
    names=list(vectors);cg=[[images.inner(vectors[x],vectors[y])/den for y in names] for x in names]
    return dict(schema='CT10PressureAttributionV1',status='pass',native_reconstruction=error,frozen_momentum_agreement=frozen_error,component_images=table,component_Gram_order=names,component_normalized_Gram=cg,total_images={k:images.measure(v) for k,v in full.items()},axis_images=axes,bands=rows,band_normalized_Gram=gram,signed_component_sum_error=sumerr,signed_band_sum_error=banderr,tighter_repeat_difference=unc,pressure_solves=images.receipts,selection=selection),dict(total={k:v.tolist() for k,v in full.items()},components={k:v.tolist() for k,v in vectors.items()},bands=[v.tolist() for v in band_vectors])


def main():
    a=argparse.ArgumentParser(description=__doc__)
    for key in ('ct8','ct9','runtime','contract','output'):a.add_argument('--'+key,type=Path,required=True)
    a=a.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    pre=sealed(a.ct8/'preflight-v6');base=sealed(a.ct9/'equivalence-v2');c=pre['case'];member=base['members']['0.001'];folder=Path(member['directory']);native=n.read(folder/'case.evolution.json');require(sha(folder/'case.evolution.json')==member['native_report_sha256'],'changed baseline')
    initial_path=a.ct9/'equivalence-v2/zero-0.001/initial.json';initial=n.read(initial_path);attempt=next(x for x in native['pimple_attempts'] if x['scheduled_complete']);sp=Path(attempt['flow_snapshot']);require(sha(sp)==attempt['flow_snapshot_sha256'],'changed snapshot');snap=n.read(sp)
    require(snap['input']==initial['observation']['state'],'initial frozen snapshot differs from actual projected state')
    m=RequestManifest(a.output,dict(schema='CT10SignedDiagnosisV1',contract=sha(a.contract),case=c,source=sha(Path(__file__)),diagnostic_runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),governing_runtime=runtime_identity(a.ct9/'build-v8/runtime/IncNavierStokesSolver'),baseline_report=sha(folder/'case.evolution.json'),initial=sha(initial_path),snapshot=sha(sp),state=exact_vector_hash(snap['input'])))
    dest=a.output/'initial';dest.mkdir();shutil.copy2(folder/'case.xml',dest/'case.xml')
    req=replay(initial,snap['input'],0.,[],export_ct10_momentum=True,ct10_export_points=True)
    packet,e=probe(a.runtime,dest,req,'loads');require(e['exit_code']==0,'native CT10 load export failed:'+str(packet.get('failure')));require(packet['observation']['state']==snap['input'],'diagnostic changed coefficients')
    p=packet['ct10_momentum'];independent=tensors_and_loads(packet,initial,c);write_json(a.output/'independent.json',independent)
    result,vectors=attribution(p,PressureImages(snap,initial,c,pre['weights'],pre['global_count']),c)
    result.update(independent=independent,execution=e,packet_sha256=sha(dest/'loads.json'),baseline_native_report_sha256=sha(folder/'case.evolution.json'),snapshot_sha256=sha(sp),initial_sha256=sha(initial_path),actual_initial_sha256=exact_vector_hash(snap['input']),native_counts=dict(stored=len(snap['input']),free_velocity=snap['velocity_size'],free_pressure=len(snap['rows'])-snap['velocity_size']),native_point_count=len(p['volume_points']),local_sparse_entries=sum(len(t['ids']) for q in p['local_sparse'] for t in q['terms'].values()))
    write_json(a.output/'pressure-vectors.json',vectors);write_json(a.output/'selection.json',result['selection']);write_json(a.output/'results.json',result);m.seal(a.output,['results.json','selection.json','independent.json','pressure-vectors.json']);print('CT10 diagnosis complete',result['selection'],flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
