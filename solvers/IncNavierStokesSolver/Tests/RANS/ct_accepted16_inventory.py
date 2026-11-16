#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh batched full A95 inventory; original scalar witnesses and fallbacks."""
import argparse,copy,ctypes,hashlib,math,os,resource,subprocess,sys,time
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
import ct_accepted3_metrics as metrics
from ct_accepted15_python import PreparedRows
from ct_accepted15_audits import bundle
from ct_accepted16_blas_audit import ExactAuditRows,exact_equal
from ct_accepted16_blas import BLAS,BLAS_SHA
from ct_accepted_open import require,run
from ct_accepted_feasibility import write_json
from rans_trace_index import TraceIndex
from rans_request_manifest import RequestManifest,sha,runtime_identity

ORIGINAL_INVENTORY=metrics.inventory


class BatchInventoryRows(ExactAuditRows):
    policy='CT16-FreshBatchedA95InventoryV1'
    # This is an algorithm-selection guard, not a changed admissibility limit.
    # Difficult near-boundary states use the complete original scalar route.
    guard=2e-10
    def __init__(self,*args):
        super().__init__(*args);self.saved_inventories=[];self.inventory_calls=0;self.scalar_fallbacks=0;self.scalar_witness_evaluations=0;self.last_arrays=None
        f64=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
        self.lib.ct16_norm9.argtypes=[ctypes.c_int64,f64,f64];self.lib.ct16_norm9.restype=ctypes.c_int
    def norm9(self,m):
        m=np.ascontiguousarray(m);out=np.empty(len(m));require(self.lib.ct16_norm9(len(m),m,out)==0,'bounded nine-entry norm rejected');return out
    def invariants(self,G,k,e):
        div=np.trace(G,axis1=1,axis2=2);I=np.eye(3)[None,:,:]
        S=(G+G.transpose(0,2,1))/2-I*div[:,None,None]/3;O=(G-G.transpose(0,2,1))/2
        ev=np.linalg.eigvalsh(S);sn=self.norm9(S);on=self.norm9(O);tau=k/e
        if np.any(sn==0) or not all(np.isfinite(a).all() for a in (div,ev,sn,on,tau)):return None
        C=1/(6.5+3*ev[:,-1]/sn*np.hypot(sn,on)*tau)
        rad=1-9*C*C*sn*sn*tau*tau
        return div,S,O,ev,sn,on,tau,C,rad
    def inventory(self,rows,coords,state,offsets,model,sides=None):
        if not (rows is self.initial['test_basis'] and coords is self.initial['observation']['grid'] and tuple(offsets)==self.offsets and model=='A95' and (sides is None or sides is self.initial['observation']['grid_sides'])):
            return ORIGINAL_INVENTORY(rows,coords,state,offsets,model,sides)
        self.inventory_calls+=1;self.last_arrays=None
        def original():
            self.scalar_fallbacks+=1
            return ORIGINAL_INVENTORY(rows,coords,state,offsets,model,sides)
        s=np.ascontiguousarray(state,dtype=np.float64)
        if s.shape!=(self.state_size,):return original()
        count=self.grid_size;v=np.empty((count,5));g=np.empty((count,5,3))
        require(self.lib.ct16_plan_grid(self.handle,count,s,len(s),v,g)==0,'inventory actual field input rejected')
        outer=np.flatnonzero(self.present[:count]);inner=np.flatnonzero(~self.present[:count].astype(bool))
        if not len(outer) or not np.isfinite(v).all() or not np.isfinite(g).all() or np.any(v[inner,2]<0) or np.any(v[outer,2]<=0) or np.any(v[outer,3]<=0):return original()
        k=v[outer,2];e=v[outer,3];G=np.zeros((len(outer),3,3));G[:,0,:]=g[outer,0,:];G[:,2,:]=g[outer,1,:]
        try:
            with np.errstate(all='ignore'):
                inv=self.invariants(G,k,e)
                if inv is None:return original()
                div,S,O,ev,sn,on,tau,C,rad=inv
                # The original scalar oracle uses math.hypot here, while its
                # independent invariant uses NumPy hypot. Preserve both paths.
                hm=np.asarray([math.hypot(float(a),float(b)) for a,b in zip(sn,on,strict=True)])
                cmu=1/(6.5+(3*ev[:,-1]/sn)*hm*tau);f=1-9*cmu*cmu*sn*sn*tau*tau
                if not np.isfinite(rad).all() or not np.isfinite(f).all() or np.any(rad<=self.guard) or np.any(f<=self.guard):return original()
                coefficient=2*np.sqrt(f)/(1+6*sn*on*tau*tau)*k*tau*tau
                Q=coefficient[:,None,None]*(O@S-S@O);nut=cmu*k*tau
                R=2*k[:,None,None]*np.eye(3)[None,:,:]/3-2*nut[:,None,None]*S+Q
                covariance=np.linalg.eigvalsh(R)[:,0]
                production=-np.sum(R*G,axis=(1,2));stress=-R+np.eye(3)[None,:,:]*(2*k/3-v[outer,4])[:,None,None]
                planar=G.copy();planar[:,0,0]-=div/2;planar[:,2,2]-=div/2
                diagnostic=self.invariants(planar,k,e)
                if diagnostic is None or not all(np.isfinite(a).all() for a in (Q,R,nut,covariance,production,stress)) or np.any(abs(covariance)<=self.guard*np.maximum(1.,abs(k))):return original()
        except (ValueError,FloatingPointError,np.linalg.LinAlgError,OverflowError):return original()
        self.last_arrays=dict(indices=outer,k=k,epsilon=e,raw_gradient=G,raw_divergence=div,strain_eigenvalues=ev,strain_norm=sn,rotation_norm=on,tau=tau,Cmu=C,radicand=rad,R=R,Q=Q,production=production,total_stress=stress,diagnostic_planar_radicand=diagnostic[-1])
        def loc(i):return dict(coordinate=coords[i],geometry=rows[i]['velocity_k']['geometry'],side=sides[i] if sides is not None else 0,sample_index=int(i))
        minimum_inner=None
        if len(inner):
            i=int(inner[np.argmin(v[inner,2])]);minimum_inner=dict(loc(i),value=float(v[i,2]))
        minima={}
        for name,field in [('k',2),('epsilon',3)]:
            i=int(outer[np.argmin(v[outer,field])]);minima[name]=dict(loc(i),value=float(v[i,field]))
        # Re-evaluate every close contender through the ORIGINAL full scalar
        # point function and preserve original observation-order tie-breaking.
        # No batch estimate is itself published as a limiting witness.
        witnesses={}
        def scalar(j):
            i=int(outer[j])
            if i not in witnesses:witnesses[i]=metrics.point(rows[i],s,offsets,model);self.scalar_witness_evaluations+=1
            return i,witnesses[i]
        rad_witness=None;cov_witness=None
        bands={name:np.flatnonzero(a<=float(np.min(a))+self.guard*max(1.,float(np.max(abs(a))))) for name,a in [('rad',rad),('cov',covariance)]}
        for j in bands['rad']:
            i,q=scalar(int(j))
            if q['status']!='admissible':return original()
            if rad_witness is None or q['radicand']<rad_witness['radicand']:rad_witness=dict(q,**loc(i))
        for j in bands['cov']:
            i,q=scalar(int(j))
            if q['status']!='admissible' or 'R' not in q:return original()
            value=float(np.linalg.eigvalsh(q['R'])[0])
            if cov_witness is None or value<cov_witness['value']:cov_witness=dict(loc(i),value=value)
        minima['covariance']=cov_witness
        # Python's original chronological scalar reduction is inexpensive and
        # avoids changing its division/RMS and signed divergence conventions.
        square=0.;maximum=0.
        for value in div:
            value=float(value);square+=value**2;maximum=max(maximum,abs(value))
        return dict(count=len(outer),domain_failure_count=0,first_failures=[],minimum_radicand=rad_witness,minimum_inner_k=minimum_inner,minima=minima,divergence_unweighted_RMS=float(np.sqrt(square/len(outer))),divergence_maximum=maximum)
    def __enter__(self):
        super().__enter__()
        for module in list(sys.modules.values()):
            if module and getattr(module,'__name__','').startswith('ct_accepted') and getattr(module,'inventory',None) is ORIGINAL_INVENTORY:self.saved_inventories.append(module);module.inventory=self.inventory
        return self
    def __exit__(self,*args):
        try:
            for module in self.saved_inventories:module.inventory=ORIGINAL_INVENTORY
        finally:super().__exit__(*args)


def controls(library):
    vk=dict(ids=[0,1],phi=[1.,0.],derivative=[[0.,1.],[0.,0.],[1.,0.]],geometry=7)
    row=dict(velocity_k=vk,epsilon=copy.deepcopy(vk),pressure=copy.deepcopy(vk))
    packet=dict(observation=dict(grid=[[0.,0.,1.]],grid_sides=[0],grid_regions=[1],offsets=list(range(0,13,2)),state=[1.,0.,0.,0.,1.,0.,1.,0.,0.,0.,.25,0.],faces=[]),test_basis=[row])
    o=packet['observation'];base=np.asarray(o['state']);cases=[('admissible',base.copy(),False)]
    for name,field,value in [('negative_k',4,-1.),('zero_epsilon',6,0.),('negative_epsilon',6,-1.)]:
        s=base.copy();s[field]=value;cases.append((name,s,True))
    s=base.copy();s[3]=-1.;cases.append(('pure_rotation',s,True))
    boundary=6.5/(math.sqrt(18)-3)
    for sign in (-1,1):
        s=base.copy();s[0]=0.;s[1]=boundary*(1+sign*1e-12);s[2]=-s[1];cases.append(('radicand_boundary_'+str(sign),s,True))
    results={}
    for name,state,fallback in cases:
        before=ORIGINAL_INVENTORY(packet['test_basis'],o['grid'],state,o['offsets'],'A95',o['grid_sides'])
        with BatchInventoryRows(packet,1,library) as rows:
            after=rows.inventory(packet['test_basis'],o['grid'],state,o['offsets'],'A95',o['grid_sides'])
            require(exact_equal(before,after),'fresh batch negative/forward boundary changed:'+name)
            require(bool(rows.scalar_fallbacks)==fallback,'batch algorithm-selection path not exercised:'+name)
            results[name]=dict(status='pass',original_scalar_fallback=fallback,domain_failure_count=after['domain_failure_count'])
    require(metrics.inventory is ORIGINAL_INVENTORY,'inventory hook remained after control')
    return dict(status='pass',cases=results,thresholds_unchanged=True,no_clipping=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));n.reference=lru_cache(None)(n.reference)
    source=Path(__file__).with_name('ct16_independent_rows.cpp');initial=a.ct14/'static-v1/target/frozen.json';histories=[a.ct14/'pair-v2'/dt/'case.evolution.json' for dt in ('0.00025','0.000125')]
    require(sha(BLAS)==BLAS_SHA,'actual NumPy BLAS bytes changed')
    m=RequestManifest(a.output,dict(schema='CT16FreshBatchedA95QualificationV1',policy=BatchInventoryRows.policy,source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=Path(__file__).resolve().parents[4],text=True).strip(),source_sha256={f:sha(Path(__file__).with_name(f)) for f in ['ct_accepted16_inventory.py','ct_accepted16_blas_audit.py','ct16_independent_rows.cpp','ct_accepted3_metrics.py','ct_accepted2_diagnosis.py','validate_a93_a95_e94e95.py']},initial_sha256=sha(initial),case_sha256=sha(a.ct14/'static-v1/results.json'),actual_history_sha256=[sha(f) for f in histories],reported_bundles_binary64_required=True,internal_tensor_scaled_tolerance=2e-10,physical_advance=False))
    library=a.output/'CT16IndependentRows.so';execution=run(['c++','-std=c++17','-O3','-fPIC','-shared','-ffp-contract=off','-fno-fast-math',str(source),str(BLAS),'-o',str(library)],a.output,'build.log');require(execution['exit_code']==0,'bounded A95 helper build failed')
    ini=n.read(initial);c=n.read(a.ct14/'static-v1/results.json')['case'];weights=n.inventory(c)[2];count=len(weights);states=[ini['observation']['state']]+[n.read(f)['state'] for f in histories]
    for q,region in zip(ini['observation']['grid'][:count],ini['observation']['grid_regions'][:count],strict=True):n.reference(n.col(c),q[2],bool(region))
    out=dict(schema='CT16FreshBatchedA95QualificationV1',status='incomplete',policy=BatchInventoryRows.policy,build_execution=execution,library=runtime_identity(library),controls=controls(library),states=[],physical_advance=False)
    for i,state in enumerate(states):
        cpu=time.process_time();wall=time.monotonic()
        with ExactAuditRows(ini,count,library),TraceIndex(ini,c['H']):
            before=bundle(ini,state,c,weights,count);wc=time.process_time();ww=time.monotonic();again=bundle(ini,state,c,weights,count)
            require(exact_equal(before,again),'original repeated full audit changed');old_warm=dict(CPU=time.process_time()-wc,wall=time.monotonic()-ww)
        old=dict(CPU=time.process_time()-cpu-old_warm['CPU'],wall=time.monotonic()-wall-old_warm['wall'],warm_bundle=old_warm)
        cpu=time.process_time();wall=time.monotonic()
        with BatchInventoryRows(ini,count,library) as rows,TraceIndex(ini,c['H']):
            prepared=dict(CPU=time.process_time()-cpu,wall=time.monotonic()-wall);after=bundle(ini,state,c,weights,count);wc=time.process_time();ww=time.monotonic();again=bundle(ini,state,c,weights,count)
            require(exact_equal(after,again),'fresh batched A95 repeated audit changed');new_warm=dict(CPU=time.process_time()-wc,wall=time.monotonic()-ww)
            require(rows.inventory_calls>0 and rows.last_arrays is not None and rows.scalar_fallbacks==0,'selected A95 batch path not exercised')
            vc=time.process_time();vw=time.monotonic();arrays=rows.last_arrays;s=np.asarray(state);errors={};worst={}
            es=n.col(c).A/(c['H']-c['displacement']);gs=c['ustar']/c['H']
            scales=dict(k=c['K0'],epsilon=es,raw_gradient=gs,raw_divergence=gs,strain_eigenvalues=gs,strain_norm=gs,rotation_norm=gs,tau=c['K0']/es,Cmu=1.,radicand=1.,R=c['K0'],Q=c['K0'],production=es,total_stress=c['K0'])
            for j,index in enumerate(arrays['indices']):
                reference=metrics.point(ini['test_basis'][int(index)],s,ini['observation']['offsets'],'A95');require(reference['status']=='admissible','original selected A95 state failed')
                for name,scale in scales.items():
                    error=float(np.max(abs(np.asarray(arrays[name][j])-reference[name]),initial=0))/scale
                    if error>errors.get(name,-1.):errors[name]=error;worst[name]=int(index)
                diagnostic=reference['diagnostic_planar_trace_removed']
                require(diagnostic['status']=='admissible','original selected planar diagnostic domain')
                error=abs(float(arrays['diagnostic_planar_radicand'][j])-diagnostic['radicand'])
                if error>errors.get('diagnostic_planar_radicand',-1.):errors['diagnostic_planar_radicand']=error;worst['diagnostic_planar_radicand']=int(index)
            verification=dict(CPU=time.process_time()-vc,wall=time.monotonic()-vw,rows=len(arrays['indices']),maximum_scaled_errors=errors,worst_indices=worst)
            write_json(a.output/('all-current-tensor-comparison-'+str(i)+'.json'),verification)
            require(max(errors.values())<=2e-10,'current batched full signed tensor/invariant differs beyond original scaled bound')
        new=dict(CPU=time.process_time()-cpu-new_warm['CPU']-verification['CPU'],wall=time.monotonic()-wall-new_warm['wall']-verification['wall'],warm_bundle=new_warm,immutable_preparation=prepared,plan_bytes=rows.plan_bytes,inventory_calls=rows.inventory_calls,scalar_witness_evaluations=rows.scalar_witness_evaluations,scalar_fallbacks=rows.scalar_fallbacks,qualification_only_tensor_comparison=verification)
        write_json(a.output/('original-bundle-'+str(i)+'.json'),before);write_json(a.output/('candidate-bundle-'+str(i)+'.json'),after)
        require(exact_equal(before,after) and after['metrics']['status']=='pass','fresh A95 bundle/decision/original witness differs')
        out['states'].append(dict(index=i,status='pass',every_reported_observation_decision_witness_exact_equal=True,original=old,batched=new,bundle_sha256=sha(a.output/('candidate-bundle-'+str(i)+'.json'))));write_json(a.output/'progress.json',out);print('CT16 fresh A95 inventory',i,old_warm['CPU'],new_warm['CPU'],flush=True)
    mapped={}
    for line in Path('/proc/self/maps').read_text().splitlines():
        f=line.split(maxsplit=5)
        if len(f)==6 and f[5].startswith('/') and ('/lib' in f[5] or '.so' in f[5]):file=Path(f[5]);require(file.is_file(),'unverifiable loaded inventory library');mapped[str(file)]=sha(file)
    out['Python_runtime']=dict(executable_sha256=sha(Path(sys.executable).resolve()),version=sys.version,numpy_version=np.__version__,actual_mapped_library_sha256=mapped)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
