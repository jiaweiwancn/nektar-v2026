#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One selected initial shear-image contrast, no predictor factor or evolution."""
import argparse,copy,gc,inspect,os,resource,shutil
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.sparse import block_diag,diags
from scipy.sparse.linalg import splu,LinearOperator,gmres
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted10_diagnose import matrix
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def maximum(v):return float(np.max(abs(np.asarray(v)),initial=0.))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('exports','runtime','contract','ct11','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1'
    exported=sealed(a.exports);target=sealed(a.ct11/'target-v1');baseline_member=target['members']['0.001']
    baseline_report=Path(baseline_member['directory'])/'case.evolution.json'
    require(sha(baseline_report)==baseline_member['native_report_sha256'],'changed baseline initial-state history')
    initial=n.read(baseline_report)['accepted_states'][0]
    require(exact_vector_hash(initial)==baseline_member['actual_initial_sha256'],'foreign baseline actual initialization')
    m=RequestManifest(a.output,dict(schema='CT13InitialShearImageV1',source=sha(Path(__file__)),contract=sha(a.contract),exports=sha(a.exports/'results.json'),CT11=sha(a.ct11/'target-v1/results.json'),baseline_report_sha256=sha(baseline_report),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),observer_policy='NativeAuditReceiptV1'))
    out={};vectors={}
    for label,system in [('baseline','original-initial'),('selected','candidate-alpha1000')]:
        source=Path(exported['systems'][system]['packet']);require(sha(source)==exported['systems'][system]['packet_sha256'],'foreign no-factor packet')
        packet=n.read(source);f=packet['pre_factor_frozen'];nv=f['velocity_size'];rows=np.asarray(f['rows'],int)
        require(f['schema']=='PreFactorFrozenCSRExportV2' and f['whole_command_sparse_factorizations']==0,'unqualified export startup')
        if label=='baseline':require(f['input']==initial and exact_vector_hash(f['input'])==baseline_member['actual_initial_sha256'],'first flow input differs from actual time-zero initialization')
        dest=a.output/label;dest.mkdir();shutil.copy2(source.parent/'case.xml',dest/'case.xml')
        req=dict(operation='replay',identity=f['identity'],descriptor=f['descriptor'],state=f['input'],reference=f['fixed_reference'],fixed_boundary=f['fixed_boundary'],time=0.,grid=[],export_faces_only=True,export_ct10_momentum=True,ct10_export_points=False)
        q,e=probe(a.runtime,dest,req,'loads');require(e['exit_code']==0,'initial native attribution failed:'+str(q.get('failure')))
        mom=q['ct10_momentum'];require(mom['maximum_reconstruction_error']<=2e-10 and mom['maximum_input_lift_difference']==0,'native initial effective-load reconstruction')
        M=matrix(f['free_velocity_mass']);G=matrix(f['G']);D=matrix(f['D']);C=matrix(f['C']);require(not np.any(C.data),'initial image requires actual C0')
        mass=splu(M.tocsc());pc=splu((D@diags(1/M.diagonal())@G).tocsc())
        require(max(mass.L.nnz+mass.U.nnz,pc.L.nnz+pc.U.nnz)<=40000000,'read-only mass/preconditioner fill guard')
        Mex=M.astype(np.longdouble);Gex=G.astype(np.longdouble);Dex=D.astype(np.longdouble)
        def direct(rhs):
            rhs=np.asarray(rhs);x=mass.solve(rhs)
            for k in range(4):
                r=np.asarray(Mex@x.astype(np.longdouble),float)-rhs
                if maximum(r)<=2e-12*maximum(rhs)+1e-25:return x
                require(k<3,'consistent mass true residual');x-=mass.solve(r)
        def action(x):return np.asarray(Dex@direct(np.asarray(Gex@np.asarray(x,dtype=np.longdouble),float)).astype(np.longdouble),float)
        W=LinearOperator(C.shape,matvec=action);P=LinearOperator(C.shape,matvec=pc.solve)
        pm=f['pressure_physical_mass'];masses=[matrix(x['mass']) for x in pm];constants=[np.asarray(x['constant']) for x in pm]
        require(len(rows[nv:])==sum(len(x) for x in constants) and np.array_equal(rows[nv:],np.arange(f['offsets'][4],f['offsets'][6])),'absolute physical pressure rows/gauge changed')
        volumes=[float(c@(A@c)) for c,A in zip(constants,masses,strict=True)]
        expected=[55*.1875,55*(20-.1875)];require(max(abs(x-y) for x,y in zip(volumes,expected,strict=True))<2e-10,'pressure mass/common physical volume')
        Mp=block_diag(masses,format='csr');us2=.07426793792441279**2
        def measure(x):
            at=0;stats={}
            for name,A,c,v in zip(('inner','outer'),masses,constants,volumes,strict=True):
                z=x[at:at+len(c)];mean=float(c@(A@z))/v;square=float(z@(A@z))/v
                stats[name]=dict(weighted_RMS=np.sqrt(max(0.,square))/us2,signed_mean=mean/us2,variation_RMS=np.sqrt(max(0.,square-mean*mean))/us2);at+=len(c)
            stats['whole']=dict(weighted_RMS=float(np.sqrt(max(0.,x@(Mp@x)/sum(volumes))))/us2,signed_mean=sum(stats[name]['signed_mean']*v for name,v in zip(('inner','outer'),volumes,strict=True))/sum(volumes))
            return stats
        checks=[]
        def solve(rhs,tol):
            mag=float(np.linalg.norm(rhs))
            if mag==0:return np.zeros(C.shape[0])
            ticks=[];kw={'rtol' if 'rtol' in inspect.signature(gmres).parameters else 'tol':tol}
            x,info=gmres(W,rhs/mag,M=P,**kw,atol=1e-27/mag,restart=80,maxiter=12,callback=ticks.append,callback_type='pr_norm');x*=mag
            r=action(x)-rhs;checks.append(dict(relative_request=tol,info=int(info),iterations=len(ticks),true_residual_maximum=maximum(r),rhs_maximum=maximum(rhs)))
            require(maximum(r)<=2e-10*maximum(rhs)+1e-25,'initial image unscaled true action residual')
            return x
        images={};measurements={}
        for kind in ('inner_shear','total_representation'):
            keys=[k for k in mom['actual'] if kind=='total_representation' or k in ('inner/volume/linear_shear','inner/volume/molecular_shear')]
            require(len(keys)>=2,'missing declared initial shear components')
            load=sum((np.asarray(mom['actual'][k])-mom['analytic'][k] for k in keys),np.zeros(mom['momentum_size']))[rows[:nv]]
            rhs=-np.asarray(Dex@direct(load).astype(np.longdouble),float);v=solve(rhs,2e-12);tight=solve(rhs,2e-14)
            # Modified_A P2 tensor basis has nine local functions, each <=1.
            repeat_upper=10*maximum(v-tight)/us2;require(repeat_upper<1e-8,'initial image physical numerical uncertainty')
            images[kind]=v;measurements[kind]=dict(pressure=measure(v),repeat_physical_max_upper=repeat_upper,free_load_maximum=maximum(load),components=keys)
        den=float(images['total_representation']@(Mp@images['total_representation']));alignment=float(images['inner_shear']@(Mp@images['total_representation']))/den
        np.savez(a.output/(label+'-images.npz'),**images)
        out[label]=dict(native_execution=e,native_packet_sha256=sha(dest/'loads.json'),no_factor_packet_sha256=sha(source),state_sha256=exact_vector_hash(f['input']),row_reconstruction=mom['maximum_reconstruction_error'],pressure_volumes=volumes,measurements=measurements,signed_shear_alignment=alignment,solve_checks=checks,pressure_image_file_sha256=sha(a.output/(label+'-images.npz')),scope='all native pressure modes,exact physical pressure mass; no mean removed,coefficient-position comparison,predictor factor or physical step')
        vectors[label]=measurements['inner_shear']['pressure']['whole']['weighted_RMS'];write_json(a.output/'progress.json',out)
        del W,P,action,direct,Mex,Gex,Dex,mass,pc,q;gc.collect()
    ratio=vectors['selected']/vectors['baseline'];uncertainty=sum(r['measurements']['inner_shear']['repeat_physical_max_upper'] for r in out.values())
    result=dict(schema='CT13InitialShearImageV2',status='diagnostic_pass',cases=out,selected_to_baseline_shear_RMS_ratio=ratio,initial_image_magnitude_decrease=ratio<1,resolved_initial_image_reduction=vectors['baseline']-vectors['selected']>max(1e-8,10*uncertainty),comparison_uncertainty_floor=max(1e-8,10*uncertainty),interpretation='static pressure-image sensitivity only; not evolved pressure uncertainty or spatial convergence; new fine modes retained; physical continuation remains independently gated')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT13 initial shear image ratio',ratio,flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
