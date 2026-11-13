#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Direct element-to-CSR physical consistent-mass diagnostic; default unchanged."""
import argparse,os,time,resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from scipy.sparse import coo_matrix,block_diag
from scipy.sparse.linalg import splu
from scipy.linalg import lu_factor,lu_solve
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import sha,require
from continue_cls_epawf_pimple import atomic

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--native',type=Path,required=True);p.add_argument('--context',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False);d=load(a.native);f=d['context'];root=Path(d['array_directory'])
    def arr(name):v=f['files'][name];require(sha(root/v['file'])==v['sha256'],'native mass input checksum');return np.memmap(root/v['file'],dtype=v['dtype'],mode='r',shape=tuple(v['shape']))
    local=arr('element_mass');n=f['velocity_field_size'];free=np.array(f['free_rows'],int);began=time.perf_counter();rows=[];cols=[];values=[]
    for record,element in zip(local,f['elements']):
        count=element['count'];ids=record[:count].astype(int);matrix=record[count:].reshape(count,count);rows.extend(np.repeat(ids,count));cols.extend(np.tile(ids,count));values.extend(matrix.ravel())
    field=coo_matrix((values,(rows,cols)),shape=(n,n)).tocsr();field.sum_duplicates();field.eliminate_zeros();full=block_diag([field,field],format='csr');mass=full[free][:,free].tocsc();assembly=time.perf_counter()-began;pre_dense_peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    factor_times=[]
    for _ in range(3):t=time.perf_counter();factor=splu(mass);factor_times.append(time.perf_counter()-t)
    gradients=np.array(arr('gradients'));Q=np.array(arr('continuity_Q'));rhs=np.vstack([np.sin(np.arange(len(free))*.17),gradients])
    solve_times=[]
    for _ in range(3):t=time.perf_counter();sol=factor.solve(rhs.T);solve_times.append(time.perf_counter()-t)
    inverse_residual=float(np.max(abs(mass@sol-rhs.T))/max(1,float(np.max(abs(rhs)))))
    # Dense reference is loaded AFTER direct CSR construction/factor/action.
    dense=np.array(arr('mass_free'));field_dense=np.array(arr('mass_field'));t=time.perf_counter();dense_factor=lu_factor(dense);dense_seconds=time.perf_counter()-t;t=time.perf_counter();reference=lu_solve(dense_factor,rhs.T);dense_solve=time.perf_counter()-t
    matrix_error=float(np.max(abs(field.toarray()-field_dense))/np.max(abs(field_dense)));row_scaled_error=float(np.max(np.max(abs(mass.toarray()-dense),axis=1)/np.maximum(1e-300,np.max(abs(dense),axis=1))));inverse_per_rhs=(np.max(abs(sol-reference),axis=0)/np.maximum(1,np.max(abs(reference),axis=0))).tolist();free_error=float(np.max(abs(mass.toarray()-dense))/np.max(abs(dense)));inverse_error=float(np.max(abs(sol-reference))/max(1,float(np.max(abs(reference)))))
    dual=np.einsum('ij,ji->i',rhs,sol);ref_dual=np.einsum('ij,ji->i',rhs,reference);dual_error=float(np.max(abs(dual-ref_dual))/max(1,float(np.max(abs(ref_dual)))))
    iq=factor.solve(Q.T);S=Q@iq;mu=lu_solve(lu_factor(S),Q@sol[:,1:]);response=sol[:,1:]-iq@mu;native=arr('global_responses');response_error=float(np.max(abs(response.T-native))/max(1,float(np.max(abs(native)))));response_per_rhs=(np.max(abs(response.T-native),axis=1)/np.maximum(1,np.max(abs(native),axis=1))).tolist();Q_raw=float(np.max(abs(Q@response)));response_norm=np.sqrt(np.maximum(0,np.einsum('ni,ij,nj->n',response.T,dense,response.T)));Q_error=float(np.max(abs(Q@response)/np.maximum(1,response_norm)))
    state=np.cos(np.arange(2*n)*.13);fixed=np.setdiff1d(np.arange(2*n),free);load_full=full@state;eliminated=load_full[free]-full[free][:,fixed]@state[fixed];elimination_error=float(np.max(abs(eliminated-mass@state[free])))
    require(max(matrix_error,free_error,row_scaled_error)<1e-11 and max(inverse_error,dual_error,response_error,max(inverse_per_rhs),max(response_per_rhs))<1e-9 and Q_error<1e-9 and inverse_residual<1e-9,'sparse mass equivalence gate failed')
    def storage(A):return A.data.nbytes+A.indices.nbytes+A.indptr.nbytes
    out=dict(version='DirectLocalPhysicalMassCSR-V1',passed=True,production_default_changed=False,physical_field_size=n,free_velocity_rows=len(free),element_count=len(local),nnz=mass.nnz,field_nnz=field.nnz,factor_L_nnz=factor.L.nnz,factor_U_nnz=factor.U.nnz,factor_fill_ratio=(factor.L.nnz+factor.U.nnz)/mass.nnz,CSR_bytes=storage(mass),factor_bytes=storage(factor.L)+storage(factor.U),dense_free_bytes=dense.nbytes,direct_assembly_seconds=assembly,sparse_factor_seconds=factor_times,sparse_multi_rhs_solve_seconds=solve_times,dense_factor_seconds=dense_seconds,dense_solve_seconds=dense_solve,matrix_relative_max=matrix_error,row_scaled_matrix_error=row_scaled_error,inverse_per_rhs=inverse_per_rhs,native_response_per_rhs=response_per_rhs,free_matrix_relative_max=free_error,inverse_relative_max=inverse_error,inverse_residual=inverse_residual,dual_relative_error=dual_error,native_constrained_response_relative_error=response_error,Q_response_mass_unit_max=Q_error,Q_response_unscaled_max=Q_raw,essential_elimination_max_error=elimination_error,peak_before_dense_reference_KiB=pre_dense_peak,peak_total_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,source='native local signed-CG physical element masses; sparse factor never receives the exported dense reference',limits='diagnostic only; pressure Schur dense; production 3000-slot guard unchanged')
    atomic(a.output/'sparse-mass.json',out)
if __name__=='__main__':main()
