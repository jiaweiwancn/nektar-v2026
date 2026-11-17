#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare native moving-cutoff rows with the frozen independent witness."""
import argparse,json,math,os,sys,time
from pathlib import Path
sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import physical_oracle as po
import physical_cache_history_check as history
import history_recovery_check as hr
require,read = po.require,history.read


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('manifest','table','out'):parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args();began=time.monotonic();require(not args.out.exists(),'operator assessment output is immutable')
    out=dict(schema='WJH27CutoffPanelOperatorCheckV1',status='FAIL',operator_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,
        manifest=hr.identity(args.manifest),table=hr.identity(args.table),checker=hr.identity(__file__),workers=1)
    try:
        require(sys.platform.startswith('linux'),'authoritative Linux operator assessment required')
        manifest=read(args.manifest);root=hr.resolve(args.manifest.resolve().parent,manifest['root'])
        require(manifest['schema']=='WJH27CutoffPanelOperatorV1'and hr.identity(args.table)['sha256']==po.TABLE_SHA,'wrong operator/table binding')
        paths={key:history.input_path(manifest[key],root)[0]for key in ('native','witness','request')}
        native,witness,request=(read(paths[key])for key in ('native','witness','request'))
        require(native['schema']=='WJH27NativeCutoffPanelRowsV1'and native['source_identity']==manifest['source_identity']and
            native['reference_identity']==po.TABLE_SHA and native['policy']=='CubicCutoffFacePanelsV1','native operator identity differs')
        require(request['source']['sha256']==witness['native']['sha256']and witness['prospective_operator_limit']==1e-10 and
            witness['prospective_reference_limit']==1e-11,'frozen independent operator witness changed')
        references={row['owner']:row for row in witness['owners']};require({row['owner']for row in native['cells']}==set(references),'original signed owner inventory differs')
        result={}
        for cell in native['cells']:
            reference=references[cell['owner']];xl,xr,rl,rr=reference['bounds'];area=(xr-xl)*(rr-rl)
            require(abs(cell['area']-area)<=1e-12*area and cell['points']>0,'native positive area does not close')
            require(max(reference['reference_errors'].values())<=1e-11,'independent reference remains unresolved')
            errors={}
            for key in ('k_volume','Ew_load'):
                a,b=np.asarray(cell[key]),np.asarray(reference['reference_rows'][key]);require(a.shape==b.shape==(16,)and np.all(np.isfinite(a)),'native row inventory differs')
                errors[key]=float(np.max(abs(a-b)))/reference['scale']
                require(errors[key]<=1e-10,'native crossing row exceeds fixed1e-10 action limit')
            result[str(cell['owner'])]=dict(errors=errors,points=cell['points'],area=cell['area'],original_bounds=reference['bounds'],scale=reference['scale'])
        ref=po.Reference(read(args.table));derivative=np.polynomial.chebyshev.chebder(ref.coefficients)*2/(ref.hi-ref.lo);maximum=0.;inverse=0.
        require(len(native['inverse_cases'])==21 and native['out_of_range_inverse_rejected'] is True and native['log_slope_bounds'][1]<0,'bounded inverse controls failed')
        for row in native['inverse_cases']:
            K=row['K'];c=float(ref.cutoff(np.asarray([K]))[0]);d=c*np.polynomial.chebyshev.chebval((2*K-ref.lo-ref.hi)/(ref.hi-ref.lo),derivative)
            inverse=max(inverse,abs(c-row['radius'])/ref.H);maximum=max(maximum,ref.Kin*abs(d-row['derivative'])/ref.H)
            require(d<0 and ref.lo<=K<=ref.hi,'inverse left its certified decreasing branch')
        require(inverse<=1e-10 and maximum<=1e-8,'inverse residual/derivative agreement failed')
        out.update(status='PASS_NATIVE_CUTOFF_PANEL_OPERATOR_ONLY',operator_gate_pass=True,cells=result,inverse_residual_over_H_MAX=inverse,
            Kin_times_derivative_difference_over_H_MAX=maximum,native_slope_bounds=native['log_slope_bounds'],independent_slope_bounds=ref.log_cutoff_slope_bounds,
            sources={key:hr.identity(path)for key,path in paths.items()},scope='Actual local conservative k/Ew rows on two retained original owners; not a physical trajectory or global uniqueness claim')
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out['elapsed_seconds']=time.monotonic()-began
    with args.out.open('x')as stream:json.dump(out,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({key:out[key]for key in ('status','reason','elapsed_seconds')if key in out}))
    if not out['operator_gate_pass']:raise SystemExit(2)


if __name__=='__main__':main()
