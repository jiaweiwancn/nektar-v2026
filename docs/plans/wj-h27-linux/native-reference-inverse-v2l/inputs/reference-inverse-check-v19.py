#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent source19 inverse/constructor input controls; no PDE claim."""
import argparse,hashlib,json,math,os,resource,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import numpy as np
from numpy.polynomial import chebyshev as cheb
from scipy.optimize import brentq


def require(condition,reason):
    if not condition:raise ValueError(reason)


def identity(path):
    p=Path(path);b=p.read_bytes();return dict(path=str(p.resolve()),bytes=len(b),sha256=hashlib.sha256(b).hexdigest())


def run(table_path,native_path,manifest_path,out):
    table=json.loads(table_path.read_text());d=json.loads(native_path.read_text());m=json.loads(manifest_path.read_text())
    require(identity(table_path)['sha256']=='b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'==d['table_sha256'],'qualified reference differs')
    require(d['schema']=='WJH27NativeReferenceInverseProbeV1'and d['status']=='RECORDED'and
        d['source_identity']=='56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'and
        d['exact_endpoints_and_return_query']is True,'native source/query guard differs')
    require(m['schema']=='WJH27ReferenceInverseMutationsV1','wrong mutation manifest')
    Kin,H=table['parameters']['K'],table['parameters']['H'];lo,hi=table['K_interval'];c=np.asarray(table['cutoff_table']['coefficients'])
    def cutoff(K):return H*math.exp(float(cheb.chebval((2*K-lo-hi)/(hi-lo),c)))
    derivative=cheb.chebder(c)*2/(hi-lo);tail=float(np.sum(abs(derivative[1:])))
    reserve=32*(len(derivative)+1)*np.finfo(float).eps*float(np.sum(abs(derivative)))
    enclosure=np.array([derivative[0]-tail-reserve,derivative[0]+tail+reserve]);bounds=np.asarray(d['log_cutoff_slope_bounds'])
    require(enclosure[1]<0 and bounds[1]<0 and np.max(abs(enclosure-bounds))<=1e-10*abs(enclosure[0]),'independent slope enclosure differs')
    extrema=cheb.chebroots(cheb.chebder(derivative));points=[-1.,1.]+[float(z.real)for z in extrema if abs(z.imag)<1e-12 and -1<z.real<1]
    slopes=cheb.chebval(np.asarray(points),derivative)
    require(np.all(slopes>=bounds[0])and np.all(slopes<=bounds[1]),'slope enclosure misses polynomial extrema')
    require(len(d['roundtrips'])==len(d['inverses'])==129,'whole-interval query coverage differs')
    rows=[];maximum=dict(cutoff_over_H=0.,K_over_Kin=0.,inverse_residual_over_H=0.,derivative_relative=0.)
    for row in d['roundtrips']:
        K=row['K'];expected=cutoff(K);slope=expected*float(cheb.chebval((2*K-lo-hi)/(hi-lo),derivative))
        errors=dict(cutoff_over_H=abs(row['c']-expected)/H,K_over_Kin=abs(row['inverse']-K)/Kin,
            inverse_residual_over_H=abs(cutoff(row['inverse'])-row['c'])/H,derivative_relative=abs(row['dc']-slope)/abs(slope))
        require(all(v<=1e-10 for v in errors.values()),'native roundtrip/derivative differs')
        for key,value in errors.items():maximum[key]=max(maximum[key],value)
    rr=d['radius_range'];require(rr[0]<rr[1]and abs(rr[0]-cutoff(hi))/H<=1e-10 and abs(rr[1]-cutoff(lo))/H<=1e-10,'inverse interval differs')
    for index,row in enumerate(d['inverses']):
        r=row['r'];K=hi if index==0 else lo if index==128 else brentq(lambda k:cutoff(k)-r,lo,hi,xtol=1e-15*Kin,rtol=4*np.finfo(float).eps)
        errors=dict(K_over_Kin=abs(row['K']-K)/Kin,inverse_residual_over_H=abs(cutoff(row['K'])-r)/H)
        require(all(v<=1e-10 for v in errors.values()),'independent direct inverse differs')
        for key,value in errors.items():maximum[key]=max(maximum[key],value)
        rows.append(dict(radius=r,native_K=row['K'],independent_K=K,errors=errors))
    require(d['inverses'][0]['K']==hi and d['inverses'][-1]['K']==lo,'endpoints are not exact')
    expected={'below inverse interval','above inverse interval','negative radius','zero radius','NaN radius','positive infinity radius','negative infinity radius'}
    require({v['name']for v in d['rejections']}==expected and len(d['rejections'])==7,'native rejection coverage differs')
    mutations={v['name']:v for v in d['mutations']};controls=[]
    require(len(mutations)==len(m['mutations'])==4,'mutation coverage differs')
    for item in m['mutations']:
        p=Path(item['file']['path']);require(identity(p)==item['file'],'synthetic table identity differs')
        native=mutations[item['name']];require('immutable table SHA256 differs'in native['default_identity']['reason']and
            item['expected_reason']in native['constructor_guard']['reason'],'native constructor rejected at the wrong guard')
        mutant=json.loads(p.read_text());poly=np.asarray(mutant['cutoff_table']['coefficients']);dm=cheb.chebder(poly)*2/(hi-lo)
        evidence=dict(name=item['name'],file=item['file'],native=native)
        if item['name']=='two_inverse_roots':
            level=poly[0];shift=poly.copy();shift[0]-=level;roots=cheb.chebroots(shift)
            roots=sorted(float(z.real)for z in roots if abs(z.imag)<1e-12 and -1<=z.real<=1)
            require(len(roots)==2 and roots[1]-roots[0]>1,'synthetic inverse is not actually ambiguous')
            evidence.update(log_radius_level=float(level),two_distinct_normalized_roots=roots,
                physical_K_roots=[(lo+hi+(hi-lo)*x)/2 for x in roots])
        elif item['name']=='constant_cutoff':require(np.all(dm==0),'constant synthetic inverse is not constant')
        elif item['name']=='increasing_cutoff':require(dm[0]>0 and np.all(dm[1:]==0),'synthetic inverse is not increasing')
        elif item['name']=='stale_derivative':
            supplied=np.asarray(mutant['cutoff_table']['logcutoff_dK_coefficients']);require(np.max(abs(dm-supplied))>1e-3,'derivative mutation is not inconsistent')
        else:raise ValueError('unadmitted mutation')
        controls.append(evidence)
    out.update(status='PASS_REFERENCE_INVERSE_GUARDS_ONLY',inverse_guard_gate_pass=True,maximum_scaled_errors=maximum,
        direct_inverse_rows=rows,native_roundtrip_count=129,independent_radius_count=129,non_root_or_invalid_radius_rejections=d['rejections'],
        slope_enclosure=enclosure.tolist(),actual_derivative_extrema=[float(np.min(slopes)),float(np.max(slopes))],mutation_controls=controls,
        scope='Source19 immutable lookup/inverse and synthetic constructor guards only. Synthetic polynomials are not physical references or a global root/nonexistence proof; no PDE, physical-pilot or production-cost closure.')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in('table','native','manifest','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();require(not a.out.exists(),'immutable inverse report exists');began=time.monotonic()
    out=dict(schema='WJH27ReferenceInverseCheckV1',status='FAIL',inverse_guard_gate_pass=False,full_H27_qualification=False,TASK3_READY=False,
        table=identity(a.table),native=identity(a.native),manifest=identity(a.manifest),checker=identity(__file__))
    try:
        require(sys.platform.startswith('linux'),'authoritative Linux check required');run(a.table,a.native,a.manifest,out)
    except Exception as error:out.update(exception=type(error).__name__,reason=str(error))
    out.update(elapsed_seconds=time.monotonic()-began,maximum_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    with a.out.open('x')as f:json.dump(out,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps({k:out[k]for k in('status','reason','elapsed_seconds')if k in out}))
    if not out['inverse_guard_gate_pass']:raise SystemExit(2)


if __name__=='__main__':main()
