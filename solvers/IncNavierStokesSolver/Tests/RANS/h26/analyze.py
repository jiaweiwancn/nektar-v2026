#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Executed numerical gates, disjoint cost screen, and compact evidence index."""
import argparse,hashlib,json,math,os,statistics,struct
from pathlib import Path
from run import THREADS
os.environ.update(THREADS)
import numpy as np

def load(p):return json.loads(p.read_text())
def identity(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def orders(values):return [math.log(a/b,2)for a,b in zip(values,values[1:])]
def state(path):
    data=path.read_bytes();payload=data[:-64];assert hashlib.sha256(payload).hexdigest()==data[-64:].decode();at=payload.index(b'\n')+1;at=payload.index(b'\n',at)+1;step,time,dt=struct.unpack_from('Qdd',payload,at);at+=24;vectors=[]
    for _ in range(4):
        n=struct.unpack_from('Q',payload,at)[0];at+=8;vectors.append(np.frombuffer(payload,dtype='d',count=n,offset=at).copy());at+=8*n
    assert at==len(payload);return dict(step=step,time=time,dt=dt,vectors=vectors)
def numeric(root,version):
    science=root/'science';suffix=f'-v{version}.json';result={};operator=[]
    for c,n in ((0.,0),(0.,1),(.1,0),(.1,1)):
        data=load(science/f'independent-c{c}-n{n}{suffix}');assert data['status']=='PASS_TINY_NATIVE_INDEPENDENT_ASSEMBLY';operator.append(data)
    result['operator_direct_traction_work']=operator
    result['maximum_operator_or_traction_discrepancy']=max(v['relative_error']for d in operator for v in d['checks'].values());assert result['maximum_operator_or_traction_discrepancy']<=1e-10
    result['maximum_direct_solution_error']=max(d['direct_solution_relative_error']for d in operator);assert result['maximum_direct_solution_error']<=1e-8
    result['maximum_independent_momentum_residual']=max(d['independent_momentum_residual']for d in operator);result['maximum_independent_continuity_residual']=max(d['independent_continuity_residual']for d in operator)
    assert max(result['maximum_independent_momentum_residual'],result['maximum_independent_continuity_residual'])<=1e-9
    result['adapter']=load(science/f'adapter{suffix}');result['nonlinear']=load(science/f'nonlinear{suffix}');result['continuous_forcing']=load(science/f'independent-formulas{suffix}')
    assert all(result[k]['status']=='PASS'for k in ('adapter','nonlinear','continuous_forcing'))
    assert result['nonlinear']['continuity_forcing']==0 and result['nonlinear']['projection_continuity']<=1e-9
    temporal=[load(science/f'temporal-poly-{i}{suffix}')for i in range(3)];result['temporal']=dict(rows=[dict(dt=d['dt'],errors=d['L2_errors'])for d in temporal],velocity_orders=orders([d['L2_errors'][0]for d in temporal]),pressure_orders=orders([d['L2_errors'][1]for d in temporal]));assert min(result['temporal']['velocity_orders']+result['temporal']['pressure_orders'])>=1.8
    result['spatial']={};result['WJ_C0_collar_traction']={}
    for curve in (0.,.1):
        rows=[load(science/f'spatial-c{curve}-n{n}{suffix}')for n in (2,4,8)]
        r=dict(rows=[dict(n=n,L2=d['L2_errors'],H1=d['H1_velocity_seminorm_error'],traction=d['traction_errors_L2_and_weak_moments'],strong_divergence=d['strong_divergence'],momentum=d['max_momentum_residual'],continuity=d['max_continuity_residual'])for n,d in zip((2,4,8),rows)])
        r.update(velocity_L2_orders=orders([d['L2_errors'][0]for d in rows]),pressure_L2_orders=orders([d['L2_errors'][1]for d in rows]),velocity_H1_orders=orders([d['H1_velocity_seminorm_error']for d in rows]),traction_L2_orders=orders([d['traction_errors_L2_and_weak_moments'][0]for d in rows]),traction_weak_orders=orders([d['traction_errors_L2_and_weak_moments'][1]for d in rows]))
        assert r['velocity_L2_orders'][-1]>=4 and r['pressure_L2_orders'][-1]>=3 and r['velocity_H1_orders'][-1]>=3
        assert min(r['traction_L2_orders'][-1],r['traction_weak_orders'][-1])>=2.8
        result['spatial'][str(curve)]=r
        wall=[load(science/f'wall-trace-c{curve}-n{n}{suffix}')for n in (2,4,8,16)]
        wall_orders=orders([d['L2_traction_jump_over_k']for d in wall]);weak_orders=orders([d['weak_traction_moments_over_k']for d in wall]);assert wall_orders[-1]>=2.8 and weak_orders[-1]>=2.8
        result['WJ_C0_collar_traction'][str(curve)]=dict(rows=wall,L2_orders=wall_orders,weak_orders=weak_orders)
    whole=load(science/f'restart-whole{suffix}');resumed=load(science/f'restart-resumed{suffix}');a=state(Path(whole['final_checkpoint']));b=state(Path(resumed['final_checkpoint']));assert a['step']==b['step']and a['time']==b['time']
    differences=[float(np.linalg.norm(x-y))for x,y in zip(a['vectors'],b['vectors'])];assert max(differences)==0
    rejection=(root/'operations'/f'restart-reject-height-v{version}'/'command.log').read_text();assert 'checkpoint scientific identity'in rejection
    result['restart']=dict(completed_time=a['time'],step=a['step'],all_four_history_vector_differences=differences,restore_seconds=resumed['restart_seconds'],changed_height_rejected=True)
    tight=load(science/f'tolerance-tight{suffix}');normal=temporal[-1];u=state(Path(normal['final_checkpoint']))['vectors'][0];v=state(Path(tight['final_checkpoint']))['vectors'][0];basis=load(science/f'algebra-c0.0-n0{suffix}');n=basis['velocity_size'];p=basis['pressure_size']
    from oracle import native_sparse
    M=native_sparse(basis['mass']);Q=native_sparse(basis['pressure_mass']);delta=u-v;du=math.sqrt(sum(delta[c*n:(c+1)*n]@M@delta[c*n:(c+1)*n]for c in range(2)));dp=math.sqrt(delta[2*n:]@Q@delta[2*n:]);assert max(du/normal['L2_errors'][0],dp/normal['L2_errors'][1])<1e-4
    result['tolerance_tightening']=dict(tolerances=[1e-11,1e-13],velocity_L2_change=du,pressure_L2_change=dp,changes_over_discretization_errors=[du/normal['L2_errors'][0],dp/normal['L2_errors'][1]])
    constant={}
    for model in ('nu','constant'):
        constant[model]={}
        for scheme in ('S','N'):
            rows=[load(science/f'constant-{model}-dt{i}-{scheme}{suffix}')for i in range(3)];uorders=orders([d['L2_errors'][0]for d in rows]);porders=orders([d['L2_errors'][1]for d in rows]);assert min(uorders)>=1.8 and min(porders)>0
            constant[model][scheme]=dict(rows=[dict(dt=d['dt'],errors=d['L2_errors'])for d in rows],velocity_orders=uorders,pressure_orders=porders)
    result['stock_native_constant_limits']=constant
    for data in temporal+[load(science/f'spatial-c{c}-n{n}{suffix}')for c in (0.,.1)for n in (2,4,8)]:assert max(data['max_momentum_residual'],data['max_continuity_residual'])<=1e-9
    result['status']='PASS';return result
def benchmark(root,version):
    science=root/'science';rows={};fixed_start=[];fixed_final=[];external_fixed=[];checkpoint=[];all_new=[];total_cost=np.zeros(6);total_wall=0;peak=0;workers=0
    config=load(science/f'matrix-v{version}/config.json')
    for name in ('Snu','Smu','Nnu','Nmu','Nspace','Ntime'):
        blocks=[]
        for r in range(3):
            label=f'block-{r}-{name}-v{version}';data=load(science/f'{label}.json');receipt=load(root/'operations'/label/'receipt.json');assert receipt['status']=='PASS';samples=data['ordinary_samples'];assert len(samples)==400 and sum(samples)>3 and data['dt']==.005 and data['time']==2.1
            costs=data['disjoint_seconds'];closure=abs(sum(costs)-sum(samples))/sum(samples);assert closure<1e-8 and min(costs)>-1e-8
            peak=max(peak,receipt['peak_tree_rss_bytes']);workers=max(workers,receipt['peak_host_workers']);record=dict(repetition=r,accepted=420,retained=400,interval=[.1,2.1],mean=statistics.mean(samples),median=statistics.median(samples),minimum=min(samples),maximum=max(samples),stdev=statistics.stdev(samples),p05=float(np.quantile(samples,.05)),p95=float(np.quantile(samples,.95)),warmed_seconds=sum(samples),disjoint_seconds=costs,attribution_fraction=1-closure,receipt=identity(root/'operations'/label/'receipt.json'),raw=identity(science/f'{label}.json'))
            if name.startswith('N'):
                records=data['records'][20:];its=[x['iterations']for x in records];record.update(linear_iterations_mean=statistics.mean(its),linear_iterations_range=[min(its),max(its)],preconditioner_rebuilds=data['preconditioner_rebuilds'],minimum_AWF_feedback=min(x['AWF_feedback']for x in records),minimum_pressure_correction=min(x['pressure_correction']for x in records),modulation_range=data['coefficient_modulation_range'],momentum=data['max_momentum_residual'],continuity=data['max_continuity_residual'])
                assert max(record['momentum'],record['continuity'])<=1e-9
                all_new.append(data)
                if name in ('Nspace','Ntime'):
                    assert record['minimum_AWF_feedback']>1e-6 and record['minimum_pressure_correction']>1e-6
                    fixed_start.append(data['startup_seconds']);fixed_final.append(data['final_output_and_assessment_seconds']);checkpoint+=data['checkpoint_seconds'];external_fixed.append(receipt['seconds']-sum(x['seconds']for x in data['records'])-sum(data['checkpoint_seconds']));total_cost+=costs;total_wall+=sum(samples)
                if name=='Ntime':assert min(data['coefficient_modulation_range'])<=.901 and max(data['coefficient_modulation_range'])>=1.099
            else:record.update(linear_iterations='direct pressure Poisson + two direct component Helmholtz solves per IMEX stage',full_strain_target_difference=data['max_full_strain_target_difference'],weak_continuity=data['max_weak_continuity'])
            record['MMS_errors']=data['L2_errors'];blocks.append(record)
        rows[name]=dict(blocks=blocks,largest_block_mean=max(x['mean']for x in blocks),mean_of_block_means=statistics.mean(x['mean']for x in blocks))
    flow=max(rows[x]['largest_block_mean']for x in ('Nspace','Ntime'));restart=load(science/f'restart-700-v{version}.json');restore=restart['restart_seconds'];assert restart['start_step']==420 and restart['end_step']==440 and restore>0
    transient=max(0.,max(sum(x['seconds']for x in d['records'][:20])-20*flow for d in all_new if d['minimum_AWF_feedback']>0))
    startup=max(fixed_start);final=max(fixed_final);other=max(0.,max(external_fixed)-startup-final);fixed=startup+final+other+restore+transient;one_write=max(checkpoint);count=101005//200;recurring=count*one_write;screen=1.25*(fixed+101005*flow+recurring)
    assert screen<=43200 and workers<=72 and peak<8*1024**3
    return dict(status='PASS',cases=rows,matching=config,cost_screen=dict(N_star=101005,c_flow=flow,C_startup=startup,C_cold_start_transient_excess=transient,C_final_output_and_assessment=final,C_other_fixed_conservative=other,C_actual_700_restart_read_restore=restore,C_fixed=fixed,scheduled_recovery_count=count,conservative_checkpoint_seconds_each=one_write,C_recurring=recurring,compact_summary_count=101005//20,compact_summary_cost='inside ordinary means; charged no second time',C_screen=screen,limit=43200,headroom=43200-screen,rule='1.25*(C_fixed+101005*c_flow+C_recurring); c_flow is the largest warmed Nspace/Ntime block mean'),required_cases_disjoint_fractions=(total_cost/total_wall).tolist(),attribution_fraction=float(total_cost.sum()/total_wall),peak_tree_rss_bytes=peak,peak_actual_host_workers=workers,stock_new_ratios=dict(Nnu_over_Snu=rows['Nnu']['mean_of_block_means']/rows['Snu']['mean_of_block_means'],Nmu_over_Smu=rows['Nmu']['mean_of_block_means']/rows['Smu']['mean_of_block_means']))
def publish(root,version,stage,output):
    assert not output.exists();output.mkdir();numerics=numeric(root,version);summary=dict(stage=stage,tested_source_revision=load(root/f'source-v{version}.json')['revision'],NUMERICS=numerics,FULL_COUPLED_MODEL_AND_PRODUCTION='NOT_QUALIFIED')
    if stage=='C':summary['COST']=benchmark(root,version);summary.update(H26_NUMERICS='PASS',H26_COST='PASS',H26_INTEGRATION_READINESS='READY')
    else:summary.update(H26_NUMERICS='PASS',H26_COST='TARGET_NOT_MET',cost_reason='final selected-binary matrix pending; predecessor v6 is preliminary evidence',H26_INTEGRATION_READINESS='NOT_READY')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Small results and receipts are published; large fields/operators remain
    # at the indexed Linux paths with exact source/runtime/input identities.
    patterns=[f'independent-*-v{version}.json',f'adapter-v{version}.json',f'nonlinear-v{version}.json',f'wall-trace-*-v{version}.json',f'temporal-poly-*-v{version}.json',f'spatial-c*-v{version}.json',f'restart-*-v{version}.json',f'tolerance-tight-v{version}.json',f'constant-*-v{version}.json']
    if stage=='C':patterns+=[f'block-*-v{version}.json',f'pilot-*-v{version}.json']
    files=set()
    for pattern in patterns:files.update((root/'science').glob(pattern))
    files.update((root/f'source-v{version}.json',root/'science'/f'runtime-v{version}.json',root/'accounting.json'))
    files.update((root/'inputs').glob(f'*-v{version}.json'));files.update((root/'operations').glob(f'*-v{version}/request.json'));files.update((root/'operations').glob(f'*-v{version}/receipt.json'))
    for folder in ('fixtures-v2',f'spatial-v{version}',f'matrix-v{version}',f'constant-limit-v{version}'):
        files.update((root/'science'/folder).glob('*.xml'));files.update((root/'science'/folder).glob('*.opt'))
    files.add(root/'science'/f'matrix-v{version}/config.json')
    index=[identity(p)for p in sorted(files)];(output/'linux-index.json').write_text(json.dumps(index,indent=2)+'\n')
    source=load(root/f'source-v{version}.json');runtime=load(root/'science'/f'runtime-v{version}.json');small=dict(source_manifest=identity(root/f'source-v{version}.json'),source_revision=source['revision'],runtime_manifest=identity(root/'science'/f'runtime-v{version}.json'),compiler=runtime['compiler'],configure=runtime['configure'],build=runtime['build'],executables=runtime['executables'],loader=runtime['loader'],library_path=runtime['libraries'],stock_object_count=len(runtime['qualified_stock_objects']))
    (output/'identities.json').write_text(json.dumps(small,indent=2)+'\n')
    if stage=='C':
        samples={}
        for name in ('Snu','Smu','Nnu','Nmu','Nspace','Ntime'):
            for repetition in range(3):
                label=f'block-{repetition}-{name}-v{version}';data=load(root/'science'/f'{label}.json');item=dict(seconds=data['ordinary_samples'])
                if name.startswith('N'):item.update(iterations=[r['iterations']for r in data['records'][20:]],rebuilds=[r['rebuilds']for r in data['records'][20:]],AWF_feedback=[r['AWF_feedback']for r in data['records'][20:]],pressure_correction=[r['pressure_correction']for r in data['records'][20:]])
                samples[label]=item
        (output/'raw-timing-samples.json').write_text(json.dumps(samples,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in summary.items()if k!='NUMERICS'and k!='COST'},indent=2));
    if stage=='C':print(json.dumps(summary['COST']['cost_screen'],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--version',type=int,required=True);p.add_argument('--stage',choices=('B','C'),required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();publish(a.root,a.version,a.stage,a.output)
