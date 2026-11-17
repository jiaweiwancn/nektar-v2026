from pathlib import Path
import os,sys,json,math,shutil,subprocess,time
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/cutoff-panels-v19';t=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27';sys.path.insert(0,str(t));sys.path.insert(0,str(r/'tools'))
from resources import atomic,identity
import implicit_algorithm_check as fields
import physical_cache_history_check as history
read=history.read
assert read(o/'operator-independent.json')['operator_gate_pass']
assert read(o/'regression-independent.json')['foundation_gate_pass']
assert read(o/'q16-panels-independent.json')['quadrature_foundation_gate_pass']
assert read(o/'radial-zero-independent.json')['four_binary_histories_exact']
native=read(o/'reference-r2.result.json');a=native['state'];source=read(r/'source-v19.json');assert a['source_identity']==source['expected_source_identity']
traces={}
for name in ('r2','h2','coarse'):
 p=r/'science'/('trace-geometry-'+name+'-v18.json');j=read(p);assert j['horizontal_trace_gate_pass']and j['native_face_gate_pass']
 traces[name]=dict(file=identity(p),**{key:value for key,value in j.items()if key!='rows'})
exact={}
for name,oldname in [('poly','poly-cache'),('hist','hist-cache'),('physical','physical-cache')]:
 old=fields.load(r/'science/coupled-fixtures-v15'/(oldname+'.result.json'),'CoupledJFNKV1');new=fields.load(r/'science/htrace-regression-v18'/(name+'.result.json'),'CoupledJFNKV1')
 hashes={key:dict(old=history.transient.state_hash(old['states'][key]),new=history.transient.state_hash(new['states'][key]))for key in fields.STATE_NAMES+fields.CACHE_NAMES}
 exact[name]=dict(histories=hashes,all_seven_byte_exact=all(v['old']==v['new']for v in hashes.values()))
 assert exact[name]['all_seven_byte_exact']
reports={}
for order in (32,48):
 p=o/('independent'+str(order)+'.json');j=read(p);assert 'exception'not in j and j['core_steady_equations_status']=='MET'and j['boundary_checks']['status']=='MET'
 reports[str(order)]=dict(file=identity(p),current_equations=j['core_steady_equations_status'],boundaries=j['boundary_checks']['status'],accounting=j['core']['accounting_identity_status'],
  residuals=j['core']['residuals'],physical_development=j['physical_development_status'],physical_envelopes=j['core']['physical_envelope'],
  physical_station_envelopes=j['profiles']['norms'],physical_ground_stations=j['profiles']['ground_stations'],original_owner_flux_jumps=j['owner_interfaces']['norms'])
j=read(o/'independent48.json');maxima={}
for row in j['profiles']['rows']:
 if row['outer']:continue
 delta=(row['values']['k']-row['retained_inlet']['k'])/.008809015145515192;key=str(row['x'])
 if key not in maxima or abs(delta)>abs(maxima[key]['scaled_k_departure']):
  maxima[key]=dict(x=row['x'],r=row['r'],z_mm=row['z_mm'],owner=row['owner'],bounds=row['bounds'],k=row['values']['k'],retained_inlet_k=row['retained_inlet']['k'],scaled_k_departure=delta,cutoff=row['values']['cutoff'],K=row['values']['K'])
reports['48']['sampled_inner_k_extrema_by_station']=maxima
account=read(r/'accounting.json');labels=[key for key in account['entries']if key.startswith(('trace-geometry-','radial-','htrace-','cutoff-','materialize-v18','materialize-v19'))]
entries={key:account['entries'][key]for key in labels};cases={}
for name in ('poly','hist','physical','reference-r2','reference-r2-reload'):
 packet=read(o/(name+'.result.json'));s=packet['state'];rows=[hr for line in(o/(name+'.steps.jsonl')).read_text().splitlines()if line.strip()for hr in[history.hr.native_json(line)]]
 cases[name]={key:packet[key]for key in ('construction_seconds','initialisation_seconds','complete_advancing_seconds','final_export_seconds','complete_step_seconds','scientific_frame_seconds','scientific_index_seconds','recovery_seconds')}
 cases[name].update(source_identity=s['source_identity'],signature=s['signature'],end_step=s['step'],end_time=s['time'],dt=s['dt'],newton_iterations=[row['nonlinear_iterations']for row in rows],
  krylov_iterations=[row['coupled_newton']['krylov_iterations']for row in rows],steady_solver=({key:value for key,value in packet['steady_solver'].items()if key not in ('best_unknowns','best_raw_residual')}if packet.get('steady_solver')else None),minimum_k=s['minimum_k'],minimum_epsilon=s['minimum_epsilon'],
  minimum_covariance_over_k=s['minimum_covariance_over_k'],minimum_principal_over_mu=s['minimum_principal_over_mu'])
v=cases['physical']['complete_step_seconds'];means=[math.fsum(v[i:i+3])/3 for i in (0,3,6,9)];N=math.ceil(505.02466521278404/.005)
conditional=dict(dt=.005,steps=N,initial_block_mean=means[0],warm_block_means=means[1:],largest_warm_mean=max(means[1:]),
 advancement_only_seconds=N*max(means[1:]),one_point_two_five_advancement_only_seconds=1.25*N*max(means[1:]),four_trajectory_advancement_only_seconds=4*N*max(means[1:]),
 complete_fixed_costs='NOT_QUALIFIED',monitoring_and_statistics_cost='NOT_QUALIFIED',long_state_and_perturbation_cost='NOT_QUALIFIED',selected_mesh_p_q_dt='NOT_QUALIFIED',production_forecast_qualified=False,
 scope='Short700-cell work-class extrapolation only; excluded components are unknown, not zero. All actual operations remain in the aggregate account, and overlapping/subset timings are not added twice.')
primary='/media/jay/新加卷1/nektar-v2026';head=subprocess.check_output(['git','-C',primary,'rev-parse','HEAD'],text=True).strip();status=subprocess.check_output(['git','-C',primary,'status','--porcelain'],text=True)
assert head=='f4b70baec10016dd098005bc5a5507b87e4a57bf'and status==' M docs/developer-guide/src/library\n'
def bytes_in(folder):return sum(p.stat().st_size for p in folder.rglob('*')if p.is_file())
source_bytes=sum(bytes_in(p)for p in r.iterdir()if p.is_dir()and p.name.startswith(('source-v','build-v','runtime-v')));science_bytes=bytes_in(r/'science');free=shutil.disk_usage(r).free
assert source_bytes<8*1024**3 and science_bytes<8*1024**3 and free>15*1024**3
resources=dict(primary=dict(path=primary,head=head,status=status,preserved=True),source_build_runtime_bytes=source_bytes,science_bytes=science_bytes,free_bytes=free,
 aggregate_closed_prior_and_current_seconds=sum(v.get('seconds',0)for v in account['entries'].values()),account_closed=account['closed'],campaign_ceiling=43200,evidence_protected=3600,
 workers_ceiling=72,max_phase_actual_host_workers=max(v.get('peak_host_workers',0)for v in entries.values()),max_phase_tree_rss_bytes=max(v.get('peak_tree_rss_bytes',0)for v in entries.values()),ranks=1,all_six_thread_limits_one=True)
result=dict(schema='WJH27CutoffPanelCheckpointV1',status='PASS_CUTOFF_INTEGRATION_FOUNDATION_ONLY',full_H27_qualification=False,TASK3_READY=False,
 active_model='FiniteNuIsotropicWallFullWJCollarV2',unapproved_v4_implemented=False,source=identity(r/'source-v19.json'),runtime=identity(r/'science/runtime-v19.json'),source_identity=a['source_identity'],
 native=identity(o/'reference-r2.result.json'),operator=identity(o/'operator-independent.json'),regression=identity(o/'regression-independent.json'),radial_comparison=identity(o/'q16-panels-independent.json'),
 cold_zero=identity(o/'radial-zero-independent.json'),source15_to18_seven_role_hashes=exact,trace_checks=traces,reports=reports,cases=cases,conditional_cost=conditional,resources=resources,
 operations=entries,mesh=dict(nx=28,nr=200,cells=5600,p_velocity_k_epsilon=3,p_pressure=2,quadratureX=8,quadratureR=12,scalar_quadrature='CubicCutoffFacePanelsV1'),physical_time_advanced_for_steady=0,
 prior_non_nested_comparisons=dict(h2_to_r2=identity(r/'science/htrace-regression-v18/comparison-h2-r2.json'),x2_vs_r2=identity(r/'science/htrace-regression-v18/comparison-x2-r2.json'),
  scope='Source18 Gauss8 candidate comparison retained with its failed current-k equation. These anisotropic diagnostic differences are not selected-source19 h/p/q accuracy qualification.'),
 physical_limits=dict(inner_k_RMS=.05,inner_k_MAX=.10),numerical_k_limits=dict(RMS=.001,MAX=.001),no_continuum_limit_or_uniqueness_claim=True,
 global_gates=dict(H27_MODEL_AND_REFERENCE='NOT_MET',H27_COUPLED_NUMERICS='NOT_MET',H27_PHYSICAL_PILOTS_AND_ACCURACY='NOT_MET',H27_RESTART_MONITORING_STORAGE='NOT_MET',H27_COMPLETE_COST='TARGET_NOT_MET',TASK3_READY='NO'),
 checker=identity(Path(__file__)))
try:
 os.environ['MPLCONFIGDIR']=str(o/'matplotlib-cache')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,2,figsize=(9,4.8),sharey=True);plotted=[]
 datasets=[('700 (14x50)',r/'science/physical-steady-independent-v11.json','#1f77b4','--'),('2800 (28x100)',r/'science/physical-reference-h2-v1/independent-h2.json','#ff7f0e','-.'),
  ('5600 x (56x100)',r/'science/physical-reference-x2-v1/independent-x2.json','#2ca02c',':'),('5600 r, panels (28x200)',o/'independent48.json','#d62728','-')]
 for label,path,color,style in datasets:
  packet=read(path);rows=[row for row in packet['profiles']['rows']if row['x']==0 and not row['outer']];assert rows
  groups={}
  for row in rows:groups.setdefault(row['owner'],[]).append(row)
  first=True
  for owner,values in sorted(groups.items()):
   values=sorted(values,key=lambda row:row['r']);height=[row['z_mm']for row in values];kval=[row['values']['k']/.008809015145515192 for row in values];departure=[(row['values']['k']-row['retained_inlet']['k'])/.008809015145515192 for row in values]
   axes[0].plot(kval,height,color=color,ls=style,lw=1.6,label=label if first else'_nolegend_');axes[1].plot(departure,height,color=color,ls=style,lw=1.6)
   if label.startswith('700'):axes[0].plot([row['retained_inlet']['k']/.008809015145515192 for row in values],height,color='0.35',lw=1.,label='Retained inlet'if first else'_nolegend_')
   first=False
  plotted.append(dict(label=label,source=identity(path),points=len(rows),owners=len(groups)))
 axes[1].axvspan(-.10,.10,color='0.5',alpha=.16);axes[1].axvline(0,color='0.5',lw=.7)
 for ax in axes:ax.set_ylim(3.5,9.5);ax.grid(True,alpha=.18)
 axes[0].set_xlabel('k / Kin');axes[0].set_ylabel('Physical elevation z (mm)');axes[1].set_xlabel('(k - retained inlet k) / Kin')
 fig.suptitle('Approved v2 inner-k profiles at x = 0',y=.96);handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
 fig.text(.5,.105,'Shading: fixed pointwise envelope +/-0.10 Kin; RMS criterion assessed separately.\nOlder meshes are diagnostic: no selected-source19 mesh-convergence claim.',ha='center',fontsize=8)
 fig.subplots_adjust(left=.09,right=.98,bottom=.25,top=.88,wspace=.12);png=o/'v2-inner-k-x0.png';svg=o/'v2-inner-k-x0.svg';fig.savefig(png,dpi=180);fig.savefig(svg,metadata={'Date':None});plt.close(fig)
 result['visualization']=dict(png=identity(png),svg=identity(svg),inputs=plotted,owners_preserved=True,scope='Original-owner profile samples; no averaging or new acceptance test')
except (ImportError,FileNotFoundError)as error:result['visualization_unavailable']=str(error)
atomic(o/'assessment.json',result)
print(json.dumps(dict(status=result['status'],physical_development=reports['48']['physical_development'],resources=resources,conditional_cost=conditional)))
