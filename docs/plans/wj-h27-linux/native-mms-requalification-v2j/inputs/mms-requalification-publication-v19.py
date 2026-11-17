from pathlib import Path
import os,sys,json,math,shutil,tarfile,csv,subprocess
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/mms-requalification-v19';t=r/'tools-mms-requalification-v19';out=r/'evidence/native-mms-requalification-v2j';out.mkdir(parents=True);raw=[]
def read(p):return json.loads(p.read_text(),parse_int=lambda token:-0.0 if token=='-0'else int(token))
temporal=read(o/'temporal/independent.json');spatial=read(o/'spatial/independent.json');cross=read(o/'crossing/independent.json')
assert temporal['study_gate_pass']is True and spatial['spatial_MMS_gate_pass']is True and cross['foundation_gate_pass']is True
cases=read(o/'cases.json');source=cases['source_identity'];assert source=='56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
assert temporal['source_identity']==source and all(v['source_identity']==source for v in spatial['runs'].values())
account=read(r/'accounting.json');labels=[key for key in account['entries']if key.startswith(('mms-requal-','mms-cross-'))];entries={key:account['entries'][key]for key in labels};case_stats={}
for name,spec in cases['cases'].items():
 p=Path(spec['session']['path']).with_suffix('.result.json');packet=read(p);data=packet['state'];assert data['source_identity']==source and data['implicit_solver']=='CoupledJFNKV1'and data['coupled_newton']['preconditioner_policy']=='CurrentDataDriftV1'
 prefix=str(p)[:-len('.result.json')];logs=[readline for line in Path(prefix+'.steps.jsonl').read_text().splitlines()for readline in[json.loads(line)]]
 assert all(v['source_identity']==source and v['implicit_solver']=='CoupledJFNKV1'for v in logs)
 key='mms-cross-'+name.split('/')[1]+'-native-v19'if name.startswith('crossing/')else'mms-requal-'+name.split('/')[1]+'-native-v19'
 steady=packet.get('steady_solver');steady={k:v for k,v in steady.items()if k not in('best_unknowns','best_raw_residual')}if steady else None
 case_stats[name]=dict(result=identity(p),session=spec['session'],cells=len(data['elements']),dt=data['dt'],end_step=data['step'],end_time=data['time'],quadrature=[data['quadratureX'],data['quadratureR']],
  whole_native_seconds=entries[key]['seconds'],construction_seconds=packet['construction_seconds'],initialisation_seconds=packet['initialisation_seconds'],
  complete_advancing_seconds=packet['complete_advancing_seconds'],complete_step_seconds=packet['complete_step_seconds'],final_export_seconds=packet['final_export_seconds'],
  newton_iterations=[row['nonlinear_iterations']for row in logs],krylov_iterations=[row['coupled_newton']['krylov_iterations']for row in logs],
  refreshes=[row['coupled_newton']['refreshes']for row in logs],reuses=[row['coupled_newton']['reuses']for row in logs],steady_solver=steady,
  native_minima={key:data[key]for key in('minimum_k','minimum_epsilon','minimum_covariance_over_k','minimum_principal_over_mu')})
def copy(p):
 target=out/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(identity(p))
drop={'coefficients','elements','input_history','current_advection','previous_advection','free_velocity_ids','k_free_ids','epsilon_free_ids','essential_flow_lift','essential_k_lift','essential_epsilon_lift','mms_bulk_body_rhs','mms_sampled_body_rhs'}
def compact(v,parent=''):
 if isinstance(v,dict):return{k:compact(value,k)for k,value in v.items()if not(v.get('schema')=='WJH27NativeCoupledStateV2'and k in drop)and k not in('best_unknowns','best_raw_residual')}
 if isinstance(v,list):return[compact(value,parent)for value in v]
 return v
for p in sorted(o.rglob('*')):
 if not p.is_file():continue
 raw.append(identity(p))
 if p.suffix in('.bin','.fld'):continue
 if p.suffix in('.xml','.opt','.jsonl'):copy(p)
 elif p.suffix=='.json':
  if p.name.endswith(('.result.json','.steady.json')):
   value=compact(read(p));value['preserved_raw_artifact']=identity(p);value['publication_scope']='Compact native evidence; full signed coefficients/history remain on Linux.';target=out/'compact'/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);atomic(target,value)
  else:copy(p)
for p in sorted((r/'inputs').iterdir()):
 if p.is_file()and('mms-requalification'in p.name or p.name.startswith(('mms-v19-','mms-crossing-'))):copy(p)
for p in sorted(t.rglob('*')):
 if p.is_file()and(p.name in('coupled_oracle.py','mms_study.py','mms_cutoff_crossing_check.py','rough-fig3.csv','admission.md','overlay.json')):copy(p)
copy(r/'launches/mms-requalification-v19.log');copy(r/'launches/mms-crossing-v19b.log');copy(r/'science/runtime-v19.json')
for label in labels:
 for name in('request.json','receipt.json','command.log'):
  p=r/'operations'/label/name
  if p.exists():copy(p)
shutil.copy2(r/'accounting.json',out/'accounting-open-snapshot.json')
def bytes_in(folder):return sum(p.stat().st_size for p in folder.rglob('*')if p.is_file())
source_bytes=sum(bytes_in(p)for p in r.iterdir()if p.is_dir()and p.name.startswith(('source-v','build-v','runtime-v')));science_bytes=bytes_in(r/'science');free=shutil.disk_usage(r).free
assert source_bytes<8*1024**3 and science_bytes<8*1024**3 and free>15*1024**3
resource=dict(source_build_runtime_bytes=source_bytes,science_bytes=science_bytes,free_bytes=free,account_closed=account['closed'],aggregate_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in account['entries'].values()),
 phase_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in entries.values()),max_phase_host_workers=max(v.get('peak_host_workers',0)for v in entries.values()),max_phase_tree_rss_bytes=max(v.get('peak_tree_rss_bytes',0)for v in entries.values()),
 ranks=1,all_six_thread_limits_one=True,worker_ceiling=72,campaign_ceiling=43200,evidence_protected=3600)
rates={key:spatial['comparisons'][1]['rates'][key]for key in('U','k','epsilon','gradU','gradk','gradepsilon','Pi','gradPi')}
assessment=dict(schema='WJH27Source19MMSCheckpointV1',status='PASS_SOURCE19_COUPLED_MMS_FOUNDATIONS_ONLY',full_H27_qualification=False,TASK3_READY=False,
 active_model='FiniteNuIsotropicWallFullWJCollarV2',v4_implemented=False,mac_parent='b97214df4bc93624660748cae3858d6936defde6',source_identity=source,runtime=identity(r/'science/runtime-v19.json'),helpers=identity(t/'overlay.json'),
 temporal=dict(report=identity(o/'temporal/independent.json'),status=temporal['status'],pairs=temporal['temporal_gate']['pairs'],controls=temporal['comparisons'],runs={key:{k:v for k,v in value.items()if k not in('assemblies',)}for key,value in temporal['runs'].items()}),
 spatial=dict(report=identity(o/'spatial/independent.json'),status=spatial['status'],fine_pair_rates=rates,comparisons=spatial['comparisons'],runs={key:{k:v for k,v in value.items()if k not in('assemblies','steady_solver')}for key,value in spatial['runs'].items()}),
 moving_cutoff=cross,cases=case_stats,resources=resource,
 retained_failures_and_open_gates='Physical reference inner-k .13279468209091205/.44168324293315253 exceeds .05/.10. Selected physical h/p/q/dt, required-duration four-case pilots, monitoring/statistics/storage and complete cost remain unqualified; no full Task3 launch.',
 aggregate_account_scope='Same account, all earlier failed and expected-rejection charges preserved; native advancement/export timers are subsets of whole operations and are not added twice.')
with(out/'native-timings.csv').open('x',newline='')as f:
 writer=csv.writer(f);writer.writerow(['case','cells','dt','step','time','whole_native_seconds','construction_seconds','initialisation_seconds','advancement_seconds','export_seconds','Newton_iterations','Krylov_iterations'])
 for name,v in case_stats.items():writer.writerow([name,v['cells'],v['dt'],v['end_step'],v['end_time'],v['whole_native_seconds'],v['construction_seconds'],v['initialisation_seconds'],v['complete_advancing_seconds'],v['final_export_seconds'],json.dumps(v['newton_iterations']),json.dumps(v['krylov_iterations'])])
try:
 os.environ['MPLCONFIGDIR']=str(o/'matplotlib-cache')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,2,figsize=(9,4));colors=dict(U='#1f77b4',k='#d62728',epsilon='#2ca02c')
 for key,color in colors.items():
  series_labels=['dt08','dt04','dt02'];xs=[temporal['runs'][name]['dt']for name in series_labels];ys=[temporal['runs'][name]['errors'][key]['RMS']for name in series_labels]
  axes[0].loglog(xs,ys,'o-',color=color,label=key)
  series_labels=['h2','h4','h8'];xs=[.5,.25,.125];ys=[spatial['runs'][name]['errors'][key]['RMS']for name in series_labels]
  axes[1].loglog(xs,ys,'o-',color=color,label=key)
 axes[0].set_title('Advancing Polynomial MMS, T=0.4');axes[0].set_xlabel('dt');axes[1].set_title('Frozen Historical MMS, t=0');axes[1].set_xlabel('uniform x spacing')
 for ax in axes:ax.set_ylabel('reference-scaled RMS error');ax.invert_xaxis();ax.grid(True,which='both',alpha=.2);ax.legend()
 fig.suptitle('Source19 independent coupled MMS errors');fig.tight_layout();png=out/'mms-convergence.png';svg=out/'mms-convergence.svg';fig.savefig(png,dpi=180);fig.savefig(svg,metadata={'Date':None});plt.close(fig)
 assessment['visualization']=dict(png=identity(png),svg=identity(svg),scope='Authoritative independent RMS integrals; rate and contamination gates are assessed separately in the reports.')
except ImportError as error:assessment['visualization_unavailable']=str(error)
atomic(out/'assessment.json',assessment)
primary='/media/jay/新加卷1/nektar-v2026';head=subprocess.check_output(['git','-C',primary,'rev-parse','HEAD'],text=True).strip();status=subprocess.check_output(['git','-C',primary,'status','--porcelain'],text=True);assert head=='f4b70baec10016dd098005bc5a5507b87e4a57bf'and status==' M docs/developer-guide/src/library\n'
atomic(out/'preserved-checkout.json',dict(path=primary,head=head,status=status,preserved=True))
atomic(out/'inventory.json',dict(schema='WJH27Source19MMSEvidenceV1',TASK3_READY=False,account_closed=False,parent=identity(r/'evidence/native-cutoff-recovery-v2i2/inventory.json'),raw_artifacts=list({v['path']:v for v in raw}.values()),phase_labels=labels,
 files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]))
archive=r/'evidence/native-mms-requalification-v2j.tar.gz'
with tarfile.open(archive,'w:gz')as tar:tar.add(out,arcname=out.name)
print(json.dumps(dict(archive=identity(archive),status=assessment['status'],resources=resource)))
