from pathlib import Path
import os,json,shutil,tarfile,csv,subprocess
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/cutoff-owner-transition-v19';t=r/'tools-cutoff-owner-transition-v19';out=r/'evidence/native-cutoff-owner-transition-v2k';out.mkdir(parents=True);raw=[]
def read(p):return json.loads(p.read_text(),parse_int=lambda token:-0.0 if token=='-0'else int(token))
j=read(o/'independent.json');assert j['foundation_gate_pass']is True;cases=read(o/'cases.json');assert j['source_identity']==cases['source_identity']=='56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
account=read(r/'accounting.json');labels=[key for key in account['entries']if key.startswith('owner-transition-')];entries={key:account['entries'][key]for key in labels}
def copy(p):
 target=out/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(identity(p))
drop={'coefficients','elements','input_history','current_advection','previous_advection','free_velocity_ids','k_free_ids','epsilon_free_ids','essential_flow_lift','essential_k_lift','essential_epsilon_lift','mms_bulk_body_rhs','mms_sampled_body_rhs'}
def compact(v):
 if isinstance(v,dict):return{k:compact(value)for k,value in v.items()if not(v.get('schema')=='WJH27NativeCoupledStateV2'and k in drop)}
 if isinstance(v,list):return[compact(value)for value in v]
 return v
for p in sorted(o.rglob('*')):
 if not p.is_file():continue
 raw.append(identity(p))
 if p.suffix in('.bin','.fld'):continue
 if p.suffix in('.xml','.opt','.jsonl')or(p.suffix=='.json'and not p.name.endswith('.result.json')and'.science-'not in p.name):copy(p)
 elif p.name.endswith('.result.json'):
  value=compact(read(p));value['preserved_raw_artifact']=identity(p);value['publication_scope']='Compact native metadata; complete signed coefficients/histories and every scientific frame remain on Linux.'
  target=out/'compact'/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);atomic(target,value)
for p in sorted((r/'inputs').iterdir()):
 if p.is_file()and('owner-transition'in p.name or p.name.startswith('cutoff_owner_transition_check-')):copy(p)
for p in(t/'h27/cutoff_owner_transition_check.py',t/'transition-overlay.json',t/'owner-transition-admission.md',r/'science/runtime-v19.json',r/'launches/owner-transition-v19.log'):copy(p)
for label in labels:
 for name in('request.json','receipt.json','command.log'):
  p=r/'operations'/label/name
  if p.exists():copy(p)
shutil.copy2(r/'accounting.json',out/'accounting-open-snapshot.json')
native={}
for name in cases['cases']:
 packet=read(o/(name+'.result.json'));d=packet['state'];rows=[json.loads(line)for line in(o/(name+'.steps.jsonl')).read_text().splitlines()]
 native[name]=dict(result=identity(o/(name+'.result.json')),whole_seconds=entries['owner-transition-'+name+'-native-v19']['seconds'],
  construction_seconds=packet['construction_seconds'],initialisation_seconds=packet['initialisation_seconds'],advancing_seconds=packet['complete_advancing_seconds'],step_seconds=packet['complete_step_seconds'],
  export_seconds=packet['final_export_seconds'],science_frame_seconds=packet['scientific_frame_seconds'],science_index_seconds=packet['scientific_index_seconds'],recovery_seconds=packet['recovery_seconds'],
  step=d['step'],time=d['time'],dt=d['dt'],policy=d['coupled_newton']['preconditioner_policy'],newton=[row['nonlinear_iterations']for row in rows],krylov=[row['coupled_newton']['krylov_iterations']for row in rows],
  refreshes=[row['coupled_newton']['refreshes']for row in rows],reuses=[row['coupled_newton']['reuses']for row in rows],native_minima={key:d[key]for key in('minimum_k','minimum_epsilon','minimum_covariance_over_k','minimum_principal_over_mu')})
with(out/'native-timings.csv').open('x',newline='')as f:
 writer=csv.writer(f);writer.writerow(['case','whole_seconds','construction_seconds','initialisation_seconds','advancing_seconds','export_seconds','final_step','final_time','Newton','Krylov','refreshes','reuses'])
 for name,v in native.items():writer.writerow([name,v['whole_seconds'],v['construction_seconds'],v['initialisation_seconds'],v['advancing_seconds'],v['export_seconds'],v['step'],v['time'],json.dumps(v['newton']),json.dumps(v['krylov']),json.dumps(v['refreshes']),json.dumps(v['reuses'])])
def bytes_in(folder):return sum(p.stat().st_size for p in folder.rglob('*')if p.is_file())
source_bytes=sum(bytes_in(p)for p in r.iterdir()if p.is_dir()and p.name.startswith(('source-v','build-v','runtime-v')));science_bytes=bytes_in(r/'science');free=shutil.disk_usage(r).free
assert source_bytes<8*1024**3 and science_bytes<8*1024**3 and free>15*1024**3
resources=dict(source_build_runtime_bytes=source_bytes,science_bytes=science_bytes,free_bytes=free,account_closed=account['closed'],aggregate_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in account['entries'].values()),
 phase_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in entries.values()),max_phase_host_workers=max(v.get('peak_host_workers',0)for v in entries.values()),max_phase_tree_rss_bytes=max(v.get('peak_tree_rss_bytes',0)for v in entries.values()),
 ranks=1,all_six_thread_limits_one=True,worker_ceiling=72,campaign_ceiling=43200,evidence_protected=3600)
assessment=dict(schema='WJH27CutoffOwnerTransitionCheckpointV1',status='PASS_CUTOFF_OWNER_CACHE_HISTORY_ONLY',full_H27_qualification=False,TASK3_READY=False,
 active_model='FiniteNuIsotropicWallFullWJCollarV2',v4_implemented=False,mac_parent='82d11ddeb21118293737fd67d5be265000a256b2',source_identity=cases['source_identity'],runtime=identity(r/'science/runtime-v19.json'),
 independent=identity(o/'independent.json'),cache_comparison=identity(o/'cache-comparison.json'),owner_transitions=j['owner_transitions'],
 sampled_K_change_over_Kin=j['sampled_K_change_over_Kin'],sampled_cutoff_change_over_H=j['sampled_cutoff_change_over_H'],
 frame_comparison_MAX={name:max(value['scaled_MAX']for frame in rows.values()for role in frame.values()for value in role.values())for name,rows in j['frame_coefficient_comparisons'].items()},
 current_equations=j['independent_equations'],retained_recovery_imports=j['retained_recovery_imports'],first_cold_inverse_event=j['first_cold_inverse_event'],
 cases=native,resources=resources,scope='Focused current-rule owner crossing, inverse-cache equivalence and cold/rotation foundation; physical reference, selected h/p/q/dt, required-duration pilots, monitoring/storage and complete production cost remain open.')
try:
 os.environ['MPLCONFIGDIR']=str(o/'matplotlib-cache')
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 fig,axes=plt.subplots(1,2,figsize=(9,4));styles={'fresh':('#1f77b4','--'),'cached':('#d62728','-'),'cold':('#2ca02c',':')}
 for name,(color,style)in styles.items():
  rows=j['trajectories'][name];axes[0].plot([v['time']for v in rows],[v['roots'][0]['x']for v in rows],color=color,ls=style,marker='o',ms=3,label=name)
 axes[0].axhline(.5,color='0.4',lw=1,ls='--');axes[0].set_xlabel('accepted MMS time');axes[0].set_ylabel('current cutoff crossing x');axes[0].legend();axes[0].grid(True,alpha=.2)
 rows=j['trajectories']['cached'];times=[v['time']for v in rows];refresh=[v['refreshes']for v in rows];reuse=[v['reuses']for v in rows]
 axes[1].bar(times,refresh,width=.014,label='refresh');axes[1].bar(times,reuse,bottom=refresh,width=.014,label='reuse');axes[1].set_xlabel('accepted MMS time');axes[1].set_ylabel('Newton preparation bases');axes[1].legend();axes[1].grid(True,axis='y',alpha=.2)
 for transition in j['owner_transitions']:
  for ax in axes:ax.axvline(transition['after']['time'],color='0.4',lw=.8,ls=':')
 fig.suptitle('Source19 actual cutoff owner crossing and inverse reuse');fig.tight_layout();png=out/'owner-transition.png';svg=out/'owner-transition.svg';fig.savefig(png,dpi=180);fig.savefig(svg,metadata={'Date':None});plt.close(fig)
 assessment['visualization']=dict(png=identity(png),svg=identity(svg),scope='Root from one retained original radial side; both radial owners independently validated. Bars count preparation bases, not Krylov applications.')
except ImportError as error:assessment['visualization_unavailable']=str(error)
atomic(out/'assessment.json',assessment)
primary='/media/jay/新加卷1/nektar-v2026';head=subprocess.check_output(['git','-C',primary,'rev-parse','HEAD'],text=True).strip();status=subprocess.check_output(['git','-C',primary,'status','--porcelain'],text=True);assert head=='f4b70baec10016dd098005bc5a5507b87e4a57bf'and status==' M docs/developer-guide/src/library\n'
atomic(out/'preserved-checkout.json',dict(path=primary,head=head,status=status,preserved=True))
atomic(out/'inventory.json',dict(schema='WJH27CutoffOwnerTransitionEvidenceV1',TASK3_READY=False,account_closed=False,parent=identity(r/'evidence/native-mms-requalification-v2j/inventory.json'),
 raw_artifacts=list({v['path']:v for v in raw}.values()),phase_labels=labels,files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]))
archive=r/'evidence/native-cutoff-owner-transition-v2k.tar.gz'
with tarfile.open(archive,'w:gz')as tar:tar.add(out,arcname=out.name)
print(json.dumps(dict(archive=identity(archive),status=assessment['status'],resources=resources)))
