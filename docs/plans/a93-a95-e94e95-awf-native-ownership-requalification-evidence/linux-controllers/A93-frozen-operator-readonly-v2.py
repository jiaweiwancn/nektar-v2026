import pathlib,os,sys,gzip,json,shutil
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from qualify_rans_preservation import replay_request

def load(path):
 with (gzip.open(path,'rt') if str(path).endswith('.gz') else pathlib.Path(path).open()) as f:return json.load(f)
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_qualification import THREADS,write_json,exact_vector_hash
from rans_preservation_oracle import Reference,weak_balance
from diagnose_rans_balanced_spatial import run_probe
from rans_spatial_resolution import analyze,common_test,test_parent_geometries,test_effective_coordinates
from rans_duration_contract import residual_reduction
old=r.with_name('a93-a95-e94e95-nonlinear-20261110');runtime=r/'runtime-resolution-v3'
os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
common=r/'diagnosis/A93-Natural-Accepted/common';cp=load(common/'packet.json.gz')
out=r/'selected-frozen-operator-comparison/A93-v2';out.parent.mkdir(exist_ok=True)
paths={m:old/'evolved-audits'/('A93-N-Accepted-duration-0.002-0.5-'+m)/'T-0.5000' for m in ('B5','B11')}
request=dict(schema='LegacyStateCorrectedOperatorReadOnlyV1',source_operator='LegacyCoordinateToleranceV1',target_operator='GeneratedVolumeOwnerV1',common_sha256=sha(common/'packet.json.gz'),source_native_packets={m:sha(p/'packet.json.gz') for m,p in paths.items()},runtime=runtime_identity(runtime/'RANSPreservationProbe'),history='no steps, no history export or signature editing')
seal=RequestManifest(out,request,False);ref=Reference('A93','Natural','Accepted');rows={};oldrows={};files=[]
for m,p in paths.items():
 folder=old/'campaigns'/('A93-N-Accepted-duration-0.002-'+m);native=load(folder/'T-0.5000/case.evolution.json');seed=load(folder/'initial-packet.json')
 dest=out/m;dest.mkdir();shutil.copy2(p/'case.xml',dest/'case.xml')
 q=replay_request(seed,[],state=native['state'],time=native['time'],export_basis=True)
 q.update(operation='legacy_operator_replay',source_numerical_operator='LegacyCoordinateToleranceV1',target_numerical_operator='GeneratedVolumeOwnerV1',fixed_boundary=native['initialization_receipt']['live_fixed_boundary'],export_resolution=True)
 packet,seconds=run_probe(runtime,dest,q)
 assert exact_vector_hash(packet['observation']['state'])==exact_vector_hash(native['state'])
 prior=load(p/'packet.json.gz')
 assert [(x['coordinate'],x['weight']) for x in prior['observation']['volume']]==[(x['coordinate'],x['weight']) for x in packet['observation']['volume']]
 cross=out/(m+'-common');cross.mkdir();shutil.copy2(common/'case.xml',cross/'case.xml')
 req=replay_request(cp,test_effective_coordinates(packet));req.update(export_test_basis=True,grid_test_geometries=test_parent_geometries(cp,packet))
 xp,_=run_probe(runtime,cross,req)
 row=dict(status='pass',time=.5,configuration=dict(model='A93',policy='Natural',wall='Accepted'),spatial=analyze(packet,ref),common_test=common_test(cp,packet,xp,ref),analytical_weak_balance=weak_balance(ref,packet['observation']),source_native_report_sha256=sha(folder/'T-0.5000/case.evolution.json'),state_sha256=exact_vector_hash(native['state']),volume_coordinates_weights_unchanged=True,performed_time_steps=0)
 write_json(dest/'analysis.json',row);rows[m]=row;oldrows[m]=load(p/'analysis.json')
 files += [str(dest/'packet.json.gz'),str(cross/'packet.json.gz'),str(dest/'analysis.json')]
result=dict(status='read_only_comparison_completed',source_operator='LegacyCoordinateToleranceV1',target_operator='GeneratedVolumeOwnerV1',source_reduction=residual_reduction(oldrows['B5'],oldrows['B11'],'A93'),corrected_operator_on_legacy_coefficients=residual_reduction(rows['B5'],rows['B11'],'A93'),members=rows,scope='old coefficients evaluated without advancing or resuming any history; not corrected evolution qualification')
write_json(out/'results.json',result);seal.seal(out,['results.json',*files])
