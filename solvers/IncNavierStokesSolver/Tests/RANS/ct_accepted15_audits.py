# SPDX-License-Identifier: MIT
"""Unchanged complete per-state observations and scheduled current-state actions."""
import gc,os,resource,time,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
import ct_accepted8_native as n
from ct_accepted5_traction import profiles
from ct_accepted10_compare import gradients,outlet
from ct_accepted3_metrics import inventory
from ct_accepted11_target import bounded_budgets
from ct_accepted_pilot import replay,constitutive_comparison
from ct_accepted_open import probe,require
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import sha

def external_cpu():
    """Actual live CPU of this controller's guards/launchers, excluding child CPU."""
    hz=os.sysconf('SC_CLK_TCK');pid=os.getppid();owned=set();command_parents=[]
    def info(pid):
        root=Path('/proc')/str(pid);f=(root/'stat').read_text().rsplit(')',1)[1].split()
        return int(f[1]),(int(f[11])+int(f[12]))/hz,(root/'cmdline').read_bytes().split(b'\0')
    while pid>1:
        try:parent,cpu,args=info(pid)
        except (FileNotFoundError,ProcessLookupError):break
        names={Path(os.fsdecode(v)).name for v in args if v}
        if names&{'ct_accepted15_budget.py','ct_accepted3_command.py','ct_accepted15_command.py'}:
            owned.add(pid)
            if names&{'ct_accepted3_command.py','ct_accepted15_command.py'}:command_parents.append(pid)
        pid=parent
    for parent in command_parents:
        try:children=(Path('/proc')/str(parent)/'task'/str(parent)/'children').read_text().split()
        except FileNotFoundError:continue
        for child in children:
            try:_,cpu,args=info(int(child))
            except (FileNotFoundError,ProcessLookupError):continue
            if any(Path(os.fsdecode(v)).name=='guard_rans_worker_budget.py' for v in args if v):owned.add(int(child))
    return sum(info(pid)[1] for pid in owned)

def scope(function,*args,**kwargs):
    start=time.monotonic();cpu=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN);external=external_cpu()
    value=function(*args,**kwargs)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    usage=dict(CPU=time.process_time()-cpu+after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall=time.monotonic()-start,self_CPU=time.process_time()-cpu,external_monitor_CPU=max(0.,external_cpu()-external))
    return value,usage

def bundle(ini,state,c,weights,count):
    return dict(metrics=n.physical_metrics(ini,state,c,weights,count),profiles=profiles(ini,state,c),gradients=gradients(ini,state,c,weights,count),outlet=outlet(ini,state,c),inventory=inventory(ini['test_basis'],ini['observation']['grid'],np.asarray(state),ini['observation']['offsets'],'A95',ini['observation']['grid_sides']))

def complete(rt,folder,ini,report,c,weights,count,actions=False,full=True):
    value=dict(per_state=[],phase_usage=dict(per_state=[]))
    # Skip only the restart input, whose independently checked predecessor is
    # retained. Every newly accepted state gets the complete unchanged bundle.
    for t,state in zip(report['audit_times'][1:],report['accepted_states'][1:],strict=True):
        q,usage=scope(bundle,ini,state,c,weights,count);require(q['metrics']['status']=='pass','independent accepted-state physical gate')
        value['phase_usage']['per_state'].append(usage)
        value['per_state'].append(dict(time=t,**q))
    if not full:
        value['scheduled_full_audit']=False
        value['scope']='every accepted-state physical/grid/gradient/traction/mass/inventory observation; full scalar/constitutive controls are separately scheduled'
        return value
    value['scheduled_full_audit']=True
    value['budgets'],value['phase_usage']['scalar_budgets']=scope(bounded_budgets,rt,folder,ini,report,c)
    require(value['budgets']['gate']['status']=='pass','all four current consecutive-state scalar controls')
    grid,sides,_,_=n.inventory(c)
    (q,e),value['phase_usage']['constitutive']=scope(probe,rt,folder,replay(ini,report['state'],report['time'],grid,grid_sides=sides,export_constitutive=True),'constitutive-replay')
    require(e['exit_code']==0,'current constitutive replay failed')
    value['constitutive']=constitutive_comparison(q,c);value['constitutive_execution']=e
    require(value['constitutive']['status']=='pass','current constitutive gate')
    del q;gc.collect()
    if actions:value['actions'],value['phase_usage']['actions']=scope(current_actions,rt,folder,ini,report,c,count)
    return value

def current_actions(rt,folder,ini,report,c,count):
    off=ini['observation']['offsets'];rr=[];cc=[];vv=[]
    for i,(row,region) in enumerate(zip(ini['test_basis'][:count],ini['observation']['grid_regions'][:count],strict=True)):
        for j,v in zip(row['pressure']['ids'],row['pressure']['phi'],strict=True):rr.append(i);cc.append(off[4+region]+j-off[4]);vv.append(v)
    B=csr_matrix((vv,(rr,cc)),shape=(count,off[6]-off[4]));base=None;result={}
    for label,orders in (('native',None),('q48',(48,12,12,16,64)),('q64',(64,16,12,20,64))):
        dest=folder/('actions-'+label);dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
        if orders:
            for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
        write(root,dest/'case.xml')
        q,e=probe(rt,dest,replay(ini,report['state'],report['time'],[],export_faces_only=True,export_sparse_backend=True,sparse_alpha=1/report['dt']),'current')
        require(e['exit_code']==0 and not q.get('failure'),'current-state sparse action control:'+label)
        b=q['sparse_backend'];require(q['observation']['state']==report['state'] and b['pressure_C_nnz']==0,'current action state/absolute pressure convention')
        for kind,key in (('predictor','predictor_residual'),('pressure','true_pressure_residual')):require(b[key]<=2e-12*b[kind+'_rhs_maximum']+1e-25,'current original-system residual')
        require(max(b[k] for k in ('predictor_factor_nnz','consistent_mass_factor_nnz','pressure_preconditioner_factor_nnz'))<=40000000,'current produced factor fill')
        require(max(b[k] for k in ('full_native_G_action_defect','consistent_mass_correction_true_residual','splitting_identity'))<=2e-10,'current full G/consistent mass/splitting')
        if base is None:base=np.asarray(b['pressure_increment'])
        difference=float(np.max(abs(B@(np.asarray(b['pressure_increment'])-base)),initial=0))/c['ustar']**2
        require(difference<=1e-8,'current governing-quadrature absolute pressure uncertainty')
        result[label]=dict(status='pass',execution=e,packet_sha256=sha(dest/'current.json'),pressure_change_scaled_maximum=difference,maximum_factor_entries=max(b[k] for k in ('predictor_factor_nnz','consistent_mass_factor_nnz','pressure_preconditioner_factor_nnz')))
        del q;gc.collect()
    return result
