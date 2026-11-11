#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H2 scientific tables/figures from immutable native measurements."""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h2 import load, save, sha
from wj_h2_profiles import Completion, QE, US, NU, D0, TOP, wall


def csv_file(path, rows):
    with Path(path).open('w') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def boundary_summary(report):
    d=report['open_boundary'];exterior={}
    for r in d['boundary_rows']:
        v=exterior.setdefault(r['role'],dict(measure=0.,mass=0.,flux=np.zeros(4)))
        v['measure']+=r['weight'];v['mass']+=r['weight']*np.dot(r['U'],r['normal'])
        v['flux']+=r['weight']*np.array(r['flux'])
    for v in exterior.values():v['flux']=v['flux'].tolist()
    interface=d['interface'];grounds=[v for v in d['boundary_rows'] if v['role']=='ground']
    return dict(mass_flux=d['mass_flux'],exterior=exterior,pressure_gauge_unknowns=d['gauge_unknowns'],
        fictitious_top_flux_repair=d['boundary_flux_correction'],
        signed_ground_outward_Tn_x_range=[min(v['actual_traction'][0] for v in grounds),max(v['actual_traction'][0] for v in grounds)],
        positive_forward_wall_stress_mean=-sum(v['weight']*v['actual_traction'][0] for v in grounds)/8,
        interface=dict(traction_point_max=max(np.linalg.norm(v['traction_defect']) for v in interface),
            mass_point_max=max(abs(v['mass_flux_defect']) for v in interface),
            shared_k_max=max(abs(v['shared_k_defect']) for v in interface),
            opposing_k_diffusive_flux_max=max(abs(sum(v['k_flux_inner_outer'])) for v in interface),
            point_epsilon_trace_max=max(abs(np.diff(v['epsilon_trace_target'])[0]) for v in interface),
            mechanical_normal_balance_point_departure=max(abs(np.diff(v['mechanical_jump_actual_target'])[0]) for v in interface),
            K_range=[min(v['K'] for v in interface),max(v['K'] for v in interface)],
            cutoff_range=[min(v['cutoff'] for v in interface),max(v['cutoff'] for v in interface)]))


def run_summary(record):
    p=Path(record['report']);d=load(p);attempts=d['pimple_attempts'];last=d['state_audits'][-1]
    snapshots=[]
    for a in attempts:
        if a.get('k_snapshot'):
            s=load(a['k_snapshot']);assert sha(a['k_snapshot'])==a['k_snapshot_sha256']
            same=np.asarray(s['raw'],dtype='<f8').tobytes()==np.asarray(s['corrected'],dtype='<f8').tobytes()
            assert same,'material numerical correction in retained physical run'
            snapshots.append(dict(path=a['k_snapshot'],sha256=sha(a['k_snapshot']),raw_corrected_identical=same))
    counts={k:sum(a[k] for a in attempts) for k in ('momentum_predictor_solves','pressure_system_solves',
        'k_transport_solves','epsilon_transport_solves','bound_calls','bound_activations','additional_governing_equation_solves')}
    assert all(a['scheduled_complete'] for a in attempts)
    for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):
        assert counts[k]==len(attempts)
    return dict(report=str(p),report_sha256=sha(p),xml_sha256=sha(record['xml']),state_size=len(d['state']),
        dt=d['dt'],accepted_steps=d['steps'],accepted_time=d['time'],flow_throughs=d['time']/8,
        attempts=len(attempts),failure=d['failure'],solve_counts=counts,
        solve_dimensions={k:attempts[0][k] for k in ('velocity_size','pressure_gauge_size','assembled_flow_size','k_solve_size','epsilon_solve_size')} if attempts else {},
        sweep_seconds=sum(a['seconds'] for a in attempts),k_snapshots=snapshots,
        linear_residual_max={k:max(a[k] for a in attempts) for k in ('predictor_linear_residual','pressure_linear_residual',
            'corrected_continuity_gauge','k_frozen_linear_residual','epsilon_frozen_linear_residual')} if attempts else {},
        imposed_epsilon_trace_max=max((a['refreshed_trace_defect'] for a in attempts),default=0.),
        final_conditions=d['wj_audit_packets'][-1]['conditions'],
        final_root=last['wj_root_minimum_scaled_residual'],final_covariance_minima=last['wj_covariance'],
        final_principal_minima=last['wj_rank_one_xz_principal_minima_over_nu_eff'],
        raw_work_defect_and_inventory_divergence=last['wj_work_defect_raw_divergence'],
        boundary=boundary_summary(d))


def timer_summary(path):
    values={}
    for line in Path(path).read_text().splitlines():
        r=json.loads(line);v=values.setdefault(r['name'],dict(calls=0,inclusive_seconds=0.,exclusive_seconds=0.))
        v['calls']+=1;v['inclusive_seconds']+=r['wall_seconds'];v['exclusive_seconds']+=r['exclusive_wall_seconds']
    return values


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--measurements',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    root=a.root;measure=load(a.measurements/'summary.json');index=load(a.measurements/'index.json');study=load(root/'study-v1/progress.json')
    assert study['status']=='OBSTRUCTED'
    c=Completion();summary=dict(schema='WJH2ScientificSummaryV1',status='OBSTRUCTED',first_limiting_gate=study['first_limiting_gate'],
        unreached_target=8.,unperformed=['developed one-flow-through precursor','evolved missing-input companion','physical hill prediction'],
        calibration=c.summary(),comparisons=measure['comparisons'],quadrature=measure['quadrature'],
        first_incident=measure['first_incident'],measurements_sha256=sha(a.measurements/'summary.json'),runs={},profiles={},weak={})
    for key,r in study['runs'].items():summary['runs'][key]=run_summary(r)
    half=load(root/'diagnosis-v1/result.json');summary['runs']['refined-half-dt-diagnostic']=run_summary(half)
    summary['missing_input_family']=load(root/'completions-compiled-v1/summary.json')
    norms={};profile_rows=[];discrepancy=[];profiles={}
    for label,identity in measure['profiles'].items():
        assert sha(identity['native'])==identity['native_sha256'];m=load(identity['native']);d=m['open_boundary'];rows=d['stations'];profiles[label]=d
        norms[label]=m['physical_mass_dual_norms'];errors=[]
        for station in (0.,2.,4.,6.,8.):
            station_rows=[v for v in rows if abs(v['physical'][0]-station)<1e-12]
            for panel in ('a','b','c','d'):
                group=[]
                for source in [v for v in c.rows if v['panel']==panel]:
                    z=float(source['y_over_H']);candidates=[v for v in station_rows if abs(v['physical'][2]-z)<2e-11]
                    assert len(candidates)==1,'experimental observation ownership/height differs'
                    v=candidates[0];R=np.array(v['R']);assert R[0,0]>=0 and R[2,2]>=0
                    value=v['U'][0] if panel=='a' else -R[0,2] if panel=='b' else math.sqrt(R[0,0] if panel=='c' else R[2,2])/US
                    observed=float(source['value']);group.append(value-observed)
                    profile_rows.append(dict(case=label,time=d['time'],x=station,panel=panel,y_over_H=z,y_mm=z*40,
                        inner=v['inner'],geometry=v['geometry'],measured=observed,modeled=value,error=value-observed,
                        U=v['U'][0],W=v['U'][2],Rxx=R[0,0],Ryy=R[1,1],Rzz=R[2,2],Rxz=R[0,2],
                        k=v['k'],epsilon=v['epsilon'],pi=v['pi'],mechanical_pressure=v['mechanical_pressure'],
                        actual_P=v['actual_production'],counterfactual_P=v['counterfactual_production'],destruction=v['destruction']))
                error=dict(case=label,time=d['time'],x=station,panel=panel,count=len(group),bias=float(np.mean(group)),
                    rms=float(np.sqrt(np.mean(np.square(group)))),max_absolute=float(max(abs(np.array(group)))))
                discrepancy.append(error);errors.append(error)
        summary['profiles'][label]=dict(identity=identity,discrepancy=errors,boundary=boundary_summary(dict(open_boundary=d)))
    summary['physical_mass_dual_norms']=norms
    csv_file(a.output/'profiles-at-measured-heights.csv',profile_rows);csv_file(a.output/'experimental-discrepancy.csv',discrepancy)
    for label,v in index.items():
        if 'independent_weak' not in v:continue
        m=load(v['native']);packet=m['weak_checks'];ind=v['independent_weak'];integrals=dict(actual_inner_P=0.,reference_inner_P=0.,E_c=0.)
        for row in packet['samples']:
            if not row['inner']:continue
            w=row['weight'];k=row['k'];G=np.array(row['G']);p0=wall(row['K'],c.H,c.h);nt=p0['a']*max(row['chart'][2]-p0['zv'],0)
            integrals['actual_inner_P']+=w*(.5*nt*np.sum((G+G.T)**2)-2*k*np.trace(G)/3)
            integrals['reference_inner_P']+=w*(row['K']/QE)**2*nt/(NU+nt)**2
            integrals['E_c']+=w*p0['A']/max(row['chart'][2],p0['rc'])
        summary['weak'][label]=dict(independent=ind,integrals=integrals,boundary=v['independent_boundary'])
    summary['timing']={k:timer_summary(root/'operations'/('study-'+k)/'timing.jsonl') for k in ('base-pilot-off','base-pilot-on','fine-pilot','refined-pilot')}
    csv_file(a.output/'comparison-table.csv',[dict(comparison=key,field=f,**v) for key,row in measure['comparisons'].items() for f,v in row['fields'].items()])
    save(a.output/'summary.json',summary)
    # Scientific figures use immutable physical fields, never an outer law
    # extrapolated into the transported inner region.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':9,'figure.dpi':150})
    fig,axes=plt.subplots(1,4,figsize=(10,4),sharey=True)
    for axis,panel in zip(axes,('a','b','c','d')):
        source=[v for v in c.rows if v['panel']==panel]
        axis.plot([float(v['value']) for v in source],[float(v['y_over_H']) for v in source],'ko',mfc='none',ms=3,label='Fig.3 rough symbols')
        for label,x,style in [('base08',0.,'-'),('base08',4.,'--'),('base08',8.,':'),('refined-fine08',8.,'-.')]:
            selected=[v for v in profiles[label]['stations'] if abs(v['physical'][0]-x)<1e-12 and v['normal_side']>=0]
            selected.sort(key=lambda v:v['physical'][2]);R=np.array([v['R'] for v in selected])
            values=[v['U'][0] for v in selected] if panel=='a' else -R[:,0,2] if panel=='b' else np.sqrt(R[:,0,0] if panel=='c' else R[:,2,2])/US
            axis.plot(values,[v['physical'][2] for v in selected],style,lw=1,label=f'{label}, x={x:g}')
        axis.axhline(D0+c.H,color='grey',lw=.6);axis.set_ylim(D0,8);axis.grid(alpha=.2)
    for ax,text in zip(axes,('U/U_inf','-R_xz/U_inf²','sqrt(R_xx)/u_tau','sqrt(R_zz)/u_tau')):ax.set_xlabel(text)
    axes[0].set_ylabel('original physical height / H_h');axes[-1].legend(fontsize=6)
    fig.suptitle('Flowing rough-flat transient at T U_inf/H_h = 0.08; not developed/calibrated')
    fig.tight_layout();fig.savefig(a.output/'physical-profiles.png');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(10,3.3))
    for name in ('base-pilot-on','fine-pilot','refined-pilot'):
        d=load(study['runs'][name]['report']);t=[v['time'] for v in d['wj_audit_packets']]
        values=[v['conditions']['principal_in_plane']['margins'][0] for v in d['wj_audit_packets']]
        axes[0].plot(t,values,'o-',ms=3,label=name)
    d=load(half['report']);axes[0].plot([v['time'] for v in d['wj_audit_packets']],
        [v['conditions']['principal_in_plane']['margins'][0] for v in d['wj_audit_packets']],'o-',ms=3,label='refined half-dt')
    w=measure['first_incident']['witnesses']['in_plane'];axes[0].plot(.06,w['minimum_over_nu_eff'],'rx',ms=8,label='retained failed candidate')
    axes[0].axhline(0,color='black',lw=.5);axes[0].set_xlabel('T U_inf/H_h');axes[0].set_ylabel('sampled principal form / nu_eff');axes[0].legend(fontsize=6)
    for label in ('startup','base08','refined-fine08'):
        rows=profiles[label]['boundary_rows'];rows=sorted([v for v in rows if v['role']=='ground'],key=lambda v:v['physical'][0])
        axes[1].plot([v['physical'][0] for v in rows],[-v['actual_traction'][0]/US**2 for v in rows],label=label)
        f=profiles[label]['interface'];axes[2].plot([v['physical'][0] for v in f],[np.linalg.norm(v['traction_defect'])/US**2 for v in f],label=label)
    axes[1].set_xlabel('x/H_h');axes[1].set_ylabel('forward wall shear / reported u_tau²');axes[1].legend(fontsize=6)
    axes[2].set_xlabel('x/H_h');axes[2].set_ylabel('point interface traction defect / u_tau²');axes[2].legend(fontsize=6)
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(a.output/'admission-and-interface.png');plt.close(fig)
    print(json.dumps(dict(status='COMPLETE',scientific_status='OBSTRUCTED',output=str(a.output))))


if __name__=='__main__':main()
