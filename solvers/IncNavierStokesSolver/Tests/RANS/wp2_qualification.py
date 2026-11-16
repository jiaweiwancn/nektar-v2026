# SPDX-License-Identifier: MIT
"""Frozen WP2 physical budgets; absolute pressure and full outer epsilon."""
import math
import numpy as np
from run_cls_epawf_picard import require


def physical_gates(s, us, K0, complete=False):
    gates={}
    def limit(name,value,bound):
        gates[name]=dict(value=float(value),limit=bound,passed=math.isfinite(value) and value<=bound)
    n=s['unanchored_free_block_norms']
    limit('four_unforced_block_norms',max(n[:4]),1e-3)
    limit('weak_continuity_unused_gauge',max(n[4:7]),1e-9)
    if 'physical_essential_max_scaled' in s:limit('physical_prescribed_boundary_data',s['physical_essential_max_scaled'],1e-9)
    limit('strong_divergence_RMS_over_ustar_H',s['strong_divergence_RMS_over_ustar_H'],1e-3)
    limit('wall_shear',s['wall_shear_RMS_relative_error'],1e-2)
    for region,v in enumerate(s['regional_errors_inner_outer']):
        r,p=v['relative_L2'],v['pointwise_normalized_maxima'];prefix=('inner','outer')[region]
        for name,i,bound in [('velocity',0,1e-3),('W',1,1e-3),('k',2,1e-2),('nut',4,1e-2),
                              ('full_velocity_gradient',5,1e-2),('pi',8,1e-2),('mechanical_pressure',9,1e-2),
                              ('selected_shear',10,1e-2),('normal_stress',11,1e-2)]:
            limit(prefix+'_L2_'+name,r[i],bound)
        for name,i in [('U',0),('W',1),('k',2),('nut',4)]:limit(prefix+'_point_'+name,p[i],1e-2)
        if region:
            limit('outer_epsilon_L2',r[3],1e-2);limit('outer_epsilon_point',p[3],1e-2)
    for band in s['log_height_bands']:
        # Volume is column0; all remaining columns use local positive profiles.
        v=band['volume_then_local_relative_RMS_U_k_epsilon_nut_W']
        key='band_'+str(band['bounds'])
        limit(key+'_U',v[1],1e-3)
        limit(key+'_k_epsilon_nut',max(v[2:5]),1e-2);limit(key+'_W',v[5],1e-3)
    b=s['boundary_status']
    inlet=abs(next(row[3] for row in b['faces'] if row[0]==0 and row[1]==-1))
    limit('normalized_net_volume_flux',abs(b['net_volume_flux'])/inlet,1e-9)
    limit('prescribed_top_net_flux_adjustment',abs(b['top_net_flux_adjustment']),0.)
    if complete:
        interface=s['interface_raw'];limit('interface_normal_full_RMS',max(abs(interface[i]) for i in [2,3])/(us*us),1e-2)
        limit('interface_k_flux',abs(interface[5])/(us*K0),1e-2)
        d=s['scalar_domain'];gates['native_scalar_domain']=dict(passed=bool(d['admissible']))
        limit('cutoff_relative_budget',d['maximum_budget_relative_defect'],1e-10)
        cov=min(s['covariance_inner_outer'][:2]+s['independent_covariance_inner_outer'])
        gates['covariance_margin']=dict(value=cov,minimum=1e-4*K0,passed=cov>=1e-4*K0)
    else:
        cov=min(s['volume_covariance_inner_outer'])
        gates['volume_covariance_margin']=dict(value=cov,minimum=1e-4*K0,passed=cov>=1e-4*K0)
    return dict(passed=all(x['passed'] for x in gates.values()),gates=gates,
                scope='complete native and independent audit' if complete else 'per-accepted-state physical GL16/stations; pair with native sided hard audits')


def common_differences(a,b,us,K0,H):
    require(a['columns']==b['columns'] and a['weights']==b['weights'] and a['regions']==b['regions'],
            'comparison needs the same union physical partition and integration weights')
    va,vb=np.asarray(a['values']),np.asarray(b['values']);w=np.asarray(a['weights']);outer=np.asarray(a['regions'])==1
    d=va-vb;out={}
    names=a['columns'].split(',')
    for i,name in enumerate(names):
        selected=outer if name=='epsilon' else np.ones(len(w),dtype=bool)
        v=float(np.sqrt(np.dot(w[selected],d[selected,i]**2)))
        out[name]=dict(dimensional_L2=v,pointwise_max=float(np.max(np.abs(d[selected,i]))))
    out['pressure_convention']='absolute, both regional fields retained; no mean shift'
    out['U_W_combined_L2']=math.hypot(out['U']['dimensional_L2'],out['W']['dimensional_L2'])
    out['scales']=dict(ustar=us,K0=K0,H=H,stress=us*us)
    return out


def verify_report(d,us,K0,allow_recovery=False):
    results=[]
    for row in d['wp2_per_state_monitoring']:
        checks=physical_gates(row['physical_checks'],us,K0)
        # Recovery may start outside the reference-error budget, by design.
        results.append(dict(time=row['time'],**checks))
    if not allow_recovery:require(all(x['passed'] for x in results),'WP2 physical preservation budget failed')
    # New native reporting exports the immutable physical essential map. Use
    # it to check all retained coefficient vectors (including older-runtime
    # prefixes), separately from the expected post-k dependent trace lag.
    table=next((x['physical_checks'] for x in d['wp2_per_state_monitoring'] if 'physical_essential_rows' in x['physical_checks']),None)
    if table:
        ids=[int(i) for i in table['physical_essential_rows']];values=table['physical_essential_values']
        for state in d['accepted_states']:
            require(all(state[i]==v for i,v in zip(ids,values)),'physical essential coefficient changed')
    for a in d['state_audits']:
        require(a['accepted_awf_hard_admissible']==[1] and min(a['accepted_awf_covariance'][:2])>=1e-4*K0,'native hard or positive covariance margin failed')
        require(max(abs(a['interface'][i]) for i in [2,3])/(us*us)<=1e-2 and abs(a['interface'][5])/(us*K0)<=1e-2,'native interface budget failed')
    for row in d['pimple_attempts']:
        k=row['k_bound'];require(k['activations']==0 and k['integral_total']==0 and k['mass_norm_correction']==0,'equilibrium/recovery required artificial bound supply')
    return results


def quarter_budget_differences(a,b,us,K0,H):
    require(a['weights']==b['weights'] and a['regions']==b['regions'],
            'changed physical comparison grid')
    ra,rb=np.asarray(a['reference_values']),np.asarray(b['reference_values'])
    reference_roundoff=float(np.max(abs(ra-rb)/np.maximum(abs(ra),1e-15)))
    # The same analytical ReferenceK integral is split on the respective native
    # inner partitions. Record its integration roundoff, without changing the
    # fixed normalization or any physical accuracy tolerance.
    require(reference_roundoff<=2e-12,'changed immutable reference beyond integration roundoff')
    va,vb,ref=(np.asarray(x) for x in [a['values'],b['values'],a['reference_values']]);w=np.asarray(a['weights']);regions=np.asarray(a['regions'])
    names=a['columns'].split(',');gates={}
    for region in [0,1]:
        sel=regions==region;measure=float(sum(w[sel]));d=va[sel]-vb[sel];r=ref[sel];weights=w[sel]
        for i,name in enumerate(names):
            if name=='epsilon' and region==0:continue
            if i in [10,11,12,13] or name=='production':continue
            if name=='U':den=math.sqrt(float(np.dot(weights,r[:,i]**2)));budget=1e-3
            elif name=='W':den=us*math.sqrt(measure);budget=1e-3
            elif name in ['pi','mechanical_pressure','total_shear','total_normal_stress']:den=us*us*math.sqrt(measure);budget=1e-2
            else:den=math.sqrt(float(np.dot(weights,r[:,i]**2)));budget=1e-2
            value=math.sqrt(float(np.dot(weights,d[:,i]**2)))/den
            gates[str(region)+'_'+name]=dict(value=value,limit=.25*budget,passed=value<=.25*budget)
            point_scale=np.maximum(abs(r[:,i]),us) if name=='U' else np.full(len(r),us) if name=='W' else np.full(len(r),us*us) if name in ['pi','mechanical_pressure','total_shear','total_normal_stress'] else abs(r[:,i])
            point=float(np.max(abs(d[:,i])/point_scale))
            gates[str(region)+'_'+name+'_point']=dict(value=point,limit=.0025,passed=point<=.0025)
        value=math.sqrt(float(np.sum(weights[:,None]*d[:,10:14]**2)))/math.sqrt(float(np.sum(weights[:,None]*r[:,10:14]**2)))
        gates[str(region)+'_full_velocity_gradient']=dict(value=value,limit=.0025,passed=value<=.0025)
    return dict(passed=all(x['passed'] for x in gates.values()),gates=gates,reference_integration_roundoff=reference_roundoff,
                pressure='absolute; original relative regional offset retained')


def homogeneity(samples,cold,us,K0):
    require(samples['weights']==cold['weights'] and samples['regions']==cold['regions'],'homogeneity changed sampling partition')
    nz=(len(samples['physical_z_partition'])-1)*int(samples['order']);v=np.asarray(samples['values']);c=np.asarray(cold['values']);w=np.asarray(samples['weights']);ref=np.asarray(samples['reference_values'])
    require(len(v)%nz==0,'common tensor sample shape')
    shape=(len(v)//nz,nz,len(v[0]));cv=c.reshape(shape);vv=v.reshape(shape);ww=w.reshape(shape[:2]);weights_x=ww[:,0]/sum(ww[:,0])
    # Radial injection is exactly constant in x; retain and check its roundoff.
    cold_x_range=float(np.max(np.ptp(cv[:,:,:4],axis=0)))
    require(cold_x_range<=1e-12,'cold reference is not the prescribed homogeneous native inlet')
    inlet=np.sum(cv*weights_x[:,None,None],axis=0);mean=np.sum(vv*weights_x[:,None,None],axis=0)
    out={}
    for label,d in [('inlet_to_station',vv-inlet[None,:,:]),('horizontal_variation',vv-mean[None,:,:])]:
        d=d.reshape(v.shape);r=np.asarray(samples['regions']);metrics={}
        for region in [0,1]:
            selected=r==region;measure=sum(w[selected]);row={}
            for i,name in enumerate(samples['columns'].split(',')):
                if name=='epsilon' and region==0:continue
                if i>=10 or name=='production':continue
                if name=='U':scale=np.maximum(abs(ref[selected,i]),us);norm=math.sqrt(float(np.dot(w[selected],ref[selected,i]**2)));budget=1e-3
                elif name=='W':scale=np.full(sum(selected),us);norm=us*math.sqrt(measure);budget=1e-3
                elif name in ['pi','mechanical_pressure','total_shear','total_normal_stress']:scale=np.full(sum(selected),us*us);norm=us*us*math.sqrt(measure);budget=1e-2
                else:scale=abs(ref[selected,i]);norm=math.sqrt(float(np.dot(w[selected],scale**2)));budget=1e-2
                value=math.sqrt(float(np.dot(w[selected],d[selected,i]**2)))/norm
                point=float(max(abs(d[selected,i])/scale))
                row[name]=dict(relative_L2=value,pointwise_scaled_max=point,budget=budget,passed=value<=budget and point<=.01)
            metrics[str(region)]=row
        out[label]=metrics
    out['cold_x_range']=cold_x_range;out['passed']=all(v['passed'] for key in ['inlet_to_station','horizontal_variation'] for row in out[key].values() for v in row.values())
    out['scope']='common physical panels at all streamwise Gauss sections; actual declared-station maxima are also included in every-state native monitoring; cold radial injection supplies the unchanged native inlet trace'
    return out
