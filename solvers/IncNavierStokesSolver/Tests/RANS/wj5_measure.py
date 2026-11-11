# SPDX-License-Identifier: MIT
"""Passive WJ divergence decomposition and constraint-aware finite estimates."""
import math
import numpy as np
from numpy.testing import assert_allclose
from wj3_measure import FIELDS
from wj4_measure import AMP,boundary_rule,physical_norms
from run_cls_epawf_picard import require


def constraint_estimate(row,coefficient_uncertainty,prescribed_uncertainty=0.,evaluation_uncertainty=0.):
    """Actual E*P action plus fixed-data/lift/evaluation/accumulation reserve.

    coefficient_uncertainty already carries the frozen-system 100x reserve.
    Prescribed data and arithmetic of the lifting/basis evaluation retain their
    own 100x reserve. A measured zero difference never sets this estimate to zero.
    """
    values=[row['free_row_L2'],row['gross_row_L2'],row['absolute_accumulation'],
            row['prescribed_accumulation'],coefficient_uncertainty,prescribed_uncertainty,evaluation_uncertainty]
    require(all(math.isfinite(x) and x>=0 for x in values),'invalid constraint uncertainty input')
    n=max(1,len(row['row']));eps=np.finfo(float).eps
    accumulation=100*eps*(32*n+16)*(1+row['absolute_accumulation']+row['prescribed_accumulation'])
    free=row['free_row_L2']*coefficient_uncertainty
    data=100*row['gross_row_L2']*prescribed_uncertainty
    evaluation=100*evaluation_uncertainty
    return dict(total=free+data+evaluation+accumulation,free_state=free,prescribed_data_and_lift=data,
        basis_geometry_evaluation=evaluation,accumulation=accumulation,free_row_L2=row['free_row_L2'],
        safety_factor=100,scope='finite frozen-system estimate through actual elimination/lift map; not a nonlinear theorem')


def boundary_uncertainty(measurement,arithmetic,rule_settings):
    """Tighten only constrained differential observations; pressure stays free.

    One full-domain pressure alignment is already used by physical_norms. The
    retained global pressure envelope is doubled for its subtraction of the
    full-volume mean. There is no regional or boundary-only pressure gauge.
    """
    weights=boundary_rule(**{k:v for k,v in rule_settings.items() if k!='all_faces'})
    rows=measurement['observations'];values=[];max_free=np.zeros(4);proof=[]
    systems=arithmetic['systems'];floors=arithmetic['estimated_comparison_floor'];sampled=arithmetic['estimated_sampled_comparison_floor']
    for r in rows:
        x,_,z=r['position'];on_boundary=x in (0.,1.) or z==2.
        v=np.zeros(6)
        for j,name in enumerate(FIELDS[:4]):
            if j==3 and r['inner']:continue
            desc=r.get('constraint_rows',[])[j]
            require(desc is not None,'missing native constraint row')
            block='momentum' if j<2 else 'k' if j==2 else 'epsilon'
            minimum=systems[block]['physical_metric_min_eigenvalue'];require(minimum>0,'invalid physical coefficient metric')
            coefficient=floors[name]/math.sqrt(minimum)
            estimate=constraint_estimate(desc,coefficient,evaluation_uncertainty=measurement['arithmetic_evaluation_sensitivity'][j])
            v[j]=estimate['total']
            if on_boundary:max_free[j]=max(max_free[j],desc['free_row_L2'])
            if on_boundary and (desc['free_row_L2']>1e-12 or z in (1.,2.)):
                proof.append(dict(position=r['position'],geometry=r['geometry'],field=name,estimate=estimate))
        v[4]=2*sampled['reduced_pressure'];v[5]=v[4]+2*v[2]/3
        values.append(v)
    q=np.array(values);bw=[];seen=set()
    for r in rows:
        key=(r['position'][0],r['position'][2]);bw.append(weights.get(key,0.) if key not in seen else 0.);seen.add(key)
    bw=np.array(bw);outer=np.array([not r['inner'] for r in rows]);mask=np.array([r['position'][0] in (0.,1.) or r['position'][2]==2. for r in rows])
    result={}
    for j,f in enumerate(FIELDS+('velocity',)):
        active=mask & outer if f=='outer_epsilon' else mask
        a=np.hypot(q[:,0],q[:,1]) if f=='velocity' else q[:,j]
        result[f]=dict(L2=float(np.sqrt(bw[active]@(a[active]**2))),sampled_Linf=float(max(a[active])),measure=float(sum(bw[active])))
    return dict(schema='WJConstraintAwareBoundaryV1',safety_factor=100,fields=result,
        maximum_boundary_free_rows=dict(zip(FIELDS[:4],max_free.tolist())),selected_proofs=proof,
        pressure='unconstrained; retained global row estimate plus one full-domain gauge subtraction; no fixed-trace exemption',
        original_estimator='WJ-4 global sampled envelope retained separately; original UNRESOLVED evidence is unchanged')


def diagnose(m):
    d=m['passive_diagnostics'];require(d and d['governing_solves']==0,'missing passive diagnostics')
    regions=[];mode_rows=[];by_element={}
    for q in d['quadrature']:by_element.setdefault(int(q['geometry']),[]).append(q)
    for r in d['pressure_regions']:
        region=int(r['region']);mass=np.array(r['pressure_mass']);rhs=np.array(r['B_U']);constant=np.array(r['constant_coefficients'])
        # Offline diagnostic mass algebra. No governing solver or trajectory
        # field is called or changed by this pressure-space projection.
        coeff=np.linalg.solve(mass,rhs);total=projected=complement=cross=0.;independent=np.zeros(len(rhs))
        for q in d['quadrature']:
            if q['region']!=region:continue
            ids=np.array(q['pressure_ids'],int);phi=np.array(q['pressure_phi']);div=q['Gxx_Gzz_div'][2];w=q['weight'];p=phi@coeff[ids]
            independent[ids]+=w*phi*div;total+=w*div*div;projected+=w*p*p;complement+=w*(div-p)**2;cross+=w*p*(div-p)
        assert_allclose(independent,rhs,rtol=1e-10,atol=2e-13)
        assert_allclose(np.array(r['B_U'])+r['lambda_m_p'],r['actual_continuity_rows'],rtol=1e-9,atol=2e-12)
        assert_allclose(total,projected+complement,rtol=2e-9,atol=2e-23)
        assert_allclose(float(constant@rhs),d['outward_net_flux'][region],rtol=1e-8,atol=2e-12)
        regions.append(dict(region=region,divergence_L2=math.sqrt(total),Q1_projection_L2=math.sqrt(projected),
            complement_L2=math.sqrt(complement),orthogonality=cross,pressure_coefficients=coeff.tolist(),
            maximum_B_U=float(max(abs(rhs))),maximum_gauge_load=float(max(abs(np.array(r['lambda_m_p'])))),
            maximum_actual_continuity_row=float(max(abs(np.array(r['actual_continuity_rows'])))),
            integrated_divergence=d['divergence_integral'][region],outward_net_flux=d['outward_net_flux'][region]))
    for geometry,qs in by_element.items():
        lo=np.array(qs[0]['lower']);hi=np.array(qs[0]['upper']);vol=(hi[0]-lo[0])*(hi[2]-lo[2]);coeff=np.zeros((3,4,4))
        for q in qs:
            x,_,z=q['position'];tx=2*(x-lo[0])/(hi[0]-lo[0])-1;tz=2*(z-lo[2])/(hi[2]-lo[2])-1
            px=np.polynomial.legendre.legvander(tx,3).ravel();pz=np.polynomial.legendre.legvander(tz,3).ravel()
            coeff+=np.array(q['Gxx_Gzz_div'])[:,None,None]*np.outer(px,pz)[None,:,:]*q['weight']/vol
        coeff*=np.outer([1,3,5,7],[1,3,5,7])[None,:,:]
        mode_rows.append(dict(geometry=geometry,region=qs[0]['region'],lower=lo.tolist(),upper=hi.tolist(),
            fields=['Gxx','Gzz','divergence'],Legendre_x_z=coeff.tolist(),
            scope='elementwise tensor Legendre coefficients, not a discontinuous replacement pressure space'))
    rows=m['observations'];w=np.array([r['weight'] for r in rows]);x=np.array([r['position'][0] for r in rows]);z=np.array([r['position'][2] for r in rows]);inner=np.array([r['inner'] for r in rows]);G=np.array([r['raw_G'] for r in rows]);a=np.column_stack([G[:,0],G[:,8],G[:,0]+G[:,8]])
    norms={}
    for scope,mask in dict(full=np.ones(len(rows),bool),inner=inner,outer=~inner,central=(x>=.25)&(x<=.75)).items():
        norms[scope]={}
        for j,f in enumerate(('Gxx','Gzz','divergence')):
            ids=np.flatnonzero(mask);idx=int(ids[np.argmax(abs(a[mask,j]))]);r=rows[idx]
            norms[scope][f]=dict(L2=float(np.sqrt(w[mask]@(a[mask,j]**2))),sampled_Linf=float(abs(a[idx,j])),signed_at_maximum=float(a[idx,j]),
                maximum={k:r[k] for k in ('position','geometry','element','requested_x_side','requested_normal_side','face_sides_xz')})
    loads=np.sum([r['signed_momentum_loads'] for r in d['near_h_loads']],axis=0)
    forcing=d['forcing_at_H'];force=np.array([r['actual_F_x'] for r in forcing]);parts=np.array([[r['mechanical_pressure_gradient_x'],r['conservative_tangential_convection_x']] for r in forcing])
    assert_allclose(force,parts.sum(axis=1),rtol=2e-10,atol=2e-13)
    return dict(schema='WJPassiveDiagnosisV1',state_sha256=d['state_sha256'],history_sha256=d['history_sha256'],
        pressure_projection=regions,element_modes=mode_rows,strong_components=norms,
        forcing=dict(range=[float(min(force)),float(max(force))],pressure_gradient_range=[float(min(parts[:,0])),float(max(parts[:,0]))],
            conservative_convection_range=[float(min(parts[:,1])),float(max(parts[:,1]))],samples=forcing),
        near_h_signed_loads=loads.tolist(),near_h_terms=d['near_h_terms'],near_h_test=d['near_h_test'],
        normal_trace_observations=sum(r['requested_normal_side']!=0 for r in rows),
        x_trace_observations=sum(r['requested_x_side']!=0 for r in rows),governing_solves=0)
