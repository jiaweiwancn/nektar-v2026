# SPDX-License-Identifier: MIT
"""Read-only physical-term, element and common-test analysis of native packets."""
import hashlib
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.sparse import csr_matrix

FIELDS=('u','w','k','epsilon_outer','Pi_inner','Pi_outer')

def basis_matrix(rows,n,derivative=None):
    data=[];indices=[];ptr=[0]
    for row in rows:
        data.extend(row['phi'] if derivative is None else row['derivative'][derivative]);indices.extend(row['ids']);ptr.append(len(indices))
    return csr_matrix((data,indices,ptr),shape=(len(rows),n))

def scales(ref):
    return np.array([ref.us,ref.us,ref.K,ref.A/ref.H,1.,1.])*ref.us/ref.H

def analyze(packet,ref):
    o=packet['observation'];r=packet['resolution'];n=sum(o['sizes']);off=o['offsets'];size=o['sizes'];fixed={int(x[0]) for x in r['constraints']};terms={k:np.array(v) for k,v in r['terms'].items()};native={k:np.array(v) for k,v in r['native'].items()}
    names=r['physical_term_names'];physical=sum((terms[k] for k in names),np.zeros(n));bound=sum(np.abs(terms[k]) for k in names);tol=2e-11*max(1.,float(np.max(bound)))
    constrained=physical.copy()
    for i,value in r['constraints']:constrained[int(i)]=o['state'][int(i)]-value
    defects=dict(physical_sum_vs_native=float(np.max(np.abs(physical-native['spatial_unmasked']))),constraint_replacement_vs_native=float(np.max(np.abs(constrained-native['stage_constrained']))),constraint_rows_vs_export=float(np.max(np.abs(constrained-terms['total_constrained']))))
    if any(v>tol for v in defects.values()):raise ValueError('term reconstruction mismatch '+str(defects))
    masses=[np.array(m) for m in r['mass']];rates=scales(ref);volumes=[8.,8.,8.,4.,4.,4.];factors=[];freeblocks=[]
    for b in range(6):
        free=np.array([i for i in range(size[b]) if off[b]+i not in fixed],int);freeblocks.append(free);factors.append(cho_factor(masses[b][np.ix_(free,free)],lower=True))
    def norm(v,b,ids=None):
        if ids is None:ids=freeblocks[b];factor=factors[b]
        else:factor=cho_factor(masses[b][np.ix_(ids,ids)],lower=True)
        a=np.asarray(v)[off[b]+ids];sq=float(a@cho_solve(factor,a));return float(np.sqrt(max(0,sq))/(rates[b]*np.sqrt(volumes[b])))
    norms={FIELDS[b]:norm(native['spatial_unmasked'],b) for b in range(6)}
    block_terms={FIELDS[b]:{k:norm(terms[k],b) for k in names} for b in range(6)}
    support=[set(),set()]
    for p in o['volume']:support[p['region']].update(p['velocity_k']['ids'])
    shared=support[0]&support[1]
    partitions={name:np.array(sorted([i for i in freeblocks[2] if (i in shared if name=='shared' else i in support[region] and i not in shared)]),int) for name,region in [('inner',0),('outer',1),('shared',0)]}
    kpart={name:norm(native['spatial_unmasked'],2,ids) for name,ids in partitions.items()}
    local=[];sumlocal=np.zeros(n)
    for key,et in r['element_terms'].items():
        region,geometry=key.split(':');element=next(e for e in o['elements'] if e['geometry']==int(geometry));v=sum((np.array(et[k]) for k in names if k in et),np.zeros(n));sumlocal+=v
        row=dict(element=key,interval=[element['lower'][2],element['upper'][2]],bounds=[element['lower'],element['upper']],free_mass_dual={FIELDS[b]:norm(v,b) for b in range(4)},free_row_l2={},constrained_row_l2={},terms={})
        for b in range(4):
            allids=np.arange(off[b],off[b]+size[b]);free=[i for i in allids if i not in fixed];bc=[i for i in allids if i in fixed]
            row['free_row_l2'][FIELDS[b]]=float(np.linalg.norm(v[free]));row['constrained_row_l2'][FIELDS[b]]=float(np.linalg.norm(v[bc]))
            row['terms'][FIELDS[b]]={k:dict(free_mass_dual=norm(np.array(value),b),signed_free_row_sum=float(np.sum(np.asarray(value)[free])),signed_reaction_row_sum=float(np.sum(np.asarray(value)[bc]))) for k,value in et.items() if k in names}
        local.append(row)
    # Body and gauge are not volume element primitives; account for them explicitly.
    localdef=float(np.max(np.abs(sumlocal+terms['gauge']+terms['prescribed_body']-physical)))
    if localdef>tol:raise ValueError('element sum mismatch '+str(localdef))
    faces={};jumps=[]
    for face in r['one_sided_faces']:
        q=tuple(face['coordinate']);faces.setdefault(q,[]).append(face)
    for q,rows in faces.items():
        if len(rows)!=2:continue
        rows.sort(key=lambda x:x['side']);a,b=rows
        jump=dict(coordinate=q,geometries=[a['geometry'],b['geometry']],sides=[a['side'],b['side']],interface=abs(q[2]-ref.H)<1e-12,
            velocity=(np.array(a['velocity'])-b['velocity']).tolist(),k=a['k']-b['k'],Pi=a['Pi']-b['Pi'],
            traction=(np.array(a['full_traction'])-b['full_traction']).tolist(),k_flux=a['k_conservative_flux']-b['k_conservative_flux'],
            epsilon_flux=None if a['epsilon_conservative_flux'] is None or b['epsilon_conservative_flux'] is None else a['epsilon_conservative_flux']-b['epsilon_conservative_flux'],
            one_sided=rows)
        jumps.append(jump)
    support_rows={}
    for b,name in enumerate(FIELDS[:4]):
        which='epsilon' if b==3 else 'velocity_k';sp={}
        for q in o['volume']:
            if b==3 and q['region']!=1:continue
            for i in q[which]['ids']:sp.setdefault(i,set()).add(q[which]['geometry'])
        ranked=[]
        for i in freeblocks[b]:
            value=float(native['spatial_unmasked'][off[b]+i]);ranked.append(dict(row=int(off[b]+i),load=value,diagonal_mass_scaled_load=value/np.sqrt(masses[b][i,i])/rates[b],support=sorted(sp.get(i,set()))))
        support_rows[name]=sorted(ranked,key=lambda x:abs(x['diagonal_mass_scaled_load']),reverse=True)[:80]
    return dict(status='pass',alpha=0.,dominant_native_rows_with_basis_support=support_rows,accumulation_tolerance=tol,reconstruction_defects=defects,element_reconstruction_defect=localdef,
        native_full_mass_dual=norms,native_k_partition_mass_dual=kpart,term_mass_dual=block_terms,
        element_localization=sorted(local,key=lambda x:x['free_mass_dual']['u'],reverse=True),one_sided_jumps=jumps,
        constrained_reactions={FIELDS[b]:dict(rows=[int(i) for i in range(off[b],off[b]+size[b]) if i in fixed],load=[float(physical[i]) for i in range(off[b],off[b]+size[b]) if i in fixed]) for b in range(6)},
        scope='native alpha-zero terms and reactions; full native test space; no time-discrete mass or semidiscrete scalar-rate inference')

def common_test(coarse,fine,cross,ref):
    """P = M_f^{-1} B_fc on native quadrature, then restrict by P^T.

    cross is coarse native test basis evaluated at fine volume integration rows.
    Coarse/fine need not be nested if a shared test subspace is selected explicitly.
    This routine requires exact polynomial inclusion and checks it before use.
    """
    co,fo=coarse['observation'],fine['observation'];cr,fr=coarse['resolution'],fine['resolution'];rows=fo['volume'];test=cross['test_basis'];assert len(rows)==len(test)
    cf={int(i) for i,v in cr['constraints']};ff={int(i) for i,v in fr['constraints']};result={};rates=scales(ref)
    for b,name in enumerate(FIELDS):
        region=None if b in (0,1,2) else 1 if b in (3,5) else 0;which='velocity_k' if b<3 else 'epsilon' if b==3 else 'pressure'
        mask=[i for i,p in enumerate(rows) if region is None or p['region']==region]
        rf=[rows[i][which] for i in mask];rc=[test[i][which] for i in mask];w=np.array([rows[i]['weight'] for i in mask])
        nf,nc=fo['sizes'][b],co['sizes'][b];vf=basis_matrix(rf,nf);vc=basis_matrix(rc,nc)
        mf=np.array(fr['mass'][b]);mc=np.array(cr['mass'][b]);B=(vf.T@vc.multiply(w[:,None])).toarray();P=cho_solve(cho_factor(mf,lower=True),B)
        ve=float(np.max(np.abs(vf@P-vc.toarray())));ge=0.
        for a in (0,2):ge=max(ge,float(np.max(np.abs(basis_matrix(rf,nf,a)@P-basis_matrix(rc,nc,a).toarray()))))
        me=float(np.max(np.abs(P.T@mf@P-mc)));den=max(1.,float(np.max(np.abs(mc))));condition=dict(value_error=ve,gradient_error=ge,mass_congruence_error=me,mass_scale=den)
        free=np.array([j for j in range(nc) if co['offsets'][b]+j not in cf],int);fixedf=np.array([j for j in range(nf) if fo['offsets'][b]+j in ff],int)
        trace=float(np.max(np.abs(P[np.ix_(fixedf,free)]),initial=0.));condition['free_test_on_fine_fixed_rows']=trace
        if ve>2e-10 or ge>2e-6 or me>2e-11*den or trace>2e-10:raise ValueError('test inclusion/orientation/constraint check failed '+name+str(condition))
        p=P[:,free];M=mc[np.ix_(free,free)];factor=cho_factor(M,lower=True);vol=8. if b<3 else 4.;scale=rates[b]*np.sqrt(vol)
        lo=np.array(fr['native']['spatial_unmasked'])[fo['offsets'][b]:fo['offsets'][b]+nf];restricted=p.T@lo
        own=np.array(cr['native']['spatial_unmasked'])[co['offsets'][b]+free]
        norm=lambda x:float(np.sqrt(max(0,x@cho_solve(factor,x)))/scale)
        terms={}
        for k in fr['physical_term_names']:
            value=p.T@np.array(fr['terms'][k])[fo['offsets'][b]:fo['offsets'][b]+nf];terms[k]=norm(value)
        result[name]=dict(prolongation_sha256=hashlib.sha256(np.asarray(P,dtype='<f8').tobytes()).hexdigest(),identity_checks=condition,coarse_common_mass_dual=norm(own),fine_restricted_common_mass_dual=norm(restricted),restricted_load_difference=norm(restricted-own),restricted_term_mass_dual=terms,common_free_test_count=len(free),fine_full_count=nf)
    return dict(status='pass',fields=result,definition='P^T r_f on identical coarse free tests and physical coarse mass/normalization; full fine norms retained separately')
