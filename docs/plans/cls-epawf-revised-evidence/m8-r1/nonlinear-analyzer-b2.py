#!/usr/bin/env python3
"""Independent analysis of an M8-R1 failed-stage diagnostic (no evolution).

Compares the instrumented accepted prefix to the immutable historical run,
recomputes physical mass dual norms and matrix actions, separates residual
exceptions from valid non-descent, and audits finite-difference directions.
The optional fresh-J replay remains a same-stage, unqualified experiment.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback

sys.dont_write_bytecode = True
THREADS = {key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    if isinstance(value, np.ndarray):
        value = value.tolist()
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def git(source, *args):
    return subprocess.check_output(['git','-c','gc.auto=0','-C',str(source),*args], text=True).strip()


class Audit:
    def __init__(self, raw, source):
        sys.path.insert(0, str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
        from cls_epawf_m7_energy_oracle import _fields
        self.fields, self.raw = _fields(raw), raw
        self.d, self.off, self.us = raw['dimension'], raw['offsets'], raw['ustar']
        self.K0, self.H, self.kappa = raw['parameters'][5], raw['parameters'][1], raw['parameters'][7]
        self.E0 = self.us**3/(self.kappa*self.H)
        self.mass = [field.mass() for field in self.fields]
        self.checks, self.errors, self._factors = 0, {}, {}
        self.gauge_variable_scale=self.us**2

    def require(self, condition, message):
        self.checks += 1
        if not condition:
            raise ValueError(message)

    def near(self, actual, expected, name, tolerance=1e-10):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        self.require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'Invalid '+name)
        error = float(np.max(np.abs(a-b)/np.maximum(1., np.abs(b)), initial=0.))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: componentwise scaled error {error} > {tolerance}')

    def constraints(self, entries):
        result = {int(row):float(value) for row,value in entries}
        self.require(len(result) == len(entries), 'Repeated constraint rows')
        return result

    def state_scale(self, block):
        return self.us if block < self.d else self.K0 if block == self.d else self.E0 if block == self.d+1 else self.us**2

    def scales(self):
        out = np.full(len(self.raw['state']), self.us**2)
        for block,field in enumerate(self.fields):
            out[self.off[block]:self.off[block]+field.size] = self.state_scale(block)
        out[-1]=self.gauge_variable_scale
        return out

    def norms(self, residual, constraints):
        residual = np.asarray(residual, float)
        values = []
        for block,field in enumerate(self.fields):
            free = np.array([i for i in range(field.size) if self.off[block]+i not in constraints], int)
            key = (block, tuple(free))
            if key not in self._factors:
                self._factors[key] = np.linalg.cholesky(self.mass[block][np.ix_(free,free)]) if len(free) else np.zeros((0,0))
            load = residual[self.off[block]+free]
            inverse = np.linalg.solve(self._factors[key], load) if len(free) else np.zeros(0)
            rate = self.us/self.H if block >= self.d+2 else self.state_scale(block)*self.us/self.H
            values.append(float(np.linalg.norm(inverse))/(rate*math.sqrt(field.volume)))
        values.append(abs(float(residual[-1]))/self.us**2)
        bc = 0.
        for row in constraints:
            block = next(b for b in range(len(self.fields)) if self.off[b] <= row < self.off[b]+self.fields[b].size)
            bc = max(bc, abs(float(residual[row]))/self.state_scale(block))
        values.append(bc)
        return np.asarray(values)

    def state_difference(self, a, b):
        delta = np.asarray(a)-np.asarray(b)
        result = []
        for block,field in enumerate(self.fields):
            x = delta[self.off[block]:self.off[block]+field.size]
            result.append(math.sqrt(max(0., float(x@self.mass[block]@x)))/(self.state_scale(block)*math.sqrt(field.volume)))
        result.append(abs(float(delta[-1]))/self.us**2)
        return result

    def norm_record(self, residual, constraints):
        values = self.norms(residual,constraints)
        return dict(blocks=values.tolist(), maximum=float(np.max(values)),
                    maximum_constraint=values[-1], gauge=values[-2])

    def linear(self, record, residual, constraints, jacobian=None):
        step = np.asarray(record['step'])
        self.near(record['maximum_dimensionless_step'], np.max(np.abs(step/self.scales())), 'dimensionless_step')
        observed = np.asarray(record['linear_residual'])
        if jacobian is not None:
            self.near(observed, jacobian@step+residual, 'independent_dense_linear_residual')
        norms = self.norms(observed,constraints)
        self.near(norms, record['linear_scaled_block_norms'], 'independent_linear_mass_dual_norms')
        return dict(dimensionless_step=float(np.max(np.abs(step/self.scales()))),
                    linear_residual=self.norm_record(observed,constraints),
                    dense_linear_product_independently_verified=jacobian is not None)

    def invalid_context(self, context, state):
        if not context or 'chart' not in context or 'local_k' not in context:
            return dict(native=context,independently_evaluated=False)
        values=np.asarray(state,float).copy()
        for row,value in self.raw['reference_constraints']:
            row=int(row)
            if row<self.off[self.d+1]:
                values[row]=value
        point=np.asarray(context['chart']);field=self.fields[self.d]
        e=field.element(point);ids,phi,_=field.rows(e,[point])
        c=values[self.off[self.d]+ids];constant=field.constant()[ids]
        centre=np.longdouble(c[0])/np.longdouble(constant[0])
        delta=c.astype(np.longdouble)-centre*constant.astype(np.longdouble)
        k=float(centre+phi[0].astype(np.longdouble)@delta)
        ordinary=float(phi[0]@c)
        bound=64*np.finfo(float).eps*(abs(float(centre))+float(np.abs(phi[0])@np.abs(np.asarray(delta,float))))
        self.near(k,context['local_k'],'independent_native_invalid_context_k')
        result=dict(native=context,independently_evaluated=True,k=k,ordinary_dot_k=ordinary,
            arithmetic_roundoff_scale=bound,magnitude_over_roundoff_scale=abs(k)/max(bound,np.finfo(float).tiny),
            actual_geometry=e['geometry'],interpretation='Exact exported guard coordinate with independent Jacobi rows; arithmetic scale is not an interval certificate')
        if context.get('sampled_K') is not None and context.get('ray_chart') is not None:
            sample=list(context['ray_chart']);sample[2]=self.raw['parameters'][3]
            sampled=field.evaluate(values,self.off[self.d],sample)[0]
            self.near(sampled,context['sampled_K'],'independent_native_invalid_context_sampled_K')
            result['sampled_K']=sampled
        if point[2]<=self.H+1e-12:
            def inner_details(coefficients):
                actual=np.asarray(coefficients,float).copy()
                for row,value in self.raw['reference_constraints']:
                    row=int(row)
                    if row<self.off[self.d+1]:
                        actual[row]=value
                local,unused=field.evaluate(actual,self.off[self.d],point)
                sample=point.copy();sample[2]=self.raw['parameters'][3]
                K=field.evaluate(actual,self.off[self.d],sample)[0]
                G=np.zeros((3,3))
                for block,axis in enumerate((0,2) if self.d==2 else (0,1,2)):
                    G[axis]=self.fields[block].evaluate(actual,self.off[block],point)[1]
                S=G+G.T;largest=float(np.linalg.eigvalsh(S)[-1])
                nu,unusedH,h=self.raw['parameters'][:3]
                root=math.sqrt(K)
                if h==0:
                    intercept=10.7
                else:
                    x=h*root/(70*nu);m=max(.5-.4*x**.7,1-.79*x**(-.28))
                    intercept=10.7*(1-x**m)
                zv=nu*intercept/root;nt=.2295*root*max(point[2]-zv,0.)
                covariance=2*local/3*np.eye(3)-nt*S
                return dict(k=local,sampled_K=K,k_over_K=local/K,nu_A=nt,G=G.tolist(),S=S.tolist(),
                    maximum_strain_eigenvalue=largest,strain_norm=float(np.linalg.norm(S)),
                    minimum_formal_covariance=float(np.linalg.eigvalsh(covariance)[0]),
                    r_R=3*nt*largest/(2*local) if local>0 else None,
                    formal_r_R=3*nt*largest/(2*local) if local!=0 else None,
                    transported_k_domain_valid=local>=0,
                    interpretation='Inner tensor/ratio at negative k is formal algebra outside the physical model domain' if local<0 else 'Inner AWF covariance at the same physical coordinate')
            result['current_inner']=inner_details(values)
            result['same_point_initial']=inner_details(self.raw['initial'])
            result['same_point_last_accepted']=inner_details(self.raw['accepted_states'][-1])
        return result


def physical_branches(base, perturbed):
    if base.get('status') != 'evaluated_source_quadrature' or perturbed.get('status') != 'evaluated_source_quadrature':
        return dict(evaluable=False, status=perturbed.get('status'), error=perturbed.get('message'))
    result = dict(evaluable=True)
    for key in ('outer_branch_pattern','cap_branch_pattern','limiter_branch_pattern'):
        a, b = base.get(key,''), perturbed.get(key,'')
        result[key] = dict(same_length=len(a)==len(b), identical=a==b,
            differing_indices=sum(x!=y for x,y in zip(a,b)) if len(a)==len(b) else None,
            interpretation='Fixed outer points/ground rays' if key!='limiter_branch_pattern' else
                'Literal adaptive-ray pattern only; physical normal points can move with current K/cutoff')
    for key in ('minimum_invariant_branch_gap','minimum_limiter_knot_gap','minimum_normalization_knot_gap'):
        result[key] = perturbed.get(key)
    return result


def iterations(records, audit, jacobian=None):
    output, counts = [], dict(valid_non_descent=0, accepted_descent=0, unevaluable_residual_or_norm=0)
    for record in records:
        constraints = audit.constraints(record['constraints'])
        residual = np.asarray(record['residual']); state = np.asarray(record['state'])
        norms = audit.norms(residual,constraints)
        audit.near(norms,record['scaled_block_norms'],'iteration_true_mass_dual_norms')
        lin = audit.linear(record['linear_step'],residual,constraints,jacobian)
        trials, accepted = [], False
        fractions = [trial['fraction'] for trial in record['trials']]
        audit.near(fractions,[2.**(-i) for i in range(len(fractions))],'unchanged_backtracking_fractions',0.)
        audit.require(len(fractions)<=21,'Unexpected line-search detail exceeds frozen policy')
        for trial in record['trials']:
            fraction = float(trial['fraction'])
            point = state+fraction*np.asarray(record['linear_step']['step'])
            if trial['evaluation_succeeded']:
                c = audit.constraints(trial['constraints']); r = np.asarray(trial['residual'])
                n = audit.norms(r,c)
                audit.near(n,trial['scaled_block_norms'],'trial_true_mass_dual_norms')
                audit.near(max(trial['scaled_block_norms']),trial['maximum_scaled_residual'],'trial_native_maximum')
                native_descent = trial['maximum_scaled_residual'] < max(record['scaled_block_norms'])
                audit.require(trial['classification']==('accepted_descent' if native_descent else 'valid_non_descent'),
                              'Stored line-search classification disagrees with unchanged native acceptance test')
                audit.require(not accepted,'Trials follow an already accepted fraction')
                accepted |= native_descent
                counts[trial['classification']] += 1
                for row,value in c.items():
                    audit.near(r[row],point[row]-value,'actual_constraint_rows')
                trials.append(dict(fraction=fraction,classification=trial['classification'],
                    norms=n.tolist(),independent_descent_margin=float(np.max(norms)-np.max(n)),
                    decision_uncertainty='Independent norm roundoff does not redefine the native strict descent test',
                    physical=trial.get('physical'),branches=physical_branches(record.get('physical',{}),trial.get('physical',{}))))
            else:
                audit.require(bool(trial.get('type')) and bool(trial.get('message')),'Unclassified residual exception')
                counts['unevaluable_residual_or_norm'] += 1
                trials.append(dict(fraction=fraction,classification='unevaluable_residual_or_norm',
                    exception_type=trial['type'],message=trial['message'],volume_exception=trial.get('volume_exception'),
                    invalid_context=audit.invalid_context(trial.get('invalid_context'),point),
                    inference='An exception alone is not proof of negative physical k or covariance'))
        audit.require(accepted==record['accepted'],'Iteration accepted flag disagrees with trial decisions')
        if not accepted:
            audit.require(len(fractions)==21,'A rejected nonlinear search did not export every frozen fraction')
        output.append(dict(iteration=record['iteration'],accepted=accepted,
            jacobian_age=record.get('jacobian_age'),jacobian_build=record.get('jacobian_build'),
            jacobian_build_time=record.get('jacobian_build_time'),rebuild_reason=record.get('rebuild_reason'),
            refreshed_flag=record.get('refreshed_flag'),residual_norms=norms.tolist(),linear=lin,
            physical=record.get('physical'),trials=trials))
    return dict(classification_counts=counts,iterations=output)


def fixed_limiter(engine, base_rays, state):
    """Compare limiter sides at common physical coordinates, not array indices."""
    engine.state = np.asarray(state)
    pattern = []
    gap = math.inf
    theta = engine.raw['robust_parameters'][2]
    for chart,points in base_rays:
        K = engine.evaluate(engine.d, [[chart[0],0.,engine.rK]])[0][0]
        k = engine.evaluate(engine.d,points)[0]
        differences = k/K-theta
        gap = min(gap,float(np.min(np.abs(differences))))
        pattern.extend(np.where(differences<0,'L','P').tolist())
    return ''.join(pattern),gap


def normalization_domain(trace, raw, source, session, audit):
    """Locate negative transported k on the actual trace/source ray rules.

    This reconstructs mathematical coefficient samples, not a new residual,
    line search, timestep or trajectory. A resulting admissible fraction is
    only a finite normalization-domain window; descent is not inferred.
    """
    import cls_epawf_m8_r1_archive as archive_module
    engine=archive_module.Reconstruction(raw,source,session.read_text())
    state=np.asarray(trace['state']);last=trace['iterations'][-1]
    step=np.asarray(last['linear_step']['step'])
    external={int(row):float(value) for row,value in raw['reference_constraints']
              if raw['offsets'][raw['dimension']]<=row<raw['offsets'][raw['dimension']+1]}
    field=engine.fields[engine.d]

    def lifted(values):
        result=np.asarray(values,float).copy()
        for row,value in external.items():
            result[row]=value
        return result

    def sampled(values,points):
        points=np.atleast_2d(points)
        centered=np.empty(len(points));plain=np.empty(len(points));bounds=np.empty(len(points))
        geometry=np.empty(len(points),int)
        groups={}
        for i,point in enumerate(points):
            groups.setdefault(field.element(point)['geometry'],[]).append(i)
        for eid,indices in groups.items():
            e=field.by_id[eid];ids,phi,_=field.rows(e,points[indices])
            c=np.asarray(values)[engine.off[engine.d]+ids]
            constant=engine.constants[engine.d][ids]
            centre=np.longdouble(c[0])/np.longdouble(constant[0])
            delta=c.astype(np.longdouble)-centre*constant.astype(np.longdouble)
            centred=centre+phi.astype(np.longdouble)@delta
            centered[indices]=centred.astype(float);plain[indices]=phi@c
            # An explicit arithmetic scale, not a rigorous interval bound on
            # native basis evaluation. It also exposes cancellation in k.
            bounds[indices]=64*np.finfo(float).eps*(abs(float(centre))+
                np.abs(phi)@np.abs(np.asarray(delta,float)))
            geometry[indices]=eid
        return centered,plain,bounds,geometry

    trace_ground=[]
    for e in engine.fields[engine.d+1].elements:
        if abs(e['lower'][2]-engine.H)>1e-12:
            continue
        nodes,_=engine.segment(e['degree'])
        trace_ground.extend([np.array([e['lower'][0]+(x+1)*(e['upper'][0]-e['lower'][0])/2,0.,engine.H])
                             for x in nodes])
    source_ground=[point for point,_ in engine.horizontal(engine.fields,engine.d,engine.H,raw['quadrature_orders'][0])]
    native_points=[]
    for e in field.elements:
        axes=[]
        for axis,count in zip((0,2),engine.native_counts):
            nodes=np.r_[-1.,np.polynomial.legendre.Legendre.basis(count-1).deriv().roots(),1.]
            axes.append(e['lower'][axis]+(nodes+1)*(e['upper'][axis]-e['lower'][axis])/2)
        native_points.extend([[x,0.,z] for z in axes[1] for x in axes[0]])
    native_points=np.asarray(native_points)
    base_rule_points=[]

    def inspect(values,label,collect=False):
        effective=lifted(values)
        answer=dict(label=label,groups={},profile_exceptions=[],
            maximum_external_k_coefficient_correction=float(np.max(np.abs(effective-np.asarray(values)))))
        all_nonnegative=True
        for kind,grounds,baseline in [('trace_baseline',trace_ground,True),
                                     ('trace_current',trace_ground,False),
                                     ('source_current',source_ground,False)]:
            result=dict(ray_count=len(grounds),samples=0,negative_count=0,minimum_k=math.inf,
                        worst_point=None,maximum_centering_difference=0.,minimum_k_to_roundoff_scale=None)
            for chart in grounds:
                view=dict(raw,upstream_transition=False,chi=0.) if baseline else raw
                try:
                    profile,metric=engine.profile(view,engine.fields,effective,chart)
                    r,_=engine.normal_rule(profile,view,engine.fields,metric)
                    points=np.column_stack([np.full_like(r,chart[0]),np.zeros_like(r),r])
                    k,ordinary,bounds,geometry=sampled(effective,points)
                except Exception as error:
                    all_nonnegative=False
                    answer['profile_exceptions'].append(dict(group=kind,ground=chart.tolist(),
                        type=type(error).__name__,message=str(error)))
                    continue
                if collect:
                    base_rule_points.extend(points.tolist())
                result['samples']+=len(points);result['negative_count']+=int(np.count_nonzero(k<0))
                result['maximum_centering_difference']=max(result['maximum_centering_difference'],float(np.max(np.abs(k-ordinary))))
                i=int(np.argmin(k))
                if k[i]<result['minimum_k']:
                    raw_k=sampled(np.asarray(values),points[i])[0][0]
                    result['minimum_k']=float(k[i])
                    result['worst_point']=dict(coordinate=points[i].tolist(),geometry=int(geometry[i]),
                        ground=chart.tolist(),sampled_K=float(profile.k),chi=float(profile.chi),
                        centered_long_double_k=float(k[i]),ordinary_dot_k=float(ordinary[i]),raw_unlifted_k=float(raw_k),
                        arithmetic_roundoff_scale=float(bounds[i]),
                        magnitude_over_roundoff_scale=abs(float(k[i]))/max(float(bounds[i]),np.finfo(float).tiny))
                all_nonnegative &= bool(np.all(k>=0))
            if not result['samples']:
                result['minimum_k']=None
            answer['groups'][kind]=result
        k,plain,bounds,geometry=sampled(effective,native_points)
        i=int(np.argmin(k))
        answer['native_tensor_GLL']=dict(samples=len(k),negative_count=int(np.count_nonzero(k<0)),minimum_k=float(k[i]),
            worst_point=dict(coordinate=native_points[i].tolist(),geometry=int(geometry[i]),
                centered_long_double_k=float(k[i]),ordinary_dot_k=float(plain[i]),
                arithmetic_roundoff_scale=float(bounds[i]),
                magnitude_over_roundoff_scale=abs(float(k[i]))/max(float(bounds[i]),np.finfo(float).tiny)))
        answer['finite_normalization_domain_nonnegative']=all_nonnegative
        return answer

    base=inspect(state,'failed_iteration_base',True)
    trials=[]
    for trial in last['trials']:
        record=inspect(state+trial['fraction']*step,'last_step_fraction_'+str(trial['fraction']))
        record['fraction']=trial['fraction'];record['native_evaluation_succeeded']=trial['evaluation_succeeded']
        trials.append(record)
    column=trace.get('assembly_failure') or {}
    fd=[]
    if 'column' in column:
        for sign in (1,-1):
            point=state.copy();point[int(column['column'])]+=sign*column['step']
            value=inspect(point,'failed_fresh_column_'+str(sign));value.update(column=column['column'],sign=sign,step=column['step'])
            fd.append(value)
    points=np.asarray(base_rule_points)
    effective=lifted(state);effective_step=lifted(state+step)-effective
    kb=sampled(effective,points)[0];dk=sampled(effective_step,points)[0]
    decreasing=dk<0
    estimates=np.where(decreasing,kb/np.where(decreasing,-dk,1.),math.inf)
    i=int(np.argmin(estimates));bound=float(estimates[i])
    window=dict(fixed_base_rule_upper_fraction=bound,limiting_base_coordinate=points[i].tolist(),
        base_k=float(kb[i]),step_delta_k=float(dk[i]),
        interpretation='Coefficient-linear estimate on fixed base quadrature points; moving physical panels are checked separately. '
            'No full residual/descent or accepted time map is asserted.',probes=[])
    if base['finite_normalization_domain_nonnegative'] and 0<bound<min(t['fraction'] for t in last['trials']):
        lo,hi=0.,min(1.,2*bound)
        high=inspect(state+hi*step,'normalization_window_upper_probe')
        window['probes'].append(dict(fraction=hi,observation=high))
        if not high['finite_normalization_domain_nonnegative']:
            for _ in range(12):
                fraction=(lo+hi)/2
                observation=inspect(state+fraction*step,'normalization_window_bisection')
                window['probes'].append(dict(fraction=fraction,observation=observation))
                if observation['finite_normalization_domain_nonnegative']:
                    lo=fraction
                else:
                    hi=fraction
            window['independent_strict_finite_domain_bracket']=[lo,hi]
            window['normalization_only']=True
    return dict(base=base,last_iteration_trials=trials,failed_fresh_column=fd,admissible_fraction_window=window,
        independent_helper=dict(path=str(Path(archive_module.__file__).resolve()),sha256=sha(archive_module.__file__)),
        interpretation='Constraints evaluates baseline/current epsilon trace on native trace-GLL wall coordinates before '
            'the source-volume Residual try/Refresh. Their normalization points differ from source-GL wall rays. '
            'Native tensor GLL is reported separately. Negative k at these points is not inferred from a covariance loss. '
            'The roundoff quantity is an arithmetic diagnostic, not a rigorous interval certificate.')


def directional(trace, raw, source, session, audit):
    cached, fresh = np.asarray(trace['cached_jacobian']),np.asarray(trace['fresh_jacobian'])
    state, base = np.asarray(trace['state']),np.asarray(trace['residual'])
    constraints = audit.constraints(trace['constraints'])
    expected_kinds = {'velocity','free_k','epsilon_with_trace','pressure_and_gauge','joint'}
    audit.require({item['kind'] for item in trace['directional_checks']}==expected_kinds,'Missing direction families')
    import cls_epawf_m8_r1_archive as archive_module
    Reconstruction = archive_module.Reconstruction
    engine = Reconstruction(raw,source,session.read_text())
    engine.state,engine.rays = state,{}
    base_rays = []
    for chart,_ in engine.horizontal(engine.fields,engine.d,engine.H,raw['quadrature_orders'][0]):
        base_rays.append((chart,engine.ray(chart)['points']))
    base_pattern,_ = fixed_limiter(engine,base_rays,state)
    # Independently reconstruct the nonlinear projected epsilon target. This
    # uses the current local k, physical sampler, seven sensors and capped ray
    # normalization, including dependence on the epsilon trace coefficients.
    engine.state,engine.rays = state,{}
    base_trace,_ = engine.trace()
    target_rows = np.array(sorted(base_trace),int)
    audit.require(len(target_rows)>0,'No nonlinear epsilon trace rows')
    audit.near([constraints[i] for i in target_rows],[base_trace[i] for i in target_rows],
               'independent_failed_state_epsilon_trace_targets')
    results = []
    for item in trace['directional_checks']:
        direction = np.asarray(item['direction'])
        audit.require(np.max(np.abs(direction))>0,'Empty exported direction')
        bubble=trace.get('admissible_dimensionless_diagnostic_columns') and item['kind'] in ('free_k','joint')
        if bubble:
            field=audit.fields[audit.d]
            for e in field.elements:
                points,_=field.quadrature(e,6);ids,phi,_=field.rows(e,points)
                x=(points[:,0]-field.lower[0])/(field.upper[0]-field.lower[0]);r=points[:,2]/audit.H
                target=np.where(points[:,2]<=audit.H,audit.K0*16*x*(1-x)*r*(1-r),0.)
                audit.near(phi@direction[audit.off[audit.d]+ids],target,'independent_positive_k_direction_bubble')
            other=np.ones(len(direction),bool);other[audit.off[audit.d]:audit.off[audit.d+1]]=False
            audit.require(np.max(np.abs(direction[other]/audit.scales()[other]))<=1.+1e-13,'Other direction components exceed declared scales')
        else:
            audit.require(np.max(np.abs(direction/audit.scales()))<=1.+1e-13,'Unscaled direction exceeds prescribed unit amplitude')
        ca,fa = cached@direction,fresh@direction
        audit.near(ca,item['cached_action'],'independent_cached_dense_action')
        audit.near(fa,item['fresh_action'],'independent_fresh_dense_action')
        audit.near([v['step'] for v in item['levels']],[2e-4,5e-5,1e-5,2e-6,5e-7],'five_level_FD_ladder',0.)
        levels,central_previous = [],None
        for level in item['levels']:
            h = float(level['step']);record=dict(step=h,sides={},central=None)
            values = {}
            for sign,label in ((1,'plus'),(-1,'minus')):
                value=level[label]
                if not value['evaluation_succeeded']:
                    record['sides'][label]=dict(evaluable=False,type=value.get('type'),message=value.get('message'),
                        volume_exception=value.get('volume_exception'),
                        invalid_context=audit.invalid_context(value.get('invalid_context'),state+sign*h*direction))
                    continue
                point = state+sign*h*direction
                r = np.asarray(value['residual']);c = audit.constraints(value['constraints'])
                audit.require(set(c)==set(constraints),'FD perturbation changed constrained row membership')
                action = (r-base)/(sign*h)
                values[label] = r
                pattern,gap = fixed_limiter(engine,base_rays,point)
                row = dict(evaluable=True,action=audit.norm_record(action,constraints),
                    cached_error=audit.norm_record(action-ca,constraints),fresh_error=audit.norm_record(action-fa,constraints),
                    branches=physical_branches(trace.get('physical',{}),value.get('physical',{})),
                    fixed_physical_limiter_changed=sum(a!=b for a,b in zip(base_pattern,pattern)),
                    fixed_physical_limiter_minimum_gap=gap)
                for index,target in c.items():
                    audit.near(r[index],point[index]-target,'FD_full_constraint_rows')
                if item['kind'] in ('free_k','epsilon_with_trace'):
                    engine.state,engine.rays = point,{}
                    independent_trace,_ = engine.trace()
                    audit.near([c[i] for i in target_rows],[independent_trace[i] for i in target_rows],
                               'independent_perturbed_epsilon_trace_targets')
                    expected = direction[target_rows]-np.array([(independent_trace[i]-base_trace[i])/(sign*h) for i in target_rows])
                    audit.near(action[target_rows],expected,'independent_dependent_epsilon_one_sided_actions',2e-9)
                    row['independent_matching_trace_action_maximum']=float(np.max(np.abs(expected)))
                    row['omitted_trace_dependency_control']=float(np.max(np.abs(expected-direction[target_rows])))
                record['sides'][label]=row
            if 'plus' in values and 'minus' in values:
                action = (values['plus']-values['minus'])/(2*h)
                reference = max(1.,max(audit.norms(fa,constraints)))
                record['central']=dict(action=audit.norm_record(action,constraints),
                    cached_error=audit.norm_record(action-ca,constraints),fresh_error=audit.norm_record(action-fa,constraints),
                    relative_cached_error=float(max(audit.norms(action-ca,constraints))/reference),
                    relative_fresh_error=float(max(audit.norms(action-fa,constraints))/reference),
                    drift_from_previous=None if central_previous is None else audit.norm_record(action-central_previous,constraints))
                central_previous=action
            levels.append(record)
        results.append(dict(kind=item['kind'],maximum_dimensionless_direction=float(np.max(np.abs(direction/audit.scales()))),
            k_direction='independently verified positive physical bubble' if bubble else 'original coefficient direction',
            cached_fresh_difference=audit.norm_record(ca-fa,constraints),levels=levels))
    accuracy=[]
    for item in results:
        central=item['levels'][-1]['central']
        accuracy.append(dict(kind=item['kind'],finest_step=item['levels'][-1]['step'],
            finest_relative_fresh_error=None if central is None else central['relative_fresh_error'],
            diagnostic_reporting_threshold=1e-7,
            threshold_met=central is not None and central['relative_fresh_error']<=1e-7,
            interpretation='Fixed diagnostic comparison threshold; any failure remains an observation, not a tolerance change or proof of a wrong smooth tangent at a switching state'))
    return dict(families=results,independent_nonlinear_trace_reconstruction=True,accuracy_observations=accuracy,
        accuracy_failures=[item for item in accuracy if not item['threshold_met']],
        fixed_limiter_samples=len(base_pattern),trace_rows=len(target_rows),
        independent_trace_helper=dict(path=str(Path(archive_module.__file__).resolve()),sha256=sha(archive_module.__file__)),
        interpretation='Native residual samples, independently recomputed matrix actions and physical mass dual norms; '
            'free-k/epsilon trace targets are independently reconstructed using indexed Jacobi/AWF/sensor algebra. '
            'Finite differences at max/limiter/cap switches are reported as one-sided observations, not a smoothness theorem.')


def column_metadata(trace,audit):
    outputs={}
    for name in ('diagnostic_jacobian_columns','last_retry_jacobian_columns'):
        records=trace.get(name)
        if records is None:
            continue
        audit.require(len(records)==len(trace['state']),'Incomplete diagnostic Jacobian column metadata')
        counts={};rejected=[];shrunk=[]
        base_state=trace['state']
        if name=='last_retry_jacobian_columns' and trace.get('fresh_jacobian_retry',{}).get('iterations'):
            base_state=trace['fresh_jacobian_retry']['iterations'][-1]['state']
        for expected,record in enumerate(records):
            audit.require(record['column']==expected,'Diagnostic column metadata is out of order')
            stencil=record['stencil'];counts[stencil]=counts.get(stencil,0)+1
            attempts=record.get('attempts',[])
            if not attempts:
                audit.require(stencil=='lifted_constraint_identity','Missing coordinate-stencil attempts')
                continue
            audit.near([attempt['step']/attempts[0]['step'] for attempt in attempts],
                       [1.,.25,.0625,.015625,.00390625][:len(attempts)],'declared_diagnostic_FD_shrinks',1e-14)
            if len(attempts)>1:
                shrunk.append(dict(column=expected,attempt_count=len(attempts),selected_step=record.get('step'),stencil=stencil))
            for attempt in attempts:
                for sign,label in ((1,'plus'),(-1,'minus')):
                    if not attempt[label+'_valid']:
                        exception=attempt[label+'_exception'];point=np.asarray(base_state).copy();point[expected]+=sign*attempt['step']
                        rejected.append(dict(column=expected,step=attempt['step'],sign=sign,type=exception['type'],
                            message=exception['message'],invalid_context=audit.invalid_context(exception.get('invalid_context'),point)))
            if stencil=='central':
                audit.require(attempts[-1]['plus_valid'] and attempts[-1]['minus_valid'],'Central diagnostic stencil has an invalid side')
        outputs[name]=dict(columns=len(records),stencils=counts,shrunk_columns=shrunk,rejected_coordinate_probes=rejected,
            original_canonical_sha256=canonical(records),scope=trace.get('diagnostic_column_admissibility_scope'))
    return outputs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--trace',required=True,type=Path)
    parser.add_argument('--run-report',required=True,type=Path)
    parser.add_argument('--archive',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--runtime-pin',required=True,type=Path)
    parser.add_argument('--session',type=Path)
    parser.add_argument('--reuse-normalization',type=Path)
    args=parser.parse_args()
    if platform.system()!='Linux':
        raise SystemExit('Scientific nonlinear analysis must run on the designated Linux host')
    args.output.mkdir(parents=True,exist_ok=False)
    report=dict(schema=1,milestone='M8-R1',phase='B nonlinear diagnostic analysis',started_utc=stamp(),
        success=False,solver_launched=False,numerical_threads=THREADS,driver_sha256=sha(__file__),failures=[])
    out=args.output/'validation.json'
    audit=None
    try:
        report['source_commit']=git(args.source,'rev-parse','HEAD')
        report['source_clean_start']=not bool(git(args.source,'status','--porcelain'))
        report['inputs']={str(p):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in
            [args.trace,args.run_report,args.archive]+([args.runtime_pin] if args.runtime_pin else [])}
        raw=json.loads(args.run_report.read_text());old=json.loads(args.archive.read_text());trace=json.loads(args.trace.read_text())
        runtime=json.loads(args.runtime_pin.read_text())
        report['producing_runtime']=dict(source_commit=runtime['source_commit'],
            source_root=runtime['source_root'],build_root=runtime['build_root'],
            solver_sha256=runtime['runtime_artifacts']['build/IncNavierStokesSolver'],
            full_pin_path=str(args.runtime_pin),full_pin_sha256=sha(args.runtime_pin),
            interpretation='Actual instrumented producing runtime, distinct from archived AA and analysis-driver source')
        audit=Audit(raw,args.source)
        audit.gauge_variable_scale=audit.us/audit.H if trace.get('admissible_dimensionless_diagnostic_columns') else audit.us**2
        report['diagnostic_variable_scaling']=dict(gauge_multiplier_scale=audit.gauge_variable_scale,
            gauge_variable_units='inverse time' if trace.get('admissible_dimensionless_diagnostic_columns') else 'historical B1 diagnostic pressure scale',
            residual_gauge_scale=audit.us**2,
            interpretation='B2 corrects diagnostic variable units only; the native physical residual and gauge residual normalization are unchanged')
        audit.require(report['source_clean_start'],'Nonlinear analysis source must be clean/pinned')
        audit.require(trace['schema']==1 and trace['kind']=='first_failed_nonlinear_stage' and
            trace['diagnostic_unqualified'] and trace['original_policy_unchanged'],'Wrong diagnostic trace contract')
        for key in ('dimension','offsets','fields','parameters','robust_parameters','transition_controls','transition_weights'):
            audit.require(raw[key]==old[key],'Instrumented accepted-prefix equation/space mismatch: '+key)
        audit.near(raw['audit_times'],old['audit_times'],'archived_accepted_times',0.)
        audit.require(len(raw['accepted_states'])==len(old['accepted_states']),'Instrumented prefix length differs')
        prefix=[]
        for time,a,b in zip(raw['audit_times'],raw['accepted_states'],old['accepted_states']):
            difference=audit.state_difference(a,b)
            audit.require(max(difference)<=1e-10,'Instrumented prefix changed actual fields/pressures/gauge')
            prefix.append(dict(time=time,physical_state_difference=difference))
        for a,b in zip(raw['native_history'],old['native_history']):
            audit.near(a,b,'unchanged_complete_native_history',1e-10)
        audit.near(raw['native_history_times'],old['native_history_times'],'unchanged_native_history_times',0.)
        audit.near(raw['state'],old['state'],'unchanged_last_accepted_coefficients',1e-10)
        report['accepted_prefix_equivalence']=dict(epochs=prefix,all_fields_pressures_gauge=True,complete_native_history=True,
            interpretation='Instrumented failure reproduction preserves the archived accepted prefix; fresh-J retry does not advance it')
        audit.near(trace['offsets'],raw['offsets'],'diagnostic_native_layout',0.)
        audit.near(trace['dimensionless_variable_scales'],audit.scales(),'independent_variable_scales',1e-14)
        audit.near(trace['time'],old['time']+old['dt'],'same_first_failed_stage_time',1e-14)
        audit.near(trace['alpha'],1.5/old['dt'],'same_IMEX2_alpha',1e-14)
        differential=raw['offsets'][raw['dimension']+2]
        history=(4*np.asarray(old['accepted_states'][-1])[:differential]-np.asarray(old['accepted_states'][-2])[:differential])/3
        audit.near(np.asarray(trace['history'])[:differential],history,'independent_failed_stage_IMEX2_history')
        report['stage']=dict(time=trace['time'],alpha=trace['alpha'],tolerance=trace['tolerance'],
            failure=trace['failure'],failure_type=trace['failure_type'],assembly_failure=trace.get('assembly_failure'),
            entry_jacobian_builds=trace['stage_entry_builds'],entry_jacobian_age=trace['stage_entry_jacobian_age'],
            cached_jacobian_time=trace['cached_jacobian_time'],cached_jacobian_iteration=trace['cached_jacobian_iteration'],
            failed_state_sha256=canonical(trace['state']),history_sha256=canonical(trace['history']),
            cached_jacobian_sha256=canonical(trace['cached_jacobian']),physical=trace.get('physical'))
        report['original_policy']=iterations(trace['iterations'],audit)
        report['diagnostic_column_metadata']=column_metadata(trace,audit)
        if trace['base_evaluation_succeeded']:
            c=audit.constraints(trace['constraints'])
            audit.near(audit.norms(trace['residual'],c),trace['scaled_block_norms'],'failed_base_true_scaled_norms')
            report['stage']['base_norms']=audit.norm_record(trace['residual'],c)
            cached=np.asarray(trace['cached_jacobian'])
            for record in trace['iterations'][-1:]:
                report['original_policy']['iterations'][-1]['linear']=audit.linear(record['linear_step'],
                    np.asarray(record['residual']),audit.constraints(record['constraints']),cached)
            session=args.session or args.run_report.with_name('case.xml')
            report['inputs'][str(session)]=dict(bytes=session.stat().st_size,sha256=sha(session))
            if trace['iterations']:
                if args.reuse_normalization:
                    prior=json.loads(args.reuse_normalization.read_text())
                    prior_trace_path=next(Path(path) for path in prior['inputs'] if path.endswith('case.nonlinear.json'))
                    audit.require(sha(prior_trace_path)==prior['inputs'][str(prior_trace_path)]['sha256'],'Prior trace bytes changed')
                    prior_trace=json.loads(prior_trace_path.read_text())
                    for key in ('state','history','alpha','time'):
                        audit.require(trace[key]==prior_trace[key],'Normalization reuse changed failed-state inputs: '+key)
                    audit.require(trace['iterations'][-1]['linear_step']['step']==prior_trace['iterations'][-1]['linear_step']['step'],
                                  'Normalization reuse changed original Newton direction')
                    for key in ('parameters','robust_parameters','transition_controls','fields'):
                        prior_raw_path=next(Path(path) for path in prior['inputs'] if path.endswith('case.evolution.json') and '/m8-r1/' in path)
                        prior_raw=json.loads(prior_raw_path.read_text())
                        audit.require(raw[key]==prior_raw[key],'Normalization reuse changed physical inputs: '+key)
                    old_source=Path(prior['producing_runtime']['source_root'])
                    math_inputs=0
                    for path,digest in prior['imported_source_files'].items():
                        relative=Path(path).relative_to(old_source)
                        audit.require(sha(args.source/relative)==digest,'Normalization reuse changed independent numerical module: '+str(relative))
                        math_inputs+=1
                    report['normalization_domain']=prior['normalization_domain']
                    report['normalization_reuse']=dict(path=str(args.reuse_normalization),sha256=sha(args.reuse_normalization),
                        exact_base_history_original_Newton_step_and_physical_inputs_equal=True,
                        independent_numerical_modules_identical=math_inputs,
                        interpretation='No new trajectory or changed equation inferred; native B2 guard points are verified independently above')
                else:
                    report['normalization_domain']=normalization_domain(trace,raw,args.source,session,audit)
                save(out,report)
            if trace.get('fresh_jacobian') is not None and trace.get('directional_checks') is not None:
                save(out,report)
                report['directional_checks']=directional(trace,raw,args.source,session,audit)
                report['observed_Jacobian_accuracy_failures']=report['directional_checks']['accuracy_failures']
        else:
            report['stage']['base_exception']=trace['base_exception']
        retry=trace.get('fresh_jacobian_retry')
        if retry:
            audit.require(retry['diagnostic_unqualified'] and not retry['history_advanced'],'Retry promoted or advanced native history')
            audit.require(retry.get('maximum_iterations',12)<=12 and len(retry.get('iterations',[]))<=12,'Retry exceeded bounded12 iterations')
            result=iterations(retry.get('iterations',[]),audit)
            if result['iterations'] and trace.get('fresh_jacobian') is not None:
                first=retry['iterations'][0]
                result['iterations'][0]['linear']=audit.linear(first['linear_step'],np.asarray(first['residual']),
                    audit.constraints(first['constraints']),np.asarray(trace['fresh_jacobian']))
            result.update(terminal=retry['terminal'],converged=retry.get('converged',False),history_advanced=False,
                diagnostic_unqualified=True,physical=retry.get('physical'),exception_type=retry.get('type'),message=retry.get('message'))
            if 'residual' in retry:
                c=audit.constraints(retry['constraints']);norms=audit.norms(retry['residual'],c)
                audit.near(norms,retry['scaled_block_norms'],'retry_terminal_true_scaled_norms')
                result['terminal_norms']=norms.tolist()
                audit.require(retry['converged']==(max(retry['scaled_block_norms'])<trace['tolerance']),
                              'Retry convergence flag disagrees with frozen nonlinear tolerance')
            report['fresh_jacobian_same_stage_retry']=result
        report['interpretation']=dict(physical_model_qualified=False,
            accepted_prefix_already_unqualified='Phase A establishes negative inner covariance and local outer-box exit at t0',
            Newton_cause='Assess original trial classifications and same-stage fresh-J comparison; no assumed refresh causality',
            no_time_history_advance=True,no_new_trajectory_from_this_analyzer=True,
            analysis_success_scope='Observation/provenance/norm/constraint consistency; numerical Jacobian accuracy failures, if present, remain explicitly recorded')
        report['source_commit_end']=git(args.source,'rev-parse','HEAD')
        report['source_clean_end']=not bool(git(args.source,'status','--porcelain'))
        audit.require(report['source_commit_end']==report['source_commit'] and report['source_clean_end'],'Analysis source changed')
        report['imported_source_files']={str(Path(module.__file__).resolve()):sha(module.__file__) for module in list(sys.modules.values())
            if getattr(module,'__file__',None) and str(Path(module.__file__).resolve()).startswith(str(args.source))}
        report['checks'],report['errors'],report['success']=audit.checks,audit.errors,True
        report['finished_utc']=stamp();save(out,report)
        print(json.dumps(dict(success=True,checks=audit.checks,output=str(out))))
        return 0
    except Exception as error:
        report['failures'].append(dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
        if audit is not None:
            report['checks'],report['errors']=audit.checks,audit.errors
        report['finished_utc']=stamp();save(out,report)
        print(traceback.format_exc(),file=sys.stderr)
        return 1


if __name__=='__main__':
    raise SystemExit(main())
