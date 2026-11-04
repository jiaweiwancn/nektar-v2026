"""Independent complex VALUE operator for the flat R3 action probes.

Only native coefficient maps, states and declared generalized-branch labels
are inputs. No production derivative or production residual is evaluated.
Numerical invocation belongs exclusively to the authoritative Linux host.
"""
import os
import sys
from functools import lru_cache

sys.dont_write_bytecode = True
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
              'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import numpy as np
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m4 import spaces as real_spaces
from cls_epawf_m6_operator_oracle import _segment_projection

LD, CD = np.longdouble, np.clongdouble


def require(value, message):
    if not value:
        raise ValueError(message)


def maximum(a, b):
    return np.where(np.real(a) >= np.real(b), a, b)


def minimum(a, b):
    return np.where(np.real(a) <= np.real(b), a, b)


def positive(a):
    return np.where(np.real(a) > 0, a, CD(0))


def rounded_base(value):
    """Base panel coordinates are binary64; keep their complex displacement."""
    return np.asarray(np.asarray(np.real(value), dtype=float), dtype=LD)+CD(1j)*np.imag(value)


@lru_cache(maxsize=16)
def gauss(order):
    return tuple(np.asarray(value, dtype=LD) for value in np.polynomial.legendre.leggauss(order))


def basis(degree, x):
    """Explicit Modified_A p1/p3 polynomials at COMPLEX coordinates."""
    require(degree in (1, 3), 'R3 oracle is deliberately restricted to the frozen p1/p3 spaces')
    x = np.atleast_1d(np.asarray(x, dtype=CD))
    p, d = [(1-x)/2, (1+x)/2], [np.full_like(x, -.5), np.full_like(x, .5)]
    if degree == 3:
        p.extend([(1-x*x)/4, x*(1-x*x)/2])
        d.extend([-x/2, (1-3*x*x)/2])
    return np.asarray(p).T, np.asarray(d).T


class Field:
    def __init__(self, record):
        self.size, self.elements = int(record['size']), record['mesh']
        self.by_id = {int(e['geometry']):e for e in self.elements}

    def element(self, points, side=0):
        point = np.real(np.asarray(points))
        candidates = []
        for e in self.elements:
            if all(e['lower'][a]-1e-11 <= point[a] <= e['upper'][a]+1e-11 for a in (0, 2)):
                if side > 0 and point[2] >= e['upper'][2]-1e-12 or side < 0 and point[2] <= e['lower'][2]+1e-12:
                    continue
                candidates.append(e)
        require(bool(candidates), 'Independent complex row missed its base element')
        return min(candidates, key=lambda e:e['geometry'])

    def rows(self, element, points):
        points = np.atleast_2d(np.asarray(points, dtype=CD))
        x0, r0 = (LD(element['lower'][a]) for a in (0, 2))
        hx, hr = (LD(element['upper'][a])-LD(element['lower'][a]) for a in (0, 2))
        px, dx = basis(element['degree'], 2*(points[:, 0]-x0)/hx-1)
        pr, dr = basis(element['degree'], 2*(points[:, 2]-r0)/hr-1)
        signs = np.asarray(element['sign'], dtype=LD)
        shape = (len(points), -1)
        phi = (pr[:, :, None]*px[:, None, :]).reshape(shape)*signs
        gradient = np.zeros((len(points), 3, phi.shape[1]), dtype=CD)
        gradient[:, 0] = (pr[:, :, None]*dx[:, None, :]).reshape(shape)*signs*2/hx
        gradient[:, 2] = (dr[:, :, None]*px[:, None, :]).reshape(shape)*signs*2/hr
        degree = element['degree']
        constant = np.asarray([signs[i] if i % (degree+1) < 2 and i//(degree+1) < 2 else 0
            for i in range(len(signs))], dtype=LD)
        return np.asarray(element['ids'], dtype=int), phi, gradient, constant

    def values(self, state, offset, points, geometry=None, side=0):
        points = np.atleast_2d(np.asarray(points, dtype=CD))
        if geometry is None:
            groups = {}
            for i, point in enumerate(points): groups.setdefault(self.element(point, side)['geometry'], []).append(i)
            values, gradient = np.empty(len(points), dtype=CD), np.empty((len(points), 3), dtype=CD)
            for key, indices in groups.items():
                value = self.values(state, offset, points[indices], key)
                values[indices], gradient[indices] = value[3], value[4]
            return None, None, None, values, gradient
        e = self.by_id[int(geometry)]; ids, phi, grad, constant = self.rows(e, points)
        coefficients = state[offset+ids]
        centre = coefficients[0]/constant[0]
        difference = coefficients-centre*constant
        values = centre+np.sum(phi*difference[None, :], axis=1, dtype=CD)
        derivatives = np.sum(grad*difference[None, None, :], axis=2, dtype=CD)
        return ids, phi, grad, values, derivatives

    def sample(self, state, offset, x, r, side=0):
        row = self.values(state, offset, [[x, 0, r]], side=side)
        return row[3][0], row[4][0]


class Profile:
    """Independent dimensional wall primitives and moving discrete normal rule."""
    def __init__(self, oracle, K, chi, relative_step):
        self.o, self.K, self.chi, self.relative_step = oracle, CD(K), CD(chi), LD(relative_step)
        self.nu, self.H, self.rough, self.rK, self.rm, _, self.Q, self.kappa = map(LD, oracle.raw['parameters'])
        root = np.sqrt(self.K)
        if self.rough:
            x = self.rough*root/(70*self.nu)
            exponent = maximum(LD('.5')-LD('.4')*x**LD('.7'), 1-LD('.79')*x**LD('-.28'))
            intercept = -LD('10.7')*np.expm1(exponent*np.log(x))
        else: intercept = CD('10.7')
        self.a, self.zv, self.ze = LD('.2295')*root, self.nu*intercept/root, LD('5.1')*self.nu/root
        self.amplitude = (self.K/self.Q)**LD('1.5')/self.kappa
        self.rc = CD(0)
        if np.real(self.chi) > 0:
            self.rc = CD(oracle.cutoff_seed(float(np.real(self.K))))
            for _ in range(4):
                slope = self.cutoff_slope(self.rc, self.relative_step)
                require(np.isfinite(slope) and abs(slope) > 0, 'Independent discrete cutoff derivative is singular')
                self.rc -= self.budget(self.rc)/slope
            require(np.real(self.rc) > 0 and np.real(self.rc) <= self.rm*(1+1e-10), 'Independent cutoff left its admissible interval')
            one, two = self.cutoff_slope(self.rc, LD('1e-3')), self.cutoff_slope(self.rc, LD('1e-4'))
            oracle.statistics['maximum_cutoff_slope_ladder_relative_difference'] = max(
                oracle.statistics['maximum_cutoff_slope_ladder_relative_difference'], float(abs(one-two)/max(LD(1), abs(two))))
            oracle.statistics['maximum_discrete_cutoff_budget_residual'] = max(
                oracle.statistics['maximum_discrete_cutoff_budget_residual'], float(abs(self.budget(self.rc))))

    def viscosity(self, r): return self.a*positive(r-self.zv)
    def original(self, r): return self.K**LD('1.5')/(LD('2.55')*maximum(r, self.ze))
    def reference_production(self, r):
        nt = self.viscosity(r)
        return nt*(self.K/self.Q)**2/(self.nu+nt)**2
    def corrected(self, r, cutoff): return self.amplitude/maximum(r, cutoff)

    def rule(self, cutoff):
        breaks = [CD(0), CD(self.H)]
        for value in [*self.o.inner_breaks, self.rough, self.zv, self.ze, self.rm, cutoff]:
            if 0 < np.real(value) < self.H: breaks.append(CD(rounded_base(value)))
        seed = minimum(CD(self.rm), self.nu/self.a)
        if np.real(cutoff) > 0: seed = minimum(seed, cutoff)
        for _ in range(4096):
            if np.real(seed) >= self.H: break
            if np.real(seed) > 0: breaks.append(CD(rounded_base(seed)))
            seed *= 2
        else: raise ValueError('Independent geometric panel sequence did not terminate')
        # Stable base ordering fixes one selected branch at coincident moving
        # breakpoints. Complex coordinates/weights retain every displacement.
        breaks.sort(key=lambda value:float(np.real(value)))
        unique = []
        for value in breaks:
            if not unique or np.real(value) != np.real(unique[-1]): unique.append(value)
        z, w = gauss(int(self.o.raw['quadrature_orders'][1]))
        points, weights = [], []
        for lo, hi in zip(unique[:-1], unique[1:]):
            half = (hi-lo)/2
            points.extend(rounded_base(lo+(1+z)*half)); weights.extend(rounded_base(w*half))
        return np.asarray(points, dtype=CD), np.asarray(weights, dtype=CD)

    def budget(self, cutoff):
        r, w = self.rule(cutoff)
        return np.sum(w*(self.corrected(r, cutoff)-self.reference_production(r)), dtype=CD)

    def cutoff_slope(self, cutoff, relative):
        h = max(abs(np.real(cutoff))*relative, LD('1e-12'))
        return (self.budget(cutoff-2*h)-8*self.budget(cutoff-h)+8*self.budget(cutoff+h)-self.budget(cutoff+2*h))/(12*h)

    def admissibility(self):
        r, w = self.rule(CD(0))
        p = np.sum(w*self.reference_production(r), dtype=CD)
        minimum_budget = np.sum(w*self.corrected(r, CD(self.rm)), dtype=CD)
        return (p-minimum_budget)/minimum_budget, self.H*np.sqrt(self.K)/self.nu


def gate(value, lower, upper):
    if np.real(value) <= lower: return CD(1)
    if np.real(value) >= upper: return CD(0)
    t = (value-lower)/(upper-lower)
    return 1-10*t**3+15*t**4-6*t**5


def inner_tensor(raw, G, k, nt):
    """Indexed value law; no derivative API and no conjugate contractions."""
    S = G+G.T; c = 2*k/3; square = np.sum(S*S, dtype=CD)
    if raw['inner_stress_model'] == 'Linear': F = nt*S
    elif np.real(k) == 0:
        # Selected admissible k->0+ limit at nonzero strain. The joint
        # k=S=0 directional-only corner has no single complex Frechet action.
        require(np.real(square) > 0, 'Joint zero-k/zero-strain corner requires a separate admissible one-sided value check')
        F = c*S/np.sqrt(square)
    elif np.real(square) == 0 or np.real(nt*np.sqrt(square)) <= LD('.8')*np.real(c): F = nt*S
    else:
        norm = np.sqrt(square); a = nt*norm; z = (a-LD('.8')*c)/(LD('.2')*c)
        F = c*(1-LD('.2')*np.exp(-z))*S/norm
    R = c*np.eye(3, dtype=CD)-F
    return F, -np.sum(R*G, dtype=CD)


def outer_tensor(raw, G, k, epsilon, r, branch, yap_branch):
    """Independent dimensional indexed CLS + prescribed beta/Gamma value."""
    S, O = G+G.T, G-G.T; eye = np.eye(3, dtype=CD)
    ss, oo = np.sum(S*S, dtype=CD), np.sum(O*O, dtype=CD)
    selected = ss if branch == 'S' else oo
    q = CD(0) if np.real(selected) == 0 else k/epsilon*np.sqrt(selected/2)
    factor = CD(1) if np.real(q) > 8 else -np.expm1(-LD('.36')*np.exp(LD('.75')*q))
    c = LD('.3')*factor/(1+LD('.35')*q**LD('1.5')); nt = c*k*k/epsilon
    t1 = np.einsum('ik,jk->ij', S, S)-ss*eye/3
    t2 = np.einsum('ik,jk->ij', O, S)+np.einsum('jk,ik->ij', O, S)
    t3 = np.einsum('ik,jk->ij', O, O)-oo*eye/3
    t4 = np.einsum('ki,lj,kl->ij', S, O, S)+np.einsum('kj,li,kl->ij', S, O, S)
    quadratic = k*nt/epsilon*(-LD('.1')*t1+LD('.1')*t2+LD('.26')*t3)
    cubic = nt*k*k/epsilon**2*c*c*(-10*t4+5*(oo-ss)*S)
    R = 2*k*eye/3-nt*S+quadratic+cubic
    if raw['model_variant'] != 'Original':
        beta, uniform, _, _, _, scale = map(LD, raw['robust_parameters'][:6])
        flow, normal = np.asarray(raw['robust_reference_flow'], dtype=LD), np.asarray(raw['robust_reference_normal'], dtype=LD)
        ref = LD(raw['ustar'])/(LD(raw['parameters'][7])*r)*(np.outer(flow, normal)+np.outer(normal, flow))
        R += (beta-1)*cubic-(uniform+scale*LD(raw['parameters'][7])*LD(raw['ustar'])*r)*(S-ref)
    P = -np.sum(R*G, dtype=CD)
    ratio = k**LD('1.5')/(LD('2.5')*epsilon*r)
    yap = LD('.83')*epsilon**2/k*((ratio-1)*ratio**2 if yap_branch == 'A' else CD(0))
    return R, nt, P, LD('1.44')*epsilon/k*P-LD('1.92')*epsilon**2/k+yap, ss-oo, ratio-1


class ValueOracle:
    def __init__(self, raw, cutoff_relative_step='1e-4'):
        require(raw['dimension'] == 2 and raw['mapped'] is False and raw['upstream_transition'] is True,
            'Independent R3 value oracle requires the declared flat2D Upstream scope')
        require(np.finfo(LD).eps < LD('1e-18'), 'Authoritative long-double precision is required')
        require('chi' in raw, 'Native probe must export requested chi explicitly')
        self.raw, self.off = raw, np.asarray(raw['offsets'], int)
        self.fields, self.real_fields = [Field(row) for row in raw['fields']], real_spaces(raw)
        self.mass = [field.mass() for field in self.real_fields]
        self.H, self.nu = LD(raw['parameters'][1]), LD(raw['parameters'][0])
        self.state = np.asarray(raw['state'], dtype=CD); self.reference = np.asarray(raw['transition_reference'], dtype=CD)
        self.relative = LD(cutoff_relative_step)
        self.trace_ids = {int(index) for e in self.fields[3].elements if abs(e['lower'][2]-float(self.H)) < 1e-11
            for index in e['ids'][:e['degree']+1]}
        self.external = {int(row):CD(value) for row,value in raw['constraints'] if
            not (self.off[3] <= row < self.off[4] and int(row)-self.off[3] in self.trace_ids)}
        self.inner_breaks = sorted({float(e[key][2]) for e in self.fields[4].elements for key in ('lower', 'upper')})
        self.columns = sorted({(float(e['lower'][0]), float(e['upper'][0])) for e in self.fields[4].elements})
        metadata = [row for row in raw['assembly_metadata'] if 'outer_branch_pattern' in row]
        require(len(metadata) == 1, 'Missing unique native generalized outer-branch binding')
        self.outer_pattern = metadata[0]['outer_branch_pattern']; self.yap_pattern = metadata[0]['outer_Yap_branch_pattern']
        require(len(self.outer_pattern) == len(self.yap_pattern) and set(self.outer_pattern) <= {'S','R'} and
            set(self.yap_pattern) <= {'A','I'}, 'Incompatible declared outer branch tables')
        self.statistics = dict(maximum_cutoff_slope_ladder_relative_difference=0.,
            maximum_discrete_cutoff_budget_residual=0., minimum_outer_invariant_switch_distance=float('inf'),
            minimum_Yap_switch_distance=float('inf'), outer_branch_ambiguities=0, evaluations=0,
            source_limiter_active_points=0, cap_active_rays=0, active_state_gate_evaluations=0,
            active_admissibility_gate_evaluations=0, minimum_normalization_lambda=1.,
            trace_projection_calls=0)

    @lru_cache(maxsize=512)
    def cutoff_seed(self, K):
        nu,H,h,rK,rm,_,Q,kappa = self.raw['parameters']
        return float(Profiles(dict(geometry=[nu,H,h,rK,rm], metric=[0.,0.], k=K, q=Q, kappa=kappa, chi=1.)).rc)

    def sample(self, state, block, x, r, side=0):
        return self.fields[block].sample(state, self.off[block], x, r, side)

    def forcing(self, state, x):
        U, gradient = self.sample(state, 0, x, self.H, 1)
        _, pressure = self.sample(state, 5, x, self.H, 1)
        _, gk = self.sample(state, 2, x, self.H, 1)
        return pressure[0]-2*gk[0]/3, 2*U*gradient[0]

    def chi(self, state, x, K, baseline):
        if baseline or self.raw['chi'] == 0: return CD(0)
        start,end,ra,rb,_,on,off,moff,mon,roff,ron = map(LD, self.raw['transition_controls'])
        if x >= end: return CD(0)
        profile = Profile(self, K, 0, self.relative); margin,reynolds = profile.admissibility()
        if np.real(margin) <= moff or np.real(reynolds) <= roff: return CD(0)
        scale2 = K/LD(self.raw['parameters'][6])+(self.nu/self.H)**2; scale = np.sqrt(scale2)
        samples = {}
        for block in (0, 1):
            for r in (self.H, ra, rb):
                samples[(block,r)] = self.sample(state, block, x, r, 1 if r == self.H else 0)[0]-self.sample(
                    self.reference, block, x, r, 1 if r == self.H else 0)[0]
        ka = self.sample(state, 2, x, ra)[0]; ea = self.sample(state, 3, x, ra)[0]
        ea0 = self.sample(self.reference, 3, x, ra)[0]
        pressure,convection = self.forcing(state, x); p0,c0 = self.forcing(self.reference, x)
        groups = [samples[(0,self.H)]**2/scale2,
            (samples[(0,rb)]-samples[(0,ra)])**2/scale2, samples[(1,ra)]**2/scale2,
            ((ka-K)/K)**2, ((ea-ea0)/ea0)**2,
            (self.H*(pressure-p0)/scale2)**2, (self.H*(convection-c0)/scale2)**2]
        departure = sum(LD(weight)*value for weight,value in zip(self.raw['transition_weights'], groups))
        state_gate=gate(departure,on*on,off*off)
        admissibility=(1-gate(margin,moff,mon))*(1-gate(reynolds,roff,ron))
        self.statistics['active_state_gate_evaluations'] += int(0 < np.real(state_gate) < 1)
        self.statistics['active_admissibility_gate_evaluations'] += int(0 < np.real(admissibility) < 1)
        return gate(CD(x),start,end)*state_gate*admissibility

    def ray(self, state, x, baseline=False):
        key = (float(x), bool(baseline))
        if key in self.cache: return self.cache[key]
        K = self.sample(state, 2, x, LD(self.raw['parameters'][3]))[0]
        require(np.real(K) > 0, 'Independent sampled K left its positive domain')
        chi = self.chi(state,x,K,baseline); profile = Profile(self,K,chi,self.relative)
        r,w = profile.rule(profile.rc if np.real(chi) > 0 else CD(0))
        points = np.column_stack([np.full_like(r,x), np.zeros_like(r), r])
        k = self.fields[2].values(state,self.off[2],points)[3]
        ea,p0 = profile.original(r), profile.reference_production(r)
        ec = profile.corrected(r,profile.rc) if np.real(chi)>0 else ea
        pre_e,pre_d = (1-chi)*ea+chi*ec, (1-chi)*ea+chi*p0
        theta,delta,smooth = (LD(value) for value in self.raw['robust_parameters'][2:5])
        ratio = k/(theta*K); limited = np.where(np.real(ratio) >= 1, CD(1), ratio)
        require(np.min(np.real(k)) >= 0, 'Independent source k left its nonnegative domain')
        limiter = 10*limited**3-15*limited**4+6*limited**5 if smooth else limited
        budget = np.sum(w*pre_d,dtype=CD); b = np.sum(w*pre_d*limiter,dtype=CD)/budget
        m = b if np.real(b) >= delta else CD(delta); lam = b/m
        self.statistics['source_limiter_active_points'] += int(np.count_nonzero(np.real(ratio) < 1))
        self.statistics['cap_active_rays'] += int(np.real(b) < delta)
        self.statistics['minimum_normalization_lambda'] = min(self.statistics['minimum_normalization_lambda'],float(np.real(lam)))
        E,D = lam*pre_e, pre_d*limiter/m
        answer = dict(profile=profile, r=r, w=w, points=points, E=E, D=D, lam=lam,
            endpoint=lam*((1-chi)*profile.original(CD(self.H))+
                chi*(profile.corrected(CD(self.H),profile.rc) if np.real(chi)>0 else profile.original(CD(self.H)))))
        self.cache[key] = answer; return answer

    def trace(self, state):
        self.statistics['trace_projection_calls'] += 1
        baseline,change,baseline_ids = {},{},set(); start = self.raw['transition_controls'][4]
        for e in self.fields[3].elements:
            if abs(e['lower'][2]-float(self.H)) > 1e-11: continue
            degree = e['degree']; nodes,transform = _segment_projection(degree)
            x = LD(e['lower'][0])+(np.asarray(nodes,dtype=LD)+1)*(LD(e['upper'][0])-LD(e['lower'][0]))/2
            values = np.asarray([[self.ray(state,xx,True)['endpoint'], self.ray(state,xx,False)['endpoint']-
                self.ray(state,xx,True)['endpoint']] for xx in x], dtype=CD)
            coefficients = np.asarray(transform,dtype=LD)@values
            for j,index in enumerate(e['ids'][:degree+1]):
                baseline.setdefault(int(index),[]).append(LD(e['sign'][j])*coefficients[j,0])
                change.setdefault(int(index),[]).append(LD(e['sign'][j])*coefficients[j,1])
            if e['lower'][0] >= start-1e-11: baseline_ids.update(map(int,e['ids'][:degree+1]))
        return {index:sum(values)/len(values)+(0 if index in baseline_ids else sum(change[index])/len(change[index]))
            for index,values in baseline.items()}

    def value(self, state, history=None, alpha=None):
        self.statistics['evaluations'] += 1; self.cache = {}
        original = np.asarray(state,dtype=CD); state = original.copy()
        for row,value in self.external.items(): state[row] = value
        result = np.zeros(len(state),dtype=CD); off = self.off
        def add(block, row, weights, source, flux):
            ids,phi,gradient = row[:3]
            load = np.sum(phi*(weights*source)[:,None],axis=0,dtype=CD)
            load += np.sum(gradient*(weights[:,None]*flux)[:,:,None],axis=(0,1),dtype=CD)
            np.add.at(result,off[block]+ids,load)
        def integrate(identity, points, weights, inner, ray=None, indices=None, branches=None):
            blocks = (0,1,2,4) if inner else (0,1,2,3,5)
            rows = {block:self.fields[block].values(state,off[block],points,identity) for block in blocks}
            U,G = np.zeros((len(points),3),dtype=CD),np.zeros((len(points),3,3),dtype=CD)
            for block,axis in ((0,0),(1,2)): U[:,axis],G[:,axis] = rows[block][3],rows[block][4]
            k,gk = rows[2][3],rows[2][4]; pi = rows[4 if inner else 5][3]
            stress,P,nt,se = [],[],[],[]
            if inner:
                viscosity = ray['profile'].viscosity(points[:,2])
                for g,kval,n in zip(G,k,viscosity):
                    F,production = inner_tensor(self.raw,g,kval,n)
                    stress.append(self.nu*(g+g.T)+F); P.append(production)
                nt = viscosity; sk = np.asarray(P)-ray['D'][indices]; diffusion = self.nu+nt
                force = sum(self.forcing(state,float(np.real(points[0,0]))))
                body = np.zeros_like(U); body[:,0] = np.where(np.real(points[:,2]) <= self.raw['parameters'][2],force,0)
            else:
                epsilon,ge = rows[3][3],rows[3][4]
                for j,(g,kval,e,point) in enumerate(zip(G,k,epsilon,points)):
                    branch,yap = branches[j]
                    R,n,production,source,gap,ygap = outer_tensor(self.raw,g,kval,e,point[2],branch,yap)
                    magnitude = max(LD(1),abs(np.sum((g+g.T)**2)),abs(np.sum((g-g.T)**2)))
                    ambiguous = abs(np.real(gap)) <= 128*np.finfo(float).eps*magnitude
                    require(ambiguous or (np.real(gap) >= 0) == (branch == 'S'), 'Independent/native invariant branches differ away from a tie')
                    require(abs(np.real(ygap)) <= 128*np.finfo(float).eps or (np.real(ygap)>0) == (yap=='A'),
                        'Independent/native Yap branches differ away from a knot')
                    self.statistics['outer_branch_ambiguities'] += int(ambiguous)
                    self.statistics['minimum_outer_invariant_switch_distance'] = min(self.statistics['minimum_outer_invariant_switch_distance'],float(abs(np.real(gap))))
                    self.statistics['minimum_Yap_switch_distance'] = min(self.statistics['minimum_Yap_switch_distance'],float(abs(np.real(ygap))))
                    stress.append(-R+2*kval*np.eye(3,dtype=CD)/3); nt.append(n); P.append(production);se.append(source)
                nt = np.asarray(nt,dtype=CD); sk=np.asarray(P)-epsilon;diffusion=nt;body=np.zeros_like(U)
                add(3,rows[3],weights,-np.asarray(se),nt[:,None]/LD('1.3')*ge-U*epsilon[:,None])
            stress=np.asarray(stress,dtype=CD)-pi[:,None,None]*np.eye(3,dtype=CD)
            for block,axis in ((0,0),(1,2)):
                add(block,rows[block],weights,-body[:,axis],stress[:,axis]-U[:,axis,None]*U)
            add(2,rows[2],weights,-sk,diffusion[:,None]*gk-U*k[:,None])
            pressure=4 if inner else 5
            add(pressure,rows[pressure],weights,np.trace(G,axis1=1,axis2=2),np.zeros_like(U))
        z,w = gauss(int(self.raw['quadrature_orders'][0]))
        for x0,x1 in self.columns:
            for xx,ww in zip(LD(x0)+(1+z)*LD(x1-x0)/2,w*LD(x1-x0)/2):
                ray=self.ray(state,xx); groups={}
                for index,point in enumerate(ray['points']): groups.setdefault(self.fields[0].element(point)['geometry'],[]).append(index)
                for identity,ids in groups.items():
                    ids=np.asarray(ids,int);integrate(identity,ray['points'][ids],ww*ray['w'][ids],True,ray,ids)
        cursor=0
        for e in self.fields[5].elements:
            points,weights=[],[]
            # Native outer quadrature enumerates physical x fastest. Only the
            # declared tie labels depend on this ordering; value rows are ours.
            for rz,wr in zip(z,w):
                for xx,wx in zip(z,w):
                    points.append([LD(e['lower'][0])+(xx+1)*LD(e['upper'][0]-e['lower'][0])/2,0,
                        LD(e['lower'][2])+(rz+1)*LD(e['upper'][2]-e['lower'][2])/2])
                    weights.append(wx*wr*LD(e['upper'][0]-e['lower'][0])*LD(e['upper'][2]-e['lower'][2])/4)
            count=len(points); require(cursor+count <= len(self.outer_pattern),'Incomplete native outer branch inventory')
            branches=list(zip(self.outer_pattern[cursor:cursor+count],self.yap_pattern[cursor:cursor+count]));cursor+=count
            integrate(e['geometry'],np.asarray(points,dtype=CD),np.asarray(weights,dtype=CD),False,branches=branches)
        require(cursor == len(self.outer_pattern),'Unconsumed native outer branch inventory')
        volume=LD(self.real_fields[0].volume)
        for block in (4,5):
            mean=np.asarray(self.mass[block]@self.real_fields[block].constant(),dtype=LD)/volume
            ids=np.arange(self.fields[block].size)+off[block]
            result[ids]+=state[-1]*mean;result[-1]+=np.sum(state[ids]*mean,dtype=CD)
        if history is not None:
            for block in range(4):
                ids=np.arange(self.fields[block].size)+off[block]
                result[ids]+=LD(alpha)*(np.asarray(self.mass[block],dtype=LD)@(state[ids]-np.asarray(history,dtype=CD)[ids]))
        for row,target in self.external.items(): result[row]=original[row]-target
        trace=self.trace(state)
        for index,target in trace.items(): result[off[3]+index]=original[off[3]+index]-target
        require(np.all(np.isfinite(result)), 'Nonfinite independent complex value loads')
        return result
