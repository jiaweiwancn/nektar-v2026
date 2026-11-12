#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent physical L2 and coefficient-trace reconstruction, never a solve
of a governing equation. All b-d stages are reconstructions, not native exports.
"""
import ast
from functools import lru_cache
import xml.etree.ElementTree as ET
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from numpy.polynomial.legendre import Legendre
from validate_cls_epawf_m3 import basis
from wj_h4_points import mapping,vector_hash
from wj_h3_analysis import values
from wj_h5_neutral import QE,KAPPA


def expression(text):
    tree=ast.parse(text.replace('^','**'),mode='eval')
    allowed=(ast.Expression,ast.BinOp,ast.UnaryOp,ast.Compare,ast.Call,ast.Name,ast.Load,ast.Constant,
        ast.Add,ast.Sub,ast.Mult,ast.Div,ast.Pow,ast.USub,ast.UAdd,ast.Lt,ast.LtE,ast.Gt,ast.GtE,ast.Eq,ast.NotEq)
    for node in ast.walk(tree):
        if not isinstance(node,allowed):raise ValueError('unsupported XML expression node')
        if isinstance(node,ast.Name)and node.id not in ('x','y','z','t','log','exp','sqrt'):raise ValueError('unsupported XML expression name')
        if isinstance(node,ast.Call)and not(isinstance(node.func,ast.Name)and node.func.id in ('log','exp','sqrt')):raise ValueError('unsupported XML call')
    code=compile(tree,'<frozen XML Function>','eval')
    def evaluate(x,y):
        x,y=np.broadcast_arrays(x,y)
        v=eval(code,{'__builtins__':{}},dict(x=x,y=y,z=0.,t=0.,log=np.log,exp=np.exp,sqrt=np.sqrt))
        v=np.broadcast_to(np.asarray(v,dtype=float),x.shape)
        if not np.all(np.isfinite(v)):raise ValueError('nonfinite XML function evaluation')
        return v
    return evaluate


class Functions:
    def __init__(self,xml,report):
        self.H,self.d0,self.top,self.J=mapping(report);self.functions={}
        for f in ET.parse(xml).findall('./CONDITIONS/FUNCTION'):
            for e in f.findall('E'):self.functions[(f.get('NAME'),e.get('VAR'))]=expression(e.get('VALUE'))

    def raw(self,name,variable,x,rho):
        x,rho=np.broadcast_arrays(x,rho);z=self.d0+np.where(rho<=self.H,rho,self.H+(rho-self.H)*self.J)
        return self.functions[(name,variable)](x,z)

    def block(self,b,x,rho,family='InitialConditions'):
        x,rho=np.broadcast_arrays(x,rho);variable=('u','v','k','epsilon','p','p')[b]
        if b in (3,5):return self.raw(family+'Outer',variable,x,rho)
        if b==4:return self.raw(family+'Inner',variable,x,rho)
        mask=rho<=self.H if family=='InitialConditions'else rho<self.H
        v=np.empty_like(rho,dtype=float)
        for take,side in ((mask,'Inner'),(~mask,'Outer')):
            if np.any(take):v[take]=self.raw(family+side,variable,x[take],rho[take])
        return v


def local_rule(el,order=16):
    g,w=np.polynomial.legendre.leggauss(order);lo,hi=el['lower'],el['upper']
    x,r=np.meshgrid((lo[0]+hi[0])/2+(hi[0]-lo[0])/2*g,(lo[2]+hi[2])/2+(hi[2]-lo[2])/2*g)
    p=basis(el['degree'],g)[0]
    B=np.einsum('ni,mj->nmij',p,p).reshape(order*order,-1)*np.asarray(el['sign'])
    weight=(np.outer(w,w)*(hi[0]-lo[0])*(hi[2]-lo[2])/4).ravel()
    return x.ravel(),r.ravel(),B,weight


def project_field(report,b,function,order=16):
    """Assemble physical mass/load equations directly from exported CG maps."""
    space=report['fields'][b];n=space['size'];M=np.zeros((n,n));load=np.zeros(n)
    H,d0,top,J=mapping(report);panels=[]
    for el in space['mesh']:
        x,r,B,w=local_rule(el,order);w=w*(1 if el['upper'][2]<=H else J)
        ids=np.asarray(el['ids'],int);M[np.ix_(ids,ids)]+=B.T@(w[:,None]*B)
        f=function(x,r);load[ids]+=B.T@(w*f);panels.append((ids,B,w,f))
    factor=cho_factor(M,check_finite=False);c=cho_solve(factor,load,check_finite=False)
    # Refine the *integral* defect, avoiding cancellation of assembled M*c and
    # load on the hierarchical basis. This is offline projection algebra, not
    # a velocity correction or governing solve. Equations/quadrature stay fixed.
    for _ in range(2):
        correction=np.zeros(n,dtype=np.longdouble)
        for ids,B,w,f in panels:
            defect=np.asarray(f,dtype=np.longdouble)-np.asarray(B,dtype=np.longdouble)@np.asarray(c[ids],dtype=np.longdouble)
            correction[ids]+=np.asarray(B,dtype=np.longdouble).T@(np.asarray(w,dtype=np.longdouble)*defect)
        c+=cho_solve(factor,np.asarray(correction,dtype=float),check_finite=False)
    defect=M@c-load
    backward=float(np.max(abs(defect))/max(np.linalg.norm(M,np.inf)*np.linalg.norm(c,np.inf)+np.linalg.norm(load,np.inf),1e-300))
    return c,dict(order=order,backward_residual=backward,load_max=float(np.max(abs(load))),
        residual_max=float(np.max(abs(defect))),physical_mass_shape=list(M.shape),integral_defect_refinements=2)


@lru_cache(None)
def gll(n):
    p=Legendre.basis(n-1);g=np.r_[-1.,p.deriv().roots(),1.];return g,2/(n*(n-1)*p(g)**2)


def constrained_segment(function,degree=3):
    g,w=gll(max(8,degree+3));B=basis(degree,g)[0];v=np.asarray(function(g))
    c=np.zeros(degree+1);c[:2]=[v[0],v[-1]]
    if degree>1:c[2:]=np.linalg.solve(B[:,2:].T@(w[:,None]*B[:,2:]),B[:,2:].T@(w*(v-B@c)))
    return c


def boundary(report,b,function,predicate):
    result={}
    for el in report['fields'][b]['mesh']:
        lo,hi=el['lower'],el['upper'];p=el['degree'];m=p+1
        for axis in (0,2):
            for side in (0,1):
                if not predicate(axis,side,el):continue
                free=2 if axis==0 else 0
                def sample(g):
                    a=np.full_like(g,hi[axis]if side else lo[axis]);z=(lo[free]+hi[free])/2+g*(hi[free]-lo[free])/2
                    return function(a,z)if axis==0 else function(z,a)
                coefficients=constrained_segment(sample,p)
                for j,c in enumerate(coefficients):
                    tensor=side+m*j if axis==0 else j+m*side;index=el['ids'][tensor];value=el['sign'][tensor]*c
                    if index in result and abs(result[index]-value)>1e-10:raise ValueError('incompatible reconstructed boundary corner')
                    result[index]=float(value)
    return result


def reconstruct(report,xml):
    assert report['steps']==0 and report['implicit_callback_count']==0 and not report['pimple_attempts']
    functions=Functions(xml,report);H=functions.H;top=functions.top
    stages={};projection=np.zeros(len(report['state']));systems=[];offset=list(map(int,report['offsets']))
    for b in range(6):
        c,meta=project_field(report,b,lambda x,r,b=b:functions.block(b,x,r));projection[offset[b]:offset[b+1]]=c;systems.append(meta)
    stages['physical_L2']=projection
    trace_ids={e['ids'][j]for e in report['fields'][3]['mesh']if e['lower'][2]==H for j in range(e['degree']+1)}
    exterior=projection.copy();constraints={}
    for b in range(4):
        def fn(x,r,b=b):
            value=functions.block(b,x,r,'BoundaryValues')
            if b<2:value=np.where(abs(r)<1e-11,0.,value)
            if b==1:value=np.where(abs(r-top)<1e-11,0.,value)
            return value
        def predicate(axis,side,el,b=b):
            if axis==0:return not side and el['lower'][0]==0
            return b<2 and ((not side and el['lower'][2]==0)or(side and b==1 and el['upper'][2]==top))
        for index,v in boundary(report,b,fn,predicate).items():
            if b==3 and index in trace_ids:continue
            row=offset[b]+index;constraints[row]=v;exterior[row]=v
    stages['exterior_lift']=exterior
    def epsilon_target(x,r):
        matching=next(e for e in report['fields'][2]['mesh']if e['lower'][2]==H)
        K=values(report,exterior,2,matching['geometry'],x,np.full_like(x,H))[0]
        if np.any(K<=0):raise ValueError('nonpositive represented endpoint K')
        return (K/QE)**1.5/(KAPPA*H)
    eps=boundary(report,3,epsilon_target,lambda axis,side,e:axis==2 and not side and e['lower'][2]==H)
    final=exterior.copy()
    for index,v in eps.items():final[offset[3]+index]=v
    stages['dependent_epsilon_lift']=final
    scales=np.array([1.,1.,report['parameters'][5],(report['parameters'][5]/QE)**1.5/(KAPPA*H),report['ustar']**2,report['ustar']**2])
    errors=[];native=np.asarray(report['state'])
    for b in range(6):
        a,z=offset[b:b+2];error=float(np.max(abs(final[a:z]-native[a:z])))
        errors.append(dict(block=b,absolute_max=error,scale=float(scales[b]),scaled=float(error/scales[b]),passed=bool(error<=2e-11*scales[b]+1e-15)))
    lineage='uncommitted native startup candidate'if report['failure']else 'accepted native CP0'
    if report['failure']:
        ev=report['physical_range_events'][-1];assert vector_hash(ev['state'])==vector_hash(native)and ev['alpha']==0 and ev['time']==0
        assert not report['accepted_states'] and report['checkpoint']==''
    else:assert len(report['accepted_states'])==1 and vector_hash(report['accepted_states'][0])==vector_hash(native)
    return stages,dict(status='PASS'if all(v['passed']for v in errors)else 'UNRESOLVED',coefficient_errors=errors,
        physical_linear_systems=systems,exterior_constraints=constraints,dependent_epsilon_constraints=eps,
        native_lineage=lineage,stage_scope='independent reconstructed L2/exterior/epsilon stages, not native intermediate exports',
        native_state_sha256=vector_hash(native),reconstructed_final_sha256=vector_hash(final))
