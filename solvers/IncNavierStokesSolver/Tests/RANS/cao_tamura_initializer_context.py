# SPDX-License-Identifier: MIT
"""Named, versioned binary64 native initializer context. V1 stays read-only."""
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import sha,require
from search_cao_tamura_budget_h2e import condition

class NativeContext:
    def __init__(self,path,verify=True):
        self.path=Path(path);self.doc=load(self.path);self.meta=self.doc['context'];self.root=Path(self.doc['array_directory']);self.loaded={}
        require(self.meta['version'] in ('NativeOuterScalarContextV1','NativeOuterScalarContextV2'),'unknown initializer schema')
        self.counts=self.meta.get('direction_counts',dict(velocity=32,k=8,epsilon=8));nv=self.counts['velocity'];nk=self.counts['k'];ne=self.counts['epsilon']
        self.blocks=self.meta.get('observation_blocks',{})
        if not self.blocks:
            # The single compatibility adapter for historical V1. New consumers
            # never encode a direction-count-dependent numeric slice.
            for name,offset,width,shape in [('chart',0,3,[3]),('physical',3,3,[3]),('gradient',6,4,[2,2]),('velocity',10,2,[2]),('velocity_basis',18,nv*6,[nv,6]),('k_basis',18+nv*6,nk,[nk]),('epsilon_basis',18+nv*6+nk,ne,[ne])]:self.blocks[name]=dict(offset=offset,width=width,shape=shape)
            for name,offset in [('k',12),('epsilon',13),('inner_nu_t',14),('inner',15),('covariance',16),('dependent_epsilon',17)]:self.blocks[name]=dict(offset=offset,width=1,shape=[1])
        require(self.blocks['velocity_basis']['shape']==[nv,6] and self.blocks['k_basis']['shape']==[nk] and self.blocks['epsilon_basis']['shape']==[ne],'schema direction counts disagree with observation blocks')
        require(self.meta['files']['basis']['shape'][0]==nv+nk+ne and self.meta['files']['mass_gram']['shape']==[nv,nv] and self.meta['files']['scalar_gram']['shape']==[nk+ne,nk+ne],'schema matrix/direction dimensions disagree')
        if self.meta['version']=='NativeOuterScalarContextV2':
            for key in ('geometry_identity','inventory_identity','reference_sha256'):require(self.meta[key]==self.doc[key],'context identity mismatch '+key)
        self.verify=verify;self.observations=self.array('observations');self.N=len(self.observations)
        intervals=[]
        for b in self.blocks.values():intervals.extend(range(b['offset'],b['offset']+b['width']))
        require(sorted(intervals)==list(range(self.observations.shape[1])),'overlapping/incomplete observation blocks')
    def array(self,name):
        if name not in self.loaded:
            d=self.meta['files'][name];p=self.root/d['file']
            if self.verify:require(sha(p)==d['sha256'],'context checksum mismatch '+name)
            self.loaded[name]=np.memmap(p,dtype=d['dtype'],mode='r',shape=tuple(d['shape']))
        return self.loaded[name]
    def block(self,name):
        b=self.blocks[name];a=self.observations[:,b['offset']:b['offset']+b['width']]
        return a.reshape(self.N) if name in ('k','epsilon','inner_nu_t','inner','covariance','dependent_epsilon') else a.reshape((self.N,)+tuple(b['shape']))
    def gradient(self):
        g=np.zeros((self.N,3,3));a=self.block('gradient')
        for i,x in enumerate((0,2)):
            for j,z in enumerate((0,2)):g[:,x,z]=a[:,i,j]
        return g
    def transforms(self):
        G=np.array(self.array('mass_gram'));nv=self.counts['velocity'];old=min(32,nv);T0,vc=condition(G[:old,:old])
        if nv>old:
            # New directions have already been mass-orthogonalized natively.
            # Keep original coordinates rather than rotating a degenerate I.
            T=np.zeros((nv,nv));T[:old,:old]=T0;T[old:,old:]=np.eye(nv-old)
            error=float(np.max(abs(T.T@G@T-np.eye(nv))));require(error<1e-8,'enriched mass Gram is not conditioned');vc.update(rank=nv,original_rank=old,enriched_orthogonality_error=error)
        else:T=T0
        sg=np.array(self.array('scalar_gram'));nk=self.counts['k'];ne=self.counts['epsilon'];Tk,kc=condition(sg[:nk,:nk]);Te,ec=condition(sg[nk:,nk:]);Ts=np.zeros((nk+ne,Tk.shape[1]+Te.shape[1]));Ts[:nk,:Tk.shape[1]]=Tk;Ts[nk:,Tk.shape[1]:]=Te
        return T,Ts,dict(velocity=vc,k=kc,epsilon=ec)
